#!/usr/bin/env python3
"""Explicit outcome, metrics, runtime split and criteria dry run of ONE muon of run_k2a_validation.sh.

  k2a_validation_muon.py MUON_DIR --muon M --label REAL|SYNTHETIC --independent TEXT --t0 0|1 --criteria v1.json

MUON_DIR = OUT/muons/muon_M (steps.tsv, timing.tsv, far_dir.txt -> the far-mode experiment muon directory E).
Writes MUON_DIR/outcome.json; prints one line. Never raises on missing pieces: they make the outcome INCONCLUSIVE.

Outcome (first that applies):
  FAIL           numerical / integrity failure: archive changed; moment path differs (archive moments vs export
                 moments, s_compare_moments 1e-9; responses 1e-7); compiled K2A differs from the Python K2A_local
                 (1e-12) or only one of them refuses; a comparison reports numerical breakage (FAIL); the shaped
                 comparison could not be made on identical strip labels / clock; production evaluator error
  S0_FAILED      the established S0_T1 reference is not available (exporter, T1 collapse, grouping, a frozen replay
                 refusal or failure, census): cause named. K2A domain status is still reported
  OUT_OF_DOMAIN  K2A blocking statuses (named, with the blocked charge fraction); the S0 reference is still computed
                 and its availability reported; no accuracy metric is claimed
  INCONCLUSIVE   anything required missing or inconclusive
  VALID          S0_T1, S0_T0 (if enabled) and K2A all available, every consistency check passed; the accuracy is
                 REPORTED (metrics), not thresholded here. VALID_T1_ONLY: as VALID but S0_T0 disabled or failed (named)
Criteria dry run: the frozen sealed criteria v1 applied to this NON-SEALED event (reference S0_T1, ion + electron,
Tp 100 and 200 ns, both planes), labelled DIAGNOSTIC; also against S0_T0 (diagnostic only). Not a sealed verdict.
"""
import argparse
import csv
import json
import math
import re
from pathlib import Path

CMP_NAMES = ("S0T1_vs_K2A", "S0T0_vs_K2A", "S0T0_vs_S0T1")
SHAPED_KEYS = ["L2_rel", "L2_unshaped_current", "L2_unshaped_charge", "lead_index_changed", "lead_peak_rel",
               "lead_peak_time_ns", "charge_rel", "amplitude_rel", "second_over_lead", "centroid_strip",
               "strip_peak_maxabs_over_T0lead", "t_f0.05_ns", "t_f0.1_ns", "t_f0.5_ns", "mult_f0.05", "mult_f0.1",
               "hitflips_f0.05", "hitflips_f0.1"]


def steps_of(p, header):
    out = {}
    if not p.is_file():
        return out
    rows = [l.split("\t") for l in p.read_text().splitlines() if l.strip()]
    for r in rows[1:] if header else rows:
        r += [""] * (4 - len(r))
        out[r[0]] = {"status": r[1], "detail": r[2], "file": r[3]}
    return out


def timing_of(p):
    t = {}
    if p.is_file():
        for r in csv.DictReader(open(p), delimiter="\t"):
            try:
                t[r["label"]] = {"wall_s": float(r["wall_s"]), "peak_rss_kib": int(r["peak_rss_kib"]) if r["peak_rss_kib"] else None,
                                 "exit": int(r["exit"])}
            except (ValueError, KeyError):
                pass
    return t


def jload(p):
    try:
        return json.loads(Path(p).read_text())
    except (OSError, ValueError):
        return None


def chosen_wall(sg):
    """Sum over groups of the wall time of the chosen ladder attempt (the single-pass replay cost) and of all attempts."""
    chosen, total, n, steps = 0.0, 0.0, 0, {}
    for g in sorted(Path(sg).glob("g_*/timing_attempts.json")) if Path(sg).is_dir() else []:
        j = jload(g) or {}
        atts = j.get("attempts", [])
        tot = sum(a.get("wall_s") or 0.0 for a in atts)
        ch = [a for a in atts if a.get("late_step_ns") == j.get("chosen_late_step_ns")]
        chosen += (ch[0].get("wall_s") or 0.0) if ch else float("nan")
        total += tot; n += 1; steps[g.parent.name] = j.get("chosen_late_step_ns")
    return chosen, total, n, steps


def cmp_metrics(j):
    if not j:
        return None
    m = {"status": j.get("status"), "worst_L2_rel": j.get("worst_L2_rel")}
    for term, tj in (j.get("terms") or {}).items():
        for s in ("U", "V"):
            pj = tj.get(s) or {}
            for k in ("unshaped_L2_rel", "cumulative_L2_rel", "window_charge_rel_diff", "max_strip_charge_diff_over_lead",
                      "shaped_L2_rel_Tp100", "shaped_L2_rel_Tp200", "common_strips", "ref_charge_on_ref_only_strips_frac",
                      "cand_charge_on_cand_only_strips_frac"):
                if k in pj:
                    m[f"{term}.{s}.{k}"] = pj[k]
            if "lead_identical" in pj:
                m[f"{term}.{s}.lead_changed"] = int(not pj["lead_identical"])
        if "UV_charge_ratio_rel_diff" in tj:
            m[f"{term}.UV_charge_ratio_rel_diff"] = tj["UV_charge_ratio_rel_diff"]
    return m


def shaped_lines(d):
    """{Tp: {variant: line}} from per_muon_Tp*.jsonl of one t_shaped_compare run (one muon)."""
    out = {}
    for f in sorted(Path(d).glob("per_muon_Tp*.jsonl")) if Path(d).is_dir() else []:
        tp = f.stem.split("Tp")[1]
        for line in f.read_text().splitlines():
            if line.strip():
                r = json.loads(line)
                out.setdefault(tp, {})[r["variant"]] = r
    return out


def shaped_metrics(lines, ref):
    m = {}
    for tp, by in lines.items():
        for v, r in by.items():
            if v == ref:
                continue
            for k in ["w_mm", "u_mm", "UV_charge_ratio_rel", "UV_amp_ratio_rel"] + [f"{s}_{x}" for s in ("U", "V") for x in SHAPED_KEYS]:
                if k in r:
                    m.setdefault(v, {}).setdefault(tp, {})[k] = r[k]
    return m


def finite(x):
    return isinstance(x, (int, float)) and math.isfinite(x)


def criteria_dry_run(crit, lines, ref, cand):
    """Frozen v1 criteria (unchanged) applied to one event; returns per-check results and PASS / FAIL / MISSING."""
    C = crit["criteria"]
    tps = [f"{t:g}" for t in crit["electronics_provisional"]["Tp_ns"]]
    checks, missing = [], []
    for tp in tps:
        by = lines.get(tp, {})
        r, r0 = by.get(cand), by.get(ref)
        if r is None or r0 is None:
            missing.append(f"Tp{tp}"); continue
        def add(name, key, val, limit, status, note=""):
            checks.append({"criterion": name, "key": key, "Tp_ns": float(tp), "value": val, "limit": limit, "status": status, "note": note})
        for s in ("U", "V"):
            pk0 = sorted(r0[f"{s}_peaks"], reverse=True)
            k = f"{s}_lead_index_changed"; v = r.get(k)
            if v == 0:
                add("lead_strip", k, v, 0, "PASS")
            else:
                amb = len(pk0) > 1 and (pk0[0] - pk0[1]) <= 0.01 * pk0[0]
                add("lead_strip", k, v, 0, "AMBIGUOUS" if amb else "FAIL", "S0 top-2 within 1 % of the lead" if amb else "")
            thr = r0["thresholds_fCeq"][s]
            for key in C["hit_flips"]["keys"]:
                k = key.replace("{U,V}", s); f = k.split("_")[-1]; v = r.get(k)
                if v == 0:
                    add("hit_flips", k, v, 0, "PASS"); continue
                t = thr[f]; a0 = [p >= t for p in r0[f"{s}_peaks"]]; a1 = [p >= t for p in r[f"{s}_peaks"]]
                flips = [i for i, (x, y) in enumerate(zip(a0, a1)) if x != y]
                marg = bool(flips) and all(abs(r0[f"{s}_peaks"][i] - t) <= 0.02 * t for i in flips)
                add("hit_flips", k, v, 0, "MARGINAL" if marg else "FAIL", f"flipping strips at S0 peak/thr {[round(r0[f'{s}_peaks'][i] / t, 4) for i in flips]}")
            for name, lim_key, scale in (("threshold_time_f0.1", "limit_abs_ns", 1.0), ("window_charge", "limit_abs_rel", 1.0),
                                         ("second_over_lead", "limit_abs", 1.0), ("shaped_L2_ceiling", "limit", 1.0)):
                k = C[name]["key"].replace("{U,V}", s); v = r.get(k); lim = C[name][lim_key] * scale
                add(name, k, v, lim, "PASS" if finite(v) and abs(v) <= lim else "FAIL", "" if finite(v) else "nonfinite")
        for name, lim_key, scale in (("UV_charge_ratio", "limit_abs_rel", 1.0), ("w_centroid", "limit_abs_um", 1e-3),
                                     ("u_centroid", "limit_abs_um", 1e-3)):
            k = C[name]["key"]; v = r.get(k); lim = C[name][lim_key] * scale
            add(name, k, v, lim, "PASS" if finite(v) and abs(v) <= lim else "FAIL", "" if finite(v) else "nonfinite")
    if missing or not checks:
        st = "MISSING"
    else:
        st = "FAIL" if any(c["status"] == "FAIL" for c in checks) else "PASS"
    fails = [f"{c['key']}@Tp{c['Tp_ns']:g}={c['value']}" for c in checks if c["status"] == "FAIL"]
    return {"status": st, "missing": missing, "n_checks": len(checks), "failed": fails,
            "exemptions": [f"{c['key']}@Tp{c['Tp_ns']:g} {c['status']}" for c in checks if c["status"] in ("AMBIGUOUS", "MARGINAL")],
            "checks": checks}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("dir", type=Path); ap.add_argument("--muon", required=True); ap.add_argument("--label", required=True)
    ap.add_argument("--independent", required=True); ap.add_argument("--t0", type=int, choices=(0, 1), required=True)
    ap.add_argument("--criteria", type=Path, required=True)
    a = ap.parse_args()
    V = a.dir; M = a.muon
    E = Path((V / "far_dir.txt").read_text().strip()) if (V / "far_dir.txt").is_file() else None
    vs = steps_of(V / "steps.tsv", True)
    es = steps_of(E / "steps.tsv", False) if E else {}
    fails, incon, notes = [], [], []

    # integrity of the archive (both sides read it)
    for d in [V] + ([E] if E else []):
        b, f = d / "archive_sha256_before.txt", d / "archive_sha256_after.txt"
        if b.is_file() and f.is_file():
            if b.read_text().split()[0] != f.read_text().split()[0]:
                fails.append(f"archive sha256 changed ({d.name})")
        elif d == V and "disk" not in vs:
            incon.append(f"archive sha256 before/after missing in {d}")
    if "disk" in vs:
        incon.append(vs["disk"]["detail"])
    if vs.get("far_experiment", {}).get("status") == "FAIL" and not es:
        incon.append(vs["far_experiment"]["detail"])

    # S0_T1 reference
    s0_cause = None
    for k in ("export", "T1", "moments"):
        if es.get(k, {}).get("status") == "FAIL":
            s0_cause = f"{k} failed: {es[k]['detail']}"
    groups = jload(E / "groups" / "groups.json") if E else None
    if s0_cause is None and E:
        if es.get("E2_groups", {}).get("status") != "OK":
            s0_cause = f"grouping: {es.get('E2_groups', {}).get('status', 'missing')}"
        else:
            for g in sorted((groups or {}).get("groups", {})):
                st = es.get(f"E2_replay_{g}", {})
                if st.get("status") != "OK":
                    s0_cause = (f"frozen replay of {g} {'REFUSED (outside the validated S0 domain)' if st.get('status') == 'INCONCLUSIVE' else st.get('status', 'missing')}: "
                                f"{st.get('detail', '')[:200]}")
                    break
            if s0_cause is None and es.get("E2_S0_local_census", {}).get("status") != "OK":
                s0_cause = f"S0 local sum / census: {es.get('E2_S0_local_census', {}).get('status', 'missing')}"
    have_t1 = bool(E) and s0_cause is None and (E / "S0_local" / "replay" / "waveforms.csv").is_file()
    if E and s0_cause is None and not have_t1:
        s0_cause = "S0_local outputs missing"

    # S0_T0 reference
    t0 = {"enabled": bool(a.t0), "status": "DISABLED", "cause": None}
    if a.t0:
        bad = [(k, s) for k, s in vs.items() if k.startswith("T0_") and s["status"] != "OK"]
        if bad:
            t0.update(status="FAILED", cause="; ".join(f"{k} {s['status']}: {s['detail'][:150]}" for k, s in bad))
        elif vs.get("T0_S0_local_census", {}).get("status") == "OK":
            t0["status"] = "OK"
        else:
            t0.update(status="FAILED" if have_t1 else "NOT_RUN", cause="T0 local sum missing")
    have_t0 = t0["status"] == "OK"

    # production K2A, domain, consistency
    kp = vs.get("K2A_prod", {}); kres = jload(V / "K2A_prod" / "result.json")
    domain = {"status": None, "status_names": None, "blocked_charge_frac": None}
    if kres:
        q = kres.get("cloud_charge_fC") or 0.0
        domain.update(status_names=kres.get("status_names"), usable=kres.get("usable"),
                      blocked_charge_frac=(kres.get("blocked_charge_fC") or 0.0) / q if q else None,
                      shifts_differing_from_production_rule=kres.get("shifts_differing_from_computed"),
                      groups=kres.get("groups"), clouds=kres.get("clouds"))
    if vs.get("births_moments", {}).get("status") == "FAIL":
        fails.append("p3_birth_moments failed: " + vs["births_moments"]["detail"][:200])
    if kp.get("status") == "OK":
        domain["status"] = "IN_DOMAIN"
    elif kp.get("status") == "BLOCKED":
        domain["status"] = "OUT_OF_DOMAIN"
    elif kp.get("status") == "FAIL":
        fails.append("production K2A: " + kp["detail"][:200])
    elif E and not s0_cause:
        incon.append("production K2A not run")
    if vs.get("moments_archive_vs_export", {}).get("status") == "FAIL":
        fails.append("archive moments differ from the export-path moments (s_compare_moments 1e-9)")
    py = es.get("E2_K2A_local", {}).get("status"); cc = vs.get("K2A_cc_export", {}).get("status")
    if py and cc:
        if py == "OK" and cc == "OK":
            st = vs.get("identity_cc_vs_python", {}).get("status")
            if st == "FAIL":
                fails.append("compiled K2A != Python K2A_local (1e-12)")
            elif st != "PASS":
                incon.append(f"implementation identity {st or 'missing'}")
        elif py != "OK" and cc == "BLOCKED":
            notes.append(f"Python K2A_local refuses as the compiled evaluator blocks (consistent): {es['E2_K2A_local']['detail'][:150]}")
        elif cc == "FAIL":
            fails.append("compiled K2A on the export moments failed")
        else:
            fails.append(f"inconsistent domain: Python K2A_local {py}, compiled {cc}")
    if domain["status"] == "IN_DOMAIN" and cc == "BLOCKED" or domain["status"] == "OUT_OF_DOMAIN" and cc == "OK":
        fails.append(f"moment path changes the K2A domain (archive moments {kp.get('status')}, export moments {cc})")
    if domain["status"] == "IN_DOMAIN" and cc == "OK":
        st = vs.get("moment_path_archive_vs_export", {}).get("status")
        if st == "FAIL":
            fails.append("K2A on archive moments differs from K2A on export moments (1e-7)")
        elif st != "PASS":
            incon.append(f"moment-path response check {st or 'missing'}")

    # comparisons
    cmps = {n: cmp_metrics(jload(V / f"cmp_{n}.json")) for n in CMP_NAMES}
    for n in CMP_NAMES:
        st = vs.get(f"cmp_{n}", {}).get("status")
        if st == "FAIL":
            fails.append(f"comparison {n}: numerical breakage (FAIL)")
        elif st == "INCONCLUSIVE":
            incon.append(f"comparison {n} INCONCLUSIVE")
    need_cmp = ["S0T1_vs_K2A"] + (["S0T0_vs_K2A", "S0T0_vs_S0T1"] if have_t0 else [])
    shaped_runs = {"refS0T1_ionelec": "S0T1", "refS0T0_ionelec": "S0T0", "refS0T1_ion": "S0T1"}
    shp = {}
    for n, ref in shaped_runs.items():
        st = vs.get(f"shaped_{n}", {}).get("status")
        if st == "FAIL":
            fails.append(f"shaped comparison {n} failed (strip labels / clock not identical?): {vs[f'shaped_{n}']['detail'][:150]}")
        lines = shaped_lines(V / f"shaped_{n}")
        if lines:
            shp[n] = {"reference": ref, "lines": lines, "metrics": shaped_metrics(lines, ref)}

    # outcome
    if fails:
        outcome = "FAIL"
    elif s0_cause:
        outcome = "S0_FAILED"
    elif domain["status"] == "OUT_OF_DOMAIN":
        outcome = "OUT_OF_DOMAIN"
    elif incon or domain["status"] != "IN_DOMAIN":
        outcome = "INCONCLUSIVE"
    else:
        missing = [n for n in need_cmp if not cmps.get(n)] + [f"shaped_{n}" for n in (["refS0T1_ionelec", "refS0T1_ion"] + (["refS0T0_ionelec"] if have_t0 else [])) if n not in shp]
        if missing:
            outcome = "INCONCLUSIVE"; incon.append("missing: " + ", ".join(missing))
        else:
            outcome = "VALID" if have_t0 else "VALID_T1_ONLY"

    # criteria dry run (development data, DIAGNOSTIC)
    crit = jload(a.criteria)
    dry = {}
    if crit and outcome in ("VALID", "VALID_T1_ONLY"):
        for n, ref in (("refS0T1_ionelec", "S0T1"), ("refS0T0_ionelec", "S0T0")):
            if n in shp:
                dry[n] = criteria_dry_run(crit, shp[n]["lines"], ref, "K2A")
                dry[n]["label"] = "DIAGNOSTIC dry run of the frozen sealed criteria v1 on a NON-SEALED event (not a sealed verdict)"
                dry[n]["decisive_reference"] = n == "refS0T1_ionelec"

    # runtime split (measured; seconds)
    et, vt = timing_of(E / "timing.tsv") if E else {}, timing_of(V / "timing.tsv")
    w = lambda T, k: (T.get(k) or {}).get("wall_s")
    s = lambda T, pre: sum(x["wall_s"] for k, x in T.items() if k.startswith(pre))
    bj = jload(V / "prod" / "birth_moments.csv.json") or {}
    kt = (kres or {}).get("timing", {}); km = (kres or {}).get("memory", {})
    c1, l1, n1, st1 = chosen_wall(E / "S0groups") if E else (0, 0, 0, {})
    c0, l0, n0, st0 = chosen_wall(V / "T0" / "S0groups")
    rt = {"unit": "s (wall, this machine, single process each)",
          "S0": {"rkf_export_archive_read_plus_ion_drift": w(et, "export_rkf"), "T1_collapse": w(et, "T1_collapse"),
                 "group_export_T1": w(et, "group_export"), "weighting_replay_T1_chosen_attempt_sum": c1,
                 "weighting_replay_T1_ladder_total_tm": s(et, "replay_g_"), "weighting_replay_T1_attempts_sum": l1,
                 "electron_step_T1": s(et, "electron_step_"), "local_sum_T1": w(et, "local_sum"), "groups_T1": n1,
                 "late_step_chosen_T1": st1, "base_centred_S0_extra": w(et, "base_centred_S0"),
                 "export_path_moments_s_kernel_not_production": w(et, "moments_export"),
                 "k2a_local_python_not_production": w(et, "k2a_local_python"),
                 "rkf_export_peak_rss_kib": (et.get("export_rkf") or {}).get("peak_rss_kib"),
                 "group_export_T0": w(vt, "T0_group_export"), "weighting_replay_T0_chosen_attempt_sum": c0 if n0 else None,
                 "weighting_replay_T0_attempts_sum": l0 if n0 else None, "electron_step_T0": s(vt, "T0_electron_step_") if n0 else None,
                 "local_sum_T0": w(vt, "T0_local_sum"), "late_step_chosen_T0": st0},
          "K2A": {"archive_read_validate_min": bj.get("read_validate_s_min"), "moments_min": bj.get("moments_s_min"),
                  "births_process_wall_all_repeats": w(vt, "births_moments"), "births_repeat": bj.get("repeat"),
                  "births_peak_rss_kib": bj.get("peak_rss_kib"), "archive_bytes": bj.get("archive_bytes"), "rows": bj.get("rows"),
                  "kernel_verify": kt.get("kernel_verify_s"), "kernel_load": kt.get("kernel_load_s"),
                  "evaluate_cold": kt.get("evaluate_cold_s"), "evaluate_warm_median": kt.get("evaluate_warm_median_s"),
                  "evaluate_repeat": kt.get("repeat"), "repeats_bit_identical": kt.get("repeats_bit_identical"),
                  "eval_process_wall_all_repeats": w(vt, "k2a_eval"), "eval_peak_rss_kib": km.get("peak_rss_kib_end"),
                  "kernel_K_bytes": km.get("kernel_K_bytes"), "kernel_K_mapped": km.get("kernel_K_mapped")}}
    S, K = rt["S0"], rt["K2A"]
    def tot(*xs):
        return sum(xs) if all(finite(x) for x in xs) else None
    S["total_T1_single_pass"] = tot(S["rkf_export_archive_read_plus_ion_drift"], S["T1_collapse"], S["group_export_T1"],
                                    S["weighting_replay_T1_chosen_attempt_sum"], S["electron_step_T1"], S["local_sum_T1"])
    S["total_T1_with_ladder"] = tot(S["rkf_export_archive_read_plus_ion_drift"], S["T1_collapse"], S["group_export_T1"],
                                    S["weighting_replay_T1_ladder_total_tm"], S["electron_step_T1"], S["local_sum_T1"])
    S["total_T0_single_pass"] = tot(S["rkf_export_archive_read_plus_ion_drift"], S["group_export_T0"],
                                    S["weighting_replay_T0_chosen_attempt_sum"], S["electron_step_T0"], S["local_sum_T0"]) if n0 else None
    K["per_event_warm"] = tot(K["moments_min"], K["evaluate_warm_median"])
    K["per_event_warm_incl_archive_read"] = tot(K["archive_read_validate_min"], K["moments_min"], K["evaluate_warm_median"])
    K["single_event_cold_incl_kernel_open"] = tot(K["archive_read_validate_min"], K["moments_min"], K["kernel_verify"],
                                                  K["kernel_load"], K["evaluate_cold"])
    sp = {}
    for ref in ("total_T0_single_pass", "total_T1_single_pass"):
        for den in ("per_event_warm", "per_event_warm_incl_archive_read", "single_event_cold_incl_kernel_open"):
            if finite(S.get(ref)) and finite(K.get(den)) and K[den] > 0:
                sp[f"{ref}/{den}"] = S[ref] / K[den]
    rt["speedup_measured"] = sp
    rt["notes"] = ["RKF export = frozen exporter (archive read + Garfield RKF ion drift in one process; not separable without "
                   "modifying it); the archive read alone is measured by p3_birth_moments (archive_read_validate_min)",
                   "single_pass = chosen late-step attempt only (ladder retries excluded); with_ladder = every attempt",
                   "K2A per_event_warm = moments + warm kernel evaluation with the kernel already open (production: kernel "
                   "opened once per job; Fast4 supplies moments in memory, so no archive read)",
                   "p3_birth_moments and k2a_strip_eval times are minima / medians over repeats; S0 steps ran once"]

    rec = {"tool": "k2a_validation_muon.py", "label": a.label, "muon": M, "independent": a.independent, "outcome": outcome,
           "failures": fails, "inconclusive": incon, "notes": notes, "s0_t1": {"available": have_t1, "cause": s0_cause},
           "s0_t0": t0, "k2a_domain": domain, "far_dir": str(E) if E else None,
           "base_centred_S0": es.get("base_centred_S0", {}).get("status"), "far_analysis": es.get("E2_far_analysis", {}).get("status"),
           "truncation": jload(E / "truncation.json") if E else None,
           "comparisons": cmps, "shaped": {n: {"reference": x["reference"], "metrics": x["metrics"]} for n, x in shp.items()},
           "criteria_v1_dry_run": dry, "runtime": rt,
           "steps": {"validation": vs, "far": es}}
    (V / "outcome.json").write_text(json.dumps(rec, indent=1, default=str) + "\n")
    m = cmps.get("S0T1_vs_K2A") or {}
    extra = f"; S0T1 vs K2A worst rel L2 {m.get('worst_L2_rel')}" if m else ""
    dr = dry.get("refS0T1_ionelec", {}).get("status")
    print(f"{a.label} muon {M}: {outcome}" + (f" (dry-run v1 {dr})" if dr else "") + extra +
          (f"; FAIL: {fails[0]}" if fails else "") + (f"; S0: {s0_cause}" if s0_cause else "") +
          (f"; K2A {domain['status_names']}" if domain["status"] == "OUT_OF_DOMAIN" else ""))


if __name__ == "__main__":
    main()
