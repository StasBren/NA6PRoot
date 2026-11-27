#!/usr/bin/env python3
"""Explicit outcome, metrics, runtime split and criteria dry run of ONE muon of run_k2a_validation.sh.

  k2a_validation_muon.py MUON_DIR --muon M --label REAL|SYNTHETIC --independent TEXT --t0 0|1 --criteria v1.json

MUON_DIR = OUT/muons/muon_M (steps.tsv, timing.tsv, far_dir.txt -> the far-mode experiment muon directory E).
Writes MUON_DIR/outcome.json; prints one line. Never raises on missing pieces: they make the outcome INCONCLUSIVE.

Two separate statements per muon:
EXECUTION outcome (first that applies) - says that references and prediction were computed and are consistent, NOT
that the prediction is accurate:
  FAIL           numerical / integrity failure: archive changed; moment path differs (archive moments vs export
                 moments, s_compare_moments 1e-9; responses 1e-7); compiled K2A differs from the Python K2A_local
                 (1e-12) or only one of them refuses; agreeing shift assignments but non-identical responses; a
                 comparison reports numerical breakage (FAIL); variants not placeable on one strip label set; evaluator
                 error
  S0_FAILED      the established S0_T1 reference is not available (exporter, T1 collapse, grouping, a frozen replay
                 refusal or failure, census): cause named. K2A domain status is still reported
  OUT_OF_DOMAIN  the PRIMARY K2A (production grouping rule) has blocking statuses (named, with the blocked charge
                 fraction); the S0 references are still computed; no accuracy is claimed
  INCONCLUSIVE   anything required missing or inconclusive
  VALID          S0_T1, S0_T0 (if enabled), the primary K2A and the diagnostic K2Aimp all computed, every consistency
                 check passed. VALID_T1_ONLY: as VALID but S0_T0 disabled or failed (named)
ACCURACY (only for VALID executions): the frozen sealed criteria v1, unchanged, applied to the PRIMARY K2A against S0_T1
(ion + electron, Tp 100 and 200 ns, both planes): PASS / FAIL, labelled DIAGNOSTIC (non-sealed event). Shown alongside,
never decisive: primary K2A vs S0_T0, and the S0-imposed K2Aimp vs S0_T1.
The primary K2A never receives S0 information; the S0-imposed variant is a diagnostic, and the per-primary assignment
comparison (AGREE / DISAGREE with charge) states whether the two coincide.
"""
import argparse
import csv
import json
import math
import re
from pathlib import Path

CMP_NAMES = ("S0T1_vs_K2A", "S0T0_vs_K2A", "S0T1_vs_K2Aimp", "S0T0_vs_K2Aimp", "S0T0_vs_S0T1", "K2Aimp_vs_K2A")
SHAPED_RUNS = {"refS0T1_ionelec": "S0T1", "refS0T0_ionelec": "S0T0", "refS0T1_ion": "S0T1", "refK2Aimp_ionelec": "K2Aimp"}
REFERENCE_LABELS = {
    "S0_T1": "detailed reference, established E2b workflow: frozen RKF ion paths of the archived births, per-primary T1 "
             "birth time (ion-weighted mean), frozen weighting replay; electron term = APPROXIMATE instantaneous step",
    "S0_T0": "detailed reference with true per-ion birth times (uncollapsed frozen export), same groups, frozen weighting "
             "replay; electron term = the same APPROXIMATE instantaneous step (not a microscopic electron-current simulation)",
    "K2A": "PRIMARY prediction: compact moments from the archive, compiled K2A with its own production grouping rule",
    "K2Aimp": "DIAGNOSTIC: the same moments with the S0 grouping imposed (isolates kernel / interpolation error)"}
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

    # PRIMARY K2A: production grouping rule (no S0 information); its domain decides OUT_OF_DOMAIN
    kp = vs.get("K2A_prod", {}); kres = jload(V / "K2A_prod" / "result.json")
    ki = vs.get("K2A_imposed", {}); kires = jload(V / "K2A_imposed" / "result.json")
    domain = {"status": None, "status_names": None, "blocked_charge_frac": None, "grouping": "production rule (no S0 input)"}
    if kres:
        q = kres.get("cloud_charge_fC") or 0.0
        domain.update(status_names=kres.get("status_names"), usable=kres.get("usable"),
                      blocked_charge_frac=(kres.get("blocked_charge_fC") or 0.0) / q if q else None,
                      groups=kres.get("groups"), clouds=kres.get("clouds"))
    if vs.get("births_moments", {}).get("status") == "FAIL":
        fails.append("p3_birth_moments failed: " + vs["births_moments"]["detail"][:200])
    if kp.get("status") == "OK":
        domain["status"] = "IN_DOMAIN"
    elif kp.get("status") == "BLOCKED":
        domain["status"] = "OUT_OF_DOMAIN"
    elif kp.get("status") == "FAIL":
        fails.append("production K2A: " + kp["detail"][:200])
    elif E and not s0_cause and "births_moments" in vs:
        incon.append("production K2A not run")
    if vs.get("k2a_stage", {}).get("status") == "FAIL":
        fails.append(vs["k2a_stage"]["detail"])
    if vs.get("moments_archive_vs_export", {}).get("status") == "FAIL":
        fails.append("archive moments differ from the export-path moments (s_compare_moments 1e-9)")
    elif vs.get("moments_archive_vs_export", {}).get("status") == "INCONCLUSIVE":
        incon.append("archive-vs-export moments not compared")

    # DIAGNOSTIC K2Aimp (S0 grouping imposed) and the equivalence tests (Python and export moments use the S0 grouping)
    imposed = {"status": ki.get("status"), "status_names": (kires or {}).get("status_names"),
               "role": "DIAGNOSTIC: S0 grouping imposed; isolates kernel / interpolation error; not evidence for the prediction"}
    if ki.get("status") == "FAIL":
        fails.append("K2A with the S0 grouping imposed: " + ki["detail"][:200])
    elif ki.get("status") == "INCONCLUSIVE" and not s0_cause:
        incon.append(ki["detail"][:200])
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
    if ki.get("status") in ("OK", "BLOCKED") and cc in ("OK", "BLOCKED") and ki["status"] != cc:
        fails.append(f"moment path changes the K2A domain (archive moments {ki['status']}, export moments {cc})")
    if ki.get("status") == "OK" and cc == "OK":
        st = vs.get("moment_path_archive_vs_export", {}).get("status")
        if st == "FAIL":
            fails.append("K2A on archive moments differs from K2A on export moments (1e-7)")
        elif st != "PASS":
            incon.append(f"moment-path response check {st or 'missing'}")
    if kp.get("status") in ("OK", "BLOCKED") and ki.get("status") in ("OK", "BLOCKED") and kp["status"] != ki["status"]:
        notes.append(f"the grouping changes the K2A domain: production rule {kp['status']}, S0-imposed {ki['status']}")

    # per-primary shift assignment: production rule vs S0 grouping
    sa = jload(V / "shift_assignment.json") or {}
    sst = vs.get("shift_assignment", {}).get("status")
    assignment = {"status": sst or "NOT_RUN", "n_disagree": sa.get("n_disagree"), "disagree_charge_frac": sa.get("disagree_charge_frac"),
                  "only_s0": sa.get("only_s0"), "only_k2a": sa.get("only_k2a"), "near_tie_primaries": sa.get("near_tie_primaries"),
                  "disagreements": sa.get("disagreements"), "grouping_identity": vs.get("grouping_identity", {}).get("status")}
    if sst == "FAIL":
        fails.append("shift assignment comparison failed: " + vs["shift_assignment"]["detail"][:200])
    elif sst == "INCONCLUSIVE":
        incon.append("shift assignment not compared: " + vs["shift_assignment"]["detail"][:200])
    elif sst == "AGREE" and kp.get("status") == "OK" and ki.get("status") == "OK":
        if assignment["grouping_identity"] == "FAIL":
            fails.append("assignments agree but the production-rule and S0-imposed responses differ")
        elif assignment["grouping_identity"] != "PASS":
            incon.append("grouping identity check missing")
    elif sst is None and kp.get("status") == "OK" and E and (E / "groups" / "groups.json").is_file():
        incon.append("shift assignment not compared")

    # comparisons
    cmps = {n: cmp_metrics(jload(V / f"cmp_{n}.json")) for n in CMP_NAMES}
    for n in CMP_NAMES:
        st = vs.get(f"cmp_{n}", {}).get("status")
        if st == "FAIL":
            fails.append(f"comparison {n}: numerical breakage (FAIL)")
        elif st == "INCONCLUSIVE":
            incon.append(f"comparison {n} INCONCLUSIVE")
    if vs.get("label_union", {}).get("status") == "FAIL":
        fails.append("variants could not be put on one absolute-strip label set: " + vs["label_union"]["detail"][:150])
    lu = jload(V / "shaped_layout" / f"muon_{M}" / "label_union.json") or {}
    have_imp = ki.get("status") == "OK"
    need_cmp = ["S0T1_vs_K2A"] + (["S0T0_vs_K2A", "S0T0_vs_S0T1"] if have_t0 else []) + \
               (["S0T1_vs_K2Aimp", "K2Aimp_vs_K2A"] + (["S0T0_vs_K2Aimp"] if have_t0 else []) if have_imp else [])
    shp = {}
    for n, ref in SHAPED_RUNS.items():
        st = vs.get(f"shaped_{n}", {}).get("status")
        if st == "FAIL":
            fails.append(f"shaped comparison {n} failed: {vs[f'shaped_{n}']['detail'][:150]}")
        lines = shaped_lines(V / f"shaped_{n}")
        if lines:
            shp[n] = {"reference": ref, "lines": lines, "metrics": shaped_metrics(lines, ref)}
    need_shp = ["refS0T1_ionelec", "refS0T1_ion"] + (["refS0T0_ionelec"] if have_t0 else []) + (["refK2Aimp_ionelec"] if have_imp else [])

    # EXECUTION outcome (not an accuracy statement)
    if fails:
        outcome = "FAIL"
    elif s0_cause:
        outcome = "S0_FAILED"
    elif domain["status"] == "OUT_OF_DOMAIN":
        outcome = "OUT_OF_DOMAIN"
    elif incon or domain["status"] != "IN_DOMAIN":
        outcome = "INCONCLUSIVE"
    else:
        missing = [n for n in need_cmp if not cmps.get(n)] + [f"shaped_{n}" for n in need_shp if n not in shp]
        if missing:
            outcome = "INCONCLUSIVE"; incon.append("missing: " + ", ".join(missing))
        else:
            outcome = "VALID" if have_t0 else "VALID_T1_ONLY"

    # ACCURACY: frozen criteria v1 applied to the PRIMARY K2A vs S0_T1 (decisive in v1); diagnostics alongside
    crit = jload(a.criteria)
    dry = {}
    if crit and outcome in ("VALID", "VALID_T1_ONLY"):
        for n, ref, cand, decisive in (("primary_vs_S0T1", "S0T1", "K2A", True), ("primary_vs_S0T0", "S0T0", "K2A", False),
                                       ("imposed_vs_S0T1", "S0T1", "K2Aimp", False)):
            src = {"S0T1": "refS0T1_ionelec", "S0T0": "refS0T0_ionelec"}[ref]
            if src in shp and cand in next(iter(shp[src]["lines"].values()), {}):
                dry[n] = criteria_dry_run(crit, shp[src]["lines"], ref, cand)
                dry[n].update(label="DIAGNOSTIC dry run of the frozen sealed criteria v1 on a NON-SEALED event (not a sealed verdict)",
                              decisive=decisive, reference=ref, candidate=cand)
    accuracy = {"status": dry["primary_vs_S0T1"]["status"] if "primary_vs_S0T1" in dry else "NOT_EVALUATED",
                "definition": "frozen sealed criteria v1 (unchanged) on the PRIMARY K2A (production grouping) vs S0_T1, "
                              "ion + electron step, Tp 100 and 200 ns, both planes; evaluated only for VALID executions",
                "label": "DIAGNOSTIC (non-sealed event)"}

    # runtime split (measured; seconds). Avalanche generation (Garfield++) is NOT part of either side.
    et, vt = timing_of(E / "timing.tsv") if E else {}, timing_of(V / "timing.tsv")
    w = lambda T, k: (T.get(k) or {}).get("wall_s")
    s = lambda T, pre: sum(x["wall_s"] for k, x in T.items() if k.startswith(pre))
    bj = jload(V / "prod" / "birth_moments.csv.json") or {}
    kt = (kres or {}).get("timing", {}); km = (kres or {}).get("memory", {})
    c1, l1, n1, st1 = chosen_wall(E / "S0groups") if E else (0, 0, 0, {})
    c0, l0, n0, st0 = chosen_wall(V / "T0" / "S0groups")
    rt = {"unit": "s (wall, this machine, single process each)",
          "scope": "cathode-response calculation from an EXISTING Garfield++ avalanche archive; avalanche generation is "
                   "excluded on both sides and is not accelerated by K2A",
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
          "K2A": {"archive_io_read_validate_min": bj.get("read_validate_s_min"), "moments_construction_min": bj.get("moments_s_min"),
                  "births_process_wall_all_repeats": w(vt, "births_moments"), "births_repeat": bj.get("repeat"),
                  "births_peak_rss_kib": bj.get("peak_rss_kib"), "archive_bytes": bj.get("archive_bytes"), "rows": bj.get("rows"),
                  "kernel_open_verify_once": kt.get("kernel_verify_s"), "kernel_open_load_once": kt.get("kernel_load_s"),
                  "response_eval_cold": kt.get("evaluate_cold_s"), "response_eval_warm_median": kt.get("evaluate_warm_median_s"),
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
    rd = K["archive_io_read_validate_min"]
    for k in ("T0", "T1"):
        t = S[f"total_{k}_single_pass"]
        S[f"response_{k}_excl_archive_read_approx"] = t - rd if finite(t) and finite(rd) else None
    K["response_eval_only"] = K["response_eval_warm_median"]
    K["moments_plus_response_eval"] = tot(K["moments_construction_min"], K["response_eval_warm_median"])
    K["archive_io_plus_moments_plus_eval"] = tot(rd, K["moments_construction_min"], K["response_eval_warm_median"])
    K["single_event_cold_incl_kernel_open"] = tot(rd, K["moments_construction_min"], K["kernel_open_verify_once"],
                                                  K["kernel_open_load_once"], K["response_eval_cold"])
    sp = {}
    for k in ("T0", "T1"):
        pairs = {"cathode_response_eval": (f"response_{k}_excl_archive_read_approx", "response_eval_only"),
                 "response_from_births_incl_moments": (f"response_{k}_excl_archive_read_approx", "moments_plus_response_eval"),
                 "end_to_end_from_archive_incl_io": (f"total_{k}_single_pass", "archive_io_plus_moments_plus_eval"),
                 "single_event_cold_incl_kernel_open": (f"total_{k}_single_pass", "single_event_cold_incl_kernel_open")}
        for name, (num, den) in pairs.items():
            if finite(S.get(num)) and finite(K.get(den)) and K[den] > 0:
                sp[f"S0_{k}/{name}"] = S[num] / K[den]
    rt["speedup_measured"] = sp
    rt["notes"] = ["scope: both sides start from the same archived Garfield++ avalanche; avalanche generation is excluded "
                   "and no acceleration of it is claimed",
                   "cathode_response_eval = S0 response (RKF ion drift + weighting replay + electron step + grouping/sum, "
                   "archive read subtracted APPROXIMATELY using the p3_birth_moments read time, since the frozen exporter "
                   "reads the archive and drifts the ions in one process) / K2A warm kernel evaluation only",
                   "response_from_births_incl_moments adds the K2A moment construction; end_to_end_from_archive_incl_io "
                   "compares both including their archive I/O; single_event_cold adds the one-time kernel open",
                   "single_pass = chosen late-step attempt only (ladder retries excluded); with_ladder = every attempt",
                   "p3_birth_moments and k2a_strip_eval times are minima / medians over repeats; S0 steps ran once"]

    rec = {"tool": "k2a_validation_muon.py", "label": a.label, "muon": M, "independent": a.independent,
           "outcome": outcome, "outcome_meaning": "EXECUTION outcome (references and prediction computed and consistent); "
                                                  "accuracy is reported separately in 'accuracy'",
           "accuracy": accuracy, "failures": fails, "inconclusive": incon, "notes": notes,
           "references": REFERENCE_LABELS,
           "s0_t1": {"available": have_t1, "cause": s0_cause}, "s0_t0": t0, "k2a_domain": domain, "k2a_imposed": imposed,
           "shift_assignment": assignment, "label_union": {k: lu.get(k) for k in ("identical", "variants")} if lu else None,
           "far_dir": str(E) if E else None,
           "base_centred_S0": es.get("base_centred_S0", {}).get("status"), "far_analysis": es.get("E2_far_analysis", {}).get("status"),
           "truncation": jload(E / "truncation.json") if E else None,
           "comparisons": cmps, "shaped": {n: {"reference": x["reference"], "metrics": x["metrics"]} for n, x in shp.items()},
           "criteria_v1_dry_run": dry, "runtime": rt,
           "steps": {"validation": vs, "far": es}}
    (V / "outcome.json").write_text(json.dumps(rec, indent=1, default=str) + "\n")
    m = cmps.get("S0T1_vs_K2A") or {}
    extra = f"; S0T1 vs K2A worst rel L2 {m.get('worst_L2_rel')}" if m else ""
    print(f"{a.label} muon {M}: execution {outcome}; accuracy (v1 dry run, primary K2A vs S0_T1) {accuracy['status']}; "
          f"shift assignment {assignment['status']}" + extra +
          (f"; FAIL: {fails[0]}" if fails else "") + (f"; S0: {s0_cause}" if s0_cause else "") +
          (f"; K2A {domain['status_names']}" if domain["status"] == "OUT_OF_DOMAIN" else ""))

if __name__ == "__main__":
    main()
