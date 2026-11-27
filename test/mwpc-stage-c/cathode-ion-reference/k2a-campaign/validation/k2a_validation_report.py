#!/usr/bin/env python3
"""Aggregate a run_k2a_validation.sh directory: K2A vs detailed S0 on independent muons, with explicit outcomes.

  k2a_validation_report.py OUTDIR [--selection selection.csv]
  k2a_validation_report.py NEW_DIR --from OUTDIR1 OUTDIR2 ... [--selection selection.csv]   combined report of several
        runs (e.g. pilot + remainder); NEW_DIR must not exist; a slot in two runs, or runs with different labels, kernels
        or git commits, are refused; nothing in the source runs is modified

Reads OUTDIR/slots.txt, provenance.txt and OUTDIR/muons/muon_*/outcome.json (k2a_validation_muon.py). Every listed
slot appears in the tables; a slot without outcome.json is INCONCLUSIVE. Writes OUTDIR/validation_report.json,
validation_report.md and per_muon.csv.

Sections, keeping EXECUTION and ACCURACY apart:
  execution   per slot: VALID / VALID_T1_ONLY / OUT_OF_DOMAIN / S0_FAILED / FAIL / INCONCLUSIVE with causes. VALID means
              references and prediction were computed and are consistent, not that the prediction is accurate.
  intended OOD  slots selected as the out-of-domain stratum (selection reason out_of_domain) are listed in their own
              table and never enter the accuracy set, whatever their execution outcome.
  accuracy    over the accuracy set = VALID executions that are not intended-OOD: the frozen criteria v1 dry run of the
              PRIMARY K2A (own production grouping) vs S0_T1 (PASS / FAIL counts; DIAGNOSTIC, non-sealed), and the
              error distributions (n, signed mean, rms, median |x|, p90 |x|, max |x|) for
                model          S0_T1 vs K2A      primary prediction, kernel / recentering / spatial error at T1 timing
                total          S0_T0 vs K2A      primary prediction against the true-birth-time reference
                timing         S0_T0 vs S0_T1    the T1 time collapse alone (no K2A)
                model_imposed  S0_T1 vs K2Aimp   DIAGNOSTIC, S0 grouping imposed (kernel / interpolation error only)
                grouping       K2Aimp vs K2A     effect of the production grouping where it differs from S0's
              unshaped (current, cumulative charge, window charge, U/V ratio, strip charge, lead) and shaped
              (provisional CR-RC^3, Tp 100/200 ns: L2, lead, hit flips, threshold / peak time, second/lead, w, u).
  shift assignment  per slot: production rule vs S0 grouping, AGREE / DISAGREE, primaries and charge.
  features, coverage, runtime (measured; cathode-response evaluation separated from archive I/O and moment building;
              avalanche generation excluded on both sides, no acceleration of it claimed).
Exit 0 if every slot has an execution outcome and none is FAIL or INCONCLUSIVE, 1 if any FAIL, 3 otherwise. No
acceptance decision is made here.
"""
import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np

FEATURES = ["abs_theta_u", "abs_theta_w", "log_gain", "mean_abs_A", "n_wires_5pct", "max_abs_m", "clouds"]
COMPONENTS = [("model", "S0T1_vs_K2A", "refS0T1_ionelec", "K2A"), ("total", "S0T0_vs_K2A", "refS0T0_ionelec", "K2A"),
              ("timing", "S0T0_vs_S0T1", "refS0T0_ionelec", "S0T1"), ("model_imposed", "S0T1_vs_K2Aimp", "refS0T1_ionelec", "K2Aimp"),
              ("grouping", "K2Aimp_vs_K2A", "refK2Aimp_ionelec", "K2A")]
COMP_TEXT = {"model": "S0_T1 vs K2A (PRIMARY prediction, own grouping)", "total": "S0_T0 vs K2A (PRIMARY, true birth times)",
             "timing": "S0_T0 vs S0_T1 (T1 time collapse alone)", "model_imposed": "S0_T1 vs K2Aimp (DIAGNOSTIC: S0 grouping imposed)",
             "grouping": "K2Aimp vs K2A (effect of the production grouping)"}
UNSHAPED = ["unshaped_L2_rel", "cumulative_L2_rel", "window_charge_rel_diff", "max_strip_charge_diff_over_lead", "lead_changed",
            "shaped_L2_rel_Tp100", "shaped_L2_rel_Tp200", "ref_charge_on_ref_only_strips_frac", "cand_charge_on_cand_only_strips_frac"]
SHAPED = ["L2_rel", "lead_index_changed", "hitflips_f0.05", "hitflips_f0.1", "t_f0.1_ns", "lead_peak_time_ns", "lead_peak_rel",
          "charge_rel", "second_over_lead", "strip_peak_maxabs_over_T0lead"]
HEADLINE = [("model", "unshaped", "ion_plus_electron_step.U.unshaped_L2_rel"), ("model", "unshaped", "ion_plus_electron_step.V.unshaped_L2_rel"),
            ("model", "shaped", "Tp100.U_L2_rel"), ("model", "shaped", "Tp100.w_mm"), ("model", "shaped", "Tp100.u_mm"),
            ("model", "unshaped", "ion_plus_electron_step.UV_charge_ratio_rel_diff")]


def stats(x):
    x = np.asarray([v for v in x if isinstance(v, (int, float)) and math.isfinite(v)], float)
    if len(x) == 0:
        return {"n": 0}
    ax = np.abs(x)
    return {"n": int(len(x)), "mean": float(x.mean()), "rms": float(np.sqrt((x ** 2).mean())), "median_abs": float(np.median(ax)),
            "p90_abs": float(np.percentile(ax, 90)), "max_abs": float(ax.max())}


def spearman(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 4:
        return None
    rx = np.argsort(np.argsort(x[ok])); ry = np.argsort(np.argsort(y[ok]))
    if rx.std() == 0 or ry.std() == 0:
        return None
    return float(np.corrcoef(rx, ry)[0, 1])


def flat_metrics(o):
    """component -> {'unshaped': {key: val}, 'shaped': {Tp.key: val}} for one outcome."""
    out = {}
    for comp, cmpn, shpn, var in COMPONENTS:
        d = {"unshaped": {}, "shaped": {}}
        c = (o.get("comparisons") or {}).get(cmpn) or {}
        for term in ("ion_plus_electron_step", "ion_only"):
            for s in ("U", "V"):
                for k in UNSHAPED:
                    key = f"{term}.{s}.{k}"
                    if key in c:
                        d["unshaped"][key] = c[key]
            if f"{term}.UV_charge_ratio_rel_diff" in c:
                d["unshaped"][f"{term}.UV_charge_ratio_rel_diff"] = c[f"{term}.UV_charge_ratio_rel_diff"]
        sm = (((o.get("shaped") or {}).get(shpn) or {}).get("metrics") or {}).get(var) or {}
        for tp, m in sm.items():
            for k in ["w_mm", "u_mm", "UV_charge_ratio_rel"] + [f"{s}_{x}" for s in ("U", "V") for x in SHAPED]:
                if k in m:
                    d["shaped"][f"Tp{tp}.{k}"] = m[k]
        if comp in ("model", "model_imposed"):            # ion-only shaped diagnostics vs S0_T1
            si = (((o.get("shaped") or {}).get("refS0T1_ion") or {}).get("metrics") or {}).get(var) or {}
            for tp, m in si.items():
                for k in ("U_L2_rel", "V_L2_rel", "w_mm", "u_mm"):
                    if k in m:
                        d["shaped"][f"ion_only.Tp{tp}.{k}"] = m[k]
        out[comp] = d
    return out


def fmt(v, k=""):
    if v is None or (isinstance(v, float) and not math.isfinite(v)):
        return "n/a"
    if isinstance(v, int):
        return str(v)
    if k.endswith("_mm"):
        return f"{1e3 * v:.3g} um"
    return f"{v:.3g}"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out", type=Path); ap.add_argument("--selection", type=Path)
    ap.add_argument("--from", dest="src", type=Path, nargs="+", help="combine these run directories into NEW_DIR")
    a = ap.parse_args()
    O = a.out
    srcs = a.src or [O]
    if a.src:
        if O.exists():
            raise SystemExit(f"refusing: {O} exists")
    def pv(d):
        t = (d / "provenance.txt").read_text() if (d / "provenance.txt").is_file() else ""
        return {k: v for k, v in (l.split("=", 1) for l in t.splitlines() if "=" in l and not l.startswith("tool_sha256"))}
    provs = [pv(d) for d in srcs]
    for k in ("label", "kernel", "git_commit", "run", "tag", "half_strips", "late_steps", "criteria_v1_sha256"):
        vals = {p.get(k) for p in provs}
        if len(vals) > 1:
            raise SystemExit(f"refusing to combine runs with different {k}: {sorted(map(str, vals))}")
    label = provs[0].get("label", "UNKNOWN")
    slots, home = [], {}
    for d in srcs:
        for l in ((d / "slots.txt").read_text().splitlines() if (d / "slots.txt").is_file() else []):
            if l.strip():
                if l.strip() in home:
                    raise SystemExit(f"slot {l.strip()} appears in {home[l.strip()]} and {d}")
                home[l.strip()] = d; slots.append(l.strip())
    if a.src:
        O.mkdir(parents=True)
        (O / "combined_from.txt").write_text("".join(f"{d.resolve()}\n" for d in srcs))
    sel = {}
    if a.selection and a.selection.is_file():
        sel = {r["slot"]: r for r in csv.DictReader(open(a.selection)) if r.get("slot")}
    outs = {}
    for m in slots:
        p = home[m] / "muons" / f"muon_{m}" / "outcome.json"
        outs[m] = json.loads(p.read_text()) if p.is_file() else {"muon": m, "outcome": "INCONCLUSIVE", "inconclusive": ["no outcome.json"]}
    counts = {}
    for o in outs.values():
        counts[o["outcome"]] = counts.get(o["outcome"], 0) + 1
    intended = {m for m in slots if sel.get(m, {}).get("reason") == "out_of_domain"}
    valid = [m for m, o in outs.items() if o["outcome"] in ("VALID", "VALID_T1_ONLY")]
    accset = [m for m in valid if m not in intended]
    acc_counts = {}
    for m in accset:
        st = (outs[m].get("accuracy") or {}).get("status", "NOT_EVALUATED")
        acc_counts[st] = acc_counts.get(st, 0) + 1
    fm = {m: flat_metrics(outs[m]) for m in valid}

    acc = {}
    for comp, *_ in COMPONENTS:
        acc[comp] = {}
        for kind in ("unshaped", "shaped"):
            keys = sorted({k for m in accset for k in fm[m][comp][kind]})
            acc[comp][kind] = {k: stats([fm[m][comp][kind].get(k) for m in accset]) for k in keys}

    def feat(m, f):
        try:
            return float(sel[m][f])
        except (KeyError, ValueError, TypeError):
            return float("nan")
    dep = {}
    for comp, kind, key in HEADLINE:
        ys = [abs(fm[m][comp][kind].get(key, float("nan"))) for m in accset]
        dep[f"{comp}.{key}"] = {f: spearman([feat(m, f) for m in accset], ys) for f in FEATURES}
    cov = {"stations": {}, "channels": {}, "features_accuracy_set": {}, "features_selected": {}}
    for m in slots:
        r = sel.get(m, {})
        for k, d in (("station", "stations"), ("channel", "channels")):
            v = r.get(k, "?")
            cov[d].setdefault(v, {"selected": 0, "accuracy_set": 0}); cov[d][v]["selected"] += 1
            cov[d][v]["accuracy_set"] += m in accset
    for f in FEATURES:
        for name, ms in (("features_selected", slots), ("features_accuracy_set", accset)):
            xs = [feat(m, f) for m in ms]; xs = [x for x in xs if math.isfinite(x)]
            cov[name][f] = [min(xs), max(xs)] if xs else None

    rkeys_s = ["rkf_export_archive_read_plus_ion_drift", "T1_collapse", "group_export_T1", "weighting_replay_T1_chosen_attempt_sum",
               "weighting_replay_T1_ladder_total_tm", "electron_step_T1", "local_sum_T1", "total_T1_single_pass",
               "total_T1_with_ladder", "response_T1_excl_archive_read_approx", "group_export_T0",
               "weighting_replay_T0_chosen_attempt_sum", "electron_step_T0", "local_sum_T0", "total_T0_single_pass",
               "response_T0_excl_archive_read_approx", "base_centred_S0_extra"]
    rkeys_k = ["archive_io_read_validate_min", "moments_construction_min", "response_eval_warm_median", "response_eval_cold",
               "kernel_open_verify_once", "kernel_open_load_once", "moments_plus_response_eval", "archive_io_plus_moments_plus_eval",
               "single_event_cold_incl_kernel_open", "eval_peak_rss_kib", "births_peak_rss_kib"]
    have_rt = [m for m in slots if outs[m].get("runtime")]
    def rstat(side, k):
        xs = [outs[m]["runtime"][side].get(k) for m in have_rt]
        xs = np.asarray([x for x in xs if isinstance(x, (int, float)) and math.isfinite(x)], float)
        return {"n": int(len(xs)), "median": float(np.median(xs)), "min": float(xs.min()), "max": float(xs.max())} if len(xs) else {"n": 0}
    runtime = {"scope": "cathode-response calculation from an existing Garfield++ avalanche archive; avalanche generation "
                        "excluded on both sides; no acceleration of it is claimed",
               "S0": {k: rstat("S0", k) for k in rkeys_s}, "K2A": {k: rstat("K2A", k) for k in rkeys_k}, "speedup_measured": {}}
    for k in sorted({k for m in have_rt for k in outs[m]["runtime"].get("speedup_measured", {})}):
        xs = np.asarray([outs[m]["runtime"]["speedup_measured"][k] for m in have_rt if k in outs[m]["runtime"]["speedup_measured"]], float)
        runtime["speedup_measured"][k] = {"n": int(len(xs)), "median": float(np.median(xs)), "min": float(xs.min()), "max": float(xs.max())}

    dry = {}
    for n in ("primary_vs_S0T1", "primary_vs_S0T0", "imposed_vs_S0T1"):
        rows = {m: (outs[m].get("criteria_v1_dry_run") or {}).get(n) for m in accset}
        dry[n] = {"PASS": sorted(m for m, r in rows.items() if r and r["status"] == "PASS"),
                  "FAIL": {m: r["failed"] for m, r in rows.items() if r and r["status"] == "FAIL"},
                  "MISSING": sorted(m for m, r in rows.items() if not r or r["status"] == "MISSING"),
                  "exemptions": {m: r["exemptions"] for m, r in rows.items() if r and r.get("exemptions")}}

    cols = ["slot", "selection_reason", "intended_ood", "execution", "accuracy_v1_primary", "independent", "channel", "station"] + FEATURES + [
        "k2a_domain", "k2a_status_names", "blocked_charge_frac", "k2a_imposed", "shift_assignment", "n_primaries_disagree",
        "disagree_charge_frac", "s0_t1", "s0_t0", "base_centred_S0", "far_analysis",
        "model_U_unshaped_L2", "model_V_unshaped_L2", "model_U_shaped_L2_Tp100", "model_V_shaped_L2_Tp100", "model_U_charge_rel",
        "model_V_charge_rel", "model_UV_ratio_rel", "model_w_um_Tp100", "model_u_um_Tp100", "model_U_t_f0.1_ns_Tp100",
        "total_U_shaped_L2_Tp100", "total_V_shaped_L2_Tp100", "total_w_um_Tp100", "total_u_um_Tp100",
        "timing_U_shaped_L2_Tp100", "timing_w_um_Tp100", "imposed_U_shaped_L2_Tp100", "imposed_w_um_Tp100",
        "grouping_U_shaped_L2_Tp100", "grouping_w_um_Tp100", "grouping_U_charge_only_in_one_grouping", "accuracy_v1_primary_vs_S0T0", "accuracy_v1_imposed_vs_S0T1",
        "S0_T0_response_s", "K2A_response_eval_s", "K2A_moments_s", "K2A_archive_io_s",
        "speedup_cathode_response_eval_T0", "speedup_end_to_end_from_archive_T0", "cause"]
    lines_csv = []
    for m in slots:
        o = outs[m]; r = sel.get(m, {}); f = fm.get(m)
        def g(comp, kind, key, scale=1.0):
            v = f[comp][kind].get(key) if f else None
            return "" if v is None else f"{scale * v:.6g}"
        rt = o.get("runtime") or {}; S, K, SP = rt.get("S0", {}), rt.get("K2A", {}), rt.get("speedup_measured", {})
        dom = o.get("k2a_domain") or {}; sa = o.get("shift_assignment") or {}; dr = o.get("criteria_v1_dry_run") or {}
        cause = "; ".join((o.get("failures") or []) + ([o["s0_t1"]["cause"]] if (o.get("s0_t1") or {}).get("cause") else []) + (o.get("inconclusive") or []))
        lines_csv.append({"slot": m, "selection_reason": r.get("reason", ""), "intended_ood": "yes" if m in intended else "",
                          "execution": o["outcome"],
                          "accuracy_v1_primary": (o.get("accuracy") or {}).get("status", "") if m not in intended else "not counted (intended OOD)",
                          "independent": o.get("independent", ""), "channel": r.get("channel", ""),
                          "station": r.get("station", ""), **{k: r.get(k, "") for k in FEATURES},
                          "k2a_domain": dom.get("status") or "", "k2a_status_names": dom.get("status_names") or "",
                          "blocked_charge_frac": "" if dom.get("blocked_charge_frac") is None else f"{dom['blocked_charge_frac']:.4g}",
                          "k2a_imposed": (o.get("k2a_imposed") or {}).get("status") or "", "shift_assignment": sa.get("status", ""),
                          "n_primaries_disagree": "" if sa.get("n_disagree") is None else sa["n_disagree"],
                          "disagree_charge_frac": "" if sa.get("disagree_charge_frac") is None else f"{sa['disagree_charge_frac']:.4g}",
                          "s0_t1": "OK" if (o.get("s0_t1") or {}).get("available") else "NO", "s0_t0": (o.get("s0_t0") or {}).get("status", ""),
                          "base_centred_S0": o.get("base_centred_S0") or "", "far_analysis": o.get("far_analysis") or "",
                          "model_U_unshaped_L2": g("model", "unshaped", "ion_plus_electron_step.U.unshaped_L2_rel"),
                          "model_V_unshaped_L2": g("model", "unshaped", "ion_plus_electron_step.V.unshaped_L2_rel"),
                          "model_U_shaped_L2_Tp100": g("model", "shaped", "Tp100.U_L2_rel"), "model_V_shaped_L2_Tp100": g("model", "shaped", "Tp100.V_L2_rel"),
                          "model_U_charge_rel": g("model", "shaped", "Tp100.U_charge_rel"), "model_V_charge_rel": g("model", "shaped", "Tp100.V_charge_rel"),
                          "model_UV_ratio_rel": g("model", "shaped", "Tp100.UV_charge_ratio_rel"),
                          "model_w_um_Tp100": g("model", "shaped", "Tp100.w_mm", 1e3), "model_u_um_Tp100": g("model", "shaped", "Tp100.u_mm", 1e3),
                          "model_U_t_f0.1_ns_Tp100": g("model", "shaped", "Tp100.U_t_f0.1_ns"),
                          "total_U_shaped_L2_Tp100": g("total", "shaped", "Tp100.U_L2_rel"), "total_V_shaped_L2_Tp100": g("total", "shaped", "Tp100.V_L2_rel"),
                          "total_w_um_Tp100": g("total", "shaped", "Tp100.w_mm", 1e3), "total_u_um_Tp100": g("total", "shaped", "Tp100.u_mm", 1e3),
                          "timing_U_shaped_L2_Tp100": g("timing", "shaped", "Tp100.U_L2_rel"), "timing_w_um_Tp100": g("timing", "shaped", "Tp100.w_mm", 1e3),
                          "imposed_U_shaped_L2_Tp100": g("model_imposed", "shaped", "Tp100.U_L2_rel"), "imposed_w_um_Tp100": g("model_imposed", "shaped", "Tp100.w_mm", 1e3),
                          "grouping_U_shaped_L2_Tp100": g("grouping", "shaped", "Tp100.U_L2_rel"), "grouping_w_um_Tp100": g("grouping", "shaped", "Tp100.w_mm", 1e3),
                          "grouping_U_charge_only_in_one_grouping": g("grouping", "unshaped", "ion_plus_electron_step.U.ref_charge_on_ref_only_strips_frac"),
                          "accuracy_v1_primary_vs_S0T0": (dr.get("primary_vs_S0T0") or {}).get("status", ""),
                          "accuracy_v1_imposed_vs_S0T1": (dr.get("imposed_vs_S0T1") or {}).get("status", ""),
                          "S0_T0_response_s": fmt(S.get("total_T0_single_pass")), "K2A_response_eval_s": fmt(K.get("response_eval_warm_median")),
                          "K2A_moments_s": fmt(K.get("moments_construction_min")), "K2A_archive_io_s": fmt(K.get("archive_io_read_validate_min")),
                          "speedup_cathode_response_eval_T0": fmt(SP.get("S0_T0/cathode_response_eval")),
                          "speedup_end_to_end_from_archive_T0": fmt(SP.get("S0_T0/end_to_end_from_archive_incl_io")), "cause": cause})
    with open(O / "per_muon.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols); w.writeheader(); w.writerows(lines_csv)

    bad = counts.get("FAIL", 0); inc = counts.get("INCONCLUSIVE", 0)
    status = "FAIL" if bad else ("INCONCLUSIVE" if inc or not slots else "COMPLETE")
    refs = next((o.get("references") for o in outs.values() if o.get("references")), None)
    rep = {"tool": "k2a_validation_report.py", "label": label, "dir": str(O), "sources": [str(d) for d in srcs],
           "provenance": provs[0], "status": status, "status_meaning": "COMPLETE = every slot executed without FAIL / INCONCLUSIVE; "
           "says nothing about accuracy", "slots": slots, "references": refs,
           "execution_counts": counts, "intended_ood": sorted(intended), "accuracy_set": accset,
           "accuracy_v1_primary_counts": acc_counts,
           "outcomes": {m: {k: outs[m].get(k) for k in ("outcome", "accuracy", "independent", "failures", "inconclusive", "notes",
                                                         "s0_t1", "s0_t0", "k2a_domain", "k2a_imposed", "shift_assignment")} for m in slots},
           "accuracy": acc, "feature_dependence_spearman": dep, "coverage": cov, "runtime": runtime, "criteria_v1_dry_run": dry,
           "meaning": "errors of the PRIMARY K2A (own grouping) against the detailed S0 are REPORTED with their distribution over "
                      "the accuracy set; the S0-imposed K2Aimp is diagnostic only; the v1 dry run is DIAGNOSTIC (non-sealed); "
                      "no acceptance decision is made here"}
    (O / "validation_report.json").write_text(json.dumps(rep, indent=1, default=str) + "\n")

    L = [f"# K2A vs detailed S0 validation ({label})", "", f"`{O}`", "",
         "Labels: " + ("SYNTHETIC test data (toy transport and kernel), NOT NA60+ physics." if label != "REAL" else
                       "REAL archived 2100 V Garfield++ K256 ion births (non-sealed), frozen S0 workflow, validated kernel."), "",
         "## Summary: execution is not accuracy", "",
         f"- Execution ({len(slots)} slots): " + ", ".join(f"{k} {v}" for k, v in sorted(counts.items())) +
         f". Run status **{status}** (all slots executed without FAIL/INCONCLUSIVE; says nothing about accuracy).",
         f"- Intended out-of-domain events (selection stratum, never counted as predictions): " +
         (", ".join(f"{m} -> {outs[m]['outcome']}" for m in sorted(intended)) if intended else "none"),
         f"- Accuracy set (VALID, not intended-OOD): {len(accset)} events. Frozen criteria v1 dry run, PRIMARY K2A vs S0_T1 "
         f"(DIAGNOSTIC, non-sealed): " + (", ".join(f"{k} {v}" for k, v in sorted(acc_counts.items())) or "none evaluated"), "",
         "## References and predictions", ""]
    for k, v in (refs or {}).items():
        L.append(f"- **{k}**: {v}")
    L += ["", "## Execution per slot", "",
          "| slot | selection | execution | accuracy v1 (primary vs S0_T1) | K2A domain | K2Aimp | shift assignment | S0_T1 | S0_T0 | cause / notes |",
          "|---|---|---|---|---|---|---|---|---|---|"]
    for r in lines_csv:
        sa = r["shift_assignment"] + (f" ({r['n_primaries_disagree']} prim., {r['disagree_charge_frac']} of charge)" if r["shift_assignment"] == "DISAGREE" else "")
        L.append(f"| {r['slot']} | {r['selection_reason']} | {r['execution']} | {r['accuracy_v1_primary']} | {r['k2a_domain']} "
                 f"{r['k2a_status_names'] if r['k2a_domain'] == 'OUT_OF_DOMAIN' else ''} | {r['k2a_imposed']} | {sa} | {r['s0_t1']} | {r['s0_t0']} | {r['cause'][:200]} |")
    dis = [(m, outs[m]["shift_assignment"]) for m in slots if (outs[m].get("shift_assignment") or {}).get("status") == "DISAGREE"]
    L += ["", "## Shift assignment: production rule vs S0 grouping", ""]
    if not dis:
        L.append("Every evaluated slot: AGREE (production-rule and S0-imposed K2A responses bit-identical) or not evaluated; see table.")
    for m, sa in dis:
        L.append(f"- {m}: {sa.get('n_disagree')} primaries, {fmt(sa.get('disagree_charge_frac'))} of the cloud charge; "
                 + "; ".join(f"p{d['primary']} S0 {d['s0_m']} vs K2A {d['k2a_m']} (Dbar {d['Dbar_U_mm']:.4g}/{d['Dbar_V_mm']:.4g} mm, "
                             f"margin {min(d['rule_margin_U'], d['rule_margin_V']):.2g})" for d in (sa.get("disagreements") or [])[:5])
                 + (f"; only S0 {sa['only_s0']}" if sa.get("only_s0") else "") + (f"; only K2A {sa['only_k2a']}" if sa.get("only_k2a") else "")
                 + "; effect on waveforms and observables: component 'grouping' below (common-strip differences and the charge on "
                   "strips only one grouping computes) and per_muon.csv grouping_*")
    if intended:
        L += ["", "## Intended out-of-domain events (reported separately, not predictions)", "",
              "| slot | execution | K2A domain | blocking statuses | blocked charge frac | S0_T1 | S0_T0 |", "|---|---|---|---|---|---|---|"]
        for r in lines_csv:
            if r["intended_ood"]:
                L.append(f"| {r['slot']} | {r['execution']} | {r['k2a_domain']} | {r['k2a_status_names']} | {r['blocked_charge_frac']} | {r['s0_t1']} | {r['s0_t0']} |")
        L.append("")
        L.append("An intended-OOD slot that executes as VALID (in domain after recentering) is also kept out of the accuracy set.")
    L += ["", f"## Error distributions over the accuracy set (n = {len(accset)})", "",
          "Signal ion + electron step unless prefixed ion_only. Shaped: provisional CR-RC^3 stand-in (not VMM3a). "
          "The electron term of every variant is the approximate instantaneous step.", ""]
    for comp, *_ in COMPONENTS:
        for kind in ("unshaped", "shaped"):
            tab = acc[comp][kind]
            if not tab:
                continue
            L += [f"### {COMP_TEXT[comp]} - {kind}", "", "| metric | n | mean | rms | median abs | p90 abs | max abs |", "|---|---|---|---|---|---|---|"]
            for k, st in tab.items():
                if st.get("n"):
                    L.append(f"| {k} | {st['n']} | {fmt(st['mean'], k)} | {fmt(st['rms'], k)} | {fmt(st['median_abs'], k)} | {fmt(st['p90_abs'], k)} | {fmt(st['max_abs'], k)} |")
            L.append("")
    L += ["## Feature dependence (Spearman rank correlation of |metric|, accuracy set; descriptive, small n)", "",
          "| metric | " + " | ".join(FEATURES) + " |", "|---" * (len(FEATURES) + 1) + "|"]
    for k, d in dep.items():
        L.append(f"| {k} | " + " | ".join("n/a" if d[f] is None else f"{d[f]:+.2f}" for f in FEATURES) + " |")
    L += ["", "## Coverage", "", "Stations (selected / accuracy set): " + ", ".join(f"{k}: {v['selected']}/{v['accuracy_set']}" for k, v in sorted(cov["stations"].items())),
          "", "Channels (selected / accuracy set): " + ", ".join(f"{k}: {v['selected']}/{v['accuracy_set']}" for k, v in sorted(cov["channels"].items())), "",
          "| feature | selected range | accuracy-set range |", "|---|---|---|"]
    for f in FEATURES:
        L.append(f"| {f} | {cov['features_selected'][f]} | {cov['features_accuracy_set'][f]} |")
    L += ["", "## Runtime (measured, seconds per event; median [min, max])", "",
          "**Scope:** both sides start from the same archived Garfield++ avalanche. Avalanche generation is excluded on both "
          "sides; nothing here measures or claims an acceleration of Garfield++ avalanche generation.", "",
          "| S0 step | n | median | min | max |", "|---|---|---|---|---|"]
    for k, st in runtime["S0"].items():
        if st.get("n"):
            L.append(f"| {k} | {st['n']} | {st['median']:.4g} | {st['min']:.4g} | {st['max']:.4g} |")
    L += ["", "| K2A step | n | median | min | max |", "|---|---|---|---|---|"]
    for k, st in runtime["K2A"].items():
        if st.get("n"):
            L.append(f"| {k} | {st['n']} | {st['median']:.4g} | {st['min']:.4g} | {st['max']:.4g} |")
    L += ["", "| measured speedup (S0 / K2A) | n | median | min | max |", "|---|---|---|---|---|"]
    for k, st in runtime["speedup_measured"].items():
        L.append(f"| {k} | {st['n']} | {st['median']:.4g} | {st['min']:.4g} | {st['max']:.4g} |")
    L += ["", "- **cathode_response_eval**: S0 response (RKF ion drift + weighting replay + electron step + grouping/sum, archive "
          "read subtracted approximately with the p3_birth_moments read time: the frozen exporter reads and drifts in one process) "
          "over the K2A warm kernel evaluation only (kernel already open, moments given: the Fast4 -> K2A production case).",
          "- **response_from_births_incl_moments**: the same numerator over K2A moment construction + evaluation.",
          "- **end_to_end_from_archive_incl_io**: both sides including their archive I/O.",
          "- **single_event_cold_incl_kernel_open**: adds the one-time kernel verification and load.",
          "- Single pass = chosen late-step attempt; ladder retries listed separately.", "",
          "## DIAGNOSTIC dry run of the frozen sealed criteria v1 (accuracy set; non-sealed; not a sealed verdict)", ""]
    for n, d in dry.items():
        L.append(f"- {n}: PASS {len(d['PASS'])}, FAIL {len(d['FAIL'])}, missing {len(d['MISSING'])}"
                 + (" (decisive definition)" if n == "primary_vs_S0T1" else " (information only)"))
        for m, fl in d["FAIL"].items():
            L.append(f"  - {m}: " + "; ".join(fl[:6]) + (" ..." if len(fl) > 6 else ""))
        for m, ex in d["exemptions"].items():
            L.append(f"  - {m} exemptions: " + "; ".join(ex[:6]))
    (O / "validation_report.md").write_text("\n".join(L) + "\n")
    print(f"{label} K2A validation report: execution {status} (" + ", ".join(f"{k} {v}" for k, v in sorted(counts.items()))
          + f"); accuracy set {len(accset)}: v1 dry run " + (", ".join(f"{k} {v}" for k, v in sorted(acc_counts.items())) or "none")
          + f"; intended OOD {sorted(intended)}; {O / 'validation_report.md'}")
    return 1 if bad else (3 if status == "INCONCLUSIVE" else 0)

if __name__ == "__main__":
    raise SystemExit(main())
