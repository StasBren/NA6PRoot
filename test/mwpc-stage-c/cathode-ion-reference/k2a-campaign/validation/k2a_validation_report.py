#!/usr/bin/env python3
"""Aggregate a run_k2a_validation.sh directory: K2A vs detailed S0 on independent muons, with explicit outcomes.

  k2a_validation_report.py OUTDIR [--selection selection.csv]

Reads OUTDIR/slots.txt, provenance.txt and OUTDIR/muons/muon_*/outcome.json (k2a_validation_muon.py). Every listed
slot appears in the tables; a slot without outcome.json is INCONCLUSIVE. Writes OUTDIR/validation_report.json,
validation_report.md and per_muon.csv.

Sections: outcomes (VALID / VALID_T1_ONLY / OUT_OF_DOMAIN / S0_FAILED / FAIL / INCONCLUSIVE, with causes); domain
coverage (stations, channels, feature ranges of the evaluated events vs the selection); accuracy distributions over
the VALID events (n, signed mean, rms, median |x|, p90 |x|, max |x|) for three error components
  model   S0_T1 vs K2A   (spatial / kernel / recentering error at fixed T1 timing)
  total   S0_T0 vs K2A   (K2A against the true-birth-time reference)
  timing  S0_T0 vs S0_T1 (the T1 time-collapse part alone, no K2A involved)
unshaped (current, cumulative charge, window charge, U/V ratio, strip charge, lead) and shaped (provisional CR-RC^3,
Tp 100/200 ns: shaped L2, lead, hit flips, threshold time, peak time, second/lead, w and u centroids);
dependence on the event features (Spearman rank correlation, descriptive only for n < 20);
runtime split and MEASURED speedup per event; the DIAGNOSTIC dry run of the frozen sealed criteria v1 (non-sealed
events; reference S0_T1 decisive as in v1, S0_T0 shown for information). Exit 0 if every slot has an outcome and none
is FAIL or INCONCLUSIVE, 1 if any FAIL, 3 otherwise. Accuracy is reported, not accepted or rejected here.
"""
import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np

FEATURES = ["abs_theta_u", "abs_theta_w", "log_gain", "mean_abs_A", "n_wires_5pct", "max_abs_m", "clouds"]
COMPONENTS = [("model", "S0T1_vs_K2A", "refS0T1_ionelec", "K2A"), ("total", "S0T0_vs_K2A", "refS0T0_ionelec", "K2A"),
              ("timing", "S0T0_vs_S0T1", "refS0T0_ionelec", "S0T1")]
UNSHAPED = ["unshaped_L2_rel", "cumulative_L2_rel", "window_charge_rel_diff", "max_strip_charge_diff_over_lead", "lead_changed",
            "shaped_L2_rel_Tp100", "shaped_L2_rel_Tp200"]
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
        if comp == "model":                               # ion-only shaped diagnostics vs S0_T1
            si = (((o.get("shaped") or {}).get("refS0T1_ion") or {}).get("metrics") or {}).get("K2A") or {}
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
    a = ap.parse_args()
    O = a.out
    prov = (O / "provenance.txt").read_text() if (O / "provenance.txt").is_file() else ""
    label = next((l.split("=", 1)[1] for l in prov.splitlines() if l.startswith("label=")), "UNKNOWN")
    slots = [l.strip() for l in (O / "slots.txt").read_text().splitlines() if l.strip()] if (O / "slots.txt").is_file() else []
    sel = {}
    if a.selection and a.selection.is_file():
        sel = {r["slot"]: r for r in csv.DictReader(open(a.selection)) if r.get("slot")}
    outs = {}
    for m in slots:
        p = O / "muons" / f"muon_{m}" / "outcome.json"
        outs[m] = json.loads(p.read_text()) if p.is_file() else {"muon": m, "outcome": "INCONCLUSIVE", "inconclusive": ["no outcome.json"]}
    counts = {}
    for o in outs.values():
        counts[o["outcome"]] = counts.get(o["outcome"], 0) + 1
    valid = [m for m, o in outs.items() if o["outcome"] in ("VALID", "VALID_T1_ONLY")]
    fm = {m: flat_metrics(outs[m]) for m in valid}

    # accuracy distributions over VALID events
    acc = {}
    for comp, *_ in COMPONENTS:
        acc[comp] = {}
        for kind in ("unshaped", "shaped"):
            keys = sorted({k for m in valid for k in fm[m][comp][kind]})
            acc[comp][kind] = {k: stats([fm[m][comp][kind].get(k) for m in valid]) for k in keys}

    # features
    def feat(m, f):
        try:
            return float(sel[m][f])
        except (KeyError, ValueError, TypeError):
            return float("nan")
    dep = {}
    for comp, kind, key in HEADLINE:
        ys = [abs(fm[m][comp][kind].get(key, float("nan"))) for m in valid]
        dep[f"{comp}.{key}"] = {f: spearman([feat(m, f) for m in valid], ys) for f in FEATURES}
    cov = {"stations": {}, "channels": {}, "features_evaluated": {}, "features_selected": {}}
    for m in slots:
        r = sel.get(m, {})
        for k, d in (("station", "stations"), ("channel", "channels")):
            v = r.get(k, "?")
            cov[d].setdefault(v, {"selected": 0, "valid": 0}); cov[d][v]["selected"] += 1
            cov[d][v]["valid"] += m in valid
    for f in FEATURES:
        for name, ms in (("features_selected", slots), ("features_evaluated", valid)):
            xs = [feat(m, f) for m in ms]; xs = [x for x in xs if math.isfinite(x)]
            cov[name][f] = [min(xs), max(xs)] if xs else None

    # runtime
    rkeys_s = ["rkf_export_archive_read_plus_ion_drift", "T1_collapse", "group_export_T1", "weighting_replay_T1_chosen_attempt_sum",
               "weighting_replay_T1_ladder_total_tm", "electron_step_T1", "local_sum_T1", "total_T1_single_pass",
               "total_T1_with_ladder", "group_export_T0", "weighting_replay_T0_chosen_attempt_sum", "electron_step_T0",
               "local_sum_T0", "total_T0_single_pass", "base_centred_S0_extra"]
    rkeys_k = ["archive_read_validate_min", "moments_min", "kernel_verify", "kernel_load", "evaluate_cold", "evaluate_warm_median",
               "per_event_warm", "per_event_warm_incl_archive_read", "single_event_cold_incl_kernel_open", "eval_peak_rss_kib",
               "births_peak_rss_kib"]
    have_rt = [m for m in slots if outs[m].get("runtime")]
    def rstat(side, k):
        xs = [outs[m]["runtime"][side].get(k) for m in have_rt]
        xs = np.asarray([x for x in xs if isinstance(x, (int, float)) and math.isfinite(x)], float)
        return {"n": int(len(xs)), "median": float(np.median(xs)), "min": float(xs.min()), "max": float(xs.max()), "sum": float(xs.sum())} if len(xs) else {"n": 0}
    runtime = {"S0": {k: rstat("S0", k) for k in rkeys_s}, "K2A": {k: rstat("K2A", k) for k in rkeys_k}}
    spk = sorted({k for m in have_rt for k in outs[m]["runtime"].get("speedup_measured", {})})
    runtime["speedup_measured"] = {}
    for k in spk:
        xs = np.asarray([outs[m]["runtime"]["speedup_measured"][k] for m in have_rt if k in outs[m]["runtime"]["speedup_measured"]], float)
        runtime["speedup_measured"][k] = {"n": int(len(xs)), "median": float(np.median(xs)), "min": float(xs.min()), "max": float(xs.max())}

    # criteria dry run
    dry = {}
    for n in ("refS0T1_ionelec", "refS0T0_ionelec"):
        rows = {m: outs[m].get("criteria_v1_dry_run", {}).get(n) for m in valid}
        dry[n] = {"PASS": sorted(m for m, r in rows.items() if r and r["status"] == "PASS"),
                  "FAIL": {m: r["failed"] for m, r in rows.items() if r and r["status"] == "FAIL"},
                  "MISSING": sorted(m for m, r in rows.items() if not r or r["status"] == "MISSING"),
                  "exemptions": {m: r["exemptions"] for m, r in rows.items() if r and r.get("exemptions")}}

    # per-muon csv
    cols = ["slot", "outcome", "independent", "channel", "station"] + FEATURES + [
        "k2a_domain", "k2a_status_names", "blocked_charge_frac", "s0_t1", "s0_t0", "base_centred_S0", "far_analysis",
        "model_U_unshaped_L2", "model_V_unshaped_L2", "model_U_shaped_L2_Tp100", "model_V_shaped_L2_Tp100", "model_U_charge_rel",
        "model_V_charge_rel", "model_UV_ratio_rel", "model_w_um_Tp100", "model_u_um_Tp100", "model_U_t_f0.1_ns_Tp100",
        "total_U_shaped_L2_Tp100", "total_V_shaped_L2_Tp100", "total_w_um_Tp100", "total_u_um_Tp100",
        "timing_U_shaped_L2_Tp100", "timing_w_um_Tp100", "dryrun_v1_S0T1", "dryrun_v1_S0T0",
        "S0_T0_s", "S0_T1_s", "K2A_warm_s", "K2A_warm_incl_archive_s", "speedup_T0_warm", "speedup_T0_incl_archive", "cause"]
    lines_csv = []
    for m in slots:
        o = outs[m]; r = sel.get(m, {}); f = fm.get(m)
        def g(comp, kind, key, scale=1.0):
            v = f[comp][kind].get(key) if f else None
            return "" if v is None else f"{scale * v:.6g}"
        rt = o.get("runtime") or {}; S, K, SP = rt.get("S0", {}), rt.get("K2A", {}), rt.get("speedup_measured", {})
        dom = o.get("k2a_domain") or {}
        cause = "; ".join((o.get("failures") or []) + ([o["s0_t1"]["cause"]] if (o.get("s0_t1") or {}).get("cause") else []) + (o.get("inconclusive") or []))
        lines_csv.append({"slot": m, "outcome": o["outcome"], "independent": o.get("independent", ""), "channel": r.get("channel", ""),
                          "station": r.get("station", ""), **{k: r.get(k, "") for k in FEATURES},
                          "k2a_domain": dom.get("status") or "", "k2a_status_names": dom.get("status_names") or "",
                          "blocked_charge_frac": "" if dom.get("blocked_charge_frac") is None else f"{dom['blocked_charge_frac']:.4g}",
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
                          "dryrun_v1_S0T1": ((o.get("criteria_v1_dry_run") or {}).get("refS0T1_ionelec") or {}).get("status", ""),
                          "dryrun_v1_S0T0": ((o.get("criteria_v1_dry_run") or {}).get("refS0T0_ionelec") or {}).get("status", ""),
                          "S0_T0_s": fmt(S.get("total_T0_single_pass")), "S0_T1_s": fmt(S.get("total_T1_single_pass")),
                          "K2A_warm_s": fmt(K.get("per_event_warm")), "K2A_warm_incl_archive_s": fmt(K.get("per_event_warm_incl_archive_read")),
                          "speedup_T0_warm": fmt(SP.get("total_T0_single_pass/per_event_warm")),
                          "speedup_T0_incl_archive": fmt(SP.get("total_T0_single_pass/per_event_warm_incl_archive_read")), "cause": cause})
    with open(O / "per_muon.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols); w.writeheader(); w.writerows(lines_csv)

    bad = counts.get("FAIL", 0); inc = counts.get("INCONCLUSIVE", 0)
    status = "FAIL" if bad else ("INCONCLUSIVE" if inc or not slots else "COMPLETE")
    rep = {"tool": "k2a_validation_report.py", "label": label, "dir": str(O), "status": status, "slots": slots,
           "outcome_counts": counts, "valid": valid, "outcomes": {m: {k: outs[m].get(k) for k in ("outcome", "independent", "failures",
           "inconclusive", "notes", "s0_t1", "s0_t0", "k2a_domain")} for m in slots},
           "accuracy": acc, "feature_dependence_spearman": dep, "coverage": cov, "runtime": runtime, "criteria_v1_dry_run": dry,
           "meaning": "accuracy of K2A against the detailed S0 on these events is REPORTED with its distribution; no acceptance "
                      "decision is made here; the dry run of the frozen v1 criteria is DIAGNOSTIC (non-sealed events)"}
    (O / "validation_report.json").write_text(json.dumps(rep, indent=1, default=str) + "\n")

    # markdown
    L = [f"# K2A vs detailed S0 validation ({label})", "",
         f"`{O}`  -  status **{status}**  -  {len(slots)} slots: " + ", ".join(f"{k} {v}" for k, v in sorted(counts.items())), "",
         "Labels: " + ("SYNTHETIC test data (toy transport and kernel), NOT NA60+ physics." if label != "REAL" else
                       "REAL archived 2100 V Garfield++ K256 ion births (non-sealed), frozen S0 workflow, validated kernel."), "",
         "## Outcomes", "", "| slot | outcome | independent | K2A domain | S0_T1 | S0_T0 | cause / notes |", "|---|---|---|---|---|---|---|"]
    for r in lines_csv:
        L.append(f"| {r['slot']} | {r['outcome']} | {r['independent']} | {r['k2a_domain']} {r['k2a_status_names'] if r['k2a_domain'] == 'OUT_OF_DOMAIN' else ''} "
                 f"| {r['s0_t1']} | {r['s0_t0']} | {r['cause'][:200]} |")
    L += ["", "## Accuracy over the VALID events", "",
          "Components: **model** = S0_T1 vs K2A; **total** = S0_T0 vs K2A; **timing** = S0_T0 vs S0_T1 (T1 collapse only). "
          "Signal ion + electron step unless prefixed ion_only. Shaped: provisional CR-RC^3 stand-in (not VMM3a).", ""]
    for comp, *_ in COMPONENTS:
        for kind in ("unshaped", "shaped"):
            tab = acc[comp][kind]
            if not tab:
                continue
            L += [f"### {comp} - {kind}", "", "| metric | n | mean | rms | median abs | p90 abs | max abs |", "|---|---|---|---|---|---|---|"]
            for k, s in tab.items():
                if s.get("n"):
                    L.append(f"| {k} | {s['n']} | {fmt(s['mean'], k)} | {fmt(s['rms'], k)} | {fmt(s['median_abs'], k)} | {fmt(s['p90_abs'], k)} | {fmt(s['max_abs'], k)} |")
            L.append("")
    L += ["## Feature dependence (Spearman rank correlation of |metric|; descriptive, small n)", "",
          "| metric | " + " | ".join(FEATURES) + " |", "|---" * (len(FEATURES) + 1) + "|"]
    for k, d in dep.items():
        L.append(f"| {k} | " + " | ".join("n/a" if d[f] is None else f"{d[f]:+.2f}" for f in FEATURES) + " |")
    L += ["", "## Coverage", "", "Stations (selected / VALID): " + ", ".join(f"{k}: {v['selected']}/{v['valid']}" for k, v in sorted(cov["stations"].items())),
          "", "Channels (selected / VALID): " + ", ".join(f"{k}: {v['selected']}/{v['valid']}" for k, v in sorted(cov["channels"].items())), "",
          "| feature | selected range | VALID range |", "|---|---|---|"]
    for f in FEATURES:
        L.append(f"| {f} | {cov['features_selected'][f]} | {cov['features_evaluated'][f]} |")
    L += ["", "## Runtime (measured, seconds per event; median [min, max] over events)", "", "| S0 step | n | median | min | max |", "|---|---|---|---|---|"]
    for k, s in runtime["S0"].items():
        if s.get("n"):
            L.append(f"| {k} | {s['n']} | {s['median']:.4g} | {s['min']:.4g} | {s['max']:.4g} |")
    L += ["", "| K2A step | n | median | min | max |", "|---|---|---|---|---|"]
    for k, s in runtime["K2A"].items():
        if s.get("n"):
            L.append(f"| {k} | {s['n']} | {s['median']:.4g} | {s['min']:.4g} | {s['max']:.4g} |")
    L += ["", "| measured speedup (S0 / K2A) | n | median | min | max |", "|---|---|---|---|---|"]
    for k, s in runtime["speedup_measured"].items():
        L.append(f"| {k} | {s['n']} | {s['median']:.4g} | {s['min']:.4g} | {s['max']:.4g} |")
    L += ["", "RKF export = archive read + Garfield RKF ion drift (frozen exporter, one process). Single pass = chosen late-step "
          "attempt; ladder retries listed separately. K2A per-event warm = moments + kernel evaluation with the kernel open "
          "(Fast4 production gives moments in memory: no archive read).", "",
          "## DIAGNOSTIC dry run of the frozen sealed criteria v1 (non-sealed events; not a sealed verdict)", ""]
    for n, d in dry.items():
        L.append(f"- reference {n}: PASS {len(d['PASS'])}, FAIL {len(d['FAIL'])}, missing {len(d['MISSING'])}"
                 + ("" if n == "refS0T1_ionelec" else " (information only; v1 is defined against S0_T1)"))
        for m, fl in d["FAIL"].items():
            L.append(f"  - {m}: " + "; ".join(fl[:6]) + (" ..." if len(fl) > 6 else ""))
        for m, ex in d["exemptions"].items():
            L.append(f"  - {m} exemptions: " + "; ".join(ex[:6]))
    (O / "validation_report.md").write_text("\n".join(L) + "\n")
    print(f"{label} K2A validation report: {status}; " + ", ".join(f"{k} {v}" for k, v in sorted(counts.items())) + f"; {O / 'validation_report.md'}")
    return 1 if bad else (3 if status == "INCONCLUSIVE" else 0)


if __name__ == "__main__":
    raise SystemExit(main())
