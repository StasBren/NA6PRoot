#!/usr/bin/env python3
"""Aggregate a recenter_experiment_* directory into results.json / results.md with explicit acceptance.

  k2a_recenter_aggregate.py EXPERIMENT_DIR

Reads env_check.json, thresholds.json and muon_*/steps.tsv (+ the cmp_*.json each step names).
E1 acceptance (per muon): both relabellings PASS. Any FAIL -> E1 FAIL; otherwise any INCONCLUSIVE -> E1 INCONCLUSIVE.
E2 (per muon): FAIL if any step FAILs (numerical breakage, census, identity violated, missing data), INCONCLUSIVE if a
step is inconclusive, WARN if any documented warning threshold is exceeded, else OK. Reference integrity: the
sha256 of every reference file read (T1_export, S0, cloud_moments.csv, K2Akern_TAG) before == after; else FAIL.
Key metrics of each comparison are tabulated (rel L2 unshaped/shaped, window charge, U/V ratio, lead strip,
missing-strip charge), with the non-OK checks listed verbatim.
Exit: 0 when E1 PASS and E2 OK/WARN for every muon and references intact; 1 otherwise.
"""
import json
import sys
from pathlib import Path

ORDER = {"PASS": 0, "OK": 0, "WARN": 1, "INCONCLUSIVE": 2, "FAIL": 3}


def worst(sts, default="INCONCLUSIVE"):
    return max(sts, key=lambda s: ORDER.get(s, 3)) if sts else default


def fmt(x, spec=".3g"):
    return "n/a" if x is None else (f"{x:{spec}}" if isinstance(x, (int, float)) else str(x))


def main():
    X = Path(sys.argv[1])
    env = json.loads((X / "env_check.json").read_text()) if (X / "env_check.json").exists() else {"status": "missing"}
    th = json.loads((X / "thresholds.json").read_text())
    res = {"experiment": str(X), "environment": env, "thresholds": th, "muons": {}}
    L = [f"# Recentering experiment results: {X.name}", "",
         f"environment: {env.get('status')} (LUT {str(env.get('lut_sha256', ''))[:16]}, replay {str(env.get('replay_sha256', ''))[:16]}, "
         f"kernel {str(env.get('kernel_id', ''))[:16]}, half-strips {env.get('half_strips')}); "
         f"source {(X / 'SOURCE_GIT_SHA.txt').read_text().strip() if (X / 'SOURCE_GIT_SHA.txt').exists() else '?'}", ""]
    overall_ok = env.get("status") == "PASS"
    for E in sorted(X.glob("muon_*")):
        steps = []
        if (E / "steps.tsv").exists():
            for line in (E / "steps.tsv").read_text().splitlines():
                parts = (line.split("\t") + ["", "", "", ""])[:4]
                steps.append({"step": parts[0], "status": parts[1], "detail": parts[2], "json": parts[3]})
        e1 = [s for s in steps if s["step"].startswith("E1_")]
        e2 = [s for s in steps if not s["step"].startswith("E1_")]
        e1_status = ("PASS" if len(e1) == 2 and all(s["status"] == "PASS" for s in e1)
                     else "FAIL" if any(s["status"] == "FAIL" for s in e1) else "INCONCLUSIVE") if e1 else None
        e2_core = [s for s in e2 if s["step"] != "base_centred_S0"]
        e2_status = worst([s["status"] for s in e2_core]) if e2_core else "INCONCLUSIVE"
        b, a = E / "reference_sha256_before.txt", E / "reference_sha256_after.txt"
        integrity = (None if not b.exists() else
                     "PASS" if a.exists() and a.read_text() == b.read_text() and b.read_text().strip() else "FAIL")
        m = {"E1": e1_status, "E2": e2_status, "reference_integrity": integrity, "steps": steps, "comparisons": {}}
        L += [f"## {E.name}", "", f"E1 {e1_status or '(not run)'} | E2 {e2_status} | reference files unchanged: "
              f"{integrity or '(far event: new files only)'}", "",
              "| step | status | detail |", "|---|---|---|"]
        L += [f"| {s['step']} | {s['status']} | {s['detail'][:180]} |" for s in steps]
        rows = []
        for s in steps:
            if s["json"].endswith(".json") and Path(s["json"]).name.startswith("cmp_") and Path(s["json"]).exists():
                c = json.loads(Path(s["json"]).read_text())
                m["comparisons"][s["step"]] = {k: c.get(k) for k in ("status", "kind", "worst_L2_rel", "error")}
                m["comparisons"][s["step"]]["non_ok_checks"] = [x for x in c.get("checks", []) if x["status"] != "OK"]
                for term, T in c.get("terms", {}).items():
                    for side in ("U", "V"):
                        o = T.get(side, {})
                        rows.append(f"| {s['step']} | {term} | {side} | {o.get('common_strips')} | {fmt(o.get('unshaped_L2_rel'))} | "
                                    f"{fmt(o.get('shaped_L2_rel_Tp100'))} | {fmt(o.get('shaped_L2_rel_Tp200'))} | "
                                    f"{fmt(o.get('window_charge_rel_diff'), '+.3g')} | {fmt(T.get('UV_charge_ratio_rel_diff'), '+.3g')} | "
                                    f"{'same' if o.get('lead_identical') else str(o.get('lead_ref')) + '->' + str(o.get('lead_cand'))} | "
                                    f"{fmt(o.get('max_strip_charge_diff_over_lead'))} | {fmt(o.get('ref_charge_on_ref_only_strips_frac'))} / "
                                    f"{fmt(o.get('cand_charge_on_cand_only_strips_frac'))} |")
        if rows:
            L += ["", "| comparison | term | plane | common strips | unshaped L2 | shaped L2 100 | shaped L2 200 | ΔQ/Q | Δ(U/V) | lead | max strip ΔQ/lead | charge on ref-only / cand-only strips |",
                  "|---|---|---|---|---|---|---|---|---|---|---|---|"] + rows
        nonok = [(k, x) for k, v in m["comparisons"].items() for x in v["non_ok_checks"]]
        if nonok:
            L += ["", "Non-OK checks:", ""] + [f"* {k}: {x['check']} = {x['value']} (limit {x['limit']}) -> {x['status']}" for k, x in nonok[:60]]
        L.append("")
        res["muons"][E.name] = m
        ok = (e1_status in (None, "PASS")) and e2_status in ("OK", "WARN", "PASS") and integrity in (None, "PASS")
        overall_ok &= ok
    res["overall"] = "PASS" if overall_ok and res["muons"] else "NOT PASSED"
    L.insert(2, f"**overall: {res['overall']}** (E1 must PASS; E2 must have no FAIL/INCONCLUSIVE; WARN listed for inspection)")
    (X / "results.json").write_text(json.dumps(res, indent=1, default=str) + "\n")
    (X / "results.md").write_text("\n".join(L) + "\n")
    print("\n".join(L[:3]))
    for k, v in res["muons"].items():
        print(f"  {k}: E1 {v['E1']}  E2 {v['E2']}  references {v['reference_integrity']}")
    return 0 if res["overall"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
