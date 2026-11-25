#!/usr/bin/env python3
"""Aggregate a recenter_experiment_* directory into results.json / results.md with explicit COMPLETENESS and acceptance.

  k2a_recenter_aggregate.py EXPERIMENT_DIR

Inputs: env_check.json (mode, muons, environment status), thresholds.json, muon_*/steps.tsv and the cmp_*.json each
comparison step names. Required content per mode (anything missing -> INCONCLUSIVE, never PASS):
  refs  E1:  exactly E1_shift_1 and E1_shift_-1, both PASS (FAIL if either FAILs, else INCONCLUSIVE)
        E2:  E2_groups, E2_S0_local_census, E2_K2A_local, E2_K2Alocal_vs_S0local, E2_S0local_vs_S0, E2_K2Alocal_vs_S0,
             E2_K2Abase_vs_S0, E2_K2Alocal_vs_K2Abase
        reference integrity: sha256 lists before and after present, non-empty and identical (else FAIL / INCONCLUSIVE)
  far   E1:  not applicable (any E1 step present -> INCONCLUSIVE: wrong mode)
        base_centred_S0: OK, or EXPECTED_REFUSAL (the frozen base-centred replay refuses far-D events; NOT a failure);
             FAIL -> FAIL; anything else -> INCONCLUSIVE
        E2:  E2_groups, E2_S0_local_census, E2_K2A_local, E2_K2Alocal_vs_S0local, E2_truncation, E2_far_analysis
             (with truncation.json / far_analysis.json whose status equals the step's), and, only if base_centred_S0
             is OK, E2_S0local_vs_S0 and E2_K2Alocal_vs_S0; archive sha256 before == after; one E2_replay_<g>
             step per group of groups.json, OK, with a timing_attempts.json (late-step ladder) whose status is OK
Every comparison step must point to an existing cmp_*.json whose own status equals the step's status (missing report ->
INCONCLUSIVE; disagreement -> FAIL). Every muon listed in env_check.json must have a muon_* directory with steps.tsv.
Duplicate step names -> INCONCLUSIVE. Per muon: FAIL if anything FAILs, else INCONCLUSIVE if anything is missing or
inconclusive, else PASS (E2 WARNs are allowed and listed). Overall: FAIL > INCONCLUSIVE > PASS; the environment check
must be PASS. Exit: 0 only for overall PASS; 1 FAIL; 3 INCONCLUSIVE.
"""
import json
import sys
from pathlib import Path

ORDER = {"PASS": 0, "OK": 0, "WARN": 1, "EXPECTED_REFUSAL": 1, "INCONCLUSIVE": 2, "FAIL": 3}
E2_REFS = ["E2_groups", "E2_S0_local_census", "E2_K2A_local", "E2_K2Alocal_vs_S0local", "E2_S0local_vs_S0",
           "E2_K2Alocal_vs_S0", "E2_K2Abase_vs_S0", "E2_K2Alocal_vs_K2Abase"]
E2_FAR = ["E2_groups", "E2_S0_local_census", "E2_K2A_local", "E2_K2Alocal_vs_S0local", "E2_truncation", "E2_far_analysis"]
JSON_STEPS = {"E2_truncation": "truncation.json", "E2_far_analysis": "far_analysis.json"}
E2_FAR_WITH_BASE = ["E2_S0local_vs_S0", "E2_K2Alocal_vs_S0"]
E1_REQ = ["E1_shift_1", "E1_shift_-1"]


def fmt(x, spec=".3g"):
    return "n/a" if x is None else (f"{x:{spec}}" if isinstance(x, (int, float)) and not isinstance(x, bool) else str(x))


def reduce3(sts):
    """FAIL > INCONCLUSIVE > PASS (OK/PASS/WARN count as PASS)."""
    if "FAIL" in sts:
        return "FAIL"
    if "INCONCLUSIVE" in sts or not sts:
        return "INCONCLUSIVE"
    return "PASS"


def evaluate_muon(E, mode):
    problems = []          # (status, text)
    steps = []
    if not (E / "steps.tsv").exists():
        return {"status": "INCONCLUSIVE", "problems": [("INCONCLUSIVE", "steps.tsv missing")], "steps": [], "comparisons": {}}
    for line in (E / "steps.tsv").read_text().splitlines():
        if not line.strip():
            continue
        parts = (line.split("\t") + ["", "", "", ""])[:4]
        steps.append({"step": parts[0], "status": parts[1], "detail": parts[2], "json": parts[3]})
    names = [s["step"] for s in steps]
    for n in sorted({n for n in names if names.count(n) > 1}):
        problems.append(("INCONCLUSIVE", f"duplicate step {n}"))
    by = {s["step"]: s for s in steps}
    for s in steps:
        if s["status"] not in ORDER:
            problems.append(("FAIL", f"{s['step']}: unknown status {s['status']!r}"))
        elif s["status"] == "EXPECTED_REFUSAL" and not (mode == "far" and s["step"] == "base_centred_S0"):
            problems.append(("FAIL", f"{s['step']}: EXPECTED_REFUSAL is only valid for base_centred_S0 in far mode"))
        elif s["status"] in ("FAIL", "INCONCLUSIVE"):
            problems.append((s["status"], f"{s['step']}: {s['status']} - {s['detail'][:200]}"))
    # E1
    e1_steps = [n for n in names if n.startswith("E1_")]
    if mode == "refs":
        for n in E1_REQ:
            if n not in by:
                problems.append(("INCONCLUSIVE", f"required E1 step {n} missing"))
        extra = sorted(set(e1_steps) - set(E1_REQ))
        if extra:
            problems.append(("INCONCLUSIVE", f"unexpected E1 steps {extra}"))
        e1 = ("FAIL" if any(by[n]["status"] == "FAIL" for n in E1_REQ if n in by)
              else "PASS" if all(n in by and by[n]["status"] == "PASS" for n in E1_REQ) and not extra else "INCONCLUSIVE")
        if e1 == "INCONCLUSIVE" and all(n in by for n in E1_REQ):
            problems.append(("INCONCLUSIVE", "E1 not PASS on both relabellings: " + ", ".join(f"{n} {by[n]['status']}" for n in E1_REQ)))
    else:
        e1 = "NOT_APPLICABLE"
        if e1_steps:
            problems.append(("INCONCLUSIVE", f"E1 steps present in far mode: {e1_steps}"))
    # E2 completeness
    if mode == "refs":
        req = list(E2_REFS)
    else:
        req = list(E2_FAR)
        b = by.get("base_centred_S0")
        if b is None:
            problems.append(("INCONCLUSIVE", "base_centred_S0 step missing"))
        elif b["status"] == "OK":
            req += E2_FAR_WITH_BASE
        elif b["status"] != "EXPECTED_REFUSAL":
            problems.append((b["status"] if b["status"] == "FAIL" else "INCONCLUSIVE", f"base_centred_S0: {b['status']}"))
        for n, f in JSON_STEPS.items():
            if n in by:
                jp = E / f
                if not jp.exists():
                    problems.append(("INCONCLUSIVE", f"{n}: report {f} missing"))
                else:
                    try:
                        st = json.loads(jp.read_text()).get("status")
                    except json.JSONDecodeError:
                        st = None
                    if st != by[n]["status"]:
                        problems.append(("FAIL", f"{n}: step status {by[n]['status']} != {f} status {st}"))
        ab, aa = E / "archive_sha256_before.txt", E / "archive_sha256_after.txt"
        if not (ab.exists() and aa.exists() and ab.read_text().strip()):
            problems.append(("INCONCLUSIVE", "archive sha256 before/after missing"))
        elif ab.read_text() != aa.read_text():
            problems.append(("FAIL", "archive file changed during the experiment"))
    if mode == "far":                                 # every group: a ladder replay step whose record says OK
        gj = E / "groups" / "groups.json"
        if gj.exists():
            for g in json.loads(gj.read_text())["groups"]:
                n = f"E2_replay_{g}"
                req.append(n)
                ta = E / "S0groups" / g / "timing_attempts.json"
                if n in by and by[n]["status"] == "OK":
                    if not ta.exists():
                        problems.append(("INCONCLUSIVE", f"{n}: timing_attempts.json missing"))
                    elif json.loads(ta.read_text()).get("status") != "OK":
                        problems.append(("FAIL", f"{n}: step OK but timing ladder status {json.loads(ta.read_text()).get('status')}"))
        elif "E2_groups" in by:
            problems.append(("INCONCLUSIVE", "groups/groups.json missing"))
    for n in req:
        if n not in by:
            problems.append(("INCONCLUSIVE", f"required E2 step {n} missing"))
    # comparison reports
    comps = {}
    for s in steps:
        is_cmp = s["step"].startswith("E1_") or "_vs_" in s["step"]
        if not is_cmp:
            continue
        jp = Path(s["json"]) if s["json"] else None
        if jp is None or jp.name != f"cmp_{s['step']}.json" or not jp.exists():
            if s["status"] not in ("FAIL", "INCONCLUSIVE") or s["step"].startswith("E2"):
                problems.append(("INCONCLUSIVE", f"{s['step']}: comparison report missing ({s['json'] or 'none'})"))
            continue
        try:
            c = json.loads(jp.read_text())
        except json.JSONDecodeError:
            problems.append(("FAIL", f"{s['step']}: unreadable report {jp}")); continue
        if c.get("status") != s["status"]:
            problems.append(("FAIL", f"{s['step']}: step status {s['status']} != report status {c.get('status')}"))
        comps[s["step"]] = c
    # integrity
    b, a = E / "reference_sha256_before.txt", E / "reference_sha256_after.txt"
    if mode == "refs":
        if not (b.exists() and a.exists()):
            integrity = "INCONCLUSIVE"; problems.append(("INCONCLUSIVE", "reference sha256 lists missing"))
        elif not b.read_text().strip():
            integrity = "INCONCLUSIVE"; problems.append(("INCONCLUSIVE", "reference sha256 list empty"))
        elif a.read_text() != b.read_text():
            integrity = "FAIL"; problems.append(("FAIL", "reference files changed during the experiment"))
        else:
            integrity = "PASS"
    else:
        integrity = "NOT_APPLICABLE"
    e2_sts = [by[n]["status"] for n in req if n in by]
    e2 = ("FAIL" if "FAIL" in e2_sts else "INCONCLUSIVE" if ("INCONCLUSIVE" in e2_sts or len(e2_sts) < len(req))
          else "WARN" if "WARN" in e2_sts else "OK")
    status = reduce3([p[0] for p in problems] + ([] if e1 in ("PASS", "NOT_APPLICABLE") else [e1])
                     + ([] if e2 in ("OK", "WARN") else [e2]) + ([] if integrity in ("PASS", "NOT_APPLICABLE") else [integrity])
                     + ["PASS"])
    return {"status": status, "E1": e1, "E2": e2, "reference_integrity": integrity, "required_E2": req,
            "problems": problems, "steps": steps, "comparisons": comps}


def main():
    X = Path(sys.argv[1])
    envp = X / "env_check.json"
    env = json.loads(envp.read_text()) if envp.exists() else {"status": "missing"}
    mode = env.get("mode")
    th = json.loads((X / "thresholds.json").read_text()) if (X / "thresholds.json").exists() else None
    res = {"experiment": str(X), "mode": mode, "environment": env, "thresholds": th, "muons": {}, "problems": []}
    if env.get("status") != "PASS":
        res["problems"].append(("FAIL", "environment check not PASS"))
    if mode not in ("refs", "far"):
        res["problems"].append(("INCONCLUSIVE", f"unknown mode {mode!r}"))
    if th is None:
        res["problems"].append(("INCONCLUSIVE", "thresholds.json missing"))
    expected = sorted(env.get("provenance", {}).keys())
    present = sorted(p.name[5:] for p in X.glob("muon_*") if p.is_dir())
    if not expected:
        res["problems"].append(("INCONCLUSIVE", "no muons listed in env_check.json"))
    for m in expected:
        if m not in present:
            res["problems"].append(("INCONCLUSIVE", f"muon {m}: no output directory"))
    for m in present:
        if m not in expected:
            res["problems"].append(("INCONCLUSIVE", f"muon {m}: output not listed in env_check.json"))
    L = [f"# Recentering experiment results: {X.name} (mode {mode})", "", "",
         f"environment: {env.get('status')} (LUT {str(env.get('lut_sha256', ''))[:16]}, replay {str(env.get('replay_sha256', ''))[:16]}, "
         f"kernel {str(env.get('kernel_id', ''))[:16]}, half-strips {env.get('half_strips')}); "
         f"source {(X / 'SOURCE_GIT_SHA.txt').read_text().strip() if (X / 'SOURCE_GIT_SHA.txt').exists() else '?'}", ""]
    statuses = [p[0] for p in res["problems"]]
    for m in present:
        E = X / f"muon_{m}"
        r = evaluate_muon(E, mode) if mode in ("refs", "far") else {"status": "INCONCLUSIVE", "problems": [], "steps": [], "comparisons": {}}
        statuses.append(r["status"])
        res["muons"][E.name] = {k: v for k, v in r.items() if k != "comparisons"}
        res["muons"][E.name]["comparisons"] = {k: {"status": c.get("status"), "kind": c.get("kind"), "worst_L2_rel": c.get("worst_L2_rel"),
                                                   "error": c.get("error"), "non_ok_checks": [x for x in c.get("checks", []) if x["status"] != "OK"]}
                                               for k, c in r["comparisons"].items()}
        L += [f"## {E.name}: {r['status']}", "", f"E1 {r.get('E1')} | E2 {r.get('E2')} | reference files unchanged: {r.get('reference_integrity')}", ""]
        if r["problems"]:
            L += ["Problems:", ""] + [f"* {s}: {t}" for s, t in r["problems"]] + [""]
        L += ["| step | status | detail |", "|---|---|---|"] + [f"| {s['step']} | {s['status']} | {s['detail'][:180]} |" for s in r["steps"]]
        rows = []
        for name, c in r["comparisons"].items():
            for term, T in c.get("terms", {}).items():
                for side in ("U", "V"):
                    o = T.get(side, {})
                    rows.append(f"| {name} | {term} | {side} | {o.get('common_strips')} | {fmt(o.get('unshaped_L2_rel'))} | "
                                f"{fmt(o.get('shaped_L2_rel_Tp100'))} | {fmt(o.get('shaped_L2_rel_Tp200'))} | "
                                f"{fmt(o.get('window_charge_rel_diff'), '+.3g')} | {fmt(T.get('UV_charge_ratio_rel_diff'), '+.3g')} | "
                                f"{'same' if o.get('lead_identical') else str(o.get('lead_ref')) + '->' + str(o.get('lead_cand'))} | "
                                f"{fmt(o.get('max_strip_charge_diff_over_lead'))} | {fmt(o.get('ref_charge_on_ref_only_strips_frac'))} / "
                                f"{fmt(o.get('cand_charge_on_cand_only_strips_frac'))} |")
        if rows:
            L += ["", "| comparison | term | plane | common strips | unshaped L2 | shaped L2 100 | shaped L2 200 | ΔQ/Q | Δ(U/V) | lead | max strip ΔQ/lead | charge on ref-only / cand-only strips |",
                  "|---|---|---|---|---|---|---|---|---|---|---|---|"] + rows
        nonok = [(k, x) for k, v in res["muons"][E.name]["comparisons"].items() for x in v["non_ok_checks"]]
        if nonok:
            L += ["", "Non-OK checks:", ""] + [f"* {k}: {x['check']} = {x['value']} (limit {x['limit']}) -> {x['status']}" for k, x in nonok[:80]]
        L.append("")
    res["overall"] = reduce3(statuses) if present else "INCONCLUSIVE"
    L[2] = (f"**overall: {res['overall']}** (refs: E1 both relabellings PASS, all E2 steps and reports present, no FAIL/"
            "INCONCLUSIVE, references unchanged; WARN listed for inspection)")
    if res["problems"]:
        L[3:3] = ["Experiment-level problems:", ""] + [f"* {s}: {t}" for s, t in res["problems"]] + [""]
    (X / "results.json").write_text(json.dumps(res, indent=1, default=str) + "\n")
    (X / "results.md").write_text("\n".join(L) + "\n")
    print(L[0]); print(L[2])
    for p in res["problems"]:
        print(f"  {p[0]}: {p[1]}")
    for k, v in res["muons"].items():
        print(f"  {k}: {v['status']}  E1 {v.get('E1')}  E2 {v.get('E2')}  references {v.get('reference_integrity')}"
              + (f"  ({len(v['problems'])} problems)" if v["problems"] else ""))
    return {"PASS": 0, "FAIL": 1, "INCONCLUSIVE": 3}[res["overall"]]


if __name__ == "__main__":
    raise SystemExit(main())
