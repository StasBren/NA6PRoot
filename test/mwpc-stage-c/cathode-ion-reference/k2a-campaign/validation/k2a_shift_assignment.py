#!/usr/bin/env python3
"""Compare the per-primary strip shifts of the K2A production rule with the S0 grouping.

  k2a_shift_assignment.py S0_GROUPS_JSON K2A_RULE_DIR OUT.json

S0_GROUPS_JSON  k2a_group_export.py groups.json ("primary_shift": primary -> [m_U, m_V], from the export's first-vertex
                positions, ion-weighted)
K2A_RULE_DIR    k2a_strip_eval output of a PRODUCTION-RULE run (no --shifts): shifts.csv (primary, m_U, m_V, computed_*,
                imposed, Dbar_U_mm, Dbar_V_mm, margin_U, margin_V) and clouds.csv (primary, N_fC)
Reports every primary on which they differ (both assignments, Dbar, rule margins, the primary's charge and charge
fraction), primaries known to only one side, and near-tie primaries (rule margin < 1e-6 strip pitch fraction).
Exit 0 AGREE, 2 DISAGREE (reported, not an error), 3 INCONCLUSIVE (inputs missing), 1 FAIL (malformed input, or the
K2A run was not a production-rule run).
"""
import argparse
import csv
import json
import sys
from pathlib import Path


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("groups", type=Path); ap.add_argument("k2a", type=Path); ap.add_argument("out", type=Path)
    a = ap.parse_args()
    rep = {"tool": "k2a_shift_assignment.py", "s0_groups": str(a.groups), "k2a_rule_dir": str(a.k2a)}
    def done(status, msg, code):
        rep["status"] = status; rep["summary"] = msg
        a.out.write_text(json.dumps(rep, indent=1) + "\n"); print(f"{status} shift assignment: {msg}"); return code
    for p in (a.groups, a.k2a / "shifts.csv", a.k2a / "clouds.csv"):
        if not p.is_file():
            return done("INCONCLUSIVE", f"missing {p}", 3)
    try:
        s0 = {int(k): (int(v[0]), int(v[1])) for k, v in json.loads(a.groups.read_text())["primary_shift"].items()}
        rows = list(csv.DictReader(open(a.k2a / "shifts.csv")))
        k2a = {int(r["primary"]): r for r in rows}
        q = {}
        for r in csv.DictReader(open(a.k2a / "clouds.csv")):
            q[int(r["primary"])] = q.get(int(r["primary"]), 0.0) + float(r["N_fC"])
    except (KeyError, ValueError, IndexError) as e:
        return done("FAIL", f"malformed input: {e!r}", 1)
    if any(int(r["imposed"]) for r in rows):
        return done("FAIL", "K2A run has imposed shifts: not a production-rule run", 1)
    qt = sum(q.values()) or float("nan")
    dis, near = [], []
    for p in sorted(set(s0) & set(k2a)):
        r = k2a[p]; mk = (int(r["m_U"]), int(r["m_V"]))
        mu, mv = float(r["margin_U"]), float(r["margin_V"])
        if min(mu, mv) < 1e-6:
            near.append(p)
        if mk != s0[p]:
            dis.append({"primary": p, "s0_m": list(s0[p]), "k2a_m": list(mk), "Dbar_U_mm": float(r["Dbar_U_mm"]),
                        "Dbar_V_mm": float(r["Dbar_V_mm"]), "rule_margin_U": mu, "rule_margin_V": mv,
                        "charge_fC": q.get(p, 0.0), "charge_frac": q.get(p, 0.0) / qt})
    only_s0, only_k2a = sorted(set(s0) - set(k2a)), sorted(set(k2a) - set(s0))
    rep.update(primaries_s0=len(s0), primaries_k2a=len(k2a), primaries_common=len(set(s0) & set(k2a)),
               only_s0=only_s0, only_k2a=only_k2a, only_k2a_charge_frac=sum(q.get(p, 0.0) for p in only_k2a) / qt,
               disagreements=dis, n_disagree=len(dis), disagree_charge_frac=sum(d["charge_frac"] for d in dis),
               near_tie_primaries=near, total_cloud_charge_fC=qt)
    if dis or only_s0 or only_k2a:
        return done("DISAGREE", f"{len(dis)} of {len(set(s0) & set(k2a))} common primaries assigned differently "
                    f"({100 * rep['disagree_charge_frac']:.3g} % of the cloud charge); only S0 {len(only_s0)}, only K2A {len(only_k2a)}; "
                    f"near ties {len(near)}", 2)
    return done("AGREE", f"all {len(s0)} primaries have the same (m_U, m_V); near ties {len(near)}", 0)


if __name__ == "__main__":
    sys.exit(main())
