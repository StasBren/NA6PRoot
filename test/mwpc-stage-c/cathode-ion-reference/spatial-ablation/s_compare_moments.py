#!/usr/bin/env python3
"""Per-cloud comparison of two cloud-moment CSVs (e.g. birth-derived vs T1_export-derived).

  s_compare_moments.py REFERENCE.csv CANDIDATE.csv [--report OUT.md] [--tol 1e-9]
Clouds are matched on (primary, wire); the sets must be identical. For every column reports max |diff| and RMS diff
(relative for N_fC and r_eff_um, absolute otherwise) and lists the worst clouds. PASS if every max difference is
below --tol (default 1e-9: rounding-level agreement; the two paths use different summation orders) and base_U/base_V
are identical. Exception: t_rms_ns (diagnostic only, not a kernel input) is sqrt(<t^2> - <t>^2) of a T1-collapsed
cloud, i.e. sqrt of a cancellation at the 1e-16 t^2 level, so it is only required to be < 1e-5 ns in both files.
"""
import argparse
import csv
from pathlib import Path

import numpy as np

COLS = [("N_fC", "rel"), ("r_eff_um", "rel"), ("t_ns", "abs"), ("t_rms_ns", "abs"), ("D_U_mm", "abs"), ("D_V_mm", "abs"),
        ("C1", "abs"), ("S1", "abs"), ("C2", "abs"), ("S2", "abs"), ("A", "abs"), ("min_clearance_um", "abs")]


def load(p):
    rows = list(csv.DictReader(open(p, newline="")))
    return {(int(r["primary"]), int(r["wire"])): r for r in rows}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("reference", type=Path); ap.add_argument("candidate", type=Path)
    ap.add_argument("--report", type=Path); ap.add_argument("--tol", type=float, default=1e-9)
    a = ap.parse_args()
    R, C = load(a.reference), load(a.candidate)
    L = [f"# Cloud-moment comparison", "", f"reference `{a.reference}`, candidate `{a.candidate}`", ""]
    ok = set(R) == set(C)
    L.append(f"clouds: reference {len(R)}, candidate {len(C)}, identical keys: {ok}")
    if not ok:
        L.append(f"only in reference: {sorted(set(R) - set(C))[:20]}; only in candidate: {sorted(set(C) - set(R))[:20]}")
    keys = sorted(set(R) & set(C))
    base_ok = all(R[k]["base_U"] == C[k]["base_U"] and R[k]["base_V"] == C[k]["base_V"] for k in keys)
    L += [f"strip frame base_U/base_V identical: {base_ok}", "", "| quantity | kind | max |diff| | RMS diff | worst cloud |",
          "|---|---|---|---|---|"]
    ok &= base_ok
    for c, kind in COLS:
        r = np.array([float(R[k][c]) for k in keys]); x = np.array([float(C[k][c]) for k in keys])
        d = (x - r) / np.where(r != 0, np.abs(r), 1) if kind == "rel" else x - r
        i = int(np.argmax(np.abs(d)))
        L.append(f"| {c} | {kind} | {np.abs(d).max():.3g} | {np.sqrt(np.mean(d ** 2)):.3g} | {keys[i]} |")
        if c == "t_rms_ns":
            ok &= bool(np.abs(r).max() < 1e-5 and np.abs(x).max() < 1e-5)
        else:
            ok &= bool(np.abs(d).max() < a.tol)
    L += ["", ("PASS" if ok else "FAIL") + f" cloud moments agree (tol {a.tol:g})"]
    txt = "\n".join(L) + "\n"
    if a.report:
        a.report.write_text(txt)
    print(txt)
    raise SystemExit(0 if ok else 1)


if __name__ == "__main__":
    main()
