#!/usr/bin/env python3
"""Compare two replay-format variant directories at full numerical precision (ion term and electron step).

  s_compare_waveforms.py REF_DIR CAND_DIR [--tol 1e-9] [--report OUT.md]
REF_DIR / CAND_DIR contain replay/waveforms.csv and prompt/prompt_waveforms.csv (e.g. K2Akern_TAG from the T1 export
moments and K2Akern_TAG_births from the birth-record moments). Reports, per file and plane: identical labels and time
grid, max |diff| (fC/ns), max |diff| / max |ref|, and relative L2. PASS if every relative L2 < --tol.
"""
import argparse
from pathlib import Path

import numpy as np


def load(p):
    with open(p) as f:
        h = f.readline().strip().split(",")
    return h, np.loadtxt(p, delimiter=",", skiprows=1, ndmin=2)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("ref", type=Path); ap.add_argument("cand", type=Path)
    ap.add_argument("--tol", type=float, default=1e-9); ap.add_argument("--report", type=Path)
    a = ap.parse_args()
    ok, L = True, [f"# Waveform comparison `{a.cand}` vs `{a.ref}`", "",
                   "| file | plane | labels+grid | max abs diff [fC/ns] | max diff / max ref | rel L2 |", "|---|---|---|---|---|---|"]
    for sub, fn in (("replay", "waveforms.csv"), ("prompt", "prompt_waveforms.csv")):
        hr, R = load(a.ref / sub / fn); hc, C = load(a.cand / sub / fn)
        same = hr == hc and R.shape == C.shape and np.array_equal(R[:, 0], C[:, 0])
        ok &= same
        if not same:
            L.append(f"| {fn} | - | MISMATCH | | | |"); continue
        for s in ("U", "V"):
            cols = [i for i, l in enumerate(hr) if l.startswith(s + "_")]
            r, c = R[:, cols], C[:, cols]
            l2 = np.linalg.norm(c - r) / max(np.linalg.norm(r), 1e-300)
            ok &= bool(l2 < a.tol)
            L.append(f"| {fn} | {s} | identical | {np.abs(c - r).max():.3g} | {np.abs(c - r).max() / np.abs(r).max():.3g} | {l2:.3g} |")
    L += ["", ("PASS" if ok else "FAIL") + f" waveforms agree (rel L2 < {a.tol:g})"]
    txt = "\n".join(L) + "\n"
    if a.report:
        a.report.write_text(txt)
    print(txt)
    raise SystemExit(0 if ok else 1)


if __name__ == "__main__":
    main()
