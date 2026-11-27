#!/usr/bin/env python3
"""Read-only early-time decomposition of the K2A residual with the EXISTING S-ladder variants (no simulation).

  k2a_closure_anatomy.py NEW_DIR --sladder RUN --tag TAG [--variants G,SR,S3,S3A,SW,S3R,K2A,K2Akern]
  k2a_closure_anatomy.py NEW_DIR --layout DIR [--variants ...]      DIR/<muon>/S0, DIR/<muon>/<variant>[_TAG] (e.g. a
                                                                    per-group synthesis made by k2a_group_synth.sh)

The S-ladder variants of s_synthesize.py are superpositions of EXACT single-ion Green responses (frozen RKF path +
frozen weighting LUT), each with ONE spatial simplification, so their residuals separate the error sources:
  G       vs S0  interpolation of the Green table only (bilinear log r, theta; linear D; tau as K2A), ion positions exact
  SR      vs G   radial closure: every ion of a cloud at the cloud's r_eff = exp<log r>
  S3 / S3A vs G  angular closure: harmonics n <= 2 (+ hemisphere term)
  SW      vs G   longitudinal closure: each cloud at its mean D
  K2A     vs G   all closures together (the K2A description evaluated by synthesis)
  K2Akern vs K2A the harmonic-kernel implementation of the same description (numerical identity expected)
For each pair, ion only, per plane: L2 and the share of the squared residual in early windows
0-10-20-50-100-200-500-1000-7000 ns, the best first-order time shift and the high-frequency share (from
k2a_residual_diagnostics.anatomy). Decomposition is descriptive (errors need not add); it says which simplification
produces the early residual. Writes NEW_DIR/closure_anatomy.json and closure_anatomy.md.
"""
import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "validation"))
import k2a_residual_diagnostics as rd  # noqa: E402
import k2a_time_grid as tg             # noqa: E402

WIN = [0, 10, 20, 50, 100, 200, 500, 1000, 7000]
PAIRS = [("G", "S0"), ("SR", "G"), ("S3", "G"), ("S3A", "G"), ("SW", "G"), ("S3R", "G"), ("K2A", "G"), ("K2A", "S0"),
         ("K2Akern", "K2A"), ("K2Akern", "S0")]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out", type=Path); ap.add_argument("--sladder", type=Path); ap.add_argument("--tag", default="")
    ap.add_argument("--layout", type=Path); ap.add_argument("--variants", default="G,SR,S3,S3A,SW,S3R,K2A,K2Akern")
    ap.add_argument("--window-ns", type=float, default=7000.0)
    a = ap.parse_args()
    if a.out.exists():
        raise SystemExit(f"refusing: {a.out} exists")
    root = a.sladder or a.layout
    if root is None:
        raise SystemExit("give --sladder RUN --tag TAG or --layout DIR")
    want = a.variants.split(",")
    h100 = rd.shaper(100.0)
    res = {"tool": "k2a_closure_anatomy.py", "root": str(root), "tag": a.tag, "windows_ns": WIN, "muons": {}}
    for M in sorted(p for p in root.iterdir() if p.is_dir() and (p / "S0" / "replay" / "waveforms.csv").is_file()):
        var = {"S0": M / "S0"}
        for v in want:
            for d in (M / f"{v}_{a.tag}", M / v):
                if (d / "replay" / "waveforms.csv").is_file():
                    var[v] = d; break
        for v, d in list(var.items()):                          # synthesis variants may lack a prompt file: ion only here
            if not (d / "prompt" / "prompt_waveforms.csv").is_file():
                var.pop(v)
        ok, grid, nb = tg.check(list(var.items()), a.window_ns)
        rec = {"variants": {k: str(v) for k, v in var.items()}, "time_grid": grid["status"], "pairs": {}}
        if ok:
            data = {k: rd.load(d) for k, d in var.items()}
            t = data["S0"][1][:nb]; dt = float(t[1] - t[0])
            for c, r in PAIRS:
                if c not in data or r not in data:
                    continue
                lr, _, ir, _ = data[r]; lc, _, ic, _ = data[c]
                rec["pairs"][f"{c}_vs_{r}"] = {}
                for s in ("U", "V"):
                    com = [l for l in lr if l.startswith(s + "_") and l in lc]
                    R0 = ir[:nb][:, [lr.index(l) for l in com]]; C0 = ic[:nb][:, [lc.index(l) for l in com]]
                    rec["pairs"][f"{c}_vs_{r}"][s] = rd.anatomy(R0, C0, t, dt, h100, WIN)
        else:
            rec["failures"] = grid["failures"]
        res["muons"][M.name] = rec
    (a.out).mkdir(parents=True)
    (a.out / "closure_anatomy.json").write_text(json.dumps(res, indent=1, default=str) + "\n")
    pct = lambda v: "n/a" if v is None or not math.isfinite(v) else f"{100 * v:.3g} %"
    L = ["# Early-time closure anatomy from existing S-ladder variants (ion only; read-only)", "", f"root `{root}` tag `{a.tag}`", "",
         "| muon | pair | plane | L2 | share 0-10 ns | 10-20 | 20-50 | 50-100 | >100 ns | best shift ns | L2 after shift | high-freq share | shaped/unshaped |",
         "|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for m, rec in res["muons"].items():
        for pn, pv in rec["pairs"].items():
            for s, x in pv.items():
                if not x:
                    continue
                w = x["window_share"]
                late = sum(v for k, v in w.items() if int(k.split("-")[0]) >= 100)
                L.append(f"| {m} | {pn} | {s} | {pct(x['L2_rel'])} | {pct(w['0-10'])} | {pct(w['10-20'])} | {pct(w['20-50'])} | {pct(w['50-100'])} | "
                         f"{pct(late)} | {x['best_time_shift_ns']:+.3g} | {pct(x['L2_after_time_shift'])} | {pct(x['high_freq_share_50ns'])} | "
                         f"{x['shaped_over_unshaped_L2']:.3g} |")
        if rec.get("failures"):
            L.append(f"| {m} | time grid FAIL | | {'; '.join(rec['failures'][:2])} |")
    L += ["", "Reading: early error in G vs S0 = Green-table interpolation; in SR vs G = single-radius closure; in S3/S3A vs G = "
          "angular truncation; in SW vs G = longitudinal collapse. K2Akern vs K2A should be ~0 (implementation identity)."]
    (a.out / "closure_anatomy.md").write_text("\n".join(L) + "\n")
    print(f"closure anatomy: {len(res['muons'])} muons -> {a.out / 'closure_anatomy.md'}")


if __name__ == "__main__":
    main()
