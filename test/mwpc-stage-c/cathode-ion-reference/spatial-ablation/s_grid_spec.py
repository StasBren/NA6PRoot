#!/usr/bin/env python3
"""Choose the Green-grid radii, wires and D range from the actual ion births of one or more path exports.

  s_grid_spec.py EXPORT [EXPORT ...] [--nr 14] [--half-strips 5]
Prints the wires present, the ion-birth radius distribution (charge weighted) and the strip-column range, and the
arguments for p3_green_grid_paths (--wires, --r-um) and s_green_table.py (--d-max-mm). Radii: log-spaced from just
inside the smallest to just outside the largest birth radius.
"""
import argparse
import csv
from pathlib import Path

import numpy as np

R_WIRE_UM = 15.0
PITCH = 1.7


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("exports", type=Path, nargs="+")
    ap.add_argument("--nr", type=int, default=14)
    ap.add_argument("--half-strips", type=int, default=5)
    a = ap.parse_args()
    R, Q, K, DU = [], [], [], []
    for d in a.exports:
        seal = dict(l.split("=", 1) for l in (d / "complete.txt").read_text().split() if "=" in l)
        rows = list(csv.DictReader(open(d / "ions.csv", newline="")))
        V = np.memmap(d / "vertices.f64", dtype="<f8", mode="r", shape=(int(seal["vertices"]), 5))
        p0 = np.asarray(V[np.array([int(r["first_vertex"]) for r in rows])])
        k = np.rint(p0[:, 1] / 3).astype(int)
        R.append(1000 * np.hypot(p0[:, 1] - 3 * k, p0[:, 2])); K.append(k); DU.append(p0[:, 3:5])
        Q.append(np.array([float(r["weight"]) for r in rows]))
    r, q, k, dcol = np.concatenate(R), np.concatenate(Q), np.concatenate(K), np.concatenate(DU)
    o = np.argsort(r); cq = np.cumsum(q[o]) / q.sum()
    qs = {p: float(r[o][np.searchsorted(cq, p)]) for p in (0.001, 0.01, 0.5, 0.99, 0.999)}
    lo, hi = max(R_WIRE_UM + 0.05, 0.995 * r.min()), 1.005 * r.max()
    radii = np.geomspace(lo, hi, a.nr)
    dmax = float(np.abs(dcol).max() + a.half_strips * PITCH + 0.5)
    wires = sorted(set(k.tolist()))
    print(f"exports={len(a.exports)} rows={len(r)} wires={wires} charge by wire="
          + ", ".join(f"{w_}:{q[k == w_].sum() / q.sum():.4f}" for w_ in wires))
    print(f"birth radius um: min {r.min():.3f}  charge quantiles {qs}  max {r.max():.3f}")
    print(f"strip column |D| max {np.abs(dcol).max():.3f} mm -> --d-max-mm {dmax:.2f}")
    print(f"--wires {','.join(map(str, wires))} --r-um {','.join(f'{x:.4f}' for x in radii)}")
    print(f"log-radius step factor {radii[1] / radii[0]:.4f}")


if __name__ == "__main__":
    main()
