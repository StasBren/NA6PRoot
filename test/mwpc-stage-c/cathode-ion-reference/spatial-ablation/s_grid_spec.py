#!/usr/bin/env python3
"""Choose the Green-grid radii, wires and D range from the actual ion births of one or more path exports.

  s_grid_spec.py EXPORT [EXPORT ...] [--nr 14] [--half-strips 5]
Prints the wires present, the ion-birth radius distribution (charge weighted) and the strip-column range, and the
arguments for p3_green_grid_paths (--wires, --r-um) and s_green_table.py (--d-max-mm).

Radii: log-spaced from r_lo to r_hi with
  r_lo = r_min - min(0.005 r_min, 0.5 (r_min - R_wire))   -> never above the innermost ion, never on/inside the anode
  r_hi = 1.005 r_max
so every observed birth radius lies inside the grid (no extrapolation) and the innermost grid start point keeps half
the innermost ion's clearance from the wire surface. Radii are printed at full precision (repr), so rounding cannot
move a node onto the wire or above r_min. Coverage is re-checked on the printed values; the script exits non-zero if an
ion is inside/on the anode or outside the printed grid.
(The previous rule, r_lo = max(R_wire + 0.05 um, 0.995 r_min), excluded ions within 0.05 um of the wire surface.)
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
    ap.add_argument("--strip-pitch-mm", type=float, default=PITCH, help="readout strip pitch (geometry study; default frozen 1.7)")
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
    qs = {p: round(float(r[o][np.searchsorted(cq, p)]), 6) for p in (0.001, 0.01, 0.5, 0.99, 0.999)}
    clear = r - R_WIRE_UM
    if clear.min() <= 0:
        raise SystemExit(f"FAIL {int(np.sum(clear <= 0))} ion births on or inside the anode (min clearance {clear.min():.3g} um)")
    lo = r.min() - min(0.005 * r.min(), 0.5 * clear.min())
    hi = 1.005 * r.max()
    radii = np.geomspace(lo, hi, a.nr)
    radii[0], radii[-1] = lo, hi                       # exact end points (geomspace rounding)
    rstr = ",".join(repr(float(x)) for x in radii)
    back = np.array([float(x) for x in rstr.split(",")])
    out = (r < back[0]) | (r > back[-1])
    dmax = float(np.abs(dcol).max() + a.half_strips * a.strip_pitch_mm + 0.5)
    wires = sorted(set(k.tolist()))
    print(f"exports={len(a.exports)} rows={len(r)} wires={wires} charge by wire="
          + ", ".join(f"{w_}:{q[k == w_].sum() / q.sum():.4f}" for w_ in wires))
    print(f"birth radius um: min {r.min():.9g}  charge quantiles {qs}  max {r.max():.6g}")
    cq_ = {p: float(f"{np.quantile(clear, p):.6g}") for p in (0.0, 0.01, 0.1, 0.5)}
    print(f"clearance from wire surface um: min {clear.min():.6g}  quantiles {cq_}; ions with clearance < 0.05 um: "
          f"{int(np.sum(clear < 0.05))} ({q[clear < 0.05].sum() / q.sum():.4%} of charge)")
    print(f"grid: r_lo {float(back[0])!r} um (clearance {back[0] - R_WIRE_UM:.6g} um), r_hi {float(back[-1])!r} um, "
          f"{a.nr} log-spaced nodes, step factor {radii[1] / radii[0]:.4f}")
    hist = np.histogram(r, bins=back, weights=q)[0] / q.sum()
    print("charge fraction per radial interval: " + ", ".join(f"[{back[i]:.4g},{back[i+1]:.4g}) {hist[i]:.4f}" for i in range(len(hist))))
    if out.any() or back[0] <= R_WIRE_UM:
        raise SystemExit(f"FAIL radial coverage: {int(out.sum())} ions outside [{float(back[0])!r}, {float(back[-1])!r}] or grid on the anode")
    print(f"PASS radial coverage: all {len(r)} ion births inside the grid; innermost node {back[0] - R_WIRE_UM:.6g} um above the wire")
    print(f"strip column |D| max {np.abs(dcol).max():.3f} mm -> --d-max-mm {dmax:.2f}")
    print(f"--wires {','.join(map(str, wires))} --r-um {rstr}")

if __name__ == "__main__":
    main()
