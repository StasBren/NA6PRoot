#!/usr/bin/env python3
"""SYNTHETIC path export in the exporter's format, for pipeline tests of the T-ladder tools ONLY (not NA60+ data).

Ion paths are straight lines from a birth point near an anode to the cathode with a 1/r-field time law
r^2 = r0^2 (1 + (t - t_b)/tau0); no RKF, no real mobility. Birth times per primary mimic the V0 numbers only roughly.
"""
import argparse
from pathlib import Path

import numpy as np

R_WIRE = 0.015


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out", type=Path)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--primaries", type=int, default=50)
    ap.add_argument("--ions", type=int, default=24)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    a.out.mkdir(parents=True, exist_ok=False)
    rows, blocks, off = [], [], 0
    w_track = rng.uniform(-0.8, 0.8)
    tmu = rng.normal(26, 9)
    for p in range(a.primaries):
        wire = rng.choice([-1, 0, 1], p=[0.2, 0.6, 0.2])
        tp = max(1.0, tmu + rng.normal(0, 10.6))
        gain = rng.lognormal(np.log(3e4), 0.9)
        side = rng.choice([-1, 1])
        for i in range(a.ions):
            th = rng.normal(side * np.pi / 2, 1.0)
            r0 = R_WIRE + rng.exponential(0.01) + 0.002
            x0 = np.array([3 * wire + r0 * np.cos(th), r0 * np.sin(th)])
            vend = 3.0 if np.sin(th) >= 0 else -3.0
            uend = 3 * wire + 1.2 * np.cos(th)
            xe = np.array([uend, vend])
            f = np.linspace(0, 1, 160) ** 2
            pts = x0[None, :] + f[:, None] * (xe - x0)[None, :]
            rho = np.hypot(pts[:, 0] - 3 * np.rint(pts[:, 0] / 3), pts[:, 1])
            if np.any(rho <= R_WIRE):
                continue
            dist = np.hypot(pts[:, 0] - 3 * wire, pts[:, 1])
            dist = np.maximum.accumulate(dist) + 1e-9 * np.arange(len(dist))
            tb = tp + rng.normal(0, 0.09)
            t = tb + 3.0 * ((dist / dist[0]) ** 2 - 1.0)
            wv = w_track + rng.normal(0, 0.05)
            blk = np.column_stack([t, pts[:, 0], pts[:, 1], np.full(len(t), wv), np.full(len(t), wv)])
            blk[-1, 2] = vend
            rows.append((len(rows), p, gain / a.ions, t[0], t[-1], off, len(t)))
            blocks.append(blk)
            off += len(t)
    V = np.concatenate(blocks)
    V.astype("<f8").tofile(a.out / "vertices.f64")
    wsum = sum(r[2] for r in rows)
    # make weights sum to an integer exact-ion census
    exact = int(round(wsum))
    scale = exact / wsum
    with open(a.out / "ions.csv", "w") as f:
        f.write("row,primary,weight,birth_ns,end_ns,first_vertex,vertices,upper_index,lower_index\n")
        for r in rows:
            f.write(f"{r[0]},{r[1]},{float(r[2]*scale)!r},{float(r[3])!r},{float(r[4])!r},{r[5]},{r[6]},0,0\n")
    (a.out / "complete.txt").write_text(f"rows={len(rows)}\nexact_ions={exact}\nweight_sum={exact}\nvertices={len(V)}\nstep_mm=0.005\n")
    print(f"SYNTHETIC export {a.out}: rows={len(rows)} vertices={len(V)} exact_ions={exact}")


if __name__ == "__main__":
    main()
