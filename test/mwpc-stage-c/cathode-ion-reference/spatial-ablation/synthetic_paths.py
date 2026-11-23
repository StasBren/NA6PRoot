#!/usr/bin/env python3
"""SYNTHETIC ion transport for pipeline tests of the S-ladder tools ONLY (not NA60+ physics, no Garfield).

A deterministic toy 'field': the path depends only on the start point (as for RKF in a static field). From (wire k,
r0, theta) the ion moves radially to r = 0.4 mm with r^2 = r0^2 (1 + t/tau0), then straight to the cathode at
u = 3k + 1.6 cos(theta), v = +-3 mm, continuing the same time law in distance from the wire.
  synthetic_paths.py grid  OUT --wires -1,0,1 --r-um 15.5,...  --ntheta 72     (same format as p3_green_grid_paths)
  synthetic_paths.py muon  OUT --seed 1                                         (T1-like timing: one time per primary)
"""
import argparse
from pathlib import Path

import numpy as np

TAU0 = 3.0
R1 = 0.4


def path(k, r0_mm, th, t0):
    c, s = np.cos(th), np.sin(th)
    if abs(s) < 1e-6:
        s = 1e-6
    rr = np.geomspace(r0_mm, R1, 60)
    p1 = np.column_stack([3 * k + rr * c, rr * s])
    end = np.array([3 * k + 1.6 * c, 3.0 * np.sign(s)])
    f = np.linspace(0, 1, 61)[1:] ** 1.5
    p2 = p1[-1] + f[:, None] * (end - p1[-1])
    pts = np.vstack([p1, p2])
    dist = np.hypot(pts[:, 0] - 3 * k, pts[:, 1])
    dist = np.maximum.accumulate(dist) + 1e-12 * np.arange(len(dist))
    t = t0 + TAU0 * ((dist / r0_mm) ** 2 - 1.0)
    t[0] = t0
    pts[-1, 1] = 3.0 * np.sign(s)
    return t, pts


def write(out, recs, extra=None):
    out.mkdir(parents=True, exist_ok=False)
    blocks, off = [], 0
    with open(out / "ions.csv", "w") as f:
        f.write("row,primary,weight,birth_ns,end_ns,first_vertex,vertices,upper_index,lower_index\n")
        for i, (prim, wt, t, pts, wU, wV) in enumerate(recs):
            blk = np.column_stack([t, pts, np.full(len(t), wU), np.full(len(t), wV)])
            f.write(f"{i},{prim},{float(wt)!r},{float(t[0])!r},{float(t[-1])!r},{off},{len(t)},0,0\n")
            blocks.append(blk); off += len(t)
    np.concatenate(blocks).astype("<f8").tofile(out / "vertices.f64")
    tot = sum(r[1] for r in recs)
    (out / "complete.txt").write_text(f"rows={len(recs)}\nexact_ions={int(round(tot))}\nweight_sum={tot!r}\nvertices={off}\nstep_mm=0.005\n")
    if extra:
        (out / "grid.csv").write_text(extra)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["grid", "muon"])
    ap.add_argument("out", type=Path)
    ap.add_argument("--wires", default="-1,0,1")
    ap.add_argument("--r-um", default="")
    ap.add_argument("--ntheta", type=int, default=72)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    if a.mode == "grid":
        recs, g = [], ["row,wire,r_um,theta_deg,u_mm,v_mm,end_u_mm,end_v_mm,end_ns"]
        for k in [int(x) for x in a.wires.split(",")]:
            for r in [float(x) for x in a.r_um.split(",")]:
                for m in range(a.ntheta):
                    th = (m + 0.5) * 2 * np.pi / a.ntheta
                    t, pts = path(k, r / 1000, th, 0.0)
                    g.append(f"{len(recs)},{k},{r!r},{th * 180 / np.pi!r},{pts[0,0]!r},{pts[0,1]!r},{pts[-1,0]!r},{pts[-1,1]!r},{t[-1]!r}")
                    recs.append((len(recs), 1.0, t, pts, 0.0, 0.0))
        write(a.out, recs, "\n".join(g) + "\n")
    else:
        rng = np.random.default_rng(a.seed)
        recs = []
        tmu = rng.normal(26, 9)
        wtrack = rng.uniform(-0.6, 0.6)
        for p in range(40):
            k = int(rng.choice([-1, 0, 1], p=[0.2, 0.6, 0.2]))
            tp = max(1.0, tmu + rng.normal(0, 10.6))
            side = rng.choice([-1, 1])
            n = 24
            gain = rng.lognormal(np.log(3e4), 0.9)
            for i in range(n):
                th = np.mod(rng.normal(side * np.pi / 2, 1.0), 2 * np.pi)
                r0 = 0.0152 + rng.exponential(0.012)
                t, pts = path(k, r0, th, tp)
                wv = wtrack + rng.normal(0, 0.05)
                recs.append((p, gain / n, t, pts, wv, wv))
        tot = sum(r[1] for r in recs)
        scale = round(tot) / tot
        recs = [(p, wt * scale, t, pts, wU, wV) for p, wt, t, pts, wU, wV in recs]
        write(a.out, recs)
    print(f"SYNTHETIC {a.mode} export {a.out}: rows={len(recs)}")


if __name__ == "__main__":
    main()
