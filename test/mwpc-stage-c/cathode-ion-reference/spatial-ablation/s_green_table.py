#!/usr/bin/env python3
"""Single-ion Green's-function table from the grid paths of p3_green_grid_paths, with the FROZEN weighting LUT.

For every grid start point (wire k, radius r, azimuth theta), plane side s (U, V) and strip-relative longitudinal
coordinate D (the exporter's strip column minus strip offset x pitch, i.e. the `w` passed to WeightingLUT.upper in
the frozen replay's query()):
    G[s,k,r,theta,D,tau] = phi_s(x(t0+tau), D) - phi_s(x(t0), D)     ion term per unit charge, frozen convention (+q)
    PHI0[s,k,r,theta,D]  = phi_s(x(t0), D)                            integral of the avalanche-electron term
with x(t) the frozen RKF path (piecewise-linear in its vertices, exactly as query() interpolates it), x mirrored to
-x for side V as in query(). Exact symmetry used: the field is 2D (no motion along the wire), so a w translation of
an ion is a shift of D; the path itself is not changed. Nothing in the table moves a charge along a path that was
computed for a different start point.

Output: a directory with G.npy (float32 memmap), PHI0.npy, axes.json. Uses the replay's own read_export(), path_for()
and the LUT class (imported from --replay-dir, unchanged).
"""
import argparse
import csv
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np


def sha256(p, chunk=1 << 22):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(chunk), b""):
            h.update(b)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("table", type=Path, help="spectral weighting LUT .npz")
    ap.add_argument("grid", type=Path, help="output directory of p3_green_grid_paths")
    ap.add_argument("out", type=Path)
    ap.add_argument("--replay-dir", type=Path, required=True)
    ap.add_argument("--d-max-mm", type=float, default=10.54, help="|D| range; must cover |strip col| + half_strips x pitch")
    ap.add_argument("--d-sub", type=int, default=10, help="D nodes per strip pitch (pitch/d-sub spacing)")
    ap.add_argument("--tau-max-ns", type=float, default=8000.0)
    ap.add_argument("--ntau", type=int, default=200)
    a = ap.parse_args()
    sys.path.insert(0, str(a.replay_dir.resolve()))
    import p3_replay_full_spectral_muon as rp  # noqa: E402  frozen replay module
    from p3_spectral_weighting_lut import P, A, WeightingLUT  # noqa: E402

    rows, vertices, weights, exact = rp.read_export(a.grid)
    grid = list(csv.DictReader(open(a.grid / "grid.csv", newline="")))
    if len(grid) != len(rows):
        raise SystemExit("grid.csv does not match ions.csv")
    wires = sorted({int(g["wire"]) for g in grid})
    radii = sorted({float(g["r_um"]) for g in grid})
    thetas = sorted({round(float(g["theta_deg"]), 9) for g in grid})
    nk, nr, nt = len(wires), len(radii), len(thetas)
    if nk * nr * nt != len(rows):
        raise SystemExit("grid is not a full wire x r x theta product")
    dD = P / a.d_sub
    nD = int(np.floor(a.d_max_mm / dD))
    D = np.arange(-nD, nD + 1) * dD
    tau = np.r_[0.0, np.geomspace(0.01, a.tau_max_ns, a.ntau - 1)]
    lut = WeightingLUT(a.table)
    a.out.mkdir(parents=True, exist_ok=False)
    shape = (2, nk, nr, nt, len(D), len(tau))
    G = np.lib.format.open_memmap(a.out / "G.npy", mode="w+", dtype=np.float32, shape=shape)
    PHI0 = np.zeros(shape[:-1], dtype=np.float64)
    kidx = {k: i for i, k in enumerate(wires)}
    ridx = {r: i for i, r in enumerate(radii)}
    tidx = {t: i for i, t in enumerate(thetas)}
    umax = 0.0
    t0w = time.perf_counter()
    for j, (r, g) in enumerate(zip(rows, grid)):
        p = rp.path_for(r, vertices)
        x = np.column_stack([np.interp(p[0, 0] + tau, p[:, 0], p[:, i]) for i in (1, 2)])
        umax = max(umax, float(np.max(np.abs(p[:, 1]))))
        ik, ir, it = kidx[int(g["wire"])], ridx[float(g["r_um"])], tidx[round(float(g["theta_deg"]), 9)]
        for side in range(2):
            xx = x if side == 0 else -x
            phi = lut.upper(np.repeat(xx, len(D), axis=0), np.tile(D, len(tau))).reshape(len(tau), len(D))
            if phi.min() < -1e-5 or phi.max() > 1 + 1e-5:
                raise SystemExit(f"LUT maximum-principle excursion at grid row {j}")
            G[side, ik, ir, it] = (phi - phi[0]).T.astype(np.float32)
            PHI0[side, ik, ir, it] = phi[0]
        if (j + 1) % 200 == 0:
            print(f"grid paths {j + 1}/{len(rows)}  {time.perf_counter() - t0w:.0f} s", flush=True)
    if a.d_max_mm + A * umax > 12.0 - 1e-6:
        raise SystemExit("D range + stereo shift exceeds the LUT phase range")
    G.flush()
    np.save(a.out / "PHI0.npy", PHI0)
    axes = {"wires": wires, "r_um": radii, "theta_deg": thetas, "D_mm": D.tolist(), "dD_mm": dD, "d_sub": a.d_sub,
            "tau_ns": tau.tolist(), "pitch_mm": P, "stereo_tan": A, "shape": list(shape),
            "lut": str(a.table), "lut_sha256": sha256(a.table), "grid": str(a.grid),
            "grid_ions_sha256": sha256(a.grid / "ions.csv"), "grid_vertices_sha256": sha256(a.grid / "vertices.f64"),
            "grid_seal": (a.grid / "complete.txt").read_text(), "max_abs_u_mm": umax,
            "convention": "G = phi(x(t0+tau)) - phi(x(t0)), PHI0 = phi(x(t0)); x -> -x for side V; per unit charge",
            "build_s": time.perf_counter() - t0w}
    (a.out / "axes.json").write_text(json.dumps(axes, indent=1) + "\n")
    print(f"PASS green table {shape}  {G.nbytes / 2**20:.0f} MiB  {axes['build_s']:.0f} s")


if __name__ == "__main__":
    main()
