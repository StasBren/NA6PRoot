#!/usr/bin/env python3
"""Independent numerical checks of a weighting LUT for its readout geometry (no archive, no transport).

  geom_lut_check.py GEOMETRY.json NUMW_DIR LUT.npz OUT.json [--refined-L 96 --refined-N 32] [--seed 7]

At held-out points (random, off every table node; near each anode 1.02 R .. 0.7 mm, and in the bulk) and phases:
  1. interpolation:  LUT.upper  vs the direct spectral solution with the LUT's own settings (L 48, kmax 6, N 24,
                     +4 images; strip pilot Correction, no table) -> table interpolation error only
  2. truncation:     direct solution with refined settings (default L 96 mm, N 32 sources) vs the LUT settings
  3. bounds:         0 <= phi <= 1 (maximum principle), phi = 0 on the anode surface, phi = indicator on the strip plane
  4. partition:      sum of contiguous strips j = -J..J (all inside the LUT phase range) vs the whole-plane (k = 0)
                     potential. The whole-plane potential does not depend on segmentation, so this is a geometry-
                     independent invariant; the finite-J tail is reported, compare it with the baseline LUT's value.
Prints PASS/FAIL against the tolerances that the frozen baseline LUT meets (--tol-*), writes OUT.json.
"""
import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import weighting_geometry as wg  # noqa: E402


def points(rng, R, PITCH, GAP, n_near=40, n_bulk=60):
    x = []
    for j in (0, 1, 2, 3, 4):
        r = R * np.exp(rng.uniform(np.log(1.02), np.log(.7 / R), n_near))
        th = rng.uniform(0, 2 * np.pi, n_near)
        x.append(np.column_stack([j * PITCH + r * np.cos(th), r * np.sin(th)]))
    u = rng.uniform(-13.4, 13.4, n_bulk); v = rng.uniform(-GAP + .02, GAP - .02, n_bulk)
    keep = np.hypot(u - PITCH * np.rint(u / PITCH), v) > .8
    x.append(np.column_stack([u[keep], v[keep]]))
    return np.vstack(x)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("geometry"); ap.add_argument("numw", type=Path); ap.add_argument("lut", type=Path)
    ap.add_argument("out", type=Path)
    ap.add_argument("--refined-L", type=float, default=96.); ap.add_argument("--refined-N", type=int, default=32)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--tol-interp", type=float, default=5e-5, help="max |LUT - direct(LUT settings)|; frozen baseline LUT: 1.5e-5")
    ap.add_argument("--tol-trunc", type=float, default=1e-6, help="max |direct(LUT) - direct(refined)|; frozen baseline: 3.6e-8")
    a = ap.parse_args()
    geom = wg.Geometry.load(a.geometry)
    mods = wg.apply(geom, a.numw)
    meta = wg.check_table(a.lut, geom)
    L = mods["p3_spectral_weighting_lut"]; S = mods["p3_fourier_strip_reconstruction_pilot"]
    M = mods["p3_fourier_mfs_mode_pilot"]
    A, P, R, PITCH, GAP = geom.stereo_tan, geom.strip_pitch_mm, M.R, M.PITCH, M.GAP
    rng = np.random.default_rng(a.seed)
    x = points(rng, R, PITCH, GAP)
    phases = np.array([-3.3, -1.1, -.6 * P, -.2, .07, .5 * P - .03, .5 * P + .03, 1.9, 4.4])
    lut = L.WeightingLUT(a.lut)
    t0 = time.perf_counter()
    direct = S.Correction(x, L.L, L.KMAX, L.N, L.IMAGE_EXTRA)
    t1 = time.perf_counter()
    refined = S.Correction(x, a.refined_L, L.KMAX, a.refined_N, L.IMAGE_EXTRA)
    t2 = time.perf_counter()
    res = {"geometry": geom.name, "geometry_id": geom.geometry_id, "lut_sha256": wg.sha256(a.lut), "points": len(x),
           "phases": phases.tolist(), "direct_s": t1 - t0, "refined_s": t2 - t1}
    pl, pd, pr = [], [], []
    for s in phases:
        w = A * x[:, 0] + s
        pl.append(lut.upper(x, w)); pd.append(direct.phi(w)); pr.append(refined.phi(w))
    pl, pd, pr = map(np.array, (pl, pd, pr))
    res["interp_max_abs"] = float(np.max(abs(pl - pd)))
    res["interp_rel_L2"] = float(np.linalg.norm(pl - pd) / np.linalg.norm(pd))
    res["trunc_max_abs"] = float(np.max(abs(pd - pr)))
    res["trunc_rel_L2"] = float(np.linalg.norm(pd - pr) / np.linalg.norm(pr))
    res["lut_min"], res["lut_max"] = float(pl.min()), float(pl.max())
    # boundaries
    th = np.linspace(.1, 2 * np.pi, 17, endpoint=False)
    xs = np.column_stack([PITCH * 2 + R * 1.0000001 * np.cos(th), R * 1.0000001 * np.sin(th)])
    res["anode_surface_max_abs"] = float(max(np.max(abs(lut.upper(xs, A * xs[:, 0] + s))) for s in phases))
    xw = np.column_stack([np.linspace(-12, 12, 25), np.full(25, GAP)])
    ind = lambda e: (abs(e) < P / 2).astype(float)
    res["strip_plane_max_abs"] = float(max(np.max(abs(lut.upper(xw, A * xw[:, 0] + s) - ind(np.full(25, s))))
                                           for s in phases if abs(abs(s) - P / 2) > 1e-6))
    # partition of unity -> whole-plane potential
    J = int(np.floor((L.PHASE_MAX - P / 2) / P))
    f0, _ = M.fit(0., L.N, 0)
    xin = x[np.abs(x[:, 0]) <= 13.4]
    tot = sum(lut.upper(xin, A * xin[:, 0] - j * P) for j in range(-J, J + 1))
    phi_plane = f0(xin).real
    res["partition_J"] = J
    res["partition_tail_max_abs"] = float(np.max(abs(phi_plane - tot)))
    res["partition_tail_min_signed"] = float(np.min(phi_plane - tot))
    res["lut_meta"] = meta
    ok = (res["interp_max_abs"] <= a.tol_interp and res["trunc_max_abs"] <= a.tol_trunc and
          res["lut_min"] >= -1e-5 and res["lut_max"] <= 1 + 1e-5 and res["anode_surface_max_abs"] <= 1e-6)
    res["tolerances"] = {"interp": a.tol_interp, "trunc": a.tol_trunc}
    res["pass"] = bool(ok)
    a.out.write_text(json.dumps(res, indent=1) + "\n")
    for k in ("interp_max_abs", "interp_rel_L2", "trunc_max_abs", "trunc_rel_L2", "lut_min", "lut_max",
              "anode_surface_max_abs", "strip_plane_max_abs", "partition_J", "partition_tail_max_abs"):
        print(f"  {k:26s} {res[k]:.6g}")
    print(("PASS" if ok else "FAIL") + f" LUT check {geom.name} ({len(x)} points x {len(phases)} phases)")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
