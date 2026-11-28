#!/usr/bin/env python3
"""Single-ion response table for the K2A explorer app.

G[r, theta, D, tau] = phi_U(x(t0 + tau), D) - phi_U(x(t0), D)   (ion term per unit charge, U side, frozen sign convention)
PHI0[r, theta, D]   = phi_U(x(t0), D)                            (prompt electron step)
computed with the independent field/drift of single_ion.py and the frozen spectral weighting LUT.
Symmetries used by the app (checked here numerically):
  G_U(theta, D) = G_U(180 - theta, -D)      (u -> -u, w -> -w leaves the cell and the stereo phase w - A u unchanged)
  G_V(theta, D) = G_U(theta + 180, D)       (the frozen convention: side V is the U table at x -> -x)
so only theta in (-90, 90) degrees is stored.
"""
import base64
import json
import sys

import numpy as np

sys.path.insert(0, "/home/claude/pres/an")
import single_ion as SI  # noqa: E402
from p3_spectral_weighting_lut import WeightingLUT  # noqa: E402

LUT = WeightingLUT(SI.FROZEN / "lut_frozen_default.npz")
R = np.geomspace(15.2, 70.0, 10)                       # um
TH = np.r_[np.arange(-87.5, 90, 5.0)]                  # deg, 36 nodes
D = np.round(np.arange(-25, 26) * 0.17, 6)             # mm, 51 nodes
TAU = np.r_[0.0, np.geomspace(1.0, 7000.0, 47)]        # ns


def phi_path(r_um, th_deg, side=0, d=D, tau=TAU):
    r = r_um * 1e-3
    s = SI.drift(r * np.cos(np.radians(th_deg)), r * np.sin(np.radians(th_deg)))
    x = s.sol(np.minimum(tau, s.t[-1])).T
    xx = x if side == 0 else -x
    out = np.empty((len(d), len(tau)))
    for j, dd in enumerate(d):
        out[j] = LUT.upper(xx, np.full(len(tau), dd))
    return out


def main():
    G = np.zeros((len(R), len(TH), len(D), len(TAU)), np.float64)
    P0 = np.zeros((len(R), len(TH), len(D)))
    for i, r in enumerate(R):
        for j, th in enumerate(TH):
            ph = phi_path(r, th)
            G[i, j] = ph - ph[:, :1]
            P0[i, j] = ph[:, 0]
        print("r", r, flush=True)
    # symmetry checks on a few points
    chk = {}
    for r, th in ((R[1], 37.5), (R[4], -62.5)):
        a = phi_path(r, th); b = phi_path(r, 180 - th, d=-D)
        chk[f"mirror r={r:.1f} th={th}"] = float(np.max(np.abs(a - b)))
        c = phi_path(r, th, side=1); e = phi_path(r, th + 180, side=0)
        chk[f"V=U(th+180) r={r:.1f} th={th}"] = float(np.max(np.abs(c - e)))
    print(chk)
    # int16 quantisation per (r, theta, D) row
    sc = np.maximum(np.max(np.abs(G), axis=3), 1e-12)
    q = np.round(G / sc[..., None] * 32767).astype(np.int16)
    err = np.max(np.abs(q / 32767 * sc[..., None] - G))
    data = {"R_um": R.tolist(), "TH_deg": TH.tolist(), "D_mm": D.tolist(), "TAU_ns": TAU.tolist(),
            "shape": list(G.shape), "G_i16_b64": base64.b64encode(q.tobytes()).decode(),
            "G_scale": sc.astype(np.float32).ravel().tolist(), "PHI0": np.round(P0, 7).ravel().tolist(),
            "quant_max_abs_err": float(err), "symmetry_checks": chk,
            "provenance": "independent field/drift (single_ion.py, Ar+ in Ar mobility table) + frozen spectral LUT sha256 ed458136e3e2…"}
    open("/home/claude/pres/app/table.json", "w").write(json.dumps(data))
    print("quant err", err, "size MB", len(data["G_i16_b64"]) / 1e6)


if __name__ == "__main__":
    main()
