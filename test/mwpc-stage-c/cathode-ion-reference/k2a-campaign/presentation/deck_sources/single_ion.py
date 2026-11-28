#!/usr/bin/env python3
"""Independent single-ion calculation (no Garfield++) and cross-check against the frozen K2A kernel.

Field   : own implementation of the thin-wire analytic field of the frozen model (9 anodes, 30 um diameter, 3 mm
          pitch, 2100 V, grounded cathode planes at v = +-3 mm), exact 2D Green function of the parallel-plate gap.
Mobility: Garfield++ data file IonMobility_Ar+_Ar.txt (E/N [Td], K [cm2/(V s)]), scaled to 293.15 K, 760 Torr as
          MediumGas does (K * T/273.15 * 760/p); same provisional choice as the reference.
Drift   : dx/dt = mu(|E|) E, scipy RK45 (rtol 1e-9), until the ion reaches a cathode.
Weight  : frozen spectral weighting LUT (p3_spectral_weighting_lut.py, built here from the frozen code).
Check   : K0(tau) = mean over the 72 grid azimuths of phi(x(t0+tau)) - phi(x(t0)) at D = 0, compared with the frozen
          kernel slice built from Garfield++ DriftLineRKF paths (K_slice_centralD.npz).
"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

FROZEN = Path("/tmp/claude-0/-home-claude/6e19f127-ea1b-521a-80a5-970a9287006d/scratchpad/frozen")
sys.path.insert(0, str(FROZEN))
from p3_spectral_weighting_lut import WeightingLUT  # noqa: E402

UP = Path("/home/claude/pres/data/upload_20261128_143502")
OUT = Path("/home/claude/pres/an/single_ion")
OUT.mkdir(exist_ok=True)

# ---- geometry (mm) ---------------------------------------------------------------------------------------------
GAP, PITCH, RW, HV = 3.0, 3.0, 0.015, 2100.0
H = 2 * GAP
UW = np.arange(-4, 5) * PITCH
EPS0 = 8.8541878128e-12 * 1e-3          # F/mm


def g0(u, v, us, vs):
    """Dirichlet Green function of -Laplace in the gap |v| < GAP (unit line charge -> potential * eps0)."""
    dx = u - us
    y, ys = v + GAP, vs + GAP
    d = np.sinh(np.pi * dx / (2 * H)) ** 2
    return np.log((d + np.sin(np.pi * (y + ys) / (2 * H)) ** 2) / (d + np.sin(np.pi * (y - ys) / (2 * H)) ** 2)) / (4 * np.pi)


def grad_g0(u, v, us, vs):
    """Analytic gradient of g0 with respect to the field point (u, v)."""
    a = np.pi / (2 * H)
    dx = u - us
    y, ys = v + GAP, vs + GAP
    d = np.sinh(a * dx) ** 2
    dd = 2 * np.sinh(a * dx) * np.cosh(a * dx) * a
    sp, sm = np.sin(a * (y + ys)), np.sin(a * (y - ys))
    cp, cm = np.cos(a * (y + ys)), np.cos(a * (y - ys))
    np_, nm = d + sp ** 2, d + sm ** 2
    du = (dd / np_ - dd / nm) / (4 * np.pi)
    dv = ((2 * sp * cp * a) / np_ - (2 * sm * cm * a) / nm) / (4 * np.pi)
    return du, dv


# wire charges (thin-wire approximation, as an analytic-field component does): potential HV on each wire surface
M = np.array([[g0(UW[i] + (RW if i == j else 0.0), 0.0, UW[j], 0.0) for j in range(9)] for i in range(9)])
LAM = np.linalg.solve(M, np.full(9, HV))           # = q / eps0  [V]


def efield(u, v):
    eu = ev = 0.0
    for j in range(9):
        du, dv = grad_g0(u, v, UW[j], 0.0)
        eu = eu - LAM[j] * du
        ev = ev - LAM[j] * dv
    return eu, ev                                     # V/mm


def potential(u, v):
    return sum(LAM[j] * g0(u, v, UW[j], 0.0) for j in range(9))


# ---- mobility --------------------------------------------------------------------------------------------------
EN_TD = np.array([0, 8, 10, 12, 15, 20, 25, 30, 40, 50, 60, 80, 100, 120, 150, 200, 250, 300, 400, 500, 600, 800,
                  1000, 1200, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 6000, 7000, 8000, 9000], float)
K_CM2 = np.array([1.53, 1.53, 1.53, 1.53, 1.52, 1.51, 1.49, 1.47, 1.44, 1.41, 1.38, 1.32, 1.27, 1.22, 1.16, 1.06,
                  0.99, 0.95, 0.85, 0.78, 0.72, 0.63, 0.56, 0.51, 0.46, 0.40, 0.36, 0.33, 0.30, 0.28, 0.26, 0.25,
                  0.23, 0.21, 0.20, 0.19])
T_K, P_TORR = 293.15, 760.0
NDENS = P_TORR * 133.322368 / (1.380649e-23 * T_K) * 1e-6     # cm^-3
MOB_SCALE = T_K / 273.15 * 760.0 / P_TORR


def mobility_mm2_per_V_ns(emag_V_per_mm):
    en = emag_V_per_mm * 10.0 / NDENS * 1e17          # Td
    k = np.interp(en, EN_TD, K_CM2) * MOB_SCALE       # cm2/(V s)
    return k * 1e2 * 1e-9                             # mm2/(V ns)


def drift(u0, v0, tmax=2.0e5):
    def rhs(t, x):
        eu, ev = efield(x[0], x[1])
        mu = mobility_mm2_per_V_ns(np.hypot(eu, ev))
        return [mu * eu, mu * ev]

    def hit(t, x):
        return GAP - abs(x[1]) - 1e-9
    hit.terminal = True
    s = solve_ivp(rhs, (0.0, tmax), [u0, v0], method="RK45", rtol=1e-9, atol=1e-12, events=hit, dense_output=True,
                  max_step=50.0)
    return s


def main():
    lut = WeightingLUT(FROZEN / "lut_frozen_default.npz")
    ks = np.load(UP / "kernel/K_slice_centralD.npz")
    tau = ks["tau_ns"]; r_nodes = ks["r_um"]; Kref = ks["K"]          # (basis, side, wire, r, tau)
    theta = (np.arange(72) + 0.5) * 5.0
    res = {"check": {}, "notes": __doc__}
    # sanity: wire surface potential and field at the surface
    eu, ev = efield(RW, 0.0)
    res["E_surface_kV_per_cm"] = float(np.hypot(eu, ev) / 100.0)
    res["phi_surface_V"] = float(potential(RW, 0.0))
    res["phi_cathode_V"] = float(potential(0.7, GAP - 1e-12))
    for ir in (0, 1, 6):
        r = r_nodes[ir] * 1e-3
        Gsides = np.zeros((2, 72, len(tau)))
        tend = np.zeros(72)
        for it, th in enumerate(theta):
            u0, v0 = r * np.cos(np.radians(th)), r * np.sin(np.radians(th))
            s = drift(u0, v0)
            tend[it] = s.t[-1]
            tt = np.minimum(tau, s.t[-1])
            x = s.sol(tt).T
            for side in range(2):
                xx = x if side == 0 else -x
                phi = lut.upper(xx, np.zeros(len(tt)))
                Gsides[side, it] = phi - phi[0]
        K0 = Gsides.mean(axis=1)
        ref = Kref[0, :, 1, ir, :]
        out = {}
        for side, nm in enumerate("UV"):
            diff = K0[side] - ref[side]
            out[nm] = {"rel_L2": float(np.linalg.norm(diff) / np.linalg.norm(ref[side])),
                       "max_abs": float(np.max(np.abs(diff))), "ref_at_8us": float(ref[side, -1]),
                       "own_at_8us": float(K0[side, -1])}
        out["drift_time_us"] = [float(tend.min() / 1e3), float(np.median(tend) / 1e3), float(tend.max() / 1e3)]
        res["check"][f"r={r_nodes[ir]:.2f}um"] = out
        np.savez(OUT / f"k0_check_r{ir}.npz", tau=tau, own=K0, ref=ref, tend=tend)
        print(f"r={r_nodes[ir]:.2f} um", json.dumps(out), flush=True)
    (OUT / "check.json").write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
