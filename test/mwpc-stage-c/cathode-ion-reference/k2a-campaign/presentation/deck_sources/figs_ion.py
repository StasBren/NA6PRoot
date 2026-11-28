#!/usr/bin/env python3
"""One-ion figure and independent cross-check figure (no Garfield++; see single_ion.py)."""
import json

import matplotlib.pyplot as plt
import numpy as np

import single_ion as SI
import style as S
from p3_spectral_weighting_lut import WeightingLUT

P = 1.7
OUT = SI.OUT


def one_ion():
    lut = WeightingLUT(SI.FROZEN / "lut_frozen_default.npz")
    r0, th = 0.01776, 87.5
    u0, v0 = r0 * np.cos(np.radians(th)), r0 * np.sin(np.radians(th))
    s = SI.drift(u0, v0)
    tend = s.t[-1]
    t = np.r_[0.0, np.geomspace(0.01, tend, 1500)]
    t[-1] = tend
    x = s.sol(t).T
    strips = {"U": [0, 1, 2], "V": [0, 1]}
    Q = {}
    for side, nm in enumerate("UV"):
        xx = x if side == 0 else -x
        for n in range(-7, 8):
            phi = lut.upper(xx, np.full(len(t), n * P))
            Q[(nm, n)] = phi                      # phi(x(t)); ion term = phi - phi[0], electron step = phi[0]
    info = {"tend_us": tend / 1e3, "r0_um": r0 * 1e3, "theta_deg": th}
    for nm in "UV":
        plane0 = sum(Q[(nm, n)][0] for n in range(-7, 8))
        plane_end = sum(Q[(nm, n)][-1] for n in range(-7, 8))
        info[f"{nm}_plane_phi0"] = float(plane0)
        info[f"{nm}_plane_phi_end"] = float(plane_end)
        info[f"{nm}0_phi0"] = float(Q[(nm, 0)][0])
        info[f"{nm}0_end"] = float(Q[(nm, 0)][-1])
        i7 = np.searchsorted(t, 7000.0)
        info[f"{nm}0_ion_at_7us"] = float(Q[(nm, 0)][i7] - Q[(nm, 0)][0])
        info[f"{nm}1_ion_max"] = float(np.max(Q[(nm, 1)] - Q[(nm, 1)][0]))
        info[f"{nm}1_end"] = float(Q[(nm, 1)][-1] - Q[(nm, 1)][0])
    # radius vs time marks
    rr = np.hypot(x[:, 0], x[:, 1])
    marks = {}
    for tm in (1.0, 10.0, 100.0, 1000.0, 7000.0):
        marks[tm] = float(np.interp(tm, t, x[:, 1]))
    info["v_at_mm"] = {str(k): v for k, v in marks.items()}
    (OUT / "one_ion.json").write_text(json.dumps(info, indent=1))

    fig = plt.figure(figsize=(12.6, 4.6))
    gs = fig.add_gridspec(1, 3, width_ratios=[0.82, 1.25, 1.25], wspace=0.34)
    a0 = fig.add_subplot(gs[0]); a1 = fig.add_subplot(gs[1]); a2 = fig.add_subplot(gs[2])
    # (a) cell cross-section with the ion path and field lines
    for th2 in np.arange(2.5, 360, 15.0):
        s2 = SI.drift(r0 * np.cos(np.radians(th2)), r0 * np.sin(np.radians(th2)))
        a0.plot(s2.y[0], s2.y[1], color="#C9CED4", lw=0.7, zorder=1)
    a0.plot(x[:, 0], x[:, 1], color=S.K2C, lw=2.2, zorder=3)
    for tm, vv in marks.items():
        if tm < 1000:
            continue
        a0.plot([0], [vv], "o", ms=4, color=S.K2C, zorder=4)
        a0.text(0.12, vv, f"{tm / 1000:g} µs", fontsize=8.5, va="center", color=S.INK2)
    a0.text(0.12, 2.75, f"arrives at {tend / 1e3:.0f} µs", fontsize=8.5, va="center", color=S.K2C)
    a0.text(-1.5, -1.05, f"after 100 ns the ion is\nonly {1e3 * marks[100.0]:.0f} µm from the wire axis", fontsize=8, color=S.INK2, va="top")
    for k in (-1, 0, 1):
        a0.plot([3 * k], [0], "o", ms=5, color=S.INK, zorder=5)
    a0.axhline(3, color=S.UC, lw=3); a0.axhline(-3, color=S.VC, lw=3)
    a0.text(-1.5, 3.1, "U cathode (strips)", color=S.UC, fontsize=9, va="bottom")
    a0.text(-1.5, -3.1, "V cathode (strips)", color=S.VC, fontsize=9, va="top")
    a0.text(0.12, -0.25, "anode", fontsize=8, color=S.INK, va="top")
    a0.set_xlim(-1.6, 1.6); a0.set_ylim(-3.7, 3.7)
    a0.set_xlabel("u, across the wires  [mm]"); a0.set_ylabel("v, towards the cathodes  [mm]")
    a0.set_aspect("auto")
    # (b) induced charge on strips, ion term, log time
    tl = t[1:] / 1e3
    def ion(nm, n):
        return (Q[(nm, n)] - Q[(nm, n)][0])[1:]
    a1.plot(tl, ion("U", 0), color=S.UC, lw=2.0, label="U strip above the ion")
    a1.plot(tl, ion("U", 1), color=S.UC, lw=1.4, ls="--", label="U neighbour (±1.7 mm)")
    a1.plot(tl, ion("V", 0), color=S.VC, lw=2.0, label="V strip below")
    a1.plot(tl, ion("V", 1), color=S.VC, lw=1.4, ls="--", label="V neighbour")
    a1.axvspan(1e-3, 7.0, color="#F1F2EF", zorder=0)
    a1.text(6.0, 0.97, "0–7 µs:\nK2A window", fontsize=8.5, color=S.INK2, va="top", ha="right")
    a1.axhline(0, color=S.MUTED, lw=0.7)
    a1.set_xscale("log"); a1.set_xlim(1e-3, tend / 1e3 * 1.3)
    a1.set_xlabel("time after ion birth"); a1.set_ylabel("induced charge, ion term  [e]")
    a1.legend(loc="upper left", fontsize=8.5, bbox_to_anchor=(0.0, 0.99))
    S.label(a1, x=0.04, y=0.66, size=9.5, extra=["one Ar⁺ ion, 2100 V", f"born {r0 * 1e3:.1f} µm from the wire axis", "(wire radius 15 µm)"])
    a1.set_xticks([1e-3, 1e-2, 1e-1, 1, 1e1]); a1.set_xticklabels(["1 ns", "10 ns", "100 ns", "1 µs", "10 µs"])
    # (c) planes: ion + electron
    for nm, c in (("U", S.UC), ("V", S.VC)):
        pl = sum(Q[(nm, n)] for n in range(-7, 8))
        a2.plot(tl, (pl - pl[0])[1:], color=c, lw=2.0, label=f"{nm} plane, ion term")
        a2.plot(tl, pl[1:], color=c, lw=1.2, ls=":", label=f"{nm} plane, ion + electrons")
    a2.axvspan(1e-3, 7.0, color="#F1F2EF", zorder=0)
    a2.axhline(0, color=S.MUTED, lw=0.7); a2.axhline(1, color=S.MUTED, lw=0.7, ls="--")
    a2.set_xscale("log"); a2.set_xlim(1e-3, tend / 1e3 * 1.3)
    a2.set_xticks([1e-3, 1e-2, 1e-1, 1, 1e1]); a2.set_xticklabels(["1 ns", "10 ns", "100 ns", "1 µs", "10 µs"])
    a2.set_xlabel("time after ion birth"); a2.set_ylabel("induced charge, whole plane  [e]")
    a2.legend(loc="upper left", fontsize=8.5)
    a2.text(0.04, 0.62, "electrons: a prompt step of\n" + f"{100 * info['U_plane_phi0']:.1f} % (U), {100 * info['V_plane_phi0']:.1f} % (V) of e\n" + "ion + electrons → 1 e on U, 0 on V", transform=a2.transAxes, fontsize=8.5, color=S.INK2, va="top")
    S.save(fig, "one_ion")
    return info


def crosscheck():
    fig, axs = plt.subplots(2, 1, figsize=(6.2, 4.6), sharex=True, gridspec_kw={"height_ratios": [2.2, 1], "hspace": 0.08})
    res = json.loads((OUT / "check.json").read_text())
    for ir, ls in ((1, "-"), (6, "--")):
        d = np.load(OUT / f"k0_check_r{ir}.npz")
        tau = d["tau"][1:] / 1e3
        for side, c in ((0, S.UC),):
            ref, own = d["ref"][side, 1:], d["own"][side, 1:]
            r_um = float(json.load(open(SI.UP / "kernel/green_axes.json"))["r_um"][ir])
            axs[0].plot(tau, ref, color=S.S0C, lw=1.8, ls=ls, label=f"Garfield++ table, r₀ = {r_um:.0f} µm")
            axs[0].plot(tau, own, color=S.K2C, lw=1.4, ls=":", label=f"independent Python, r₀ = {r_um:.0f} µm")
            m = np.abs(ref) > 0.02 * np.abs(ref).max()
            axs[1].plot(tau[m], 100 * (own[m] - ref[m]) / ref[m], color=S.K2C, lw=1.2, ls=ls)
    axs[0].set_xscale("log"); axs[0].set_ylabel("induced charge, strip above wire  [e]", fontsize=10.5)
    axs[0].legend(fontsize=8, loc="upper left", frameon=False)
    S.label(axs[0], x=0.97, y=0.35, ha="right", size=9.5, extra=["ion term, angle-averaged", "U strip at D = 0"])
    axs[1].axhline(0, color=S.MUTED, lw=0.7)
    axs[1].set_ylabel("diff.  [%]"); axs[1].set_xlabel("time after ion birth  [µs]")
    axs[1].set_xlim(1e-3, 8.0)
    axs[1].set_xticks([1e-3, 1e-2, 1e-1, 1]); axs[1].set_xticklabels(["1 ns", "10 ns", "100 ns", "1 µs"])
    axs[1].set_xlabel("time after ion birth")
    S.save(fig, "crosscheck")
    return res


if __name__ == "__main__":
    print(json.dumps(one_ion(), indent=1))
    crosscheck()
