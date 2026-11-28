#!/usr/bin/env python3
"""Ten-parameter figures.

cloud_anatomy : two ILLUSTRATIVE ion clouds of equal charge (schematic ion positions), their ten K2A moments computed
                with the exact definitions of cathode_engine.moments_from_births / BirthMoments.hh, and the strip-frame
                coordinates D_U, D_V in a top view.
param_effects : effect of each parameter on the induced strip charge (ion term, 0-7 us), computed with the independent
                field/drift of single_ion.py and the frozen weighting LUT; and the full ion-by-ion signal of a cloud vs
                its ten-number (K2A-formula) reconstruction at the same r_eff (no table interpolation involved).
"""
import json

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle, Rectangle

import single_ion as SI
import style as S
from p3_spectral_weighting_lut import WeightingLUT

RW = 15.0
P = 1.7
rng = np.random.default_rng(7)
TAU = np.r_[0.0, np.geomspace(1.0, 7000.0, 160)]


def make_clouds(n=60):
    th1 = np.radians(rng.normal(75, 13, n))
    r1 = RW + rng.lognormal(np.log(2.3), 0.35, n)
    lobe = rng.random(n) < 0.6
    th2 = np.radians(np.where(lobe, rng.normal(245, 22, n), rng.normal(140, 18, n)))
    r2 = RW + rng.lognormal(np.log(5.0), 0.55, n)
    return (r1, th1), (r2, th2)


def moments(r, th, q=None):
    q = np.ones_like(r) if q is None else q
    m = lambda x: np.sum(q * x) / np.sum(q)
    return {"r_eff": float(np.exp(m(np.log(r)))), "C1": m(np.cos(th)), "S1": m(np.sin(th)), "C2": m(np.cos(2 * th)),
            "S2": m(np.sin(2 * th)), "A": m(np.sign(np.sin(th)))}


def density(mo, th):
    return (1 + 2 * (mo["C1"] * np.cos(th) + mo["S1"] * np.sin(th) + mo["C2"] * np.cos(2 * th) + mo["S2"] * np.sin(2 * th))) / (2 * np.pi)


def anatomy(c1, c2):
    fig = plt.figure(figsize=(13.0, 5.6))
    gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 1.05], wspace=0.28, bottom=0.27)
    out = {}
    for i, ((r, th), title) in enumerate(((c1, "cloud 1: compact, upper side"), (c2, "cloud 2: wide, wrapped, mostly lower side"))):
        ax = fig.add_subplot(gs[i])
        ax.add_patch(Rectangle((-45, 0), 90, 45, color=S.UC, alpha=0.06, lw=0))
        ax.add_patch(Rectangle((-45, -45), 90, 45, color=S.VC, alpha=0.08, lw=0))
        ax.text(-42, 41, "→ U cathode", fontsize=8, color=S.UC, va="top")
        ax.text(-42, -41, "→ V cathode", fontsize=8, color=S.VC, va="bottom")
        ax.add_patch(Circle((0, 0), RW, color="#9AA0A6", zorder=2))
        ax.text(0, 0, "anode\n30 µm", ha="center", va="center", fontsize=7.5, color="white", zorder=3)
        x, y = r * np.cos(th), r * np.sin(th)
        ax.plot(x, y, "o", ms=3.2, color=S.K2C, alpha=0.75, zorder=4, label="ions (schematic)")
        mo = moments(r, th)
        out[f"cloud{i + 1}"] = mo
        ax.add_patch(Circle((0, 0), mo["r_eff"], fill=False, ls="--", lw=1.2, color=S.INK, zorder=5))
        # truncated angular density (n <= 2) drawn as a ring modulation
        tt = np.linspace(0, 2 * np.pi, 400)
        rho = density(mo, tt)
        rr = 28 + 14 * np.clip(rho, 0, None) / rho.max()
        ax.plot(rr * np.cos(tt), rr * np.sin(tt), color=S.INK2, lw=1.3, zorder=3, label="angular shape kept (n ≤ 2)")
        ax.plot(30 * np.cos(tt), 30 * np.sin(tt), color=S.GRID, lw=0.8, zorder=2)
        # first-harmonic vector
        ax.annotate("", xy=(30 * mo["C1"], 30 * mo["S1"]), xytext=(0, 0), arrowprops=dict(arrowstyle="-|>", color=S.INK, lw=1.6), zorder=6)
        # one ion: r and theta
        j = int(np.argmax(r * (np.abs(np.degrees(th) % 360 - (75 if i == 0 else 245)) < 15)))
        ax.plot([0, x[j]], [0, y[j]], color=S.INK, lw=0.8, zorder=6)
        ax.text(x[j] * 1.08 + (-12 if i == 0 else 1), y[j] * 1.08, "r, θ", fontsize=8.5, zorder=7)
        ax.set_xlim(-45, 45); ax.set_ylim(-45, 45); ax.set_aspect("equal")
        ax.set_xlabel("u − 3k  [µm]"); ax.set_ylabel("v  [µm]")
        ax.set_title(title, fontsize=10.5, loc="left")
        txt = (f"r_eff = {mo['r_eff']:.1f} µm (dashed)\n"
               f"C₁ = {mo['C1']:+.2f}, S₁ = {mo['S1']:+.2f} (arrow ×30 µm)\n"
               f"C₂ = {mo['C2']:+.2f}, S₂ = {mo['S2']:+.2f}\nA = {mo['A']:+.2f}")
        ax.text(0.0, -0.2, txt, transform=ax.transAxes, fontsize=8.5, ha="left", va="top")
        if i == 0:
            fig.legend(loc="lower right", fontsize=8.5, bbox_to_anchor=(0.98, 0.02), ncol=1, frameon=False,
                       handles=ax.get_legend_handles_labels()[0] + [plt.Line2D([], [], color=S.INK, ls="--", lw=1.2), plt.Line2D([], [], color=S.INK, lw=1.6, marker=">")],
                       labels=["ions (schematic positions, 60 per cloud)", "angular shape kept by K2A (n ≤ 2), schematic scale", "r_eff = exp⟨ln r⟩", "(C₁, S₁): mean direction"])
    # top view: strip frame
    ax = fig.add_subplot(gs[2])
    for k in (-1, 0, 1):
        ax.axvline(3 * k, color="#9AA0A6", lw=1.6)
    ax.text(0.12, -3.85, "anode wires (along w)", fontsize=8, color=S.INK2)
    uu = np.linspace(-4.5, 4.5, 2)
    for n in range(-3, 4):
        for sgn, c in ((+1, S.UC), (-1, S.VC)):
            lo = (n - 0.5) * P
            ax.plot(uu, lo + sgn * 0.1 * uu, color=c, lw=0.9, alpha=0.85)
    u0, w0 = 0.4, 0.55
    ax.plot([u0], [w0], "o", ms=8, color=S.K2C, zorder=5)
    ax.text(u0 + 0.2, w0 + 0.15, "cloud", color=S.K2C, fontsize=9)
    U0, V0 = w0 - 0.1 * u0, w0 + 0.1 * u0
    ax.annotate("", xy=(u0, w0), xytext=(u0, 0.1 * u0), arrowprops=dict(arrowstyle="<->", color=S.UC, lw=1.3))
    ax.text(u0 - 0.12, 0.3, "D_U", color=S.UC, fontsize=9, ha="right", va="center")
    ax.annotate("", xy=(u0 + 0.9, w0), xytext=(u0 + 0.9, -0.1 * (u0 + 0.9)), arrowprops=dict(arrowstyle="<->", color=S.VC, lw=1.3))
    ax.text(u0 + 1.02, 0.35, "D_V", color=S.VC, fontsize=9, va="center")
    ax.plot(uu, 0.1 * uu, color=S.UC, lw=0.8, ls=":")
    ax.plot(uu, -0.1 * uu, color=S.VC, lw=0.8, ls=":")
    ax.set_xlim(-4.5, 4.5); ax.set_ylim(-4.2, 4.2)
    ax.set_xlabel("u, across the wires  [mm]"); ax.set_ylabel("w, along the wires  [mm]")
    ax.set_title("top view: strip boundaries, U (blue) and V (orange)", fontsize=10.5, loc="left")
    ax.text(-4.3, 3.95, "U = w − 0.1u, V = w + 0.1u\npitch 1.7 mm in w, tilt ±5.7°\n(dotted: strip centre lines)", fontsize=8, va="top",
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec=S.GRID))
    S.save(fig, "cloud_anatomy")
    return out


class Paths:
    def __init__(self):
        self.lut = WeightingLUT(SI.FROZEN / "lut_frozen_default.npz")
        self.cache = {}

    def phi(self, r_um, th, D):
        """phi_s(x(t0 + TAU), D) for one ion start point; returns (2 sides, len(D), len(TAU))."""
        key = (round(r_um, 6), round(th, 9))
        if key not in self.cache:
            r = r_um * 1e-3
            s = SI.drift(r * np.cos(th), r * np.sin(th))
            self.cache[key] = s.sol(np.minimum(TAU, s.t[-1])).T
        x = self.cache[key]
        out = np.empty((2, len(D), len(TAU)))
        for side in range(2):
            xx = x if side == 0 else -x
            for j, d in enumerate(D):
                out[side, j] = self.lut.upper(xx, np.full(len(TAU), d))
        return out

    def cloud(self, r, th, D, w=None):
        w = np.ones_like(r) if w is None else w
        acc = 0
        for ri, ti, wi in zip(r, th, w):
            p = self.phi(ri, ti, D)
            acc = acc + wi * (p - p[:, :, :1])
        return acc / np.sum(w)


def k2a_reconstruction(P_, mo, D, M=72):
    th = (np.arange(M) + 0.5) * 2 * np.pi / M
    sm = np.sign(np.sin(th))
    sp = sm - sm.mean()
    for n in (1, 2):
        for b in (np.cos(n * th), np.sin(n * th)):
            sp = sp - b * (sp @ b) / (b @ b)
    a = (1 + 2 * sum(mo[f"C{n}"] * np.cos(n * th) + mo[f"S{n}"] * np.sin(n * th) for n in (1, 2))) / M
    beta = (mo["A"] - a @ sm) / (sp @ sm)
    a = a + beta * sp
    return P_.cloud(np.full(M, mo["r_eff"]), th, D, a * M), float(beta)


def effects(c1, c2, mo):
    P_ = Paths()
    D = np.array([0.0, P, 0.85])
    tl = TAU[1:] / 1e3
    fig, axs = plt.subplots(1, 4, figsize=(15.0, 4.2), gridspec_kw={"wspace": 0.32})
    r1, th1 = c1
    # (a) radius: cloud-1 shape, ion offsets above the wire rescaled so r_eff = 16, 17.5, 20, 25 um
    ax = axs[0]
    off = r1 - RW
    curves = {}
    for reff in (16.0, 17.5, 20.0, 25.0):
        lo, hi = 1e-3, 50.0
        for _ in range(60):
            s_ = np.sqrt(lo * hi)
            re = np.exp(np.mean(np.log(RW + s_ * off)))
            lo, hi = (s_, hi) if re < reff else (lo, s_)
        curves[reff] = P_.cloud(RW + s_ * off, th1, D[:1])[0, 0, 1:]
    for reff, ls in ((16.0, ":"), (17.5, "-"), (20.0, "--"), (25.0, "-.")):
        ax.plot(tl, curves[reff] / curves[17.5], color=S.UC, ls=ls, lw=1.6, label=f"r_eff = {reff:g} µm")
    ax.set_ylim(0.3, 1.35)
    ax.text(0.96, 0.96, "same angles; ratio to\nr_eff = 17.5 µm", transform=ax.transAxes, ha="right", va="top", fontsize=8.5, color=S.INK2)
    ax.set_title("r_eff: matters early, not late", fontsize=10.5, loc="left")
    ax.set_ylabel("U strip above, relative", fontsize=10.5)
    # (b) angle / hemisphere: cloud 1 (upper), mirror (lower), uniform ring
    ax = axs[1]
    ring = ((np.arange(72) + 0.5) * 2 * np.pi / 72)
    for (rr, tt, lab, ls) in ((r1, th1, "upper side (A = +1)", "-"), (r1, -th1, "lower side (A = −1)", "--"),
                              (np.full(72, np.exp(np.mean(np.log(r1)))), ring, "uniform ring (A = 0)", ":")):
        g = P_.cloud(rr, tt, D[:1])
        ax.plot(tl, g[0, 0, 1:], color=S.UC, ls=ls, lw=1.6, label=lab)
    ax.text(0.04, 0.66, "U strip above the wire;\nthe V strip is the mirror image", transform=ax.transAxes, fontsize=8.5, color=S.INK2, va="top")
    ax.set_ylabel("induced charge per ion  [e]", fontsize=10.5)
    ax.set_title("angles and A: which plane", fontsize=10.5, loc="left")
    # (c) D: strip offset
    ax = axs[2]
    g = P_.cloud(r1, th1, D)
    for j, (lab, ls) in enumerate((("D = 0 (strip centre)", "-"), ("D = 1.7 mm (next strip)", "--"), ("D = 0.85 mm (strip edge)", ":"))):
        ax.plot(tl, g[0, j, 1:], color=S.UC, ls=ls, lw=1.6, label=lab)
    ax.set_title("D_U, D_V: which strips share", fontsize=10.5, loc="left")
    # (d) cloud 2: ion by ion vs ten numbers (same r_eff node: compression only)
    ax = axs[3]
    r2, th2 = c2
    full = P_.cloud(r2, th2, D[:1])
    rec, beta = k2a_reconstruction(P_, mo["cloud2"], D[:1])
    for side, c in ((0, S.UC), (1, S.VC)):
        ax.plot(tl, full[side, 0, 1:], color=S.S0C, lw=1.8, label="ion by ion (60 ions)" if side == 0 else None)
        ax.plot(tl, rec[side, 0, 1:], color=K2C if False else S.K2C, lw=1.4, ls="--", label="ten numbers (K2A formula)" if side == 0 else None)
    ax.text(0.04, 0.72, "upper pair: V strip below\nlower pair: U strip above", transform=ax.transAxes, fontsize=8.5, color=S.INK2, va="top")
    ax.set_title("cloud 2: all ions vs ten numbers", fontsize=10.5, loc="left")
    e = {}
    for side, nm in ((0, "U"), (1, "V")):
        f, rr_ = full[side, 0, 1:], rec[side, 0, 1:]
        e[nm] = {"rel_L2": float(np.linalg.norm(rr_ - f) / np.linalg.norm(f)), "at_7us_rel": float(rr_[-1] / f[-1] - 1)}
    e["beta"] = beta
    for ax in axs:
        ax.set_xscale("log")
        ax.legend(loc="lower center" if ax is axs[0] else "upper left", fontsize=7.8, ncol=2 if ax is axs[0] else 1); ax.set_xlim(1e-3, 7.0)
        ax.set_xticks([1e-3, 1e-2, 1e-1, 1]); ax.set_xticklabels(["1 ns", "10 ns", "100 ns", "1 µs"])
        ax.set_xlabel("time after birth (0–7 µs)")
    S.label(axs[3], x=0.97, y=0.15, ha="right", size=9, extra=[f"difference: U {100 * e['U']['rel_L2']:.2f} %, V {100 * e['V']['rel_L2']:.2f} % (L2)"])
    S.save(fig, "param_effects")
    return e


if __name__ == "__main__":
    c1, c2 = make_clouds()
    mo = anatomy(c1, c2)
    e = effects(c1, c2, mo)
    (SI.OUT / "clouds.json").write_text(json.dumps({"moments": mo, "cloud2_reconstruction": e}, indent=1))
    print(json.dumps({"moments": mo, "cloud2_reconstruction": e}, indent=1))
