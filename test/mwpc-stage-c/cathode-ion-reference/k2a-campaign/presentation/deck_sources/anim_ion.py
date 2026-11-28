#!/usr/bin/env python3
"""Animated GIFs: one ion and an (illustrative) avalanche ion cloud drifting to the cathodes, with the induced charge and
current on strips. Same independent field/drift (single_ion.py) and frozen weighting LUT as the slide-2 figure.
Time runs logarithmically (each frame a fixed step in log t), because the signal grows about equally per decade."""
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.patches import Circle, Rectangle

import single_ion as SI
import style as S
from p3_spectral_weighting_lut import WeightingLUT

P = 1.7
RW = 0.015
LUT = WeightingLUT(SI.FROZEN / "lut_frozen_default.npz")
LU, LV = "#8FB0E3", "#F0B98A"     # lighter U / V for neighbour strips


def ion_signals(r_um, th, tgrid, strips=(0, 1)):
    """phi(x(t)) for U and V strips at offsets n*P; returns positions (len(t), 2) and dict[(side, n)] -> phi(t)."""
    r = r_um * 1e-3
    s = SI.drift(r * np.cos(th), r * np.sin(th))
    tt = np.minimum(tgrid, s.t[-1])
    x = s.sol(tt).T
    out = {}
    for side, nm in ((0, "U"), (1, "V")):
        xx = x if side == 0 else -x
        for n in strips:
            out[(nm, n)] = LUT.upper(xx, np.full(len(tt), -n * P))
    return x, out, s.t[-1]


def current(t, q):
    return np.gradient(q, t)


def log_frames(t0, t1, n):
    return np.geomspace(t0, t1, n)


def fmt_t(t):
    return f"{t:.3g} ns" if t < 1000 else f"{t / 1000:.3g} µs"


def build(kind):
    rng = np.random.default_rng(3)
    if kind == "one":
        ions = [(17.76, np.radians(87.5))]
    else:
        n = 60
        up = rng.random(n) < 0.8
        th = np.where(up, rng.normal(80, 22, n), rng.normal(235, 25, n))
        r = 15 + rng.lognormal(np.log(2.4), 0.4, n)
        ions = list(zip(r, np.radians(th)))
    tgrid = np.r_[0.0, np.geomspace(0.01, 1.2e5, 900)]
    X, PH, TE = [], [], []
    for r, th in ions:
        x, ph, te = ion_signals(r, th, tgrid)
        X.append(x); PH.append(ph); TE.append(te)
    tend = max(TE)
    m = tgrid <= tend * 1.02
    tg = tgrid[m]
    X = np.array([x[m] for x in X])                                   # (ions, t, 2)
    keys = list(PH[0].keys())
    Q = {k: np.mean([ph[k][m] - ph[k][0] for ph in PH], axis=0) for k in keys}   # ion term per ion [e]
    E = {k: np.mean([ph[k][0] for ph in PH]) for k in keys}                      # electron step per ion [e]
    I = {k: current(tg, Q[k]) for k in keys}                                    # e/ns
    return ions, tg, X, Q, E, I, tend


def animate(kind, fname, nframes=125, hold=18):
    ions, tg, X, Q, E, I, tend = build(kind)
    frames = log_frames(0.01, tend, nframes)
    frames = np.r_[frames, np.full(hold, tend)]
    fig = plt.figure(figsize=(13.0, 6.2), dpi=100)
    gs = fig.add_gridspec(2, 3, width_ratios=[1.35, 0.75, 1.5], height_ratios=[1, 1], wspace=0.38, hspace=0.42,
                          left=0.06, right=0.98, top=0.9, bottom=0.1)
    ax_cell = fig.add_subplot(gs[:, 0]); ax_zoom = fig.add_subplot(gs[0, 1]); ax_info = fig.add_subplot(gs[1, 1])
    ax_q = fig.add_subplot(gs[0, 2]); ax_i = fig.add_subplot(gs[1, 2])
    title = "One positive ion drifting to the cathode" if kind == "one" else "An avalanche: 60 positive ions drifting to the cathodes"
    fig.suptitle(title, x=0.06, ha="left", fontsize=15, fontweight="bold")
    # ---- cell view
    if kind == "one":
        # view along u: (w, v) plane, strips visible side by side
        ax_cell.set_xlim(-4.5, 4.5); ax_cell.set_ylim(-3.6, 3.6)
        ax_cell.axhline(0, color="#9AA0A6", lw=2.0)
        ax_cell.text(-4.35, 0.12, "anode wire (along w)", fontsize=8.5, color=S.INK2)
        strips = {}
        for side, v0, col, lcol in (("U", 3.0, S.UC, LU), ("V", -3.0, S.VC, LV)):
            for n in (-2, -1, 0, 1, 2):
                rect = Rectangle((n * P - P / 2 + 0.03, v0 - 0.12 if side == "U" else v0 - 0.0), P - 0.06, 0.12,
                                 color=col if n == 0 else lcol, alpha=0.35 if abs(n) < 2 else 0.12, lw=0)
                ax_cell.add_patch(rect)
                strips[(side, n)] = rect
                ax_cell.text(n * P, v0 + (0.12 if side == "U" else -0.14), f"{side}{n:+d}" if n else f"{side}0",
                             ha="center", va="bottom" if side == "U" else "top", fontsize=8.5, color=col)
        ax_cell.set_xlabel("w, along the wire  [mm]"); ax_cell.set_ylabel("v  [mm]")
        ax_cell.set_title("seen along the wire: the strips side by side", fontsize=10, loc="left")
        trail, = ax_cell.plot([], [], color=S.K2C, lw=1.2, alpha=0.6)
        dots, = ax_cell.plot([], [], "o", color=S.K2C, ms=7)
        cell_xy = lambda k: (np.zeros(1), X[:, k, 1])
        trail_xy = lambda k: (np.zeros(k + 1), X[0, :k + 1, 1])
    else:
        ax_cell.set_xlim(-3.6, 3.6); ax_cell.set_ylim(-3.6, 3.6)
        ax_cell.axhline(3, color=S.UC, lw=4); ax_cell.axhline(-3, color=S.VC, lw=4)
        ax_cell.text(-3.45, 3.15, "U cathode (strips)", fontsize=9, color=S.UC)
        ax_cell.text(-3.45, -3.15, "V cathode (strips)", fontsize=9, color=S.VC, va="top")
        for k in (-1, 0, 1):
            ax_cell.plot([3 * k], [0], "o", color=S.INK, ms=5)
        ax_cell.set_xlabel("u, across the wires  [mm]"); ax_cell.set_ylabel("v  [mm]")
        ax_cell.set_title("cross-section: ions fan out along the field lines", fontsize=10, loc="left")
        strips = {}
        dots, = ax_cell.plot([], [], "o", color=S.K2C, ms=3.5, alpha=0.8)
        trail = None
        cell_xy = lambda k: (X[:, k, 0], X[:, k, 1])
    ax_cell.set_aspect("equal")
    # ---- zoom near the wire (u, v)
    ax_zoom.add_patch(Circle((0, 0), RW * 1e3, color="#9AA0A6"))
    zdots, = ax_zoom.plot([], [], "o", color=S.K2C, ms=4 if kind == "one" else 2.5)
    edots, = ax_zoom.plot([], [], "o", color="#2E7D38", ms=4 if kind == "one" else 3)
    ax_zoom.set_aspect("equal"); ax_zoom.set_xlim(-80, 80); ax_zoom.set_ylim(-80, 80)
    ax_zoom.set_xlabel("u  [µm]", fontsize=9); ax_zoom.set_ylabel("v  [µm]", fontsize=9)
    ax_zoom.set_title("auto-zoom around the wire", fontsize=10, loc="left")
    ax_zoom.tick_params(labelsize=8)
    ax_info.axis("off")
    clock = ax_info.text(0.0, 0.95, "", fontsize=17, fontweight="bold", va="top", transform=ax_info.transAxes)
    info = ax_info.text(0.0, 0.62, "", fontsize=9.5, va="top", color=S.INK2, transform=ax_info.transAxes)
    # ---- charge and current
    sel = [(("U", 0), S.UC, "-", "U strip above (U0)"), (("U", 1), LU, "-", "U neighbour (U±1)"),
           (("V", 0), S.VC, "-", "V strip below (V0)"), (("V", 1), LV, "-", "V neighbour (V±1)")]
    full_q = {k: Q[k] + E[k] for k, *_ in sel}       # with the prompt electron step
    ymax = max(np.max(np.abs(full_q[k])) for k, *_ in sel)
    lq, li_pos, li_neg = {}, {}, {}
    for k, c, ls, lab in sel:
        lq[k], = ax_q.plot([], [], color=c, lw=2.0, label=lab)
        li_pos[k], = ax_i.plot([], [], color=c, lw=1.8)
        li_neg[k], = ax_i.plot([], [], color=c, lw=1.8, ls="--")
    for a in (ax_q, ax_i):
        a.set_xscale("log"); a.set_xlim(0.01, tend * 1.3)
        a.set_xticks([0.01, 0.1, 1, 10, 100, 1e3, 1e4]); a.set_xticklabels(["10 ps", "", "1 ns", "", "100 ns", "", "10 µs"])
        a.axvspan(0.01, 7000, color="#F1F2EF", zorder=0)
    ax_q.set_ylim(-0.08 * ymax, 1.08 * ymax); ax_q.axhline(0, color=S.MUTED, lw=0.7)
    ax_q.set_ylabel("induced charge  [e per ion]", fontsize=10)
    ax_q.legend(loc="upper left", fontsize=8.5)
    ax_q.text(6000, 1.02 * ymax, "0–7 µs: K2A window", fontsize=8, color=S.INK2, ha="right", va="top")
    imax = max(np.max(np.abs(I[k][tg > 0.01])) for k, *_ in sel)
    ax_i.set_yscale("log"); ax_i.set_ylim(imax * 1e-7, imax * 3)
    ax_i.set_ylabel("|current|  [e/ns per ion]", fontsize=10); ax_i.set_xlabel("time after the avalanche (log scale)")
    ax_i.plot([], [], color=S.INK2, ls="--", label="dashed: negative current")
    ax_i.legend(loc="upper right", fontsize=8.5)
    tnow_q = ax_q.axvline(0.3, color=S.K2C, lw=1.0)
    tnow_i = ax_i.axvline(0.3, color=S.K2C, lw=1.0)
    S.label(ax_i, x=0.03, y=0.13, size=9, extra=None)
    e0 = sum(E[(s, 0)] for s in "UV") / 2

    def upd(fi):
        tf = frames[fi]
        k = int(np.searchsorted(tg, tf))
        k = min(k, len(tg) - 1)
        xs, ys = cell_xy(k)
        dots.set_data(xs, ys)
        if trail is not None:
            trail.set_data(*trail_xy(k))
        zdots.set_data(X[:, k, 0] * 1e3, X[:, k, 1] * 1e3)
        hw = max(40.0, 1.25 * 1e3 * np.max(np.hypot(X[:, k, 0], X[:, k, 1])))
        ax_zoom.set_xlim(-hw, hw); ax_zoom.set_ylim(-hw, hw)
        if tf > 0.1:
            edots.set_data([], [])
        else:
            r0 = np.hypot(X[:, 0, 0], X[:, 0, 1]) * 1e3
            re = r0 - (r0 - 15.3) * (tf / 0.1)
            edots.set_data(X[:, 0, 0] * 1e3 * re / r0, X[:, 0, 1] * 1e3 * re / r0)
        for kk, *_ in sel:
            q = Q[kk][: k + 1] + E[kk]
            lq[kk].set_data(np.maximum(tg[: k + 1], 0.01), q)
            i = I[kk][: k + 1]
            tt = np.maximum(tg[: k + 1], 0.01)
            li_pos[kk].set_data(tt, np.where(i > 0, i, np.nan))
            li_neg[kk].set_data(tt, np.where(i < 0, -i, np.nan))
        for (side, n), rect in strips.items():
            key = (side, abs(n)) if (side, abs(n)) in Q else None
            if key:
                rect.set_alpha(0.15 + 0.85 * np.clip((Q[key][k] + E[key]) / max(ymax, 1e-9), 0, 1))
        tnow_q.set_xdata([tf, tf]); tnow_i.set_xdata([tf, tf])
        clock.set_text(f"t = {fmt_t(tf)}")
        vmax = np.max(np.abs(X[:, k, 1]))
        txt = (f"ion distance from the wire centre: {vmax * 1e3:.0f} µm" if vmax < 1 else f"ion distance from the wire centre: {vmax:.2f} mm") if kind == "one" \
            else f"furthest ion from the wire plane: {vmax:.2f} mm" if vmax >= 0.1 else f"furthest ion from the wire plane: {vmax * 1e3:.0f} µm"
        txt += f"\nU0 charge so far: {Q[('U', 0)][k] + E[('U', 0)]:.3f} e per ion"
        txt += "\n\ngreen dots: avalanche electrons;\nthey reach the wire surface\nwithin ≈ 0.1 ns: a prompt step\nof ≈ " + f"{100 * E[('U', 0)]:.1f} % of e on U0"
        if fi >= len(frames) - hold:
            txt += "\n\nall ions have reached the cathodes"
        info.set_text(txt)
        return []

    anim = FuncAnimation(fig, upd, frames=len(frames), blit=False)
    anim.save(fname, writer=PillowWriter(fps=10), dpi=100)
    upd(len(frames) - 1)
    fig.savefig(fname.replace(".gif", "_last.png"), dpi=100)
    plt.close(fig)


if __name__ == "__main__":
    kind = sys.argv[1]
    animate(kind, f"/home/claude/pres/fig/anim_{kind}.gif")
    print("done", kind)
