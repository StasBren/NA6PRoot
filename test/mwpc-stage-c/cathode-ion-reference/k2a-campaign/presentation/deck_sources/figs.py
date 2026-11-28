#!/usr/bin/env python3
"""All data figures of the K2A deck. Each function documents its source files."""
import glob
import json
import re
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, FancyArrowPatch, Rectangle, Wedge

import style as S

UP = "/home/claude/pres/data/upload_20261128_143502"
A = pd.read_pickle("/home/claude/pres/an/acc.pkl")
P = pd.read_pickle("/home/claude/pres/an/pop.pkl")
ST = pd.read_pickle("/home/claude/pres/an/stress.pkl")
EX = f"{UP}/examples"


def aligned(m, variant):
    """S0T1 labels; K2A (K2A_full) columns reordered onto them (zero where a strip is absent)."""
    l0, t0, c0 = S.load_variant(f"{EX}/muon_{m}/S0T1")
    lk, tk, ck = S.load_variant(f"{EX}/muon_{m}/{variant}")
    K = np.column_stack([ck[:, lk.index(x)] if x in lk else np.zeros(len(ck)) for x in l0])
    nb = min(len(c0), len(K))
    return l0, t0[:nb], c0[:nb], K[:nb]


# ---------------------------------------------------------------- slide 2: headline
def fig_headline():
    fig, ax = plt.subplots(figsize=(6.1, 4.5))
    b = np.linspace(0, 1.1, 111)
    for p, c, ls in (("U", S.UC, "-"), ("V", S.VC, "--")):
        x = 100 * A[f"sh100:{p}_L2_rel"].values
        ax.hist(x, bins=b, histtype="step", color=c, ls=ls, lw=1.8, label=f"{p} plane: {S.ms(x, ' %')}")
    ax.axvline(1.0, color=S.K2C, lw=1.4, ls=":")
    ax.text(0.975, ax.get_ylim()[1] * 0.30, "frozen\ncriterion v1:\nL2 ≤ 1 %", ha="right", va="center", fontsize=9.5, color=S.K2C)
    ax.set_xlabel("Shaped relative L2, K2A vs S0_T1  [%]")
    ax.set_ylabel("Events / 0.01 %")
    ax.set_xlim(0, 1.1)
    S.label(ax, x=0.36, y=0.95, extra=["192 independent muons, 2100 V, B = 0", "CR-RC³ stand-in, Tp = 100 ns, 0–6 µs"])
    ax.legend(loc="upper right", bbox_to_anchor=(1.0, 0.70), fontsize=9.5)
    ax.set_ylim(0, ax.get_ylim()[1] * 1.05)
    S.save(fig, "headline_L2")


# ---------------------------------------------------------------- slide 3: signal formation
def fig_cell():
    """Schematic of the reference cell (drawn from the frozen constants: pitch 3 mm, gap 3 mm, 30 µm wires)."""
    fig, ax = plt.subplots(figsize=(5.6, 4.4))
    ax.plot(3 * np.arange(-4, 5), np.zeros(9), "o", ms=7, color=S.INK, zorder=5)
    for v in (3, -3):
        ax.plot([-13.5, 13.5], [v, v], color=S.INK, lw=2.4)
        for x0 in np.arange(-13.5, 13.5, 1.7):
            ax.plot([x0 + 0.1, x0 + 1.6], [v + (0.22 if v > 0 else -0.22)] * 2, color=S.UC if v > 0 else S.VC, lw=3, solid_capstyle="butt")
    ax.text(-13.5, 3.75, "U strips, pitch 1.7 mm  (U = w − 0.1 u)", fontsize=9.5, color=S.UC)
    ax.text(-13.5, -4.35, "V strips, pitch 1.7 mm  (V = w + 0.1 u)", fontsize=9.5, color=S.VC)
    for (x0, y0, x1, y1, rad) in ((0.2, 0.25, 1.6, 2.8, 0.25), (-0.15, -0.25, -1.3, -2.8, 0.25)):
        ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=13, color=S.K2C, lw=1.8,
                                     connectionstyle=f"arc3,rad={rad}"))
    ax.text(2.2, 1.55, "ions drift to the cathodes\n≈ 40–100 µs (illustration)", color=S.K2C, fontsize=9.5, ha="left", va="center")
    ax.annotate("", xy=(6, -0.9), xytext=(3, -0.9), arrowprops=dict(arrowstyle="<->", color=S.INK2, lw=1))
    ax.text(4.5, -1.15, "3 mm", ha="center", va="top", fontsize=9.5, color=S.INK2)
    ax.annotate("", xy=(-11.2, 3), xytext=(-11.2, 0), arrowprops=dict(arrowstyle="<->", color=S.INK2, lw=1))
    ax.text(-10.9, 1.5, "3 mm", ha="left", va="center", fontsize=9.5, color=S.INK2)
    ax.text(-3.0, -0.75, "anodes (30 µm)", ha="center", fontsize=9, color=S.INK2, va="top")
    ax.set_xlim(-14, 14); ax.set_ylim(-4.8, 4.4)
    ax.set_xlabel("u, across the wires  [mm]"); ax.set_ylabel("v  [mm]")
    ax.set_yticks([-3, 0, 3])
    ax.minorticks_off(); ax.tick_params(top=False, right=False)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    S.save(fig, "cell_schematic")


def fig_signal():
    """Example 0219 (median-L2 event): S0_T1 full-length replay (100 µs), ion term, collecting strip and a neighbour."""
    m = "0219"
    lab, t, ion = S.load_variant(f"{EX}/muon_{m}/S0T1_full", prompt=False)
    _, _, pr = S.TS.load_wave(f"{EX}/muon_{m}/S0T1_full/prompt/prompt_waveforms.csv")
    q = ion.sum(0) * 10
    U = [j for j, l in enumerate(lab) if l.startswith("U_")]
    lead = U[int(np.argmax(np.abs(q[U])))]
    nb_ = U[int(np.argsort(np.abs(q[U]))[-2])]
    nice = lambda l: l.replace("_", " ").replace("-", "−")
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(6.4, 4.4), gridspec_kw=dict(wspace=0.42))
    n = 100
    a1.plot(t[:n], ion[:n, lead], color=S.S0C, lw=1.6, label=f"collecting {nice(lab[lead])}")
    a1.plot(t[:n], ion[:n, nb_], color=S.UC, lw=1.4, ls="--", label=f"neighbour {nice(lab[nb_])}")
    a1.axhline(0, color=S.GRID, lw=0.8)
    a1.set_xlabel("Time  [ns]"); a1.set_ylabel("Ion current  [fC/ns]")
    a1.set_xlim(0, 1000); a1.set_ylim(-0.02, 0.22)
    a1.legend(loc="center right", fontsize=8.5)
    a1.text(0.97, 0.62, "fast current while the ions\nleave the anode region", transform=a1.transAxes, ha="right", fontsize=8.5, color=S.INK2)
    for idx, c, ls, sc in ((lead, S.S0C, "-", 1), (nb_, S.UC, "--", 1)):
        a2.plot(t / 1000, np.cumsum(ion[:, idx]) * 10 * sc, color=c, ls=ls, lw=1.6)
    a2.axvspan(0.005, 7.0, color="#F3E1E4", lw=0)
    a2.axhline(0, color=S.GRID, lw=0.8)
    frac7 = np.cumsum(ion[:, lead])[699] / ion[:, lead].sum()
    a2.text(0.0065, 150, f"0–7 µs window\n(compared with S0):\n{100 * frac7:.0f} % of the final\ncollecting-strip charge", fontsize=8.5, color=S.K2C, va="top")
    a2.text(9, 45, "neighbour: bipolar,\nnet ≈ 0 at 100 µs", fontsize=8.5, color=S.UC, ha="right")
    a2.set_xscale("log"); a2.set_xlim(0.005, 100); a2.set_ylim(-10, 160)
    a2.set_xlabel("Time  [µs]"); a2.set_ylabel("Cumulative induced charge  [fC]")
    S.label(a1, x=0.04, y=0.96, size=9.5, extra=[f"S0_T1, muon {m}", "ion term"])
    S.save(fig, "signal_anatomy")
    return {"slot": m, "lead": lab[lead], "neighbour": lab[nb_], "frac_at_7us": float(frac7), "lead_charge_fC": float(q[lead]),
            "neighbour_charge_fC": float(q[nb_]), "neighbour_prompt_fC": float(pr[:, nb_].sum() * 10)}


# ---------------------------------------------------------------- slide 4: S0 cost
def fig_s0cost():
    comp = [("RKF ion drift\n(+ archive read)", "rtS0:rkf_export_archive_read_plus_ion_drift"),
            ("weighting-field replay", "rtS0:weighting_replay_T0_chosen_attempt_sum"),
            ("electron step", "rtS0:electron_step_T0"),
            ("S0_T0 total", "rtS0:total_T0_single_pass")]
    fig, ax = plt.subplots(figsize=(6.1, 4.5))
    data = [A[c].clip(lower=1e-2).values for _, c in comp]
    bp = ax.boxplot(data, vert=False, whis=(5, 95), showfliers=True, widths=0.5, patch_artist=True,
                    flierprops=dict(marker="o", ms=3, mfc="none", mec=S.MUTED),
                    medianprops=dict(color=S.K2C, lw=2), boxprops=dict(facecolor="#EEF0EC", edgecolor=S.INK),
                    whiskerprops=dict(color=S.INK), capprops=dict(color=S.INK))
    ax.set_yticks(range(1, len(comp) + 1), [n for n, _ in comp])
    ax.set_xscale("log"); ax.set_xlim(0.1, 3e4); ax.set_ylim(4.6, 0.4)
    ax.set_xlabel("Wall time per event  [s]   (one CPU process)")
    for i, d in enumerate(data):
        ax.text(0.13, i + 1, f"{np.median(d):.3g} s", va="center", ha="left", fontsize=10, color=S.K2C)
    ax.invert_yaxis()
    S.label(ax, x=0.97, y=0.97, ha="right", extra=["192 events, median in red", "box 25–75 %, whiskers 5–95 %"])
    ax.tick_params(axis="y", which="minor", left=False, right=False)
    S.save(fig, "s0_cost")


# ---------------------------------------------------------------- slide 5: moments
def acc_moments():
    rows = []
    for s in A.slot:
        f = glob.glob(f"{UP}/preflight/moments_nonsealed/*/muon_{s}.moments.csv")[0]
        m = pd.read_csv(f); m["slot"] = s; rows.append(m)
    return pd.concat(rows)


def fig_moments():
    M = acc_moments()
    fig = plt.figure(figsize=(6.6, 4.5))
    gs = fig.add_gridspec(2, 2, hspace=0.5, wspace=0.32)
    ax0 = fig.add_subplot(gs[:, 0])
    # definition sketch around one anode
    ax0.add_patch(Circle((0, 0), 0.15, color=S.INK, zorder=4))
    ax0.add_patch(Wedge((0, 0), 1.6, 0, 180, color="#F3E1E4", zorder=0))
    ax0.add_patch(Wedge((0, 0), 1.6, 180, 360, color="#E3E9F3", zorder=0))
    ax0.add_patch(Circle((0, 0), 0.9, fill=False, ls="--", color=S.INK2, lw=1.2))
    ax0.plot([0, 0.9 * np.cos(0.9)], [0, 0.9 * np.sin(0.9)], color=S.INK, lw=1.2)
    ax0.text(0.22, 0.55, r"$r_{\rm eff}$", fontsize=11)
    th = np.linspace(0, 0.9, 30); ax0.plot(0.42 * np.cos(th), 0.42 * np.sin(th), color=S.K2C, lw=1.2)
    ax0.text(0.5, 0.12, "θ", color=S.K2C, fontsize=12)
    ax0.plot([-1.7, 1.7], [0, 0], color=S.MUTED, lw=0.8)
    ax0.text(-1.55, 1.35, "v > 0: +1", fontsize=9.5, color=S.K2C)
    ax0.text(-1.55, -1.5, "v < 0: −1", fontsize=9.5, color=S.UC)
    ax0.text(0, -2.05, "A = ⟨sign sin θ⟩ (charge-weighted)", ha="center", fontsize=9, color=S.INK2)
    ax0.set_xlim(-1.8, 1.8); ax0.set_ylim(-2.25, 1.8); ax0.set_aspect("equal"); ax0.axis("off")
    q = M.N_fC.values
    a1 = fig.add_subplot(gs[0, 1]); a2 = fig.add_subplot(gs[1, 1])
    a1.hist(M.r_eff_um, bins=np.linspace(16.5, 19, 51), weights=None, histtype="step", color=S.INK, lw=1.4)
    a1.set_xlabel(r"$r_{\rm eff}=\exp\langle\ln r\rangle$  [µm]"); a1.set_ylabel("Clouds")
    S.label(a1, x=0.03, y=0.94, size=8.5)
    a1.text(0.97, 0.94, f"{len(M)} clouds\nmean {M.r_eff_um.mean():.2f} µm\nRMS {M.r_eff_um.std(ddof=0):.2f} µm",
            transform=a1.transAxes, ha="right", va="top", fontsize=8.5)
    a1.set_ylim(0, a1.get_ylim()[1] * 1.7)
    a2.hist(M.A, bins=np.linspace(-1, 1, 41), histtype="step", color=S.INK, lw=1.4, weights=q, density=True)
    a2.set_xlabel("Hemisphere asymmetry A"); a2.set_ylabel("Charge-weighted\ndensity")
    a2.text(0.5, 0.93, f"mean {np.average(M.A, weights=q):+.2f}, RMS {np.sqrt(np.average((M.A - np.average(M.A, weights=q)) ** 2, weights=q)):.2f}",
            transform=a2.transAxes, ha="center", va="top", fontsize=8.5)
    a2.set_ylim(0, a2.get_ylim()[1] * 1.3)
    S.label(a2, x=0.5, y=1.25, size=9) if False else None
    S.save(fig, "moments")
    h1 = np.hypot(M.C1, M.S1)
    return {"n_clouds": int(len(M)), "reff_mean": float(M.r_eff_um.mean()), "reff_rms": float(M.r_eff_um.std(ddof=0)),
            "reff_min": float(M.r_eff_um.min()), "reff_max": float(M.r_eff_um.max()),
            "A_mean_q": float(np.average(M.A, weights=q)), "frac_absA_gt_0.9": float(np.average(np.abs(M.A) > 0.9, weights=q)),
            "h1_mean": float(h1.mean())}


# ---------------------------------------------------------------- slide 6: approximation budget (S-ladder closure)
def closure_table(path):
    rows = []
    for l in open(path):
        if l.startswith("| muon_") or l.startswith("| g_"):
            c = [x.strip() for x in l.strip().strip("|").split("|")]
            rows.append({"muon": c[0], "pair": c[1], "plane": c[2], "L2": float(c[3].rstrip(" %")),
                         "s0_10": float(c[4].rstrip(" %")), "s10_20": float(c[5].rstrip(" %")), "s20_50": float(c[6].rstrip(" %")),
                         "s50_100": float(c[7].rstrip(" %")), "s100": float(c[8].rstrip(" %")), "shift": float(c[9]),
                         "L2_shift": float(c[10].rstrip(" %"))})
    return pd.DataFrame(rows)


def fig_closure():
    C = closure_table(f"{UP}/probes/closure_sladder/closure_anatomy.md")
    order = [("G_vs_S0", "Green-table nodes +\ninterpolation (G vs S0)"), ("SW_vs_G", "one D per cloud\n(SW vs G)"),
             ("S3A_vs_G", "harmonics n ≤ 2 + hemisphere\n(S3A vs G)"), ("S3_vs_G", "harmonics n ≤ 2 only\n(S3 vs G)"),
             ("SR_vs_G", "one radius r_eff\n(SR vs G)"), ("K2A_vs_S0", "all combined\n(K2A vs S0)")]
    fig, ax = plt.subplots(figsize=(6.1, 4.5))
    y = np.arange(len(order))
    res = {}
    for i, (p, lab) in enumerate(order):
        v = C[C.pair == p].L2.values
        res[p] = {"mean": float(v.mean()), "min": float(v.min()), "max": float(v.max()), "n": int(len(v))}
        col = S.K2C if p == "K2A_vs_S0" else S.INK
        ax.barh(i, v.mean(), color=col if p == "K2A_vs_S0" else "#D5D8DC", edgecolor=col, height=0.55)
        ax.plot(v, np.full_like(v, i) + np.linspace(-0.15, 0.15, len(v)), "o", ms=3.2, mfc="white", mec=col)
        ax.text(max(v.max(), v.mean()) + 0.03, i, f"{v.mean():.3f} %", va="center", fontsize=9.5, color=S.INK2)
    ax.set_yticks(y, [l for _, l in order], fontsize=9.5)
    ax.invert_yaxis(); ax.set_xlim(0, 0.85)
    ax.set_xlabel("Unshaped relative L2, ion term, 0–7 µs  [%]")
    ax.tick_params(axis="y", which="both", left=False, right=False)
    S.label(ax, x=0.97, y=0.97, ha="right", size=10, extra=["5 reference muons × 2 planes;", "bar = mean, points = muon-planes"])
    S.save(fig, "closure_budget")
    return res


# ---------------------------------------------------------------- slide 9: waveforms
def fig_waveforms(m, name, tp=100.0, tmax=1500):
    lab, t, c0, ck = aligned(m, "K2A_full")
    Y0 = S.shaped(c0, 10.0, tp); YK = S.shaped(ck, 10.0, tp)
    ts = np.arange(len(Y0)) * 1.0
    fig, axs = plt.subplots(2, 2, figsize=(11.6, 4.5), sharex=True, gridspec_kw=dict(height_ratios=[2.3, 1], hspace=0.07, wspace=0.17))
    out = {}
    for j, side in enumerate("UV"):
        cols = [i for i, l in enumerate(lab) if l.startswith(side + "_")]
        sgn = 1.0 if np.abs(Y0[:, cols].max()) >= np.abs(Y0[:, cols].min()) else -1.0
        pk = (sgn * Y0[:, cols]).max(0)
        top = [cols[i] for i in np.argsort(pk)[::-1][:3]]
        lead = float((sgn * Y0[:, top[0]]).max())
        a, r = axs[0, j], axs[1, j]
        for k, cidx in enumerate(top):
            lw = 1.8 if k == 0 else 1.2
            a.plot(ts, Y0[:, cidx], color=S.S0C, lw=lw, alpha=1 if k == 0 else 0.75)
            a.plot(ts, YK[:, cidx], color=S.K2C, lw=lw, ls="--", alpha=1 if k == 0 else 0.85)
            r.plot(ts, 100 * (YK[:, cidx] - Y0[:, cidx]) / lead, color=[S.INK, S.UC, S.VC][k], lw=1.1,
                   label=lab[cidx].replace("_", " ").replace("-", "−") + (" (lead)" if k == 0 else ""))
        a.axhline(0, color=S.GRID, lw=0.8); r.axhline(0, color=S.GRID, lw=0.8)
        l2 = np.linalg.norm(YK[:, cols] - Y0[:, cols]) / np.linalg.norm(Y0[:, cols])
        out[side] = float(l2)
        a.set_ylabel("Shaped signal  [fC-eq.]" if j == 0 else "")
        r.set_ylabel("Δ / lead peak [%]" if j == 0 else "", fontsize=10)
        r.set_xlabel("Time  [ns]"); a.set_xlim(0, tmax)
        a.text(0.97, 0.93, f"{side} plane · plane L2 = {100 * l2:.3f} %", transform=a.transAxes, ha="right", va="top", fontsize=10)
        r.legend(loc="upper right", ncol=3, fontsize=8.5, handlelength=1.2)
        lim = max(0.25, 1.25 * np.abs(100 * (YK[:, top] - Y0[:, top]) / lead).max())
        r.set_ylim(-lim, lim)
    S.label(axs[0, 1], x=0.97, y=0.80, ha="right", extra=[f"muon {m}, CR-RC³ Tp = {tp:g} ns", "three largest strips per plane"])
    axs[0, 0].legend([Line2D([], [], color=S.S0C, lw=1.8), Line2D([], [], color=S.K2C, lw=1.8, ls="--")],
                     ["S0_T1 (ion + electron step)", "K2A"], loc="center right", fontsize=9.5)
    S.save(fig, name)
    return out


# ---------------------------------------------------------------- slide 10: distributions + complexity
def fig_dists():
    fig, axs = plt.subplots(2, 2, figsize=(6.4, 4.6), gridspec_kw=dict(hspace=0.55, wspace=0.28))
    panels = [("sh100:{p}_L2_rel", "Shaped L2, Tp = 100 ns  [%]", np.linspace(0.1, 0.32, 45)),
              ("sh200:{p}_L2_rel", "Shaped L2, Tp = 200 ns  [%]", np.linspace(0.05, 0.4, 45)),
              ("m:ion_plus_electron_step.{p}.unshaped_L2_rel", "Unshaped current L2  [%]", np.linspace(0.3, 0.85, 45)),
              ("m:ion_plus_electron_step.{p}.window_charge_rel_diff", "Window charge ΔQ/Q  [%]", np.linspace(-0.65, 0.65, 53))]
    for ax, (key, xl, b) in zip(axs.flat, panels):
        for p, c, ls in (("U", S.UC, "-"), ("V", S.VC, "--")):
            x = 100 * A[key.format(p=p)].values
            ax.hist(x, bins=b, histtype="step", color=c, ls=ls, lw=1.5, label=f"{p}: {S.ms(x, '', '{:.3f}')}")
        ax.set_xlabel(xl, fontsize=10); ax.set_ylabel("Events", fontsize=10)
        ax.legend(loc="upper left", fontsize=7.6, handlelength=1.4, borderaxespad=0.2)
        ax.set_ylim(0, ax.get_ylim()[1] * 1.6)
        ax.tick_params(labelsize=9)
    fig.text(0.98, 0.955, "192 events, K2A vs S0_T1, ion + electron step", ha="right", va="bottom", fontsize=9, color=S.INK2)
    fig.text(0.02, 0.955, "NA60+", ha="left", va="bottom", fontsize=10, fontweight="bold", fontstyle="italic", family="DejaVu Sans")
    fig.text(0.115, 0.955, "Simulation", ha="left", va="bottom", fontsize=10)
    fig.subplots_adjust(top=0.93)
    S.save(fig, "dists")


def fig_complexity():
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(5.0, 4.6), gridspec_kw=dict(hspace=0.55))
    y = 100 * A["sh100:U_L2_rel"].values; yv = 100 * A["sh100:V_L2_rel"].values
    x = A["k2a_clouds"].values
    a1.plot(x, y, "o", ms=3.2, mfc="none", mec=S.UC, alpha=0.7, label="U")
    a1.plot(x + 0.25, yv, "s", ms=3.0, mfc="none", mec=S.VC, alpha=0.7, label="V")
    bins = [6, 12, 15, 18, 21, 27]
    for lo, hi in zip(bins[:-1], bins[1:]):
        sel = (x >= lo) & (x < hi)
        if sel.sum() > 2:
            v = np.r_[y[sel], yv[sel]]
            a1.errorbar((lo + hi - 1) / 2, v.mean(), yerr=v.std(ddof=0), xerr=(hi - lo) / 2, color=S.INK, marker="D", ms=4, capsize=0, lw=1.4)
    a1.set_xlabel("Clouds (primary × wire) per event"); a1.set_ylabel("Shaped L2  [%]")
    a1.set_ylim(0.1, 0.36)
    a1.legend(loc="upper right", ncol=2, fontsize=8.5, handlelength=1)
    a1.text(0.02, 0.93, "◆ mean ± RMS of U and V", transform=a1.transAxes, fontsize=8.5, va="top", color=S.INK2)
    # station / channel
    for k, (ch, mk) in enumerate((("jpsi", "o"), ("phi", "s"), ("omega", "^"))):
        for st in range(6):
            sel = (A.station_i == st) & (A.channel_i == ch)
            v = np.r_[100 * A.loc[sel, "sh100:U_L2_rel"], 100 * A.loc[sel, "sh100:V_L2_rel"]]
            if len(v):
                a2.errorbar(st + (k - 1) * 0.2, v.mean(), yerr=v.std(ddof=0), marker=mk, ms=5, mfc="white" if k else S.INK,
                            mec=S.INK, color=S.INK, lw=1.0, capsize=0, label=("J/ψ" if ch == "jpsi" else "φ" if ch == "phi" else "ω") if st == 0 else None)
    a2.set_xticks(range(6)); a2.set_xlabel("Muon station"); a2.set_ylabel("Shaped L2  [%]")
    a2.set_ylim(0.1, 0.36); a2.set_xlim(-0.5, 5.5)
    fig.text(0.02, 0.955, "NA60+", ha="left", va="bottom", fontsize=10, fontweight="bold", fontstyle="italic", family="DejaVu Sans")
    fig.text(0.155, 0.955, "Simulation · shaped L2, Tp = 100 ns", ha="left", va="bottom", fontsize=10)
    a2.tick_params(axis="x", which="minor", bottom=False, top=False)
    a2.legend(loc="upper right", ncol=3, fontsize=8.5, handlelength=1, columnspacing=0.8)
    a2.text(0.02, 0.93, "mean ± RMS (U and V)", transform=a2.transAxes, fontsize=8.5, va="top", color=S.INK2)
    fig.subplots_adjust(top=0.93)
    S.save(fig, "complexity")
    from scipy.stats import spearmanr
    return {"spearman_clouds_U": float(spearmanr(x, y).statistic), "spearman_thetau_U": float(spearmanr(A.abs_theta_u_pop if "abs_theta_u_pop" in A else A.theta_u_deg.abs(), y).statistic)}


# ---------------------------------------------------------------- slide 11: observables
def fig_observables():
    fig, axs = plt.subplots(2, 2, figsize=(6.9, 4.6), gridspec_kw=dict(hspace=0.55, wspace=0.3))
    a = axs[0, 0]; x = 1e3 * A["sh100:u_mm"].values
    a.hist(x, bins=np.linspace(-1.6, 1.6, 41), histtype="step", color=S.INK, lw=1.5)
    a.set_xlabel("Δu centroid  [µm]", fontsize=10); a.set_ylabel("Events", fontsize=10)
    a.text(0.03, 0.95, f"{S.ms(x, ' µm', '{:.2f}')}\ncriterion |Δu| ≤ 50 µm", transform=a.transAxes, va="top", fontsize=7.8)
    a.set_ylim(0, a.get_ylim()[1] * 1.55)
    a = axs[0, 1]; x = 1e3 * A["sh100:w_mm"].values
    a.hist(x, bins=np.linspace(-0.16, 0.16, 41), histtype="step", color=S.INK, lw=1.5)
    a.set_xlabel("Δw centroid  [µm]", fontsize=10); a.set_ylabel("Events", fontsize=10)
    a.text(0.03, 0.95, f"{S.ms(x, ' µm', '{:.3f}')}\ncriterion |Δw| ≤ 10 µm", transform=a.transAxes, va="top", fontsize=7.8)
    a.set_ylim(0, a.get_ylim()[1] * 1.55)
    a = axs[1, 0]
    for p, c, ls in (("U", S.UC, "-"), ("V", S.VC, "--")):
        x = A[f"sh100:{p}_t_f0.1_ns"].values
        a.hist(x, bins=np.linspace(0.06, 0.14, 33), histtype="step", color=c, ls=ls, lw=1.5, label=f"{p}: {S.ms(x, ' ns', '{:.3f}')}")
    a.set_xlabel("Δt at 10 % threshold  [ns]", fontsize=10); a.set_ylabel("Events", fontsize=10)
    a.legend(loc="upper left", fontsize=7.6, handlelength=1.3); a.set_ylim(0, a.get_ylim()[1] * 1.65)
    a.text(0.97, 0.66, "criterion |Δt| ≤ 1 ns", transform=a.transAxes, ha="right", fontsize=7.8)
    a = axs[1, 1]
    for p, c, ls in (("U", S.UC, "-"), ("V", S.VC, "--")):
        x = 100 * A[f"sh100:{p}_lead_peak_rel"].values
        a.hist(x, bins=np.linspace(-0.24, -0.08, 33), histtype="step", color=c, ls=ls, lw=1.5, label=f"{p}: {S.ms(x, ' %', '{:.3f}')}")
    a.set_xlabel("Lead-strip peak, (K2A−S0)/S0  [%]", fontsize=10); a.set_ylabel("Events", fontsize=10)
    a.legend(loc="upper left", fontsize=7.6, handlelength=1.3); a.set_ylim(0, a.get_ylim()[1] * 1.65)
    for ax in axs.flat:
        ax.tick_params(labelsize=9)
    fig.text(0.02, 0.955, "NA60+", ha="left", va="bottom", fontsize=10, fontweight="bold", fontstyle="italic", family="DejaVu Sans")
    fig.text(0.115, 0.955, "Simulation", ha="left", va="bottom", fontsize=10)
    fig.text(0.98, 0.955, "192 events, K2A − S0_T1, Tp = 100 ns", ha="right", va="bottom", fontsize=9, color=S.INK2)
    fig.subplots_adjust(top=0.93)
    S.save(fig, "observables")


# ---------------------------------------------------------------- slide 12: runtime
def fig_runtime():
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11.6, 4.4), gridspec_kw=dict(wspace=0.28, width_ratios=[1, 1.05]))
    comp = [("S0_T0\ntotal", "rtS0:total_T0_single_pass", S.S0C), ("K2A\narchive read", "rtK2A:archive_io_read_validate_min", S.K2C),
            ("K2A\nmoments", "rtK2A:moments_construction_min", S.K2C), ("K2A\nevaluation", "rtK2A:response_eval_warm_median", S.K2C),
            ("kernel open\n(once/job)", "rtK2A:kernel_open_verify_once", S.MUTED)]
    data = [A[c].values for _, c, _ in comp]
    bp = a1.boxplot(data, whis=(5, 95), widths=0.5, patch_artist=True, showfliers=True,
                    flierprops=dict(marker="o", ms=3, mfc="none", mec=S.MUTED), medianprops=dict(color="white", lw=2))
    for b, (_, _, c) in zip(bp["boxes"], comp):
        b.set_facecolor(c); b.set_edgecolor(c)
    for k in ("whiskers", "caps"):
        for i, l in enumerate(bp[k]):
            l.set_color(comp[i // 2][2])
    a1.set_yscale("log"); a1.set_ylim(1e-4, 3e4)
    a1.set_xticks(range(1, 6), [n for n, _, _ in comp], fontsize=9)
    a1.tick_params(axis="x", which="both", bottom=False, top=False)
    a1.set_ylabel("Wall time per event  [s]")
    for i, d in enumerate(data):
        a1.text(i + 1, np.quantile(d, 0.95) * 3, f"{np.median(d):.3g} s", ha="center", fontsize=8.5, color=S.INK2)
    S.label(a1, x=0.97, y=0.95, ha="right", extra=["192 events, median shown"])
    x = A["k2a_clouds"].values
    for key, lab, mk, c in (("sp:S0_T0/cathode_response_eval", "response evaluation only", "o", S.K2C),
                            ("sp:S0_T0/response_from_births_incl_moments", "+ moment construction", "s", S.INK2),
                            ("sp:S0_T0/end_to_end_from_archive_incl_io", "+ archive I/O (end to end)", "^", S.INK)):
        v = A[key].values
        a2.plot(x, v, mk, ms=3.4, mfc="none", mec=c, alpha=0.8, label=f"{lab}: median {np.median(v):,.0f}×")
    a2.set_yscale("log"); a2.set_ylim(1e3, 3e6)
    a2.set_xlabel("Clouds per event"); a2.set_ylabel("Speedup  S0_T0 / K2A")
    a2.legend(loc="upper left", fontsize=9)
    a2.text(0.97, 0.04, "avalanche generation excluded on both sides\nisolated high points: event 0262 (12,000 s RKF wall time)",
            transform=a2.transAxes, ha="right", fontsize=8.5, color=S.INK2)
    S.save(fig, "runtime")


# ---------------------------------------------------------------- slide 13: selection bias + domain
def fig_bias():
    pp = P[P.m_clouds.notna()].copy()
    acc = pp.slot.isin(A.slot)
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.0, 3.6), gridspec_kw=dict(wspace=0.32))
    b = np.logspace(np.log10(1000), np.log10(1e5), 36)
    a1.hist(pp.ion_rows_saved, bins=b, histtype="stepfilled", color="#E4E6E9", edgecolor=S.MUTED, label=f"all non-sealed ({len(pp)})")
    a1.hist(pp.ion_rows_saved[acc], bins=b, histtype="step", color=S.K2C, lw=1.8, label=f"evaluated ({acc.sum()})")
    for _, r in ST.iterrows():
        a1.axvline(r.ion_rows_saved, color=S.INK, ls=":", lw=1.2)
    a1.text(ST.ion_rows_saved.min() * 0.93, a1.get_ylim()[1] * 0.72, "stress\nevents", fontsize=8, color=S.INK, ha="right")
    a1.set_xscale("log"); a1.set_xlabel("Archived ion samples per event (K256 rows)", fontsize=10); a1.set_ylabel("Events")
    a1.legend(loc="upper left", fontsize=8.5)
    a1.set_ylim(0, a1.get_ylim()[1] * 1.35)
    cb = np.logspace(np.log10(5), np.log10(800), 36)
    a2.hist(pp.m_clouds, bins=cb, histtype="stepfilled", color="#E4E6E9", edgecolor=S.MUTED, density=False)
    a2.hist(pp.m_clouds[acc], bins=cb, histtype="step", color=S.K2C, lw=1.8)
    a2.set_ylim(0, a2.get_ylim()[1] * 1.3)
    a2.set_xscale("log"); a2.set_xlabel("Clouds per event", fontsize=10); a2.set_ylabel("Events")
    frac = float((pp.m_clouds <= A.k2a_clouds.max()).mean())
    a2.text(0.97, 0.95, f"evaluated: ≤ {int(A.k2a_clouds.max())} clouds\n({100 * frac:.0f} % of the\npopulation lies there)",
            transform=a2.transAxes, ha="right", va="top", fontsize=8.5)
    for ax in (a1, a2):
        ax.tick_params(labelsize=9)
    fig.text(0.02, 1.0, "NA60+", ha="left", va="bottom", fontsize=10, fontweight="bold", fontstyle="italic", family="DejaVu Sans")
    fig.text(0.11, 1.0, "Simulation · non-sealed 2100 V preflight population", ha="left", va="bottom", fontsize=10)
    S.save(fig, "selection_bias")
    return {"frac_pop_clouds_le_max": frac, "acc_clouds_max": int(A.k2a_clouds.max()), "pop_clouds_median": float(pp.m_clouds.median())}


# ---------------------------------------------------------------- appendix
def fig_kernel():
    z = np.load(f"{UP}/kernel/K_slice_centralD.npz")
    K, tau, r = z["K"], z["tau_ns"], z["r_um"]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11.6, 4.2), gridspec_kw=dict(wspace=0.25))
    cm = plt.get_cmap("Greys")
    sel = [0, 3, 6, 9, 13]
    for k, i in enumerate(sel):
        a1.plot(tau[1:], K[0, 0, 1, i, 1:], color=cm(0.35 + 0.6 * k / (len(sel) - 1)), lw=1.6, label=f"r = {r[i]:.1f} µm")
    a1.set_xscale("log"); a1.set_xlabel("τ since birth  [ns]"); a1.set_ylabel("K₀: induced charge per unit ion charge")
    a1.legend(loc="center left", fontsize=9)
    S.label(a1, x=0.04, y=0.95, extra=["U side, wire 0, D = 0 mm (strip centre)"])
    for k, tt in enumerate((10, 100, 1000, 7000)):
        it = int(np.argmin(np.abs(tau - tt)))
        a2.plot(r, K[0, 0, 1, :, it] / K[0, 0, 1, 0, it], "o-", ms=3.5, color=[S.INK, S.UC, S.VC, S.K2C][k], lw=1.3, label=f"τ = {tau[it]:,.0f} ns" if tau[it] >= 100 else f"τ = {tau[it]:.1f} ns")
    a2.set_xscale("log"); a2.set_xlabel("Birth radius r  [µm]"); a2.set_ylabel("K₀(r) / K₀(r = 15 µm)")
    from matplotlib.ticker import FixedLocator, NullFormatter, ScalarFormatter
    a2.xaxis.set_major_locator(FixedLocator([15, 20, 30, 50, 100])); a2.xaxis.set_major_formatter(ScalarFormatter())
    a2.xaxis.set_minor_formatter(NullFormatter())
    a2.legend(loc="lower left", fontsize=9)
    a2.text(0.97, 0.95, "early response depends strongly on r:\nthe single-radius closure matters most\nin the first tens of ns",
            transform=a2.transAxes, ha="right", va="top", fontsize=9, color=S.INK2)
    S.save(fig, "kernel")


def fig_horizon():
    fig, ax = plt.subplots(figsize=(6.1, 4.2))
    res = {}
    for m, c, ls in (("0219", S.INK, "-"), ("0581", S.UC, "--"), ("0637", S.VC, "-.")):
        lab, t, ion = S.load_variant(f"{EX}/muon_{m}/S0T1_full", prompt=False)
        for side in ("U",):
            cols = [j for j, l in enumerate(lab) if l.startswith(side + "_")]
            q = ion[:, cols].sum(0)
            lead = cols[int(np.argmax(np.abs(q)))]
            cq = np.cumsum(ion[:, lead]) / ion[:, lead].sum()
            ax.plot(t / 1000, cq, color=c, ls=ls, lw=1.6, label=f"muon {m}, strip {lab[lead].replace('_', ' ').replace('-', '−')}")
            res[m] = float(cq[699])
    ax.axvline(7, color=S.K2C, ls=":", lw=1.3); ax.axvline(6, color=S.MUTED, ls=":", lw=1.0)
    ax.text(7.6, 0.5, "7 µs:\nK2A window", color=S.K2C, fontsize=9)
    ax.text(5.6, 0.5, "6 µs: shaped\ncomparison", color=S.INK2, fontsize=8.5, ha="right", va="bottom")
    ax.set_xscale("log"); ax.set_xlim(0.01, 100); ax.set_ylim(0, 1.05)
    ax.set_xlabel("Time  [µs]"); ax.set_ylabel("Collecting-strip charge / final charge")
    ax.legend(loc="upper left", fontsize=9, bbox_to_anchor=(0, 0.88)); S.label(ax, x=0.03, y=0.97, extra=None)
    ax.text(0.97, 0.03, "S0_T1, ion term, 100 µs replay", transform=ax.transAxes, ha="right", fontsize=9, color=S.INK2)
    S.save(fig, "horizon")
    return res


def fig_early():
    C = closure_table(f"{UP}/probes/closure_0637/closure_anatomy.md")
    fig, ax = plt.subplots(figsize=(6.1, 4.2))
    pairs = [("G_vs_S0", "G vs S0"), ("SR_vs_G", "SR vs G"), ("S3A_vs_G", "S3A vs G"), ("K2A_vs_S0", "K2A vs S0")]
    w = ["s0_10", "s10_20", "s20_50", "s50_100", "s100"]; wl = ["0–10", "10–20", "20–50", "50–100", ">100 ns"]
    cols = ["#5A6470", S.K2C, "#E8A0AB", "#C9CED4", "#EEF0EC"]
    y = np.arange(len(pairs))
    res = {}
    for i, (p, lab) in enumerate(pairs):
        d = C[C.pair == p]
        L2 = d.L2.mean(); sh = d[w].mean().values
        left = 0
        for k, s in enumerate(sh):
            ax.barh(i, L2 * s / 100, left=left, color=cols[k], edgecolor="white", height=0.55, label=wl[k] if i == 0 else None)
            left += L2 * s / 100
        ax.text(L2 + 0.02, i, f"{L2:.2f} %", va="center", fontsize=9.5)
        res[p] = {"L2": float(L2), "share_10_20": float(sh[1])}
    ax.set_yticks(y, [l for _, l in pairs]); ax.invert_yaxis()
    ax.set_xlabel("Unshaped L2, ion term, muon 0637 (mean over groups)  [%]"); ax.set_xlim(0, 1.45)
    ax.tick_params(axis="y", which="both", left=False, right=False)
    ax.legend(title="share of L2² by time window", loc="upper center", bbox_to_anchor=(0.5, -0.2), fontsize=8.5, title_fontsize=8.5, ncol=5,
              handlelength=1.2, columnspacing=0.8)
    S.label(ax, x=0.97, y=0.97, ha="right", size=10)
    S.save(fig, "early_0637")
    return res


def fig_station_table():
    pass


if __name__ == "__main__" and len(sys.argv) == 1:
    info = {}
    fig_headline(); fig_cell()
    info["signal"] = fig_signal(); fig_s0cost(); info["moments"] = fig_moments(); info["closure"] = fig_closure()
    info["wave_typical"] = fig_waveforms("0219", "waves_typical")
    info["wave_worst"] = fig_waveforms("0581", "waves_worst")
    fig_dists(); info["complexity"] = fig_complexity(); fig_observables(); fig_runtime(); info["bias"] = fig_bias()
    fig_kernel(); info["horizon"] = fig_horizon(); info["early"] = fig_early()
    json.dump(info, open("/home/claude/pres/an/fig_info.json", "w"), indent=1)
    print(json.dumps(info, indent=1))


def criteria_table():
    rows = []
    for s in A.slot:
        o = json.load(open(f"{UP}/main/muons/muon_{s}/outcome.json"))
        for c in o["criteria_v1_dry_run"]["primary_vs_S0T1"]["checks"]:
            rows.append({"slot": s, "crit": c["criterion"], "key": c["key"], "Tp": c["Tp_ns"], "value": c["value"], "limit": c["limit"],
                         "status": c["status"], "note": c.get("note", "")})
    return pd.DataFrame(rows)


def fig_margin():
    C = criteria_table()
    cont = [("shaped_L2_ceiling", "shaped L2 ≤ 1 %"), ("window_charge", "window charge |ΔQ|/Q ≤ 0.5 %"),
            ("UV_charge_ratio", "U/V charge ratio |Δ| ≤ 1 %"), ("threshold_time_f0.1", "time at 10 % threshold |Δt| ≤ 1 ns"),
            ("second_over_lead", "second/lead |Δ| ≤ 0.005"), ("u_centroid", "u centroid |Δu| ≤ 50 µm"), ("w_centroid", "w centroid |Δw| ≤ 10 µm")]
    fig, ax = plt.subplots(figsize=(7.2, 4.5))
    out = {}
    for i, (k, lab) in enumerate(cont):
        d = C[C.crit == k].copy(); d["r"] = d.value.abs() / d.limit
        per = d.groupby("slot").r.max().values
        worst = d.loc[d.r.idxmax()]
        out[k] = {"max_ratio": float(per.max()), "median_ratio": float(np.median(per)), "worst_slot": worst.slot, "worst_key": worst.key,
                  "worst_Tp": float(worst.Tp), "worst_value": float(worst.value), "n_above_0.8": int((per > 0.8).sum()), "n_above_0.5": int((per > 0.5).sum())}
        y = i + np.random.default_rng(i).uniform(-0.22, 0.22, len(per))
        ax.plot(per, y, "o", ms=2.6, mfc="none", mec=S.K2C if per.max() > 0.8 else S.INK2, alpha=0.6)
        ax.plot([np.median(per)] * 2, [i - 0.3, i + 0.3], color=S.INK, lw=2)
        ax.text(1.03, i, f"max {100 * per.max():.0f} %", va="center", fontsize=9.5, color=S.K2C if per.max() > 0.8 else S.INK2)
    ax.axvline(1.0, color=S.K2C, lw=1.4, ls=":")
    ax.set_yticks(range(len(cont)), [l for _, l in cont], fontsize=9.5)
    ax.invert_yaxis(); ax.set_xlim(0, 1.22)
    ax.set_xlabel("Fraction of the frozen v1 limit used  (per event: worst plane and Tp)")
    ax.tick_params(axis="y", which="both", left=False, right=False)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0], ["0", "25 %", "50 %", "75 %", "100 %"])
    S.label(ax, x=0.40, y=0.30, size=9.5, extra=None)
    ax.text(0.40, 0.235, "192 events · bar = median\nred: criteria above 80 % of the limit", transform=ax.transAxes, fontsize=8.5, color=S.INK2, va="top")
    S.save(fig, "criteria_margin")
    disc = {}
    for k in ("lead_strip", "hit_flips"):
        d = C[C.crit == k]
        disc[k] = {"decisions": int(len(d)), "nonzero": int((d.value != 0).sum()),
                   "exempt": d.loc[d.value != 0, ["slot", "key", "Tp", "status", "note"]].values.tolist()}
    out["discrete"] = disc
    return out


if __name__ == "__main__" and len(sys.argv) > 1 and sys.argv[1] == "margin":
    r = fig_margin()
    json.dump(r, open("/home/claude/pres/an/margin.json", "w"), indent=1, default=str)
    print(json.dumps(r, indent=1, default=str))
