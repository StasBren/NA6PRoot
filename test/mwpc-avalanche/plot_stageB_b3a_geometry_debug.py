#!/usr/bin/env python3
"""
Static geometry-debug plot for Stage B3a.

This plot is deliberately independent of Garfield output. It visualises the
pure geometry of one tilted strip family in the (u,w) readout plane and makes
explicit where strips k=-1,0,+1 cross one selected anode wire.

Coordinate convention:
    u : across anode wires, in the detector plane
    w : along anode wires, in the detector plane
    v : normal to the cathode planes

For a +alpha strip family:
    x_alpha = w cos(alpha) - u sin(alpha)

Strip k is centred where
    x_alpha = k p

so on a wire at u=u_wire:
    w_k = u_wire tan(alpha) + k p / cos(alpha).

This is useful for understanding why a point that looks "under strip -1" in
ordinary w can still couple most strongly to strip 0: the relevant coordinate
is x_alpha, not w alone.
"""

import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def strip_center_w(u_mm, tan_alpha, pitch_mm, k):
    alpha = np.arctan(tan_alpha)
    return u_mm * tan_alpha + k * pitch_mm / np.cos(alpha)


def projected_x(w_mm, u_mm, tan_alpha):
    alpha = np.arctan(tan_alpha)
    return w_mm * np.cos(alpha) - u_mm * np.sin(alpha)


def main():
    p = argparse.ArgumentParser(
        description="Plot a static Stage-B3a tilted-strip geometry diagnostic."
    )
    p.add_argument("--tan-alpha", type=float, default=0.10)
    p.add_argument("--wire-u-mm", type=float, default=4.0)
    p.add_argument("--w0-mm", type=float, default=0.52)
    p.add_argument("--strip-pitch-mm", type=float, default=1.7)
    p.add_argument("--wire-pitch-mm", type=float, default=4.0)
    p.add_argument(
        "--output",
        default="stageB_b3a_geometry_debug.png",
    )
    args = p.parse_args()

    if args.strip_pitch_mm <= 0 or args.wire_pitch_mm <= 0:
        raise ValueError("Pitch values must be positive.")

    alpha = np.arctan(args.tan_alpha)
    alpha_deg = np.degrees(alpha)

    u_wire = args.wire_u_mm
    w0 = args.w0_mm
    pstrip = args.strip_pitch_mm

    xalpha = projected_x(w0, u_wire, args.tan_alpha)
    eta = xalpha / pstrip

    focus_ks = np.array([-1, 0, 1], dtype=int)
    focus_w = np.array(
        [
            strip_center_w(u_wire, args.tan_alpha, pstrip, k)
            for k in focus_ks
        ],
        dtype=float,
    )

    nearest_index = int(
        focus_ks[np.argmin(np.abs(focus_w - w0))]
    )

    wire_indices = np.arange(-2, 3)
    wire_positions = wire_indices * args.wire_pitch_mm

    umin = float(wire_positions.min() - 1.2)
    umax = float(wire_positions.max() + 1.2)

    w_pad = 0.45
    wmin = min(float(focus_w.min()), w0) - w_pad
    wmax = max(float(focus_w.max()), w0) + w_pad

    fig, (ax_top, ax_1d) = plt.subplots(
        1,
        2,
        figsize=(13.5, 5.8),
        gridspec_kw={"width_ratios": [1.12, 0.88]},
    )

    fig.suptitle(
        rf"Stage B3a geometry: why the relevant coordinate is $x_\alpha$, not $w$ alone"
        "\n"
        rf"$\tan\alpha={args.tan_alpha:g}$ "
        rf"($\alpha={alpha_deg:.2f}^\circ$), "
        rf"$u_{{\rm wire}}={u_wire:g}$ mm, "
        rf"$p={pstrip:g}$ mm",
        fontsize=15,
    )

    # ------------------------------------------------------------------
    # Left: top view in (u,w).
    # ------------------------------------------------------------------
    u_grid = np.linspace(umin, umax, 500)

    for m, u in zip(wire_indices, wire_positions):
        selected = np.isclose(u, u_wire)
        ax_top.plot(
            [u, u],
            [wmin, wmax],
            color="tab:blue" if selected else "0.65",
            linewidth=2.8 if selected else 1.0,
            alpha=0.95 if selected else 0.45,
            zorder=2,
        )
        ax_top.text(
            u,
            wmax - 0.04 * (wmax - wmin),
            rf"$m={m}$",
            ha="center",
            va="top",
            fontsize=9,
        )

    # Strip centre lines.
    for k in range(-3, 4):
        w_line = (
            u_grid * args.tan_alpha
            + k * pstrip / np.cos(alpha)
        )
        focus = k in (-1, 0, 1)
        ax_top.plot(
            u_grid,
            w_line,
            linestyle="--",
            linewidth=1.8 if focus else 0.9,
            alpha=0.72 if focus else 0.38,
            color="tab:red" if focus else "0.55",
            zorder=1,
        )

    # Intersections of k=-1,0,+1 with selected wire.
    for k, wk in zip(focus_ks, focus_w):
        ax_top.plot(
            [u_wire],
            [wk],
            marker="o",
            color="black",
            markersize=5.5,
            linestyle="none",
            zorder=5,
        )
        ax_top.annotate(
            rf"$k={k:+d}$",
            xy=(u_wire, wk),
            xytext=(7, 4),
            textcoords="offset points",
            fontsize=10,
        )

    ax_top.plot(
        [u_wire],
        [w0],
        marker="*",
        markersize=20,
        color="tab:green",
        markeredgecolor="black",
        markeredgewidth=0.5,
        linestyle="none",
        zorder=7,
    )
    ax_top.axhline(w0, color="0.7", linewidth=0.8, alpha=0.55)

    w_center = strip_center_w(
        u_wire, args.tan_alpha, pstrip, 0
    )
    ax_top.plot(
        [u_wire],
        [w_center],
        marker="x",
        markersize=10,
        markeredgewidth=2,
        color="tab:red",
        linestyle="none",
        zorder=6,
    )

    ax_top.annotate(
        r"strip 0 centre on this wire",
        xy=(u_wire, w_center),
        xytext=(20, 18),
        textcoords="offset points",
        arrowprops=dict(arrowstyle="->", linewidth=0.9),
        fontsize=9,
    )

    info = (
        rf"$w_0={w0:+.3f}$ mm" + "\n"
        + rf"$w_{{\rm center}}=u\tan\alpha={w_center:+.3f}$ mm" + "\n"
        + rf"$x_\alpha=w_0\cos\alpha-u\sin\alpha={xalpha:+.3f}$ mm" + "\n"
        + rf"$x_\alpha/p={eta:+.3f}$" + "\n"
        + rf"nearest local strip centre: $k={nearest_index:+d}$"
    )
    ax_top.text(
        0.03,
        0.97,
        info,
        transform=ax_top.transAxes,
        ha="left",
        va="top",
        fontsize=10,
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.92),
    )

    ax_top.set_xlim(umin, umax)
    ax_top.set_ylim(wmin, wmax)
    ax_top.set_xlabel(r"$u$ across anode wires [mm]")
    ax_top.set_ylabel(r"$w$ along anode wires [mm]")
    ax_top.set_title(
        r"Top view: tilted strip centres crossing the selected wire"
    )
    ax_top.grid(alpha=0.18)

    # ------------------------------------------------------------------
    # Right: reduce the geometry to one selected wire.
    # ------------------------------------------------------------------
    ax_1d.set_title(
        r"Same geometry reduced to the selected wire"
    )

    y_av = 0.56
    y_centres = 0.30

    ax_1d.axhline(
        y_centres,
        color="0.75",
        linewidth=1.0,
    )

    for k, wk in zip(focus_ks, focus_w):
        ax_1d.axvline(
            wk,
            linestyle="--",
            linewidth=1.6,
            alpha=0.72,
            color="tab:red",
        )
        ax_1d.plot(
            [wk],
            [y_centres],
            marker="o",
            color="black",
            markersize=6,
        )
        ax_1d.text(
            wk,
            y_centres - 0.08,
            rf"$k={k:+d}$" + "\n" + rf"$w_k={wk:+.3f}$",
            ha="center",
            va="top",
            fontsize=9,
        )

    ax_1d.plot(
        [w0],
        [y_av],
        marker="*",
        markersize=22,
        color="tab:green",
        markeredgecolor="black",
        markeredgewidth=0.5,
        linestyle="none",
        label="avalanche",
    )
    ax_1d.axvline(
        w0,
        color="tab:green",
        linewidth=1.4,
        alpha=0.65,
    )

    # Explicit distance from avalanche to strip-0 centre.
    ax_1d.annotate(
        "",
        xy=(w_center, 0.76),
        xytext=(w0, 0.76),
        arrowprops=dict(arrowstyle="<->", linewidth=1.2),
    )
    ax_1d.text(
        0.5 * (w_center + w0),
        0.80,
        rf"$w_0-w_{{k=0}}={w0-w_center:+.3f}$ mm",
        ha="center",
        va="bottom",
        fontsize=10,
    )

    ax_1d.text(
        0.03,
        0.97,
        (
            r"On a tilted family, the strip index is set by"
            "\n"
            r"$x_\alpha=w\cos\alpha-u\sin\alpha$,"
            "\n"
            r"not by $w$ alone."
        ),
        transform=ax_1d.transAxes,
        ha="left",
        va="top",
        fontsize=10,
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.92),
    )

    ax_1d.set_xlim(wmin, wmax)
    ax_1d.set_ylim(0.0, 1.0)
    ax_1d.set_yticks([])
    ax_1d.set_xlabel(
        r"$w$ position along the selected anode wire [mm]"
    )
    ax_1d.grid(axis="x", alpha=0.20)

    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.90))
    fig.savefig(args.output, dpi=170)
    plt.close(fig)

    print("\n=== STAGE B3a GEOMETRY DEBUG ===")
    print(f"tan(alpha)             : {args.tan_alpha:g}")
    print(f"alpha                  : {alpha_deg:.4f} deg")
    print(f"selected wire u        : {u_wire:g} mm")
    print(f"avalanche w0           : {w0:g} mm")
    print(f"strip-0 centre         : {w_center:+.6f} mm")
    print(f"x_alpha                : {xalpha:+.6f} mm")
    print(f"x_alpha / pitch        : {eta:+.6f}")
    for k, wk in zip(focus_ks, focus_w):
        print(
            f"strip k={k:+d} centre      : {wk:+.6f} mm"
            f"  (w0-wk={w0-wk:+.6f} mm)"
        )
    print(f"nearest local strip    : {nearest_index:+d}")
    print(f"Wrote                  : {args.output}")


if __name__ == "__main__":
    main()
