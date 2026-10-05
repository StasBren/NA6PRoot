#!/usr/bin/env python3
"""
Animate the Stage-B3a response of ONE tilted strip family while the avalanche
moves along one selected anode wire.

Input:
    stageB_tilt_scan.csv

No Garfield rerun is needed. The script uses the existing B3a scan, selects one
tilt angle and one wire, averages over the stored avalanche ensemble at each
simulated w0 point, and animates the response.

Coordinate convention used here:
    u : across anode wires, in the detector plane
    w : along anode wires, in the detector plane
    v : normal to the cathode planes

For a +alpha strip family, the coordinate perpendicular to the strips is

    x_alpha = w cos(alpha) - u sin(alpha).

Panels:
  1) Top view (u,w): parallel anode wires, tilted strips, moving avalanche.
  2) Early three-strip signal sharing A_-1, A_0, A_+1.
  3) Three-strip response versus the along-wire position w0.
  4) Left-right asymmetry versus projected coordinate x_alpha.

The scan stores event-by-event early induced charge. We use

    A_k = |Q_k(T_obs)|

and normalize the three central strips to

    A_3 = A_-1 + A_0 + A_+1.

For visual smoothness the GIF linearly interpolates between the simulated w0
points. Circular markers in the response curves are the actual simulated
positions; interpolation is display-only.
"""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.patches import Polygon
import numpy as np
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(
        description="Animate one Stage-B3a tilted strip family."
    )
    p.add_argument(
        "--input",
        default="stageB_tilt_scan.csv",
        help="Stage-B3a scan CSV.",
    )
    p.add_argument(
        "--tan-alpha",
        type=float,
        default=0.10,
        help="Tilt to select (default: tan(alpha)=0.10).",
    )
    p.add_argument(
        "--wire-index",
        type=int,
        default=1,
        help=(
            "Selected anode wire index m. Default +1 deliberately shows the "
            "tilt-induced response shift; m=0 looks almost like alpha=0."
        ),
    )
    p.add_argument(
        "--strip-pitch-mm",
        type=float,
        default=1.7,
        help="Strip pitch used in the B3a scan (default: 1.7 mm).",
    )
    p.add_argument(
        "--strip-width-mm",
        type=float,
        default=1.7,
        help="Strip width used in the B3a scan (default: 1.7 mm).",
    )
    p.add_argument(
        "--wire-pitch-mm",
        type=float,
        default=4.0,
        help="Anode-wire pitch used in the B3a scan (default: 4 mm).",
    )
    p.add_argument(
        "--duration-s",
        type=float,
        default=10.0,
        help="Total GIF playback duration including return sweep.",
    )
    p.add_argument(
        "--fps",
        type=int,
        default=20,
        help="GIF frame rate.",
    )
    p.add_argument(
        "--output",
        default=None,
        help="Output GIF name.",
    )
    return p.parse_args()


def nearest_value(values, target):
    values = np.asarray(sorted(set(values)), dtype=float)
    return float(values[np.argmin(np.abs(values - target))])


def ensemble_mean(df):
    cols = [
        "alpha_deg",
        "wire_u_mm",
        "w_shift_mm",
        "projected_x_mm",
        "projected_eta",
        "observation_ns",
        "three_strip_capture_fraction",
        "neighbor_fraction_three_strip",
        "central_fraction_three_strip",
        "left_right_asymmetry",
        "A_strip_-2_fC",
        "A_strip_-1_fC",
        "A_strip_0_fC",
        "A_strip_1_fC",
        "A_strip_2_fC",
    ]
    return (
        df.groupby("xi_w", as_index=False)[cols]
        .mean()
        .sort_values("xi_w")
        .reset_index(drop=True)
    )


def add_strip_band(ax, k, pitch, width, alpha, umin, umax):
    """
    Draw one ideal strip band in the (u,w) plane.

    Strip centre line:
        x_alpha = w cos(alpha) - u sin(alpha) = k p

    therefore
        w(u) = u tan(alpha) + k p / cos(alpha).

    The physical strip width is measured perpendicular to the strip, so the
    band boundaries are x_alpha = k p +/- width/2.
    """
    ca = np.cos(alpha)
    ta = np.tan(alpha)

    x_lo = k * pitch - 0.5 * width
    x_hi = k * pitch + 0.5 * width

    w_lo_umin = umin * ta + x_lo / ca
    w_lo_umax = umax * ta + x_lo / ca
    w_hi_umin = umin * ta + x_hi / ca
    w_hi_umax = umax * ta + x_hi / ca

    poly = Polygon(
        [
            (umin, w_lo_umin),
            (umax, w_lo_umax),
            (umax, w_hi_umax),
            (umin, w_hi_umin),
        ],
        closed=True,
        alpha=0.18,
        linewidth=1.0,
    )
    ax.add_patch(poly)

    # Dashed strip centre.
    us = np.array([umin, umax])
    ws = us * ta + (k * pitch) / ca
    ax.plot(us, ws, linestyle="--", linewidth=0.9, alpha=0.55)


def main():
    args = parse_args()

    if args.duration_s <= 0:
        raise ValueError("--duration-s must be positive")
    if args.fps <= 0:
        raise ValueError("--fps must be positive")
    if args.strip_pitch_mm <= 0 or args.strip_width_mm <= 0:
        raise ValueError("Strip pitch/width must be positive")
    if args.strip_width_mm > args.strip_pitch_mm + 1e-12:
        raise ValueError("Strip width must not exceed strip pitch")

    df = pd.read_csv(args.input)

    required = {
        "tan_alpha",
        "alpha_deg",
        "wire_index",
        "wire_u_mm",
        "xi_w",
        "w_shift_mm",
        "projected_x_mm",
        "projected_eta",
        "observation_ns",
        "left_right_asymmetry",
        "A_strip_-2_fC",
        "A_strip_-1_fC",
        "A_strip_0_fC",
        "A_strip_1_fC",
        "A_strip_2_fC",
    }
    missing = required - set(df.columns)
    if missing:
        raise KeyError(f"Missing columns in {args.input}: {sorted(missing)}")

    tan_alpha = nearest_value(df["tan_alpha"].unique(), args.tan_alpha)
    available_wires = np.asarray(
        sorted(set(df.loc[np.isclose(df["tan_alpha"], tan_alpha), "wire_index"])),
        dtype=int,
    )
    if len(available_wires) == 0:
        raise RuntimeError("No wire indices found for selected tilt.")

    wire_index = int(
        available_wires[np.argmin(np.abs(available_wires - args.wire_index))]
    )

    sel = df[
        np.isclose(df["tan_alpha"], tan_alpha)
        & (df["wire_index"].astype(int) == wire_index)
    ].copy()

    if sel.empty:
        raise RuntimeError(
            f"No rows found for tan(alpha)={tan_alpha:g}, wire m={wire_index}"
        )

    mean = ensemble_mean(sel)

    alpha_deg = float(mean["alpha_deg"].iloc[0])
    alpha = np.deg2rad(alpha_deg)
    wire_u = float(mean["wire_u_mm"].iloc[0])
    obs_ns = float(mean["observation_ns"].iloc[0])

    xi_data = mean["xi_w"].to_numpy(dtype=float)
    w_data = mean["w_shift_mm"].to_numpy(dtype=float)
    xalpha_data = mean["projected_x_mm"].to_numpy(dtype=float)

    a_m1 = mean["A_strip_-1_fC"].to_numpy(dtype=float)
    a_0 = mean["A_strip_0_fC"].to_numpy(dtype=float)
    a_p1 = mean["A_strip_1_fC"].to_numpy(dtype=float)
    a3 = a_m1 + a_0 + a_p1

    if np.any(a3 <= 0):
        raise RuntimeError("Non-positive three-strip signal encountered.")

    f_m1 = 100.0 * a_m1 / a3
    f_0 = 100.0 * a_0 / a3
    f_p1 = 100.0 * a_p1 / a3
    asym = mean["left_right_asymmetry"].to_numpy(dtype=float)

    # Smooth ping-pong sweep.
    n_total = max(8, int(round(args.duration_s * args.fps)))
    n_forward = max(4, n_total // 2 + 1)

    w_forward = np.linspace(float(w_data.min()), float(w_data.max()), n_forward)
    w_backward = w_forward[-2:0:-1]
    w_path = np.concatenate([w_forward, w_backward])

    phase_src = np.linspace(0.0, 1.0, len(w_path))
    phase_dst = np.linspace(0.0, 1.0, n_total, endpoint=False)
    w_frames = np.interp(phase_dst, phase_src, w_path)

    # Since x_alpha is linear in w at fixed u, direct interpolation in w is safe.
    fm1_frames = np.interp(w_frames, w_data, f_m1)
    f0_frames = np.interp(w_frames, w_data, f_0)
    fp1_frames = np.interp(w_frames, w_data, f_p1)
    asym_frames = np.interp(w_frames, w_data, asym)
    a3_frames = np.interp(w_frames, w_data, a3)

    xalpha_frames = (
        w_frames * np.cos(alpha) - wire_u * np.sin(alpha)
    )
    eta_frames = xalpha_frames / args.strip_pitch_mm

    expected_center = wire_u * tan_alpha

    # ------------------------------------------------------------------
    # Figure layout.
    # ------------------------------------------------------------------
    fig = plt.figure(figsize=(13.2, 8.0))
    gs = fig.add_gridspec(
        2,
        2,
        width_ratios=(1.0, 1.20),
        height_ratios=(1.06, 0.94),
    )

    ax_geom = fig.add_subplot(gs[0, 0])
    ax_bars = fig.add_subplot(gs[1, 0])
    ax_resp = fig.add_subplot(gs[0, 1])
    ax_asym = fig.add_subplot(gs[1, 1])

    fig.suptitle(
        rf"Stage B3a: response of one tilted strip family"
        "\n"
        rf"$\tan\alpha={tan_alpha:g}$ ($\alpha={alpha_deg:.2f}^\circ$), "
        rf"wire $m={wire_index}$ at $u={wire_u:g}$ mm, "
        rf"$T_{{\mathrm{{obs}}}}={obs_ns:g}$ ns",
        fontsize=14,
    )

    # ------------------------------------------------------------------
    # Panel 1: true top view in (u,w).
    # Wires run along w and are separated along u.
    # ------------------------------------------------------------------
    wire_indices_draw = np.arange(-2, 3)
    u_wires = wire_indices_draw * args.wire_pitch_mm

    umin = float(u_wires.min() - 1.2)
    umax = float(u_wires.max() + 1.2)

    w_margin = 2.0 * args.strip_pitch_mm
    wmin = min(float(w_data.min()), expected_center) - w_margin
    wmax = max(float(w_data.max()), expected_center) + w_margin

    # Draw tilted strips k=-3..+3.
    for k in range(-3, 4):
        add_strip_band(
            ax_geom,
            k,
            args.strip_pitch_mm,
            args.strip_width_mm,
            alpha,
            umin,
            umax,
        )

    # Draw anode wires. Each wire is one line parallel to w.
    for m, u in zip(wire_indices_draw, u_wires):
        lw = 2.2 if m == wire_index else 1.0
        alpha_line = 0.95 if m == wire_index else 0.45
        ax_geom.plot(
            [u, u],
            [wmin, wmax],
            linewidth=lw,
            alpha=alpha_line,
        )
        ax_geom.text(
            u,
            wmax - 0.05 * (wmax - wmin),
            rf"$m={m}$",
            ha="center",
            va="top",
            fontsize=8,
        )

    avalanche_marker, = ax_geom.plot(
        [wire_u],
        [w_frames[0]],
        marker="*",
        markersize=18,
        linestyle="none",
        zorder=6,
    )

    # Response-center point x_alpha=0 on the selected wire.
    center_marker, = ax_geom.plot(
        [wire_u],
        [expected_center],
        marker="x",
        markersize=9,
        linestyle="none",
        zorder=5,
    )

    motion_line, = ax_geom.plot(
        [wire_u, wire_u],
        [w_data.min(), w_data.max()],
        linestyle=":",
        linewidth=1.3,
        alpha=0.65,
    )

    geom_info = ax_geom.text(
        0.03,
        0.97,
        "",
        transform=ax_geom.transAxes,
        ha="left",
        va="top",
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.88),
    )

    ax_geom.set_xlim(umin, umax)
    ax_geom.set_ylim(wmin, wmax)
    ax_geom.set_xlabel(r"$u$ across anode wires [mm]")
    ax_geom.set_ylabel(r"$w$ along anode wires [mm]")
    ax_geom.set_title(
        r"Top view: wires $\parallel w$, tilted strips in the $(u,w)$ plane"
    )
    ax_geom.grid(alpha=0.18)

    # ------------------------------------------------------------------
    # Panel 2: current frame's three-strip early sharing.
    # ------------------------------------------------------------------
    bar_ks = np.array([-1, 0, 1])
    bars = ax_bars.bar(
        bar_ks,
        [fm1_frames[0], f0_frames[0], fp1_frames[0]],
    )

    ax_bars.set_xlim(-1.75, 1.75)
    ax_bars.set_ylim(
        0.0,
        max(70.0, 1.12 * np.max([f_m1, f_0, f_p1])),
    )
    ax_bars.set_xticks(bar_ks)
    ax_bars.set_xticklabels([r"$-1$", r"$0$", r"$+1$"])
    ax_bars.set_xlabel("strip index")
    ax_bars.set_ylabel(r"$A_k/A_3$ [%]")
    ax_bars.set_title(
        r"Early signal sharing, "
        r"$A_k=|Q_k(T_{\mathrm{obs}})|$"
    )
    ax_bars.grid(axis="y", alpha=0.22)

    # ------------------------------------------------------------------
    # Panel 3: spatial response versus ordinary along-wire coordinate w0.
    # ------------------------------------------------------------------
    ax_resp.plot(w_data, f_m1, marker="o", label=r"strip $-1$")
    ax_resp.plot(w_data, f_0, marker="o", label=r"strip $0$")
    ax_resp.plot(w_data, f_p1, marker="o", label=r"strip $+1$")

    resp_m1, = ax_resp.plot(
        [w_frames[0]], [fm1_frames[0]],
        marker="o", markersize=9, linestyle="none",
    )
    resp_0, = ax_resp.plot(
        [w_frames[0]], [f0_frames[0]],
        marker="o", markersize=9, linestyle="none",
    )
    resp_p1, = ax_resp.plot(
        [w_frames[0]], [fp1_frames[0]],
        marker="o", markersize=9, linestyle="none",
    )
    resp_cursor = ax_resp.axvline(
        w_frames[0], linestyle="--", linewidth=1.1
    )
    ax_resp.axvline(
        expected_center,
        linestyle=":",
        linewidth=1.2,
        label=rf"$x_\alpha=0$: $w=u\tan\alpha={expected_center:+.2f}$ mm",
    )

    ax_resp.set_xlabel(r"along-wire avalanche position $w_0$ [mm]")
    ax_resp.set_ylabel("fraction of three-strip early signal [%]")
    ax_resp.set_title("Three-strip response on the selected anode wire")
    ax_resp.legend(fontsize=8)
    ax_resp.grid(alpha=0.22)

    # ------------------------------------------------------------------
    # Panel 4: asymmetry versus projected coordinate x_alpha.
    # This is the natural coordinate of the tilted strip family.
    # ------------------------------------------------------------------
    order_x = np.argsort(xalpha_data)
    ax_asym.plot(
        xalpha_data[order_x],
        asym[order_x],
        marker="o",
    )
    asym_marker, = ax_asym.plot(
        [xalpha_frames[0]],
        [asym_frames[0]],
        marker="o",
        markersize=9,
        linestyle="none",
    )
    asym_cursor = ax_asym.axvline(
        xalpha_frames[0],
        linestyle="--",
        linewidth=1.1,
    )
    ax_asym.axhline(0.0, linewidth=0.9)
    ax_asym.axvline(0.0, linewidth=0.9)

    ax_asym.set_xlabel(
        r"projected coordinate "
        r"$x_\alpha=w_0\cos\alpha-u_{\mathrm{wire}}\sin\alpha$ [mm]"
    )
    ax_asym.set_ylabel(
        r"$R=(A_{+1}-A_{-1})/(A_{-1}+A_0+A_{+1})$"
    )
    ax_asym.set_title(
        r"Left-right asymmetry in the strip family's natural coordinate"
    )
    ax_asym.grid(alpha=0.22)

    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.93))

    def update(i):
        w0 = float(w_frames[i])
        xalpha = float(xalpha_frames[i])

        # Top-view geometry.
        avalanche_marker.set_data([wire_u], [w0])
        geom_info.set_text(
            rf"$w_0={w0:+.3f}$ mm" + "\n"
            + rf"$x_\alpha={xalpha:+.3f}$ mm" + "\n"
            + rf"$x_\alpha/p={eta_frames[i]:+.3f}$" + "\n"
            + rf"$R={asym_frames[i]:+.3f}$"
        )

        # Bars.
        for bar, value in zip(
            bars,
            [fm1_frames[i], f0_frames[i], fp1_frames[i]],
        ):
            bar.set_height(float(value))

        # Response moving markers.
        resp_m1.set_data([w0], [fm1_frames[i]])
        resp_0.set_data([w0], [f0_frames[i]])
        resp_p1.set_data([w0], [fp1_frames[i]])
        resp_cursor.set_xdata([w0, w0])

        # Asymmetry in projected coordinate.
        asym_marker.set_data([xalpha], [asym_frames[i]])
        asym_cursor.set_xdata([xalpha, xalpha])

        return (
            avalanche_marker,
            center_marker,
            motion_line,
            geom_info,
            *bars,
            resp_m1,
            resp_0,
            resp_p1,
            resp_cursor,
            asym_marker,
            asym_cursor,
        )

    ani = FuncAnimation(
        fig,
        update,
        frames=len(w_frames),
        interval=1000.0 / args.fps,
        blit=False,
        repeat=True,
    )

    if args.output:
        output = Path(args.output)
    else:
        output = Path(
            f"stageB_tilt_position_response_tan{tan_alpha:g}"
            f"_wire{wire_index:+d}.gif"
        )

    ani.save(output, writer=PillowWriter(fps=args.fps))
    plt.close(fig)

    print("\n=== STAGE B3a TILTED-FAMILY POSITION-RESPONSE GIF ===")
    print(f"input                  : {args.input}")
    print(f"selected tan(alpha)    : {tan_alpha:g}")
    print(f"alpha                  : {alpha_deg:.4f} deg")
    print(f"selected wire index    : {wire_index}")
    print(f"selected wire u        : {wire_u:g} mm")
    print(f"strip pitch / width    : {args.strip_pitch_mm:g} / {args.strip_width_mm:g} mm")
    print(f"observation time       : {obs_ns:g} ns")
    print(f"expected R=0 center    : w = u tan(alpha) = {expected_center:+.4f} mm")
    print(f"simulated w0 points    : {len(w_data)}")
    print(f"ensemble events/w0     : approximately {len(sel) // len(w_data)}")
    print(f"animation frames       : {len(w_frames)}")
    print(f"duration               : {args.duration_s:g} s")
    print(f"fps                    : {args.fps}")
    print("NOTE: animation is linearly interpolated between simulated w0 points.")
    print(f"Wrote                  : {output}")


if __name__ == "__main__":
    main()
