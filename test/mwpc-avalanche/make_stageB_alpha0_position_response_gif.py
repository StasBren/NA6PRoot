#!/usr/bin/env python3
"""
Animate the Stage-B2 alpha=0 strip response while the stored microscopic
avalanche is translated across one strip pitch.

Input:
    stageB_strip_geometry_scan.csv

This script does NOT run Garfield again. It uses the already generated Stage-B2
scan, selects one (pitch, width) geometry, averages over the stored avalanche
ensemble at each simulated xi = w0 / p, and animates the response versus w0.

Panels:
  1) 1D readout geometry and moving avalanche position w0.
  2) Early signal sharing among strips -1, 0, +1.
  3) Full three-strip response curves with the current w0 highlighted.
  4) Left-right asymmetry
         R = (A_{+1} - A_{-1}) / (A_{-1} + A_0 + A_{+1}).

Here
    A_k = |Q_k(T_obs)|
is the magnitude of the accumulated induced charge at the fixed observation
time stored in the Stage-B2 file (100 ns in the standard scan).

For visual smoothness the animation linearly interpolates between the simulated
xi points. The circular markers in the response panels are the actual simulated
positions; interpolation is display-only.
"""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(
        description="Animate Stage-B2 alpha=0 response across one strip pitch."
    )
    p.add_argument(
        "--input",
        default="stageB_strip_geometry_scan.csv",
        help="Stage-B2 scan CSV.",
    )
    p.add_argument(
        "--pitch-mm",
        type=float,
        default=1.7,
        help="Strip pitch to select (default: 1.7 mm).",
    )
    p.add_argument(
        "--width-mm",
        type=float,
        default=1.7,
        help="Strip width to select (default: 1.7 mm).",
    )
    p.add_argument(
        "--duration-s",
        type=float,
        default=10.0,
        help="Total GIF playback duration including return sweep (default: 10 s).",
    )
    p.add_argument(
        "--fps",
        type=int,
        default=20,
        help="GIF frame rate (default: 20 fps).",
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
        "true_w_mm",
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
        df.groupby("xi", as_index=False)[cols]
        .mean()
        .sort_values("xi")
        .reset_index(drop=True)
    )


def main():
    args = parse_args()

    if args.duration_s <= 0:
        raise ValueError("--duration-s must be positive")
    if args.fps <= 0:
        raise ValueError("--fps must be positive")

    df = pd.read_csv(args.input)

    required = {
        "strip_pitch_mm",
        "strip_width_mm",
        "xi",
        "true_w_mm",
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

    pitch = nearest_value(df["strip_pitch_mm"].unique(), args.pitch_mm)
    width = nearest_value(df["strip_width_mm"].unique(), args.width_mm)

    sel = df[
        np.isclose(df["strip_pitch_mm"], pitch)
        & np.isclose(df["strip_width_mm"], width)
    ].copy()

    if sel.empty:
        raise RuntimeError(
            f"No rows found for pitch={pitch:g} mm, width={width:g} mm"
        )

    mean = ensemble_mean(sel)

    xi_data = mean["xi"].to_numpy(dtype=float)
    w_data = mean["true_w_mm"].to_numpy(dtype=float)
    obs_ns = float(mean["observation_ns"].iloc[0])

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

    # A smooth ping-pong sweep avoids a discontinuous jump at GIF repeat.
    n_total = max(8, int(round(args.duration_s * args.fps)))
    n_forward = max(4, n_total // 2 + 1)

    xi_forward = np.linspace(float(xi_data.min()), float(xi_data.max()), n_forward)
    xi_backward = xi_forward[-2:0:-1]
    xi_frames = np.concatenate([xi_forward, xi_backward])

    # Adjust frame count by resampling the ping-pong coordinate to the requested
    # total count. This preserves a smooth loop even when n_total is odd.
    phase_src = np.linspace(0.0, 1.0, len(xi_frames))
    phase_dst = np.linspace(0.0, 1.0, n_total, endpoint=False)
    xi_frames = np.interp(phase_dst, phase_src, xi_frames)

    w_frames = xi_frames * pitch
    fm1_frames = np.interp(xi_frames, xi_data, f_m1)
    f0_frames = np.interp(xi_frames, xi_data, f_0)
    fp1_frames = np.interp(xi_frames, xi_data, f_p1)
    asym_frames = np.interp(xi_frames, xi_data, asym)
    a3_frames = np.interp(xi_frames, xi_data, a3)

    # ------------------------------------------------------------------
    # Figure layout.
    # ------------------------------------------------------------------
    fig = plt.figure(figsize=(12.5, 7.5))
    gs = fig.add_gridspec(
        2,
        2,
        width_ratios=(1.0, 1.28),
        height_ratios=(0.92, 1.08),
    )

    ax_geom = fig.add_subplot(gs[0, 0])
    ax_bars = fig.add_subplot(gs[1, 0])
    ax_resp = fig.add_subplot(gs[0, 1])
    ax_asym = fig.add_subplot(gs[1, 1])

    fig.suptitle(
        rf"Stage B2: early strip response while $w_0$ moves across one pitch"
        "\n"
        rf"$\alpha=0$, $p={pitch:g}$ mm, $s={width:g}$ mm, "
        rf"$T_{{\rm obs}}={obs_ns:g}$ ns",
        fontsize=14,
    )

    # ------------------------------------------------------------------
    # Panel 1: readout geometry in the sensitive w coordinate.
    # ------------------------------------------------------------------
    ks = np.array([-2, -1, 0, 1, 2], dtype=int)
    strip_centers = ks * pitch

    geom_y = 0.0
    strip_h = 0.22
    for k, center in zip(ks, strip_centers):
        rect = Rectangle(
            (center - 0.5 * width, geom_y - 0.5 * strip_h),
            width,
            strip_h,
            alpha=0.22,
            linewidth=1.2,
        )
        ax_geom.add_patch(rect)
        ax_geom.text(
            center,
            geom_y - 0.28,
            rf"$k={k:+d}$" if k != 0 else r"$k=0$",
            ha="center",
            va="top",
            fontsize=9,
        )

    ax_geom.axhline(geom_y, linewidth=1.0)
    ax_geom.axvline(0.0, linestyle="--", linewidth=1.0, alpha=0.65)

    avalanche_marker, = ax_geom.plot(
        [w_frames[0]],
        [0.75],
        marker="*",
        markersize=18,
        linestyle="none",
    )
    avalanche_drop = ax_geom.vlines(
        w_frames[0],
        geom_y + 0.12,
        0.68,
        linestyles="--",
        linewidth=1.2,
    )

    w0_arrow = ax_geom.annotate(
        "",
        xy=(w_frames[0], 0.52),
        xytext=(0.0, 0.52),
        arrowprops=dict(arrowstyle="<->", linewidth=1.2),
    )
    w0_text = ax_geom.text(
        0.5 * w_frames[0],
        0.57,
        "",
        ha="center",
        va="bottom",
    )

    info = ax_geom.text(
        0.03,
        0.95,
        "",
        transform=ax_geom.transAxes,
        ha="left",
        va="top",
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.85),
    )

    x_extent = 2.75 * pitch
    ax_geom.set_xlim(-x_extent, x_extent)
    ax_geom.set_ylim(-0.52, 1.05)
    ax_geom.set_yticks([])
    ax_geom.set_xlabel(r"sensitive coordinate $w$ [mm]")
    ax_geom.set_title("Readout geometry and avalanche position")
    ax_geom.grid(axis="x", alpha=0.22)

    # ------------------------------------------------------------------
    # Panel 2: instantaneous frame's three-strip early signal sharing.
    # A_k = |Q_k(T_obs)|, normalized to A_-1 + A_0 + A_+1.
    # ------------------------------------------------------------------
    bar_ks = np.array([-1, 0, 1])
    bar_values = np.array([fm1_frames[0], f0_frames[0], fp1_frames[0]])
    bars = ax_bars.bar(bar_ks, bar_values)

    ax_bars.set_xlim(-1.75, 1.75)
    ax_bars.set_ylim(0.0, max(70.0, 1.12 * np.max([f_m1, f_0, f_p1])))
    ax_bars.set_xticks(bar_ks)
    ax_bars.set_xticklabels([r"$-1$", r"$0$", r"$+1$"])
    ax_bars.set_xlabel("strip index")
    ax_bars.set_ylabel(r"$A_k/A_3$ [%]")
    ax_bars.set_title(
        r"Early signal sharing: $A_k=|Q_k(T_{m obs})|$, "
        r"$A_3=A_{-1}+A_0+A_{+1}$"
    )
    ax_bars.grid(axis="y", alpha=0.22)

    # ------------------------------------------------------------------
    # Panel 3: full response curves. Markers are actual simulated xi points.
    # ------------------------------------------------------------------
    ax_resp.plot(xi_data, f_m1, marker="o", label=r"strip $-1$")
    ax_resp.plot(xi_data, f_0, marker="o", label=r"strip $0$")
    ax_resp.plot(xi_data, f_p1, marker="o", label=r"strip $+1$")

    resp_m1, = ax_resp.plot(
        [xi_frames[0]], [fm1_frames[0]], marker="o", markersize=9, linestyle="none"
    )
    resp_0, = ax_resp.plot(
        [xi_frames[0]], [f0_frames[0]], marker="o", markersize=9, linestyle="none"
    )
    resp_p1, = ax_resp.plot(
        [xi_frames[0]], [fp1_frames[0]], marker="o", markersize=9, linestyle="none"
    )
    resp_cursor = ax_resp.axvline(xi_frames[0], linestyle="--", linewidth=1.1)

    ax_resp.set_xlabel(r"normalized avalanche position $\xi=w_0/p$")
    ax_resp.set_ylabel("fraction of three-strip early signal [%]")
    ax_resp.set_title("Three-strip spatial response")
    ax_resp.legend()
    ax_resp.grid(alpha=0.22)

    # ------------------------------------------------------------------
    # Panel 4: left-right asymmetry.
    # ------------------------------------------------------------------
    ax_asym.plot(xi_data, asym, marker="o")
    asym_marker, = ax_asym.plot(
        [xi_frames[0]],
        [asym_frames[0]],
        marker="o",
        markersize=9,
        linestyle="none",
    )
    asym_cursor = ax_asym.axvline(xi_frames[0], linestyle="--", linewidth=1.1)
    ax_asym.axhline(0.0, linewidth=0.9)

    ax_asym.set_xlabel(r"normalized avalanche position $\xi=w_0/p$")
    ax_asym.set_ylabel(
        r"$R=(A_{+1}-A_{-1})/(A_{-1}+A_0+A_{+1})$"
    )
    ax_asym.set_title("Left-right asymmetry: position-sensitive observable")
    ax_asym.grid(alpha=0.22)

    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.93))

    def update(i):
        xi = float(xi_frames[i])
        w0 = float(w_frames[i])

        # Geometry.
        avalanche_marker.set_data([w0], [0.75])

        # vlines returns a LineCollection; update its one segment.
        avalanche_drop.set_segments([[(w0, geom_y + 0.12), (w0, 0.68)]])

        w0_arrow.xy = (w0, 0.52)
        w0_arrow.set_position((0.0, 0.52))
        w0_text.set_position((0.5 * w0, 0.57))
        w0_text.set_text(rf"$w_0={w0:+.3f}$ mm")

        info.set_text(
            rf"$\xi=w_0/p={xi:+.3f}$" + "\n"
            + rf"$A_3={a3_frames[i]:.4g}$ fC" + "\n"
            + rf"$R={asym_frames[i]:+.3f}$"
        )

        # Bars.
        values = [fm1_frames[i], f0_frames[i], fp1_frames[i]]
        for bar, value in zip(bars, values):
            bar.set_height(float(value))

        # Response moving markers.
        resp_m1.set_data([xi], [fm1_frames[i]])
        resp_0.set_data([xi], [f0_frames[i]])
        resp_p1.set_data([xi], [fp1_frames[i]])
        resp_cursor.set_xdata([xi, xi])

        # Asymmetry moving marker.
        asym_marker.set_data([xi], [asym_frames[i]])
        asym_cursor.set_xdata([xi, xi])

        return (
            avalanche_marker,
            avalanche_drop,
            w0_arrow,
            w0_text,
            info,
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
        frames=len(xi_frames),
        interval=1000.0 / args.fps,
        blit=False,
        repeat=True,
    )

    if args.output:
        output = Path(args.output)
    else:
        output = Path(
            f"stageB_alpha0_position_response_p{pitch:g}_s{width:g}.gif"
        )

    ani.save(output, writer=PillowWriter(fps=args.fps))
    plt.close(fig)

    print("\n=== STAGE B2 alpha=0 POSITION-RESPONSE GIF ===")
    print(f"input                  : {args.input}")
    print(f"selected pitch         : {pitch:g} mm")
    print(f"selected strip width   : {width:g} mm")
    print(f"observation time       : {obs_ns:g} ns")
    print(f"simulated xi points    : {len(xi_data)}")
    print(f"ensemble events/xi     : approximately {len(sel) // len(xi_data)}")
    print(f"animation frames       : {len(xi_frames)}")
    print(f"duration               : {args.duration_s:g} s")
    print(f"fps                    : {args.fps}")
    print("NOTE: animation is linearly interpolated between simulated xi points.")
    print(f"Wrote                  : {output}")


if __name__ == "__main__":
    main()
