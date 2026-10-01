#!/usr/bin/env python3

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.patches import Rectangle


def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Create a GIF of one positive ion drifting from the anode wire "
            "to the cathode while the induced strip current evolves in physical time."
        )
    )
    p.add_argument(
        "--prefix",
        default="phaseB_single_ion_2mm",
        help="Input prefix (default: phaseB_single_ion_2mm)",
    )
    p.add_argument(
        "--strip",
        default="strip_0",
        help="Strip current to display (default: strip_0)",
    )
    p.add_argument(
        "--strip-pitch-mm",
        type=float,
        default=1.7,
        help="Ideal strip pitch used for the drawing (default: 1.7 mm)",
    )
    p.add_argument(
        "--strip-width-mm",
        type=float,
        default=1.7,
        help="Ideal strip width used for the drawing (default: 1.7 mm)",
    )
    p.add_argument(
        "--half-strips",
        type=int,
        default=4,
        help="Number of strips drawn on each side of strip 0 (default: 4)",
    )
    p.add_argument(
        "--duration-s",
        type=float,
        default=12.0,
        help=(
            "Playback duration of the full physical drift (default: 12 s). "
            "Changing this only slows/speeds playback; frame times remain uniformly "
            "spaced in physical detector time."
        ),
    )
    p.add_argument(
        "--fps",
        type=int,
        default=20,
        help="GIF frame rate (default: 20 fps)",
    )
    p.add_argument(
        "--output",
        default=None,
        help="Output GIF name (default: <prefix>_<strip>_drift_current.gif)",
    )
    return p.parse_args()


def main():
    args = parse_args()

    if args.duration_s <= 0 or args.fps <= 0:
        raise ValueError("--duration-s and --fps must be positive.")
    if args.strip_pitch_mm <= 0 or args.strip_width_mm <= 0:
        raise ValueError("Strip pitch/width must be positive.")
    if args.strip_width_mm > args.strip_pitch_mm:
        raise ValueError("Strip width cannot exceed strip pitch in this ideal model.")
    if args.half_strips < 1:
        raise ValueError("--half-strips must be at least 1.")

    prefix = Path(args.prefix)
    drift_file = Path(str(prefix) + "_driftline.csv")
    wave_file = Path(str(prefix) + "_waveforms.csv")

    drift = pd.read_csv(drift_file)
    wave = pd.read_csv(wave_file)

    current_col = f"{args.strip}_fC_per_ns"
    if current_col not in wave.columns:
        raise KeyError(
            f"Column {current_col!r} not found in {wave_file}. "
            f"Available columns: {list(wave.columns)}"
        )

    required_drift = {"time_ns", "u_mm", "v_mm", "w_mm"}
    missing = required_drift.difference(drift.columns)
    if missing:
        raise KeyError(f"Missing driftline columns: {sorted(missing)}")

    td_ns = drift["time_ns"].to_numpy(dtype=float)
    u_mm = drift["u_mm"].to_numpy(dtype=float)
    v_mm = drift["v_mm"].to_numpy(dtype=float)
    w_mm = drift["w_mm"].to_numpy(dtype=float)

    tw_ns = wave["time_ns"].to_numpy(dtype=float)
    current = wave[current_col].to_numpy(dtype=float)

    if len(td_ns) < 2 or len(tw_ns) < 2:
        raise RuntimeError("Need at least two drift points and two waveform bins.")

    t_start_ns = max(float(td_ns[0]), 0.0)
    t_end_ns = float(td_ns[-1])
    if t_end_ns <= t_start_ns:
        raise RuntimeError("Invalid drift time range.")

    # Uniform physical-time sampling. Playback duration controls only how slowly
    # the detector's microsecond-scale motion is shown to the viewer.
    n_frames = max(2, int(round(args.duration_s * args.fps)))
    frame_t_ns = np.linspace(t_start_ns, t_end_ns, n_frames)

    frame_u = np.interp(frame_t_ns, td_ns, u_mm)
    frame_v = np.interp(frame_t_ns, td_ns, v_mm)
    frame_w = np.interp(frame_t_ns, td_ns, w_mm)
    frame_i = np.interp(frame_t_ns, tw_ns, current)

    # Draw the ideal Stage-B0/1 geometry in the local (w,v) plane.
    # For the 2 mm minus-side run, the readout cathode is reached at negative v.
    v_wire = 0.0
    v_cathode = float(v_mm[-1])
    gap_mm = abs(v_cathode - v_wire)

    strip_indices = np.arange(-args.half_strips, args.half_strips + 1)
    strip_centers = strip_indices * args.strip_pitch_mm

    w_extent = max(
        abs(strip_centers[0] - 0.5 * args.strip_width_mm),
        abs(strip_centers[-1] + 0.5 * args.strip_width_mm),
        2.5,
    )

    # Current axis only needs the physical drift interval.
    wave_mask = (tw_ns >= t_start_ns) & (tw_ns <= t_end_ns + 0.5 * np.median(np.diff(tw_ns)))
    tw_plot_us = tw_ns[wave_mask] / 1000.0
    current_plot = current[wave_mask]

    i_min = float(np.min(current_plot))
    i_max = float(np.max(current_plot))
    i_span = max(i_max - i_min, max(abs(i_max), abs(i_min)) * 0.1, 1e-30)

    fig = plt.figure(figsize=(11, 6))
    gs = fig.add_gridspec(2, 2, width_ratios=(1.0, 1.45), height_ratios=(1.0, 0.7))
    ax_geom = fig.add_subplot(gs[:, 0])
    ax_i = fig.add_subplot(gs[0, 1])
    ax_q = fig.add_subplot(gs[1, 1], sharex=ax_i)

    # -----------------------------
    # Left: chamber geometry.
    # -----------------------------
    ax_geom.set_xlim(-w_extent, w_extent)
    pad_v = 0.15 * gap_mm
    if v_cathode < 0:
        ax_geom.set_ylim(v_cathode - pad_v, v_wire + pad_v)
    else:
        ax_geom.set_ylim(v_wire - pad_v, v_cathode + pad_v)

    ax_geom.axhline(v_wire, linewidth=2)
    ax_geom.text(
        0.98,
        v_wire,
        "anode wire plane  v = 0",
        ha="right",
        va="bottom",
        transform=ax_geom.get_yaxis_transform(),
    )

    cathode_thickness = max(0.04 * gap_mm, 0.04)
    for k, center in zip(strip_indices, strip_centers):
        x0 = center - 0.5 * args.strip_width_mm
        y0 = v_cathode - 0.5 * cathode_thickness
        rect = Rectangle(
            (x0, y0),
            args.strip_width_mm,
            cathode_thickness,
            fill=(k == 0),
            alpha=0.35 if k == 0 else 0.12,
            linewidth=1.0,
        )
        ax_geom.add_patch(rect)
        if abs(k) <= 2:
            ax_geom.text(
                center,
                v_cathode - np.sign(v_cathode if v_cathode != 0 else -1) * 0.08 * gap_mm,
                f"{k}",
                ha="center",
                va="center",
                fontsize=9,
            )

    ax_geom.axhline(v_cathode, linewidth=1.2, linestyle="--")
    ax_geom.text(
        0.98,
        v_cathode,
        f"readout cathode  v = {v_cathode:.2f} mm",
        ha="right",
        va="top" if v_cathode < 0 else "bottom",
        transform=ax_geom.get_yaxis_transform(),
    )

    # Full trajectory as a faint reference curve plus animated ion marker.
    ax_geom.plot(w_mm, v_mm, linewidth=1.2, alpha=0.35)
    ion_marker, = ax_geom.plot([frame_w[0]], [frame_v[0]], marker="o", markersize=10)
    trail, = ax_geom.plot([], [], linewidth=2.0, alpha=0.65)

    ax_geom.set_xlabel("w [mm]  (along anode wire)")
    ax_geom.set_ylabel("v [mm]  (chamber normal)")
    ax_geom.set_title("Ion drift in the ideal 2 mm strip geometry")
    ax_geom.grid(alpha=0.25)

    time_text = ax_geom.text(
        0.04,
        0.96,
        "",
        transform=ax_geom.transAxes,
        ha="left",
        va="top",
        fontsize=11,
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.75),
    )

    # -----------------------------
    # Right top: induced current.
    # -----------------------------
    ax_i.plot(tw_plot_us, current_plot, linewidth=1.6)
    current_marker, = ax_i.plot([frame_t_ns[0] / 1000.0], [frame_i[0]], marker="o", markersize=7)
    current_cursor = ax_i.axvline(frame_t_ns[0] / 1000.0, linestyle="--", linewidth=1)

    ax_i.set_xlim(t_start_ns / 1000.0, t_end_ns / 1000.0)
    ax_i.set_ylim(i_min - 0.08 * i_span, i_max + 0.08 * i_span)
    ax_i.set_ylabel(f"{args.strip} current [fC/ns]")
    ax_i.set_title("Shockley–Ramo induced current")
    ax_i.grid(alpha=0.25)

    # -----------------------------
    # Right bottom: cumulative charge.
    # -----------------------------
    dt_wave_ns = float(np.median(np.diff(tw_ns)))
    q_cum = np.cumsum(current) * dt_wave_ns
    q_plot = q_cum[wave_mask]
    ax_q.plot(tw_plot_us, q_plot, linewidth=1.6)
    frame_q = np.interp(frame_t_ns, tw_ns, q_cum)
    q_marker, = ax_q.plot([frame_t_ns[0] / 1000.0], [frame_q[0]], marker="o", markersize=7)
    q_cursor = ax_q.axvline(frame_t_ns[0] / 1000.0, linestyle="--", linewidth=1)

    ax_q.set_xlabel("physical detector time [us]")
    ax_q.set_ylabel("cumulative Q [fC]")
    ax_q.set_title("Integrated induced charge")
    ax_q.grid(alpha=0.25)

    fig.suptitle(
        "Single positive ion: drift toward the 2 mm cathode and strip-0 signal",
        fontsize=14,
    )
    fig.tight_layout()

    def update(frame):
        t_ns = frame_t_ns[frame]
        t_us = t_ns / 1000.0

        ion_marker.set_data([frame_w[frame]], [frame_v[frame]])
        trail.set_data(frame_w[: frame + 1], frame_v[: frame + 1])

        current_marker.set_data([t_us], [frame_i[frame]])
        current_cursor.set_xdata([t_us, t_us])

        q_marker.set_data([t_us], [frame_q[frame]])
        q_cursor.set_xdata([t_us, t_us])

        time_text.set_text(
            f"physical time = {t_us:.2f} us\n"
            f"v = {frame_v[frame]:.3f} mm\n"
            f"I = {frame_i[frame]:.3e} fC/ns"
        )

        return (
            ion_marker,
            trail,
            current_marker,
            current_cursor,
            q_marker,
            q_cursor,
            time_text,
        )

    animation = FuncAnimation(
        fig,
        update,
        frames=n_frames,
        interval=1000.0 / args.fps,
        blit=False,
        repeat=True,
    )

    if args.output:
        out_gif = Path(args.output)
    else:
        out_gif = Path(str(prefix) + f"_{args.strip}_drift_current.gif")

    writer = PillowWriter(fps=args.fps)
    animation.save(out_gif, writer=writer)
    plt.close(fig)

    physical_duration_us = (t_end_ns - t_start_ns) / 1000.0
    slowdown = args.duration_s / (physical_duration_us * 1e-6)

    print(f"Driftline       : {drift_file}")
    print(f"Waveform        : {wave_file}")
    print(f"Strip           : {args.strip}")
    print(f"Physical drift  : {physical_duration_us:.6f} us")
    print(f"GIF duration    : {args.duration_s:.3f} s")
    print(f"Frames / fps    : {n_frames} / {args.fps}")
    print(f"Playback slowdown: {slowdown:.3e} x")
    print(f"Wrote           : {out_gif}")


if __name__ == "__main__":
    main()
