#!/usr/bin/env python3

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.patches import Rectangle


def strip_label(k: int) -> str:
    if k < 0:
        return f"strip_m{-k}"
    if k > 0:
        return f"strip_p{k}"
    return "strip_0"


def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Animate one positive ion in physical detector time together with "
            "induced currents and cumulative induced charge on neighbouring strips."
        )
    )
    p.add_argument(
        "--prefix",
        default="phaseB_single_ion_2mm",
        help="Input prefix (default: phaseB_single_ion_2mm)",
    )
    p.add_argument(
        "--shown-half-strips",
        type=int,
        default=2,
        help="Show currents/sharing for strips -N..+N (default: 2)",
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
        "--duration-s",
        type=float,
        default=12.0,
        help="Playback duration of the full physical drift (default: 12 s)",
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
        help="Output GIF name",
    )
    return p.parse_args()


def main():
    args = parse_args()
    if args.shown_half_strips < 1:
        raise ValueError("--shown-half-strips must be >= 1")
    if args.duration_s <= 0 or args.fps <= 0:
        raise ValueError("--duration-s and --fps must be positive")

    prefix = Path(args.prefix)
    drift_file = Path(str(prefix) + "_driftline.csv")
    wave_file = Path(str(prefix) + "_waveforms.csv")
    drift = pd.read_csv(drift_file)
    wave = pd.read_csv(wave_file)

    ks = np.arange(-args.shown_half_strips, args.shown_half_strips + 1)
    labels = [strip_label(int(k)) for k in ks]

    for label in labels:
        if f"{label}_fC_per_ns" not in wave.columns:
            raise KeyError(f"Missing waveform column {label}_fC_per_ns")

    td_ns = drift["time_ns"].to_numpy(dtype=float)
    vd_mm = drift["v_mm"].to_numpy(dtype=float)
    wd_mm = drift["w_mm"].to_numpy(dtype=float)

    tw_ns = wave["time_ns"].to_numpy(dtype=float)
    t0_ns = max(float(td_ns[0]), 0.0)
    t1_ns = float(td_ns[-1])
    if t1_ns <= t0_ns:
        raise RuntimeError("Invalid drift time range")

    n_frames = max(2, int(round(args.duration_s * args.fps)))
    frame_t_ns = np.linspace(t0_ns, t1_ns, n_frames)
    frame_t_us = frame_t_ns / 1000.0

    frame_v = np.interp(frame_t_ns, td_ns, vd_mm)
    frame_w = np.interp(frame_t_ns, td_ns, wd_mm)

    frame_currents = {}
    for label in labels:
        frame_currents[label] = np.interp(
            frame_t_ns,
            tw_ns,
            wave[f"{label}_fC_per_ns"].to_numpy(dtype=float),
        )

    # Cumulative induced charge on each shown strip.
    # This is the detector-level quantity Q_k(t) = integral_0^t I_k(t') dt'.
    dt_wave_ns = float(np.median(np.diff(tw_ns)))
    q_cum = {}
    frame_q = {}
    for label in labels:
        current = wave[f"{label}_fC_per_ns"].to_numpy(dtype=float)
        q_cum[label] = np.cumsum(current) * dt_wave_ns
        frame_q[label] = np.interp(frame_t_ns, tw_ns, q_cum[label])

    central_current = wave["strip_0_fC_per_ns"].to_numpy(dtype=float)
    q0_cum = q_cum["strip_0"]
    frame_q0 = frame_q["strip_0"]

    physical_mask = (tw_ns >= t0_ns) & (tw_ns <= t1_ns + 0.5 * dt_wave_ns)
    tw_plot_us = tw_ns[physical_mask] / 1000.0

    fig = plt.figure(figsize=(12, 7))
    gs = fig.add_gridspec(2, 2, width_ratios=(1.0, 1.35), height_ratios=(1.0, 0.85))
    ax_geom = fig.add_subplot(gs[0, 0])
    ax_qbars = fig.add_subplot(gs[1, 0])
    ax_current = fig.add_subplot(gs[0, 1])
    ax_q = fig.add_subplot(gs[1, 1], sharex=ax_current)

    # --------------------------------------------------------------
    # Geometry panel.
    # --------------------------------------------------------------
    v_wire = 0.0
    v_cathode = float(vd_mm[-1])
    gap_mm = abs(v_cathode - v_wire)
    strip_centers = ks * args.strip_pitch_mm
    w_extent = max(abs(strip_centers[0]), abs(strip_centers[-1])) + args.strip_pitch_mm

    ax_geom.set_xlim(-w_extent, w_extent)
    pad_v = 0.15 * gap_mm
    if v_cathode < 0:
        ax_geom.set_ylim(v_cathode - pad_v, v_wire + pad_v)
    else:
        ax_geom.set_ylim(v_wire - pad_v, v_cathode + pad_v)

    ax_geom.axhline(v_wire, linewidth=2)
    ax_geom.text(
        0.98, v_wire, "anode wire plane",
        ha="right", va="bottom", transform=ax_geom.get_yaxis_transform()
    )

    strip_rects = []
    cathode_thickness = max(0.05 * gap_mm, 0.05)
    for k, center in zip(ks, strip_centers):
        rect = Rectangle(
            (
                center - 0.5 * args.strip_width_mm,
                v_cathode - 0.5 * cathode_thickness,
            ),
            args.strip_width_mm,
            cathode_thickness,
            alpha=0.18,
            linewidth=1.0,
        )
        ax_geom.add_patch(rect)
        strip_rects.append(rect)
        ax_geom.text(
            center,
            v_cathode - np.sign(v_cathode if v_cathode != 0 else -1) * 0.10 * gap_mm,
            f"{int(k):+d}",
            ha="center",
            va="center",
            fontsize=9,
        )

    ax_geom.axhline(v_cathode, linestyle="--", linewidth=1.2)
    ax_geom.plot(wd_mm, vd_mm, linewidth=1.0, alpha=0.25)
    ion_marker, = ax_geom.plot([frame_w[0]], [frame_v[0]], marker="o", markersize=10)
    trail, = ax_geom.plot([], [], linewidth=2.0, alpha=0.6)

    time_text = ax_geom.text(
        0.04,
        0.96,
        "",
        transform=ax_geom.transAxes,
        ha="left",
        va="top",
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.8),
    )

    ax_geom.set_xlabel("w [mm]")
    ax_geom.set_ylabel("v [mm]")
    ax_geom.set_title("Ion position and readout strips")
    ax_geom.grid(alpha=0.25)

    # --------------------------------------------------------------
    # Cumulative induced charge on neighbouring strips.
    # --------------------------------------------------------------
    qbar_values0 = np.array([frame_q[label][0] for label in labels], dtype=float)
    bars = ax_qbars.bar(ks, qbar_values0)

    all_q_values = np.concatenate(
        [np.asarray(frame_q[label], dtype=float) for label in labels]
    )
    qmin = float(np.min(all_q_values))
    qmax = float(np.max(all_q_values))
    qspan = max(qmax - qmin, max(abs(qmin), abs(qmax)) * 0.1, 1e-30)
    ax_qbars.set_ylim(qmin - 0.10 * qspan, qmax + 0.10 * qspan)
    ax_qbars.axhline(0.0, linewidth=0.8)
    ax_qbars.set_xticks(ks)
    ax_qbars.set_xlabel("strip index k")
    ax_qbars.set_ylabel(r"$Q_k(t)$ [fC]")
    ax_qbars.set_title("Cumulative induced charge on neighbouring strips")
    ax_qbars.grid(axis="y", alpha=0.25)

    # --------------------------------------------------------------
    # Neighbour-strip current panel.
    # --------------------------------------------------------------
    current_lines = {}
    current_markers = {}
    all_current_values = []

    for label, k in zip(labels, ks):
        values = wave.loc[physical_mask, f"{label}_fC_per_ns"].to_numpy(dtype=float)
        line, = ax_current.plot(
            tw_plot_us,
            values,
            linewidth=1.4,
            label=f"k={int(k):+d}",
        )
        marker, = ax_current.plot(
            [frame_t_us[0]],
            [frame_currents[label][0]],
            marker="o",
            markersize=5,
            linestyle="none",
        )
        current_lines[label] = line
        current_markers[label] = marker
        all_current_values.append(values)

    all_current_values = np.concatenate(all_current_values)
    imin = float(np.min(all_current_values))
    imax = float(np.max(all_current_values))
    ispan = max(imax - imin, max(abs(imin), abs(imax)) * 0.1, 1e-30)
    ax_current.set_ylim(imin - 0.08 * ispan, imax + 0.08 * ispan)

    current_cursor = ax_current.axvline(frame_t_us[0], linestyle="--", linewidth=1)
    ax_current.set_ylabel("induced current [fC/ns]")
    ax_current.set_title("Currents on central and neighbouring strips")
    ax_current.legend(ncol=3, fontsize=8)
    ax_current.grid(alpha=0.25)

    # --------------------------------------------------------------
    # Cumulative central-strip charge.
    # --------------------------------------------------------------
    q_plot = q0_cum[physical_mask]
    ax_q.plot(tw_plot_us, q_plot, linewidth=1.5)
    q_marker, = ax_q.plot([frame_t_us[0]], [frame_q0[0]], marker="o", markersize=6)
    q_cursor = ax_q.axvline(frame_t_us[0], linestyle="--", linewidth=1)
    ax_q.set_xlabel("physical detector time [us]")
    ax_q.set_ylabel(r"$Q_0(t)$ [fC]")
    ax_q.set_title("Cumulative induced charge on strip 0")
    ax_q.grid(alpha=0.25)

    fig.suptitle(
        "Single positive ion: drift, induced currents and strip charges",
        fontsize=14,
    )
    fig.tight_layout()

    def update(frame):
        t_us = frame_t_us[frame]

        ion_marker.set_data([frame_w[frame]], [frame_v[frame]])
        trail.set_data(frame_w[: frame + 1], frame_v[: frame + 1])

        for bar, label in zip(bars, labels):
            bar.set_height(float(frame_q[label][frame]))

        for label in labels:
            current_markers[label].set_data(
                [t_us], [frame_currents[label][frame]]
            )

        current_cursor.set_xdata([t_us, t_us])
        q_marker.set_data([t_us], [frame_q0[frame]])
        q_cursor.set_xdata([t_us, t_us])

        time_text.set_text(
            f"physical time = {t_us:.2f} us\n"
            f"v = {frame_v[frame]:.3f} mm\n"
            f"Q0 = {frame_q0[frame]:.3e} fC"
        )

        return (
            ion_marker,
            trail,
            *bars,
            *current_markers.values(),
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
        out_gif = Path(str(prefix) + "_strip_sharing.gif")

    animation.save(out_gif, writer=PillowWriter(fps=args.fps))
    plt.close(fig)

    physical_duration_us = (t1_ns - t0_ns) / 1000.0
    slowdown = args.duration_s / (physical_duration_us * 1e-6)

    print(f"Driftline        : {drift_file}")
    print(f"Waveforms        : {wave_file}")
    print(f"Physical drift   : {physical_duration_us:.6f} us")
    print(f"Shown strips     : {int(ks[0])} .. {int(ks[-1])}")
    print(f"GIF duration     : {args.duration_s:.3f} s")
    print(f"Playback slowdown: {slowdown:.3e} x")
    print(f"Wrote            : {out_gif}")


if __name__ == "__main__":
    main()
