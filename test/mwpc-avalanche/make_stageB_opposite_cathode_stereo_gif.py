#!/usr/bin/env python3
"""Animate actual Garfield electron/ion trajectories in the B3c.2 geometry."""

import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--trajectories", required=True)
    p.add_argument("--output", default="stageB_opposite_cathode_stereo.gif")
    p.add_argument("--gap-minus-mm", type=float, default=2.0)
    p.add_argument("--gap-plus-mm", type=float, default=4.0)
    p.add_argument("--wire-pitch-mm", type=float, default=4.0)
    p.add_argument("--strip-pitch-mm", type=float, default=1.7)
    p.add_argument("--strip-width-mm", type=float, default=1.7)
    p.add_argument("--tan-alpha", type=float, default=0.10)
    p.add_argument("--frames", type=int, default=90)
    p.add_argument("--fps", type=int, default=12)
    return p.parse_args()


def sample_track_at_time(track, t):
    tt = track["t_ns"].to_numpy()
    if t < tt[0]:
        return None
    if t >= tt[-1]:
        row = track.iloc[-1]
        return row["u_mm"], row["v_mm"], row["w_mm"]

    j = np.searchsorted(tt, t)
    a = track.iloc[j - 1]
    b = track.iloc[j]
    f = (t - a["t_ns"]) / max(b["t_ns"] - a["t_ns"], 1e-12)
    return (
        a["u_mm"] + f * (b["u_mm"] - a["u_mm"]),
        a["v_mm"] + f * (b["v_mm"] - a["v_mm"]),
        a["w_mm"] + f * (b["w_mm"] - a["w_mm"]),
    )


def main():
    args = parse_args()
    df = pd.read_csv(args.trajectories)
    if df.empty:
        raise RuntimeError("Trajectory CSV is empty.")

    alpha = np.arctan(args.tan_alpha)

    e_tracks = {
        tid: g.sort_values("t_ns")
        for tid, g in df[df["carrier"] == "electron"].groupby("track_id")
    }
    i_tracks = {
        tid: g.sort_values("t_ns")
        for tid, g in df[df["carrier"] == "ion"].groupby("track_id")
    }

    t_e_max = (
        max(g["t_ns"].max() for g in e_tracks.values())
        if e_tracks else 0.0
    )
    t_i_max = (
        max(g["t_ns"].max() for g in i_tracks.values())
        if i_tracks else 0.0
    )
    t_max = max(t_e_max, t_i_max)

    # Use a nonlinear time grid: resolve the ns electron burst and still show
    # the much slower ion drift to the cathodes.
    n_fast = max(15, args.frames // 3)
    n_slow = args.frames - n_fast
    fast_end = min(max(20.0, 1.2 * t_e_max), t_max)
    t_fast = np.linspace(0., fast_end, n_fast, endpoint=False)
    if n_slow > 0 and t_max > fast_end:
        t_slow = np.geomspace(max(fast_end, 1e-3), t_max, n_slow)
        times = np.r_[t_fast, t_slow]
    else:
        times = np.linspace(0., t_max, args.frames)

    fig = plt.figure(figsize=(12.5, 6.2))
    gs = fig.add_gridspec(1, 2, width_ratios=[0.95, 1.15], wspace=0.24)
    ax_side = fig.add_subplot(gs[0, 0])
    ax_top = fig.add_subplot(gs[0, 1])

    # Side view: u-v.
    ax_side.axhline(-args.gap_minus_mm, lw=2)
    ax_side.axhline(0., lw=1.5)
    ax_side.axhline(+args.gap_plus_mm, lw=2)
    ax_side.set_xlim(-5., 5.)
    ax_side.set_ylim(-args.gap_minus_mm - 0.5, args.gap_plus_mm + 0.5)
    ax_side.set_xlabel("u across wires [mm]")
    ax_side.set_ylabel("v chamber normal [mm]")
    ax_side.set_title("Actual Garfield trajectories: side view")

    # Anode wires seen in u-v cross-section.
    wire_us = np.arange(-4., 4.01, args.wire_pitch_mm)
    ax_side.scatter(wire_us, np.zeros_like(wire_us), s=45, zorder=4)

    side_e = ax_side.scatter([], [], s=8, label="electrons")
    side_i = ax_side.scatter([], [], s=8, label="positive ions")
    ax_side.legend(loc="upper right")

    # Top view: u-w; show both stereo families on opposite cathodes.
    ax_top.set_xlim(-9., 9.)
    ax_top.set_ylim(-5., 5.)
    ax_top.set_xlabel("u across wires [mm]")
    ax_top.set_ylabel("w along wires [mm]")
    ax_top.set_title("Top view: +alpha on v=-2 mm, -alpha on v=+4 mm")

    ugrid = np.linspace(-9., 9., 400)
    for k in range(-3, 4):
        # +alpha family natural coordinate:
        # x+ = w cos(a) - u sin(a) = k p
        w_plus = (k * args.strip_pitch_mm + ugrid * np.sin(alpha)) / np.cos(alpha)
        # -alpha:
        # x- = w cos(a) + u sin(a) = k p
        w_minus = (k * args.strip_pitch_mm - ugrid * np.sin(alpha)) / np.cos(alpha)
        ax_top.plot(ugrid, w_plus, lw=1.2, alpha=0.55)
        ax_top.plot(ugrid, w_minus, lw=1.2, ls="--", alpha=0.55)

    # Anode wires run along w at fixed u.
    for u in np.arange(-8., 8.01, args.wire_pitch_mm):
        ax_top.axvline(u, lw=0.8, alpha=0.25)

    top_e = ax_top.scatter([], [], s=8)
    top_i = ax_top.scatter([], [], s=8)

    time_text = fig.text(
        0.5, 0.97, "", ha="center", va="top", fontsize=13
    )

    def update(frame):
        t = float(times[frame])

        epos = []
        for g in e_tracks.values():
            p = sample_track_at_time(g, t)
            if p is not None:
                epos.append(p)

        ipos = []
        for g in i_tracks.values():
            p = sample_track_at_time(g, t)
            if p is not None:
                ipos.append(p)

        epos = np.asarray(epos) if epos else np.empty((0, 3))
        ipos = np.asarray(ipos) if ipos else np.empty((0, 3))

        side_e.set_offsets(
            np.c_[epos[:, 0], epos[:, 1]]
            if len(epos) else np.empty((0, 2))
        )
        side_i.set_offsets(
            np.c_[ipos[:, 0], ipos[:, 1]]
            if len(ipos) else np.empty((0, 2))
        )
        top_e.set_offsets(
            np.c_[epos[:, 0], epos[:, 2]]
            if len(epos) else np.empty((0, 2))
        )
        top_i.set_offsets(
            np.c_[ipos[:, 0], ipos[:, 2]]
            if len(ipos) else np.empty((0, 2))
        )

        if t < 1000.:
            time_text.set_text(f"physical time = {t:.2f} ns")
        else:
            time_text.set_text(f"physical time = {t/1000.:.2f} us")

        return side_e, side_i, top_e, top_i, time_text

    anim = FuncAnimation(
        fig, update, frames=len(times),
        interval=1000. / args.fps, blit=False
    )
    anim.save(args.output, writer=PillowWriter(fps=args.fps))
    plt.close(fig)
    print(f"Saved {args.output}")


if __name__ == "__main__":
    main()
