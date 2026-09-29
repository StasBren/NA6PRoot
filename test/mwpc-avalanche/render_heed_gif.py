#!/usr/bin/env python3
"""
Render one real Heed + Garfield microscopic muon event as an animated GIF.

Input files are produced by mwpc_heed_visualization:
  <prefix>_meta.txt
  <prefix>_clusters.csv
  <prefix>_electrons.csv
  <prefix>_paths.csv

The animation has three stages:
  1. the muon crosses the gas gap,
  2. Heed ionisation clusters and conduction-electron seeds appear,
  3. actual Garfield microscopic drift trajectories are revealed in time.

Only the u-v projection is drawn. The CSV keeps w as well, so a 3-D or
w-v visualisation can be added later without rerunning Heed/Garfield.
"""

import argparse
import bisect
import csv
import math
import os
import shutil
import tempfile
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from PIL import Image


def read_meta(path):
    meta = {}
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or "=" not in line:
                continue
            key, value = line.split("=", 1)
            meta[key] = value
    return meta


def read_csv(path):
    with open(path, "r", encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))


def f(row, key):
    return float(row[key])


def i(row, key):
    return int(row[key])


def group_paths(rows):
    paths = {}
    for row in rows:
        eid = i(row, "electron")
        paths.setdefault(eid, []).append(
            {
                "u": f(row, "u_mm"),
                "v": f(row, "v_mm"),
                "w": f(row, "w_mm"),
                "t": f(row, "t_ns"),
                "status": i(row, "status"),
                "point": i(row, "point"),
            }
        )
    for pts in paths.values():
        pts.sort(key=lambda p: p["point"])
    return paths


def status_name(status):
    names = {
        -1: "left drift area",
        -3: "error",
        -5: "outside drift medium",
        -7: "attachment",
        -16: "below transport cut",
        -17: "outside time window",
    }
    return names.get(status, "collected / other")


def visible_wire_positions(meta):
    """Return only the wires bracketing the projected muon path in u."""
    pitch = float(meta["wire_pitch_mm"])
    u_start = float(meta.get("muon_u_start_mm", meta["muon_u_mm"]))
    u_end = float(meta.get("muon_u_end_mm", u_start))

    umin = min(u_start, u_end)
    umax = max(u_start, u_end)

    left_index = math.floor(umin / pitch)
    right_index = math.ceil(umax / pitch)

    # A normal-incidence track entirely inside one cell should still show
    # both bounding wires.
    if right_index == left_index:
        right_index += 1

    return [k * pitch for k in range(left_index, right_index + 1)]


def draw_geometry(ax, meta):
    pitch = float(meta["wire_pitch_mm"])
    diameter_mm = float(meta["wire_diameter_um"]) * 1.0e-3
    gap_minus = float(meta["gap_minus_mm"])
    gap_plus = float(meta["gap_plus_mm"])

    # Cathode planes.
    ax.axhline(+gap_plus, lw=3, alpha=0.8)
    ax.axhline(-gap_minus, lw=3, alpha=0.8)

    # Presentation view: show only the two wires relevant to this cell.
    # The Garfield field calculation itself still contains the full 9-wire
    # local array; this is only a visual simplification.
    wires = visible_wire_positions(meta)
    radius = max(0.07, 0.5 * diameter_mm)
    for x in wires:
        ax.add_patch(Circle((x, 0.0), radius, zorder=5))

    ax.text(
        0.01,
        0.97,
        "cathode",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=10,
    )
    ax.text(
        0.01,
        0.05,
        "cathode",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=10,
    )
    ax.text(
        0.985,
        0.51,
        "anode wires",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=10,
    )


def render_frame(
    out_path,
    frame,
    n_muon,
    n_cluster,
    n_drift,
    meta,
    clusters,
    electrons,
    paths,
    tmin,
    tmax,
):
    fig, ax = plt.subplots(figsize=(12.8, 7.2))
    draw_geometry(ax, meta)

    pitch = float(meta["wire_pitch_mm"])
    gap_minus = float(meta["gap_minus_mm"])
    gap_plus = float(meta["gap_plus_mm"])
    u_start = float(meta.get("muon_u_start_mm", meta["muon_u_mm"]))
    u_end = float(meta.get("muon_u_end_mm", u_start))
    v_start = float(meta["muon_v_start_mm"])
    v_end = float(meta["muon_v_end_mm"])

    wires = visible_wire_positions(meta)
    left_wire = min(wires)
    right_wire = max(wires)
    margin = max(0.35 * pitch, 1.0)
    ax.set_xlim(left_wire - margin, right_wire + margin)
    ax.set_ylim(-gap_minus - 0.7, gap_plus + 0.7)
    ax.set_xlabel("u  — across anode wires [mm]")
    ax.set_ylabel("v  — chamber normal [mm]")
    ax.set_title("One muon in the MWPC: Heed ionisation → Garfield electron drift")

    # Stage 1: muon moves through the chamber.
    if frame < n_muon:
        frac = frame / max(1, n_muon - 1)
        v_now = v_start + frac * (v_end - v_start)
        u_now = u_start + frac * (u_end - u_start)
        ax.plot([u_start, u_now], [v_start, v_now], lw=2.5, label="muon track")
        ax.scatter([u_now], [v_now], s=60, marker="v", zorder=8)
        stage = "1  Muon crossing"

    # Stage 2: freeze complete muon track and reveal Heed clusters.
    elif frame < n_muon + n_cluster:
        ax.plot([u_start, u_end], [v_start, v_end], lw=2.0, alpha=0.8, label="muon track")
        frac = (frame - n_muon) / max(1, n_cluster - 1)
        nshow = max(1, int(math.ceil(frac * len(clusters))))

        shown = clusters[:nshow]
        if shown:
            xs = [f(r, "u_mm") for r in shown]
            ys = [f(r, "v_mm") for r in shown]
            sizes = [18.0 + 7.0 * math.sqrt(max(1, i(r, "electrons"))) for r in shown]
            ax.scatter(xs, ys, s=sizes, alpha=0.85, zorder=7, label="Heed clusters")

        # Show only the presentation-selected conduction-electron seeds.
        # This keeps a single large Heed cluster from becoming an opaque
        # orange cloud while preserving the full cluster multiplicity through
        # the blue marker size.
        shown_clusters = {i(r, "cluster") for r in shown}
        eseeds = [
            r for r in electrons
            if i(r, "cluster") in shown_clusters and i(r, "drifted") == 1
        ]
        if eseeds:
            ax.scatter(
                [f(r, "u_mm") for r in eseeds],
                [f(r, "v_mm") for r in eseeds],
                s=10,
                alpha=0.55,
                zorder=6,
                label="selected conduction e⁻ seeds",
            )
        stage = "2  Heed ionisation clusters"

    # Stage 3: reveal real Garfield microscopic trajectories in physical time.
    else:
        ax.plot([u_start, u_end], [v_start, v_end], lw=1.5, alpha=0.35, label="muon track")

        if clusters:
            ax.scatter(
                [f(r, "u_mm") for r in clusters],
                [f(r, "v_mm") for r in clusters],
                s=[14.0 + 5.0 * math.sqrt(max(1, i(r, "electrons"))) for r in clusters],
                alpha=0.45,
                zorder=5,
                label="Heed clusters",
            )

        if electrons:
            visual_seeds = [r for r in electrons if i(r, "drifted") == 1]
            ax.scatter(
                [f(r, "u_mm") for r in visual_seeds],
                [f(r, "v_mm") for r in visual_seeds],
                s=8,
                alpha=0.25,
                zorder=4,
                label="electron start",
            )

        frac = (frame - n_muon - n_cluster) / max(1, n_drift - 1)
        t_now = tmin + frac * (tmax - tmin)

        current_u = []
        current_v = []
        attached_u = []
        attached_v = []

        for eid, pts in paths.items():
            times = [p["t"] for p in pts]
            if t_now < times[0]:
                continue
            j = bisect.bisect_right(times, t_now)
            shown = pts[: max(1, j)]

            ax.plot(
                [p["u"] for p in shown],
                [p["v"] for p in shown],
                lw=0.9,
                alpha=0.55,
                zorder=3,
            )

            pnow = shown[-1]
            if j < len(pts):
                current_u.append(pnow["u"])
                current_v.append(pnow["v"])
            else:
                if pts[-1]["status"] == -7:
                    attached_u.append(pnow["u"])
                    attached_v.append(pnow["v"])
                else:
                    current_u.append(pnow["u"])
                    current_v.append(pnow["v"])

        if current_u:
            ax.scatter(current_u, current_v, s=16, zorder=9, label="drifting e⁻")
        if attached_u:
            ax.scatter(
                attached_u,
                attached_v,
                s=24,
                marker="x",
                zorder=9,
                label="attached e⁻",
            )

        stage = f"3  Garfield microscopic drift   t = {t_now - tmin:.1f} ns"

    # Information panel.
    n_clusters = meta.get("n_clusters", "?")
    n_e = meta.get("n_conduction_electrons", "?")
    n_vis = meta.get("n_visualised_electrons", "?")
    dE = float(meta.get("total_energy_transfer_eV", "nan"))
    info = (
        f"Ar/CO₂ 70/30   |   μ⁻ {float(meta['momentum_GeV']):g} GeV/c\n"
        f"Heed: {n_clusters} clusters, {n_e} conduction e⁻, "
        f"ΣΔE = {dE / 1000.0:.2f} keV\n"
        f"Garfield trajectories shown: {n_vis}\n"
        f"sampling: {meta.get('sample_mode', 'unknown')}   |   "
        f"theta = {float(meta.get('theta_deg', 0.0)):g}°   |   "
        f"displayed wires: " + ", ".join(f"{x:g}" for x in wires) + " mm"
    )
    ax.text(
        0.015,
        0.015,
        info,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=10,
        bbox=dict(boxstyle="round,pad=0.4", alpha=0.9),
    )
    ax.text(
        0.5,
        0.955,
        stage,
        transform=ax.transAxes,
        ha="center",
        va="top",
        fontsize=13,
        fontweight="bold",
    )

    # Deduplicate legend labels.
    handles, labels = ax.get_legend_handles_labels()
    seen = set()
    h2, l2 = [], []
    for h, lab in zip(handles, labels):
        if lab in seen:
            continue
        seen.add(lab)
        h2.append(h)
        l2.append(lab)
    if h2:
        ax.legend(h2, l2, loc="upper right", framealpha=0.9, fontsize=9)

    ax.grid(alpha=0.15)
    fig.tight_layout()
    fig.savefig(out_path, dpi=120)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--prefix", default="heed_vis")
    parser.add_argument("--output", default="heed_muon_drift.gif")
    parser.add_argument("--fps", type=float, default=18.0)
    parser.add_argument("--muon-frames", type=int, default=18)
    parser.add_argument("--cluster-frames", type=int, default=24)
    parser.add_argument("--drift-frames", type=int, default=70)
    parser.add_argument("--hold-frames", type=int, default=18)
    parser.add_argument("--keep-frames", action="store_true")
    args = parser.parse_args()

    prefix = Path(args.prefix)
    meta = read_meta(str(prefix) + "_meta.txt")
    clusters = read_csv(str(prefix) + "_clusters.csv")
    electrons = read_csv(str(prefix) + "_electrons.csv")
    path_rows = read_csv(str(prefix) + "_paths.csv")
    paths = group_paths(path_rows)

    if not paths:
        raise RuntimeError(
            "No Garfield drift paths were found. "
            "Check that the C++ event dump completed successfully."
        )

    all_times = [p["t"] for pts in paths.values() for p in pts]
    tmin = min(all_times)
    tmax = max(all_times)

    total_frames = args.muon_frames + args.cluster_frames + args.drift_frames

    if args.keep_frames:
        frame_dir = Path("heed_gif_frames")
        frame_dir.mkdir(parents=True, exist_ok=True)
        cleanup = False
    else:
        frame_dir = Path(tempfile.mkdtemp(prefix="heed_gif_frames_"))
        cleanup = True

    frame_paths = []
    try:
        for frame in range(total_frames):
            out = frame_dir / f"frame_{frame:04d}.png"
            render_frame(
                out,
                frame,
                args.muon_frames,
                args.cluster_frames,
                args.drift_frames,
                meta,
                clusters,
                electrons,
                paths,
                tmin,
                tmax,
            )
            frame_paths.append(out)

        # Hold the final frame for readability in the presentation.
        for j in range(args.hold_frames):
            out = frame_dir / f"frame_{total_frames + j:04d}.png"
            shutil.copyfile(frame_paths[-1], out)
            frame_paths.append(out)

        images = [Image.open(p).convert("RGB") for p in frame_paths]
        duration_ms = int(round(1000.0 / args.fps))
        images[0].save(
            args.output,
            save_all=True,
            append_images=images[1:],
            duration=duration_ms,
            loop=0,
            optimize=False,
        )
        for img in images:
            img.close()

        print(f"Wrote {args.output}")
        print(
            f"Frames: {len(frame_paths)}, fps: {args.fps:g}, "
            f"duration: {len(frame_paths) / args.fps:.1f} s"
        )
        if args.keep_frames:
            print(f"PNG frames kept in {frame_dir}/")
    finally:
        if cleanup:
            shutil.rmtree(frame_dir, ignore_errors=True)


if __name__ == "__main__":
    main()
