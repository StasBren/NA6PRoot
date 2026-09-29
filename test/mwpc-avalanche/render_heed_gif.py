#!/usr/bin/env python3
"""
Render one real Heed + Garfield microscopic muon event as an animated GIF.

Inputs are produced by mwpc_heed_visualization:
  <prefix>_meta.txt
  <prefix>_clusters.csv
  <prefix>_electrons.csv
  <prefix>_paths.csv
  <prefix>_field.csv

Stages:
  1. muon crossing,
  2. Heed ionisation clusters / selected conduction-electron seeds,
  3. actual Garfield microscopic electron drift trajectories.

The optional equipotential contours are sampled from the same
ComponentAnalyticField used to transport the electrons.
"""

import argparse
import bisect
import csv
import math
import shutil
import tempfile
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import numpy as np
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


def load_field_grid(rows):
    us = np.array(sorted({float(r["u_mm"]) for r in rows}))
    vs = np.array(sorted({float(r["v_mm"]) for r in rows}))
    nu = len(us)
    nv = len(vs)

    potential = np.full((nv, nu), np.nan)
    status = np.full((nv, nu), -999, dtype=int)

    u_index = {x: j for j, x in enumerate(us)}
    v_index = {x: j for j, x in enumerate(vs)}

    for row in rows:
        iu = u_index[float(row["u_mm"])]
        iv = v_index[float(row["v_mm"])]
        potential[iv, iu] = float(row["potential_V"])
        status[iv, iu] = int(row["status"])

    return us, vs, potential, status


def default_wire_positions(meta):
    """Wires bracketing the projected muon path in u."""
    pitch = float(meta["wire_pitch_mm"])
    u_start = float(meta.get("muon_u_start_mm", meta["muon_u_mm"]))
    u_end = float(meta.get("muon_u_end_mm", u_start))

    umin = min(u_start, u_end)
    umax = max(u_start, u_end)

    left_index = math.floor(umin / pitch)
    right_index = math.ceil(umax / pitch)
    if right_index == left_index:
        right_index += 1

    return [k * pitch for k in range(left_index, right_index + 1)]


def visible_wire_positions(meta, xmin, xmax):
    pitch = float(meta["wire_pitch_mm"])
    kmin = math.floor(xmin / pitch) - 1
    kmax = math.ceil(xmax / pitch) + 1
    return [
        k * pitch
        for k in range(kmin, kmax + 1)
        if xmin - 0.15 * pitch <= k * pitch <= xmax + 0.15 * pitch
    ]


def draw_contours(ax, meta, field, xmin, xmax, ymin, ymax):
    if field is None:
        return

    us, vs, potential, status = field
    U, V = np.meshgrid(us, vs)

    mask = (
        (U < xmin)
        | (U > xmax)
        | (V < ymin)
        | (V > ymax)
        | (status != 0)
        | ~np.isfinite(potential)
    )
    phi = np.ma.masked_where(mask, potential)

    hv = float(meta["anode_voltage_V"])
    levels = (
        np.array([0.02, 0.05, 0.10, 0.20, 0.35,
                  0.50, 0.65, 0.80, 0.90, 0.96])
        * hv
    )
    try:
        ax.contour(
            U,
            V,
            phi,
            levels=levels,
            linewidths=0.75,
            alpha=0.42,
            zorder=1,
        )
    except ValueError:
        pass


def draw_geometry(ax, meta, xmin, xmax):
    pitch = float(meta["wire_pitch_mm"])
    diameter_mm = float(meta["wire_diameter_um"]) * 1.0e-3
    gap_minus = float(meta["gap_minus_mm"])
    gap_plus = float(meta["gap_plus_mm"])

    ax.axhline(+gap_plus, lw=3, alpha=0.8, zorder=2)
    ax.axhline(-gap_minus, lw=3, alpha=0.8, zorder=2)

    wires = visible_wire_positions(meta, xmin, xmax)
    radius = max(0.07, 0.5 * diameter_mm)
    for x in wires:
        ax.add_patch(Circle((x, 0.0), radius, zorder=5))

    ax.text(
        0.01, 0.97, "cathode",
        transform=ax.transAxes,
        ha="left", va="top", fontsize=10,
    )
    ax.text(
        0.01, 0.05, "cathode",
        transform=ax.transAxes,
        ha="left", va="bottom", fontsize=10,
    )
    ax.text(
        0.985, 0.51, "anode wires",
        transform=ax.transAxes,
        ha="right", va="bottom", fontsize=10,
    )
    return wires


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
    field,
    background,
    u_min_mm,
    u_max_mm,
):
    fig, ax = plt.subplots(figsize=(12.8, 7.2))

    pitch = float(meta["wire_pitch_mm"])
    gap_minus = float(meta["gap_minus_mm"])
    gap_plus = float(meta["gap_plus_mm"])
    u_start = float(meta.get("muon_u_start_mm", meta["muon_u_mm"]))
    u_end = float(meta.get("muon_u_end_mm", u_start))
    v_start = float(meta["muon_v_start_mm"])
    v_end = float(meta["muon_v_end_mm"])

    defaults = default_wire_positions(meta)
    margin = max(0.35 * pitch, 1.0)
    xmin = min(defaults) - margin if u_min_mm is None else u_min_mm
    xmax = max(defaults) + margin if u_max_mm is None else u_max_mm

    if xmin >= xmax:
        raise RuntimeError("--u-min-mm must be smaller than --u-max-mm.")

    ymin = -gap_minus - 0.7
    ymax = gap_plus + 0.7

    if background == "contours":
        draw_contours(ax, meta, field, xmin, xmax, -gap_minus, gap_plus)

    wires = draw_geometry(ax, meta, xmin, xmax)

    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_xlabel("u — across anode wires [mm]")
    ax.set_ylabel("v — chamber normal [mm]")
    ax.set_title("One muon in the MWPC: Heed ionisation → Garfield electron drift")

    if background == "contours":
        ax.text(
            0.985, 0.985,
            "Garfield equipotential contours",
            transform=ax.transAxes,
            ha="right", va="top",
            fontsize=9, alpha=0.8,
        )

    if frame < n_muon:
        frac = frame / max(1, n_muon - 1)
        v_now = v_start + frac * (v_end - v_start)
        u_now = u_start + frac * (u_end - u_start)
        ax.plot(
            [u_start, u_now], [v_start, v_now],
            lw=2.5, label="muon track",
        )
        ax.scatter([u_now], [v_now], s=60, marker="v", zorder=8)
        stage = "1  Muon crossing"

    elif frame < n_muon + n_cluster:
        ax.plot(
            [u_start, u_end], [v_start, v_end],
            lw=2.0, alpha=0.8, label="muon track",
        )
        frac = (frame - n_muon) / max(1, n_cluster - 1)
        nshow = max(1, int(math.ceil(frac * len(clusters))))

        shown = clusters[:nshow]
        if shown:
            xs = [f(r, "u_mm") for r in shown]
            ys = [f(r, "v_mm") for r in shown]
            sizes = [
                18.0 + 7.0 * math.sqrt(max(1, i(r, "electrons")))
                for r in shown
            ]
            ax.scatter(
                xs, ys,
                s=sizes, alpha=0.85, zorder=7,
                label="Heed clusters",
            )

        shown_clusters = {i(r, "cluster") for r in shown}
        eseeds = [
            r for r in electrons
            if i(r, "cluster") in shown_clusters
            and i(r, "drifted") == 1
        ]
        if eseeds:
            ax.scatter(
                [f(r, "u_mm") for r in eseeds],
                [f(r, "v_mm") for r in eseeds],
                s=10, alpha=0.55, zorder=6,
                label="selected conduction e⁻ seeds",
            )
        stage = "2  Heed ionisation clusters"

    else:
        ax.plot(
            [u_start, u_end], [v_start, v_end],
            lw=1.5, alpha=0.35, label="muon track",
        )

        if clusters:
            ax.scatter(
                [f(r, "u_mm") for r in clusters],
                [f(r, "v_mm") for r in clusters],
                s=[
                    14.0 + 5.0 * math.sqrt(max(1, i(r, "electrons")))
                    for r in clusters
                ],
                alpha=0.45, zorder=5,
                label="Heed clusters",
            )

        visual_seeds = [r for r in electrons if i(r, "drifted") == 1]
        if visual_seeds:
            ax.scatter(
                [f(r, "u_mm") for r in visual_seeds],
                [f(r, "v_mm") for r in visual_seeds],
                s=8, alpha=0.25, zorder=4,
                label="electron start",
            )

        frac = (frame - n_muon - n_cluster) / max(1, n_drift - 1)
        t_now = tmin + frac * (tmax - tmin)

        current_u = []
        current_v = []
        attached_u = []
        attached_v = []

        for pts in paths.values():
            times = [p["t"] for p in pts]
            if t_now < times[0]:
                continue

            j = bisect.bisect_right(times, t_now)
            shown = pts[:max(1, j)]

            ax.plot(
                [p["u"] for p in shown],
                [p["v"] for p in shown],
                lw=0.9, alpha=0.55, zorder=3,
            )

            pnow = shown[-1]
            if j < len(pts):
                current_u.append(pnow["u"])
                current_v.append(pnow["v"])
            elif pts[-1]["status"] == -7:
                attached_u.append(pnow["u"])
                attached_v.append(pnow["v"])
            else:
                current_u.append(pnow["u"])
                current_v.append(pnow["v"])

        if current_u:
            ax.scatter(
                current_u, current_v,
                s=16, zorder=9, label="drifting e⁻",
            )
        if attached_u:
            ax.scatter(
                attached_u, attached_v,
                s=24, marker="x", zorder=9,
                label="attached e⁻",
            )

        stage = (
            f"3  Garfield microscopic drift   "
            f"t = {t_now - tmin:.1f} ns"
        )

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
        0.015, 0.015, info,
        transform=ax.transAxes,
        ha="left", va="bottom",
        fontsize=10,
        bbox=dict(boxstyle="round,pad=0.4", alpha=0.9),
    )
    ax.text(
        0.5, 0.955, stage,
        transform=ax.transAxes,
        ha="center", va="top",
        fontsize=13, fontweight="bold",
    )

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
    parser.add_argument(
        "--u-min-mm", type=float, default=None,
        help="Explicit minimum displayed local-u coordinate [mm]",
    )
    parser.add_argument(
        "--u-max-mm", type=float, default=None,
        help="Explicit maximum displayed local-u coordinate [mm]",
    )
    parser.add_argument(
        "--background",
        choices=["none", "contours"],
        default="none",
        help="Optional Garfield electrostatic background",
    )
    parser.add_argument("--keep-frames", action="store_true")
    args = parser.parse_args()

    prefix = Path(args.prefix)
    meta = read_meta(str(prefix) + "_meta.txt")
    clusters = read_csv(str(prefix) + "_clusters.csv")
    electrons = read_csv(str(prefix) + "_electrons.csv")
    path_rows = read_csv(str(prefix) + "_paths.csv")
    paths = group_paths(path_rows)

    field = None
    field_path = Path(str(prefix) + "_field.csv")
    if field_path.exists():
        field = load_field_grid(read_csv(field_path))
    elif args.background == "contours":
        raise RuntimeError(
            f"{field_path} is missing. Rerun mwpc_heed_visualization "
            "after pulling/rebuilding the updated C++ code."
        )

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
                field,
                args.background,
                args.u_min_mm,
                args.u_max_mm,
            )
            frame_paths.append(out)

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
