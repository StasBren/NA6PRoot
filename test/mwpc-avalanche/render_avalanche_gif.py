#!/usr/bin/env python3
"""
Render one real Garfield++ microscopic avalanche as a GIF.

Inputs are produced by mwpc_phase_a_avalanche_visualization:
  <prefix>_meta.txt
  <prefix>_field.csv       actual ComponentAnalyticField E and potential
  <prefix>_endpoints.csv   all avalanche electron start/end points
  <prefix>_paths.csv       stored microscopic drift-line points + Garfield time

The animation uses the actual Garfield time stored along each drift line.
The background can show equipotential contours, electric-field arrows, both,
or neither.
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
    out = {}
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or "=" not in line:
                continue
            key, value = line.split("=", 1)
            out[key] = value
    return out


def read_csv(path):
    with open(path, "r", encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))


def group_paths(rows):
    paths = {}
    for row in rows:
        eid = int(row["electron"])
        paths.setdefault(eid, []).append(
            {
                "point": int(row["point"]),
                "u": float(row["u_mm"]),
                "v": float(row["v_mm"]),
                "w": float(row["w_mm"]),
                "t": float(row["t_ns"]),
                "status": int(row["status"]),
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
    ex = np.full((nv, nu), np.nan)
    ey = np.full((nv, nu), np.nan)
    status = np.full((nv, nu), -999, dtype=int)

    u_index = {x: i for i, x in enumerate(us)}
    v_index = {x: i for i, x in enumerate(vs)}

    for r in rows:
        iu = u_index[float(r["u_mm"])]
        iv = v_index[float(r["v_mm"])]
        potential[iv, iu] = float(r["potential_V"])
        ex[iv, iu] = float(r["ex_Vcm"])
        ey[iv, iu] = float(r["ey_Vcm"])
        status[iv, iu] = int(r["status"])

    return us, vs, potential, ex, ey, status


def nearest_wire(seed_u, pitch):
    return round(seed_u / pitch) * pitch


def clip_paths_near_wire(paths, wire_u, radius_mm):
    """Keep only the part of each trajectory inside a circle around the wire.

    One point immediately before first entry is retained when available so the
    approach to the avalanche region remains visually continuous. Trajectories
    that never enter the requested radius are omitted entirely.
    """
    if radius_mm is None:
        return paths

    clipped = {}
    r2max = radius_mm * radius_mm

    for eid, pts in paths.items():
        inside_indices = [
            j for j, p in enumerate(pts)
            if (p["u"] - wire_u) ** 2 + p["v"] ** 2 <= r2max
        ]
        if not inside_indices:
            continue

        first = inside_indices[0]
        last = inside_indices[-1]

        # Keep one point just before entry to make the incoming segment visible.
        first_keep = max(0, first - 1)
        clipped[eid] = pts[first_keep:last + 1]

    return clipped


def draw_background(ax, meta, field, background, half_width, v_half_window):
    us, vs, potential, ex, ey, status = field

    seed_u = float(meta["seed_u_mm"])
    pitch = float(meta["wire_pitch_mm"])
    wire_u = nearest_wire(seed_u, pitch)
    hv = float(meta["anode_voltage_V"])

    xmin = wire_u - half_width
    xmax = wire_u + half_width

    if v_half_window is None:
        ymin = max(float(vs.min()), -half_width)
        ymax = min(float(vs.max()), +half_width)
    else:
        ymin = max(float(vs.min()), -v_half_window)
        ymax = min(float(vs.max()), +v_half_window)

    U, V = np.meshgrid(us, vs)

    inside = (
        (U >= xmin) & (U <= xmax) &
        (V >= ymin) & (V <= ymax)
    )

    good = (status == 0) & inside & np.isfinite(potential)
    phi = np.ma.masked_where(~good, potential)

    if background in ("contours", "both"):
        fracs = np.array([0.02, 0.05, 0.10, 0.20, 0.35,
                          0.50, 0.65, 0.80, 0.90, 0.96])
        levels = fracs * hv
        try:
            ax.contour(U, V, phi, levels=levels,
                       linewidths=0.75, alpha=0.60)
        except ValueError:
            pass

    if background in ("field", "both"):
        # Show field direction, not magnitude; otherwise the near-wire values
        # visually dominate everything.
        step_u = max(1, len(us) // 24)
        step_v = max(1, len(vs) // 18)
        uu = U[::step_v, ::step_u]
        vv = V[::step_v, ::step_u]
        qx = ex[::step_v, ::step_u]
        qy = ey[::step_v, ::step_u]
        qs = status[::step_v, ::step_u]

        mag = np.hypot(qx, qy)
        mask = (
            (uu >= xmin) & (uu <= xmax) &
            (vv >= ymin) & (vv <= ymax) &
            (qs == 0) & np.isfinite(mag) & (mag > 0)
        )
        qx = np.where(mask, qx / np.where(mag > 0, mag, 1), np.nan)
        qy = np.where(mask, qy / np.where(mag > 0, mag, 1), np.nan)
        ax.quiver(uu, vv, qx, qy, angles="xy",
                  scale_units="xy", scale=7.5,
                  width=0.0022, alpha=0.35)

    # Cathodes are usually outside this avalanche zoom, but draw them if visible.
    gap_minus = float(meta["gap_minus_mm"])
    gap_plus = float(meta["gap_plus_mm"])
    if ymin <= -gap_minus <= ymax:
        ax.axhline(-gap_minus, linewidth=2.0, alpha=0.8)
    if ymin <= +gap_plus <= ymax:
        ax.axhline(+gap_plus, linewidth=2.0, alpha=0.8)

    # Relevant anode wire(s). Usually only the nearest wire is inside the zoom.
    wire_radius = 0.5 * float(meta["wire_diameter_um"]) * 1.e-3
    for k in range(-2, 3):
        xw = wire_u + k * pitch
        if xmin <= xw <= xmax:
            ax.add_patch(Circle(
                (xw, 0.0),
                wire_radius,
                fill=False,
                linewidth=2.2,
                zorder=20,
            ))

    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    return wire_u, xmin, xmax, ymin, ymax


def render_frame(out_path, tcut, meta, field, paths,
                 background, half_width, v_half_window,
                 path_radius_mm):
    fig, ax = plt.subplots(figsize=(10.0, 6.0))

    wire_u, xmin, xmax, ymin, ymax = draw_background(
        ax, meta, field, background, half_width, v_half_window
    )

    seed_u = float(meta["seed_u_mm"])
    seed_v = float(meta["seed_v_mm"])

    # For an ultra-tight avalanche view the original seed can be far outside
    # the plotted region. Only draw it if it actually lies in the frame.
    if xmin <= seed_u <= xmax and ymin <= seed_v <= ymax:
        ax.scatter([seed_u], [seed_v], s=35, marker="x",
                   zorder=25, label="initial electron")

    current_u = []
    current_v = []
    finished_u = []
    finished_v = []
    attached_u = []
    attached_v = []
    born = 0

    for eid, pts in paths.items():
        times = [p["t"] for p in pts]
        if not times or tcut < times[0]:
            continue

        born += 1
        j = bisect.bisect_right(times, tcut)
        shown = pts[:max(1, j)]

        ax.plot(
            [p["u"] for p in shown],
            [p["v"] for p in shown],
            linewidth=0.75,
            alpha=0.35,
            zorder=5,
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
                finished_u.append(pnow["u"])
                finished_v.append(pnow["v"])

    if current_u:
        ax.scatter(current_u, current_v, s=12,
                   zorder=15, label="drifting avalanche e⁻")
    if finished_u:
        ax.scatter(finished_u, finished_v, s=7,
                   alpha=0.45, zorder=10, label="finished e⁻ trajectories")
    if attached_u:
        ax.scatter(attached_u, attached_v, s=18, marker="x",
                   alpha=0.65, zorder=14, label="attached e⁻")

    total_e = int(float(meta["avalanche_electrons"]))
    stored = int(float(meta["stored_drift_lines"]))
    t0 = float(meta["t_min_ns"])
    dt = tcut - t0

    ax.set_xlabel("u — across wires [mm]")
    ax.set_ylabel("v — chamber normal [mm]")
    ax.set_title("Microscopic avalanche near an MWPC anode wire")

    info = (
        f"Garfield++ AvalancheMicroscopic\n"
        f"Ar/CO₂ 70/30, V_anode = {float(meta['anode_voltage_V']):g} V\n"
        f"avalanche size = {total_e} e⁻\n"
        f"stored trajectories = {stored}\n"
        f"branches started by frame = {born}\n"
        f"Δt = {dt:.3f} ns"
        + (f"\nshown path radius ≤ {path_radius_mm:g} mm"
           if path_radius_mm is not None else "")
    )
    ax.text(
        0.015, 0.02, info,
        transform=ax.transAxes,
        ha="left", va="bottom",
        fontsize=9.5,
        bbox=dict(boxstyle="round,pad=0.35", alpha=0.92),
        zorder=30,
    )

    bg_label = {
        "contours": "Garfield equipotential contours",
        "field": "Garfield electric-field direction",
        "both": "Garfield potential + field direction",
        "none": "no field overlay",
    }[background]
    ax.text(
        0.985, 0.985, bg_label,
        transform=ax.transAxes,
        ha="right", va="top",
        fontsize=9,
        alpha=0.8,
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
        ax.legend(h2, l2, loc="upper left",
                  bbox_to_anchor=(1.01, 1.0),
                  frameon=False, fontsize=8.5)

    ax.grid(alpha=0.12)
    fig.tight_layout()
    fig.savefig(out_path, dpi=130, bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--prefix", default="avalanche_vis")
    parser.add_argument("--output", default="avalanche.gif")
    parser.add_argument(
        "--background",
        choices=["contours", "field", "both", "none"],
        default="contours",
    )
    parser.add_argument("--half-width-mm", type=float, default=1.25)
    parser.add_argument("--v-half-window-mm", type=float, default=None)
    parser.add_argument(
        "--path-radius-mm", type=float, default=None,
        help="Show only trajectory points within this radius of the nearest anode wire")
    parser.add_argument(
        "--late-time-power", type=float, default=1.7,
        help=">1 concentrates more animation frames near the late avalanche")
    parser.add_argument("--frames", type=int, default=80)
    parser.add_argument("--hold-frames", type=int, default=18)
    parser.add_argument("--fps", type=float, default=18.0)
    parser.add_argument("--keep-frames", action="store_true")
    args = parser.parse_args()

    prefix = Path(args.prefix)
    meta = read_meta(str(prefix) + "_meta.txt")
    field_rows = read_csv(str(prefix) + "_field.csv")
    path_rows = read_csv(str(prefix) + "_paths.csv")

    field = load_field_grid(field_rows)
    paths = group_paths(path_rows)

    seed_u = float(meta["seed_u_mm"])
    pitch = float(meta["wire_pitch_mm"])
    wire_u = nearest_wire(seed_u, pitch)

    # For the close-up avalanche movie, discard the long pre-avalanche drift
    # and keep only trajectory portions close to the selected anode wire.
    paths = clip_paths_near_wire(paths, wire_u, args.path_radius_mm)

    if not paths:
        raise RuntimeError(
            "No microscopic avalanche paths remain in the requested view. "
            "Increase --path-radius-mm or check the C++ run output."
        )

    all_times = [p["t"] for pts in paths.values() for p in pts]
    tmin = min(all_times)
    tmax = max(all_times)

    # Start the movie only when a stored trajectory reaches the requested
    # near-wire region. Use a nonlinear time scan with dense sampling near the
    # late avalanche, where multiplication becomes visually rapid.
    x = np.linspace(0.0, 1.0, args.frames)
    power = max(1.0, args.late_time_power)
    fractions = 1.0 - (1.0 - x) ** power
    times = tmin + fractions * (tmax - tmin)

    if args.keep_frames:
        frame_dir = Path("avalanche_gif_frames")
        frame_dir.mkdir(parents=True, exist_ok=True)
        cleanup = False
    else:
        frame_dir = Path(tempfile.mkdtemp(prefix="avalanche_gif_frames_"))
        cleanup = True

    frame_paths = []
    try:
        for iframe, tcut in enumerate(times):
            out = frame_dir / f"frame_{iframe:04d}.png"
            render_frame(
                out,
                float(tcut),
                meta,
                field,
                paths,
                args.background,
                args.half_width_mm,
                args.v_half_window_mm,
                args.path_radius_mm,
            )
            frame_paths.append(out)

        for j in range(args.hold_frames):
            out = frame_dir / f"frame_{len(frame_paths):04d}.png"
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
            f"frames={len(frame_paths)}, fps={args.fps:g}, "
            f"duration={len(frame_paths) / args.fps:.1f} s"
        )
        if args.keep_frames:
            print(f"PNG frames kept in {frame_dir}/")
    finally:
        if cleanup:
            shutil.rmtree(frame_dir, ignore_errors=True)


if __name__ == "__main__":
    main()
