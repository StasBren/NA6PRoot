#!/usr/bin/env python3
"""PNG-only PI atlas of fixed-grid anode collection in checkerboard MWPCs.

Products: 18 station overviews (3 conditions in physical XY), 894 full-resolution
chamber profiles (3 conditions), 18 station contact sheets and nine diagnostic
figures (multiplicity, |theta_u|, phase robustness).

No new simulation. All displayed wire activity is relative to ALL gas-volume
crossings of its own module, including overlapping checkerboard chamber hits.
"""
from __future__ import annotations
import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.collections import LineCollection
from matplotlib.patches import Rectangle
from matplotlib.cm import ScalarMappable

from mwpc_checkerboard_projection import (
    SOURCES, NX, NY, GAS_X_CM, GAS_Y_CM, MS0_GAS_Y_CM,
    build_chambers,
)
from pi_fixed_wire_occupancy import CONDITIONS, METRICS, PHASES, positions, n_wires

LABELS = {"jpsi40-bore": "J/psi", "phi40-bore": "phi", "omega40-bore": "omega"}
CAPTIONS = ("Wire >=1 collected e-", "Wire >=5 collected e-",
            "Wire >=1 collected e- by 100 ns")
GEOM = build_chambers()
BY_DET = {ch.det_id: ch for station in GEOM for ch in station}


def load_wire_csv(path):
    arrays = np.full((3, 298, 128, 3), np.nan, dtype=float)
    counts = np.zeros((3, 298, 3), dtype=np.int64)
    exp = np.full((3, 298), -1, dtype=int)
    seen = np.zeros((3, 298, 128), dtype=bool)
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        for row in reader:
            src = SOURCES.index(row["source"])
            d, j = int(row["det_id"]), int(row["wire_index"])
            station = int(row["station"])
            if d not in BY_DET or BY_DET[d].station != station or not (0 <= j < n_wires(station)):
                raise ValueError(f"Bad physical wire identity DET{d} wire{j}")
            expected_y = BY_DET[d].y + positions(station)[j] / 10.
            if abs(float(row["wire_y_global_cm"]) - expected_y) > 1e-6:
                raise ValueError(f"Wire global Y does not match geometry: DET{d} wire{j}")
            n = int(row["chamber_crossings"])
            if seen[src, d, j] or (exp[src, d] >= 0 and exp[src, d] != n):
                raise ValueError("Duplicate wire or inconsistent chamber denominator")
            seen[src, d, j] = True
            exp[src, d] = n
            for k, col in enumerate(("n_ge1e", "n_ge5e", "n_le100ns")):
                c = int(row[col])
                if c < 0 or c > n:
                    raise ValueError(f"Invalid per-wire count DET{d}, {col}: {c}/{n}")
                counts[src, d, k] += c
                arrays[src, d, j, k] = c / n if n else np.nan
    for src in range(3):
        for d, ch in BY_DET.items():
            if exp[src, d] < 0 or not np.all(seen[src, d, :n_wires(ch.station)]):
                raise ValueError(f"Missing anodes source={src} det={d}")
    return arrays, exp, counts


def read_csv(path):
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def color_norm(per_channel_values):
    finite = per_channel_values[np.isfinite(per_channel_values)]
    positive = finite[finite > 0]
    if positive.size == 0:
        return Normalize(0., 1.), False
    vmax = float(np.quantile(positive, .995))
    if vmax <= 0:
        vmax = float(np.max(positive))
    return Normalize(vmin=0., vmax=vmax, clip=False), bool(np.max(positive) > vmax)


def shared_extent():
    widths = [GAS_X_CM + (NX[s] - 1) * (GAS_X_CM - 3.) for s in range(6)]
    heights = [(MS0_GAS_Y_CM if s == 0 else GAS_Y_CM) +
               (NY[s] - 1) * ((MS0_GAS_Y_CM if s == 0 else GAS_Y_CM) - 3.)
               for s in range(6)]
    return max(widths), max(heights)


def overview(out, arr, exp, source_i, station, norm, cmap):
    src = SOURCES[source_i]
    fig, axes = plt.subplots(1, 3, figsize=(20, 7), layout="constrained")
    width, height = shared_extent()
    for k, ax in enumerate(axes):
        segments, values = [], []
        for ch in GEOM[station]:
            n = n_wires(station)
            pp = arr[source_i, ch.det_id, :n, k]
            # 8 individual anodes per visible band (32 mm nominal spacing).
            for j0 in range(0, n, 8):
                pp_group = pp[j0:j0 + 8]
                if not np.any(np.isfinite(pp_group)):
                    continue
                val = float(np.nanmean(pp_group))
                wire_y = ch.y + float(np.mean(positions(station)[j0:j0+8])) / 10
                segments.append(((ch.x - GAS_X_CM/2 + .8, wire_y),
                                 (ch.x + GAS_X_CM/2 - .8, wire_y)))
                values.append(val)
            gy = MS0_GAS_Y_CM if station == 0 else GAS_Y_CM
            ax.add_patch(Rectangle(
                (ch.x - GAS_X_CM/2, ch.y - gy/2), GAS_X_CM, gy,
                fill=False, edgecolor=".25", linewidth=.6, zorder=4,
            ))
            ax.text(ch.x, ch.y + gy/2 - 2.5, f"{ch.det_id}:{ch.layer}",
                    ha="center", va="top", fontsize=5.2, color=".3", zorder=6)
        lc = LineCollection(segments, cmap=cmap, norm=norm,
                            linewidths=2.8, zorder=3)
        lc.set_array(np.asarray(values))
        ax.add_collection(lc)
        ax.set_xlim(-width/2, width/2)
        ax.set_ylim(-height/2, height/2)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel("Global X [cm], along wires")
        ax.set_ylabel("Global Y [cm], across wires")
        ax.set_title(CAPTIONS[k], fontsize=11)
        ax.grid(False)
        fig.colorbar(lc, ax=ax, shrink=.65, label="Fraction of module crossings")
    number = int(np.sum(exp[source_i, [c.det_id for c in GEOM[station]]]))
    fig.suptitle(
        f"{LABELS[src]}  |  MS{station}  |  {number:,} chamber crossings  |  "
        "8 wires per visible stripe\n"
        "Nominal fixed 4-mm grid; 3 conditions; per-module denominators",
        fontsize=13,
    )
    path = out / "station_overview" / f"wire_overview_{src}_MS{station}.png"
    fig.savefig(path, dpi=135)
    plt.close(fig)


def details(out, arr, exp, source_i, ch):
    src = SOURCES[source_i]
    n = n_wires(ch.station)
    data = arr[source_i, ch.det_id, :n]
    exposure = int(exp[source_i, ch.det_id])
    fig, axes = plt.subplots(3, 1, figsize=(8, 7.0), sharex=True, layout="constrained")
    idx = np.arange(n)
    for k, ax in enumerate(axes):
        values = data[:, k]
        if exposure > 0:
            ax.step(idx, values, where="mid", linewidth=1.1)
            ax.fill_between(idx, values, step="mid", alpha=.2)
            ax.set_ylim(0, max(.01, float(np.nanmax(values))*1.17))
        else:
            ax.text(.5, .5, "No projected gas crossings", transform=ax.transAxes,
                    ha="center", va="center")
        ax.set_ylabel("P per crossing")
        ax.set_title(CAPTIONS[k], fontsize=10)
        ax.grid(alpha=.18)
    axes[-1].set_xlabel("Physical anode index along global Y (0 = bottom wire)")
    fig.suptitle(
        f"{LABELS[src]} | MS{ch.station} | DET{ch.det_id:03d} "
        f"(row={ch.row}, col={ch.col}, {ch.layer}) | {exposure:,} crossings\n"
        f"{n} nominal wires at 4 mm; no avalanche/electronics threshold",
        fontsize=11,
    )
    dest = out / "per_chamber" / src / f"MS{ch.station}"
    dest.mkdir(parents=True, exist_ok=True)
    fig.savefig(dest / f"wire_DET{ch.det_id:03d}.png", dpi=120)
    plt.close(fig)


def contact_sheet(out, arr, exp, source_i, station, norm, cmap):
    src = SOURCES[source_i]
    nx, ny = NX[station], NY[station]
    fig, axes = plt.subplots(ny, nx, figsize=(nx * 2.5, ny * 1.65),
                              squeeze=False, layout="constrained")
    for ch in GEOM[station]:
        ax = axes[ny - 1 - ch.row, ch.col]
        n = n_wires(station)
        vals = arr[source_i, ch.det_id, :n, :]
        ax.imshow(
            np.ma.masked_invalid(vals), origin="lower", interpolation="nearest",
            aspect="auto", norm=norm, cmap=cmap,
        )
        ax.set_title(f"{ch.det_id}:{ch.layer}  N={exp[source_i,ch.det_id]}",
                     fontsize=7, pad=1)
        ax.set_xticks([0, 1, 2])
        ax.set_xticklabels(["1e", "5e", "100ns"], fontsize=6)
        ax.set_yticks([])
    fig.suptitle(
        f"{LABELS[src]} MS{station} | fixed physical wire index bottom to top\n"
        "Three columns per chamber: >=1e, >=5e, <=100 ns",
        fontsize=14,
    )
    sm = ScalarMappable(norm=norm, cmap=cmap)
    fig.colorbar(sm, ax=list(axes.ravel()), shrink=.55,
                 label="P per chamber crossing")
    path = out / "contact_sheets" / f"wire_contacts_{src}_MS{station}.png"
    fig.savefig(path, dpi=120)
    plt.close(fig)


def plot_multiplicity(out, rows, source):
    subset = [r for r in rows if r["source"] == source
              and float(r["phase_mm"]) == 0]
    fig, axes = plt.subplots(2, 3, figsize=(17, 8), layout="constrained")
    labels = ["0", "1", "2", "3", "4", "5+"]
    categories = ("geom", "ge1e", "ge5e", "le100ns")
    for st, ax in enumerate(axes.flat):
        total = [r for r in subset if int(r["station"]) == st]
        for mi, name in enumerate(categories):
            selected = [r for r in total if r["condition"] == name]
            if not selected:
                continue
            row = selected[0]
            n = int(row["crossings"])
            if n == 0:
                continue
            probs = np.asarray([int(row["N"+k]) for k in
                                ("0", "1", "2", "3", "4", "5plus")]) / n
            ax.bar(np.arange(6) + (mi-1.5)*.19, probs, width=.18,
                   label=name)
        ax.set_xticks(np.arange(6), labels)
        ax.set_title(f"MS{st}")
        ax.set_yscale("symlog", linthresh=.001)
        ax.set_ylim(0, 1.1)
        ax.set_ylabel("Fraction of gas crossings")
        ax.set_xlabel("Number of active wires")
        if st == 0: ax.legend(fontsize=8)
    fig.suptitle(
        f"{LABELS[source]} | wire multiplicity: geometric vs Garfield thresholds",
        fontsize=14
    )
    fig.savefig(out / "physics_summary" / f"multiplicity_{source}.png", dpi=140)
    plt.close(fig)


def plot_angles(out, rows, source):
    subset = [r for r in rows if r["source"] == source
              and float(r["phase_mm"]) == 0]
    fig, axes = plt.subplots(2, 3, figsize=(17, 8), layout="constrained")
    for st, ax in enumerate(axes.flat):
        st_rows = [r for r in subset if int(r["station"]) == st]
        for metric in METRICS:
            vals = sorted(
                [r for r in st_rows if r["condition"] == metric
                 and int(r["crossings"]) >= 30],
                key=lambda r: float(r["theta_abs_low_deg"]),
            )
            xx = [(.5*(float(r["theta_abs_low_deg"]) +
                         float(r["theta_abs_high_deg"]))) for r in vals]
            yy = [float(r["P_ge2"]) for r in vals]
            ax.plot(xx, yy, "o-", markersize=3, linewidth=1.2, label=metric)
        ax.set_title(f"MS{st}")
        ax.set_xlim(0, 32)
        ax.set_ylim(0, 1)
        ax.set_xlabel(r"$|\theta_u|$ [deg]")
        ax.set_ylabel("P(>=2 anodes per gas crossing)")
        ax.grid(alpha=.2)
        if st == 0: ax.legend(fontsize=8)
    fig.suptitle(
        f"{LABELS[source]} | geometric vs collected multiplicity vs across-wire angle",
        fontsize=14
    )
    fig.savefig(out / "physics_summary" / f"angle_dependence_{source}.png",
                dpi=140)
    plt.close(fig)


def plot_phase(out, rows, source):
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), layout="constrained")
    for ax, metric in zip(axes, ("ge1e", "ge5e")):
        for phase in PHASES:
            rr = sorted(
                [r for r in rows if r["source"] == source
                 and r["condition"] == metric and float(r["phase_mm"]) == phase],
                key=lambda r: int(r["station"]),
            )
            ax.plot([int(r["station"]) for r in rr],
                    [float(r["P_ge2"]) for r in rr],
                    "o-", label=f"grid phase {phase:+g} mm")
        ax.set_xticks(list(range(6)), [f"MS{k}" for k in range(6)])
        ax.set_ylim(bottom=0)
        ax.set_ylabel("P(>=2 active wires)")
        ax.set_title(metric)
        ax.grid(alpha=.2)
        ax.legend(fontsize=8)
    fig.suptitle(
        f"{LABELS[source]} | physical anode-grid origin sensitivity "
        "(same ionisation and Garfield random quantiles)",
        fontsize=12,
    )
    fig.savefig(out / "physics_summary" / f"phase_sensitivity_{source}.png",
                dpi=140)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--input", type=Path,
        default=Path("test_runs/mwpc_stage_c/pi_fixed_wire_maps"),
    )
    ap.add_argument(
        "--no-details", action="store_true",
        help="Only for quick smoke rendering, no 894 chamber files",
    )
    args = ap.parse_args()
    arr, exp, counts = load_wire_csv(args.input / "pi_fixed_wire_counts.csv")
    mult_rows = read_csv(args.input / "pi_fixed_wire_station_multiplicity.csv")
    angle_rows = read_csv(args.input / "pi_fixed_wire_angle_response.csv")
    for folder in ("station_overview", "contact_sheets", "physics_summary",
                   "per_chamber"):
        (args.input / folder).mkdir(exist_ok=True)
    cmap = plt.get_cmap("magma").copy()
    cmap.set_bad("white")
    for i, src in enumerate(SOURCES):
        norm, clipped = color_norm(arr[i])
        print(
            f"{src}: shared within-channel P color max {norm.vmax:.5f}, "
            f"high-end clipping={clipped}", flush=True,
        )
        for st in range(6):
            overview(args.input, arr, exp, i, st, norm, cmap)
            contact_sheet(args.input, arr, exp, i, st, norm, cmap)
            if not args.no_details:
                for ch in GEOM[st]:
                    details(args.input, arr, exp, i, ch)
            print(f"  Finished MS{st} for {src}", flush=True)
        plot_multiplicity(args.input, mult_rows, src)
        plot_angles(args.input, angle_rows, src)
        plot_phase(args.input, mult_rows, src)
    print(
        f"PNG atlas complete: 18 station overviews, 18 contact sheets, "
        f"{0 if args.no_details else 894} per-chamber detail images, "
        f"9 physics summaries. {args.input}",
        flush=True,
    )


if __name__ == "__main__":
    main()
