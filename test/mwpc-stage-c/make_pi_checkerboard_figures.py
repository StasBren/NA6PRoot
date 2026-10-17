#!/usr/bin/env python3
"""Make PI-ready NA60+/DiCE checkerboard maps and accepted-muon angle spectra.

Inputs: Stage-C reconstructed gas crossings, one row per MODULE crossing.
Source production was historical GEANT4 reference-plane trajectories; results
are straight-line geometrical projections, not direct GEANT4 chamber hits.
Requires numpy and matplotlib, but not ROOT, Garfield or pandas.
"""
from __future__ import annotations

import argparse
import csv
from collections import Counter
from itertools import groupby
from math import ceil, isfinite
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.colors import LogNorm, Normalize, PowerNorm

from mwpc_checkerboard_projection import (
    SOURCES, NX, NY, GAS_X_CM, GAS_Y_CM, MS0_GAS_Y_CM,
    OVERLAP_CM, build_chambers,
)

CHANNEL_LABELS = {
    "jpsi40-bore": "J/psi", "phi40-bore": "phi", "omega40-bore": "omega",
}
ANGLE_VARIABLES = ("theta_u_deg", "theta_w_deg")
ANGLE_LABELS = {
    "theta_u_deg": r"$\theta_u$ (across wires)",
    "theta_w_deg": r"$\theta_w$ (along wires)",
}
WIRE_PITCH_CM = 0.4
GEOMETRY = build_chambers()


def footprint(station):
    gas_y = MS0_GAS_Y_CM if station == 0 else GAS_Y_CM
    width = GAS_X_CM + (NX[station] - 1) * (GAS_X_CM - OVERLAP_CM)
    height = gas_y + (NY[station] - 1) * (gas_y - OVERLAP_CM)
    return width, height


def allocate(bin_cm, xy_range="shared"):
    """Common centred physical XY frame for all stations, unless overridden."""
    if xy_range not in ("shared", "station"):
        raise ValueError(f"Unknown xy_range: {xy_range}")
    shared_width = max(footprint(st)[0] for st in range(6))
    shared_height = max(footprint(st)[1] for st in range(6))
    data = {}
    for station in range(6):
        width, height = (
            (shared_width, shared_height) if xy_range == "shared"
            else footprint(station)
        )
        nxb, nyb = ceil(width / bin_cm), ceil(height / bin_cm)
        data[station] = {
            "chamber_counts": np.zeros((NY[station], NX[station]), dtype=np.int64),
            "spatial": np.zeros((nyb, nxb), dtype=np.int64),
            "xedges": np.linspace(-width / 2, width / 2, nxb + 1),
            "yedges": np.linspace(-height / 2, height / 2, nyb + 1),
            "angles": {field: [] for field in ANGLE_VARIABLES},
            "n_muons": 0, "n_crossings": 0, "n_multi": 0,
            "multiplicity": Counter(), "det_count": Counter(),
        }
    return data


def track_key(row):
    return (row["job"], row["event"], row["track"], row["station"])


def consume_crossings(path, entries, limit_groups=0):
    """A source CSV is ordered by reference track. Module hits for a muon
    and station are contiguous, as emitted by mwpc_checkerboard_projection.py."""
    required = {
        "job", "event", "track", "station", "det_id", "row", "col", "layer",
        "x_in_cm", "y_in_cm", "x_out_cm", "y_out_cm",
        "theta_u_deg", "theta_w_deg",
    }
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path}: missing columns {sorted(missing)}")
        processed = 0
        for key, records in groupby(reader, track_key):
            if limit_groups and processed >= limit_groups:
                break
            processed += 1
            station = int(key[3])
            if station not in entries:
                raise ValueError(f"Invalid station {station}")
            data = entries[station]
            det_ids = set()
            first = True
            for row in records:
                if first:
                    for field in ANGLE_VARIABLES:
                        theta = float(row[field])
                        if not isfinite(theta):
                            raise ValueError(f"Invalid angle {field}: {theta}")
                        data["angles"][field].append(theta)
                    data["n_muons"] += 1
                    first = False

                det_id = int(row["det_id"])
                if det_id in det_ids:
                    raise ValueError(
                        f"Repeated det_id={det_id} for reference state {key}"
                    )
                det_ids.add(det_id)
                rr, cc = int(row["row"]), int(row["col"])
                if not (0 <= rr < NY[station] and 0 <= cc < NX[station]):
                    raise ValueError(f"Invalid module row/col: {row}")
                chamber = GEOMETRY[station][rr * NX[station] + cc]
                if chamber.det_id != det_id or chamber.layer != row["layer"]:
                    raise ValueError(f"Geometry mismatch for MS{station} det_id={det_id}")

                data["chamber_counts"][rr, cc] += 1
                data["det_count"][det_id] += 1
                data["n_crossings"] += 1

                # IMPORTANT: one spatial entry for EACH gas-volume crossing.
                # A muon crossing two chambers contributes two spatial entries.
                x = .5 * (float(row["x_in_cm"]) + float(row["x_out_cm"]))
                y = .5 * (float(row["y_in_cm"]) + float(row["y_out_cm"]))
                xe, ye = data["xedges"], data["yedges"]
                if not (xe[0] - 1e-6 <= x <= xe[-1] + 1e-6 and
                        ye[0] - 1e-6 <= y <= ye[-1] + 1e-6):
                    raise ValueError(
                        f"Midpoint outside MS{station} footprint: ({x}, {y})"
                    )
                ix = min(len(xe) - 2, max(
                    0, int(np.searchsorted(xe, x, side="right") - 1)))
                iy = min(len(ye) - 2, max(
                    0, int(np.searchsorted(ye, y, side="right") - 1)))
                data["spatial"][iy, ix] += 1

            mult = len(det_ids)
            if not mult:
                raise AssertionError(f"Empty muon/station group: {key}")
            data["multiplicity"][mult] += 1
            data["n_multi"] += int(mult >= 2)
    return processed


def validate(source, entries, input_dir, scenario, complete):
    """Internal sum rules and independent check against Stage-C projector."""
    for station, data in entries.items():
        hits = data["n_crossings"]
        assert int(data["chamber_counts"].sum()) == hits
        assert int(data["spatial"].sum()) == hits
        assert sum(k * v for k, v in data["multiplicity"].items()) == hits
        assert sum(data["multiplicity"].values()) == data["n_muons"]
        assert all(len(data["angles"][a]) == data["n_muons"]
                   for a in ANGLE_VARIABLES)

    if not complete:
        return

    station_file = input_dir / "stations.csv"
    if not station_file.is_file():
        raise FileNotFoundError(
            f"Missing full-run cross-check table: {station_file}"
        )
    with station_file.open(newline="", encoding="utf-8") as stream:
        rows = [r for r in csv.DictReader(stream)
                if r["source"] == source
                and float(r["scenario_ms2_z_cm"]) == scenario]
    if len(rows) != 6:
        raise ValueError(f"Expected six station rows for {source}")
    for row in rows:
        station = int(row["station"])
        data = entries[station]
        expected_hits = int(row["chamber_crossings_total"])
        expected_muons = int(row["projected_muons_with_gas"])
        if (data["n_crossings"] != expected_hits
                or data["n_muons"] != expected_muons):
            raise ValueError(
                f"MS{station} {source}: stations.csv says "
                f"{expected_hits} hits/{expected_muons} muons, "
                f"but reconstruction has "
                f"{data['n_crossings']} hits/{data['n_muons']} muons"
            )

    module_file = input_dir / "modules.csv"
    if not module_file.is_file():
        raise FileNotFoundError(
            f"Missing full-run cross-check table: {module_file}"
        )
    with module_file.open(newline="", encoding="utf-8") as stream:
        rows = [r for r in csv.DictReader(stream)
                if r["source"] == source
                and float(r["scenario_ms2_z_cm"]) == scenario]
    if len(rows) != 298:
        raise ValueError(f"Expected 298 module rows for {source}")
    for row in rows:
        st, det = int(row["station"]), int(row["det_id"])
        if int(row["muon_crossings"]) != entries[st]["det_count"][det]:
            raise ValueError(
                f"MS{st}, module {det}: crossings differ from modules.csv"
            )


def build_color_norm(values, mode="sqrt", upper_quantile=0.995):
    """Normalize raw occupancy without altering its counts.

    The percentile is calculated over strictly positive bins, so empty area
    does not dilute the high-end cutoff. Saturated bins are explicitly
    signalled by the colorbar's 'extend=max' indicator.
    """
    arr = np.asarray(values, dtype=float)
    positive = arr[np.isfinite(arr) & (arr > 0)]
    if positive.size == 0:
        return Normalize(vmin=0.0, vmax=1.0), False

    maximum = float(positive.max())
    upper = (float(np.quantile(positive, upper_quantile))
             if upper_quantile < 1.0 else maximum)
    upper = max(upper, float(positive.min()))
    clipped = maximum > upper + max(1.e-12, abs(maximum) * 1.e-12)

    if mode == "linear":
        return Normalize(vmin=0.0, vmax=upper), clipped
    if mode == "sqrt":
        return PowerNorm(gamma=0.5, vmin=0.0, vmax=upper), clipped
    if mode == "log":
        lower = float(positive.min())
        if upper <= lower:
            upper = lower * (1.0 + 1.e-6)
        return LogNorm(vmin=lower, vmax=upper), clipped
    raise ValueError(f"Unknown color normalization mode: {mode}")


def occupancy_cmap():
    cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad("white")  # genuinely empty bins have no hits
    return cmap


def plot_spatial(ax, data, station, mode="sqrt", upper_quantile=0.995):
    xe, ye = data["xedges"], data["yedges"]
    pixel_area = (xe[1] - xe[0]) * (ye[1] - ye[0])
    values = data["spatial"] / pixel_area
    norm, saturated = build_color_norm(values, mode, upper_quantile)
    image = ax.imshow(
        np.ma.masked_less_equal(values, 0),
        origin="lower", interpolation="nearest",
        extent=[xe[0], xe[-1], ye[0], ye[-1]], aspect="equal",
        norm=norm, cmap=occupancy_cmap(),
    )
    for ch in GEOMETRY[station]:
        half_y = (MS0_GAS_Y_CM if station == 0 else GAS_Y_CM) / 2
        ax.add_patch(Rectangle(
            (ch.x - GAS_X_CM / 2, ch.y - half_y), GAS_X_CM, 2 * half_y,
            fill=False, linewidth=.45, edgecolor="white", alpha=.28,
        ))
    ax.set_xlabel("Global X (cm), along wires")
    ax.set_ylabel("Global Y (cm), across wires")
    ax.set_title(
        f"Spatial hit density (per chamber crossing; {mode} scale)",
        fontsize=11,
    )
    return image, saturated


def plot_modules(ax, data, station, mode="sqrt", upper_quantile=0.995):
    """Per-chamber crossing counts positioned at their exact physical XY centres.

    Each colored tile represents one module centre (Voronoi-style), not the
    full gas-volume polygon: neighboring gas volumes overlap by 3 cm and
    may independently register the same muon.
    """
    counts = data["chamber_counts"]
    norm, saturated = build_color_norm(counts, mode, upper_quantile)
    pitch_x = GAS_X_CM - OVERLAP_CM
    pitch_y = (MS0_GAS_Y_CM if station == 0 else GAS_Y_CM) - OVERLAP_CM
    x_boundaries = (np.arange(NX[station] + 1) - NX[station] / 2) * pitch_x
    y_boundaries = (np.arange(NY[station] + 1) - NY[station] / 2) * pitch_y
    image = ax.pcolormesh(
        x_boundaries, y_boundaries, np.ma.masked_less_equal(counts, 0),
        norm=norm, cmap=occupancy_cmap(), shading="flat",
    )
    width, height = footprint(station)
    ax.add_patch(Rectangle(
        (-width / 2, -height / 2), width, height,
        fill=False, linewidth=0.7, edgecolor="0.5",
    ))
    for ch in GEOMETRY[station]:
        count = counts[ch.row, ch.col]
        intensity = float(norm(float(count))) if count else 0.0
        label_color = (
            "0.35" if count == 0
            else "white" if intensity < .63 else "black"
        )
        ax.text(
            ch.x, ch.y, ch.layer, ha="center", va="center",
            fontsize=6.2 if station >= 4 else 8,
            color=label_color, alpha=.85,
        )
    ax.set_xlim(float(data["xedges"][0]), float(data["xedges"][-1]))
    ax.set_ylim(float(data["yedges"][0]), float(data["yedges"][-1]))
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("Global X (cm), along wires")
    ax.set_ylabel("Global Y (cm), across wires")
    ax.set_title(
        f"Crossings / chamber; A/B/C/D = stagger layer ({mode} scale)",
        fontsize=11,
    )
    return image, saturated


def make_occupancy_figures(
    output_dir, source, entries, scenario, mode="sqrt", upper_quantile=0.995,
):
    folder = output_dir / "occupancy"
    folder.mkdir(parents=True, exist_ok=True)
    for station, data in entries.items():
        fig, axs = plt.subplots(
            1, 2, figsize=(14.4, 6.5), layout="constrained"
        )
        image, clipped = plot_spatial(
            axs[0], data, station, mode, upper_quantile
        )
        fig.colorbar(
            image, ax=axs[0], shrink=.8, label="Crossings / cm$^2$",
            extend="max" if clipped else "neither",
        )
        image, clipped_modules = plot_modules(
            axs[1], data, station, mode, upper_quantile
        )
        fig.colorbar(
            image, ax=axs[1], shrink=.8, label="Crossings / module",
            extend="max" if clipped_modules else "neither",
        )
        n_muon = data["n_muons"]
        n_hits = data["n_crossings"]
        p_multi = data["n_multi"] / n_muon if n_muon else 0.
        limit_label = (
            f"positive-bin P{100 * upper_quantile:g}" if upper_quantile < 1.0
            else "full range"
        )
        fig.suptitle(
            f"{CHANNEL_LABELS[source]} | MS{station} | "
            f"298-module checkerboard | MS2 at {scenario:.0f} cm\n"
            f"{n_hits:,} chamber crossings from {n_muon:,} distinct muons | "
            f"multi-module muons: {p_multi:.1%} | "
            f"{mode} colors, {limit_label}",
            fontsize=12,
        )
        stem = f"hitmap_{source}_MS{station}_ms2_{int(scenario)}"
        fig.savefig(folder / (stem + ".png"), dpi=175)
        plt.close(fig)


def make_angle_figures(output_dir, source, entries, scenario):
    """Six station panels; one value per unique accepted muon/station."""
    folder = output_dir / "angles"
    folder.mkdir(parents=True, exist_ok=True)
    rows = []
    for field in ANGLE_VARIABLES:
        xlim = 40 if field == "theta_u_deg" else 60
        edges = np.linspace(
            -xlim, xlim, 161 if xlim == 40 else 241
        )
        heights, means, stds, q95s, n_out = [], [], [], [], []
        for station in range(6):
            angles = np.asarray(entries[station]["angles"][field])
            hist, _ = np.histogram(angles, bins=edges)
            density = (
                hist / (len(angles) * (edges[1] - edges[0]))
                if len(angles) else hist
            )
            heights.append(density)
            means.append(float(angles.mean()) if len(angles) else np.nan)
            stds.append(float(angles.std()) if len(angles) else np.nan)
            q95s.append(float(np.quantile(np.abs(angles), .95))
                        if len(angles) else np.nan)
            n_out.append(int(np.count_nonzero(
                (angles < edges[0]) | (angles > edges[-1])
            )))

        ymax = max((float(v.max()) for v in heights), default=1.) * 1.15
        if ymax <= 0:
            ymax = 1.
        fig, axs = plt.subplots(
            2, 3, figsize=(14, 8.7),
            sharex=True, sharey=True, layout="constrained",
        )
        for station, ax in enumerate(axs.flat):
            ax.stairs(heights[station], edges, fill=True, alpha=.65)
            ax.axvline(0, color="0.35", lw=.9, ls="--")
            ax.set_ylim(0, ymax)
            ax.set_title(
                f"MS{station} — N={entries[station]['n_muons']:,}"
            )
            ax.text(
                .96, .94,
                f"mean={means[station]:+.2f}°\n"
                f"σ={stds[station]:.2f}°\n"
                f"P95 |θ|={q95s[station]:.2f}°",
                transform=ax.transAxes, ha="right", va="top",
                fontsize=9,
                bbox={"facecolor": "white", "alpha": .85,
                      "edgecolor": "none"},
            )
            if station >= 3:
                ax.set_xlabel(ANGLE_LABELS[field] + " [deg]")
            if station % 3 == 0:
                ax.set_ylabel("Density (1/degree)")
            rows.append({
                "source": source, "ms2_z_cm": scenario,
                "station": station, "axis": field,
                "unique_muons": entries[station]["n_muons"],
                "mean_deg": means[station], "std_deg": stds[station],
                "p95_abs_deg": q95s[station],
                "out_of_display_range": n_out[station],
            })
        fig.suptitle(
            f"{CHANNEL_LABELS[source]} | accepted muon incidence angles | "
            f"{ANGLE_LABELS[field]}\n"
            "ONE entry per muon/station "
            "(multi-chamber tracks counted once)",
            fontsize=14,
        )
        stem = (
            f"angles_{source}_{field.replace('_deg', '')}"
            f"_ms2_{int(scenario)}"
        )
        fig.savefig(folder / (stem + ".png"), dpi=175)
        plt.close(fig)
    return rows


def write_csv(path, rows):
    if not rows:
        raise ValueError(f"Cannot write empty summary {path}")
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)


def run(args):
    args.out.mkdir(parents=True, exist_ok=True)
    # Remove only PDFs produced by older versions of this script.
    # Keeping them would misleadingly suggest current figures have PDF peers.
    for directory, pattern in (
        (args.out / "occupancy", "hitmap_*_MS*_ms2_*.pdf"),
        (args.out / "angles", "angles_*_ms2_*.pdf"),
    ):
        if directory.is_dir():
            for old_pdf in directory.glob(pattern):
                old_pdf.unlink()
    station_rows, angle_rows = [], []
    full_run = args.limit_per_source == 0
    for source in SOURCES:
        filename = (
            args.input / f"crossings_{source}_ms2_{int(args.ms2_z)}.csv"
        )
        if not filename.is_file():
            raise FileNotFoundError(
                f"Missing {filename}; run "
                "run_c_checkerboard_projection.sh full first"
            )
        entries = allocate(args.bin_cm, args.xy_range)
        print(f"Reading {filename} ...", flush=True)
        n_groups = consume_crossings(
            filename, entries, args.limit_per_source
        )
        validate(
            source, entries, args.input, args.ms2_z, full_run
        )
        for st, data in entries.items():
            n_muon, n_hits = data["n_muons"], data["n_crossings"]
            row = {
                "source": source, "ms2_z_cm": args.ms2_z,
                "station": st,
                "unique_muon_station_crossings": n_muon,
                "chamber_crossings": n_hits,
                "extra_chamber_crossings": n_hits - n_muon,
                "muons_with_2plus_chambers": data["n_multi"],
                "multi_chamber_muon_fraction":
                    data["n_multi"] / n_muon if n_muon else np.nan,
            }
            for k in (1, 2, 3, 4):
                row[f"muons_with_{k}_chambers"] = data["multiplicity"][k]
            row["muons_with_5plus_chambers"] = sum(
                v for k, v in data["multiplicity"].items() if k >= 5
            )
            station_rows.append(row)
            print(
                f"  MS{st}: {n_muon:,} unique muons, "
                f"{n_hits:,} chamber crossings, "
                f"{data['n_multi']/n_muon:.2%} multi-chamber",
                flush=True,
            )
        make_occupancy_figures(
            args.out, source, entries, args.ms2_z,
            args.occupancy_scale, args.occupancy_clip_quantile,
        )
        angle_rows.extend(
            make_angle_figures(args.out, source, entries, args.ms2_z)
        )
        print(
            f"  Processed {n_groups:,} contributing "
            "muon/station groups",
            flush=True,
        )

    write_csv(
        args.out / "pi_station_crossing_summary.csv", station_rows
    )
    write_csv(
        args.out / "pi_muon_angles_summary.csv", angle_rows
    )
    (args.out / "READ_ME_FIRST.txt").write_text(
        "NA60+/DiCE MWPC checkerboard: PI geometry/angle figures\n"
        "===================================================\n\n"
        "18 occupancy figures = 3 dimuon channels x 6 stations. PNG only.\n"
        f"XY axis framing: {args.xy_range}; same physical XY limits across "
        "all six stations when set to shared.\n"
        "Chamber counts use tiles centered on real XY positions; "
        "the physical sensitive gas rectangles can overlap.\n"
        f"Occupancy color scale: {args.occupancy_scale}, "
        f"upper quantile of positive bins: {args.occupancy_clip_quantile}.\n"
        "Bins above the displayed color maximum are saturated and marked "
        "by a colorbar extension. The counts themselves are unchanged.\n"
        "Empty bins are white. Colorbar labels retain original physical "
        "counts, not transformed intensities.\n"
        "  LEFT: spatial density of gas-crossing midpoints. "
        "ONE entry PER CHAMBER CROSSING.\n"
        "  RIGHT: counts per chamber; A/B/C/D are checkerboard Z levels.\n"
        "6 angle figures in PNG = 3 channels x theta_u/theta_w, "
        "each with MS0..MS5.\n"
        "  Angles use ONE entry per accepted muon/station, "
        "even for multi-chamber crossings.\n\n"
        f"MS2 baseline z = {args.ms2_z:.0f} cm.\n"
        "Standard gas 66.56x51.20 cm, narrow MS0 66.56x30.72 cm.\n"
        "Module grid: 4x8, 4x5, 6x7, 6x7, 8x9, 9x10 = 298 modules.\n"
        "Active overlap X/Y=3 cm; stagger A/B/C/D = -6,-2,+2,+6 cm.\n"
        "Wires along X; 4 mm wire pitch in reference Garfield model.\n"
        "Approx. wires: standard 128, MS0 77 (rounded).\n"
        "Approx. total: 36,512 wires, pending edge/wire-origin design.\n\n"
        "LIMITATIONS: Historical GEANT4 states projected by straight lines "
        "onto the new gas boxes.\n"
        "NOT direct gas hits from the updated geometry; "
        "not real Edep or readout efficiency.\n"
        "Each source simulated separately; no physical luminosity "
        "normalization or unbiased efficiency.\n"
        "Normal full run cross-checks stations.csv and modules.csv.\n"
        "Smoke-prefix mode skips checks against full production tables.\n",
        encoding="utf-8",
    )
    print(
        f"Saved 18 PNG dual-panel occupancy figures, "
        f"6 angle figures, 2 CSV summaries: {args.out}"
    )


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--input", type=Path,
        default=Path("test_runs/mwpc_stage_c/checkerboard_full"),
    )
    ap.add_argument(
        "--out", type=Path,
        default=Path("test_runs/mwpc_stage_c/pi_checkerboard_figures"),
    )
    ap.add_argument(
        "--ms2-z", type=float, choices=[530., 550.], default=530.,
    )
    ap.add_argument(
        "--bin-cm", type=float, default=2.,
        help="Spatial 2D hitmap bin width [cm]",
    )
    ap.add_argument(
        "--xy-range", choices=("shared", "station"), default="shared",
        help="Global physical X/Y axes: shared across stations (default) "
             "or fitted individually to each station",
    )
    ap.add_argument(
        "--limit-per-source", type=int, default=0,
        help="Muon/station groups per channel: smoke only",
    )
    ap.add_argument(
        "--occupancy-scale", choices=("sqrt", "linear", "log"),
        default="sqrt",
        help="Color normalization for occupancy plots (default: sqrt)",
    )
    ap.add_argument(
        "--occupancy-clip-quantile", type=float, default=0.995,
        help="Upper color-limit quantile of positive bins (default: .995); "
             "1.0 means no clipping",
    )
    args = ap.parse_args()
    if args.bin_cm <= 0 or args.limit_per_source < 0:
        ap.error(
            "--bin-cm must be positive and "
            "--limit-per-source must be non-negative"
        )
    if not 0 < args.occupancy_clip_quantile <= 1:
        ap.error("--occupancy-clip-quantile must be in (0, 1]")
    run(args)


if __name__ == "__main__":
    main()
