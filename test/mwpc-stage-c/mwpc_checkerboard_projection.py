#!/usr/bin/env python3
"""Project historical muon reference states onto the 298-chamber NA60+/DiCE MWPC layout.

Inputs: full bore_muon_crossings.csv files (J/psi, phi, omega), not the balanced
36k response sample. Output is geometric straight-line *projection*, NOT
new Geant4 transport, gas Edep, fired-wire probabilities, or physics rates.
No ROOT/Garfield dependency; matplotlib is optional.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

SOURCES = ("jpsi40-bore", "phi40-bore", "omega40-bore")
NX = (4, 4, 6, 6, 8, 9)
NY = (8, 5, 7, 7, 9, 10)
Z_STATION = (300., 340., None, 590., 810., 850.)
MS2_SCENARIOS = (530., 550.)
GAS_X_CM = 66.56  # installed X: along wires
GAS_Y_CM = 51.20  # installed Y: across wires
MS0_GAS_Y_CM = 30.72
OVERLAP_CM = 3.0
GAS_THICKNESS_CM = 0.6
STAGGER_STEP_CM = 4.0
LAYERS = "ABCD"

REQUIRED = {
    "job", "event", "track", "pdg", "station", "x_cm", "y_cm",
    "z_cm", "px_GeV", "py_GeV", "pz_GeV",
}
CROSSING_FIELDS = (
    "source", "scenario_ms2_z_cm", "job", "event", "track", "pdg",
    "station", "det_id", "row", "col", "layer",
    "chamber_center_x_cm", "chamber_center_y_cm", "gas_center_z_cm",
    "reference_z_cm", "reference_to_gas_dz_cm",
    "x_in_cm", "y_in_cm", "z_in_cm",
    "x_out_cm", "y_out_cm", "z_out_cm",
    "gas_path_mm", "du_mm", "dw_mm", "theta_u_deg", "theta_w_deg", "p_GeV",
)


@dataclass(frozen=True)
class Chamber:
    station: int
    det_id: int
    row: int
    col: int
    layer: str
    x: float
    y: float
    z_offset: float
    half_x: float
    half_y: float


def build_chambers():
    """Replicate NA6PMuonSpec::createGeometry row-major ID and ABCD assignment."""
    groups = []
    det_id = 0
    for station, (nx, ny) in enumerate(zip(NX, NY)):
        gas_y = MS0_GAS_Y_CM if station == 0 else GAS_Y_CM
        pitch_x = GAS_X_CM - OVERLAP_CM
        pitch_y = gas_y - OVERLAP_CM
        group = []
        for row in range(ny):
            for col in range(nx):
                q = 2 * (row & 1) + (col & 1)
                group.append(Chamber(
                    station, det_id, row, col, LAYERS[q],
                    (col - (nx - 1) / 2.) * pitch_x,
                    (row - (ny - 1) / 2.) * pitch_y,
                    (q - 1.5) * STAGGER_STEP_CM,
                    GAS_X_CM / 2., gas_y / 2.,
                ))
                det_id += 1
        groups.append(group)
    assert det_id == 298
    return groups


def intersect_gas(x0, y0, z0, sx, sy, ch, station_z):
    """Return downstream-z interval of ray inside full sensitive gas box."""
    zc = station_z + ch.z_offset
    zmin = zc - GAS_THICKNESS_CM / 2.
    zmax = zc + GAS_THICKNESS_CM / 2.
    for pos, slope, centre, half in (
        (x0, sx, ch.x, ch.half_x),
        (y0, sy, ch.y, ch.half_y),
    ):
        low, high = centre - half, centre + half
        if abs(slope) < 1e-14:
            if pos < low or pos > high:
                return None
            continue
        a = z0 + (low - pos) / slope
        b = z0 + (high - pos) / slope
        zmin = max(zmin, min(a, b))
        zmax = min(zmax, max(a, b))
        if zmax - zmin <= 1.e-10:
            return None
    return (zmin, zmax) if zmax - zmin > 1.e-10 else None


def candidate_chambers(group, station, station_z, x0, y0, z0, sx, sy):
    """Avoid testing all chambers; conservative bounding box for all z layers."""
    low = station_z - 1.5 * STAGGER_STEP_CM - GAS_THICKNESS_CM / 2.
    high = station_z + 1.5 * STAGGER_STEP_CM + GAS_THICKNESS_CM / 2.
    xa, xb = x0 + sx * (low - z0), x0 + sx * (high - z0)
    ya, yb = y0 + sy * (low - z0), y0 + sy * (high - z0)
    gy = MS0_GAS_Y_CM if station == 0 else GAS_Y_CM
    px, py = GAS_X_CM - OVERLAP_CM, gy - OVERLAP_CM
    nx, ny = NX[station], NY[station]
    eps = 1.e-9
    cmin = max(0, math.ceil((min(xa, xb) - GAS_X_CM / 2.) / px + (nx - 1) / 2. - eps))
    cmax = min(nx - 1, math.floor((max(xa, xb) + GAS_X_CM / 2.) / px + (nx - 1) / 2. + eps))
    rmin = max(0, math.ceil((min(ya, yb) - gy / 2.) / py + (ny - 1) / 2. - eps))
    rmax = min(ny - 1, math.floor((max(ya, yb) + gy / 2.) / py + (ny - 1) / 2. + eps))
    for row in range(rmin, rmax + 1):
        for col in range(cmin, cmax + 1):
            yield group[row * nx + col]


def save_table(path, rows):
    if not rows:
        raise ValueError("Cannot save an empty table: " + str(path))
    with path.open("w", newline="", encoding="utf-8") as stream:
        w = csv.DictWriter(stream, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--base", type=Path,
        default=Path.home() /
        "na60-dice-mwpc-simulation/geometry/station-layout-study/full-setup/"
        "mnp33-bore-runs/20260912-202112-0wrxhdg_",
    )
    parser.add_argument(
        "--out", type=Path,
        default=Path("test_runs/mwpc_stage_c/c_checkerboard_projection"),
    )
    parser.add_argument(
        "--limit-per-source", type=int, default=0,
        help="Limit INPUT REFERENCE ROWS per source for smoke test (0 = all).",
    )
    parser.add_argument("--write-crossings", action="store_true",
                        help="Write potentially large per-track/per-module CSVs.")
    parser.add_argument("--no-plots", action="store_true")
    args = parser.parse_args()
    if args.limit_per_source < 0:
        parser.error("--limit-per-source must be non-negative")

    inputs = {s: args.base / "analysis" / s / "bore_muon_crossings.csv"
              for s in SOURCES}
    missing = [str(p) for p in inputs.values() if not p.is_file()]
    if missing:
        parser.error("Missing historical crossing CSVs:\n" + "\n".join(missing))

    args.out.mkdir(parents=True, exist_ok=True)
    groups = build_chambers()
    module_count = Counter()
    module_length = Counter()
    module_tu = Counter()
    module_tw = Counter()
    station_stats = {}
    skipped = Counter()
    opened = {}
    writers = {}

    if args.write_crossings:
        for source in SOURCES:
            for z in MS2_SCENARIOS:
                path = args.out / f"crossings_{source}_ms2_{int(z)}.csv"
                stream = path.open("w", newline="", encoding="utf-8")
                opened[(source, z)] = stream
                writer = csv.DictWriter(stream, fieldnames=CROSSING_FIELDS)
                writer.writeheader()
                writers[(source, z)] = writer

    try:
        for source, path in inputs.items():
            print(f"Reading {source}: {path}", flush=True)
            total = 0
            with path.open(newline="", encoding="utf-8") as stream:
                reader = csv.DictReader(stream)
                absent = REQUIRED - set(reader.fieldnames or [])
                if absent:
                    raise ValueError(f"{path}: missing columns {sorted(absent)}")
                for row in reader:
                    if args.limit_per_source and total >= args.limit_per_source:
                        break
                    total += 1
                    try:
                        st = int(row["station"])
                        pdg = int(row["pdg"])
                        x0, y0, z0 = (float(row[k]) for k in
                                      ("x_cm", "y_cm", "z_cm"))
                        px, py, pz = (float(row[k]) for k in
                                      ("px_GeV", "py_GeV", "pz_GeV"))
                    except (TypeError, ValueError):
                        skipped[(source, "bad_row")] += 1
                        continue
                    if st not in range(6) or abs(pdg) != 13:
                        skipped[(source, "not_forward_muon_station")] += 1
                        continue
                    if not all(map(math.isfinite, (x0, y0, z0, px, py, pz))) or pz <= 0:
                        skipped[(source, "invalid_or_nonforward")] += 1
                        continue

                    sx, sy = px / pz, py / pz
                    tu = math.degrees(math.atan2(py, pz))
                    tw = math.degrees(math.atan2(px, pz))
                    pabs = math.sqrt(px * px + py * py + pz * pz)
                    length_scale = 10. * math.sqrt(1. + sx * sx + sy * sy)

                    for ms2z in MS2_SCENARIOS:
                        zs = ms2z if st == 2 else Z_STATION[st]
                        stats = station_stats.setdefault((source, ms2z, st), Counter())
                        stats["reference_states"] += 1
                        multiplicity = 0
                        for ch in candidate_chambers(
                            groups[st], st, zs, x0, y0, z0, sx, sy,
                        ):
                            slab = intersect_gas(x0, y0, z0, sx, sy, ch, zs)
                            if slab is None:
                                continue
                            zin, zout = slab
                            path_mm = (zout - zin) * length_scale
                            if path_mm <= 1.e-8:
                                continue
                            multiplicity += 1
                            key = (source, ms2z, ch.det_id)
                            module_count[key] += 1
                            module_length[key] += path_mm
                            module_tu[key] += tu
                            module_tw[key] += tw
                            if args.write_crossings:
                                writers[(source, ms2z)].writerow({
                                    "source": source, "scenario_ms2_z_cm": ms2z,
                                    "job": row["job"], "event": row["event"],
                                    "track": row["track"], "pdg": pdg,
                                    "station": st, "det_id": ch.det_id,
                                    "row": ch.row, "col": ch.col, "layer": ch.layer,
                                    "chamber_center_x_cm": ch.x,
                                    "chamber_center_y_cm": ch.y,
                                    "gas_center_z_cm": zs + ch.z_offset,
                                    "reference_z_cm": z0,
                                    "reference_to_gas_dz_cm": zs + ch.z_offset - z0,
                                    "x_in_cm": x0 + sx * (zin - z0),
                                    "y_in_cm": y0 + sy * (zin - z0),
                                    "z_in_cm": zin,
                                    "x_out_cm": x0 + sx * (zout - z0),
                                    "y_out_cm": y0 + sy * (zout - z0),
                                    "z_out_cm": zout,
                                    "gas_path_mm": path_mm,
                                    "du_mm": 10. * sy * (zout - zin),
                                    "dw_mm": 10. * sx * (zout - zin),
                                    "theta_u_deg": tu, "theta_w_deg": tw,
                                    "p_GeV": pabs,
                                })
                        stats["station_with_hit"] += (multiplicity > 0)
                        stats["station_without_hit"] += (multiplicity == 0)
                        stats["chamber_crossings"] += multiplicity
                        stats[f"chamber_multiplicity_{min(multiplicity, 5)}"] += 1
            print(f"  processed reference rows: {total:,}", flush=True)
    finally:
        for stream in opened.values():
            stream.close()

    modules = []
    for source in SOURCES:
        for ms2z in MS2_SCENARIOS:
            for group in groups:
                for ch in group:
                    key = (source, ms2z, ch.det_id)
                    n = module_count[key]
                    ref = station_stats.get(
                        (source, ms2z, ch.station), Counter()
                    )["reference_states"]
                    zs = ms2z if ch.station == 2 else Z_STATION[ch.station]
                    modules.append({
                        "source": source, "scenario_ms2_z_cm": ms2z,
                        "station": ch.station, "det_id": ch.det_id,
                        "row": ch.row, "col": ch.col, "layer": ch.layer,
                        "gas_center_x_cm": ch.x, "gas_center_y_cm": ch.y,
                        "gas_center_z_cm": zs + ch.z_offset,
                        "gas_width_x_cm": GAS_X_CM,
                        "gas_height_y_cm": 2. * ch.half_y,
                        "muon_crossings": n,
                        "fraction_of_reference_states": n / ref if ref else math.nan,
                        "mean_gas_path_mm": module_length[key] / n if n else math.nan,
                        "mean_theta_u_deg": module_tu[key] / n if n else math.nan,
                        "mean_theta_w_deg": module_tw[key] / n if n else math.nan,
                    })
    stations = []
    for source in SOURCES:
        for ms2z in MS2_SCENARIOS:
            for st in range(6):
                cnt = station_stats.get((source, ms2z, st), Counter())
                ref = cnt["reference_states"]
                covered = cnt["station_with_hit"]
                stations.append({
                    "source": source, "scenario_ms2_z_cm": ms2z,
                    "station": st,
                    "station_center_z_cm": ms2z if st == 2 else Z_STATION[st],
                    "reference_states": ref,
                    "projected_muons_with_gas": covered,
                    "projected_muons_without_gas": cnt["station_without_hit"],
                    "muons_fraction_with_gas": covered / ref if ref else math.nan,
                    "chamber_crossings_total": cnt["chamber_crossings"],
                    "extra_chamber_crossings_due_to_overlap":
                        cnt["chamber_crossings"] - covered,
                    **{
                        (f"muons_in_{k}_chambers" if k < 5 else
                         "muons_in_5plus_chambers"):
                        cnt[f"chamber_multiplicity_{k}"]
                        for k in range(6)
                    },
                })
    save_table(args.out / "modules.csv", modules)
    save_table(args.out / "stations.csv", stations)
    lookup = {(r["source"], r["scenario_ms2_z_cm"], r["det_id"]): r
              for r in modules}
    comparison = []
    for source in SOURCES:
        for ch in groups[2]:
            a = lookup[(source, 530., ch.det_id)]
            b = lookup[(source, 550., ch.det_id)]
            comparison.append({
                "source": source, "station": 2, "det_id": ch.det_id,
                "row": ch.row, "col": ch.col,
                "n_ms2_530": a["muon_crossings"],
                "n_ms2_550": b["muon_crossings"],
                "delta_550_minus_530":
                    b["muon_crossings"] - a["muon_crossings"],
            })
    save_table(args.out / "ms2_comparison.csv", comparison)
    metadata = {
        "inputs": {k: str(v) for k, v in inputs.items()},
        "reference_station_z_cm": [300, 340, 550, 590, 810, 850],
        "ms2_scenarios_cm": list(MS2_SCENARIOS),
        "station_NX": list(NX), "station_NY": list(NY),
        "chambers": sum(x * y for x, y in zip(NX, NY)),
        "gas_x_cm": GAS_X_CM, "gas_y_cm": GAS_Y_CM,
        "ms0_gas_y_cm": MS0_GAS_Y_CM,
        "gas_thickness_cm": GAS_THICKNESS_CM, "overlap_cm": OVERLAP_CM,
        "stagger_step_cm": STAGGER_STEP_CM,
        "checkerboard": "q=2*(row%2)+(col%2); z(A,B,C,D)=-6,-2,+2,+6 cm",
        "source_filter": "PDG +-13, forward pz>0, historical reference-plane muons",
        "limit_per_source": args.limit_per_source,
        "skipped": {f"{a}|{b}": n for (a, b), n in skipped.items()},
        "limitations": [
            "Straight-line gas-box projections, not new Geant4 hits.",
            "No energy deposition, thresholds or wire signal model.",
            "Old transported geometries differ from 298-chamber setup.",
            "MS2=530 extrapolates around 20 cm backwards from z=550 reference.",
            "Per-module totals double-count distinct tracks in gas-overlap regions.",
            "Generated dimuon sources are not luminosity-normalized PbPb fluence.",
            "Smoke-test --limit-per-source is a file prefix, not a random sample.",
        ],
    }
    (args.out / "metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
    )

    if not args.no_plots:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            for source in SOURCES:
                for ms2z in MS2_SCENARIOS:
                    fig, axes = plt.subplots(
                        2, 3, figsize=(14, 8), constrained_layout=True
                    )
                    for st, ax in enumerate(axes.flat):
                        arr = [
                            [
                                module_count[(source, ms2z, ch.det_id)]
                                for ch in groups[st][r * NX[st]:(r + 1) * NX[st]]
                            ]
                            for r in range(NY[st])
                        ]
                        im = ax.imshow(arr, origin="lower", aspect="auto")
                        ax.set_title(f"MS{st}, chamber crossings={sum(map(sum, arr)):,}")
                        ax.set_xlabel("column: X (along wires)")
                        ax.set_ylabel("row: Y (across wires)")
                        fig.colorbar(im, ax=ax, shrink=0.7)
                    fig.suptitle(f"{source} | MS2 z={ms2z:.0f} cm")
                    fig.savefig(
                        args.out / f"map_{source}_ms2_{int(ms2z)}.png", dpi=150
                    )
                    plt.close(fig)
        except ImportError:
            print("matplotlib unavailable: CSVs saved; no maps")

    print(f"Saved to: {args.out}")
    print("MS2: source, z, reference states, unique muons in gas, gas crossings")
    for row in stations:
        if row["station"] == 2:
            print(
                f'{row["source"]:13s} {row["scenario_ms2_z_cm"]:.0f}'
                f' {row["reference_states"]:8,d}'
                f' {row["projected_muons_with_gas"]:8,d}'
                f' {row["chamber_crossings_total"]:8,d}'
            )
    print("CAUTION: projected geometry, not Geant4 Edep, hit efficiency or fired wires.")


if __name__ == "__main__":
    main()
