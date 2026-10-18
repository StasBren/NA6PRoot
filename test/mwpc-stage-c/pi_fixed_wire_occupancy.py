#!/usr/bin/env python3
"""Per-anode occupancy on fixed nominal wire grids for all 298 MWPC chambers.

Inputs: full checkerboard projected crossings, empirical Heed library, saved
Garfield single-electron transport-box responses. No new Geant4 or Garfield run.

Each chamber has a fixed 4-mm pitch grid centred inside its active gas:
77 MS0 wires, 128 other wires. For a given gas crossing a Poisson Heed cluster
population is sampled along its ACTUAL 3-D gas segment, with empirical Heed
cluster electron multiplicities. All primary electrons sample Garfield transport
responses by their initial (u,v) box. A collected endpoint falling outside the
finite nominal wire grid is counted as an edge loss, not put onto a fake wire.

The simulation includes three coherent mechanical grid phases (nominal, -1,
+1 mm), sharing underlying ionisation and random transport quantiles to make
phase comparisons less noisy. Only nominal-grid per-wire counts are exported.

Interpretation: "fired" means >=1 or >=5 *collected primary electrons*,
or >=1 collected within 100 ns, NOT a measured electronics threshold or
simulated amplified signal. The 6-mm projected gas is mapped to Garfield's
[-2,+4] mm gap with the existing +1-mm coordinate translation. This is an
explicit approximation, not an electrostatic equivalence.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import time
from pathlib import Path

import numpy as np

from mwpc_checkerboard_projection import (
    SOURCES, NX, NY, GAS_Y_CM, MS0_GAS_Y_CM, build_chambers,
)
from run_pi_per_chamber_response import read_garfield, read_heed

PITCH = 4.0
PHASES = (-1.0, 0.0, 1.0)
CONDITIONS = ("ge1e", "ge5e", "le100ns")
METRICS = ("geom", "ge1e", "ge5e", "le100ns")
ANGLE_EDGES = (0., 2., 4., 6., 8., 10., 12., 15., 20., 30., 90.)
CHAMBERS = build_chambers()
LOOKUP = {c.det_id: c for group in CHAMBERS for c in group}
assert len(LOOKUP) == 298


def n_wires(station):
    return 77 if station == 0 else 128


def positions(station, phase=0.):
    n = n_wires(station)
    return (np.arange(n, dtype=float) - (n - 1.) / 2.) * PITCH + phase


def nominal_wire_index(u_mm, station, phase):
    """Nearest fixed wire on infinite periodic grid; may fall outside gas."""
    return np.floor(
        (np.asarray(u_mm) - phase) / PITCH + (n_wires(station) - 1) / 2. + .5
    ).astype(np.int64)


def cluster_segment(row, rng, lam, marks):
    """One ionisation realization reused at all tested physical wire phases."""
    path = float(row["gas_path_mm"])
    if not (math.isfinite(path) and 0 < path < 1000):
        raise ValueError(f"Bad gas path: {path}")
    ncl = int(rng.poisson(lam * path))
    if ncl == 0:
        return np.empty(0), np.empty(0), 0, path
    s = rng.random(ncl)
    mult = marks[rng.integers(len(marks), size=ncl)]
    frac = np.repeat(s, mult)
    yin, yout = float(row["y_in_cm"]), float(row["y_out_cm"])
    zin, zout = float(row["z_in_cm"]), float(row["z_out_cm"])
    zc = float(row["gas_center_z_cm"])
    cy = float(row["chamber_center_y_cm"])
    u_mm = 10. * (yin + frac * (yout - yin) - cy)
    v_mm = 10. * (zin + frac * (zout - zin) - zc) + 1.
    if v_mm.min() < -2.001 or v_mm.max() > 4.001:
        raise ValueError(f"Gas/transport gap mismatch: {v_mm.min()}..{v_mm.max()}")
    return u_mm, v_mm, ncl, path


def electronic_response(u_mm, v_mm, station, phase, lookup, random01):
    """Count electrons deposited by each finite anode (no gain/electronics).

    Return wire-index arrays, per-wire electron counts, <=100ns counts and
    separately the number of collected electrons landing beyond finite wires.
    """
    offsets, sizes, response = lookup
    n = n_wires(station)
    if not u_mm.size:
        return np.zeros(n, dtype=np.int32), np.zeros(n, dtype=np.int32), 0

    closest = nominal_wire_index(u_mm, station, phase)
    local = u_mm - (closest - (n - 1.) / 2.) * PITCH - phase
    iu = np.clip(np.floor((local + 2.) / .2).astype(np.int64), 0, 19)
    iv = np.clip(np.floor((v_mm + 2.) / .2).astype(np.int64), 0, 29)
    ibox = 30 * iu + iv
    chosen = offsets[ibox] + np.minimum(
        (random01 * sizes[ibox]).astype(np.int64), sizes[ibox] - 1
    )
    trials = response[chosen]
    collected = trials[:, 0] > .5
    dest = closest + trials[:, 1].astype(np.int64)
    inside = collected & (dest >= 0) & (dest < n)
    edge_loss = int(np.count_nonzero(collected & ~inside))
    counts = np.bincount(dest[inside], minlength=n).astype(np.int32)
    early = inside & np.isfinite(trials[:, 2]) & (trials[:, 2] >= 0.) & (trials[:, 2] <= 100.)
    early_counts = np.bincount(dest[early], minlength=n).astype(np.int32)
    return counts, early_counts, edge_loss


def geometric_cells(row, station, phase):
    """Number of finite wire cells intersected by straight gas segment."""
    offset = float(row["chamber_center_y_cm"])
    a = nominal_wire_index(10. * (float(row["y_in_cm"]) - offset), station, phase)
    b = nominal_wire_index(10. * (float(row["y_out_cm"]) - offset), station, phase)
    start, stop = sorted((int(a), int(b)))
    return max(0, min(n_wires(station) - 1, stop) - max(0, start) + 1)


def load_expected(path):
    expected = {(src, d): 0 for src in SOURCES for d in LOOKUP}
    with path.open(newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            if float(r["scenario_ms2_z_cm"]) != 530.:
                continue
            key = (r["source"], int(r["det_id"]))
            if key not in expected:
                raise ValueError(f"Unexpected chamber in modules.csv: {key}")
            if int(r["station"]) != LOOKUP[key[1]].station:
                raise ValueError(f"Station/det_id mismatch in modules.csv: {key}")
            expected[key] = int(r["muon_crossings"])
    if len(expected) != 3 * 298:
        raise ValueError("Incomplete expected-count table")
    return expected


def calc(args):
    out = args.out
    out.mkdir(parents=True, exist_ok=True)
    modules = args.crossings / "modules.csv"
    for input_path in (args.heed_tracks, args.heed_clusters, args.garfield, modules):
        if not input_path.is_file():
            raise FileNotFoundError(f"Required input missing: {input_path}")
    expected = load_expected(modules)
    lam, marks = read_heed(args.heed_tracks, args.heed_clusters)
    response = read_garfield(args.garfield)
    print(f"Heed density {lam:.6f} clusters/mm, marks {len(marks):,}", flush=True)

    # Arrays for every source/module/wire; zero-exposure module kept as NaN.
    wire_counts = np.zeros((3, 298, 128, len(CONDITIONS)), dtype=np.int64)
    hits = np.zeros((3, 298), dtype=np.int64)
    edge_electrons = np.zeros((3, 298), dtype=np.int64)
    # [source, station, phase, metric, wire multiplicity category (0..4,5+)]
    mult = np.zeros((3, 6, 3, len(METRICS), 6), dtype=np.int64)
    # [source, exact physical detector, grid phase, criterion, 0..4/5+]
    # This directly implements PI's P(0/1/2/3/4/5+ anodes | ONE gas hit)
    # for EVERY real det_id, unlike the older station-only averages.
    module_mult = np.zeros(
        (3, 298, len(PHASES), len(METRICS), 6), dtype=np.int64
    )
    # per-theta bin [source, station, theta bin, phase, metric] -> (#>=2)
    angle_num = np.zeros((3, 6, len(ANGLE_EDGES) - 1, 3, len(METRICS)), dtype=np.int64)
    angle_den = np.zeros((3, 6, len(ANGLE_EDGES) - 1), dtype=np.int64)
    t0 = time.monotonic()
    rng = np.random.default_rng(args.seed)

    for source_i, source in enumerate(SOURCES):
        input_path = args.crossings / f"crossings_{source}_ms2_530.csv"
        if not input_path.is_file():
            raise FileNotFoundError(f"Missing full projected crossing CSV: {input_path}")
        print(f"Processing {source}: {input_path}", flush=True)
        with input_path.open(newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            needed = {
                "source", "station", "det_id", "y_in_cm", "y_out_cm",
                "z_in_cm", "z_out_cm", "gas_center_z_cm",
                "chamber_center_y_cm", "gas_path_mm", "theta_u_deg",
            }
            missing = needed - set(reader.fieldnames or [])
            if missing:
                raise ValueError(f"Missing CSV columns: {sorted(missing)}")
            for index, row in enumerate(reader, 1):
                st, det = int(row["station"]), int(row["det_id"])
                if row["source"] != source or LOOKUP[det].station != st:
                    raise ValueError(f"Wrong source/station/det_id: {row}")
                ch = LOOKUP[det]
                if abs(float(row["chamber_center_y_cm"]) - ch.y) > .0001:
                    raise ValueError(f"Incorrect chamber centre: {det}")
                u, v, ncl, path = cluster_segment(row, rng, lam, marks)
                q = rng.random(u.size)
                theta = abs(float(row["theta_u_deg"]))
                itheta = int(np.searchsorted(ANGLE_EDGES, theta, side="right") - 1)
                itheta = max(0, min(len(ANGLE_EDGES) - 2, itheta))
                hits[source_i, det] += 1
                angle_den[source_i, st, itheta] += 1
                for iphase, phase in enumerate(PHASES):
                    geom = geometric_cells(row, st, phase)
                    counts, early, lost = electronic_response(
                        u, v, st, phase, response, q
                    )
                    values = (
                        geom,
                        int(np.count_nonzero(counts >= 1)),
                        int(np.count_nonzero(counts >= 5)),
                        int(np.count_nonzero(early >= 1)),
                    )
                    for metric_i, value in enumerate(values):
                        mult[source_i, st, iphase, metric_i, min(value, 5)] += 1
                        module_mult[source_i, det, iphase, metric_i,
                                    min(value, 5)] += 1
                        angle_num[source_i, st, itheta, iphase, metric_i] += value >= 2
                    if iphase == 1:
                        n = n_wires(st)
                        wire_counts[source_i, det, :n, 0] += (counts >= 1)
                        wire_counts[source_i, det, :n, 1] += (counts >= 5)
                        wire_counts[source_i, det, :n, 2] += (early >= 1)
                        edge_electrons[source_i, det] += lost
                if args.max_per_source and index >= args.max_per_source:
                    break
                if index % 100000 == 0:
                    print(
                        f"  {source}: {index:,} crossings in "
                        f"{(time.monotonic() - t0)/60:.1f} min", flush=True
                    )
        if not args.max_per_source:
            for det in LOOKUP:
                if hits[source_i, det] != expected[source, det]:
                    raise ValueError(
                        f"Wrong chamber crossing count {source} DET{det}: "
                        f"{hits[source_i, det]} vs {expected[source, det]}"
                    )
            print(
                f"PASS {source}: every one of 298 module counters agrees",
                flush=True
            )

    with (out / "pi_fixed_wire_counts.csv").open(
        "w", newline="", encoding="utf-8"
    ) as f:
        fields = (
            "source", "station", "det_id", "row", "col", "layer",
            "wire_index", "wire_y_global_cm", "n_wires_in_chamber",
            "chamber_crossings", "reference_expected_crossings",
            "n_ge1e", "n_ge5e", "n_le100ns",
            "P_ge1e", "P_ge5e", "P_le100ns",
            "edge_collected_electrons_lost_from_module",
        )
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for isrc, src in enumerate(SOURCES):
            for det, ch in sorted(LOOKUP.items()):
                n = n_wires(ch.station)
                p = positions(ch.station)
                exposure = int(hits[isrc, det])
                for j in range(n):
                    counts = wire_counts[isrc, det, j, :]
                    writer.writerow({
                        "source": src, "station": ch.station, "det_id": det,
                        "row": ch.row, "col": ch.col, "layer": ch.layer,
                        "wire_index": j,
                        "wire_y_global_cm": ch.y + p[j]/10,
                        "n_wires_in_chamber": n, "chamber_crossings": exposure,
                        "reference_expected_crossings": expected[src, det],
                        "n_ge1e": int(counts[0]),
                        "n_ge5e": int(counts[1]),
                        "n_le100ns": int(counts[2]),
                        "P_ge1e": counts[0]/exposure if exposure else math.nan,
                        "P_ge5e": counts[1]/exposure if exposure else math.nan,
                        "P_le100ns": counts[2]/exposure if exposure else math.nan,
                        "edge_collected_electrons_lost_from_module":
                            int(edge_electrons[isrc, det]),
                    })
    with (out / "pi_fixed_wire_station_multiplicity.csv").open(
        "w", newline="", encoding="utf-8"
    ) as f:
        fields = ("source", "station", "phase_mm", "condition", "crossings",
                  "N0", "N1", "N2", "N3", "N4", "N5plus", "P_ge2", "mean_nwire")
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for isrc, src in enumerate(SOURCES):
            for st in range(6):
                for ip, phase in enumerate(PHASES):
                    for im, metric in enumerate(METRICS):
                        counts = mult[isrc, st, ip, im]
                        n = int(counts.sum())
                        writer.writerow({
                            "source": src, "station": st, "phase_mm": phase,
                            "condition": metric, "crossings": n,
                            **{k: int(counts[j]) for j, k in enumerate(
                                ("N0", "N1", "N2", "N3", "N4", "N5plus")
                            )},
                            "P_ge2": float(counts[2:].sum()/n) if n else math.nan,
                            "mean_nwire": (
                                float(np.dot(counts, np.arange(6))/n)
                                if n else math.nan
                            ),
                        })
    # Exact finite-wire multiplicity per physical gas chamber. Keep all
    # zero-hit modules in output with undefined (NaN) conditional P values.
    # This table is derived from the SAME per-crossing fold as wire_counts,
    # unlike the historical random-phase per-module report.
    with (out / "pi_fixed_wire_module_multiplicity.csv").open(
        "w", newline="", encoding="utf-8"
    ) as f:
        fields = (
            "source", "station", "det_id", "row", "col", "layer",
            "phase_mm", "condition", "crossings", "reference_expected_crossings",
            "N0", "N1", "N2", "N3", "N4", "N5plus",
            "P0", "P1", "P2", "P3", "P4", "P5plus",
            "P_ge2", "low_stats_lt100",
        )
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for isrc, src in enumerate(SOURCES):
            for det, ch in sorted(LOOKUP.items()):
                nhits = int(hits[isrc, det])
                for ip, phase in enumerate(PHASES):
                    for im, metric in enumerate(METRICS):
                        bins = module_mult[isrc, det, ip, im]
                        if int(bins.sum()) != nhits:
                            raise ValueError(
                                f"Module multiplicity sum mismatch "
                                f"{src} DET{det} phase={phase} {metric}"
                            )
                        if nhits and int(bins[1:].sum()) > nhits:
                            raise AssertionError("Invalid chamber multiplicity")
                        writer.writerow({
                            "source": src, "station": ch.station,
                            "det_id": det, "row": ch.row,
                            "col": ch.col, "layer": ch.layer,
                            "phase_mm": phase, "condition": metric,
                            "crossings": nhits,
                            "reference_expected_crossings": expected[src, det],
                            **{f"N{k}": int(bins[j]) for j,k in enumerate(
                                ("0","1","2","3","4","5plus")
                            )},
                            **{f"P{k}": (bins[j]/nhits if nhits else math.nan)
                               for j,k in enumerate(
                                   ("0","1","2","3","4","5plus")
                               )},
                            "P_ge2": bins[2:].sum()/nhits if nhits else math.nan,
                            "low_stats_lt100": int(nhits < 100),
                        })

    with (out / "pi_fixed_wire_angle_response.csv").open(
        "w", newline="", encoding="utf-8"
    ) as f:
        fields = ("source", "station", "theta_abs_low_deg",
                  "theta_abs_high_deg", "phase_mm", "condition",
                  "crossings", "n_ge2", "P_ge2")
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for isrc, src in enumerate(SOURCES):
            for st in range(6):
                for it in range(len(ANGLE_EDGES)-1):
                    n = int(angle_den[isrc, st, it])
                    for ip, phase in enumerate(PHASES):
                        for im, metric in enumerate(METRICS):
                            c = int(angle_num[isrc, st, it, ip, im])
                            writer.writerow({
                                "source": src, "station": st,
                                "theta_abs_low_deg": ANGLE_EDGES[it],
                                "theta_abs_high_deg": ANGLE_EDGES[it + 1],
                                "phase_mm": phase, "condition": metric,
                                "crossings": n, "n_ge2": c,
                                "P_ge2": c/n if n else math.nan,
                            })

    (out / "pi_fixed_wire_settings.json").write_text(
        json.dumps({
            "crossings": str(args.crossings.resolve()),
            "garfield_map": str(args.garfield.resolve()),
            "heed_library": [str(args.heed_tracks.resolve()),
                             str(args.heed_clusters.resolve())],
            "ms2_z_cm": 530, "pitch_mm": PITCH,
            "phase_mm": list(PHASES), "nominal_wires_ms0": 77,
            "nominal_wires_standard": 128, "seed": args.seed,
            "full_statistics": args.max_per_source == 0,
            "limits": [
                "original GEANT plane states projected as straight lines",
                "not updated Geant4 gas hits or Edep",
                "no avalanche gain, electronics threshold or induced signal",
                "finite nominal wire grid is an unconfirmed engineering approximation",
                "garfield B=0 T, symmetric 6mm gas mapped to asymmetric 2+4 mm",
                "concurrent gas volumes are treated as independent ionization",
                "100 ns means individual primary-electron drift, not shaped pulse",
                "wires falling beyond nominal outermost anodes counted as edge loss",
                "phase shifts change wire origins but not nominal number of wires",
            ],
        }, indent=2) + "\n",
        encoding="utf-8",
    )
    print(
        f"Saved per-anode counts, per-module 0/1/2/3/4/5+ "
        f"multiplicity, station and angle summaries under {out}",
        flush=True,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--crossings", type=Path,
        default=Path("test_runs/mwpc_stage_c/checkerboard_full"),
    )
    parser.add_argument(
        "--heed-tracks", type=Path,
        default=Path(
            "test_runs/mwpc_stage_c/c_heed_library_2p4/heed_2p4_tracks.csv"
        ),
    )
    parser.add_argument(
        "--heed-clusters", type=Path,
        default=Path(
            "test_runs/mwpc_stage_c/c_heed_library_2p4/heed_2p4_clusters.csv"
        ),
    )
    parser.add_argument(
        "--garfield", type=Path, default=None,
    )
    parser.add_argument(
        "--out", type=Path,
        default=Path("test_runs/mwpc_stage_c/pi_fixed_wire_maps"),
    )
    parser.add_argument("--seed", type=int, default=314159)
    parser.add_argument(
        "--max-per-source", type=int, default=0,
        help="Smoke mode only; 0 processes all 2.49M chamber crossings",
    )
    args = parser.parse_args()
    if args.max_per_source < 0:
        parser.error("--max-per-source must be nonnegative")
    if args.garfield is None:
        hybrid = Path(
            "test_runs/mwpc_stage_c/c_transport_box_map_hybrid_boundary_2p4/"
            "c_transport_box_electrons.csv"
        )
        baseline = Path(
            "test_runs/mwpc_stage_c/c_transport_box_map_gm2p0_gp4p0/"
            "c_transport_box_electrons.csv"
        )
        args.garfield = hybrid if hybrid.is_file() else baseline
    calc(args)


if __name__ == "__main__":
    main()
