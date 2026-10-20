#!/usr/bin/env python3
"""PI: validate centered Garfield 4mm LUT sample quality + in-memory lookup rate.

Uses the already frozen B=0, 3+3mm p=4mm native Garfield LUT; no new physics
simulations. Uncertainty in each box is from independently seeded microscopic
Garfield trials, NOT a muon-event uncertainty. Timing deliberately excludes
Geant4, Heed, disk I/O, Garfield construction and Python interpreter startup.
"""
from __future__ import annotations
import argparse
import csv
import json
import math
import resource
import time
from pathlib import Path

import numpy as np
from pi_centered_electron_sampler import load_centered_map
from pi_stagec_performance_report import stats

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
BASE = REPO / "test_runs/mwpc_stage_c"
MAP_DIR = BASE / "c_transport_box_map_gm3p0_gp3p0"


def assess(pool, bootstrap=200):
    rows = []
    for bid in range(pool.shape[0]):
        hit = pool[bid]
        collected = hit[:, 0] == 1
        valid_time = np.isfinite(hit[:, 2]) & (hit[:, 2] >= 0)
        good = hit[collected & valid_time, 2]
        measured = stats(good, seed=180019+bid, replicates=bootstrap)
        # Symmetric splitting into independent sample-index halves diagnoses
        # sensitivity to finite trials; no fabricated event correlations.
        half = len(hit) // 2
        first = hit[:half]
        second = hit[half:]
        ta = first[(first[:, 0] == 1) &
                   np.isfinite(first[:, 2]) & (first[:, 2] >= 0), 2]
        tb = second[(second[:, 0] == 1) &
                    np.isfinite(second[:, 2]) & (second[:, 2] >= 0), 2]
        mean_delta = (float(abs(np.mean(ta)-np.mean(tb)))
                      if len(ta) and len(tb) else None)
        rows.append({
            "box_id": bid, "n_trials": len(hit),
            "n_collected": int(np.count_nonzero(collected)),
            "n_collected_with_known_time": len(good),
            "n_collected_unknown_time":
                int(np.count_nonzero(collected))-len(good),
            "collection_fraction": float(np.mean(collected)),
            "mean_drift_ns": measured["mean"],
            "rms_drift_ns": measured["rms"],
            "mean_stat_error_ns": measured["mean_error"],
            "rms_stat_error_ns": measured["rms_error"],
            "abs_split_half_mean_delta_ns": mean_delta,
            "insufficient_collected_n_lt_10": len(good) < 10,
        })
    return rows


def microbenchmark(pool, n_electrons=300000, batch=15000):
    if n_electrons <= 0 or batch <= 0:
        raise ValueError("Positive sample and batch sizes required")
    rng = np.random.default_rng(20261020)
    lut_ns_per_lookup = []
    geom_ns_per_lookup = []
    observed = 0
    while observed < n_electrons:
        n = min(batch, n_electrons-observed)
        u = rng.uniform(-1.999, 1.999, n)
        v = rng.uniform(-2.999, 2.999, n)
        sample = rng.integers(0, pool.shape[1], n)
        start = time.perf_counter()
        iu = np.clip(np.floor((u+2.)/.2).astype(np.int64), 0, 19)
        iv = np.clip(np.floor((v+3.)/.2).astype(np.int64), 0, 29)
        response = pool[iu*30+iv, sample]
        # Consume response to include random access (collection/wire/time).
        n_collected = int(np.count_nonzero(response[:,0]))
        time_lut = time.perf_counter() - start
        if n_collected < 0:
            raise AssertionError("Unreachable")
        start = time.perf_counter()
        # GEOMETRY ONLY: nearest nominal anode, no drift/loss/time physics.
        wire = np.floor(u/4. + 0.5).astype(np.int64)
        n_geometric = int(np.count_nonzero(wire == 0))
        time_geometric = time.perf_counter() - start
        if n_geometric != n:
            raise AssertionError("Unexpected nominal centered wire cell")
        lut_ns_per_lookup.append(time_lut/n*1.e9)
        geom_ns_per_lookup.append(time_geometric/n*1.e9)
        observed += n
    return {
        "n_synthetic_positions": observed,
        "batches": len(lut_ns_per_lookup),
        "basis": ("Nanoseconds per electron, vectorized in-memory NumPy "
                  "lookup ONLY, excluding sampling of starting points, "
                  "Heed, Geant4, full Garfield, disk I/O."),
        "lut_nanoseconds_per_electron":
            stats(lut_ns_per_lookup, seed=20261020, replicates=1000),
        "geometry_nanoseconds_per_electron":
            stats(geom_ns_per_lookup, seed=20261021, replicates=1000),
        "geometry_limitations":
            "Nearest wire ID only: NOT a physics-equivalent replacement "
            "for Garfield drift, collection, time, loss, or field-wire effects.",
        "peak_process_maxrss_kb_linux": resource.getrusage(
            resource.RUSAGE_SELF).ru_maxrss,
    }


def run(lut, manifest, out, n_electrons=300000, batch=15000, bootstrap=200):
    started = time.perf_counter()
    pool, meta = load_centered_map(lut, manifest)
    seconds_loading = time.perf_counter() - started
    if pool.ndim != 3 or pool.shape[0] != 600 or pool.shape[2] != 3:
        raise ValueError("Unexpected centered lookup array shape")
    out.mkdir(parents=True, exist_ok=True)
    rows = assess(pool, bootstrap=bootstrap)
    with (out / "pi_centered_lut_box_quality.csv").open(
            "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    performance = microbenchmark(pool, n_electrons, batch)
    report = {
        "geometry_id": meta["geometry_id"],
        "map_sha256": meta["response_csv_sha256"],
        "electron_trials_per_box": pool.shape[1],
        "number_of_boxes": pool.shape[0],
        "load_wall_seconds": seconds_loading,
        "boxes_with_less_than_10_collected_trials":
            sum(row["insufficient_collected_n_lt_10"] for row in rows),
        "total_garfield_trials": sum(row["n_trials"] for row in rows),
        "total_collected_trials": sum(row["n_collected"] for row in rows),
        "total_collected_with_known_time":
            sum(row["n_collected_with_known_time"] for row in rows),
        "mean_collection_fraction_uniform_boxes":
            sum(row["collection_fraction"] for row in rows)/len(rows),
        "in_memory_benchmark": performance,
        "statistics": (
            "Box time means/RMS and bootstrap errors condition on collection. "
            "Samples are Garfield single-electron trials, not muon events. "
            "A uniform-box mean is NOT a real muon-angle-weighted mean. "
            "Validation against NEW independent direct Garfield trajectories "
            "is still required to measure LUT interpolation bias."),
    }
    path = out / "pi_centered_lut_benchmark.json"
    path.write_text(json.dumps(report, indent=2) + "\n",
                    encoding="utf-8")
    print(f"PASS: {len(rows)} LUT boxes, {report['total_garfield_trials']} "
          f"single-electron trials; summary: {path}")
    print(f"Low-stat collected bins (<10): "
          f"{report['boxes_with_less_than_10_collected_trials']}")
    print("NOTICE: this is LUT precision + lookup throughput, not direct "
          "Garfield physics equivalence or full pipeline throughput.")
    return report


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--map", type=Path,
                   default=MAP_DIR/"c_transport_box_electrons.csv")
    p.add_argument("--manifest", type=Path,
                   default=MAP_DIR/"pi_centered_transport_manifest.json")
    p.add_argument("--out", type=Path,
                   default=BASE/"pi_centered_lut_benchmark")
    p.add_argument("--electrons", type=int, default=300000)
    p.add_argument("--batch", type=int, default=15000)
    p.add_argument("--bootstrap", type=int, default=200)
    args = p.parse_args()
    if args.bootstrap < 2:
        p.error("bootstrap must be >= 2")
    run(args.map, args.manifest, args.out, args.electrons,
        args.batch, args.bootstrap)


if __name__ == "__main__":
    main()
