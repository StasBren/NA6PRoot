#!/usr/bin/env python3
"""Read-only, SHA-audited incidence-angle and >=3-anode diagnostics.

Processes the EXISTING completed donor-weighted 100k-per-sample event CSVs.
No Geant4, Heed, Garfield, LUT sampling, or stochastic reweighting.
Outputs new files in a separate directory; frozen comparison is untouched.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path

CHANNELS = ("Jpsi", "Phi", "Omega")
STATIONS = (0, 3, 5)
EDGES = (0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 90, 180)
REQUIRED = {
    "source", "chunk", "event", "track_id", "station", "det_id",
    "theta_u_in_deg", "theta_w_in_deg", "path_length_mm",
    "n_primaries", "n_source_cells", "geometry_n_positive", "lut_n_positive",
    "geometry_n_gt5pct", "lut_n_gt5pct",
    "geometry_q2_fraction", "lut_q2_fraction",
    "geometry_edge_charge", "lut_edge_charge",
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(2**20), b""):
            h.update(chunk)
    return h.hexdigest()


def angle_bin(theta: float) -> tuple[int, int]:
    if not math.isfinite(theta) or not 0 <= theta <= 180:
        raise ValueError(f"Invalid absolute theta_u: {theta}")
    for a, b in zip(EDGES, EDGES[1:]):
        if a <= theta < b or (theta == 180 and b == 180):
            return a, b
    raise AssertionError(f"No angle bin for {theta}")


def add(bucket: dict, row: dict) -> None:
    bucket["n"] += 1
    for method in ("geometry", "lut"):
        n5 = int(row[f"{method}_n_gt5pct"])
        bucket[f"{method}_ge2"] += n5 >= 2
        bucket[f"{method}_ge3"] += n5 >= 3
        bucket[f"{method}_ge4"] += n5 >= 4
    bucket["source_cells_ge3"] += int(row["n_source_cells"]) >= 3


def zero_bucket() -> dict:
    return dict(n=0, geometry_ge2=0, geometry_ge3=0, geometry_ge4=0,
                lut_ge2=0, lut_ge3=0, lut_ge4=0, source_cells_ge3=0)


def write_rows(path: Path, headers: list[str], rows: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as f:
        out = csv.DictWriter(f, fieldnames=headers, extrasaction="ignore")
        out.writeheader()
        out.writerows(rows)


def diagnose(base: Path, output: Path) -> None:
    if output.exists():
        raise FileExistsError(f"Output already exists; refuse overwrite: {output}")
    if not base.is_dir():
        raise FileNotFoundError(base)
    raw = defaultdict(zero_bucket)
    samples = []
    high = []
    third = []
    input_hashes = {}
    n_all = 0
    for channel in CHANNELS:
        for station in STATIONS:
            sample_id = f"{channel}_MS{station}"
            event_path = base / f"events_{sample_id}.csv.gz"
            summary_path = base / f"summary_{sample_id}.json"
            if not event_path.is_file() or not summary_path.is_file():
                raise FileNotFoundError(f"Missing frozen full-run output for {sample_id}")
            summary = json.loads(summary_path.read_text(encoding="utf-8"))
            if (summary["channel"] != channel or summary["station"] != station or
                    summary["n_crossings"] != 100000):
                raise ValueError(f"Expected frozen full 100k summary: {sample_id}")
            input_hashes[event_path.name] = sha256(event_path)
            input_hashes[summary_path.name] = sha256(summary_path)
            counts = zero_bucket()
            with gzip.open(event_path, "rt", encoding="utf-8", newline="") as f:
                rd = csv.DictReader(f)
                if not REQUIRED.issubset(rd.fieldnames or ()):
                    raise ValueError(f"Missing per-muon columns in {event_path}")
                for row in rd:
                    if int(row["station"]) != station:
                        raise ValueError("Per-event station mismatch")
                    theta = abs(float(row["theta_u_in_deg"]))
                    lo, hi = angle_bin(theta)
                    add(raw[("POOLED", -1, lo, hi)], row)
                    add(raw[(channel, station, lo, hi)], row)
                    add(counts, row)
                    if theta >= 30 or (int(row["geometry_n_gt5pct"]) >= 3 or
                                       int(row["lut_n_gt5pct"]) >= 3):
                        entry = {"channel": channel, "station": f"MS{station}",
                                 "abs_theta_u_deg": f"{theta:.7g}",
                                 **{k: row[k] for k in sorted(REQUIRED)}}
                        if theta >= 30:
                            high.append(entry)
                        if (int(row["geometry_n_gt5pct"]) >= 3 or
                                int(row["lut_n_gt5pct"]) >= 3):
                            third.append(entry)
            if counts["n"] != 100000:
                raise ValueError(f"{sample_id}: {counts['n']} records, expected 100000")
            for multiplicity in (2, 3):
                check = summary["threshold_comparison"][f"ge{multiplicity}_anodes_gt5pct"]
                for method in ("geometry", "lut"):
                    observed = counts[f"{method}_ge{multiplicity}"] / counts["n"]
                    if abs(observed - check[f"{method}_rate"]) > 1e-12:
                        raise ValueError(f"{sample_id}: per-event/summary mismatch for {method} >= {multiplicity}")
            samples.append({"channel": channel, "station": f"MS{station}",
                            "n": counts["n"], **{k: v for k, v in counts.items() if k != "n"}})
            n_all += counts["n"]
    if n_all != 900000:
        raise ValueError("Expected 900,000 paired physical muon crossings")
    standard = [b for b in zip(EDGES, EDGES[1:]) if b[1] <= 20]
    denom = sum(raw[("POOLED", -1, *b)]["n"] for b in standard)
    common_weights = {b: raw[("POOLED", -1, *b)]["n"] / denom for b in standard}
    standardized = []
    for channel in CHANNELS:
        for station in STATIONS:
            rec = {"channel": channel, "station": f"MS{station}"}
            for method in ("geometry", "lut"):
                rec[f"{method}_common_angle_0_20"] = sum(
                    common_weights[b] * raw[(channel, station, *b)][f"{method}_ge2"] /
                    raw[(channel, station, *b)]["n"] for b in standard)
            standardized.append(rec)

    bin_rows = []
    for (channel, station, lo, hi), counts in sorted(raw.items()):
        n = counts["n"]
        bin_rows.append({
            "channel": channel, "station": "pooled" if station < 0 else f"MS{station}",
            "theta_lo_deg": lo, "theta_hi_deg": hi, "n": n,
            "geometry_ge2_5pct": counts["geometry_ge2"],
            "lut_ge2_5pct": counts["lut_ge2"],
            "geometry_ge3_5pct": counts["geometry_ge3"],
            "lut_ge3_5pct": counts["lut_ge3"],
            "geometry_ge4_5pct": counts["geometry_ge4"],
            "lut_ge4_5pct": counts["lut_ge4"],
            "n_source_cells_ge3": counts["source_cells_ge3"],
            "geometry_ge2_rate": counts["geometry_ge2"] / n,
            "lut_ge2_rate": counts["lut_ge2"] / n,
            "geometry_ge3_rate": counts["geometry_ge3"] / n,
            "lut_ge3_rate": counts["lut_ge3"] / n,
        })
    output.mkdir(parents=True)
    write_rows(output / "angle_bins_fine.csv", list(bin_rows[0]), bin_rows)
    write_rows(output / "standardized_0_20.csv", list(standardized[0]), standardized)
    write_rows(output / "sample_counts.csv", list(samples[0]), samples)
    columns = ["channel", "station", "abs_theta_u_deg"] + sorted(REQUIRED)
    write_rows(output / "high_angle_muons_ge30.csv", columns, high)
    write_rows(output / "three_plus_anodes_5pct.csv", columns, third)
    manifest = {
        "protocol": "NA60_LUT500_100K_POSTRUN_ANGLE_TAIL_V1",
        "source": "Read-only existing full 9x100k paired output, NO physics rerun",
        "n_crossings": n_all,
        "n_high_angle_ge30": len(high),
        "n_three_plus_union": len(third),
        "physical_anode_significance": "each >5% of total collected anode charge",
        "standardization": "Pooled empirical angle-bin weights, restricted to 0<=|theta_u|<20 degrees; descriptive, not causal",
        "input_sha256": input_hashes,
        "outputs": ["angle_bins_fine.csv", "standardized_0_20.csv",
                    "sample_counts.csv", "high_angle_muons_ge30.csv",
                    "three_plus_anodes_5pct.csv"],
    }
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"PASS: {n_all:,} frozen paired crossings; 3+ union {len(third)}; |theta_u|>=30: {len(high)}")
    print("Fine bins:", output / "angle_bins_fine.csv")
    print("Angle standardized:", output / "standardized_0_20.csv")
    print("Rare 3+ events:", output / "three_plus_anodes_5pct.csv")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--base", type=Path, default=Path.home() /
                   "na6p/studies/lut500_vs_geometry_100k/comparison_100000_donor")
    p.add_argument("--out", type=Path)
    a = p.parse_args()
    diagnose(a.base, a.out or a.base / "angle_tail_diagnostics_v1")


if __name__ == "__main__":
    main()
