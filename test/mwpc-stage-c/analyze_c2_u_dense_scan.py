#!/usr/bin/env python3
"""
Summarise a long C2 dense across-wire scan using only the Python standard
library.  The intent is to leave a compact point-by-point file ready for later
physics analysis after an unattended run.
"""

from __future__ import annotations

import argparse
import csv
import math
import statistics
from collections import defaultdict
from pathlib import Path


def f(row, key):
    try:
        return float(row[key])
    except (KeyError, ValueError, TypeError):
        return math.nan


def mean(xs):
    vals = [x for x in xs if math.isfinite(x)]
    return statistics.fmean(vals) if vals else math.nan


def stdev(xs):
    vals = [x for x in xs if math.isfinite(x)]
    return statistics.stdev(vals) if len(vals) >= 2 else 0.0 if vals else math.nan


def frac(rows, predicate):
    if not rows:
        return math.nan
    return sum(1 for r in rows if predicate(r)) / len(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, type=Path)
    ap.add_argument("--output", required=True, type=Path)
    args = ap.parse_args()

    with args.input.open(newline="", encoding="utf-8") as src:
        rows = list(csv.DictReader(src))
    if not rows:
        raise SystemExit("No rows found in combined event summary.")

    groups = defaultdict(list)
    for row in rows:
        groups[f(row, "u0_mm")].append(row)

    fields = [
        "u0_mm",
        "events",
        "valid_fraction",
        "mean_u_reco_mm",
        "std_u_reco_mm",
        "mean_u_reco_minus_u0_mm",
        "std_u_reco_minus_u0_mm",
        "mean_w_reco_mm",
        "std_w_reco_mm",
        "mean_gain_u_mm",
        "std_gain_u_mm",
        "mean_gain_w_mm",
        "std_gain_w_mm",
        "mean_wire_cog_u_mm",
        "std_wire_cog_u_mm",
        "central_wire_fraction",
        "negative_neighbor_fraction",
        "positive_neighbor_fraction",
        "multiwire_fraction",
        "mean_minus_coverage",
        "mean_plus_coverage",
        "mean_primary_electrons",
        "mean_avalanche_electrons",
    ]

    out_rows = []
    for u0 in sorted(groups):
        grows = groups[u0]
        valid = [r for r in grows if r.get("readout_valid") == "1"]

        u_reco = [f(r, "u_stereo_mm") for r in valid]
        w_reco = [f(r, "w_stereo_mm") for r in valid]
        gain_u = [f(r, "gain_weighted_seed_u_mm") for r in valid]
        gain_w = [f(r, "gain_weighted_seed_w_mm") for r in valid]
        wire_u = [f(r, "wire_cog_u_mm") for r in valid]

        # For 4-mm wire pitch, |wire_cog_u| < 1 mm is an unambiguous central
        # wire assignment even if a rare two-wire event yields a weighted CoG.
        central = lambda r: abs(f(r, "wire_cog_u_mm")) < 1.0
        neg = lambda r: f(r, "wire_cog_u_mm") <= -1.0
        pos = lambda r: f(r, "wire_cog_u_mm") >= +1.0
        multi = lambda r: int(float(r.get("active_wires", "0"))) > 1

        out_rows.append({
            "u0_mm": u0,
            "events": len(grows),
            "valid_fraction": len(valid) / len(grows),
            "mean_u_reco_mm": mean(u_reco),
            "std_u_reco_mm": stdev(u_reco),
            "mean_u_reco_minus_u0_mm": mean([x - u0 for x in u_reco]),
            "std_u_reco_minus_u0_mm": stdev([x - u0 for x in u_reco]),
            "mean_w_reco_mm": mean(w_reco),
            "std_w_reco_mm": stdev(w_reco),
            "mean_gain_u_mm": mean(gain_u),
            "std_gain_u_mm": stdev(gain_u),
            "mean_gain_w_mm": mean(gain_w),
            "std_gain_w_mm": stdev(gain_w),
            "mean_wire_cog_u_mm": mean(wire_u),
            "std_wire_cog_u_mm": stdev(wire_u),
            "central_wire_fraction": frac(valid, central),
            "negative_neighbor_fraction": frac(valid, neg),
            "positive_neighbor_fraction": frac(valid, pos),
            "multiwire_fraction": frac(valid, multi),
            "mean_minus_coverage": mean([f(r, "minus_coverage") for r in valid]),
            "mean_plus_coverage": mean([f(r, "plus_coverage") for r in valid]),
            "mean_primary_electrons": mean([f(r, "primary_electrons") for r in valid]),
            "mean_avalanche_electrons": mean([f(r, "avalanche_electrons") for r in valid]),
        })

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as dst:
        writer = csv.DictWriter(dst, fieldnames=fields)
        writer.writeheader()
        writer.writerows(out_rows)

    n_valid = sum(1 for r in rows if r.get("readout_valid") == "1")
    print("Dense C2 scan summary")
    print(f"  events          : {len(rows)}")
    print(f"  valid events    : {n_valid}")
    print(f"  u0 points       : {len(out_rows)}")
    print(f"  wrote           : {args.output}")


if __name__ == "__main__":
    main()
