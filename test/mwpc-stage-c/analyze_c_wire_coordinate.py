#!/usr/bin/env python3
"""
Summarise Stage-C wire-only response scans.

Input is a combined event-summary CSV written by run_c_wire_boundary_scan.sh.
The analysis characterises the coordinate measured by the anode wires:
  - dominant-wire probability;
  - multi-wire probability;
  - neighbour-charge sharing;
  - binary wire-coordinate bias/RMS;
  - charge-weighted wire-CoG bias/RMS;
  - approximate 16/50/84% transition positions at each wire boundary.

No cathode-strip information is used.
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


def i(row, key):
    return int(float(row[key]))


def mean(xs):
    vals = [x for x in xs if math.isfinite(x)]
    return statistics.fmean(vals) if vals else math.nan


def stdev(xs):
    vals = [x for x in xs if math.isfinite(x)]
    return statistics.stdev(vals) if len(vals) >= 2 else 0.0 if vals else math.nan


def rms(xs):
    vals = [x for x in xs if math.isfinite(x)]
    return math.sqrt(statistics.fmean(x * x for x in vals)) if vals else math.nan


def interp_crossing(points, target):
    """Linear interpolation of x where monotonic-ish y crosses target."""
    pts = sorted((x, y) for x, y in points if math.isfinite(x) and math.isfinite(y))
    for (x0, y0), (x1, y1) in zip(pts, pts[1:]):
        if (y0 - target) == 0:
            return x0
        if (y0 - target) * (y1 - target) <= 0 and y1 != y0:
            return x0 + (target - y0) * (x1 - x0) / (y1 - y0)
    return math.nan


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, type=Path)
    ap.add_argument("--output", required=True, type=Path)
    ap.add_argument("--transition-output", type=Path)
    ap.add_argument("--wire-pitch-mm", type=float, default=4.0)
    args = ap.parse_args()

    with args.input.open(newline="", encoding="utf-8") as src:
        rows = list(csv.DictReader(src))
    if not rows:
        raise SystemExit("No rows found.")

    groups = defaultdict(list)
    for row in rows:
        groups[f(row, "u0_mm")].append(row)

    fields = [
        "u0_mm",
        "events",
        "valid_fraction",
        "mean_active_wires",
        "multiwire_fraction",
        "dominant_wire_minus1_fraction",
        "dominant_wire_0_fraction",
        "dominant_wire_plus1_fraction",
        "mean_dominant_fraction",
        "mean_wire_cog_u_mm",
        "std_wire_cog_u_mm",
        "mean_binary_residual_mm",
        "std_binary_residual_mm",
        "rms_binary_residual_mm",
        "mean_wire_cog_residual_mm",
        "std_wire_cog_residual_mm",
        "rms_wire_cog_residual_mm",
        "mean_primary_electrons",
        "mean_avalanche_electrons",
    ]

    out_rows = []
    for u0 in sorted(groups):
        grows = groups[u0]
        valid = [
            r for r in grows
            if i(r, "collected_electrons") > 0
            and math.isfinite(f(r, "dominant_wire_u_mm"))
        ]

        def frac(pred):
            return (
                sum(1 for r in valid if pred(r)) / len(valid)
                if valid else math.nan
            )

        binary_res = [f(r, "binary_residual_mm") for r in valid]
        cog_res = [f(r, "wire_cog_residual_mm") for r in valid]
        wire_cog = [f(r, "wire_cog_u_mm") for r in valid]

        out_rows.append({
            "u0_mm": u0,
            "events": len(grows),
            "valid_fraction": len(valid) / len(grows),
            "mean_active_wires": mean([f(r, "active_wires") for r in valid]),
            "multiwire_fraction": frac(lambda r: i(r, "active_wires") > 1),
            "dominant_wire_minus1_fraction": frac(
                lambda r: i(r, "dominant_wire_index") == -1
            ),
            "dominant_wire_0_fraction": frac(
                lambda r: i(r, "dominant_wire_index") == 0
            ),
            "dominant_wire_plus1_fraction": frac(
                lambda r: i(r, "dominant_wire_index") == 1
            ),
            "mean_dominant_fraction": mean(
                [f(r, "dominant_wire_fraction") for r in valid]
            ),
            "mean_wire_cog_u_mm": mean(wire_cog),
            "std_wire_cog_u_mm": stdev(wire_cog),
            "mean_binary_residual_mm": mean(binary_res),
            "std_binary_residual_mm": stdev(binary_res),
            "rms_binary_residual_mm": rms(binary_res),
            "mean_wire_cog_residual_mm": mean(cog_res),
            "std_wire_cog_residual_mm": stdev(cog_res),
            "rms_wire_cog_residual_mm": rms(cog_res),
            "mean_primary_electrons": mean(
                [f(r, "primary_electrons") for r in valid]
            ),
            "mean_avalanche_electrons": mean(
                [f(r, "avalanche_electrons") for r in valid]
            ),
        })

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as dst:
        w = csv.DictWriter(dst, fieldnames=fields)
        w.writeheader()
        w.writerows(out_rows)

    trans_rows = []
    # Positive boundary: central wire 0 -> +1 around +pitch/2.
    pos = [
        (r["u0_mm"], r["dominant_wire_plus1_fraction"])
        for r in out_rows if r["u0_mm"] > 0
    ]
    # Negative boundary: -1 -> 0 as u increases.  Use neighbour (-1)
    # probability, which falls from ~1 to ~0; reflect x to |u| so it grows.
    neg = [
        (-r["u0_mm"], r["dominant_wire_minus1_fraction"])
        for r in out_rows if r["u0_mm"] < 0
    ]

    for side, pts in (("positive", pos), ("negative_abs", neg)):
        # For negative_abs, sort by |u|. Probability of -1 should increase
        # with |u|, matching the positive convention.
        pts = sorted(pts)
        u16 = interp_crossing(pts, 0.16)
        u50 = interp_crossing(pts, 0.50)
        u84 = interp_crossing(pts, 0.84)
        trans_rows.append({
            "side": side,
            "u16_mm": u16,
            "u50_mm": u50,
            "u84_mm": u84,
            "width_16_84_mm": (
                u84 - u16
                if math.isfinite(u16) and math.isfinite(u84)
                else math.nan
            ),
            "geometric_boundary_mm": 0.5 * args.wire_pitch_mm,
            "u50_minus_geometric_boundary_mm": (
                u50 - 0.5 * args.wire_pitch_mm
                if math.isfinite(u50) else math.nan
            ),
        })

    transition_output = (
        args.transition_output
        if args.transition_output
        else args.output.with_name("c_wire_transition_summary.csv")
    )
    with transition_output.open("w", newline="", encoding="utf-8") as dst:
        fields2 = [
            "side", "u16_mm", "u50_mm", "u84_mm", "width_16_84_mm",
            "geometric_boundary_mm", "u50_minus_geometric_boundary_mm",
        ]
        w = csv.DictWriter(dst, fieldnames=fields2)
        w.writeheader()
        w.writerows(trans_rows)

    print("Stage-C wire-coordinate summary")
    print(f"  events      : {len(rows)}")
    print(f"  u0 points   : {len(out_rows)}")
    print(f"  summary     : {args.output}")
    print(f"  transitions : {transition_output}")
    for r in trans_rows:
        print(
            f"  {r['side']:12s}: "
            f"u50={r['u50_mm']:.4f} mm, "
            f"width16-84={r['width_16_84_mm']:.4f} mm"
        )


if __name__ == "__main__":
    main()
