#!/usr/bin/env python3
"""
Summarise the fast PI cluster->Garfield transport response.

Primary questions:
  * how often more than one wire is active;
  * binary-wire and wire-CoG resolution versus u0 and charge-weighted truth;
  * whether dominant-wire arrival asymmetry carries sub-wire information;
  * channel/station dependence from the realistic MNP33 crossing sample.
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


def finite(vals):
    return [x for x in vals if math.isfinite(x)]


def mean(vals):
    x = finite(vals)
    return statistics.fmean(x) if x else math.nan


def stdev(vals):
    x = finite(vals)
    return statistics.stdev(x) if len(x) > 1 else (0.0 if x else math.nan)


def rms(vals):
    x = finite(vals)
    return math.sqrt(statistics.fmean(v * v for v in x)) if x else math.nan


def corr(a, b):
    pairs = [
        (x, y)
        for x, y in zip(a, b)
        if math.isfinite(x) and math.isfinite(y)
    ]
    if len(pairs) < 2:
        return math.nan
    xa = [x for x, _ in pairs]
    xb = [y for _, y in pairs]
    ma = statistics.fmean(xa)
    mb = statistics.fmean(xb)
    va = sum((x - ma) ** 2 for x in xa)
    vb = sum((y - mb) ** 2 for y in xb)
    if va <= 0 or vb <= 0:
        return math.nan
    return sum((x - ma) * (y - mb) for x, y in pairs) / math.sqrt(va * vb)


def quantile(vals, q):
    x = sorted(finite(vals))
    if not x:
        return math.nan
    if len(x) == 1:
        return x[0]
    z = (len(x) - 1) * q
    lo = int(math.floor(z))
    hi = int(math.ceil(z))
    if lo == hi:
        return x[lo]
    return x[lo] + (z - lo) * (x[hi] - x[lo])


def summary_row(label, rows):
    valid = [r for r in rows if i(r, "collected_electrons") > 0]
    u0 = [f(r, "u0_cell_mm") for r in valid]
    a = [f(r, "dominant_wire_mean_cos_phi") for r in valid]
    alr = [f(r, "dominant_wire_A_LR") for r in valid]
    rb = [f(r, "binary_residual_u0_mm") for r in valid]
    rbc = [f(r, "binary_residual_charge_truth_mm") for r in valid]
    rc = [f(r, "wire_cog_residual_u0_mm") for r in valid]
    rcc = [f(r, "wire_cog_residual_charge_truth_mm") for r in valid]
    return {
        **label,
        "events": len(rows),
        "valid_events": len(valid),
        "valid_fraction": len(valid) / len(rows) if rows else math.nan,
        "mean_primary_electrons": mean(
            [f(r, "primary_electrons") for r in rows]
        ),
        "mean_collection_fraction": mean(
            [f(r, "collection_fraction") for r in valid]
        ),
        "multiwire_fraction": (
            sum(i(r, "active_wires") > 1 for r in valid) / len(valid)
            if valid
            else math.nan
        ),
        "mean_dominant_wire_fraction": mean(
            [f(r, "dominant_wire_fraction") for r in valid]
        ),
        "binary_bias_u0_mm": mean(rb),
        "binary_rms_u0_mm": rms(rb),
        "binary_bias_charge_truth_mm": mean(rbc),
        "binary_rms_charge_truth_mm": rms(rbc),
        "wire_cog_bias_u0_mm": mean(rc),
        "wire_cog_rms_u0_mm": rms(rc),
        "wire_cog_bias_charge_truth_mm": mean(rcc),
        "wire_cog_rms_charge_truth_mm": rms(rcc),
        "mean_drift_time_ns": mean(
            [f(r, "mean_drift_time_ns") for r in valid]
        ),
        "mean_fraction_before_25ns": mean(
            [f(r, "fraction_collected_before_25ns") for r in valid]
        ),
        "mean_fraction_before_50ns": mean(
            [f(r, "fraction_collected_before_50ns") for r in valid]
        ),
        "mean_fraction_before_100ns": mean(
            [f(r, "fraction_collected_before_100ns") for r in valid]
        ),
        "mean_fraction_before_150ns": mean(
            [f(r, "fraction_collected_before_150ns") for r in valid]
        ),
        "mean_fraction_before_200ns": mean(
            [f(r, "fraction_collected_before_200ns") for r in valid]
        ),
        "q05_earliest_drift_time_ns": quantile(
            [f(r, "earliest_drift_time_ns") for r in valid], 0.05
        ),
        "q50_earliest_drift_time_ns": quantile(
            [f(r, "earliest_drift_time_ns") for r in valid], 0.50
        ),
        "q95_earliest_drift_time_ns": quantile(
            [f(r, "earliest_drift_time_ns") for r in valid], 0.95
        ),
        "corr_u0_mean_cos_phi": corr(u0, a),
        "corr_u0_A_LR": corr(u0, alr),
        "mean_abs_cos_phi": mean([abs(x) for x in a]),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--events", required=True, type=Path)
    ap.add_argument("--output-summary", required=True, type=Path)
    ap.add_argument("--output-u-bins", required=True, type=Path)
    ap.add_argument("--u-bin-mm", type=float, default=0.2)
    args = ap.parse_args()

    with args.events.open(newline="", encoding="utf-8") as fobj:
        rows = list(csv.DictReader(fobj))
    if not rows:
        raise SystemExit("No event rows.")

    groups = defaultdict(list)
    for r in rows:
        groups[(r["channel"], i(r, "station"))].append(r)

    fields = [
        "channel", "station", "events", "valid_events", "valid_fraction",
        "mean_primary_electrons", "mean_collection_fraction",
        "multiwire_fraction", "mean_dominant_wire_fraction",
        "binary_bias_u0_mm", "binary_rms_u0_mm",
        "binary_bias_charge_truth_mm", "binary_rms_charge_truth_mm",
        "wire_cog_bias_u0_mm", "wire_cog_rms_u0_mm",
        "wire_cog_bias_charge_truth_mm", "wire_cog_rms_charge_truth_mm",
        "mean_drift_time_ns",
        "mean_fraction_before_25ns", "mean_fraction_before_50ns",
        "mean_fraction_before_100ns", "mean_fraction_before_150ns",
        "mean_fraction_before_200ns",
        "q05_earliest_drift_time_ns", "q50_earliest_drift_time_ns",
        "q95_earliest_drift_time_ns",
        "corr_u0_mean_cos_phi", "corr_u0_A_LR", "mean_abs_cos_phi",
    ]

    out = [summary_row({"channel": "ALL", "station": -1}, rows)]
    for (channel, station), grows in sorted(groups.items()):
        out.append(
            summary_row(
                {"channel": channel, "station": station},
                grows,
            )
        )

    args.output_summary.parent.mkdir(parents=True, exist_ok=True)
    with args.output_summary.open("w", newline="", encoding="utf-8") as fobj:
        w = csv.DictWriter(fobj, fieldnames=fields)
        w.writeheader()
        w.writerows(out)

    # Bin within one wire pitch to expose the response curve versus true phase.
    u_groups = defaultdict(list)
    for r in rows:
        u = f(r, "u0_cell_mm")
        if not math.isfinite(u):
            continue
        ib = math.floor((u + 2.0) / args.u_bin_mm)
        u_groups[ib].append(r)

    ubin_fields = [
        "u_bin", "u_center_mm", "events", "valid_events",
        "mean_active_wires", "multiwire_fraction",
        "mean_dominant_wire_index", "mean_wire_cog_u_mm",
        "mean_cos_phi", "std_cos_phi", "mean_A_LR", "std_A_LR",
        "binary_bias_u0_mm", "binary_rms_u0_mm",
        "wire_cog_bias_u0_mm", "wire_cog_rms_u0_mm",
        "mean_collection_fraction", "mean_drift_time_ns",
        "mean_fraction_before_25ns", "mean_fraction_before_50ns",
        "mean_fraction_before_100ns", "mean_fraction_before_150ns",
        "mean_fraction_before_200ns",
    ]
    ubin_rows = []
    for ib, grows in sorted(u_groups.items()):
        valid = [r for r in grows if i(r, "collected_electrons") > 0]
        cosv = [f(r, "dominant_wire_mean_cos_phi") for r in valid]
        alr = [f(r, "dominant_wire_A_LR") for r in valid]
        rb = [f(r, "binary_residual_u0_mm") for r in valid]
        rc = [f(r, "wire_cog_residual_u0_mm") for r in valid]
        ubin_rows.append({
            "u_bin": ib,
            "u_center_mm": -2.0 + (ib + 0.5) * args.u_bin_mm,
            "events": len(grows),
            "valid_events": len(valid),
            "mean_active_wires": mean(
                [f(r, "active_wires") for r in valid]
            ),
            "multiwire_fraction": (
                sum(i(r, "active_wires") > 1 for r in valid) / len(valid)
                if valid
                else math.nan
            ),
            "mean_dominant_wire_index": mean(
                [f(r, "dominant_wire_index") for r in valid]
            ),
            "mean_wire_cog_u_mm": mean(
                [f(r, "wire_cog_u_mm") for r in valid]
            ),
            "mean_cos_phi": mean(cosv),
            "std_cos_phi": stdev(cosv),
            "mean_A_LR": mean(alr),
            "std_A_LR": stdev(alr),
            "binary_bias_u0_mm": mean(rb),
            "binary_rms_u0_mm": rms(rb),
            "wire_cog_bias_u0_mm": mean(rc),
            "wire_cog_rms_u0_mm": rms(rc),
            "mean_collection_fraction": mean(
                [f(r, "collection_fraction") for r in valid]
            ),
            "mean_drift_time_ns": mean(
                [f(r, "mean_drift_time_ns") for r in valid]
            ),
            "mean_fraction_before_25ns": mean(
                [f(r, "fraction_collected_before_25ns") for r in valid]
            ),
            "mean_fraction_before_50ns": mean(
                [f(r, "fraction_collected_before_50ns") for r in valid]
            ),
            "mean_fraction_before_100ns": mean(
                [f(r, "fraction_collected_before_100ns") for r in valid]
            ),
            "mean_fraction_before_150ns": mean(
                [f(r, "fraction_collected_before_150ns") for r in valid]
            ),
            "mean_fraction_before_200ns": mean(
                [f(r, "fraction_collected_before_200ns") for r in valid]
            ),
        })

    with args.output_u_bins.open("w", newline="", encoding="utf-8") as fobj:
        w = csv.DictWriter(fobj, fieldnames=ubin_fields)
        w.writeheader()
        w.writerows(ubin_rows)

    allrow = out[0]
    print("PI fast-response summary")
    print(f"  events                    : {allrow['events']}")
    print(f"  valid fraction            : {allrow['valid_fraction']:.6f}")
    print(f"  multiwire fraction        : {allrow['multiwire_fraction']:.6f}")
    print(
        f"  binary RMS vs u0 [mm]     : "
        f"{allrow['binary_rms_u0_mm']:.6f}"
    )
    print(
        f"  wire-CoG RMS vs u0 [mm]   : "
        f"{allrow['wire_cog_rms_u0_mm']:.6f}"
    )
    print(
        f"  corr(u0,<cos phi>)        : "
        f"{allrow['corr_u0_mean_cos_phi']:.6f}"
    )
    print(
        f"  corr(u0,A_LR)             : "
        f"{allrow['corr_u0_A_LR']:.6f}"
    )
    print(f"  wrote summary             : {args.output_summary}")
    print(f"  wrote u bins              : {args.output_u_bins}")


if __name__ == "__main__":
    main()
