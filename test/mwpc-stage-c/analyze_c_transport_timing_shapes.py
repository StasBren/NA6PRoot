#!/usr/bin/env python3
"""
Characterise drift-time distributions in the completed 2+4 mm transport map.

The PI hypothesis to test is whether the microscopic arrival-time distribution
from a small (u,v) box can be represented by a Gaussian-like packet and how
much charge/electron count arrives within electronics-relevant windows such as
100 ns.

Analysis is conditional on collected electrons and, where requested, on the
destination wire.  This matters near the separatrix where a single box can feed
two neighbouring wires and the unconditional timing distribution can become a
mixture rather than a single Gaussian.

Outputs:
  *_boxes.csv              one row per box, all collected destinations combined
  *_box_wires.csv          one row per box and destination wire
  *_representative.csv     a small automatically selected subset for inspection
"""

from __future__ import annotations

import argparse
import csv
import math
import statistics
from collections import defaultdict
from pathlib import Path


def ff(r, k):
    try:
        return float(r[k])
    except (KeyError, TypeError, ValueError):
        return math.nan


def ii(r, k):
    return int(float(r[k]))


def qtile(x, q):
    x = sorted(v for v in x if math.isfinite(v))
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


def moments(times):
    x = [v for v in times if math.isfinite(v)]
    n = len(x)
    if not n:
        return {
            "n": 0, "mean": math.nan, "std": math.nan,
            "skew": math.nan, "excess_kurtosis": math.nan,
            "q16": math.nan, "q50": math.nan, "q84": math.nan,
            "frac25": math.nan, "frac50": math.nan, "frac100": math.nan,
            "frac150": math.nan, "frac200": math.nan,
            "within1": math.nan, "within2": math.nan,
        }
    mu = statistics.fmean(x)
    if n > 1:
        s = statistics.stdev(x)
    else:
        s = 0.0

    if s > 0:
        m2 = sum((v - mu) ** 2 for v in x) / n
        m3 = sum((v - mu) ** 3 for v in x) / n
        m4 = sum((v - mu) ** 4 for v in x) / n
        skew = m3 / (m2 ** 1.5) if m2 > 0 else math.nan
        excess = m4 / (m2 * m2) - 3.0 if m2 > 0 else math.nan
        within1 = sum(abs(v - mu) <= s for v in x) / n
        within2 = sum(abs(v - mu) <= 2 * s for v in x) / n
    else:
        skew = 0.0
        excess = 0.0
        within1 = 1.0
        within2 = 1.0

    return {
        "n": n,
        "mean": mu,
        "std": s,
        "skew": skew,
        "excess_kurtosis": excess,
        "q16": qtile(x, 0.16),
        "q50": qtile(x, 0.50),
        "q84": qtile(x, 0.84),
        "frac25": sum(v <= 25 for v in x) / n,
        "frac50": sum(v <= 50 for v in x) / n,
        "frac100": sum(v <= 100 for v in x) / n,
        "frac150": sum(v <= 150 for v in x) / n,
        "frac200": sum(v <= 200 for v in x) / n,
        "within1": within1,
        "within2": within2,
    }


def row_from(group, wire_label):
    first = group[0]
    collected = [r for r in group if ii(r, "collected") == 1]
    if wire_label != "ALL":
        collected = [r for r in collected if ii(r, "wire_index") == wire_label]
    t = [ff(r, "drift_time_ns") for r in collected]
    m = moments(t)
    return {
        "iu": ii(first, "iu"),
        "iv": ii(first, "iv"),
        "u_center_mm": ff(first, "u_center_mm"),
        "v_center_mm": ff(first, "v_center_mm"),
        "wire_index": wire_label,
        "n_collected": m["n"],
        "mean_ns": m["mean"],
        "std_ns": m["std"],
        "skewness": m["skew"],
        "excess_kurtosis": m["excess_kurtosis"],
        "q16_ns": m["q16"],
        "median_ns": m["q50"],
        "q84_ns": m["q84"],
        "fraction_before_25ns": m["frac25"],
        "fraction_before_50ns": m["frac50"],
        "fraction_before_100ns": m["frac100"],
        "fraction_before_150ns": m["frac150"],
        "fraction_before_200ns": m["frac200"],
        "fraction_within_mean_1sigma": m["within1"],
        "fraction_within_mean_2sigma": m["within2"],
    }


FIELDS = [
    "iu", "iv", "u_center_mm", "v_center_mm", "wire_index",
    "n_collected", "mean_ns", "std_ns", "skewness", "excess_kurtosis",
    "q16_ns", "median_ns", "q84_ns",
    "fraction_before_25ns", "fraction_before_50ns",
    "fraction_before_100ns", "fraction_before_150ns",
    "fraction_before_200ns",
    "fraction_within_mean_1sigma", "fraction_within_mean_2sigma",
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, type=Path)
    ap.add_argument("--output-prefix", required=True, type=Path)
    args = ap.parse_args()

    groups = defaultdict(list)
    with args.input.open(newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            groups[(ii(r, "iu"), ii(r, "iv"))].append(r)

    if not groups:
        raise SystemExit("No transport rows.")

    all_rows = []
    wire_rows = []
    for key in sorted(groups):
        g = groups[key]
        all_rows.append(row_from(g, "ALL"))
        wires = sorted({
            ii(r, "wire_index")
            for r in g
            if ii(r, "collected") == 1
        })
        for w in wires:
            wire_rows.append(row_from(g, w))

    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    boxes_path = Path(str(args.output_prefix) + "_boxes.csv")
    wires_path = Path(str(args.output_prefix) + "_box_wires.csv")
    rep_path = Path(str(args.output_prefix) + "_representative.csv")

    for path, rows in ((boxes_path, all_rows), (wires_path, wire_rows)):
        with path.open("w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=FIELDS)
            w.writeheader()
            w.writerows(rows)

    # Select interpretable cases:
    # central-u boxes near anode and both cathodes, plus cell-boundary boxes.
    targets = [
        (0.1, 0.1, "central near anode"),
        (0.1, -1.9, "central near - cathode"),
        (0.1, 3.9, "central near + cathode"),
        (-1.9, 0.1, "- boundary near anode"),
        (-1.9, -1.9, "- boundary near - cathode"),
        (-1.9, 3.9, "- boundary near + cathode"),
        (1.9, 0.1, "+ boundary near anode"),
        (1.9, -1.9, "+ boundary near - cathode"),
        (1.9, 3.9, "+ boundary near + cathode"),
    ]

    reps = []
    for tu, tv, label in targets:
        candidates = sorted(
            all_rows,
            key=lambda r: (
                abs(r["u_center_mm"] - tu) +
                abs(r["v_center_mm"] - tv)
            ),
        )
        if not candidates:
            continue
        base = dict(candidates[0])
        base["label"] = label
        reps.append(base)

    rep_fields = ["label"] + FIELDS
    with rep_path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=rep_fields)
        w.writeheader()
        w.writerows(reps)

    good = [r for r in all_rows if r["n_collected"] >= 20]
    if good:
        mean_f100 = statistics.fmean(
            r["fraction_before_100ns"] for r in good
        )
        mean_abs_skew = statistics.fmean(abs(r["skewness"]) for r in good)
        mean_abs_kurt = statistics.fmean(
            abs(r["excess_kurtosis"]) for r in good
        )
        mean_cov1 = statistics.fmean(
            r["fraction_within_mean_1sigma"] for r in good
        )
        mean_cov2 = statistics.fmean(
            r["fraction_within_mean_2sigma"] for r in good
        )
    else:
        mean_f100 = mean_abs_skew = mean_abs_kurt = math.nan
        mean_cov1 = mean_cov2 = math.nan

    print("Transport timing-shape diagnostics")
    print(f"  boxes                         : {len(all_rows)}")
    print(f"  conditional box-wire rows     : {len(wire_rows)}")
    print(f"  boxes with >=20 collected e-  : {len(good)}")
    print(f"  mean fraction <=100 ns        : {mean_f100:.6f}")
    print(f"  mean |skewness|               : {mean_abs_skew:.6f}")
    print(f"  mean |excess kurtosis|        : {mean_abs_kurt:.6f}")
    print(
        "  Gaussian reference coverages : "
        "68.27% within 1 sigma, 95.45% within 2 sigma"
    )
    print(f"  observed mean 1-sigma coverage: {mean_cov1:.6f}")
    print(f"  observed mean 2-sigma coverage: {mean_cov2:.6f}")
    print(f"  wrote                         : {boxes_path}")
    print(f"  wrote                         : {wires_path}")
    print(f"  wrote                         : {rep_path}")


if __name__ == "__main__":
    main()
