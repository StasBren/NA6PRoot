#!/usr/bin/env python3
"""
Summarise the Stage-C single-electron transport box map.

For each (u,v) gas box, calculate:
  * collection efficiency;
  * probability to end on wire -1, 0, +1;
  * drift-time mean/std/median and central 68% interval;
  * hit-side asymmetry A_LR = (N(+u side)-N(-u side))/N_collected;
  * angular moment <cos(phi_hit)> (smooth left/right asymmetry);
  * <sin(phi_hit)> (top/bottom component);
  * incidence-angle mean/std.

The resulting CSV is the transport response kernel that can later be weighted
by a GEANT cluster distribution.
"""

from __future__ import annotations

import argparse
import csv
import math
import statistics
from collections import defaultdict
from pathlib import Path


def ff(row, key):
    try:
        return float(row[key])
    except (KeyError, TypeError, ValueError):
        return math.nan


def ii(row, key):
    try:
        return int(float(row[key]))
    except (KeyError, TypeError, ValueError):
        return 0


def finite(xs):
    return [x for x in xs if math.isfinite(x)]


def mean(xs):
    v = finite(xs)
    return statistics.fmean(v) if v else math.nan


def stdev(xs):
    v = finite(xs)
    return statistics.stdev(v) if len(v) > 1 else (0.0 if v else math.nan)


def quantile(xs, q):
    v = sorted(finite(xs))
    if not v:
        return math.nan
    if len(v) == 1:
        return v[0]
    x = (len(v) - 1) * q
    lo = int(math.floor(x))
    hi = int(math.ceil(x))
    if lo == hi:
        return v[lo]
    return v[lo] + (x - lo) * (v[hi] - v[lo])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, type=Path)
    ap.add_argument("--output", required=True, type=Path)
    args = ap.parse_args()

    with args.input.open(newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    if not rows:
        raise SystemExit("No electron rows found.")

    groups = defaultdict(list)
    for r in rows:
        groups[(ii(r, "iu"), ii(r, "iv"))].append(r)

    fields = [
        "iu", "iv", "box_id", "u_center_mm", "v_center_mm",
        "electrons_launched", "electrons_collected", "collection_fraction",
        "p_wire_minus1", "p_wire_0", "p_wire_plus1", "p_other_wire",
        "mean_drift_time_ns", "std_drift_time_ns",
        "q16_drift_time_ns", "median_drift_time_ns", "q84_drift_time_ns",
        "A_LR", "mean_cos_phi_hit", "std_cos_phi_hit",
        "mean_sin_phi_hit", "std_sin_phi_hit",
        "mean_phi_hit_deg", "std_phi_hit_deg",
        "mean_incidence_normal_deg", "std_incidence_normal_deg",
        "A_LR_wire_minus1", "A_LR_wire_0", "A_LR_wire_plus1",
        "mean_cos_phi_wire_minus1", "mean_cos_phi_wire_0",
        "mean_cos_phi_wire_plus1",
        "mean_incidence_wire_minus1_deg",
        "mean_incidence_wire_0_deg",
        "mean_incidence_wire_plus1_deg",
    ]

    out_rows = []
    for (iu, iv), grows in sorted(groups.items()):
        collected = [r for r in grows if ii(r, "collected") == 1]
        n = len(grows)
        nc = len(collected)

        def p_wire(k):
            return (
                sum(ii(r, "wire_index") == k for r in collected) / nc
                if nc else math.nan
            )

        p_m1 = p_wire(-1)
        p_0 = p_wire(0)
        p_p1 = p_wire(+1)
        if nc:
            n_known = sum(
                ii(r, "wire_index") in (-1, 0, 1)
                for r in collected
            )
            p_other = 1.0 - n_known / nc
        else:
            p_other = math.nan

        sides = [ii(r, "hit_side_u") for r in collected]
        n_plus = sum(s > 0 for s in sides)
        n_minus = sum(s < 0 for s in sides)
        denom = n_plus + n_minus
        a_lr = (
            (n_plus - n_minus) / denom
            if denom else math.nan
        )

        drift = [ff(r, "drift_time_ns") for r in collected]
        cosphi = [ff(r, "cos_phi_hit") for r in collected]
        sinphi = [ff(r, "sin_phi_hit") for r in collected]
        phi = [ff(r, "phi_hit_deg") for r in collected]
        inc = [ff(r, "incidence_normal_deg") for r in collected]

        conditional = {}
        for wire in (-1, 0, +1):
            wr = [r for r in collected if ii(r, "wire_index") == wire]
            wsides = [ii(r, "hit_side_u") for r in wr]
            wp = sum(s > 0 for s in wsides)
            wm = sum(s < 0 for s in wsides)
            conditional[wire] = {
                "A": (wp - wm) / (wp + wm) if (wp + wm) else math.nan,
                "cos": mean([ff(r, "cos_phi_hit") for r in wr]),
                "inc": mean([ff(r, "incidence_normal_deg") for r in wr]),
            }

        first = grows[0]
        out_rows.append({
            "iu": iu,
            "iv": iv,
            "box_id": ii(first, "box_id"),
            "u_center_mm": ff(first, "u_center_mm"),
            "v_center_mm": ff(first, "v_center_mm"),
            "electrons_launched": n,
            "electrons_collected": nc,
            "collection_fraction": nc / n if n else math.nan,
            "p_wire_minus1": p_m1,
            "p_wire_0": p_0,
            "p_wire_plus1": p_p1,
            "p_other_wire": p_other,
            "mean_drift_time_ns": mean(drift),
            "std_drift_time_ns": stdev(drift),
            "q16_drift_time_ns": quantile(drift, 0.16),
            "median_drift_time_ns": quantile(drift, 0.50),
            "q84_drift_time_ns": quantile(drift, 0.84),
            "A_LR": a_lr,
            "mean_cos_phi_hit": mean(cosphi),
            "std_cos_phi_hit": stdev(cosphi),
            "mean_sin_phi_hit": mean(sinphi),
            "std_sin_phi_hit": stdev(sinphi),
            "mean_phi_hit_deg": mean(phi),
            "std_phi_hit_deg": stdev(phi),
            "mean_incidence_normal_deg": mean(inc),
            "std_incidence_normal_deg": stdev(inc),
            "A_LR_wire_minus1": conditional[-1]["A"],
            "A_LR_wire_0": conditional[0]["A"],
            "A_LR_wire_plus1": conditional[+1]["A"],
            "mean_cos_phi_wire_minus1": conditional[-1]["cos"],
            "mean_cos_phi_wire_0": conditional[0]["cos"],
            "mean_cos_phi_wire_plus1": conditional[+1]["cos"],
            "mean_incidence_wire_minus1_deg": conditional[-1]["inc"],
            "mean_incidence_wire_0_deg": conditional[0]["inc"],
            "mean_incidence_wire_plus1_deg": conditional[+1]["inc"],
        })

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(out_rows)

    print("Stage-C transport box-map summary")
    print(f"  raw electrons : {len(rows)}")
    print(f"  boxes         : {len(out_rows)}")
    print(
        f"  collected     : "
        f"{sum(r['electrons_collected'] for r in out_rows)}"
    )
    print(f"  wrote         : {args.output}")


if __name__ == "__main__":
    main()
