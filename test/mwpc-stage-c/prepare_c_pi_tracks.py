#!/usr/bin/env python3
"""
Prepare realistic MWPC track-crossing distributions for the PI response-map
study from the large J/psi, phi and omega MNP33 bore scans.

Input files are the existing bore_muon_crossings.csv files with columns:
  job,event,track,pdg,station,x_cm,y_cm,z_cm,
  px_GeV,py_GeV,pz_GeV,p_GeV,theta_x_deg,theta_y_deg

Coordinate mapping used by the MWPC Stage-C convention:
  global X -> local w (along wires)
  global Y -> local u (across wires)
  global Z -> local v (chamber normal)

Therefore:
  theta_w = atan(px/pz) = theta_x
  theta_u = atan(py/pz) = theta_y

The absolute local wire-cell phase is intentionally NOT inferred from y_cm.
The current chamber layout does not define a frozen wire-origin phase for the
old bore-crossing files.  A later synthetic-cluster step samples u0 uniformly
inside one 4-mm wire pitch (or can use an explicitly supplied phase).
"""

from __future__ import annotations

import argparse
import csv
import math
import random
import statistics
from collections import defaultdict
from pathlib import Path


CHANNELS = ("jpsi40-bore", "phi40-bore", "omega40-bore")


def ff(row, key):
    try:
        return float(row[key])
    except (KeyError, TypeError, ValueError):
        return math.nan


def ii(row, key):
    return int(float(row[key]))


def quantile(values, q):
    x = sorted(v for v in values if math.isfinite(v))
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


def mean(values):
    x = [v for v in values if math.isfinite(v)]
    return statistics.fmean(x) if x else math.nan


def stdev(values):
    x = [v for v in values if math.isfinite(v)]
    return statistics.stdev(x) if len(x) > 1 else (0.0 if x else math.nan)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--base",
        required=True,
        type=Path,
        help=".../mnp33-bore-runs/<run-id> directory",
    )
    ap.add_argument(
        "--output-dir",
        type=Path,
        default=Path("test_runs/mwpc_stage_c/c_pi_tracks"),
    )
    ap.add_argument(
        "--sample-per-channel-station",
        type=int,
        default=2000,
        help="Rows retained per channel/station in compact sample; 0 = all",
    )
    ap.add_argument("--seed", type=int, default=180001)
    args = ap.parse_args()

    rng = random.Random(args.seed)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    groups = defaultdict(list)
    total_rows = 0

    for channel in CHANNELS:
        path = args.base / "analysis" / channel / "bore_muon_crossings.csv"
        if not path.exists():
            raise SystemExit(f"Missing input: {path}")

        with path.open(newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            required = {
                "job", "event", "track", "pdg", "station",
                "x_cm", "y_cm", "z_cm",
                "px_GeV", "py_GeV", "pz_GeV", "p_GeV",
                "theta_x_deg", "theta_y_deg",
            }
            missing = required - set(reader.fieldnames or [])
            if missing:
                raise SystemExit(
                    f"{path} missing columns: {sorted(missing)}"
                )

            for row in reader:
                station = ii(row, "station")
                px = ff(row, "px_GeV")
                py = ff(row, "py_GeV")
                pz = ff(row, "pz_GeV")
                p = ff(row, "p_GeV")

                theta_w = math.degrees(math.atan2(px, pz))
                theta_u = math.degrees(math.atan2(py, pz))
                # Keep both recomputed slopes and the old precomputed columns
                # so a convention mismatch is immediately visible.
                rec = {
                    "channel": channel,
                    "job": row["job"],
                    "event": ii(row, "event"),
                    "track": ii(row, "track"),
                    "pdg": ii(row, "pdg"),
                    "station": station,
                    "x_cm": ff(row, "x_cm"),
                    "y_cm": ff(row, "y_cm"),
                    "z_cm": ff(row, "z_cm"),
                    "px_GeV": px,
                    "py_GeV": py,
                    "pz_GeV": pz,
                    "p_GeV": p,
                    "theta_u_deg": theta_u,
                    "theta_w_deg": theta_w,
                    "theta_y_input_deg": ff(row, "theta_y_deg"),
                    "theta_x_input_deg": ff(row, "theta_x_deg"),
                    "tan_theta_u": py / pz if pz != 0 else math.nan,
                    "tan_theta_w": px / pz if pz != 0 else math.nan,
                }
                groups[(channel, station)].append(rec)
                total_rows += 1

    summary_fields = [
        "channel", "station", "crossings",
        "mean_p_GeV", "std_p_GeV", "q05_p_GeV",
        "q50_p_GeV", "q95_p_GeV",
        "mean_theta_u_deg", "std_theta_u_deg",
        "q05_theta_u_deg", "q50_theta_u_deg", "q95_theta_u_deg",
        "mean_theta_w_deg", "std_theta_w_deg",
        "q05_theta_w_deg", "q50_theta_w_deg", "q95_theta_w_deg",
        "mean_abs_theta_u_deg", "mean_abs_theta_w_deg",
        "max_abs_theta_input_recalc_diff_deg",
    ]
    summary_rows = []

    for (channel, station), rows in sorted(groups.items()):
        p = [r["p_GeV"] for r in rows]
        tu = [r["theta_u_deg"] for r in rows]
        tw = [r["theta_w_deg"] for r in rows]
        diffs = [
            max(
                abs(r["theta_u_deg"] - r["theta_y_input_deg"]),
                abs(r["theta_w_deg"] - r["theta_x_input_deg"]),
            )
            for r in rows
        ]
        summary_rows.append({
            "channel": channel,
            "station": station,
            "crossings": len(rows),
            "mean_p_GeV": mean(p),
            "std_p_GeV": stdev(p),
            "q05_p_GeV": quantile(p, 0.05),
            "q50_p_GeV": quantile(p, 0.50),
            "q95_p_GeV": quantile(p, 0.95),
            "mean_theta_u_deg": mean(tu),
            "std_theta_u_deg": stdev(tu),
            "q05_theta_u_deg": quantile(tu, 0.05),
            "q50_theta_u_deg": quantile(tu, 0.50),
            "q95_theta_u_deg": quantile(tu, 0.95),
            "mean_theta_w_deg": mean(tw),
            "std_theta_w_deg": stdev(tw),
            "q05_theta_w_deg": quantile(tw, 0.05),
            "q50_theta_w_deg": quantile(tw, 0.50),
            "q95_theta_w_deg": quantile(tw, 0.95),
            "mean_abs_theta_u_deg": mean(abs(x) for x in tu),
            "mean_abs_theta_w_deg": mean(abs(x) for x in tw),
            "max_abs_theta_input_recalc_diff_deg": max(diffs) if diffs else math.nan,
        })

    summary_path = args.output_dir / "c_pi_track_distribution_summary.csv"
    with summary_path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=summary_fields)
        w.writeheader()
        w.writerows(summary_rows)

    sample_fields = [
        "channel", "job", "event", "track", "pdg", "station",
        "x_cm", "y_cm", "z_cm",
        "px_GeV", "py_GeV", "pz_GeV", "p_GeV",
        "theta_u_deg", "theta_w_deg", "tan_theta_u", "tan_theta_w",
    ]
    sample_path = args.output_dir / "c_pi_track_sample.csv"
    n_sample = 0
    with sample_path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=sample_fields)
        w.writeheader()
        for key in sorted(groups):
            rows = groups[key]
            if args.sample_per_channel_station > 0:
                take = min(args.sample_per_channel_station, len(rows))
                selected = rng.sample(rows, take)
            else:
                selected = rows
            for r in selected:
                w.writerow({k: r[k] for k in sample_fields})
                n_sample += 1

    print("PI realistic MWPC track preparation")
    print(f"  source crossings : {total_rows}")
    print(f"  channel/stations : {len(groups)}")
    print(f"  compact sample   : {n_sample}")
    print(f"  summary          : {summary_path}")
    print(f"  sample           : {sample_path}")
    print()
    print("Coordinate convention:")
    print("  theta_u = theta_y = atan(py/pz)  [across wires]")
    print("  theta_w = theta_x = atan(px/pz)  [along wires]")
    print("  absolute y_cm is NOT folded into a wire cell in this step.")


if __name__ == "__main__":
    main()
