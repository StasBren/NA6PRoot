#!/usr/bin/env python3
"""Reconstruct Stage-C Garfield (u,v) electron drift-time heat maps.

Reads completed single-electron transport-box sample (hybrid or baseline).
For each launch box, reports collection probability, mean/median collected
electron drift time, and fraction of ALL launched electrons collected by 100ns.
Time statistics exclude lost electrons by definition; missing bins are white.
"""
from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def analyze(path):
    time_groups = defaultdict(list)
    started = np.zeros((30, 20), dtype=np.int64)
    collected = np.zeros((30, 20), dtype=np.int64)
    in100 = np.zeros((30, 20), dtype=np.int64)
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        required = {"iu", "iv", "u_center_mm", "v_center_mm",
                    "collected", "drift_time_ns"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(f"Missing fields {required-set(reader.fieldnames or [])}")
        for row in reader:
            iu, iv = int(float(row["iu"])), int(float(row["iv"]))
            if not (0 <= iu < 20 and 0 <= iv < 30):
                raise ValueError(f"Unexpected u/v box {iu},{iv}")
            if abs(float(row["u_center_mm"]) - (-1.9 + .2*iu)) > 1e-5:
                raise ValueError("Wrong u centre")
            if abs(float(row["v_center_mm"]) - (-1.9 + .2*iv)) > 1e-5:
                raise ValueError("Wrong v centre")
            started[iv, iu] += 1
            if int(float(row["collected"])) == 1:
                collected[iv, iu] += 1
                t = float(row["drift_time_ns"])
                if math.isfinite(t) and t >= 0.:
                    time_groups[(iv, iu)].append(t)
                    if t <= 100.:
                        in100[iv, iu] += 1
    if np.any(started == 0):
        raise ValueError(f"Incomplete Garfield map: {np.count_nonzero(started==0)} boxes empty")
    mean = np.full((30, 20), np.nan)
    med = np.full((30, 20), np.nan)
    for (iv, iu), times in time_groups.items():
        mean[iv, iu] = np.mean(times)
        med[iv, iu] = np.median(times)
    efficiency = collected / started
    in100_fraction = in100 / started
    return started, collected, mean, med, efficiency, in100_fraction


def save(args):
    inputs = analyze(args.input)
    started, collected, mean, med, efficiency, in100 = inputs
    args.out.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 2, figsize=(13, 10), layout="constrained")
    specs = (
        ("Mean collected-electron drift time", mean, "ns", "viridis"),
        ("Median collected-electron drift time", med, "ns", "viridis"),
        ("Electron collection fraction", efficiency, "fraction", "magma"),
        ("Collected by 100 ns / all launched", in100, "fraction", "magma"),
    )
    for ax, (name, values, label, color) in zip(axes.flat, specs):
        cmap = plt.get_cmap(color).copy()
        cmap.set_bad("white")
        image = ax.imshow(
            np.ma.masked_invalid(values), origin="lower", aspect="equal",
            interpolation="nearest", cmap=cmap, extent=[-2, 2, -2, 4],
            vmin=0 if label == "fraction" else None,
            vmax=1 if label == "fraction" else None,
        )
        fig.colorbar(image, ax=ax, label=label, shrink=.78)
        ax.plot(0, 0, marker="o", ms=4, color="white",
                mec="black", zorder=3)
        ax.set_xlabel("Local u across wires [mm]")
        ax.set_ylabel("Local v, chamber normal [mm]")
        ax.set_title(name)
    fig.suptitle(
        "Stage C Garfield 20 x 30 box grid | B = 0 T\n"
        "Ar/CO2 70:30 | 6-mm 2+4 gap | 4-mm anode pitch | 1800 V\n"
        f"{int(started.sum()):,} launched electrons",
        fontsize=12,
    )
    target = args.out / "stage_c_garfield_drift_time_grid.png"
    fig.savefig(target, dpi=165)
    plt.close(fig)

    summary = args.out / "stage_c_garfield_drift_time_grid.csv"
    with summary.open("w", newline="", encoding="utf-8") as stream:
        fields = ("iu", "iv", "u_center_mm", "v_center_mm", "N_trials",
                  "N_collected", "mean_drift_time_ns", "median_drift_time_ns",
                  "P_collected", "P_collected_by_100ns")
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for iv in range(30):
            for iu in range(20):
                writer.writerow({
                    "iu": iu, "iv": iv,
                    "u_center_mm": -1.9 + .2*iu,
                    "v_center_mm": -1.9 + .2*iv,
                    "N_trials": int(started[iv, iu]),
                    "N_collected": int(collected[iv, iu]),
                    "mean_drift_time_ns": mean[iv, iu],
                    "median_drift_time_ns": med[iv, iu],
                    "P_collected": efficiency[iv, iu],
                    "P_collected_by_100ns": in100[iv, iu],
                })
    print(f"Drift-time box map from {args.input}:", flush=True)
    print(f"  trials={int(started.sum()):,}, collected={int(collected.sum()):,}", flush=True)
    print(f"  heat map: {target}", flush=True)
    print(f"  box CSV : {summary}", flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", type=Path, required=True)
    ap.add_argument(
        "--out", type=Path,
        default=Path("test_runs/mwpc_stage_c/c_drift_grid_figures"),
    )
    args = ap.parse_args()
    save(args)


if __name__ == "__main__":
    main()
