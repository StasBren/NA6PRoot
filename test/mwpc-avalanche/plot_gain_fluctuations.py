#!/usr/bin/env python3
"""
Plot the single-electron avalanche-gain distribution produced by
mwpc_phase_a_gain.

The plotted gain proxy is:
    G_eff = number of avalanche electrons collected on an anode wire

This is intentionally a Phase-A response diagnostic, not a calibrated detector
gain measurement. Zero-collected events are reported separately because they
mix pre-avalanche loss/attachment with true multiplication statistics.
"""

import argparse
import csv
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def read_gain_csv(path):
    gains = []
    attached = []
    produced = []
    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            gains.append(int(row["collected_on_wire"]))
            attached.append(int(row["attached"]))
            produced.append(int(row["avalanche_electrons"]))
    return np.asarray(gains), np.asarray(attached), np.asarray(produced)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("csv_file")
    parser.add_argument(
        "--output", default="gain_distribution.png",
        help="Output PNG file")
    parser.add_argument(
        "--bins", type=int, default=60,
        help="Histogram bins over the displayed positive-gain range")
    parser.add_argument(
        "--upper-percentile", type=float, default=99.5,
        help="Upper x-axis percentile for positive-gain events")
    args = parser.parse_args()

    gains, attached, produced = read_gain_csv(args.csv_file)
    n = len(gains)
    if n == 0:
        raise RuntimeError("Input CSV contains no events.")

    positive = gains[gains > 0]
    zero_fraction = np.mean(gains == 0)

    mean_all = np.mean(gains)
    sigma_all = np.std(gains, ddof=1) if n > 1 else 0.0
    cv_all = sigma_all / mean_all if mean_all > 0 else float("nan")
    median_all = np.median(gains)

    if len(positive) == 0:
        raise RuntimeError("No positive-gain events were found.")

    mean_pos = np.mean(positive)
    median_pos = np.median(positive)
    q10_pos, q90_pos = np.quantile(positive, [0.10, 0.90])
    xmax = np.percentile(positive, args.upper_percentile)
    xmax = max(xmax, 1.0)

    shown = positive[positive <= xmax]

    fig, ax = plt.subplots(figsize=(9.0, 5.6))
    ax.hist(
        shown,
        bins=args.bins,
        range=(0, xmax),
        histtype="stepfilled",
        alpha=0.72,
        label="Positive collected gain",
    )

    ax.axvline(mean_pos, linestyle="--", linewidth=1.8,
               label=f"Positive mean = {mean_pos:.0f}")
    ax.axvline(median_pos, linestyle=":", linewidth=2.0,
               label=f"Positive median = {median_pos:.0f}")

    ax.set_xlabel("Collected avalanche electrons per initial electron")
    ax.set_ylabel("Events")
    ax.set_title("Single-electron avalanche gain fluctuations")
    ax.grid(alpha=0.18)
    ax.legend()

    note = (
        f"N = {n}\n"
        f"all-event mean = {mean_all:.0f}\n"
        f"all-event σ/mean = {cv_all:.2f}\n"
        f"zero-collected fraction = {100.0 * zero_fraction:.1f}%\n"
        f"positive 10–90% interval = {q10_pos:.0f}–{q90_pos:.0f}\n"
        f"x-axis truncated at {args.upper_percentile:g}th percentile"
    )
    ax.text(
        0.985, 0.965, note,
        transform=ax.transAxes,
        ha="right", va="top",
        fontsize=9.5,
        bbox=dict(boxstyle="round,pad=0.4", alpha=0.9),
    )

    fig.tight_layout()
    fig.savefig(args.output, dpi=220)
    plt.close(fig)

    print(f"Events                    : {n}")
    print(f"Mean gain, all events     : {mean_all:.3f}")
    print(f"Sigma, all events         : {sigma_all:.3f}")
    print(f"Sigma / mean              : {cv_all:.4f}")
    print(f"Median, all events        : {median_all:.3f}")
    print(f"Zero-collected fraction   : {zero_fraction:.5f}")
    print(f"Positive-gain mean        : {mean_pos:.3f}")
    print(f"Positive-gain median      : {median_pos:.3f}")
    print(f"Positive-gain 10% quantile: {q10_pos:.3f}")
    print(f"Positive-gain 90% quantile: {q90_pos:.3f}")
    print(f"Wrote                     : {args.output}")


if __name__ == "__main__":
    main()
