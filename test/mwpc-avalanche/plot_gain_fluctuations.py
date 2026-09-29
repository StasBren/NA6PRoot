#!/usr/bin/env python3
"""
Plot the single-electron avalanche-gain distribution produced by
mwpc_phase_a_gain.

The plotted gain proxy is:
    G_eff = number of avalanche electrons collected on an anode wire

The histogram is normalised to the total number of positive-gain events, so
each bin shows the fraction of positive-gain trials in that gain interval.
Zero-collected events are reported separately because they mix pre-avalanche
loss/attachment with multiplication statistics.

This is intentionally a Phase-A response diagnostic, not a calibrated detector
gain measurement.
"""

import argparse
import csv

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
    xmax = max(np.percentile(positive, args.upper_percentile), 1.0)

    # Keep the high-gain tail out of the visible range for readability, but
    # normalise with respect to ALL positive-gain events. Therefore the sum of
    # the visible bin fractions is slightly below one when the upper tail is
    # truncated.
    shown = positive[positive <= xmax]
    weights = np.full(
        shown.shape, 1.0 / float(len(positive)), dtype=float
    )

    fig, ax = plt.subplots(figsize=(11.2, 5.8))

    ax.hist(
        shown,
        bins=args.bins,
        range=(0, xmax),
        weights=weights,
        histtype="stepfilled",
        alpha=0.72,
        label="Positive collected gain",
    )

    ax.axvline(
        mean_pos,
        linestyle="--",
        linewidth=1.8,
        label=f"Positive mean = {mean_pos:.0f}",
    )
    ax.axvline(
        median_pos,
        linestyle=":",
        linewidth=2.0,
        label=f"Positive median = {median_pos:.0f}",
    )

    ax.set_xlabel("Collected avalanche electrons per initial electron")
    ax.set_ylabel("Fraction of positive-gain events per bin")
    ax.set_title("Single-electron avalanche gain fluctuations")
    ax.grid(alpha=0.18)

    # Reserve a clean right-hand margin for annotations so that neither the
    # legend nor the numerical summary covers the histogram.
    fig.subplots_adjust(right=0.72)

    ax.legend(
        loc="upper left",
        bbox_to_anchor=(1.02, 1.00),
        borderaxespad=0.0,
        frameon=False,
    )

    note = (
        f"N = {n}\n"
        f"all-event mean = {mean_all:.0f}\n"
        f"all-event σ/mean = {cv_all:.2f}\n"
        f"zero-collected fraction = {100.0 * zero_fraction:.1f}%\n\n"
        f"positive-gain events:\n"
        f"mean = {mean_pos:.0f}\n"
        f"median = {median_pos:.0f}\n"
        f"10–90% interval = {q10_pos:.0f}–{q90_pos:.0f}\n\n"
        f"x-axis truncated at\n"
        f"{args.upper_percentile:g}th percentile"
    )
    ax.text(
        1.02,
        0.68,
        note,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=9.5,
        bbox=dict(boxstyle="round,pad=0.45", facecolor="white", alpha=0.95),
    )

    fig.savefig(args.output, dpi=220, bbox_inches="tight")
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
