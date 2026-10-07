#!/usr/bin/env python3
"""Compare microscopic avalanche multiplication timing with prompt strip current.

The real-avalanche executable already exports:
  <prefix>_electron_endpoints.csv
  <prefix>_electron_waveforms.csv

Each row of electron_endpoints is one microscopic electron trajectory.
Its t0_ns is the trajectory start time; for secondary trajectories this is
where/when a new avalanche electron is created.  Histogramming t0 therefore
provides a direct diagnostic of the avalanche multiplication chronology.

This plot deliberately does NOT claim that dN/dt must be identical to current:
Shockley-Ramo current also depends on electron velocity and strip weighting
field.  We use the timing comparison only to test whether multiple prompt
current peaks coincide with bursts of microscopic avalanche production.
"""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument(
        "--prefix",
        default="stageB_posscan_p0p40",
        help="Prefix produced by mwpc_phase_b_real_avalanche.",
    )
    p.add_argument(
        "--bin-ns",
        type=float,
        default=0.02,
        help="Histogram bin width for electron trajectory start times.",
    )
    p.add_argument(
        "--threshold-fraction",
        type=float,
        default=0.01,
        help="Fraction of the peak prompt current used to define the zoom window.",
    )
    p.add_argument(
        "--output",
        default="stageB_avalanche_burst_structure.png",
    )
    return p.parse_args()


def strip_label(i):
    if i < 0:
        return f"strip_m{-i}"
    if i > 0:
        return f"strip_p{i}"
    return "strip_0"


def main():
    args = parse_args()
    prefix = Path(args.prefix)

    endpoints_file = str(prefix) + "_electron_endpoints.csv"
    wave_file = str(prefix) + "_electron_waveforms.csv"

    endpoints = pd.read_csv(endpoints_file)
    wave = pd.read_csv(wave_file)

    if endpoints.empty:
        raise RuntimeError(f"{endpoints_file} is empty.")
    if wave.empty:
        raise RuntimeError(f"{wave_file} is empty.")
    if args.bin_ns <= 0:
        raise ValueError("--bin-ns must be positive.")

    t0 = endpoints["t0_ns"].to_numpy(dtype=float)
    t = wave["time_ns"].to_numpy(dtype=float)

    key = [-1, 0, 1]
    currents = {}
    for k in key:
        col = f"{strip_label(k)}_electron_fC_per_ns"
        if col not in wave.columns:
            raise KeyError(f"Missing waveform column {col!r}")
        currents[k] = wave[col].to_numpy(dtype=float)

    # Define the displayed interval from the actual prompt current.
    envelope = np.max(np.vstack([np.abs(currents[k]) for k in key]), axis=0)
    peak_current = float(np.max(envelope))
    if peak_current <= 0:
        raise RuntimeError("No non-zero prompt electron current found.")

    active = np.flatnonzero(
        envelope >= args.threshold_fraction * peak_current
    )
    if not len(active):
        raise RuntimeError("Could not identify the avalanche-current interval.")

    dt_wave = float(np.median(np.diff(t)))
    margin = max(0.08, 5.0 * dt_wave)
    xmin = max(float(t[0]), float(t[active[0]]) - margin)
    xmax = min(float(t[-1]), float(t[active[-1]]) + margin)

    # Histogram trajectory starts on the same scale.  The first trajectory is
    # the seed electron; one count is negligible for a developed avalanche but
    # we keep it for a literal chronology of all Garfield electron trajectories.
    start = np.floor(min(xmin, float(np.min(t0))) / args.bin_ns) * args.bin_ns
    stop = np.ceil(max(xmax, float(np.max(t0))) / args.bin_ns) * args.bin_ns
    edges = np.arange(start, stop + 1.01 * args.bin_ns, args.bin_ns)
    counts, edges = np.histogram(t0, bins=edges)
    centers = 0.5 * (edges[:-1] + edges[1:])

    in_zoom = (centers >= xmin) & (centers <= xmax)
    current_zoom = (t >= xmin) & (t <= xmax)

    # A smoothed creation-rate guide suppresses single-bin discreteness while
    # retaining sub-ns burst structure.
    kernel_bins = max(1, int(round(0.04 / args.bin_ns)))
    kernel = np.ones(kernel_bins, dtype=float) / kernel_bins
    smooth_counts = np.convolve(counts.astype(float), kernel, mode="same")

    fig, axes = plt.subplots(
        2, 1, figsize=(9.6, 7.0), sharex=True,
        gridspec_kw={"height_ratios": [1.0, 1.35]},
    )
    ax_n, ax_i = axes

    ax_n.bar(
        centers[in_zoom],
        counts[in_zoom],
        width=0.92 * args.bin_ns,
        alpha=0.45,
        label=r"electron trajectory starts per bin",
    )
    ax_n.plot(
        centers[in_zoom],
        smooth_counts[in_zoom],
        linewidth=2.0,
        label=r"smoothed creation-rate guide",
    )
    ax_n.set_ylabel("new electron trajectories / bin")
    ax_n.set_title(
        "Microscopic avalanche multiplication occurs in short bursts"
    )
    ax_n.grid(alpha=0.2)
    ax_n.legend()

    styles = {-1: "--", 0: "-", 1: ":"}
    for k in key:
        ax_i.plot(
            t[current_zoom],
            currents[k][current_zoom],
            linestyle=styles[k],
            label=f"strip {k}",
        )
    ax_i.set_xlabel("time [ns]")
    ax_i.set_ylabel("electron current [fC/ns]")
    ax_i.set_title(
        "Prompt Shockley-Ramo current from the same microscopic avalanche"
    )
    ax_i.grid(alpha=0.2)
    ax_i.legend()

    ax_i.set_xlim(xmin, xmax)

    fig.suptitle(
        "Stage B diagnostic: avalanche production timing vs prompt strip current"
    )
    fig.tight_layout()
    fig.savefig(args.output, dpi=220)
    plt.close(fig)

    # Also write a compact aligned table for inspection / future slide notes.
    hist = pd.DataFrame(
        {
            "time_ns": centers,
            "electron_trajectory_starts": counts,
            "smoothed_starts": smooth_counts,
        }
    )
    hist = hist[(hist["time_ns"] >= xmin) & (hist["time_ns"] <= xmax)]
    csv_out = str(Path(args.output).with_suffix("")) + "_histogram.csv"
    hist.to_csv(csv_out, index=False)

    # Print the strongest multiplication bins and central-strip current bins.
    n_top = min(5, len(hist))
    top_birth = hist.nlargest(n_top, "smoothed_starts")

    mask = current_zoom
    tmp = pd.DataFrame(
        {
            "time_ns": t[mask],
            "abs_strip0_current": np.abs(currents[0][mask]),
        }
    )
    top_current = tmp.nlargest(min(5, len(tmp)), "abs_strip0_current")

    print("\n=== STAGE B: AVALANCHE BURST STRUCTURE ===")
    print(f"prefix                    : {args.prefix}")
    print(f"electron trajectories     : {len(endpoints)}")
    print(f"histogram bin width       : {args.bin_ns:g} ns")
    print(f"display interval          : {xmin:.4f} .. {xmax:.4f} ns")
    print("\nStrongest electron-production bins (smoothed):")
    print(top_birth.to_string(index=False))
    print("\nStrongest strip-0 prompt-current bins:")
    print(top_current.to_string(index=False))
    print("\nInterpretation:")
    print(
        "  Coincident timing supports the interpretation that multiple prompt "
        "current peaks reflect bursty development of one microscopic avalanche. "
        "Exact amplitudes need not match because i = -q v dot E_w."
    )
    print("\nWrote:")
    print(" ", args.output)
    print(" ", csv_out)


if __name__ == "__main__":
    main()
