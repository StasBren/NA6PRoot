#!/usr/bin/env python3
"""
Diagnostic plots for one Stage-B synthetic avalanche-ion cloud run.

Inputs are produced by mwpc_phase_b_ion_cloud:
  <prefix>_waveforms.csv
  <prefix>_summary.csv
  <prefix>_ions.csv

The plots deliberately stay at detector-signal level:
  * strip induced currents I_k(t),
  * full integrated induced charge Q_k = integral I_k dt,
  * spatial charge-sharing profile across strips,
  * start/end w distribution of the synthetic ion cloud.

No frontend-electronics shaping is applied here.
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--prefix", default="stageB_cloud_2mm")
    p.add_argument(
        "--current-strips",
        default="-2,-1,0,1,2",
        help="Comma-separated strip indices to draw (default: -2,-1,0,1,2)",
    )
    p.add_argument(
        "--output-prefix",
        default=None,
        help="Default: <prefix>_diagnostics",
    )
    return p.parse_args()


def label_for_strip(i):
    if i < 0:
        return f"strip_m{-i}"
    if i > 0:
        return f"strip_p{i}"
    return "strip_0"


def weighted_rms(x, weights):
    weights = np.asarray(weights, dtype=float)
    x = np.asarray(x, dtype=float)
    total = np.sum(weights)
    if total <= 0:
        return np.nan, np.nan
    w = weights / total
    mean = np.sum(w * x)
    rms = np.sqrt(np.sum(w * (x - mean) ** 2))
    return mean, rms


def main():
    args = parse_args()
    prefix = Path(args.prefix)
    out = (
        Path(args.output_prefix)
        if args.output_prefix
        else Path(str(prefix) + "_diagnostics")
    )

    wave = pd.read_csv(str(prefix) + "_waveforms.csv")
    summary = pd.read_csv(str(prefix) + "_summary.csv").sort_values(
        "center_w_mm"
    )
    ions = pd.read_csv(str(prefix) + "_ions.csv")

    strip_ids = [
        int(x.strip())
        for x in args.current_strips.split(",")
        if x.strip()
    ]

    t_us = wave["time_ns"].to_numpy(dtype=float) / 1000.0

    # 1) Currents on neighbouring strips.
    plt.figure(figsize=(9, 5.5))
    for i in strip_ids:
        col = label_for_strip(i) + "_fC_per_ns"
        if col not in wave.columns:
            continue
        plt.plot(t_us, wave[col], label=f"strip {i}")
    plt.axhline(0.0, linewidth=0.8)
    plt.xlabel("time [us]")
    plt.ylabel("induced ion current [fC/ns]")
    plt.title("Stage B: synthetic ion cloud - strip currents")
    plt.legend()
    plt.tight_layout()
    plt.savefig(str(out) + "_currents.png", dpi=200)
    plt.close()

    # 2) Integrated charge sharing.
    x = summary["center_w_mm"].to_numpy(dtype=float)
    q = summary["integrated_ion_signal_fC"].to_numpy(dtype=float)
    absq = np.abs(q)
    qnorm = absq / absq.sum() if absq.sum() > 0 else absq
    mean_q, rms_q = weighted_rms(x, absq)

    fig, ax = plt.subplots(figsize=(9, 5.5))
    ax.bar(x, q)
    ax.axhline(0.0, linewidth=0.8)
    ax.set_xlabel("strip center w [mm]")
    ax.set_ylabel("integrated induced charge [fC]")
    ax.set_title(
        f"Stage B: full integrated strip charge "
        f"(abs-weight RMS = {rms_q:.3f} mm)"
    )
    fig.tight_layout()
    fig.savefig(str(out) + "_integrated_charge.png", dpi=200)
    plt.close(fig)

    # 3) Positive sharing diagnostic based on |Q_k|.
    fig, ax = plt.subplots(figsize=(9, 5.5))
    ax.bar(x, 100.0 * qnorm)
    ax.set_xlabel("strip center w [mm]")
    ax.set_ylabel("|Q_k| / sum |Q_j| [%]")
    ax.set_title(
        f"Stage B: charge-sharing diagnostic "
        f"(centroid = {mean_q:.3f} mm, RMS = {rms_q:.3f} mm)"
    )
    fig.tight_layout()
    fig.savefig(str(out) + "_charge_sharing.png", dpi=200)
    plt.close(fig)

    # 4) Synthetic cloud geometry in w.
    w0 = ions["w0_mm"].to_numpy(dtype=float)
    w1 = ions["w1_mm"].to_numpy(dtype=float)
    bins = max(20, min(80, int(np.sqrt(len(w0)) * 3)))

    plt.figure(figsize=(9, 5.5))
    plt.hist(w0, bins=bins, histtype="step", density=True, label="start w")
    plt.hist(w1, bins=bins, histtype="step", density=True, label="end w")
    plt.xlabel("w [mm]")
    plt.ylabel("density [a.u.]")
    plt.title(
        "Stage B: synthetic cloud width along the wire\n"
        f"sigma_start = {np.std(w0, ddof=1):.4f} mm, "
        f"sigma_end = {np.std(w1, ddof=1):.4f} mm"
    )
    plt.legend()
    plt.tight_layout()
    plt.savefig(str(out) + "_cloud_w.png", dpi=200)
    plt.close()

    print("\n=== STAGE B: SYNTHETIC ION-CLOUD DIAGNOSTICS ===")
    print(f"prefix                       : {prefix}")
    print(f"ions                         : {len(ions)}")
    print(f"sigma_w start [mm]           : {np.std(w0, ddof=1):.6f}")
    print(f"sigma_w end [mm]             : {np.std(w1, ddof=1):.6f}")
    print(f"|Q|-weighted centroid [mm]   : {mean_q:.6f}")
    print(f"|Q|-weighted cluster RMS [mm]: {rms_q:.6f}")
    print("\nWrote:")
    for suffix in [
        "_currents.png",
        "_integrated_charge.png",
        "_charge_sharing.png",
        "_cloud_w.png",
    ]:
        print(" ", str(out) + suffix)


if __name__ == "__main__":
    main()
