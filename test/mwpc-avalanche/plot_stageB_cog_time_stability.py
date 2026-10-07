#!/usr/bin/env python3
"""Presentation plot: stability of the CoG position versus observation time.

Input is the summary produced by analyze_stageB_position_scan.py.  The same
microscopic avalanche is translated across one strip pitch, so changes of
w_hat with T isolate the time dependence of the readout response rather than
avalanche-to-avalanche fluctuations.

For each true position w0 we compare

    Delta w_hat(T) = w_hat(T) - w_hat(T_ref)

and summarize the RMS and maximum absolute shift across the scan positions.
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
        "--input",
        default="stageB_position_scan_summary.csv",
        help="Summary from analyze_stageB_position_scan.py",
    )
    p.add_argument(
        "--reference-window-ns",
        type=float,
        default=200.0,
        help="Observation window used as the position reference.",
    )
    p.add_argument(
        "--output-prefix",
        default="stageB_cog_time_stability",
    )
    return p.parse_args()


def main():
    args = parse_args()
    df = pd.read_csv(args.input)

    required = {"true_w_mm", "window_ns", "w_hat_mm"}
    missing = required - set(df.columns)
    if missing:
        raise KeyError(f"Missing columns in {args.input}: {sorted(missing)}")
    if df.empty:
        raise RuntimeError("Input summary is empty.")

    windows = np.sort(df["window_ns"].unique())
    tref = min(windows, key=lambda x: abs(x - args.reference_window_ns))

    ref = (
        df[np.isclose(df["window_ns"], tref)]
        [["true_w_mm", "w_hat_mm"]]
        .rename(columns={"w_hat_mm": "w_hat_ref_mm"})
    )
    if ref.empty:
        raise RuntimeError("Could not construct reference-window sample.")

    rows = []
    joined_by_window = {}

    for T in windows:
        d = (
            df[np.isclose(df["window_ns"], T)]
            [["true_w_mm", "w_hat_mm"]]
            .merge(ref, on="true_w_mm", how="inner")
            .sort_values("true_w_mm")
        )
        d["delta_w_hat_mm"] = d["w_hat_mm"] - d["w_hat_ref_mm"]
        joined_by_window[float(T)] = d

        delta_um = 1000.0 * d["delta_w_hat_mm"].to_numpy(dtype=float)
        rows.append(
            {
                "window_ns": float(T),
                "reference_window_ns": float(tref),
                "n_positions": len(d),
                "mean_shift_um": float(np.mean(delta_um)),
                "rms_shift_um": float(np.sqrt(np.mean(delta_um**2))),
                "max_abs_shift_um": float(np.max(np.abs(delta_um))),
            }
        )

    metrics = pd.DataFrame(rows)
    metrics.to_csv(args.output_prefix + "_metrics.csv", index=False)

    fig, axes = plt.subplots(
        1, 2, figsize=(12.0, 5.2),
        gridspec_kw={"width_ratios": [1.45, 1.0]},
    )
    ax0, ax1 = axes

    for T in windows:
        if np.isclose(T, tref):
            continue
        d = joined_by_window[float(T)]
        ax0.plot(
            d["true_w_mm"],
            1000.0 * d["delta_w_hat_mm"],
            marker="o",
            label=f"{T:g} ns",
        )

    ax0.axhline(0.0, color="black", linewidth=0.9)
    ax0.set_xlabel("true avalanche position $w_0$ [mm]")
    ax0.set_ylabel(
        rf"$\hat w(T)-\hat w({tref:g}\,\mathrm{{ns}})$ [$\mu$m]"
    )
    ax0.set_title("Position shift across one strip pitch")
    ax0.grid(alpha=0.25)
    ax0.legend(title="observation window")

    ax1.plot(
        metrics["window_ns"],
        metrics["rms_shift_um"],
        marker="o",
        label="RMS shift",
    )
    ax1.plot(
        metrics["window_ns"],
        metrics["max_abs_shift_um"],
        marker="s",
        label="maximum |shift|",
    )
    ax1.set_xscale("log")
    ax1.set_xlabel("observation window $T$ [ns]")
    ax1.set_ylabel(r"change in reconstructed position [$\mu$m]")
    ax1.set_title(f"Relative to $T={tref:g}$ ns")
    ax1.grid(alpha=0.25)
    ax1.legend()

    fig.suptitle(
        "Stage B: CoG position remains stable while the signal composition evolves"
    )
    fig.tight_layout()
    fig.savefig(args.output_prefix + ".png", dpi=220)
    plt.close(fig)

    print("\n=== STAGE B: CoG TIME-STABILITY ===")
    print(f"input                    : {args.input}")
    print(f"reference window         : {tref:g} ns")
    print(metrics.to_string(index=False))
    print("\nWrote:")
    print(" ", args.output_prefix + ".png")
    print(" ", args.output_prefix + "_metrics.csv")


if __name__ == "__main__":
    main()
