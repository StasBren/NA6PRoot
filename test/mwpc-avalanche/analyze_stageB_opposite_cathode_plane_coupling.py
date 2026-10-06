#!/usr/bin/env python3
"""Analyse Stage B3c wire-aware complete-cathode coupling."""

import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument(
        "--input",
        default="stageB_opposite_cathode_plane_coupling.csv",
    )
    p.add_argument(
        "--output-prefix",
        default="stageB_opposite_cathode_plane_coupling",
    )
    return p.parse_args()


def safe_corr(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if len(x) < 2 or np.std(x) == 0. or np.std(y) == 0.:
        return np.nan
    return float(np.corrcoef(x, y)[0, 1])


def main():
    args = parse_args()
    df = pd.read_csv(args.input)
    if df.empty:
        raise RuntimeError("Input CSV contains no usable events.")

    qminus = df["q_minus_total_fC"].to_numpy(dtype=float)
    qplus = df["q_plus_total_fC"].to_numpy(dtype=float)
    aminus = np.abs(qminus)
    aplus = np.abs(qplus)

    rho_signed = safe_corr(qminus, qplus)
    rho_abs = safe_corr(aminus, aplus)

    eps = 1.e-30
    ratio_plus_over_minus = aplus / np.maximum(aminus, eps)

    summary = pd.DataFrame([{
        "events": len(df),
        "seed_side": str(df["seed_side"].iloc[0]),
        "mean_abs_q_minus_fC": float(np.mean(aminus)),
        "mean_abs_q_plus_fC": float(np.mean(aplus)),
        "median_abs_q_minus_fC": float(np.median(aminus)),
        "median_abs_q_plus_fC": float(np.median(aplus)),
        "mean_plus_over_minus_abs_ratio": float(np.mean(ratio_plus_over_minus)),
        "median_plus_over_minus_abs_ratio": float(np.median(ratio_plus_over_minus)),
        "rho_signed_qminus_qplus": rho_signed,
        "rho_abs_qminus_qplus": rho_abs,
        "mean_minus_fraction_of_two_cathodes":
            float(np.mean(df["abs_minus_fraction_of_two_cathodes"])),
        "mean_plus_fraction_of_two_cathodes":
            float(np.mean(df["abs_plus_fraction_of_two_cathodes"])),
    }])
    summary.to_csv(args.output_prefix + "_summary.csv", index=False)

    # Plot 1: event-by-event physical coupling to the two complete cathodes.
    fig, ax = plt.subplots(figsize=(7.3, 6.2))
    ax.scatter(aminus, aplus, s=14, alpha=0.28)
    ax.set_xlabel("|Q| on minus cathode at Tobs [fC]")
    ax.set_ylabel("|Q| on plus cathode at Tobs [fC]")
    ax.set_title(
        "Stage B3c: same avalanche coupled to both physical cathodes\n"
        f"corr(|Q-|, |Q+|) = {rho_abs:.4f}"
    )
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(args.output_prefix + "_cathode_correlation.png", dpi=200)
    plt.close(fig)

    # Plot 2: mean electron/ion/total magnitudes on each side.
    labels = ["electron", "positive ion", "total"]
    minus_vals = [
        np.mean(np.abs(df["q_minus_electron_fC"])),
        np.mean(np.abs(df["q_minus_ion_fC"])),
        np.mean(np.abs(df["q_minus_total_fC"])),
    ]
    plus_vals = [
        np.mean(np.abs(df["q_plus_electron_fC"])),
        np.mean(np.abs(df["q_plus_ion_fC"])),
        np.mean(np.abs(df["q_plus_total_fC"])),
    ]

    x = np.arange(len(labels))
    width = 0.36
    fig, ax = plt.subplots(figsize=(8.0, 5.6))
    ax.bar(x - width/2, minus_vals, width, label="minus cathode")
    ax.bar(x + width/2, plus_vals, width, label="plus cathode")
    ax.set_xticks(x, labels)
    ax.set_ylabel("mean |induced charge| at Tobs [fC]")
    ax.set_title("Stage B3c: early-signal coupling to the two cathodes")
    ax.grid(axis="y", alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(args.output_prefix + "_component_comparison.png", dpi=200)
    plt.close(fig)

    # Plot 3: event-by-event fraction of the two-cathode signal.
    fig, ax = plt.subplots(figsize=(7.6, 5.4))
    ax.hist(
        100. * df["abs_minus_fraction_of_two_cathodes"],
        bins=30, alpha=0.65, label="minus fraction"
    )
    ax.hist(
        100. * df["abs_plus_fraction_of_two_cathodes"],
        bins=30, alpha=0.65, label="plus fraction"
    )
    ax.set_xlabel("fraction of |Q-| + |Q+| [%]")
    ax.set_ylabel("events")
    ax.set_title("Stage B3c: how the early cathode signal is shared between sides")
    ax.grid(axis="y", alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(args.output_prefix + "_side_fraction.png", dpi=200)
    plt.close(fig)

    print("\n=== STAGE B3c: OPPOSITE-CATHODE COMPLETE-PLANE CONTROL ===")
    print(f"events                         : {len(df)}")
    print(f"seed side                      : {df['seed_side'].iloc[0]}")
    print(
        "mean |Q-| / |Q+| [fC]         : "
        f"{np.mean(aminus):.6g} / {np.mean(aplus):.6g}"
    )
    print(
        "median |Q-| / |Q+| [fC]       : "
        f"{np.median(aminus):.6g} / {np.median(aplus):.6g}"
    )
    print(
        "median |Q+|/|Q-|               : "
        f"{np.median(ratio_plus_over_minus):.6g}"
    )
    print(f"corr signed Q-/Q+              : {rho_signed:.6f}")
    print(f"corr |Q-|/|Q+|                 : {rho_abs:.6f}")
    print(
        "mean two-cathode fractions -/+ : "
        f"{100.*np.mean(df['abs_minus_fraction_of_two_cathodes']):.2f}% / "
        f"{100.*np.mean(df['abs_plus_fraction_of_two_cathodes']):.2f}%"
    )
    print("\nWrote:")
    print(" ", args.output_prefix + "_summary.csv")
    print(" ", args.output_prefix + "_cathode_correlation.png")
    print(" ", args.output_prefix + "_component_comparison.png")
    print(" ", args.output_prefix + "_side_fraction.png")


if __name__ == "__main__":
    main()
