#!/usr/bin/env python3
"""
Analyze the Stage-B one-seed ensemble reconstruction study.

At each true position w0 and observation window T, the C++ simulation produces
many independent microscopic avalanches and a three-strip CoG estimate w_hat.

We deliberately separate three quantities:

  bias(w0, T) = <w_hat - w0>

  sigma(w0, T) = std(w_hat - w0)
               = event-by-event stochastic spread at fixed true position

  RMSE(w0, T) = sqrt(<(w_hat - w0)^2>)
              = sqrt(bias^2 + sigma_population^2)

This is the first intrinsic one-seed spatial-response estimate. It is NOT yet
the full chamber resolution because Heed primary-ionisation statistics,
electronics noise/thresholds, real strip geometry, and the magnetic field are
not yet included.
"""

import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument(
        "--input",
        default="stageB_resolution_ensemble.csv",
    )
    p.add_argument(
        "--representative-window-ns",
        type=float,
        default=100.0,
    )
    p.add_argument(
        "--output-prefix",
        default="stageB_resolution_ensemble",
    )
    return p.parse_args()


def main():
    args = parse_args()
    df = pd.read_csv(args.input)

    required = {
        "true_w_mm",
        "window_ns",
        "reconstructable",
        "w_hat_mm",
        "residual_mm",
        "random_seed",
        "garfield_ions",
    }
    missing = required - set(df.columns)
    if missing:
        raise KeyError(f"Missing columns: {sorted(missing)}")

    all_positions = np.sort(df["true_w_mm"].unique())
    all_windows = np.sort(df["window_ns"].unique())

    rows = []
    for (w0, T), g in df.groupby(["true_w_mm", "window_ns"]):
        total = len(g)
        good = g[g["reconstructable"] == 1].copy()
        n = len(good)

        efficiency = n / total if total else np.nan

        if n > 0:
            r = good["residual_mm"].to_numpy(dtype=float)
            what = good["w_hat_mm"].to_numpy(dtype=float)

            bias = float(np.mean(r))
            sigma = float(np.std(r, ddof=1)) if n > 1 else np.nan
            sigma_pop = float(np.std(r, ddof=0))
            rmse = float(np.sqrt(np.mean(r * r)))
            mean_what = float(np.mean(what))
            mean_gain = float(np.mean(good["garfield_ions"]))
            median_gain = float(np.median(good["garfield_ions"]))
        else:
            bias = sigma = sigma_pop = rmse = mean_what = np.nan
            mean_gain = median_gain = np.nan

        rows.append(
            {
                "true_w_mm": w0,
                "window_ns": T,
                "n_total": total,
                "n_reconstructed": n,
                "reconstructable_fraction": efficiency,
                "mean_w_hat_mm": mean_what,
                "bias_mm": bias,
                "sigma_mm": sigma,
                "rmse_mm": rmse,
                "mean_garfield_ions_reconstructed": mean_gain,
                "median_garfield_ions_reconstructed": median_gain,
                "rmse_identity_check_mm": (
                    np.sqrt(bias * bias + sigma_pop * sigma_pop)
                    if np.isfinite(bias) and np.isfinite(sigma_pop)
                    else np.nan
                ),
            }
        )

    summary = pd.DataFrame(rows).sort_values(["window_ns", "true_w_mm"])
    summary.to_csv(args.output_prefix + "_summary.csv", index=False)

    Tref = min(
        all_windows,
        key=lambda x: abs(x - args.representative_window_ns),
    )
    ref = summary[np.isclose(summary["window_ns"], Tref)].sort_values(
        "true_w_mm"
    )

    # ------------------------------------------------------------------
    # 1) Most intuitive plot: ensemble mean reconstruction and its spread.
    # ------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(8.0, 6.5))

    xmin = float(np.min(all_positions))
    xmax = float(np.max(all_positions))
    ax.plot(
        [xmin, xmax],
        [xmin, xmax],
        linestyle="--",
        label="ideal: mean w_hat = w0",
    )

    ax.errorbar(
        ref["true_w_mm"],
        ref["mean_w_hat_mm"],
        yerr=ref["sigma_mm"],
        marker="o",
        capsize=4,
        label=f"ensemble mean +/- sigma, T={Tref:g} ns",
    )

    ax.set_xlabel("true position w0 [mm]")
    ax.set_ylabel("reconstructed position w_hat [mm]")
    ax.set_title(
        "Stage B: ensemble position reconstruction across one strip pitch"
    )
    ax.legend()
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(
        args.output_prefix + "_reconstruction_with_spread.png",
        dpi=200,
    )
    plt.close(fig)

    # ------------------------------------------------------------------
    # 2) Separate systematic bias from stochastic resolution.
    # ------------------------------------------------------------------
    fig, axes = plt.subplots(2, 1, figsize=(9.5, 8.0), sharex=True)
    ax0, ax1 = axes

    for T in all_windows:
        d = summary[np.isclose(summary["window_ns"], T)].sort_values(
            "true_w_mm"
        )
        ax0.plot(
            d["true_w_mm"],
            1000.0 * d["bias_mm"],
            marker="o",
            label=f"{T:g} ns",
        )
        ax1.plot(
            d["true_w_mm"],
            1000.0 * d["sigma_mm"],
            marker="o",
            label=f"{T:g} ns",
        )

    ax0.axhline(0.0, linewidth=0.8)
    ax0.set_ylabel("bias <w_hat - w0> [um]")
    ax0.set_title("Systematic CoG bias")
    ax0.legend()

    ax1.set_xlabel("true position w0 [mm]")
    ax1.set_ylabel("stochastic sigma [um]")
    ax1.set_title("Event-by-event spatial spread at fixed true position")

    fig.suptitle(
        "Stage B: bias and intrinsic one-seed spatial resolution"
    )
    fig.tight_layout()
    fig.savefig(
        args.output_prefix + "_bias_and_resolution_vs_position.png",
        dpi=200,
    )
    plt.close(fig)

    # ------------------------------------------------------------------
    # 3) Reconstruction fraction. For this one-seed study it mostly reflects
    #    zero/small-gain avalanche probability, not full chamber efficiency.
    # ------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(9.5, 5.5))
    for T in all_windows:
        d = summary[np.isclose(summary["window_ns"], T)].sort_values(
            "true_w_mm"
        )
        ax.plot(
            d["true_w_mm"],
            100.0 * d["reconstructable_fraction"],
            marker="o",
            label=f"{T:g} ns",
        )

    ax.set_xlabel("true position w0 [mm]")
    ax.set_ylabel("reconstructable events [%]")
    ax.set_ylim(0.0, 105.0)
    ax.set_title(
        "Stage B: reconstructable fraction in the one-seed avalanche study"
    )
    ax.legend()
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(
        args.output_prefix + "_reconstructable_fraction.png",
        dpi=200,
    )
    plt.close(fig)

    # ------------------------------------------------------------------
    # Compact global numbers by window.
    # ------------------------------------------------------------------
    global_rows = []
    for T in all_windows:
        good = df[
            np.isclose(df["window_ns"], T)
            & (df["reconstructable"] == 1)
        ].copy()

        if len(good):
            r = good["residual_mm"].to_numpy(dtype=float)
            global_rmse = float(np.sqrt(np.mean(r * r)))
        else:
            global_rmse = np.nan

        s = summary[np.isclose(summary["window_ns"], T)]
        global_rows.append(
            {
                "window_ns": T,
                "mean_local_sigma_um": 1000.0 * float(
                    np.nanmean(s["sigma_mm"])
                ),
                "max_local_sigma_um": 1000.0 * float(
                    np.nanmax(s["sigma_mm"])
                ),
                "mean_abs_bias_um": 1000.0 * float(
                    np.nanmean(np.abs(s["bias_mm"]))
                ),
                "max_abs_bias_um": 1000.0 * float(
                    np.nanmax(np.abs(s["bias_mm"]))
                ),
                "global_uncalibrated_rmse_um": 1000.0 * global_rmse,
                "mean_reconstructable_fraction": float(
                    np.nanmean(s["reconstructable_fraction"])
                ),
            }
        )

    global_summary = pd.DataFrame(global_rows)
    global_summary.to_csv(
        args.output_prefix + "_global_summary.csv",
        index=False,
    )

    print("\n=== STAGE B: ENSEMBLE RESOLUTION SUMMARY ===")
    print(
        "Interpretation: local sigma = stochastic spread; "
        "bias = systematic CoG nonlinearity."
    )
    print(
        "This is a one-seed intrinsic study, not yet the full muon-chamber "
        "resolution.\n"
    )
    print(global_summary.to_string(index=False))

    print(f"\nRepresentative window for error-bar plot: {Tref:g} ns")
    print("\nWrote:")
    print(" ", args.output_prefix + "_summary.csv")
    print(" ", args.output_prefix + "_global_summary.csv")
    print(" ", args.output_prefix + "_reconstruction_with_spread.png")
    print(" ", args.output_prefix + "_bias_and_resolution_vs_position.png")
    print(" ", args.output_prefix + "_reconstructable_fraction.png")


if __name__ == "__main__":
    main()
