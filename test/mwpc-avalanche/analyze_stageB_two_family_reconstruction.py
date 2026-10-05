#!/usr/bin/env python3
"""
Stage B3b: ideal two-family (+alpha / -alpha) position reconstruction.

This script consumes the event-by-event CSV produced by phaseB_tilt_scan.cxx
when that scan is run with BOTH +tan(alpha) and -tan(alpha), for example

    --tan-alphas=-0.10,0.10

The physical avalanche bank is the same for both strip families.  In the ideal
geometry the two families measure the projected coordinates

    x_plus  = w cos(alpha) - u sin(alpha)
    x_minus = w cos(alpha) + u sin(alpha).

Each family supplies a left-right asymmetry

    R = (A_{+1} - A_{-1}) / (A_{-1} + A_0 + A_{+1}),

which is calibrated to the corresponding projected coordinate x.  The two
projected coordinates can then be combined:

    w = (x_plus + x_minus) / (2 cos(alpha))
    u = (x_minus - x_plus) / (2 sin(alpha)).

To avoid evaluating the reconstruction on exactly the same events used to
derive the calibration, events are split by event_index parity:

    even event_index -> calibration
    odd  event_index -> validation.

The calibration is intentionally simple and transparent: for each strip
family we fit x(R) with a cubic polynomial.  This is a proof-of-principle
Stage-B reconstruction, not yet the final detector calibration.

Outputs:
    <prefix>_summary.png
    <prefix>_validation_events.csv
    <prefix>_truth_summary.csv
    <prefix>_metrics.csv
"""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(
        description="Reconstruct (u,w) from two opposite tilted strip families."
    )
    p.add_argument(
        "--input",
        default="stageB_two_family_scan.csv",
        help="Event-level Stage-B tilt scan containing +/- tan(alpha).",
    )
    p.add_argument(
        "--tan-alpha",
        type=float,
        default=0.10,
        help="Positive |tan(alpha)| to use (default 0.10).",
    )
    p.add_argument(
        "--poly-degree",
        type=int,
        default=3,
        help="Degree of x(R) calibration polynomial (default 3).",
    )
    p.add_argument(
        "--output-prefix",
        default="stageB_two_family_reconstruction",
        help="Prefix for plots and CSV outputs.",
    )
    return p.parse_args()


def nearest(values, target):
    values = np.asarray(sorted(set(values)), dtype=float)
    if len(values) == 0:
        raise RuntimeError("No candidate values are available.")
    return float(values[np.argmin(np.abs(values - target))])


def polyfit_x_from_r(df, degree):
    r = df["left_right_asymmetry"].to_numpy(dtype=float)
    x = df["projected_x_mm"].to_numpy(dtype=float)

    if len(r) <= degree:
        raise RuntimeError(
            f"Not enough calibration rows ({len(r)}) for polynomial degree {degree}."
        )

    coeff = np.polyfit(r, x, degree)
    x_fit = np.polyval(coeff, r)
    residual = x_fit - x

    return {
        "coeff": coeff,
        "r_min": float(np.min(r)),
        "r_max": float(np.max(r)),
        "x_min": float(np.min(x)),
        "x_max": float(np.max(x)),
        "bias_mm": float(np.mean(residual)),
        "sigma_mm": float(np.std(residual, ddof=1)),
        "rmse_mm": float(np.sqrt(np.mean(residual ** 2))),
    }


def eval_calibration(r, cal):
    r = np.asarray(r, dtype=float)
    return np.polyval(cal["coeff"], r)


def family_mean_curve(df):
    return (
        df.groupby(
            ["wire_index", "xi_w", "wire_u_mm", "w_shift_mm", "projected_x_mm"],
            as_index=False,
        )["left_right_asymmetry"]
        .mean()
        .sort_values("projected_x_mm")
        .reset_index(drop=True)
    )


def coordinate_metrics(residual):
    residual = np.asarray(residual, dtype=float)
    return {
        "bias_mm": float(np.mean(residual)),
        "sigma_mm": float(np.std(residual, ddof=1)),
        "rmse_mm": float(np.sqrt(np.mean(residual ** 2))),
    }


def main():
    args = parse_args()

    if args.tan_alpha <= 0:
        raise ValueError("--tan-alpha must be positive.")
    if args.poly_degree < 1 or args.poly_degree > 7:
        raise ValueError("--poly-degree must be between 1 and 7.")

    df = pd.read_csv(args.input)

    required = {
        "tan_alpha",
        "alpha_deg",
        "wire_index",
        "wire_u_mm",
        "xi_w",
        "w_shift_mm",
        "projected_x_mm",
        "event_index",
        "random_seed",
        "observation_ns",
        "left_right_asymmetry",
        "A_strip_-1_fC",
        "A_strip_0_fC",
        "A_strip_1_fC",
    }
    missing = required - set(df.columns)
    if missing:
        raise KeyError(f"Missing required columns: {sorted(missing)}")

    positive_values = [x for x in df["tan_alpha"].unique() if x > 0]
    negative_values = [x for x in df["tan_alpha"].unique() if x < 0]

    if not positive_values or not negative_values:
        raise RuntimeError(
            "The input must contain both +tan(alpha) and -tan(alpha). "
            "Rerun mwpc_phase_b_tilt_scan with e.g. "
            "--tan-alphas=-0.10,0.10."
        )

    tan_plus = nearest(positive_values, +args.tan_alpha)
    tan_minus = nearest(negative_values, -args.tan_alpha)

    if abs(abs(tan_plus) - abs(tan_minus)) > 1.e-6:
        raise RuntimeError(
            f"Selected tilts are not opposite: {tan_plus} and {tan_minus}."
        )

    alpha = np.arctan(abs(tan_plus))
    cos_a = np.cos(alpha)
    sin_a = np.sin(alpha)
    alpha_deg = np.degrees(alpha)

    plus = df[np.isclose(df["tan_alpha"], tan_plus)].copy()
    minus = df[np.isclose(df["tan_alpha"], tan_minus)].copy()

    # Same physical avalanche event and same translated true position must be
    # paired across the two strip families.
    pair_keys = [
        "wire_index",
        "wire_u_mm",
        "xi_w",
        "w_shift_mm",
        "event_index",
        "random_seed",
        "observation_ns",
    ]

    plus_cols = pair_keys + [
        "projected_x_mm",
        "left_right_asymmetry",
        "A_strip_-1_fC",
        "A_strip_0_fC",
        "A_strip_1_fC",
    ]
    minus_cols = plus_cols

    paired = plus[plus_cols].merge(
        minus[minus_cols],
        on=pair_keys,
        how="inner",
        suffixes=("_plus", "_minus"),
        validate="one_to_one",
    )

    if paired.empty:
        raise RuntimeError("No +alpha/-alpha event pairs were found.")

    # Split by the avalanche event index, not by geometry point. This ensures
    # that one microscopic avalanche is never in both calibration and
    # validation samples.
    calibration = paired[paired["event_index"] % 2 == 0].copy()
    validation = paired[paired["event_index"] % 2 == 1].copy()

    if calibration.empty or validation.empty:
        raise RuntimeError(
            "Calibration/validation parity split produced an empty sample."
        )

    # Build family-specific calibration datasets.
    cal_plus_df = pd.DataFrame({
        "left_right_asymmetry": calibration["left_right_asymmetry_plus"],
        "projected_x_mm": calibration["projected_x_mm_plus"],
    })
    cal_minus_df = pd.DataFrame({
        "left_right_asymmetry": calibration["left_right_asymmetry_minus"],
        "projected_x_mm": calibration["projected_x_mm_minus"],
    })

    cal_plus = polyfit_x_from_r(cal_plus_df, args.poly_degree)
    cal_minus = polyfit_x_from_r(cal_minus_df, args.poly_degree)

    validation["x_plus_hat_mm"] = eval_calibration(
        validation["left_right_asymmetry_plus"], cal_plus
    )
    validation["x_minus_hat_mm"] = eval_calibration(
        validation["left_right_asymmetry_minus"], cal_minus
    )

    # Flag events that would require extrapolation beyond the calibration
    # asymmetry range. Keep them in the output CSV, but do not use them in the
    # quoted intrinsic-resolution metrics.
    validation["plus_in_calibration_range"] = (
        (validation["left_right_asymmetry_plus"] >= cal_plus["r_min"])
        & (validation["left_right_asymmetry_plus"] <= cal_plus["r_max"])
    )
    validation["minus_in_calibration_range"] = (
        (validation["left_right_asymmetry_minus"] >= cal_minus["r_min"])
        & (validation["left_right_asymmetry_minus"] <= cal_minus["r_max"])
    )
    validation["reconstructable"] = (
        validation["plus_in_calibration_range"]
        & validation["minus_in_calibration_range"]
    )

    validation["w_hat_mm"] = (
        validation["x_plus_hat_mm"] + validation["x_minus_hat_mm"]
    ) / (2.0 * cos_a)
    validation["u_hat_mm"] = (
        validation["x_minus_hat_mm"] - validation["x_plus_hat_mm"]
    ) / (2.0 * sin_a)

    validation["w_residual_mm"] = (
        validation["w_hat_mm"] - validation["w_shift_mm"]
    )
    validation["u_residual_mm"] = (
        validation["u_hat_mm"] - validation["wire_u_mm"]
    )

    good = validation[validation["reconstructable"]].copy()
    if good.empty:
        raise RuntimeError(
            "No validation events lie inside both calibration ranges."
        )

    w_metrics = coordinate_metrics(good["w_residual_mm"])
    u_metrics = coordinate_metrics(good["u_residual_mm"])

    reconstructable_fraction = len(good) / len(validation)

    # Truth-point summary: event-by-event mean and spread at each (u,w).
    truth_summary = (
        good.groupby(
            ["wire_index", "wire_u_mm", "xi_w", "w_shift_mm"],
            as_index=False,
        )
        .agg(
            n=("event_index", "size"),
            u_hat_mean_mm=("u_hat_mm", "mean"),
            u_hat_sigma_mm=("u_hat_mm", "std"),
            w_hat_mean_mm=("w_hat_mm", "mean"),
            w_hat_sigma_mm=("w_hat_mm", "std"),
            x_plus_hat_mean_mm=("x_plus_hat_mm", "mean"),
            x_minus_hat_mean_mm=("x_minus_hat_mm", "mean"),
        )
        .sort_values(["wire_u_mm", "w_shift_mm"])
        .reset_index(drop=True)
    )

    # Aggregate further for clean one-dimensional linearity plots.
    u_summary = (
        good.groupby("wire_u_mm", as_index=False)
        .agg(
            u_hat_mean_mm=("u_hat_mm", "mean"),
            u_hat_sigma_mm=("u_hat_mm", "std"),
            n=("u_hat_mm", "size"),
        )
        .sort_values("wire_u_mm")
    )
    w_summary = (
        good.groupby("w_shift_mm", as_index=False)
        .agg(
            w_hat_mean_mm=("w_hat_mm", "mean"),
            w_hat_sigma_mm=("w_hat_mm", "std"),
            n=("w_hat_mm", "size"),
        )
        .sort_values("w_shift_mm")
    )

    # Family mean curves for a visually clean calibration panel.
    plus_cal_mean = family_mean_curve(
        plus[plus["event_index"] % 2 == 0]
    )
    minus_cal_mean = family_mean_curve(
        minus[minus["event_index"] % 2 == 0]
    )

    # ------------------------------------------------------------------
    # Compact 2x2 summary figure.
    # ------------------------------------------------------------------
    fig, axes = plt.subplots(2, 2, figsize=(13.5, 9.0))
    ax_cal, ax_w, ax_u, ax_res = axes.flat

    fig.suptitle(
        rf"Stage B3b: two-family reconstruction, "
        rf"$\tan\alpha=\pm{abs(tan_plus):g}$ "
        rf"($\alpha={alpha_deg:.2f}^\circ$)",
        fontsize=16,
    )

    # Calibration panel.
    ax_cal.scatter(
        plus_cal_mean["projected_x_mm"],
        plus_cal_mean["left_right_asymmetry"],
        s=28,
        alpha=0.72,
        label=r"$+\alpha$ calibration means",
    )
    ax_cal.scatter(
        minus_cal_mean["projected_x_mm"],
        minus_cal_mean["left_right_asymmetry"],
        s=28,
        alpha=0.72,
        label=r"$-\alpha$ calibration means",
    )

    r_min = min(cal_plus["r_min"], cal_minus["r_min"])
    r_max = max(cal_plus["r_max"], cal_minus["r_max"])
    r_grid = np.linspace(r_min, r_max, 400)

    # Plot x(R) as R(x) by evaluating a dense R grid and swapping axes.
    ax_cal.plot(
        np.polyval(cal_plus["coeff"], r_grid),
        r_grid,
        linewidth=1.8,
        label=r"$+\alpha$: cubic $x(R)$",
    )
    ax_cal.plot(
        np.polyval(cal_minus["coeff"], r_grid),
        r_grid,
        linewidth=1.8,
        linestyle="--",
        label=r"$-\alpha$: cubic $x(R)$",
    )
    ax_cal.axhline(0.0, linewidth=0.8)
    ax_cal.axvline(0.0, linewidth=0.8)
    ax_cal.set_xlabel(r"projected coordinate $x_\pm$ [mm]")
    ax_cal.set_ylabel(r"left-right asymmetry $R_\pm$")
    ax_cal.set_title("One asymmetry calibration per strip family")
    ax_cal.grid(alpha=0.22)
    ax_cal.legend(fontsize=8)

    # w reconstruction.
    ax_w.errorbar(
        w_summary["w_shift_mm"],
        w_summary["w_hat_mean_mm"],
        yerr=w_summary["w_hat_sigma_mm"],
        fmt="o-",
        capsize=3,
        label="validation mean ± event spread",
    )
    lo_w = min(
        w_summary["w_shift_mm"].min(),
        w_summary["w_hat_mean_mm"].min(),
    )
    hi_w = max(
        w_summary["w_shift_mm"].max(),
        w_summary["w_hat_mean_mm"].max(),
    )
    ax_w.plot([lo_w, hi_w], [lo_w, hi_w], "--", label="ideal")
    ax_w.set_xlabel(r"true $w$ [mm]")
    ax_w.set_ylabel(r"reconstructed $\hat w$ [mm]")
    ax_w.set_title(
        rf"Along-wire reconstruction: "
        rf"$\sigma(\hat w-w)={1e3*w_metrics['sigma_mm']:.1f}\ \mu$m"
    )
    ax_w.grid(alpha=0.22)
    ax_w.legend(fontsize=8)

    # u reconstruction.
    ax_u.errorbar(
        u_summary["wire_u_mm"],
        u_summary["u_hat_mean_mm"],
        yerr=u_summary["u_hat_sigma_mm"],
        fmt="o-",
        capsize=3,
        label="validation mean ± event spread",
    )
    lo_u = min(
        u_summary["wire_u_mm"].min(),
        u_summary["u_hat_mean_mm"].min(),
    )
    hi_u = max(
        u_summary["wire_u_mm"].max(),
        u_summary["u_hat_mean_mm"].max(),
    )
    ax_u.plot([lo_u, hi_u], [lo_u, hi_u], "--", label="ideal")
    ax_u.set_xlabel(r"true $u$ [mm]")
    ax_u.set_ylabel(r"reconstructed $\hat u$ [mm]")
    ax_u.set_title(
        rf"Across-wire reconstruction: "
        rf"$\sigma(\hat u-u)={1e3*u_metrics['sigma_mm']:.1f}\ \mu$m"
    )
    ax_u.grid(alpha=0.22)
    ax_u.legend(fontsize=8)

    # Residual distributions.
    ax_res.hist(
        1e3 * good["w_residual_mm"],
        bins=30,
        histtype="step",
        linewidth=1.8,
        label=(
            rf"$\hat w-w$: bias={1e3*w_metrics['bias_mm']:.1f} $\mu$m, "
            rf"$\sigma={1e3*w_metrics['sigma_mm']:.1f}$ $\mu$m"
        ),
    )
    ax_res.hist(
        1e3 * good["u_residual_mm"],
        bins=30,
        histtype="step",
        linewidth=1.8,
        label=(
            rf"$\hat u-u$: bias={1e3*u_metrics['bias_mm']:.1f} $\mu$m, "
            rf"$\sigma={1e3*u_metrics['sigma_mm']:.1f}$ $\mu$m"
        ),
    )
    ax_res.axvline(0.0, linewidth=0.8)
    ax_res.set_xlabel("reconstruction residual [µm]")
    ax_res.set_ylabel("validation events")
    ax_res.set_title(
        rf"Intrinsic two-family residuals; usable fraction "
        rf"$={100*reconstructable_fraction:.1f}\%$"
    )
    ax_res.grid(alpha=0.20)
    ax_res.legend(fontsize=8)

    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.95))

    prefix = Path(args.output_prefix)
    summary_png = Path(str(prefix) + "_summary.png")
    validation_csv = Path(str(prefix) + "_validation_events.csv")
    truth_csv = Path(str(prefix) + "_truth_summary.csv")
    metrics_csv = Path(str(prefix) + "_metrics.csv")

    fig.savefig(summary_png, dpi=170)
    plt.close(fig)

    validation.to_csv(validation_csv, index=False)
    truth_summary.to_csv(truth_csv, index=False)

    metrics = pd.DataFrame(
        [
            {
                "quantity": "w",
                **w_metrics,
                "sigma_um": 1e3 * w_metrics["sigma_mm"],
                "bias_um": 1e3 * w_metrics["bias_mm"],
                "rmse_um": 1e3 * w_metrics["rmse_mm"],
            },
            {
                "quantity": "u",
                **u_metrics,
                "sigma_um": 1e3 * u_metrics["sigma_mm"],
                "bias_um": 1e3 * u_metrics["bias_mm"],
                "rmse_um": 1e3 * u_metrics["rmse_mm"],
            },
        ]
    )
    metrics["tan_alpha_abs"] = abs(tan_plus)
    metrics["alpha_deg_abs"] = alpha_deg
    metrics["validation_events"] = len(validation)
    metrics["reconstructable_events"] = len(good)
    metrics["reconstructable_fraction"] = reconstructable_fraction
    metrics.to_csv(metrics_csv, index=False)

    print("\n=== STAGE B3b: IDEAL TWO-FAMILY RECONSTRUCTION ===")
    print(f"input                         : {args.input}")
    print(f"selected tan(alpha)           : +{tan_plus:g}, {tan_minus:g}")
    print(f"|alpha|                       : {alpha_deg:.4f} deg")
    print(f"paired events                 : {len(paired)}")
    print(f"calibration pairs             : {len(calibration)}")
    print(f"validation pairs              : {len(validation)}")
    print(
        f"inside both calibration ranges: {len(good)} / {len(validation)} "
        f"({100*reconstructable_fraction:.2f}%)"
    )
    print("")
    print("Calibration quality x(R):")
    print(
        f"  +alpha sigma / RMSE         : "
        f"{1e3*cal_plus['sigma_mm']:.2f} / "
        f"{1e3*cal_plus['rmse_mm']:.2f} um"
    )
    print(
        f"  -alpha sigma / RMSE         : "
        f"{1e3*cal_minus['sigma_mm']:.2f} / "
        f"{1e3*cal_minus['rmse_mm']:.2f} um"
    )
    print("")
    print("Two-family validation:")
    print(
        f"  w bias / sigma / RMSE       : "
        f"{1e3*w_metrics['bias_mm']:.2f} / "
        f"{1e3*w_metrics['sigma_mm']:.2f} / "
        f"{1e3*w_metrics['rmse_mm']:.2f} um"
    )
    print(
        f"  u bias / sigma / RMSE       : "
        f"{1e3*u_metrics['bias_mm']:.2f} / "
        f"{1e3*u_metrics['sigma_mm']:.2f} / "
        f"{1e3*u_metrics['rmse_mm']:.2f} um"
    )
    print("")
    print("Interpretation:")
    print(
        "  The sum x_plus + x_minus reconstructs w, while their difference "
        "reconstructs u."
    )
    print(
        "  Because u is divided by 2 sin(alpha), a small tilt angle amplifies "
        "the projected-coordinate fluctuations in the u reconstruction."
    )
    print("")
    print(f"Wrote                         : {summary_png}")
    print(f"Wrote                         : {validation_csv}")
    print(f"Wrote                         : {truth_csv}")
    print(f"Wrote                         : {metrics_csv}")


if __name__ == "__main__":
    main()
