#!/usr/bin/env python3
"""
Analyze Stage B3a: one ideal tilted strip family.

The key geometry is

    x_alpha = w cos(alpha) - u sin(alpha)

for a +alpha strip family.  A translated avalanche on another anode wire should
therefore look shifted in raw w, but the responses should collapse when plotted
against x_alpha.

This script intentionally emphasizes that geometric statement before any
resolution or reconstruction study.
"""

import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--input", default="stageB_tilt_scan.csv")
    p.add_argument("--output-prefix", default="stageB_tilt_scan")
    p.add_argument("--target-tan-alpha", type=float, default=0.10)
    p.add_argument(
        "--fit-projected-half-range-mm",
        type=float,
        default=0.50,
        help="Use |x_alpha| below this value for local linear sensitivity fits.",
    )
    return p.parse_args()


def mean_response(df):
    group_cols = [
        "tan_alpha",
        "alpha_deg",
        "wire_index",
        "wire_u_mm",
        "xi_w",
        "w_shift_mm",
        "projected_x_mm",
        "projected_eta",
    ]
    value_cols = [
        "three_strip_capture_fraction",
        "neighbor_fraction_three_strip",
        "central_fraction_three_strip",
        "left_right_asymmetry",
        "A_strip_-2_fC",
        "A_strip_-1_fC",
        "A_strip_0_fC",
        "A_strip_1_fC",
        "A_strip_2_fC",
    ]
    return (
        df.groupby(group_cols, as_index=False)[value_cols]
        .mean()
        .sort_values(["tan_alpha", "wire_index", "w_shift_mm"])
    )


def nearest_value(values, target):
    values = np.asarray(sorted(set(values)), dtype=float)
    return float(values[np.argmin(np.abs(values - target))])


def local_linear_fit(g, xcol, ycol, half_range):
    fit = g[np.abs(g[xcol]) <= half_range].copy()
    if len(fit) < 3:
        return np.nan, np.nan
    slope, intercept = np.polyfit(
        fit[xcol].to_numpy(dtype=float),
        fit[ycol].to_numpy(dtype=float),
        1,
    )
    return float(slope), float(intercept)


def main():
    args = parse_args()
    df = pd.read_csv(args.input)

    required = {
        "tan_alpha",
        "alpha_deg",
        "wire_index",
        "wire_u_mm",
        "w_shift_mm",
        "projected_x_mm",
        "projected_eta",
        "left_right_asymmetry",
        "A_strip_-1_fC",
        "A_strip_0_fC",
        "A_strip_1_fC",
    }
    missing = required - set(df.columns)
    if missing:
        raise KeyError(f"Missing columns: {sorted(missing)}")

    mean = mean_response(df)
    mean.to_csv(args.output_prefix + "_mean_response.csv", index=False)

    tan_target = nearest_value(
        mean["tan_alpha"].unique(), args.target_tan_alpha
    )
    target = mean[np.isclose(mean["tan_alpha"], tan_target)].copy()
    alpha_target = float(target["alpha_deg"].iloc[0])

    # --------------------------------------------------------------
    # 1) Raw response versus w.
    # Tilt makes different anode wires appear horizontally displaced.
    # --------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(9.2, 6.0))
    for wire_index, g in target.groupby("wire_index"):
        g = g.sort_values("w_shift_mm")
        wire_u = float(g["wire_u_mm"].iloc[0])
        ax.plot(
            g["w_shift_mm"],
            g["left_right_asymmetry"],
            marker="o",
            label=rf"wire $m={wire_index}$, $u={wire_u:g}$ mm",
        )

    ax.axhline(0.0, linewidth=1.0)
    ax.set_xlabel(r"along-wire avalanche position $w_0$ [mm]")
    ax.set_ylabel(
        r"left-right asymmetry $R=(A_{+1}-A_{-1})/A_3$"
    )
    ax.set_title(
        rf"Stage B3a: raw response for one tilted family, "
        rf"$\tan\alpha={tan_target:g}$ ($\alpha={alpha_target:.2f}^\circ$)"
    )
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(
        args.output_prefix + "_raw_response_vs_w.png",
        dpi=200,
    )
    plt.close(fig)

    # --------------------------------------------------------------
    # 2) Projection test.
    # If the ideal tilted-strip geometry is implemented correctly, curves
    # from different wire indices should collapse against x_alpha.
    # --------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(9.2, 6.0))
    for wire_index, g in target.groupby("wire_index"):
        g = g.sort_values("projected_x_mm")
        wire_u = float(g["wire_u_mm"].iloc[0])
        ax.plot(
            g["projected_x_mm"],
            g["left_right_asymmetry"],
            marker="o",
            label=rf"wire $m={wire_index}$, $u={wire_u:g}$ mm",
        )

    ax.axhline(0.0, linewidth=1.0)
    ax.axvline(0.0, linewidth=1.0)
    ax.set_xlabel(
        r"projected coordinate "
        r"$x_\alpha=w_0\cos\alpha-u_{\rm wire}\sin\alpha$ [mm]"
    )
    ax.set_ylabel(
        r"left-right asymmetry $R=(A_{+1}-A_{-1})/A_3$"
    )
    ax.set_title(
        rf"Stage B3a: projection collapse at $\tan\alpha={tan_target:g}$"
    )
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(
        args.output_prefix + "_projected_coordinate_collapse.png",
        dpi=200,
    )
    plt.close(fig)

    # --------------------------------------------------------------
    # 3) Where does the response center move?
    #
    # x_alpha = 0 predicts exactly
    #       w_center = u_wire tan(alpha).
    #
    # Fit R(w) locally around projected x = 0 and compare.
    # --------------------------------------------------------------
    center_rows = []

    for (tan_alpha, alpha_deg, wire_index, wire_u), g in mean.groupby(
        ["tan_alpha", "alpha_deg", "wire_index", "wire_u_mm"]
    ):
        local = g[
            np.abs(g["projected_x_mm"])
            <= args.fit_projected_half_range_mm
        ].copy()

        if len(local) >= 3:
            slope_w, intercept_w = np.polyfit(
                local["w_shift_mm"].to_numpy(dtype=float),
                local["left_right_asymmetry"].to_numpy(dtype=float),
                1,
            )
            if abs(slope_w) > 1.e-12:
                w_center = -intercept_w / slope_w
            else:
                w_center = np.nan
        else:
            slope_w = np.nan
            w_center = np.nan

        expected = wire_u * tan_alpha

        center_rows.append(
            {
                "tan_alpha": tan_alpha,
                "alpha_deg": alpha_deg,
                "wire_index": wire_index,
                "wire_u_mm": wire_u,
                "measured_w_center_mm": w_center,
                "expected_w_center_mm": expected,
                "center_residual_mm": w_center - expected
                if np.isfinite(w_center)
                else np.nan,
                "local_dR_dw_per_mm": slope_w,
            }
        )

    centers = pd.DataFrame(center_rows)
    centers.to_csv(
        args.output_prefix + "_wire_center_shift_summary.csv",
        index=False,
    )

    fig, ax = plt.subplots(figsize=(9.2, 6.0))
    for tan_alpha, g in centers.groupby("tan_alpha"):
        g = g.sort_values("wire_u_mm")
        alpha_deg = float(g["alpha_deg"].iloc[0])

        ax.plot(
            g["wire_u_mm"],
            g["measured_w_center_mm"],
            marker="o",
            label=rf"measured, $\tan\alpha={tan_alpha:g}$",
        )
        ax.plot(
            g["wire_u_mm"],
            g["expected_w_center_mm"],
            linestyle="--",
            label=rf"$w_{{\rm center}}=u\tan\alpha$ "
                  rf"($\alpha={alpha_deg:.2f}^\circ$)",
        )

    ax.axhline(0.0, linewidth=1.0)
    ax.axvline(0.0, linewidth=1.0)
    ax.set_xlabel(r"anode-wire position $u_{\rm wire}$ [mm]")
    ax.set_ylabel(r"response-center position $w_{\rm center}$ [mm]")
    ax.set_title(
        r"Stage B3a: tilt maps wire position into an along-wire strip shift"
    )
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8, ncol=2)
    fig.tight_layout()
    fig.savefig(
        args.output_prefix + "_wire_center_shift.png",
        dpi=200,
    )
    plt.close(fig)

    # --------------------------------------------------------------
    # 4) Intrinsic one-family sensitivity in its own projected coordinate.
    # It should stay close to the alpha=0 response for small tilts.
    # --------------------------------------------------------------
    sensitivity_rows = []
    for (tan_alpha, alpha_deg), g in mean.groupby(
        ["tan_alpha", "alpha_deg"]
    ):
        slope_x, intercept_x = local_linear_fit(
            g,
            "projected_x_mm",
            "left_right_asymmetry",
            args.fit_projected_half_range_mm,
        )

        sensitivity_rows.append(
            {
                "tan_alpha": tan_alpha,
                "alpha_deg": alpha_deg,
                "dR_dxalpha_per_mm": slope_x,
                "intercept": intercept_x,
            }
        )

    sensitivity = pd.DataFrame(sensitivity_rows).sort_values("tan_alpha")
    sensitivity.to_csv(
        args.output_prefix + "_projected_sensitivity_summary.csv",
        index=False,
    )

    fig, ax = plt.subplots(figsize=(8.4, 5.6))
    ax.plot(
        sensitivity["tan_alpha"],
        sensitivity["dR_dxalpha_per_mm"],
        marker="o",
    )
    ax.set_xlabel(r"$\tan\alpha$")
    ax.set_ylabel(r"$dR/dx_\alpha$ [mm$^{-1}$]")
    ax.set_title(
        r"Stage B3a: one-family sensitivity in the projected coordinate"
    )
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(
        args.output_prefix + "_projected_sensitivity_vs_tilt.png",
        dpi=200,
    )
    plt.close(fig)

    print("\n=== STAGE B3a: TILT PROJECTION CHECK ===")
    print(
        f"target tan(alpha) for response plots = {tan_target:g} "
        f"(alpha = {alpha_target:.3f} deg)"
    )
    print(
        "\nGeometric prediction: x_alpha = w cos(alpha) - u sin(alpha)"
    )
    print(
        "Zero-asymmetry center should follow w_center = u_wire tan(alpha)."
    )

    print("\nMeasured response-center shifts:")
    print(
        centers[
            [
                "tan_alpha",
                "wire_index",
                "wire_u_mm",
                "measured_w_center_mm",
                "expected_w_center_mm",
                "center_residual_mm",
            ]
        ].to_string(index=False)
    )

    finite_res = centers["center_residual_mm"].to_numpy(dtype=float)
    finite_res = finite_res[np.isfinite(finite_res)]
    if len(finite_res):
        print(
            "\nRMS(measured - expected center) = "
            f"{np.sqrt(np.mean(finite_res**2)):.6f} mm"
        )

    print("\nProjected-coordinate sensitivity:")
    print(sensitivity.to_string(index=False))

    print("\nWrote:")
    print(" ", args.output_prefix + "_mean_response.csv")
    print(" ", args.output_prefix + "_raw_response_vs_w.png")
    print(" ", args.output_prefix + "_projected_coordinate_collapse.png")
    print(" ", args.output_prefix + "_wire_center_shift_summary.csv")
    print(" ", args.output_prefix + "_wire_center_shift.png")
    print(" ", args.output_prefix + "_projected_sensitivity_summary.csv")
    print(" ", args.output_prefix + "_projected_sensitivity_vs_tilt.png")


if __name__ == "__main__":
    main()
