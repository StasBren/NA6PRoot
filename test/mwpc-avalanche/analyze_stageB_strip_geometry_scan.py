#!/usr/bin/env python3
"""
Analyze the Stage-B2 alpha=0 strip pitch x strip width scan.

We focus on signal-formation observables rather than a final resolution number:

1) Neighbor sharing at xi = 0:
      F_neigh = (A_-1 + A_+1) / (A_-1 + A_0 + A_+1)

2) Left-right position sensitivity near the strip center:
      S = d/dxi [(A_+1 - A_-1) / (A_-1 + A_0 + A_+1)]

3) Three-strip capture:
      F_3 = sum_{k=-1}^{1} A_k / sum_all_strips A_k

The reference sandbox point p = 1.7 mm, width = 1.7 mm is plotted explicitly.
"""

import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--input", default="stageB_strip_geometry_scan.csv")
    p.add_argument("--output-prefix", default="stageB_strip_geometry_scan")
    p.add_argument("--baseline-pitch-mm", type=float, default=1.7)
    p.add_argument("--baseline-width-mm", type=float, default=1.7)
    p.add_argument("--fit-xi-max", type=float, default=0.30)
    return p.parse_args()


def grouped_mean(df):
    cols = [
        "total_abs_signal_all_fC",
        "three_strip_abs_signal_fC",
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
        df.groupby(
            ["strip_pitch_mm", "strip_width_mm", "fill_factor", "xi", "true_w_mm"],
            as_index=False,
        )[cols]
        .mean()
        .sort_values(["strip_pitch_mm", "strip_width_mm", "xi"])
    )


def geometry_summary(mean_df, fit_xi_max):
    rows = []
    for (pitch, width), g in mean_df.groupby(
        ["strip_pitch_mm", "strip_width_mm"]
    ):
        g = g.sort_values("xi")
        center = g.iloc[np.argmin(np.abs(g["xi"].to_numpy(dtype=float)))]

        fit = g[np.abs(g["xi"]) <= fit_xi_max].copy()
        if len(fit) >= 3:
            slope, intercept = np.polyfit(
                fit["xi"].to_numpy(dtype=float),
                fit["left_right_asymmetry"].to_numpy(dtype=float),
                1,
            )
        else:
            slope = np.nan
            intercept = np.nan

        rows.append(
            {
                "strip_pitch_mm": pitch,
                "strip_width_mm": width,
                "fill_factor": width / pitch,
                "center_neighbor_fraction":
                    float(center["neighbor_fraction_three_strip"]),
                "center_central_fraction":
                    float(center["central_fraction_three_strip"]),
                "mean_three_strip_capture_fraction":
                    float(np.mean(g["three_strip_capture_fraction"])),
                "center_total_abs_signal_fC":
                    float(center["total_abs_signal_all_fC"]),
                "asymmetry_slope_per_xi": slope,
                "asymmetry_slope_per_mm": slope / pitch if np.isfinite(slope) else np.nan,
                "asymmetry_intercept": intercept,
            }
        )

    return pd.DataFrame(rows).sort_values(
        ["strip_pitch_mm", "strip_width_mm"]
    )


def heatmap(summary, value_col, title, cbar_label, output):
    pitches = np.sort(summary["strip_pitch_mm"].unique())
    widths = np.sort(summary["strip_width_mm"].unique())

    grid = np.full((len(widths), len(pitches)), np.nan)
    for _, r in summary.iterrows():
        iy = np.where(np.isclose(widths, r["strip_width_mm"]))[0][0]
        ix = np.where(np.isclose(pitches, r["strip_pitch_mm"]))[0][0]
        grid[iy, ix] = r[value_col]

    fig, ax = plt.subplots(figsize=(8.3, 6.2))
    masked = np.ma.masked_invalid(grid)
    im = ax.imshow(
        masked,
        origin="lower",
        aspect="auto",
        interpolation="nearest",
    )

    ax.set_xticks(np.arange(len(pitches)))
    ax.set_xticklabels([f"{x:g}" for x in pitches])
    ax.set_yticks(np.arange(len(widths)))
    ax.set_yticklabels([f"{x:g}" for x in widths])
    ax.set_xlabel("strip pitch p [mm]")
    ax.set_ylabel("strip width s [mm]")
    ax.set_title(title)

    for iy in range(len(widths)):
        for ix in range(len(pitches)):
            if np.isfinite(grid[iy, ix]):
                ax.text(
                    ix, iy, f"{grid[iy, ix]:.2f}",
                    ha="center", va="center", fontsize=8
                )

    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label(cbar_label)
    fig.tight_layout()
    fig.savefig(output, dpi=200)
    plt.close(fig)


def main():
    args = parse_args()
    df = pd.read_csv(args.input)

    required = {
        "strip_pitch_mm",
        "strip_width_mm",
        "fill_factor",
        "xi",
        "true_w_mm",
        "three_strip_capture_fraction",
        "neighbor_fraction_three_strip",
        "central_fraction_three_strip",
        "left_right_asymmetry",
        "A_strip_-1_fC",
        "A_strip_0_fC",
        "A_strip_1_fC",
    }
    missing = required - set(df.columns)
    if missing:
        raise KeyError(f"Missing columns: {sorted(missing)}")

    mean_df = grouped_mean(df)
    summary = geometry_summary(mean_df, args.fit_xi_max)

    mean_df.to_csv(args.output_prefix + "_mean_response.csv", index=False)
    summary.to_csv(args.output_prefix + "_geometry_summary.csv", index=False)

    h1 = summary.copy()
    h1["center_neighbor_percent"] = 100.0 * h1["center_neighbor_fraction"]
    heatmap(
        h1,
        "center_neighbor_percent",
        "Stage B2: neighbour sharing for a centered avalanche",
        r"$F_{\mathrm{neigh}}=(A_{-1}+A_{+1})/A_3$ [%]",
        args.output_prefix + "_neighbor_sharing_heatmap.png",
    )

    heatmap(
        summary,
        "asymmetry_slope_per_xi",
        r"Stage B2: normalized sub-strip position sensitivity at $\alpha=0$",
        r"$dR/d\xi$,  $R=(A_{+1}-A_{-1})/A_3$",
        args.output_prefix + "_position_sensitivity_heatmap.png",
    )

    heatmap(
        summary,
        "asymmetry_slope_per_mm",
        r"Stage B2: physical position sensitivity at $\alpha=0$",
        r"$dR/dw$ [mm$^{-1}$],  $R=(A_{+1}-A_{-1})/A_3$",
        args.output_prefix + "_position_sensitivity_per_mm_heatmap.png",
    )

    h3 = summary.copy()
    h3["capture_percent"] = 100.0 * h3["mean_three_strip_capture_fraction"]
    heatmap(
        h3,
        "capture_percent",
        r"Stage B2: fraction of early signal contained in strips $-1,0,+1$",
        "mean three-strip capture [%]",
        args.output_prefix + "_three_strip_capture_heatmap.png",
    )

    # --------------------------------------------------------------
    # 4) Compact design trade-off view.
    #
    # x: physical position sensitivity per mm
    # y: neighbour sharing at the strip centre
    # marker size/color: three-strip capture
    #
    # A useful geometry should not maximize one axis blindly; it should retain
    # appreciable neighbour sharing and position sensitivity while keeping most
    # of the early signal inside a compact 3-strip cluster.
    # --------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(9.2, 6.5))

    x = summary["asymmetry_slope_per_mm"].to_numpy(dtype=float)
    y = 100.0 * summary["center_neighbor_fraction"].to_numpy(dtype=float)
    capture = 100.0 * summary["mean_three_strip_capture_fraction"].to_numpy(dtype=float)

    sizes = 35.0 + 18.0 * (capture - np.nanmin(capture))
    sc = ax.scatter(
        x,
        y,
        s=sizes,
        c=capture,
        alpha=0.85,
    )

    for _, r in summary.iterrows():
        ax.annotate(
            rf"$p={r['strip_pitch_mm']:g},\ s={r['strip_width_mm']:g}$",
            (
                r["asymmetry_slope_per_mm"],
                100.0 * r["center_neighbor_fraction"],
            ),
            xytext=(4, 4),
            textcoords="offset points",
            fontsize=7,
        )

    baseline = summary[
        np.isclose(summary["strip_pitch_mm"], args.baseline_pitch_mm)
        & np.isclose(summary["strip_width_mm"], args.baseline_width_mm)
    ]
    if not baseline.empty:
        b = baseline.iloc[0]
        ax.scatter(
            [b["asymmetry_slope_per_mm"]],
            [100.0 * b["center_neighbor_fraction"]],
            s=180,
            facecolors="none",
            edgecolors="black",
            linewidths=1.5,
            label=(
                rf"reference $p={args.baseline_pitch_mm:g}$ mm, "
                rf"$s={args.baseline_width_mm:g}$ mm"
            ),
        )
        ax.legend()

    cbar = fig.colorbar(sc, ax=ax)
    cbar.set_label("mean three-strip capture [%]")

    ax.set_xlabel(
        r"physical left-right sensitivity  $dR/dw$ [mm$^{-1}$], "
        r"$R=(A_{+1}-A_{-1})/A_3$"
    )
    ax.set_ylabel(r"neighbour sharing at $\xi=0$ [%]")
    ax.set_title(r"Stage B2: geometry trade-off at $\alpha=0$")
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(
        args.output_prefix + "_geometry_tradeoff.png",
        dpi=200,
    )
    plt.close(fig)

    # --------------------------------------------------------------
    # 5) Reference response curve, useful as a bridge from previous slides.
    # --------------------------------------------------------------
    distances = (
        (mean_df["strip_pitch_mm"] - args.baseline_pitch_mm) ** 2
        + (mean_df["strip_width_mm"] - args.baseline_width_mm) ** 2
    )
    i = int(np.argmin(distances.to_numpy(dtype=float)))
    bp = float(mean_df.iloc[i]["strip_pitch_mm"])
    bw = float(mean_df.iloc[i]["strip_width_mm"])

    base = mean_df[
        np.isclose(mean_df["strip_pitch_mm"], bp)
        & np.isclose(mean_df["strip_width_mm"], bw)
    ].sort_values("xi")

    a3 = (
        base["A_strip_-1_fC"]
        + base["A_strip_0_fC"]
        + base["A_strip_1_fC"]
    )

    fig, ax = plt.subplots(figsize=(9.0, 5.8))
    ax.plot(
        base["xi"],
        100.0 * base["A_strip_-1_fC"] / a3,
        marker="o",
        label=r"strip $-1$",
    )
    ax.plot(
        base["xi"],
        100.0 * base["A_strip_0_fC"] / a3,
        marker="o",
        label=r"strip $0$",
    )
    ax.plot(
        base["xi"],
        100.0 * base["A_strip_1_fC"] / a3,
        marker="o",
        label=r"strip $+1$",
    )
    ax.set_xlabel(r"normalized true position $\xi=w_0/p$")
    ax.set_ylabel("fraction of three-strip early signal [%]")
    ax.set_title(
        rf"Stage B2: reference spatial response, $p={bp:g}$ mm, $s={bw:g}$ mm"
    )
    ax.legend()
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(
        args.output_prefix + "_baseline_response_curve.png",
        dpi=200,
    )
    plt.close(fig)

    print("\n=== STAGE B2: STRIP GEOMETRY SUMMARY ===")
    print(
        "Observables are ensemble-averaged positive-avalanche early signals; "
        "no reconstruction algorithm is assumed."
    )
    print("\ncenter_neighbor_fraction: adjacent-strip sharing at xi=0")
    print("asymmetry_slope_per_xi: sensitivity per fraction of one strip pitch")
    print("asymmetry_slope_per_mm: physical sensitivity per millimetre")
    print("mean_three_strip_capture_fraction: signal retained in a 3-strip cluster")

    display = summary.copy()
    display["center_neighbor_fraction"] *= 100.0
    display["mean_three_strip_capture_fraction"] *= 100.0

    print(
        "\n"
        + display[
            [
                "strip_pitch_mm",
                "strip_width_mm",
                "fill_factor",
                "center_neighbor_fraction",
                "asymmetry_slope_per_xi",
                "asymmetry_slope_per_mm",
                "mean_three_strip_capture_fraction",
            ]
        ].to_string(index=False)
    )

    print("\nWrote:")
    print(" ", args.output_prefix + "_mean_response.csv")
    print(" ", args.output_prefix + "_geometry_summary.csv")
    print(" ", args.output_prefix + "_neighbor_sharing_heatmap.png")
    print(" ", args.output_prefix + "_position_sensitivity_heatmap.png")
    print(" ", args.output_prefix + "_position_sensitivity_per_mm_heatmap.png")
    print(" ", args.output_prefix + "_three_strip_capture_heatmap.png")
    print(" ", args.output_prefix + "_geometry_tradeoff.png")
    print(" ", args.output_prefix + "_baseline_response_curve.png")


if __name__ == "__main__":
    main()
