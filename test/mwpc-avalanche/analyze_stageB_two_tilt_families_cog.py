#!/usr/bin/env python3
"""
Stage B3b: proposal-aligned two-family reconstruction using cluster CoG.

This is the clean successor to the first calibrated-R proof of principle.

For each tilted strip family, use the induced strip amplitudes directly:

    x_CoG = sum_k A_k x_k / sum_k A_k

with strip centres x_k = k p in that family's natural projected coordinate.

The two families measure

    x_plus  = w cos(alpha) - u sin(alpha)
    x_minus = w cos(alpha) + u sin(alpha)

and therefore

    w = (x_plus + x_minus) / (2 cos(alpha))
    u = (x_minus - x_plus) / (2 sin(alpha))

No R->x calibration is used here.

The input is the Stage-B tilt scan CSV containing BOTH +alpha and -alpha.
The same microscopic Garfield avalanche bank is reused for both families, so
rows can be paired event-by-event.

Scope: ideal two-family readout.  This script still does not choose whether
the two physical strip families live on the same cathode or on opposite
cathodes; that electrostatic placement study is the next step.
"""

import argparse
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--input", default="stageB_two_tilt_families_large.csv")
    p.add_argument("--output-prefix", default="stageB_two_tilt_families_cog")
    p.add_argument("--tan-alpha", type=float, default=0.10)
    p.add_argument("--strip-pitch-mm", type=float, default=1.7)
    return p.parse_args()


def nearest_value(values, target):
    values = np.asarray(sorted(set(values)), dtype=float)
    return float(values[np.argmin(np.abs(values - target))])


def discover_amplitude_columns(df):
    found = []
    pat = re.compile(r"^A_strip_(-?[0-9]+)_fC$")
    for col in df.columns:
        m = pat.match(col)
        if m:
            found.append((int(m.group(1)), col))
    found.sort(key=lambda x: x[0])
    if not found:
        raise RuntimeError("No A_strip_k_fC columns found in input CSV.")
    strip_ids = np.array([k for k, _ in found], dtype=int)
    amp_cols = [col for _, col in found]
    return strip_ids, amp_cols


def add_cluster_cog(df, strip_ids, amp_cols, pitch_mm):
    out = df.copy()
    amps = out[amp_cols].to_numpy(dtype=float)
    strip_x = strip_ids.astype(float) * pitch_mm

    asum = np.sum(amps, axis=1)
    good = np.isfinite(asum) & (asum > 0.)

    x_cog = np.full(len(out), np.nan, dtype=float)
    x_cog[good] = (
        np.sum(amps[good] * strip_x[None, :], axis=1) / asum[good]
    )

    # Diagnostic compact-cluster capture around strongest strip.
    imax = np.argmax(amps, axis=1)
    capture3 = np.full(len(out), np.nan)
    capture5 = np.full(len(out), np.nan)
    for i, j in enumerate(imax):
        if not good[i]:
            continue
        s3 = np.sum(amps[i, max(0, j - 1):min(len(strip_ids), j + 2)])
        s5 = np.sum(amps[i, max(0, j - 2):min(len(strip_ids), j + 3)])
        capture3[i] = s3 / asum[i]
        capture5[i] = s5 / asum[i]

    out["x_cog_mm"] = x_cog
    out["cluster_sum_fC"] = asum
    out["cluster_reconstructable"] = good & np.isfinite(x_cog)
    out["capture3_fraction"] = capture3
    out["capture5_fraction"] = capture5
    return out


def metrics(residual):
    residual = np.asarray(residual, dtype=float)
    return {
        "bias": float(np.mean(residual)),
        "sigma": float(np.std(residual, ddof=1)) if len(residual) > 1 else np.nan,
        "rms": float(np.sqrt(np.mean(residual**2))),
    }


def main():
    args = parse_args()

    if args.strip_pitch_mm <= 0.:
        raise ValueError("--strip-pitch-mm must be positive")
    if args.tan_alpha <= 0.:
        raise ValueError("--tan-alpha must be positive; both +/- values are used")

    df = pd.read_csv(args.input)

    required = {
        "tan_alpha", "alpha_deg", "wire_index", "wire_u_mm",
        "xi_w", "w_shift_mm", "projected_x_mm",
        "event_index", "random_seed", "observation_ns",
    }
    missing = required - set(df.columns)
    if missing:
        raise KeyError(f"Missing columns: {sorted(missing)}")

    strip_ids, amp_cols = discover_amplitude_columns(df)
    p = args.strip_pitch_mm

    tplus = nearest_value(df["tan_alpha"].unique(), +args.tan_alpha)
    tminus = nearest_value(df["tan_alpha"].unique(), -args.tan_alpha)
    if tplus <= 0. or tminus >= 0.:
        raise RuntimeError("Input must contain both positive and negative tilt values.")

    plus = add_cluster_cog(
        df[np.isclose(df["tan_alpha"], tplus)].copy(),
        strip_ids, amp_cols, p
    )
    minus = add_cluster_cog(
        df[np.isclose(df["tan_alpha"], tminus)].copy(),
        strip_ids, amp_cols, p
    )

    alpha = abs(float(plus["alpha_deg"].iloc[0])) * np.pi / 180.
    cos_a = np.cos(alpha)
    sin_a = np.sin(alpha)

    keys = [
        "wire_index", "wire_u_mm", "xi_w", "w_shift_mm",
        "event_index", "random_seed",
    ]

    pcols = keys + [
        "projected_x_mm", "x_cog_mm", "cluster_reconstructable",
        "capture3_fraction", "capture5_fraction",
    ]
    mcols = keys + [
        "projected_x_mm", "x_cog_mm", "cluster_reconstructable",
        "capture3_fraction", "capture5_fraction",
    ]

    paired_all = plus[pcols].merge(
        minus[mcols],
        on=keys,
        suffixes=("_plus", "_minus"),
        how="inner",
    )
    total_pairs = len(paired_all)

    paired = paired_all[
        paired_all["cluster_reconstructable_plus"]
        & paired_all["cluster_reconstructable_minus"]
    ].copy()
    if paired.empty:
        raise RuntimeError("No reconstructable paired events.")

    # Per-family projected-coordinate residuals.
    paired["dx_plus_mm"] = (
        paired["x_cog_mm_plus"] - paired["projected_x_mm_plus"]
    )
    paired["dx_minus_mm"] = (
        paired["x_cog_mm_minus"] - paired["projected_x_mm_minus"]
    )

    # 2D inversion.
    paired["w_reco_mm"] = (
        paired["x_cog_mm_plus"] + paired["x_cog_mm_minus"]
    ) / (2. * cos_a)
    paired["u_reco_mm"] = (
        paired["x_cog_mm_minus"] - paired["x_cog_mm_plus"]
    ) / (2. * sin_a)

    paired["w_true_mm"] = paired["w_shift_mm"]
    paired["u_true_mm"] = paired["wire_u_mm"]

    paired["dw_mm"] = paired["w_reco_mm"] - paired["w_true_mm"]
    paired["du_mm"] = paired["u_reco_mm"] - paired["u_true_mm"]

    mxp = metrics(paired["dx_plus_mm"])
    mxm = metrics(paired["dx_minus_mm"])
    mw = metrics(paired["dw_mm"])
    mu = metrics(paired["du_mm"])

    reconstructable_fraction = len(paired) / total_pairs if total_pairs else 0.

    paired.to_csv(
        args.output_prefix + "_event_reconstruction.csv",
        index=False,
    )

    summary = pd.DataFrame([{
        "tan_alpha_magnitude": args.tan_alpha,
        "alpha_deg": alpha * 180. / np.pi,
        "strip_pitch_mm": p,
        "first_exported_strip": int(strip_ids[0]),
        "last_exported_strip": int(strip_ids[-1]),
        "n_exported_strips": int(len(strip_ids)),
        "paired_events_before_cut": int(total_pairs),
        "paired_reconstructable_events": int(len(paired)),
        "reconstructable_fraction": reconstructable_fraction,
        "xplus_bias_mm": mxp["bias"],
        "xplus_sigma_mm": mxp["sigma"],
        "xplus_rms_mm": mxp["rms"],
        "xminus_bias_mm": mxm["bias"],
        "xminus_sigma_mm": mxm["sigma"],
        "xminus_rms_mm": mxm["rms"],
        "w_bias_mm": mw["bias"],
        "w_sigma_mm": mw["sigma"],
        "w_rms_mm": mw["rms"],
        "u_bias_mm": mu["bias"],
        "u_sigma_mm": mu["sigma"],
        "u_rms_mm": mu["rms"],
        "mean_capture3_fraction": float(np.nanmean(np.r_[
            paired["capture3_fraction_plus"],
            paired["capture3_fraction_minus"],
        ])),
        "mean_capture5_fraction": float(np.nanmean(np.r_[
            paired["capture5_fraction_plus"],
            paired["capture5_fraction_minus"],
        ])),
        "projection_error_correlation_rho": float(
            np.corrcoef(
                paired["dx_plus_mm"].to_numpy(dtype=float),
                paired["dx_minus_mm"].to_numpy(dtype=float),
            )[0, 1]
        ),
        "common_mode_sigma_um": float(
            np.std(
                0.5e3 * (
                    paired["dx_plus_mm"].to_numpy(dtype=float)
                    + paired["dx_minus_mm"].to_numpy(dtype=float)
                ),
                ddof=1,
            )
        ),
        "differential_mode_sigma_um": float(
            np.std(
                0.5e3 * (
                    paired["dx_minus_mm"].to_numpy(dtype=float)
                    - paired["dx_plus_mm"].to_numpy(dtype=float)
                ),
                ddof=1,
            )
        ),
    }])
    summary.to_csv(args.output_prefix + "_summary.csv", index=False)

    # ------------------------------------------------------------
    # Plot 1: each strip family reconstructs its own projected coordinate.
    # ------------------------------------------------------------
    fig, (axp, axm) = plt.subplots(1, 2, figsize=(11.4, 5.2))

    for ax, truth_col, reco_col, title, m in [
        (
            axp, "projected_x_mm_plus", "x_cog_mm_plus",
            "+alpha family", mxp
        ),
        (
            axm, "projected_x_mm_minus", "x_cog_mm_minus",
            "-alpha family", mxm
        ),
    ]:
        truth = paired[truth_col].to_numpy()
        reco = paired[reco_col].to_numpy()
        lo = min(np.min(truth), np.min(reco))
        hi = max(np.max(truth), np.max(reco))
        ax.scatter(truth, reco, s=9, alpha=0.16)
        ax.plot([lo, hi], [lo, hi], "--", color="black", linewidth=1.3)
        ax.set_xlabel("true projected coordinate [mm]")
        ax.set_ylabel("cluster-CoG coordinate [mm]")
        ax.set_title(
            f"{title}: bias={1.e3*m['bias']:+.1f} um, "
            f"sigma={1.e3*m['sigma']:.1f} um"
        )
        ax.grid(alpha=0.25)

    fig.suptitle(
        "Stage B3b: each tilted strip family measures one projected coordinate"
    )
    fig.tight_layout(rect=(0., 0., 1., 0.93))
    fig.savefig(
        args.output_prefix + "_projected_coordinates.png",
        dpi=200,
    )
    plt.close(fig)

    # ------------------------------------------------------------
    # Plot 2: final reconstructed u and w.
    # ------------------------------------------------------------
    fig, (axw, axu) = plt.subplots(1, 2, figsize=(11.6, 5.2))

    wtrue = paired["w_true_mm"].to_numpy()
    wreco = paired["w_reco_mm"].to_numpy()
    lo = min(np.min(wtrue), np.min(wreco))
    hi = max(np.max(wtrue), np.max(wreco))
    axw.scatter(wtrue, wreco, s=10, alpha=0.18)
    axw.plot([lo, hi], [lo, hi], "--", color="black")
    axw.set_xlabel("true w [mm]")
    axw.set_ylabel("reconstructed w [mm]")
    axw.set_title(
        f"w: bias={1.e3*mw['bias']:+.1f} um, "
        f"sigma={1.e3*mw['sigma']:.1f} um"
    )
    axw.grid(alpha=0.25)

    utrue = paired["u_true_mm"].to_numpy()
    ureco = paired["u_reco_mm"].to_numpy()
    lo = min(np.min(utrue), np.min(ureco))
    hi = max(np.max(utrue), np.max(ureco))
    axu.scatter(utrue, ureco, s=10, alpha=0.18)
    axu.plot([lo, hi], [lo, hi], "--", color="black")
    axu.set_xlabel("true u [mm]")
    axu.set_ylabel("reconstructed u [mm]")
    axu.set_title(
        f"u: bias={1.e3*mu['bias']:+.1f} um, "
        f"sigma={1.e3*mu['sigma']:.1f} um"
    )
    axu.grid(alpha=0.25)

    fig.suptitle(
        f"Stage B3b: 2D reconstruction from cluster CoG, tan alpha = ±{args.tan_alpha:g}"
    )
    fig.tight_layout(rect=(0., 0., 1., 0.93))
    fig.savefig(
        args.output_prefix + "_reconstructed_vs_true.png",
        dpi=200,
    )
    plt.close(fig)

    # ------------------------------------------------------------
    # Plot 3: compact 2D true points and reconstructed event cloud.
    # ------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(8.4, 6.4))
    ax.scatter(
        paired["u_reco_mm"], paired["w_reco_mm"],
        s=9, alpha=0.18, label="cluster-CoG reconstructed events"
    )
    truth_grid = (
        paired[["u_true_mm", "w_true_mm"]]
        .drop_duplicates()
        .sort_values(["u_true_mm", "w_true_mm"])
    )
    ax.scatter(
        truth_grid["u_true_mm"], truth_grid["w_true_mm"],
        marker="x", s=70, linewidths=1.5, color="black",
        label="true scan points"
    )
    ax.set_xlabel("u [mm]")
    ax.set_ylabel("w [mm]")
    ax.set_title("Stage B3b: two cluster-CoG projections reconstruct the in-plane point")
    ax.grid(alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(
        args.output_prefix + "_2d_reconstruction.png",
        dpi=200,
    )
    plt.close(fig)

    # ------------------------------------------------------------
    # Plot 4: presentation-level resolution summary.
    # ------------------------------------------------------------
    labels = ["x(+alpha)", "x(-alpha)", "w", "u"]
    sigmas = [
        1.e3 * mxp["sigma"],
        1.e3 * mxm["sigma"],
        1.e3 * mw["sigma"],
        1.e3 * mu["sigma"],
    ]

    fig, ax = plt.subplots(figsize=(7.8, 5.6))
    bars = ax.bar(labels, sigmas)
    for bar, value in zip(bars, sigmas):
        ax.text(
            bar.get_x() + bar.get_width() / 2.,
            value + max(sigmas) * 0.02,
            f"{value:.1f}",
            ha="center", va="bottom",
        )
    ax.set_ylabel("event-to-event spread sigma [um]")
    ax.set_title("Stage B3b: cluster-CoG resolution through the two-family reconstruction")
    ax.set_ylim(0., max(sigmas) * 1.15)
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(
        args.output_prefix + "_resolution_summary.png",
        dpi=200,
    )
    plt.close(fig)

    # ------------------------------------------------------------
    # Plot 5: correlation of the two projected-coordinate errors.
    # This diagnoses why the difference coordinate u can reconstruct much
    # better than either individual projected coordinate.
    # ------------------------------------------------------------
    dplus_um = 1.e3 * paired["dx_plus_mm"].to_numpy(dtype=float)
    dminus_um = 1.e3 * paired["dx_minus_mm"].to_numpy(dtype=float)
    rho = float(np.corrcoef(dplus_um, dminus_um)[0, 1])

    common_um = 0.5 * (dplus_um + dminus_um)
    differential_um = 0.5 * (dminus_um - dplus_um)

    fig, ax = plt.subplots(figsize=(7.2, 6.5))
    ax.scatter(
        dplus_um, dminus_um,
        s=10, alpha=0.18,
        label=f"paired avalanches, rho = {rho:.4f}"
    )
    lo = min(np.min(dplus_um), np.min(dminus_um))
    hi = max(np.max(dplus_um), np.max(dminus_um))
    ax.plot([lo, hi], [lo, hi], "--", color="black",
            linewidth=1.4, label="perfect common-mode: delta_- = delta_+")
    ax.axhline(0., linewidth=0.8)
    ax.axvline(0., linewidth=0.8)
    ax.set_xlabel("error in +alpha projected coordinate, delta_+ [um]")
    ax.set_ylabel("error in -alpha projected coordinate, delta_- [um]")
    ax.set_title(
        "Stage B3b: projected-coordinate errors are strongly common-mode"
    )
    ax.grid(alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(
        args.output_prefix + "_projection_error_correlation.png",
        dpi=200,
    )
    plt.close(fig)

    # Save the common/differential decomposition explicitly for later
    # comparison with physical two-cathode readout.
    paired["delta_common_um"] = common_um
    paired["delta_differential_um"] = differential_um
    paired.to_csv(
        args.output_prefix + "_event_reconstruction.csv",
        index=False,
    )

    print("\n=== STAGE B3b: TWO-FAMILY CLUSTER-CoG RECONSTRUCTION ===")
    print(f"selected tan(alpha)      : {tplus:+g}, {tminus:+g}")
    print(f"|alpha|                  : {alpha * 180. / np.pi:.4f} deg")
    print(f"strip pitch              : {p:g} mm")
    print(
        f"exported strips          : {strip_ids[0]} .. {strip_ids[-1]} "
        f"({len(strip_ids)} total)"
    )
    print(
        f"paired events            : {len(paired)} / {total_pairs} "
        f"({100.*reconstructable_fraction:.1f}% reconstructable)"
    )
    print("\nProjected-coordinate reconstruction:")
    print(
        f"  x(+alpha) bias/sigma   : "
        f"{1.e3*mxp['bias']:+.2f} / {1.e3*mxp['sigma']:.2f} um"
    )
    print(
        f"  x(-alpha) bias/sigma   : "
        f"{1.e3*mxm['bias']:+.2f} / {1.e3*mxm['sigma']:.2f} um"
    )
    print("\n2D inversion:")
    print(
        f"  w bias/sigma/RMS       : "
        f"{1.e3*mw['bias']:+.2f} / {1.e3*mw['sigma']:.2f} / "
        f"{1.e3*mw['rms']:.2f} um"
    )
    print(
        f"  u bias/sigma/RMS       : "
        f"{1.e3*mu['bias']:+.2f} / {1.e3*mu['sigma']:.2f} / "
        f"{1.e3*mu['rms']:.2f} um"
    )
    print("\nError correlation:")
    print(f"  rho(delta_+, delta_-)   : {rho:.6f}")
    print(
        f"  sigma common-mode       : "
        f"{np.std(common_um, ddof=1):.2f} um"
    )
    print(
        f"  sigma differential-mode : "
        f"{np.std(differential_um, ddof=1):.2f} um"
    )
    print(
        "\nNOTE: no R calibration is used; each family coordinate is the "
        "cluster center of gravity over all exported strips."
    )
    print(
        "NOTE: physical placement of the two families (same cathode versus "
        "opposite cathodes) remains the next electrostatic study."
    )
    print("\nWrote:")
    print(" ", args.output_prefix + "_event_reconstruction.csv")
    print(" ", args.output_prefix + "_summary.csv")
    print(" ", args.output_prefix + "_projected_coordinates.png")
    print(" ", args.output_prefix + "_reconstructed_vs_true.png")
    print(" ", args.output_prefix + "_2d_reconstruction.png")
    print(" ", args.output_prefix + "_resolution_summary.png")
    print(" ", args.output_prefix + "_projection_error_correlation.png")


if __name__ == "__main__":
    main()
