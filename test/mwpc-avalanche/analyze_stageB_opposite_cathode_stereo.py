#!/usr/bin/env python3
"""Analyse Stage B3c.2 opposite-cathode stereo hybrid output."""

import argparse
import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument(
        "--input",
        default="stageB_opposite_cathode_stereo.csv",
    )
    p.add_argument(
        "--output-prefix",
        default="stageB_opposite_cathode_stereo",
    )
    p.add_argument("--strip-pitch-mm", type=float, default=1.7)
    return p.parse_args()


def discover(df, family):
    pat = re.compile(rf"^A_{family}_strip_(-?[0-9]+)_fC$")
    found = []
    for c in df.columns:
        m = pat.match(c)
        if m:
            found.append((int(m.group(1)), c))
    found.sort(key=lambda x: x[0])
    if not found:
        raise RuntimeError(f"No amplitudes found for family '{family}'.")
    return np.array([x[0] for x in found], dtype=int), [x[1] for x in found]


def add_cog(df, family, strip_ids, amp_cols, pitch):
    amps = df[amp_cols].to_numpy(dtype=float)
    xk = strip_ids.astype(float) * pitch
    s = np.sum(amps, axis=1)
    good = np.isfinite(s) & (s > 0.)
    cog = np.full(len(df), np.nan)
    cog[good] = np.sum(amps[good] * xk[None, :], axis=1) / s[good]
    df[f"x_{family}_cog_mm"] = cog
    df[f"{family}_reconstructable"] = good
    df[f"{family}_cluster_abs_fC"] = s
    return df


def metrics(x):
    x = np.asarray(x, dtype=float)
    return {
        "bias": float(np.mean(x)),
        "sigma": float(np.std(x, ddof=1)),
        "rms": float(np.sqrt(np.mean(x*x))),
    }


def main():
    args = parse_args()
    df = pd.read_csv(args.input)
    if df.empty:
        raise RuntimeError("Input is empty.")

    ids_m, cols_m = discover(df, "minus")
    ids_p, cols_p = discover(df, "plus")
    if not np.array_equal(ids_m, ids_p):
        raise RuntimeError("Minus and plus strip index sets differ.")

    df = add_cog(df, "minus", ids_m, cols_m, args.strip_pitch_mm)
    df = add_cog(df, "plus", ids_p, cols_p, args.strip_pitch_mm)
    df = df[df["minus_reconstructable"] & df["plus_reconstructable"]].copy()

    alpha = np.deg2rad(float(df["alpha_deg"].iloc[0]))
    ca = np.cos(alpha)
    sa = np.sin(alpha)

    # +alpha family is on minus cathode.
    df["dx_plusalpha_mm"] = (
        df["x_minus_cog_mm"] - df["x_plusalpha_true_mm"]
    )
    # -alpha family is on plus cathode.
    df["dx_minusalpha_mm"] = (
        df["x_plus_cog_mm"] - df["x_minusalpha_true_mm"]
    )

    df["w_true_mm"] = df["w_shift_mm"]
    df["u_true_mm"] = df["wire_u_mm"]
    df["w_reco_mm"] = (
        df["x_minus_cog_mm"] + df["x_plus_cog_mm"]
    ) / (2. * ca)
    df["u_reco_mm"] = (
        df["x_plus_cog_mm"] - df["x_minus_cog_mm"]
    ) / (2. * sa)

    df["dw_mm"] = df["w_reco_mm"] - df["w_true_mm"]
    df["du_mm"] = df["u_reco_mm"] - df["u_true_mm"]

    mp = metrics(df["dx_plusalpha_mm"])
    mm = metrics(df["dx_minusalpha_mm"])
    mw = metrics(df["dw_mm"])
    mu = metrics(df["du_mm"])

    dplus = 1.e3 * df["dx_plusalpha_mm"].to_numpy()
    dminus = 1.e3 * df["dx_minusalpha_mm"].to_numpy()
    rho = float(np.corrcoef(dplus, dminus)[0, 1])
    common = 0.5 * (dplus + dminus)
    differential = 0.5 * (dminus - dplus)

    qminus = np.abs(df["q_cathode_minus_full_fC"].to_numpy())
    qplus = np.abs(df["q_cathode_plus_full_fC"].to_numpy())
    rho_amp = float(np.corrcoef(qminus, qplus)[0, 1])

    df["delta_common_um"] = common
    df["delta_differential_um"] = differential
    df.to_csv(args.output_prefix + "_event_reconstruction.csv", index=False)

    summary = pd.DataFrame([{
        "events": len(df),
        "tan_alpha": float(df["tan_alpha"].iloc[0]),
        "alpha_deg": float(df["alpha_deg"].iloc[0]),
        "seed_side": str(df["seed_side"].iloc[0]),
        "x_plusalpha_bias_um": 1.e3 * mp["bias"],
        "x_plusalpha_sigma_um": 1.e3 * mp["sigma"],
        "x_minusalpha_bias_um": 1.e3 * mm["bias"],
        "x_minusalpha_sigma_um": 1.e3 * mm["sigma"],
        "w_bias_um": 1.e3 * mw["bias"],
        "w_sigma_um": 1.e3 * mw["sigma"],
        "u_bias_um": 1.e3 * mu["bias"],
        "u_sigma_um": 1.e3 * mu["sigma"],
        "projection_error_rho": rho,
        "common_mode_sigma_um": float(np.std(common, ddof=1)),
        "differential_mode_sigma_um": float(np.std(differential, ddof=1)),
        "physical_cathode_amplitude_rho": rho_amp,
        "mean_abs_q_minus_fC": float(np.mean(qminus)),
        "mean_abs_q_plus_fC": float(np.mean(qplus)),
        "median_abs_qplus_over_qminus":
            float(np.median(qplus / np.maximum(qminus, 1.e-30))),
    }])
    summary.to_csv(args.output_prefix + "_summary.csv", index=False)

    # 1. Projected coordinates.
    fig, axs = plt.subplots(1, 2, figsize=(11.5, 5.2))
    for ax, truth, reco, title, m in [
        (
            axs[0], df["x_plusalpha_true_mm"], df["x_minus_cog_mm"],
            "+alpha on minus cathode", mp
        ),
        (
            axs[1], df["x_minusalpha_true_mm"], df["x_plus_cog_mm"],
            "-alpha on plus cathode", mm
        ),
    ]:
        lo = min(truth.min(), reco.min())
        hi = max(truth.max(), reco.max())
        ax.scatter(truth, reco, s=8, alpha=0.15)
        ax.plot([lo, hi], [lo, hi], "--", color="black")
        ax.set_xlabel("true projected coordinate [mm]")
        ax.set_ylabel("cluster-CoG coordinate [mm]")
        ax.set_title(
            f"{title}\nbias={1.e3*m['bias']:+.1f} um, "
            f"sigma={1.e3*m['sigma']:.1f} um"
        )
        ax.grid(alpha=0.25)
    fig.suptitle(
        "Stage B3c.2: opposite cathodes provide the two stereo projections"
    )
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    fig.savefig(args.output_prefix + "_projected_coordinates.png", dpi=200)
    plt.close(fig)

    # 2. Physical-coordinate reconstruction.
    fig, axs = plt.subplots(1, 2, figsize=(11.5, 5.2))
    for ax, truth, reco, label, m in [
        (axs[0], df["w_true_mm"], df["w_reco_mm"], "w", mw),
        (axs[1], df["u_true_mm"], df["u_reco_mm"], "u", mu),
    ]:
        lo = min(truth.min(), reco.min())
        hi = max(truth.max(), reco.max())
        ax.scatter(truth, reco, s=8, alpha=0.15)
        ax.plot([lo, hi], [lo, hi], "--", color="black")
        ax.set_xlabel(f"true {label} [mm]")
        ax.set_ylabel(f"reconstructed {label} [mm]")
        ax.set_title(
            f"{label}: bias={1.e3*m['bias']:+.1f} um, "
            f"sigma={1.e3*m['sigma']:.1f} um"
        )
        ax.grid(alpha=0.25)
    fig.suptitle(
        "Stage B3c.2: 2D reconstruction with opposite-cathode stereo families"
    )
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    fig.savefig(args.output_prefix + "_reconstructed_vs_true.png", dpi=200)
    plt.close(fig)

    # 3. Error correlation.
    fig, ax = plt.subplots(figsize=(7.1, 6.4))
    ax.scatter(dplus, dminus, s=9, alpha=0.16,
               label=f"rho = {rho:.4f}")
    lo = min(dplus.min(), dminus.min())
    hi = max(dplus.max(), dminus.max())
    ax.plot([lo, hi], [lo, hi], "--", color="black",
            label="perfect common-mode")
    ax.axhline(0., linewidth=0.8)
    ax.axvline(0., linewidth=0.8)
    ax.set_xlabel("error in +alpha / minus-cathode coordinate [um]")
    ax.set_ylabel("error in -alpha / plus-cathode coordinate [um]")
    ax.set_title(
        "Stage B3c.2: do opposite cathodes preserve common-mode CoG fluctuations?"
    )
    ax.grid(alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(args.output_prefix + "_projection_error_correlation.png", dpi=200)
    plt.close(fig)

    # 4. Resolution + amplitude summary.
    fig, axs = plt.subplots(1, 2, figsize=(11.0, 5.0))
    labels = ["x(+a)", "x(-a)", "w", "u"]
    vals = [
        1.e3 * mp["sigma"], 1.e3 * mm["sigma"],
        1.e3 * mw["sigma"], 1.e3 * mu["sigma"],
    ]
    bars = axs[0].bar(labels, vals)
    for b, v in zip(bars, vals):
        axs[0].text(
            b.get_x() + b.get_width()/2, v + max(vals)*0.02,
            f"{v:.1f}", ha="center"
        )
    axs[0].set_ylabel("event-to-event spread sigma [um]")
    axs[0].set_title("coordinate resolution")
    axs[0].grid(axis="y", alpha=0.25)

    axs[1].scatter(qminus, qplus, s=8, alpha=0.15)
    axs[1].set_xlabel("|Q| on minus cathode [fC]")
    axs[1].set_ylabel("|Q| on plus cathode [fC]")
    axs[1].set_title(f"physical cathode coupling, rho={rho_amp:.4f}")
    axs[1].grid(alpha=0.25)

    fig.suptitle(
        "Stage B3c.2: physical coupling and stereo reconstruction"
    )
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    fig.savefig(args.output_prefix + "_summary.png", dpi=200)
    plt.close(fig)

    print("\n=== STAGE B3c.2: OPPOSITE-CATHODE STEREO ANALYSIS ===")
    print(f"events                         : {len(df)}")
    print(f"seed side                      : {df['seed_side'].iloc[0]}")
    print(f"projection error rho           : {rho:.6f}")
    print(
        f"common / differential sigma    : "
        f"{np.std(common, ddof=1):.2f} / "
        f"{np.std(differential, ddof=1):.2f} um"
    )
    print(
        f"x(+alpha) bias/sigma           : "
        f"{1.e3*mp['bias']:+.2f} / {1.e3*mp['sigma']:.2f} um"
    )
    print(
        f"x(-alpha) bias/sigma           : "
        f"{1.e3*mm['bias']:+.2f} / {1.e3*mm['sigma']:.2f} um"
    )
    print(
        f"w bias/sigma                   : "
        f"{1.e3*mw['bias']:+.2f} / {1.e3*mw['sigma']:.2f} um"
    )
    print(
        f"u bias/sigma                   : "
        f"{1.e3*mu['bias']:+.2f} / {1.e3*mu['sigma']:.2f} um"
    )
    print(f"physical cathode amplitude rho : {rho_amp:.6f}")
    print(
        f"median |Q+|/|Q-|               : "
        f"{np.median(qplus / np.maximum(qminus, 1.e-30)):.4f}"
    )
    print("\nWrote:")
    print(" ", args.output_prefix + "_event_reconstruction.csv")
    print(" ", args.output_prefix + "_summary.csv")
    print(" ", args.output_prefix + "_projected_coordinates.png")
    print(" ", args.output_prefix + "_reconstructed_vs_true.png")
    print(" ", args.output_prefix + "_projection_error_correlation.png")
    print(" ", args.output_prefix + "_summary.png")


if __name__ == "__main__":
    main()
