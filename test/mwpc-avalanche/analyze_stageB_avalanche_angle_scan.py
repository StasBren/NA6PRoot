#!/usr/bin/env python3
"""Analyse avalanche localisation and cathode coupling versus wire azimuth."""

import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--input", default="stageB_avalanche_angle_scan.csv")
    p.add_argument(
        "--output-prefix",
        default="stageB_avalanche_angle_scan",
    )
    p.add_argument("--wire-diam-um", type=float, default=30.0)
    return p.parse_args()


def circ_diff_deg(a, b):
    return (np.asarray(a) - np.asarray(b) + 180.) % 360. - 180.


def sem(x):
    x = np.asarray(x, dtype=float)
    return np.std(x, ddof=1) / np.sqrt(len(x)) if len(x) > 1 else 0.0


def main():
    args = parse_args()
    df = pd.read_csv(args.input)
    if df.empty:
        raise RuntimeError("Input CSV is empty.")

    theta_order = np.sort(df["theta_seed_deg"].unique())
    grouped = df.groupby("theta_seed_deg", sort=True)

    rows = []
    for theta, g in grouped:
        rows.append({
            "theta_seed_deg": theta,
            "events": len(g),
            "mean_birth_theta_deg":
                np.rad2deg(np.arctan2(
                    np.mean(np.sin(np.deg2rad(g["birth_mean_theta_deg"]))),
                    np.mean(np.cos(np.deg2rad(g["birth_mean_theta_deg"]))),
                )),
            "mean_abs_birth_theta_minus_seed_deg":
                float(np.mean(np.abs(g["birth_theta_minus_seed_deg"]))),
            "mean_birth_theta_sigma_deg":
                float(g["birth_theta_sigma_deg"].mean()),
            "mean_birth_radius_um":
                float(g["birth_mean_radius_um"].mean()),
            "mean_minus_fraction":
                float(g["abs_minus_fraction_of_two_cathodes"].mean()),
            "sem_minus_fraction":
                sem(g["abs_minus_fraction_of_two_cathodes"]),
            "mean_plus_fraction":
                float(g["abs_plus_fraction_of_two_cathodes"].mean()),
            "mean_abs_q_minus_fC":
                float(np.mean(np.abs(g["q_minus_total_fC"]))),
            "mean_abs_q_plus_fC":
                float(np.mean(np.abs(g["q_plus_total_fC"]))),
            "median_garfield_ions":
                float(np.median(g["garfield_ions"])),
            "mean_garfield_ions":
                float(np.mean(g["garfield_ions"])),
        })

    summary = pd.DataFrame(rows)
    summary.to_csv(args.output_prefix + "_summary_by_angle.csv", index=False)

    # 1. Does the avalanche stay localised around the incoming azimuth?
    fig, ax = plt.subplots(figsize=(7.3, 6.3))
    ax.scatter(
        df["theta_seed_deg"],
        df["birth_mean_theta_deg"],
        s=14, alpha=0.25,
        label="individual avalanches"
    )
    ax.plot([0, 360], [0, 360], "--", color="black", label="same azimuth")
    ax.set_xlim(-10, 370)
    ax.set_ylim(-190, 190)
    ax.set_xlabel("seed azimuth around wire [deg]")
    ax.set_ylabel("mean ionisation-birth azimuth [deg]")
    ax.set_title("Stage B3c: where around the wire does the avalanche form?")
    ax.grid(alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(args.output_prefix + "_birth_angle_mapping.png", dpi=200)
    plt.close(fig)

    # A cleaner circular-error plot avoids the +/-180 degree wrap.
    fig, ax = plt.subplots(figsize=(8.0, 5.2))
    vals = []
    means = []
    errs = []
    for theta in theta_order:
        g = df[df["theta_seed_deg"] == theta]
        d = circ_diff_deg(g["birth_mean_theta_deg"], theta)
        vals.append(theta)
        means.append(np.mean(d))
        errs.append(np.std(d, ddof=1) if len(d) > 1 else 0.)
    ax.errorbar(vals, means, yerr=errs, marker="o", capsize=3)
    ax.axhline(0., color="black", linestyle="--")
    ax.set_xlabel("seed azimuth around wire [deg]")
    ax.set_ylabel("mean avalanche azimuth - seed azimuth [deg]")
    ax.set_title("Stage B3c: avalanche localisation around the anode wire")
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(args.output_prefix + "_birth_angle_error.png", dpi=200)
    plt.close(fig)

    # 2. Mean ionisation-birth centroid in the u-v cross-section.
    fig, ax = plt.subplots(figsize=(6.7, 6.5))
    sc = ax.scatter(
        df["birth_mean_u_um"], df["birth_mean_v_um"],
        c=df["theta_seed_deg"], s=18, alpha=0.5
    )
    wire_r = 0.5 * args.wire_diam_um
    phi = np.linspace(0., 2. * np.pi, 300)
    ax.plot(wire_r * np.cos(phi), wire_r * np.sin(phi), color="black")
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("u [um]")
    ax.set_ylabel("v [um]")
    ax.set_title("Stage B3c: avalanche ionisation centroid around one wire")
    cb = fig.colorbar(sc, ax=ax)
    cb.set_label("seed azimuth [deg]")
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(args.output_prefix + "_birth_centroid_uv.png", dpi=200)
    plt.close(fig)

    # 3. Cathode sharing versus avalanche azimuth.
    fig, ax = plt.subplots(figsize=(8.2, 5.3))
    ax.errorbar(
        summary["theta_seed_deg"],
        100. * summary["mean_minus_fraction"],
        yerr=100. * summary["sem_minus_fraction"],
        marker="o", capsize=3, label="minus cathode"
    )
    ax.plot(
        summary["theta_seed_deg"],
        100. * summary["mean_plus_fraction"],
        marker="o", label="plus cathode"
    )
    ax.set_xlabel("seed azimuth around wire [deg]")
    ax.set_ylabel("fraction of |Q-| + |Q+| [%]")
    ax.set_title("Stage B3c: cathode signal sharing versus avalanche azimuth")
    ax.grid(alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(args.output_prefix + "_cathode_fraction_vs_angle.png", dpi=200)
    plt.close(fig)

    # 4. Absolute signal amplitudes.
    fig, ax = plt.subplots(figsize=(8.2, 5.3))
    ax.plot(
        summary["theta_seed_deg"],
        summary["mean_abs_q_minus_fC"],
        marker="o", label="minus cathode"
    )
    ax.plot(
        summary["theta_seed_deg"],
        summary["mean_abs_q_plus_fC"],
        marker="o", label="plus cathode"
    )
    ax.set_xlabel("seed azimuth around wire [deg]")
    ax.set_ylabel("mean |Q(Tobs)| [fC]")
    ax.set_title("Stage B3c: complete-cathode signal versus avalanche azimuth")
    ax.grid(alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(args.output_prefix + "_cathode_signal_vs_angle.png", dpi=200)
    plt.close(fig)

    # 5. Gain diagnostic: this tells us whether some angles preferentially
    # produce a different avalanche size in the asymmetric cell.
    fig, ax = plt.subplots(figsize=(8.2, 5.3))
    ax.plot(
        summary["theta_seed_deg"],
        summary["median_garfield_ions"],
        marker="o", label="median avalanche ions"
    )
    ax.set_xlabel("seed azimuth around wire [deg]")
    ax.set_ylabel("Garfield positive ions per avalanche")
    ax.set_title("Stage B3c diagnostic: gain versus seed azimuth")
    ax.grid(alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(args.output_prefix + "_gain_vs_angle.png", dpi=200)
    plt.close(fig)

    print("\n=== STAGE B3c: AVALANCHE ANGLE-SCAN ANALYSIS ===")
    print(f"events                  : {len(df)}")
    print(
        "mean |birth angle - seed angle| : "
        f"{np.mean(np.abs(df['birth_theta_minus_seed_deg'])):.2f} deg"
    )
    print(
        "mean intrinsic birth angular spread: "
        f"{df['birth_theta_sigma_deg'].mean():.2f} deg"
    )
    print(
        "minus-fraction range    : "
        f"{100*summary['mean_minus_fraction'].min():.2f}% .. "
        f"{100*summary['mean_minus_fraction'].max():.2f}%"
    )
    print(
        "plus-fraction range     : "
        f"{100*summary['mean_plus_fraction'].min():.2f}% .. "
        f"{100*summary['mean_plus_fraction'].max():.2f}%"
    )
    print("\nWrote:")
    for suffix in [
        "_summary_by_angle.csv",
        "_birth_angle_mapping.png",
        "_birth_angle_error.png",
        "_birth_centroid_uv.png",
        "_cathode_fraction_vs_angle.png",
        "_cathode_signal_vs_angle.png",
        "_gain_vs_angle.png",
    ]:
        print(" ", args.output_prefix + suffix)


if __name__ == "__main__":
    main()
