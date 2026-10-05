#!/usr/bin/env python3
import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import re
import numpy as np
import pandas as pd

from analyze_stageB_two_tilt_families import (
    add_local_measurement,
    calibration_curve, invert_asymmetry, nearest_value,
)

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


def metrics(residual):
    r = np.asarray(residual, dtype=float)
    return np.mean(r), np.std(r, ddof=1), np.sqrt(np.mean(r * r))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", default="stageB_two_tilt_families.csv")
    ap.add_argument("--output-prefix", default="stageB_cog_vs_calibration")
    ap.add_argument("--tan-alpha", type=float, default=0.10)
    ap.add_argument("--strip-pitch-mm", type=float, default=1.7)
    ap.add_argument("--calibration-parity", type=int, choices=(0, 1), default=0)
    ap.add_argument("--calibration-bins", type=int, default=23)
    ap.add_argument("--max-local-eta", type=float, default=0.58)
    args = ap.parse_args()

    df = pd.read_csv(args.input)
    t = nearest_value(df["tan_alpha"].unique(), args.tan_alpha)
    strip_ids, amp_cols = discover_amplitude_columns(df)

    # The shared helper currently expects at least the central five columns,
    # which remain present in the expanded CSV.
    fam = add_local_measurement(
        df[np.isclose(df["tan_alpha"], t)].copy(),
        args.strip_pitch_mm,
    )

    cal = fam[fam["event_index"] % 2 == args.calibration_parity].copy()
    val = fam[fam["event_index"] % 2 != args.calibration_parity].copy()

    curve, inv_r, inv_eta = calibration_curve(
        cal, args.calibration_bins, args.max_local_eta
    )

    val = val[val["local_reconstructable"]].copy()
    p = args.strip_pitch_mm

    # Proposal-aligned three-strip cluster CoG around the strongest strip.
    # With equal pitch this is exactly x_center + p*R.
    x_center = val["local_center_strip"].to_numpy(dtype=float) * p
    val["x_cog3_mm"] = x_center + p * val["local_asymmetry"]

    # CoG over the central five strips, kept for direct comparison with the
    # earlier study.
    ids5 = np.array([-2, -1, 0, 1, 2], dtype=int)
    cols5 = [f"A_strip_{k}_fC" for k in ids5]
    missing5 = [col for col in cols5 if col not in val.columns]
    if missing5:
        raise RuntimeError(f"Missing central-five amplitude columns: {missing5}")
    amps5 = val[cols5].to_numpy(dtype=float)
    x5 = ids5.astype(float) * p
    val["x_cog5_mm"] = (
        np.sum(amps5 * x5[None, :], axis=1) /
        np.sum(amps5, axis=1)
    )

    # Proposal-aligned cluster CoG using every strip exported by the current
    # simulation.  With --half-strips 7 this is a 15-strip centroid.
    amps_all = val[amp_cols].to_numpy(dtype=float)
    x_all = strip_ids.astype(float) * p
    val["x_cog_cluster_mm"] = (
        np.sum(amps_all * x_all[None, :], axis=1) /
        np.sum(amps_all, axis=1)
    )

    # Existing calibrated R -> local position estimator.
    eta_cal, clipped = invert_asymmetry(
        val["local_asymmetry"].to_numpy(dtype=float),
        inv_r, inv_eta
    )
    val["x_calibrated_mm"] = x_center + p * eta_cal
    val["calibration_clipped"] = clipped

    truth = val["projected_x_mm"].to_numpy(dtype=float)
    methods = {
        "3-strip CoG": val["x_cog3_mm"].to_numpy(),
        "5-strip CoG": val["x_cog5_mm"].to_numpy(),
        "cluster CoG": val["x_cog_cluster_mm"].to_numpy(),
        "calibrated 3-strip R": val["x_calibrated_mm"].to_numpy(),
    }

    rows = []
    for name, reco in methods.items():
        b, s, rms = metrics(reco - truth)
        rows.append({
            "method": name,
            "bias_um": 1.e3 * b,
            "sigma_um": 1.e3 * s,
            "rms_um": 1.e3 * rms,
        })
    pd.DataFrame(rows).to_csv(
        args.output_prefix + "_summary.csv", index=False
    )
    val.to_csv(args.output_prefix + "_events.csv", index=False)

    # 1) Show why R is already a three-strip CoG coordinate.
    fig, ax = plt.subplots(figsize=(8.8, 6.0))
    good = val[np.abs(val["local_eta_true"]) <= args.max_local_eta]
    ax.scatter(
        good["local_eta_true"], good["local_asymmetry"],
        s=12, alpha=0.18, label="validation avalanches"
    )
    ax.plot(
        curve["eta"], curve["R_monotonic"],
        linewidth=2.3, label="mean simulated response R(eta)"
    )
    q = np.linspace(-args.max_local_eta, args.max_local_eta, 200)
    ax.plot(q, q, "--", linewidth=1.7,
            label="raw 3-strip CoG: eta_CoG = R")
    ax.axhline(0., linewidth=0.8)
    ax.axvline(0., linewidth=0.8)
    ax.set_xlabel("true local projected position eta = x_local / p")
    ax.set_ylabel("measured three-strip asymmetry R")
    ax.set_title("Stage B3b: CoG and R are the same three-strip observable")
    ax.grid(alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(args.output_prefix + "_R_vs_true_position.png", dpi=200)
    plt.close(fig)

    # 2) Residuals show whether calibration actually helps.
    fig, ax = plt.subplots(figsize=(9.2, 6.0))
    for name, reco in methods.items():
        ax.scatter(truth, 1.e3 * (reco - truth),
                   s=13, alpha=0.28, label=name)
    ax.axhline(0., color="black", linewidth=1.0)
    ax.set_xlabel("true projected coordinate x [mm]")
    ax.set_ylabel("reconstructed - true [um]")
    ax.set_title("Stage B3b: raw CoG versus response-calibrated position")
    ax.grid(alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(args.output_prefix + "_residuals_vs_position.png", dpi=200)
    plt.close(fig)

    print("\n=== STAGE B3b: CoG VS CALIBRATED ESTIMATOR ===")
    print(f"selected family : tan(alpha) = {t:+g}")
    print(f"strip pitch     : {p:g} mm")
    print(f"validation rows : {len(val)}")
    print(
        f"exported strips : {strip_ids[0]} .. {strip_ids[-1]} "
        f"({len(strip_ids)} strips)"
    )
    print("\nResidual metrics:")
    for row in rows:
        print(
            f"  {row['method']:18s} "
            f"bias={row['bias_um']:+8.2f} um  "
            f"sigma={row['sigma_um']:7.2f} um  "
            f"RMS={row['rms_um']:7.2f} um"
        )
    print(
        "\nIdentity used: x_CoG3 = x_center + p*R, "
        "R=(A_right-A_left)/A3."
    )
    print(
        "The cluster CoG uses every strip exported by the simulation. "
        "Repeat with a larger --half-strips value if a convergence check "
        "of the far tails is needed."
    )


if __name__ == "__main__":
    main()
