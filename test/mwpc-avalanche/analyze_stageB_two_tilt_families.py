#!/usr/bin/env python3
"""
Stage B3b: proof-of-principle reconstruction with two ideal strip families.

This analysis uses ONE Stage-B3a-style scan containing BOTH +alpha and -alpha.
The C++ scan reuses the same microscopic Garfield avalanche bank for every
requested tilt angle, wire index and along-wire translation, so rows from the
two families can be paired event by event.

Coordinate convention:
    u : across anode wires, in the detector plane
    w : along anode wires, in the detector plane
    v : normal to the cathode planes

For the two ideal families:
    x_plus  = w cos(alpha) - u sin(alpha)
    x_minus = w cos(alpha) + u sin(alpha)

Each family measures its own projected coordinate through local three-strip
charge sharing.  We define the local asymmetry around the strongest strip

    R = (A_right - A_left) / (A_left + A_center + A_right).

A calibration subset is used to learn R -> x_local/p.  A disjoint validation
subset is then reconstructed event by event:

    x_plus  <- (+alpha family)
    x_minus <- (-alpha family)

and finally

    w = (x_plus + x_minus) / (2 cos alpha)
    u = (x_minus - x_plus) / (2 sin alpha).

IMPORTANT SCOPE:
This is an IDEAL two-family reconstruction test.  The two strip families are
treated as independent readout projections of the same avalanche.  It does NOT
yet choose whether the real detector places both families on the same cathode
or on opposite cathodes, and it does not include inter-layer electrostatic
cross-coupling, electronics thresholds/noise, or finite PCB geometry.
"""

import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


STRIP_IDS = np.array([-2, -1, 0, 1, 2], dtype=int)
A_COLS = [f"A_strip_{k}_fC" for k in STRIP_IDS]


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument(
        "--input",
        default="stageB_two_tilt_families.csv",
        help="phaseB_tilt_scan output containing both +alpha and -alpha.",
    )
    p.add_argument(
        "--output-prefix",
        default="stageB_two_tilt_families",
    )
    p.add_argument(
        "--tan-alpha",
        type=float,
        default=0.10,
        help="Use +/- this tilt magnitude.",
    )
    p.add_argument(
        "--strip-pitch-mm",
        type=float,
        default=1.7,
    )
    p.add_argument(
        "--calibration-parity",
        type=int,
        choices=(0, 1),
        default=0,
        help=(
            "event_index parity used for calibration. "
            "The opposite parity is validation."
        ),
    )
    p.add_argument(
        "--calibration-bins",
        type=int,
        default=23,
        help="Number of bins in local projected coordinate for R calibration.",
    )
    p.add_argument(
        "--max-local-eta",
        type=float,
        default=0.58,
        help="Keep calibration rows with |x_local/p| below this value.",
    )
    return p.parse_args()


def nearest_value(values, target):
    values = np.asarray(sorted(set(values)), dtype=float)
    return float(values[np.argmin(np.abs(values - target))])


def add_local_measurement(df, pitch_mm):
    """
    Convert fixed global strip amplitudes (-2..+2) into a local 3-strip
    observable centred on the strongest strip.

    This mirrors what a real cluster finder would do at first order:
    choose the strongest channel as the local centre, then use its immediate
    neighbours to form a left-right asymmetry.

    Rows whose strongest strip is +/-2 cannot be reconstructed with the current
    CSV because one required neighbour would lie outside the stored range.
    """
    out = df.copy()
    amps = out[A_COLS].to_numpy(dtype=float)

    imax = np.argmax(amps, axis=1)
    k_center = STRIP_IDS[imax]

    valid = np.isin(k_center, [-1, 0, 1])

    a_left = np.full(len(out), np.nan)
    a_center = np.full(len(out), np.nan)
    a_right = np.full(len(out), np.nan)

    for i in np.where(valid)[0]:
        j = imax[i]
        a_left[i] = amps[i, j - 1]
        a_center[i] = amps[i, j]
        a_right[i] = amps[i, j + 1]

    a3 = a_left + a_center + a_right
    r_local = (a_right - a_left) / a3

    x_true = out["projected_x_mm"].to_numpy(dtype=float)
    x_local = x_true - k_center.astype(float) * pitch_mm

    out["local_center_strip"] = k_center
    out["local_A_left_fC"] = a_left
    out["local_A_center_fC"] = a_center
    out["local_A_right_fC"] = a_right
    out["local_A3_fC"] = a3
    out["local_asymmetry"] = r_local
    out["local_x_true_mm"] = x_local
    out["local_eta_true"] = x_local / pitch_mm
    out["local_reconstructable"] = valid & np.isfinite(r_local) & (a3 > 0.)

    return out


def calibration_curve(df, bins, max_eta):
    """
    Build a monotonic shared calibration curve eta_local -> R.

    We bin event-level calibration data, average R in each eta bin, then impose
    a weak monotonicity constraint with cumulative maximum.  This avoids any
    external fitting dependency and is sufficient for the present ideal
    proof-of-principle.
    """
    good = df[
        df["local_reconstructable"]
        & np.isfinite(df["local_eta_true"])
        & np.isfinite(df["local_asymmetry"])
        & (np.abs(df["local_eta_true"]) <= max_eta)
    ].copy()

    if len(good) < 20:
        raise RuntimeError("Too few calibration rows for a stable curve.")

    edges = np.linspace(-max_eta, max_eta, bins + 1)
    good["eta_bin"] = pd.cut(
        good["local_eta_true"],
        bins=edges,
        include_lowest=True,
        labels=False,
    )

    curve = (
        good.dropna(subset=["eta_bin"])
        .groupby("eta_bin", as_index=False)
        .agg(
            eta=("local_eta_true", "mean"),
            R=("local_asymmetry", "mean"),
            R_std=("local_asymmetry", "std"),
            count=("local_asymmetry", "size"),
        )
        .sort_values("eta")
        .reset_index(drop=True)
    )

    if len(curve) < 7:
        raise RuntimeError("Calibration binning left too few populated bins.")

    # Ideal response is increasing.  Statistical fluctuations in a finite
    # calibration sample can create tiny local reversals; remove only those.
    curve["R_monotonic"] = np.maximum.accumulate(
        curve["R"].to_numpy(dtype=float)
    )

    # np.interp requires a strictly increasing x-array for robust inversion.
    r = curve["R_monotonic"].to_numpy(dtype=float)
    eta = curve["eta"].to_numpy(dtype=float)

    keep = np.ones(len(curve), dtype=bool)
    keep[1:] = np.diff(r) > 1.e-8

    inv_r = r[keep]
    inv_eta = eta[keep]

    if len(inv_r) < 5:
        raise RuntimeError("Calibration curve is not sufficiently monotonic.")

    return curve, inv_r, inv_eta


def invert_asymmetry(r_values, inv_r, inv_eta):
    r_values = np.asarray(r_values, dtype=float)
    eta = np.interp(
        r_values,
        inv_r,
        inv_eta,
        left=inv_eta[0],
        right=inv_eta[-1],
    )
    clipped = (r_values < inv_r[0]) | (r_values > inv_r[-1])
    return eta, clipped


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
    if args.calibration_bins < 9:
        raise ValueError("--calibration-bins must be >= 9")

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
        *A_COLS,
    }
    missing = required - set(df.columns)
    if missing:
        raise KeyError(f"Missing columns: {sorted(missing)}")

    tplus = nearest_value(df["tan_alpha"].unique(), +args.tan_alpha)
    tminus = nearest_value(df["tan_alpha"].unique(), -args.tan_alpha)

    if tplus <= 0. or tminus >= 0.:
        raise RuntimeError(
            "Input must contain both positive and negative tilt values."
        )

    plus = add_local_measurement(
        df[np.isclose(df["tan_alpha"], tplus)].copy(),
        args.strip_pitch_mm,
    )
    minus = add_local_measurement(
        df[np.isclose(df["tan_alpha"], tminus)].copy(),
        args.strip_pitch_mm,
    )

    alpha = abs(float(plus["alpha_deg"].iloc[0])) * np.pi / 180.
    cos_a = np.cos(alpha)
    sin_a = np.sin(alpha)

    calibration_parity = args.calibration_parity
    validation_parity = 1 - calibration_parity

    plus_cal = plus[plus["event_index"] % 2 == calibration_parity].copy()
    minus_cal = minus[minus["event_index"] % 2 == calibration_parity].copy()

    # Shared local calibration.  In the ideal geometry +alpha and -alpha
    # should have the same local response once expressed in x_local/p.
    cal_all = pd.concat(
        [
            plus_cal.assign(family="+alpha"),
            minus_cal.assign(family="-alpha"),
        ],
        ignore_index=True,
    )

    curve, inv_r, inv_eta = calibration_curve(
        cal_all,
        args.calibration_bins,
        args.max_local_eta,
    )

    # --------------------------------------------------------------
    # Validation: pair the same Garfield avalanche in the two families.
    # --------------------------------------------------------------
    plus_val = plus[plus["event_index"] % 2 == validation_parity].copy()
    minus_val = minus[minus["event_index"] % 2 == validation_parity].copy()

    keys = [
        "wire_index",
        "wire_u_mm",
        "xi_w",
        "w_shift_mm",
        "event_index",
        "random_seed",
    ]

    pcols = keys + [
        "projected_x_mm",
        "local_center_strip",
        "local_asymmetry",
        "local_reconstructable",
    ]
    mcols = keys + [
        "projected_x_mm",
        "local_center_strip",
        "local_asymmetry",
        "local_reconstructable",
    ]

    paired_all = plus_val[pcols].merge(
        minus_val[mcols],
        on=keys,
        suffixes=("_plus", "_minus"),
        how="inner",
    )

    total_pairs = len(paired_all)

    paired = paired_all[
        paired_all["local_reconstructable_plus"]
        & paired_all["local_reconstructable_minus"]
    ].copy()

    if paired.empty:
        raise RuntimeError("No validation events reconstructable in both families.")

    eta_plus, clip_plus = invert_asymmetry(
        paired["local_asymmetry_plus"].to_numpy(dtype=float),
        inv_r,
        inv_eta,
    )
    eta_minus, clip_minus = invert_asymmetry(
        paired["local_asymmetry_minus"].to_numpy(dtype=float),
        inv_r,
        inv_eta,
    )

    p = args.strip_pitch_mm

    paired["x_local_reco_plus_mm"] = eta_plus * p
    paired["x_local_reco_minus_mm"] = eta_minus * p

    paired["x_reco_plus_mm"] = (
        paired["local_center_strip_plus"].to_numpy(dtype=float) * p
        + paired["x_local_reco_plus_mm"]
    )
    paired["x_reco_minus_mm"] = (
        paired["local_center_strip_minus"].to_numpy(dtype=float) * p
        + paired["x_local_reco_minus_mm"]
    )

    paired["calibration_clipped_plus"] = clip_plus
    paired["calibration_clipped_minus"] = clip_minus

    paired["w_reco_mm"] = (
        paired["x_reco_plus_mm"] + paired["x_reco_minus_mm"]
    ) / (2. * cos_a)

    paired["u_reco_mm"] = (
        paired["x_reco_minus_mm"] - paired["x_reco_plus_mm"]
    ) / (2. * sin_a)

    paired["w_true_mm"] = paired["w_shift_mm"]
    paired["u_true_mm"] = paired["wire_u_mm"]

    paired["dw_mm"] = paired["w_reco_mm"] - paired["w_true_mm"]
    paired["du_mm"] = paired["u_reco_mm"] - paired["u_true_mm"]

    paired.to_csv(
        args.output_prefix + "_event_reconstruction.csv",
        index=False,
    )

    mw = metrics(paired["dw_mm"])
    mu = metrics(paired["du_mm"])

    reconstructable_fraction = (
        len(paired) / total_pairs if total_pairs else 0.
    )
    clipped_fraction = np.mean(
        paired["calibration_clipped_plus"]
        | paired["calibration_clipped_minus"]
    )

    # --------------------------------------------------------------
    # Plot 1: local response calibration.
    # --------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(8.8, 5.8))

    for family, g in cal_all[
        cal_all["local_reconstructable"]
        & (np.abs(cal_all["local_eta_true"]) <= args.max_local_eta)
    ].groupby("family"):
        # Bin only for display.
        edges = np.linspace(
            -args.max_local_eta,
            args.max_local_eta,
            args.calibration_bins + 1,
        )
        gg = g.copy()
        gg["bin"] = pd.cut(
            gg["local_eta_true"],
            edges,
            labels=False,
            include_lowest=True,
        )
        disp = (
            gg.dropna(subset=["bin"])
            .groupby("bin", as_index=False)
            .agg(
                eta=("local_eta_true", "mean"),
                R=("local_asymmetry", "mean"),
            )
            .sort_values("eta")
        )
        ax.plot(
            disp["eta"],
            disp["R"],
            marker="o",
            linewidth=1.2,
            label=family,
        )

    ax.plot(
        curve["eta"],
        curve["R_monotonic"],
        color="black",
        linewidth=2.2,
        label="shared calibration used for inversion",
    )
    ax.axhline(0., linewidth=0.9)
    ax.axvline(0., linewidth=0.9)
    ax.set_xlabel(
        r"local projected coordinate "
        r"$eta_{m local}=x_{m local}/p$"
    )
    ax.set_ylabel(
        r"local asymmetry "
        r"$R=(A_{m right}-A_{m left})/A_3$"
    )
    ax.set_title(
        rf"Stage B3b: +$alpha$ and -$alpha$ families share the same local response"
    )
    ax.grid(alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(
        args.output_prefix + "_local_calibration.png",
        dpi=200,
    )
    plt.close(fig)

    # --------------------------------------------------------------
    # Plot 2: reconstructed versus true u and w.
    # --------------------------------------------------------------
    fig, (axw, axu) = plt.subplots(1, 2, figsize=(11.6, 5.2))

    wmin = min(paired["w_true_mm"].min(), paired["w_reco_mm"].min())
    wmax = max(paired["w_true_mm"].max(), paired["w_reco_mm"].max())
    axw.scatter(
        paired["w_true_mm"],
        paired["w_reco_mm"],
        s=16,
        alpha=0.45,
    )
    axw.plot([wmin, wmax], [wmin, wmax], "--", color="black")
    axw.set_xlabel(r"true $w$ [mm]")
    axw.set_ylabel(r"reconstructed $w$ [mm]")
    axw.set_title(
        rf"$w$: bias={1.e3*mw['bias']:+.1f} $mu$m, "
        rf"$sigma={1.e3*mw['sigma']:.1f}$ $mu$m"
    )
    axw.grid(alpha=0.25)

    umin = min(paired["u_true_mm"].min(), paired["u_reco_mm"].min())
    umax = max(paired["u_true_mm"].max(), paired["u_reco_mm"].max())
    axu.scatter(
        paired["u_true_mm"],
        paired["u_reco_mm"],
        s=16,
        alpha=0.45,
    )
    axu.plot([umin, umax], [umin, umax], "--", color="black")
    axu.set_xlabel(r"true $u$ [mm]")
    axu.set_ylabel(r"reconstructed $u$ [mm]")
    axu.set_title(
        rf"$u$: bias={1.e3*mu['bias']:+.1f} $mu$m, "
        rf"$sigma={1.e3*mu['sigma']:.1f}$ $mu$m"
    )
    axu.grid(alpha=0.25)

    fig.suptitle(
        rf"Stage B3b: ideal two-family 2D reconstruction, "
        rf"$	analpha=pm{args.tan_alpha:g}$"
    )
    fig.tight_layout(rect=(0., 0., 1., 0.93))
    fig.savefig(
        args.output_prefix + "_reconstructed_vs_true.png",
        dpi=200,
    )
    plt.close(fig)

    # --------------------------------------------------------------
    # Plot 3: 2D true grid and reconstructed event cloud.
    # --------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(8.4, 6.4))

    ax.scatter(
        paired["u_reco_mm"],
        paired["w_reco_mm"],
        s=13,
        alpha=0.35,
        label="reconstructed events",
    )

    truth = (
        paired[["u_true_mm", "w_true_mm"]]
        .drop_duplicates()
        .sort_values(["u_true_mm", "w_true_mm"])
    )
    ax.scatter(
        truth["u_true_mm"],
        truth["w_true_mm"],
        marker="x",
        s=70,
        linewidths=1.6,
        color="black",
        label="true scan points",
    )

    ax.set_xlabel(r"$u$ [mm]")
    ax.set_ylabel(r"$w$ [mm]")
    ax.set_title(
        r"Stage B3b: two projected coordinates reconstruct the in-plane point"
    )
    ax.grid(alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(
        args.output_prefix + "_2d_reconstruction.png",
        dpi=200,
    )
    plt.close(fig)

    # --------------------------------------------------------------
    # Compact summary CSV.
    # --------------------------------------------------------------
    summary = pd.DataFrame(
        [
            {
                "tan_alpha_magnitude": args.tan_alpha,
                "alpha_deg": alpha * 180. / np.pi,
                "strip_pitch_mm": p,
                "calibration_event_parity": calibration_parity,
                "validation_event_parity": validation_parity,
                "paired_validation_events_before_local_cut": total_pairs,
                "paired_reconstructable_events": len(paired),
                "reconstructable_fraction": reconstructable_fraction,
                "calibration_clipped_fraction": clipped_fraction,
                "w_bias_mm": mw["bias"],
                "w_sigma_mm": mw["sigma"],
                "w_rms_mm": mw["rms"],
                "u_bias_mm": mu["bias"],
                "u_sigma_mm": mu["sigma"],
                "u_rms_mm": mu["rms"],
            }
        ]
    )
    summary.to_csv(
        args.output_prefix + "_summary.csv",
        index=False,
    )

    print("\n=== STAGE B3b: TWO IDEAL TILTED STRIP FAMILIES ===")
    print(
        f"selected tan(alpha)      : {tplus:+g}, {tminus:+g}"
    )
    print(
        f"|alpha|                  : {alpha * 180. / np.pi:.4f} deg"
    )
    print(f"strip pitch              : {p:g} mm")
    print(
        f"calibration / validation : event parity "
        f"{calibration_parity} / {validation_parity}"
    )
    print(
        f"validation pairs         : {len(paired)} / {total_pairs} "
        f"({100.*reconstructable_fraction:.1f}% reconstructable)"
    )
    print(
        f"calibration clipping     : {100.*clipped_fraction:.2f}%"
    )
    print("\nEvent-by-event reconstruction:")
    print(
        f"  w bias / sigma / RMS   : "
        f"{1.e3*mw['bias']:+.2f} / "
        f"{1.e3*mw['sigma']:.2f} / "
        f"{1.e3*mw['rms']:.2f} um"
    )
    print(
        f"  u bias / sigma / RMS   : "
        f"{1.e3*mu['bias']:+.2f} / "
        f"{1.e3*mu['sigma']:.2f} / "
        f"{1.e3*mu['rms']:.2f} um"
    )
    print(
        "\nNOTE: this is an ideal independent-two-family test; "
        "physical layer placement is the next step."
    )
    print("\nWrote:")
    print(" ", args.output_prefix + "_event_reconstruction.csv")
    print(" ", args.output_prefix + "_summary.csv")
    print(" ", args.output_prefix + "_local_calibration.png")
    print(" ", args.output_prefix + "_reconstructed_vs_true.png")
    print(" ", args.output_prefix + "_2d_reconstruction.png")


if __name__ == "__main__":
    main()
