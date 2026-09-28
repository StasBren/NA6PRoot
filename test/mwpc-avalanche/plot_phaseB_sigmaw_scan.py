#!/usr/bin/env python3
"""Plot a Phase-B1b scan over the synthetic avalanche width sigma_w.

Expected input prefixes are produced by mwpc_phase_b_ion_cloud, for example:
  phaseB_sigmaw_005
  phaseB_sigmaw_020
  phaseB_sigmaw_050
  phaseB_sigmaw_100

Example:
  python3 plot_phaseB_sigmaw_scan.py \
    --runs "0.05=phaseB_sigmaw_005" "0.20=phaseB_sigmaw_020" \
           "0.50=phaseB_sigmaw_050" "1.00=phaseB_sigmaw_100" \
    --output-prefix phaseB_sigmaw_scan

The label before '=' is sigma_w in mm.
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def parse_run(spec):
    if "=" not in spec:
        raise ValueError(f"Run must be SIGMA_MM=PREFIX, got {spec!r}")
    sigma, prefix = spec.split("=", 1)
    return float(sigma), Path(prefix)


def load_run(sigma_mm, prefix):
    summary = pd.read_csv(str(prefix) + "_summary.csv")
    ions = pd.read_csv(str(prefix) + "_ions.csv")
    wave = pd.read_csv(str(prefix) + "_waveforms.csv")

    summary = summary.sort_values("center_w_mm").reset_index(drop=True)

    q = summary["integrated_ion_signal_fC"].to_numpy()
    x = summary["center_w_mm"].to_numpy()
    absq = np.abs(q)
    qabs_sum = absq.sum()

    qfrac = absq / qabs_sum if qabs_sum > 0 else np.zeros_like(absq)
    mean_x = np.sum(qfrac * x)
    rms = np.sqrt(np.sum(qfrac * (x - mean_x) ** 2))

    central = summary.loc[summary["strip"] == 0].iloc[0]
    active_1pct = int(np.sum(qfrac >= 0.01))
    active_5pct = int(np.sum(qfrac >= 0.05))

    # Endpoint distribution in w: useful to connect cloud geometry to strip charge.
    w_end = ions["w1_mm"].to_numpy()
    w_start = ions["w0_mm"].to_numpy()

    # Peak-current spatial width, using peak absolute current as a positive weight.
    ip = summary["peak_abs_current_fC_per_ns"].to_numpy()
    ip_sum = ip.sum()
    ipfrac = ip / ip_sum if ip_sum > 0 else np.zeros_like(ip)
    mean_ip = np.sum(ipfrac * x)
    peak_rms = np.sqrt(np.sum(ipfrac * (x - mean_ip) ** 2))

    return {
        "sigma_mm": sigma_mm,
        "prefix": str(prefix),
        "summary": summary,
        "ions": ions,
        "wave": wave,
        "qfrac": qfrac,
        "charge_rms_mm": rms,
        "peak_rms_mm": peak_rms,
        "central_fraction": float(qfrac[summary["strip"].to_numpy() == 0][0]),
        "active_1pct": active_1pct,
        "active_5pct": active_5pct,
        "w_start_sigma_mm": float(np.std(w_start, ddof=1)),
        "w_end_sigma_mm": float(np.std(w_end, ddof=1)),
        "central_q_fC": float(central["integrated_ion_signal_fC"]),
        "central_peak_fC_per_ns": float(central["peak_abs_current_fC_per_ns"]),
        "ramo_rel_error": float(
            (central["integrated_ion_signal_fC"] - central["ramo_endpoint_fC"])
            / central["ramo_endpoint_fC"]
        ) if central["ramo_endpoint_fC"] != 0 else np.nan,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", nargs="+", required=True,
                    help="SIGMA_W_MM=PREFIX entries.")
    ap.add_argument("--output-prefix", default="phaseB_sigmaw_scan")
    args = ap.parse_args()

    runs = sorted([load_run(*parse_run(x)) for x in args.runs],
                  key=lambda r: r["sigma_mm"])

    # 1) Integrated charge sharing profile.
    plt.figure(figsize=(8, 5))
    for r in runs:
        s = r["summary"]
        plt.plot(s["center_w_mm"], 100. * r["qfrac"],
                 marker="o", label=f"sigma_w = {r['sigma_mm']:.2f} mm")
    plt.xlabel("Strip center w [mm]")
    plt.ylabel("|Integrated charge| fraction [%]")
    plt.title("Phase B1b: charge sharing vs synthetic avalanche width")
    plt.legend()
    plt.tight_layout()
    plt.savefig(args.output_prefix + "_charge_sharing.png", dpi=200)
    plt.close()

    # 2) Spatial width of the charge and peak-current profiles.
    sigmas = np.array([r["sigma_mm"] for r in runs])
    charge_rms = np.array([r["charge_rms_mm"] for r in runs])
    peak_rms = np.array([r["peak_rms_mm"] for r in runs])

    plt.figure(figsize=(8, 5))
    plt.plot(sigmas, charge_rms, marker="o", label="Integrated-charge RMS")
    plt.plot(sigmas, peak_rms, marker="o", label="Peak-current RMS")
    plt.xlabel("Input cloud sigma_w [mm]")
    plt.ylabel("Strip-profile RMS width [mm]")
    plt.title("Phase B1b: input cloud width -> readout-cluster width")
    plt.legend()
    plt.tight_layout()
    plt.savefig(args.output_prefix + "_cluster_width.png", dpi=200)
    plt.close()

    # 3) Central-strip fraction.
    central = np.array([100. * r["central_fraction"] for r in runs])

    plt.figure(figsize=(8, 5))
    plt.plot(sigmas, central, marker="o")
    plt.xlabel("Input cloud sigma_w [mm]")
    plt.ylabel("Central-strip share of |integrated charge| [%]")
    plt.title("Phase B1b: central-strip charge fraction")
    plt.tight_layout()
    plt.savefig(args.output_prefix + "_central_fraction.png", dpi=200)
    plt.close()

    # 4) Verify w-spread is preserved by this translationally invariant field.
    w0 = np.array([r["w_start_sigma_mm"] for r in runs])
    w1 = np.array([r["w_end_sigma_mm"] for r in runs])

    plt.figure(figsize=(8, 5))
    plt.plot(sigmas, w0, marker="o", label="Ion cloud at start")
    plt.plot(sigmas, w1, marker="o", label="Ion endpoints")
    plt.xlabel("Requested sigma_w [mm]")
    plt.ylabel("Measured sigma(w) [mm]")
    plt.title("Phase B1b: cloud width along the wire")
    plt.legend()
    plt.tight_layout()
    plt.savefig(args.output_prefix + "_w_width_check.png", dpi=200)
    plt.close()

    rows = []
    for r in runs:
        rows.append({
            "sigma_w_input_mm": r["sigma_mm"],
            "sigma_w_start_mm": r["w_start_sigma_mm"],
            "sigma_w_end_mm": r["w_end_sigma_mm"],
            "charge_profile_rms_mm": r["charge_rms_mm"],
            "peak_current_profile_rms_mm": r["peak_rms_mm"],
            "central_charge_fraction": r["central_fraction"],
            "strips_ge_1pct": r["active_1pct"],
            "strips_ge_5pct": r["active_5pct"],
            "central_Q_fC": r["central_q_fC"],
            "central_peak_fC_per_ns": r["central_peak_fC_per_ns"],
            "central_Ramo_rel_error": r["ramo_rel_error"],
        })

    table = pd.DataFrame(rows)
    table.to_csv(args.output_prefix + "_summary.csv", index=False)

    print(table.to_string(index=False))
    print("\nWrote:")
    for suffix in [
        "_charge_sharing.png",
        "_cluster_width.png",
        "_central_fraction.png",
        "_w_width_check.png",
        "_summary.csv",
    ]:
        print(" ", args.output_prefix + suffix)


if __name__ == "__main__":
    main()
