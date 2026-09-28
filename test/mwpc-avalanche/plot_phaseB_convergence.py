#!/usr/bin/env python3
"""Compare Phase-B1a single-ion convergence runs.

Example:
  python3 plot_phaseB_convergence.py \
    --runs "auto50=conv_auto50" "auto10=conv_auto10" \
           "step0.10=conv_step0p10" "step0.05=conv_step0p05" \
           "step0.02=conv_step0p02" \
    --output-prefix phaseB_convergence
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def parse_run(spec):
    if "=" not in spec:
        raise ValueError(f"Run must be LABEL=PREFIX, got {spec!r}")
    label, prefix = spec.split("=", 1)
    return label, Path(prefix)


def load_run(label, prefix):
    wave = pd.read_csv(str(prefix) + "_waveforms.csv")
    summary = pd.read_csv(str(prefix) + "_summary.csv")
    drift = pd.read_csv(str(prefix) + "_driftline.csv")
    central = summary.loc[summary["strip"] == 0].iloc[0]

    dt = float(np.median(np.diff(wave["time_ns"])))
    current = wave["strip_0_fC_per_ns"].to_numpy()
    cumulative = np.cumsum(current * dt)

    q_int = float(central["integrated_ion_signal_fC"])
    q_ramo = float(central["ramo_delta_phi_fC"])
    rel_err = np.nan if q_ramo == 0 else (q_int - q_ramo) / q_ramo

    positive_steps = drift.loc[drift["step_mm"] > 0, "step_mm"]
    return {
        "label": label,
        "prefix": str(prefix),
        "wave": wave,
        "summary": summary,
        "drift": drift,
        "dt_ns": dt,
        "cumulative": cumulative,
        "q_int": q_int,
        "q_ramo": q_ramo,
        "rel_err": rel_err,
        "peak": float(central["peak_abs_current_fC_per_ns"]),
        "n_points": len(drift),
        "drift_time_ns": float(drift["time_ns"].iloc[-1] - drift["time_ns"].iloc[0]),
        "max_actual_step_mm": float(positive_steps.max()) if len(positive_steps) else 0.,
        "median_actual_step_mm": float(positive_steps.median()) if len(positive_steps) else 0.,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", nargs="+", required=True,
                    help="One or more LABEL=PREFIX entries.")
    ap.add_argument("--output-prefix", default="phaseB_convergence")
    args = ap.parse_args()

    runs = [load_run(*parse_run(x)) for x in args.runs]

    # Central-strip current overlay.
    plt.figure(figsize=(8, 5))
    for r in runs:
        plt.plot(r["wave"]["time_ns"] / 1000.,
                 r["wave"]["strip_0_fC_per_ns"],
                 label=r["label"])
    plt.xlabel("Time [µs]")
    plt.ylabel("Current on strip 0 [fC/ns]")
    plt.title("Phase B1a convergence: central-strip waveform")
    plt.legend()
    plt.tight_layout()
    plt.savefig(args.output_prefix + "_waveform.png", dpi=200)
    plt.close()

    # Cumulative central-strip charge.
    plt.figure(figsize=(8, 5))
    for r in runs:
        plt.plot(r["wave"]["time_ns"] / 1000.,
                 r["cumulative"], label=r["label"])
    plt.xlabel("Time [µs]")
    plt.ylabel("Integrated charge on strip 0 [fC]")
    plt.title("Phase B1a convergence: cumulative induced charge")
    plt.legend()
    plt.tight_layout()
    plt.savefig(args.output_prefix + "_cumulative.png", dpi=200)
    plt.close()

    # Actual RKF step size along the trajectory.
    plt.figure(figsize=(8, 5))
    for r in runs:
        d = r["drift"]
        plt.plot(d["time_ns"] / 1000., d["step_mm"], label=r["label"])
    plt.xlabel("Time [µs]")
    plt.ylabel("RKF spatial step [mm]")
    plt.title("Phase B1a convergence: actual drift-line step size")
    plt.legend()
    plt.tight_layout()
    plt.savefig(args.output_prefix + "_rkf_steps.png", dpi=200)
    plt.close()

    rows = []
    for r in runs:
        rows.append({
            "run": r["label"],
            "dt_ns": r["dt_ns"],
            "drift_points": r["n_points"],
            "drift_time_us": r["drift_time_ns"] / 1000.,
            "max_actual_step_mm": r["max_actual_step_mm"],
            "median_actual_step_mm": r["median_actual_step_mm"],
            "Q0_integrated_fC": r["q_int"],
            "Q0_Ramo_fC": r["q_ramo"],
            "Q0_rel_error": r["rel_err"],
            "peak_abs_current_fC_per_ns": r["peak"],
        })

    table = pd.DataFrame(rows)
    table.to_csv(args.output_prefix + "_summary.csv", index=False)
    print(table.to_string(index=False))
    print("\nWrote:")
    print(" ", args.output_prefix + "_waveform.png")
    print(" ", args.output_prefix + "_cumulative.png")
    print(" ", args.output_prefix + "_rkf_steps.png")
    print(" ", args.output_prefix + "_summary.csv")


if __name__ == "__main__":
    main()
