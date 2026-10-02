#!/usr/bin/env python3
"""
Validate the fast truncated-ion transport used by phaseB_resolution_ensemble
against the full DriftLineRKF waveform from phaseB_real_avalanche_signal.

Both runs must use the same microscopic avalanche:
  same w0, same Garfield random seed, one attempt.

The script compares finite-window Q_k(T) on strips -1, 0, +1.
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument(
        "--ensemble",
        default="stageB_resolution_validation.csv",
    )
    p.add_argument(
        "--full-prefix",
        default="stageB_resolution_validation_full",
    )
    return p.parse_args()


def infer_dt(t):
    return float(np.median(np.diff(t)))


def integrate_to(t, y, T):
    t = np.asarray(t, dtype=float)
    y = np.asarray(y, dtype=float)
    dt = infer_dt(t)
    left = t - 0.5 * dt
    right = t + 0.5 * dt
    overlap = np.minimum(right, T) - np.maximum(left, 0.0)
    overlap = np.clip(overlap, 0.0, None)
    return float(np.sum(y * overlap))


def label(k):
    if k < 0:
        return f"strip_m{-k}"
    if k > 0:
        return f"strip_p{k}"
    return "strip_0"


def main():
    args = parse_args()

    ensemble_path = Path(args.ensemble)
    electron_path = Path(args.full_prefix + "_electron_waveforms.csv")
    ion_path = Path(args.full_prefix + "_ion_waveforms.csv")

    missing = [
        str(p)
        for p in (ensemble_path, electron_path, ion_path)
        if not p.exists()
    ]
    if missing:
        raise FileNotFoundError(
            "Validation prerequisites are missing: "
            + ", ".join(missing)
            + ". Run both the full-avalanche and truncated-ensemble "
              "validation commands before this script."
        )

    fast = pd.read_csv(ensemble_path)
    e = pd.read_csv(electron_path)
    i = pd.read_csv(ion_path)

    te = e["time_ns"].to_numpy(dtype=float)
    ti = i["time_ns"].to_numpy(dtype=float)
    e_end = te[-1] + 0.5 * infer_dt(te)

    rows = []
    for _, r in fast.iterrows():
        T = float(r["window_ns"])

        for k in (-1, 0, 1):
            lab = label(k)
            qe = integrate_to(
                te,
                e[f"{lab}_electron_fC_per_ns"].to_numpy(dtype=float),
                min(T, e_end),
            )
            qi = integrate_to(
                ti,
                i[f"{lab}_ion_fC_per_ns"].to_numpy(dtype=float),
                T,
            )
            q_full = qe + qi
            q_fast = float(r[f"Q_strip_{k}_fC"])

            scale = max(abs(q_full), 1e-30)
            rows.append(
                {
                    "window_ns": T,
                    "strip": k,
                    "Q_full_fC": q_full,
                    "Q_fast_fC": q_fast,
                    "abs_diff_fC": q_fast - q_full,
                    "relative_diff": (q_fast - q_full) / scale,
                }
            )

    out = pd.DataFrame(rows)
    print("\n=== FAST TRUNCATED ION TRANSPORT VALIDATION ===")
    print(out.to_string(index=False))
    print(
        "\nmax |relative difference| = "
        f"{np.max(np.abs(out['relative_diff'])):.6g}"
    )


if __name__ == "__main__":
    main()
