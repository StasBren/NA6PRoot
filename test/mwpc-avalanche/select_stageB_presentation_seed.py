#!/usr/bin/env python3
"""
Select a reproducible, unconditioned example avalanche for presentation.

The input is the Stage-B ensemble CSV.  At one chosen true position (default
w0 = 0), we collapse repeated rows from the different observation windows to
one row per random seed, keep successful non-zero avalanches, and select the
event whose Garfield ion count is closest to the median positive gain.

This avoids cherry-picking an unusually large avalanche while still giving a
clean, reproducible example for the signal-formation plots.
"""

import argparse
import sys

import numpy as np
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument(
        "--input",
        default="stageB_resolution_ensemble.csv",
        help="Ensemble CSV produced by mwpc_phase_b_resolution_ensemble.",
    )
    p.add_argument(
        "--true-w-mm",
        type=float,
        default=0.0,
        help="True position from which to choose the example avalanche.",
    )
    p.add_argument(
        "--seed-only",
        action="store_true",
        help="Print only the selected random seed to stdout.",
    )
    return p.parse_args()


def main():
    args = parse_args()
    df = pd.read_csv(args.input)

    required = {
        "true_w_mm",
        "random_seed",
        "avalanche_ok",
        "garfield_electrons",
        "garfield_ions",
        "recorded_ion_births",
    }
    missing = required - set(df.columns)
    if missing:
        raise KeyError(f"Missing columns: {sorted(missing)}")

    d = df[np.isclose(df["true_w_mm"], args.true_w_mm)].copy()
    if d.empty:
        raise RuntimeError(
            f"No rows found at true_w_mm={args.true_w_mm:g}."
        )

    # The same event appears once per observation window.
    d = (
        d.sort_values("random_seed")
        .drop_duplicates(subset=["random_seed"], keep="first")
    )

    good = d[
        (d["avalanche_ok"] == 1)
        & (d["garfield_ions"] > 0)
    ].copy()

    if good.empty:
        raise RuntimeError(
            "No successful non-zero avalanches were found at the requested "
            "position."
        )

    median_ions = float(np.median(good["garfield_ions"]))
    good["distance_to_median"] = np.abs(
        good["garfield_ions"] - median_ions
    )
    selected = good.sort_values(
        ["distance_to_median", "random_seed"]
    ).iloc[0]

    seed = int(selected["random_seed"])
    ne = int(selected["garfield_electrons"])
    ni = int(selected["garfield_ions"])
    nb = int(selected["recorded_ion_births"])

    if args.seed_only:
        print(seed)
        return

    print("=== STAGE B: PRESENTATION BASELINE SELECTION ===")
    print(f"true position w0        : {args.true_w_mm:g} mm")
    print(f"positive events found  : {len(good)}")
    print(f"median Garfield ions   : {median_ions:.1f}")
    print(f"selected random seed   : {seed}")
    print(f"selected e / ions      : {ne} / {ni}")
    print(f"recorded ion births    : {nb}")
    print()
    print(
        "Selection rule: successful non-zero avalanche closest to the "
        "median positive gain. No minimum-gain conditioning is applied."
    )


if __name__ == "__main__":
    main()
