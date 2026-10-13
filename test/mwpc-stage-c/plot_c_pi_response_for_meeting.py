#!/usr/bin/env python3
"""Make presentation-ready diagnostic figures for the PI MWPC response study."""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def save(fig, path):
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--u-bins", required=True)
    ap.add_argument("--transport-summary", required=True)
    ap.add_argument("--heed-tracks", required=True)
    ap.add_argument("--output-dir", required=True, type=Path)
    args = ap.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    u = pd.read_csv(args.u_bins)
    tr = pd.read_csv(args.transport_summary)
    ht = pd.read_csv(args.heed_tracks)

    # 1) Heed cluster multiplicity vs Poisson with same mean.
    counts = ht["clusters"].to_numpy(dtype=int)
    mu = float(np.mean(counts))
    kmin = max(0, int(np.min(counts)))
    kmax = int(np.percentile(counts, 99.8))
    ks = np.arange(kmin, kmax + 1)
    obs = np.array([(counts == k).mean() for k in ks])
    # Stable Poisson recursion.
    p = math.exp(-mu)
    probs = {0: p}
    for k in range(1, max(kmax, 1) + 1):
        p *= mu / k
        probs[k] = p
    pois = np.array([probs.get(int(k), 0.0) for k in ks])

    fig, ax = plt.subplots(figsize=(8.6, 5.5))
    ax.step(ks, obs, where="mid", linewidth=1.8, label="Heed")
    ax.plot(ks, pois, marker="o", markersize=3, linewidth=1.2,
            label=f"Poisson, mean={mu:.2f}")
    ax.set_xlabel("Ionisation clusters per muon track")
    ax.set_ylabel("Probability")
    ax.set_title("Heed cluster multiplicity is Poisson-like")
    ax.grid(alpha=0.2)
    ax.legend()
    save(fig, args.output_dir / "pi_01_heed_poisson.png")

    # 2) 2D response map: mean drift time.
    pivot = tr.pivot(index="v_center_mm", columns="u_center_mm",
                     values="mean_drift_time_ns").sort_index()
    fig, ax = plt.subplots(figsize=(8.8, 6.0))
    extent = [
        pivot.columns.min() - 0.1, pivot.columns.max() + 0.1,
        pivot.index.min() - 0.1, pivot.index.max() + 0.1,
    ]
    im = ax.imshow(
        pivot.to_numpy(),
        origin="lower",
        aspect="auto",
        extent=extent,
        interpolation="nearest",
    )
    cb = fig.colorbar(im, ax=ax)
    cb.set_label("Mean drift time [ns]")
    ax.set_xlabel("u across wires [mm]")
    ax.set_ylabel("v chamber normal [mm]")
    ax.set_title("Garfield transport map: mean electron drift time")
    save(fig, args.output_dir / "pi_02_transport_drift_time_map.png")

    # 3) Arrival-side observables vs true sub-wire coordinate.
    fig, ax = plt.subplots(figsize=(8.6, 5.5))
    ax.plot(u["u_center_mm"], u["mean_cos_phi"], marker="o",
            label=r"$\langle\cos\phi_{hit}\rangle$")
    ax.plot(u["u_center_mm"], u["mean_A_LR"], marker="s",
            label=r"$A_{LR}$")
    ax.axhline(0.0, linewidth=1.0)
    ax.axvline(0.0, linewidth=1.0, linestyle="--")
    ax.set_xlabel("True sub-wire coordinate $u_0$ [mm]")
    ax.set_ylabel("Arrival-side observable")
    ax.set_title("Arrival direction retains the sign of the sub-wire position")
    ax.set_ylim(-1.05, 1.05)
    ax.grid(alpha=0.2)
    ax.legend()
    save(fig, args.output_dir / "pi_03_arrival_asymmetry_vs_u.png")

    # 4) Timing windows vs coordinate.
    fig, ax = plt.subplots(figsize=(8.6, 5.5))
    for t in (25, 50, 100):
        ax.plot(
            u["u_center_mm"],
            u[f"mean_fraction_before_{t}ns"],
            marker="o",
            label=f"arrived by {t} ns",
        )
    ax.set_xlabel("True sub-wire coordinate $u_0$ [mm]")
    ax.set_ylabel("Fraction of collected electrons")
    ax.set_title("Finite time window carries information on distance from wire")
    ax.set_ylim(-0.02, 1.02)
    ax.grid(alpha=0.2)
    ax.legend()
    save(fig, args.output_dir / "pi_04_timing_windows_vs_u.png")

    # 5) Multi-wire probability.
    fig, ax = plt.subplots(figsize=(8.6, 5.5))
    ax.plot(
        u["u_center_mm"],
        100.0 * u["multiwire_fraction"],
        marker="o",
    )
    ax.set_xlabel("True sub-wire coordinate $u_0$ [mm]")
    ax.set_ylabel("Events with >1 active wire [%]")
    ax.set_title("Multi-wire response identifies the cell boundary")
    ax.set_ylim(-2, 102)
    ax.grid(alpha=0.2)
    save(fig, args.output_dir / "pi_05_multiwire_fraction_vs_u.png")

    # 6) Mean drift time vs coordinate.
    fig, ax = plt.subplots(figsize=(8.6, 5.5))
    ax.plot(
        u["u_center_mm"],
        u["mean_drift_time_ns"],
        marker="o",
    )
    ax.set_xlabel("True sub-wire coordinate $u_0$ [mm]")
    ax.set_ylabel("Mean drift time [ns]")
    ax.set_title("Mean electron drift time measures distance from the wire")
    ax.grid(alpha=0.2)
    save(fig, args.output_dir / "pi_06_mean_drift_time_vs_u.png")

    print("Wrote:")
    for p in sorted(args.output_dir.glob("pi_*.png")):
        print(f"  {p}")


if __name__ == "__main__":
    main()
