#!/usr/bin/env python3
"""Plot the completed brute-force Stage-C wire-boundary validation scan."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--points", required=True)
    ap.add_argument("--transitions", required=True)
    ap.add_argument("--output-dir", required=True, type=Path)
    args = ap.parse_args()

    p = pd.read_csv(args.points)
    t = pd.read_csv(args.transitions)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Fold the two sides to |u| and use the probability that the outer
    # neighbour wire dominates.  This makes left/right symmetry visible.
    neg = p[p["u0_mm"] < 0].copy()
    neg["abs_u_mm"] = -neg["u0_mm"]
    neg["p_outer"] = neg["dominant_wire_minus1_fraction"]

    pos = p[p["u0_mm"] > 0].copy()
    pos["abs_u_mm"] = pos["u0_mm"]
    pos["p_outer"] = pos["dominant_wire_plus1_fraction"]

    fig, ax = plt.subplots(figsize=(8.6, 5.5))
    ax.plot(
        neg["abs_u_mm"], neg["p_outer"],
        marker="o", label="negative boundary",
    )
    ax.plot(
        pos["abs_u_mm"], pos["p_outer"],
        marker="s", label="positive boundary",
    )
    ax.axvline(2.0, linestyle="--", linewidth=1.2,
               label="geometric boundary")
    ax.axhline(0.5, linestyle=":", linewidth=1.0)
    ax.set_xlabel(r"$|u_0|$ [mm]")
    ax.set_ylabel("Probability outer neighbour is dominant")
    ax.set_title("Brute-force Heed + avalanche wire-boundary transition")
    ax.set_ylim(-0.03, 1.03)
    ax.grid(alpha=0.2)
    ax.legend()
    fig.tight_layout()
    fig.savefig(args.output_dir / "pi_bruteforce_boundary_transition.png", dpi=180)
    plt.close(fig)

    neg["multi"] = neg["multiwire_fraction"]
    pos["multi"] = pos["multiwire_fraction"]
    fig, ax = plt.subplots(figsize=(8.6, 5.5))
    ax.plot(neg["abs_u_mm"], 100 * neg["multi"], marker="o",
            label="negative boundary")
    ax.plot(pos["abs_u_mm"], 100 * pos["multi"], marker="s",
            label="positive boundary")
    ax.axvline(2.0, linestyle="--", linewidth=1.2,
               label="geometric boundary")
    ax.set_xlabel(r"$|u_0|$ [mm]")
    ax.set_ylabel("Events with >1 active wire [%]")
    ax.set_title("Brute-force avalanche sharing is localised near the boundary")
    ax.set_ylim(-3, 103)
    ax.grid(alpha=0.2)
    ax.legend()
    fig.tight_layout()
    fig.savefig(args.output_dir / "pi_bruteforce_multiwire_transition.png", dpi=180)
    plt.close(fig)

    print("Transition summary:")
    print(t.to_string(index=False))
    print("Wrote:")
    print(args.output_dir / "pi_bruteforce_boundary_transition.png")
    print(args.output_dir / "pi_bruteforce_multiwire_transition.png")


if __name__ == "__main__":
    main()
