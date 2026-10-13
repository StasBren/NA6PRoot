#!/usr/bin/env python3
"""Plot representative high-stat drift-time distributions with Gaussian overlays."""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


TARGETS = [
    (0.1, 0.1, "central_near_anode"),
    (1.5, 0.1, "approaching_boundary"),
    (1.9, 0.1, "boundary_near_anode"),
    (1.9, 3.9, "boundary_near_plus_cathode"),
]


def nearest_box(df, u0, v0):
    boxes = (
        df[["iu", "iv", "u_center_mm", "v_center_mm"]]
        .drop_duplicates()
        .copy()
    )
    d2 = (boxes["u_center_mm"] - u0) ** 2 + (boxes["v_center_mm"] - v0) ** 2
    row = boxes.loc[d2.idxmin()]
    return int(row["iu"]), int(row["iv"]), float(row["u_center_mm"]), float(row["v_center_mm"])


def gaussian_pdf(x, mu, sigma):
    return np.exp(-0.5 * ((x - mu) / sigma) ** 2) / (sigma * np.sqrt(2 * np.pi))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--output-dir", required=True, type=Path)
    ap.add_argument("--bins", type=int, default=40)
    args = ap.parse_args()

    df = pd.read_csv(args.input)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    for u0, v0, tag in TARGETS:
        iu, iv, uc, vc = nearest_box(df, u0, v0)
        g = df[(df["iu"] == iu) & (df["iv"] == iv) & (df["collected"] == 1)].copy()
        if len(g) < 20:
            print(f"skip {tag}: only {len(g)} collected electrons")
            continue

        # Unconditional collected-electron timing.
        t = g["drift_time_ns"].to_numpy(dtype=float)
        t = t[np.isfinite(t)]
        mu = float(np.mean(t))
        sigma = float(np.std(t, ddof=1))

        fig, ax = plt.subplots(figsize=(8.6, 5.5))
        values, edges, _ = ax.hist(
            t,
            bins=args.bins,
            density=True,
            alpha=0.65,
            label="Garfield electrons",
        )
        if sigma > 0:
            xx = np.linspace(edges[0], edges[-1], 400)
            ax.plot(
                xx,
                gaussian_pdf(xx, mu, sigma),
                linewidth=2.0,
                label=f"Gaussian with same mean/std",
            )
        ax.set_xlabel("Drift time [ns]")
        ax.set_ylabel("Probability density")
        ax.set_title(f"Drift-time distribution at u={uc:.1f} mm, v={vc:.1f} mm")
        ax.grid(alpha=0.2)
        ax.legend()
        fig.tight_layout()
        out = args.output_dir / f"pi_timing_{tag}_all_wires.png"
        fig.savefig(out, dpi=180)
        plt.close(fig)

        # Near a separatrix, show destination-wire-conditioned distributions.
        wires = sorted(g["wire_index"].dropna().astype(int).unique())
        if len(wires) > 1:
            fig, ax = plt.subplots(figsize=(8.6, 5.5))
            for wire in wires:
                tw = g[g["wire_index"] == wire]["drift_time_ns"].to_numpy(dtype=float)
                tw = tw[np.isfinite(tw)]
                if len(tw) < 10:
                    continue
                ax.hist(
                    tw,
                    bins=args.bins,
                    density=True,
                    histtype="step",
                    linewidth=1.8,
                    label=f"wire {wire}",
                )
            ax.set_xlabel("Drift time [ns]")
            ax.set_ylabel("Probability density")
            ax.set_title(
                f"Timing conditioned on destination wire, u={uc:.1f} mm, v={vc:.1f} mm"
            )
            ax.grid(alpha=0.2)
            ax.legend()
            fig.tight_layout()
            out = args.output_dir / f"pi_timing_{tag}_by_wire.png"
            fig.savefig(out, dpi=180)
            plt.close(fig)

    print("Wrote:")
    for p in sorted(args.output_dir.glob("pi_timing_*.png")):
        print(f"  {p}")


if __name__ == "__main__":
    main()
