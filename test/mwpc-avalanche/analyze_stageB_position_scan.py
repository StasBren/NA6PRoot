#!/usr/bin/env python3
"""
Controlled Stage-B position-response scan across one strip pitch.

Each input run is one real Garfield microscopic avalanche, generated with the
same random seed but translated in true w relative to the strip pattern.

This deliberately measures RESPONSE LINEARITY / BIAS, not resolution yet.

For each observation window T:
  A_k(T) = | integral_0^T I_k(t) dt |

and for the chosen 3-strip cluster (-1,0,+1 by default),

  w_hat(T) = sum_k w_k A_k / sum_k A_k.

Outputs:
  1) strip_amplitudes_vs_true_position.png
     Shows how the three strip amplitudes exchange as the avalanche moves.

  2) reconstructed_vs_true.png
     Shows whether CoG tracks the true position.

  3) residual_vs_true.png
     Shows deterministic bias/non-linearity w_hat - w_true.

This is the bridge from Stage-B signal formation to spatial reconstruction.
"""

import argparse
import glob
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument(
        "--prefix-glob",
        default="stageB_posscan*_event_summary.csv",
        help="Glob selecting event-summary files from the position scan.",
    )
    p.add_argument(
        "--windows-ns",
        default="25,50,100,200",
        help="Comma-separated early integration windows in ns.",
    )
    p.add_argument(
        "--cluster-strips",
        default="-1,0,1",
        help="Strip indices used in CoG reconstruction.",
    )
    p.add_argument(
        "--representative-window-ns",
        type=float,
        default=100.0,
        help="Window used for the strip-amplitude plot (default: 100 ns).",
    )
    p.add_argument(
        "--output-prefix",
        default="stageB_position_scan",
    )
    return p.parse_args()


def strip_label(i):
    if i < 0:
        return f"strip_m{-i}"
    if i > 0:
        return f"strip_p{i}"
    return "strip_0"


def infer_dt(times):
    if len(times) < 2:
        raise RuntimeError("Need at least two time bins.")
    return float(np.median(np.diff(times)))


def integrate_to(times, current, T):
    times = np.asarray(times, dtype=float)
    current = np.asarray(current, dtype=float)
    dt = infer_dt(times)

    left = times - 0.5 * dt
    right = times + 0.5 * dt
    overlap = np.minimum(right, T) - np.maximum(left, 0.0)
    overlap = np.clip(overlap, 0.0, None)

    return float(np.sum(current * overlap))


def load_total_charge(prefix, strip_ids, T):
    e = pd.read_csv(prefix + "_electron_waveforms.csv")
    ion = pd.read_csv(prefix + "_ion_waveforms.csv")

    te = e["time_ns"].to_numpy(dtype=float)
    ti = ion["time_ns"].to_numpy(dtype=float)

    charges = {}
    for k in strip_ids:
        label = strip_label(int(k))
        ce = e[f"{label}_electron_fC_per_ns"].to_numpy(dtype=float)
        ci = ion[f"{label}_ion_fC_per_ns"].to_numpy(dtype=float)

        qe = integrate_to(te, ce, min(T, te[-1] + 0.5 * infer_dt(te)))
        qi = integrate_to(ti, ci, T)
        charges[int(k)] = qe + qi

    return charges


def main():
    args = parse_args()

    windows = sorted(
        float(x.strip())
        for x in args.windows_ns.split(",")
        if x.strip()
    )
    cluster_strips = [
        int(x.strip())
        for x in args.cluster_strips.split(",")
        if x.strip()
    ]
    if len(cluster_strips) < 2:
        raise ValueError("Need at least two strips for CoG reconstruction.")

    files = sorted(glob.glob(args.prefix_glob))
    if not files:
        raise FileNotFoundError(
            f"No files matched {args.prefix_glob!r}"
        )

    rows = []

    for event_file in files:
        event = pd.read_csv(event_file).iloc[0]

        if "seed_w_mm" not in event.index:
            raise KeyError(
                f"{event_file} lacks seed_w_mm; rerun with updated executable."
            )

        true_w = float(event["seed_w_mm"])
        prefix = event_file[: -len("_event_summary.csv")]

        for T in windows:
            q = load_total_charge(prefix, cluster_strips, T)

            amps = np.array(
                [abs(q[k]) for k in cluster_strips],
                dtype=float,
            )
            centers = np.array(
                [1.7 * k for k in cluster_strips],
                dtype=float,
            )

            total_amp = float(np.sum(amps))
            if total_amp <= 0.0:
                what = np.nan
            else:
                what = float(np.sum(centers * amps) / total_amp)

            row = {
                "prefix": prefix,
                "true_w_mm": true_w,
                "window_ns": T,
                "w_hat_mm": what,
                "residual_mm": what - true_w,
                "total_cluster_amplitude_fC": total_amp,
            }

            for k in cluster_strips:
                row[f"A_strip_{k}_fC"] = abs(q[k])
                row[f"fraction_strip_{k}"] = (
                    abs(q[k]) / total_amp if total_amp > 0 else np.nan
                )

            rows.append(row)

    df = pd.DataFrame(rows).sort_values(
        ["window_ns", "true_w_mm"]
    )
    df.to_csv(args.output_prefix + "_summary.csv", index=False)

    # ------------------------------------------------------------------
    # 1) Strip amplitudes at one representative early window.
    # ------------------------------------------------------------------
    Tref = min(
        windows,
        key=lambda x: abs(x - args.representative_window_ns),
    )
    dref = df[np.isclose(df["window_ns"], Tref)].sort_values("true_w_mm")

    fig, ax = plt.subplots(figsize=(9.5, 5.8))
    for k in cluster_strips:
        ax.plot(
            dref["true_w_mm"],
            dref[f"A_strip_{k}_fC"],
            marker="o",
            label=f"strip {k}",
        )
    ax.set_xlabel("true avalanche position w0 [mm]")
    ax.set_ylabel(f"|Q_k(T={Tref:g} ns)| [fC]")
    ax.set_title(
        "Stage B: three-strip response across one strip pitch"
    )
    ax.legend()
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(
        args.output_prefix + "_strip_amplitudes_vs_true_position.png",
        dpi=200,
    )
    plt.close(fig)

    # ------------------------------------------------------------------
    # 2) Reconstructed coordinate vs true coordinate.
    # ------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(7.5, 6.5))
    xmin = float(df["true_w_mm"].min())
    xmax = float(df["true_w_mm"].max())
    ax.plot([xmin, xmax], [xmin, xmax], linestyle="--", label="ideal w_hat = w0")

    for T in windows:
        d = df[np.isclose(df["window_ns"], T)].sort_values("true_w_mm")
        ax.plot(
            d["true_w_mm"],
            d["w_hat_mm"],
            marker="o",
            label=f"{T:g} ns",
        )

    ax.set_xlabel("true position w0 [mm]")
    ax.set_ylabel("CoG reconstructed position w_hat [mm]")
    ax.set_title(
        "Stage B: does three-strip charge sharing encode sub-strip position?"
    )
    ax.legend()
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(
        args.output_prefix + "_reconstructed_vs_true.png",
        dpi=200,
    )
    plt.close(fig)

    # ------------------------------------------------------------------
    # 3) Deterministic bias / nonlinearity.
    # ------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(9.5, 5.8))
    ax.axhline(0.0, linewidth=0.8)

    for T in windows:
        d = df[np.isclose(df["window_ns"], T)].sort_values("true_w_mm")
        ax.plot(
            d["true_w_mm"],
            1000.0 * d["residual_mm"],
            marker="o",
            label=f"{T:g} ns",
        )

    ax.set_xlabel("true position w0 [mm]")
    ax.set_ylabel("w_hat - w0 [um]")
    ax.set_title(
        "Stage B: CoG reconstruction bias across one strip pitch"
    )
    ax.legend()
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(
        args.output_prefix + "_residual_vs_true.png",
        dpi=200,
    )
    plt.close(fig)

    print("\n=== STAGE B: CONTROLLED POSITION RESPONSE ===")
    print(f"runs found              : {len(files)}")
    print(f"true-w range            : {xmin:.3f} .. {xmax:.3f} mm")
    print(f"CoG strips              : {cluster_strips}")
    print(f"observation windows     : {windows} ns")
    print(f"representative window   : {Tref:g} ns")

    print("\nMaximum absolute deterministic bias:")
    for T in windows:
        d = df[np.isclose(df["window_ns"], T)]
        max_bias_um = 1000.0 * float(np.nanmax(np.abs(d["residual_mm"])))
        print(f"  {T:7g} ns : {max_bias_um:9.3f} um")

    print("\nWrote:")
    print(" ", args.output_prefix + "_summary.csv")
    print(" ", args.output_prefix + "_strip_amplitudes_vs_true_position.png")
    print(" ", args.output_prefix + "_reconstructed_vs_true.png")
    print(" ", args.output_prefix + "_residual_vs_true.png")


if __name__ == "__main__":
    main()
