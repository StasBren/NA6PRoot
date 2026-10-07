#!/usr/bin/env python3
"""Synthesis plot for Stage B: microscopic current -> strip charge -> CoG.

Uses one already-generated real Garfield avalanche.  The upper panels show the
prompt electron burst and slow ion current on their natural time scales.  The
lower panels integrate the same waveforms and show how the strip amplitudes and
the reconstructed CoG position emerge as the observation window grows.

No electronics shaping, threshold, ADC, or noise is applied.
"""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument(
        "--prefix",
        default="stageB_posscan_p0p40",
        help="Prefix produced by mwpc_phase_b_real_avalanche.",
    )
    p.add_argument(
        "--cluster-strips",
        default="-2,-1,0,1,2",
        help="Strip indices used for CoG in the lower-right panel.",
    )
    p.add_argument(
        "--max-window-ns",
        type=float,
        default=5000.0,
        help="Largest integration window shown.",
    )
    p.add_argument(
        "--cog-min-fraction",
        type=float,
        default=0.01,
        help=(
            "Hide CoG points until the accumulated cluster amplitude exceeds "
            "this fraction of its final value. This avoids showing a centroid "
            "with a nearly zero denominator before the avalanche signal arrives."
        ),
    )
    p.add_argument(
        "--output",
        default="stageB_response_kernel.png",
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


def main():
    args = parse_args()
    prefix = Path(args.prefix)

    e = pd.read_csv(str(prefix) + "_electron_waveforms.csv")
    ion = pd.read_csv(str(prefix) + "_ion_waveforms.csv")
    strip_summary = pd.read_csv(str(prefix) + "_strip_summary.csv")
    event_summary = pd.read_csv(str(prefix) + "_event_summary.csv").iloc[0]

    cluster = [
        int(x.strip())
        for x in args.cluster_strips.split(",")
        if x.strip()
    ]
    if not cluster:
        raise ValueError("--cluster-strips cannot be empty.")

    centers_map = {
        int(row["strip"]): float(row["center_w_mm"])
        for _, row in strip_summary.iterrows()
    }
    missing_centers = [k for k in cluster if k not in centers_map]
    if missing_centers:
        raise KeyError(f"Missing strip centers for {missing_centers}")

    te = e["time_ns"].to_numpy(dtype=float)
    ti = ion["time_ns"].to_numpy(dtype=float)
    e_full = float(te[-1] + 0.5 * infer_dt(te))
    i_full = float(ti[-1] + 0.5 * infer_dt(ti))
    tmax = min(float(args.max_window_ns), i_full)
    if tmax <= 0:
        raise RuntimeError("No positive integration range available.")

    key = [-1, 0, 1]
    ecur = {}
    icur = {}
    for k in set(key + cluster):
        label = strip_label(k)
        ecol = f"{label}_electron_fC_per_ns"
        icol = f"{label}_ion_fC_per_ns"
        if ecol not in e.columns or icol not in ion.columns:
            raise KeyError(f"Missing waveform columns for strip {k}")
        ecur[k] = e[ecol].to_numpy(dtype=float)
        icur[k] = ion[icol].to_numpy(dtype=float)

    peak = max(float(np.max(np.abs(ecur[k]))) for k in key)
    if peak > 0:
        combined = np.max(np.vstack([np.abs(ecur[k]) for k in key]), axis=0)
        active = np.flatnonzero(combined > 0.01 * peak)
        if len(active):
            i0 = max(0, int(active[0]) - 10)
            i1 = min(len(te) - 1, int(active[-1]) + 10)
            electron_xlim = (te[i0], te[i1])
        else:
            electron_xlim = (te[0], te[-1])
    else:
        electron_xlim = (te[0], te[-1])

    tmin = max(0.2, min(1.0, tmax))
    windows = np.unique(
        np.concatenate(
            [
                np.geomspace(tmin, tmax, 180),
                np.array([25., 50., 100., 200., 500., 1000.]),
            ]
        )
    )
    windows = windows[(windows > 0.) & (windows <= tmax)]

    amp = {k: [] for k in cluster}
    xhat = []
    cluster_sum = []

    for T in windows:
        amps = []
        centers = []
        for k in cluster:
            qe = integrate_to(te, ecur[k], min(T, e_full))
            qi = integrate_to(ti, icur[k], T)
            a = abs(qe + qi)
            amp[k].append(a)
            amps.append(a)
            centers.append(centers_map[k])

        amps = np.asarray(amps, dtype=float)
        centers = np.asarray(centers, dtype=float)
        denom = float(np.sum(amps))
        cluster_sum.append(denom)
        xhat.append(
            float(np.sum(centers * amps) / denom)
            if denom > 0 else np.nan
        )

    cluster_sum = np.asarray(cluster_sum, dtype=float)
    xhat = np.asarray(xhat, dtype=float)
    final_cluster_sum = float(np.nanmax(cluster_sum))
    valid_cog = cluster_sum >= args.cog_min_fraction * final_cluster_sum
    xhat_visible = np.where(valid_cog, xhat, np.nan)

    fig, axes = plt.subplots(2, 2, figsize=(12.2, 8.4))
    ax_e, ax_i, ax_q, ax_x = axes.flat

    styles = {-1: "--", 0: "-", 1: ":"}
    for k in key:
        ax_e.plot(te, ecur[k], linestyle=styles[k], label=f"strip {k}")
    ax_e.set_xlim(*electron_xlim)
    ax_e.set_xlabel("time [ns]")
    ax_e.set_ylabel("electron current [fC/ns]")
    ax_e.set_title("1. Prompt avalanche-electron current")
    ax_e.grid(alpha=0.2)
    ax_e.legend()

    for k in key:
        ax_i.plot(ti / 1000.0, icur[k], linestyle=styles[k], label=f"strip {k}")
    ax_i.set_xlim(0.0, min(tmax, ti[-1]) / 1000.0)
    ax_i.set_xlabel(r"time [$\mu$s]")
    ax_i.set_ylabel("positive-ion current [fC/ns]")
    ax_i.set_title("2. Slow positive-ion response")
    ax_i.grid(alpha=0.2)
    ax_i.legend()

    for k in key:
        if k in amp:
            ax_q.semilogx(windows, amp[k], linestyle=styles[k], label=f"strip {k}")
    ax_q.set_xlabel("observation window $T$ [ns]")
    ax_q.set_ylabel(r"$A_k(T)=|Q_k(T)|$ [fC]")
    ax_q.set_title("3. Integrated strip amplitudes")
    ax_q.grid(alpha=0.2)
    ax_q.legend()

    true_w = float(event_summary.get("seed_w_mm", np.nan))
    ax_x.semilogx(windows, 1000.0 * xhat_visible, linewidth=2.2, label="5-strip CoG")
    if np.isfinite(true_w):
        ax_x.axhline(
            1000.0 * true_w,
            color="black",
            linestyle="--",
            label=rf"true $w_0={true_w:g}$ mm",
        )
    if np.any(valid_cog):
        conv = 1000.0 * xhat_visible[np.where(valid_cog)[0][-1]]
        ax_x.axhline(
            conv,
            color="0.35",
            linestyle=":",
            label=rf"late CoG $\approx {conv:.0f}\,\mu$m",
        )
    ax_x.set_xlabel("observation window $T$ [ns]")
    ax_x.set_ylabel(r"cluster-CoG $\hat w(T)$ [$\mu$m]")
    ax_x.set_title("4. Position stabilizes before full ion collection")
    ax_x.grid(alpha=0.2)
    ax_x.legend()

    fig.suptitle("Stage B response kernel: microscopic motion -> induced charge -> position")
    fig.tight_layout()
    fig.savefig(args.output, dpi=220)
    plt.close(fig)

    anchors = [25., 50., 100., 200.]
    print("\n=== STAGE B: RESPONSE KERNEL ===")
    print(f"prefix                   : {args.prefix}")
    print(f"CoG strips               : {cluster}")
    print(f"CoG visibility threshold : {100.0 * args.cog_min_fraction:.2f}% of max cluster amplitude")
    if np.isfinite(true_w):
        print(f"true position            : {true_w:.4f} mm")
    print("\nCoG versus observation window:")
    for T in anchors:
        if T > tmax:
            continue
        j = int(np.argmin(np.abs(windows - T)))
        print(
            f"  T={windows[j]:7.2f} ns : "
            f"w_hat={xhat[j]: .6f} mm"
            + (
                f", residual={1000.0*(xhat[j]-true_w):+.2f} um"
                if np.isfinite(true_w) else ""
            )
        )
    print("\nWrote:")
    print(" ", args.output)


if __name__ == "__main__":
    main()
