#!/usr/bin/env python3
"""
Stage-B time-window analysis for the synthetic avalanche-ion cloud.

This script stays strictly at detector-signal level. It takes the strip current
waveforms produced by mwpc_phase_b_ion_cloud and derives:

  1) peak absolute induced current on each strip,
       I_peak,k = max_t |I_k(t)|

  2) finite-window induced charge,
       Q_k(T) = integral_0^T I_k(t) dt

  3) a positive cluster-width diagnostic based on |Q_k(T)|,
       w_k(T) = |Q_k(T)| / sum_j |Q_j(T)|

  4) the corresponding |Q|-weighted centroid and RMS width vs integration time.

No frontend electronics, shaping, threshold, ADC, or VMM model is applied here.
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument(
        "--prefix",
        default="stageB_cloud_2mm",
        help="Input prefix from mwpc_phase_b_ion_cloud",
    )
    p.add_argument(
        "--windows-ns",
        default="50,100,200,500,1000,5000",
        help=(
            "Comma-separated finite integration windows in ns. "
            "The full available waveform is added automatically."
        ),
    )
    p.add_argument(
        "--output-prefix",
        default=None,
        help="Default: <prefix>_time_windows",
    )
    return p.parse_args()


def strip_label(i):
    if i < 0:
        return f"strip_m{-i}"
    if i > 0:
        return f"strip_p{i}"
    return "strip_0"


def weighted_centroid_rms(x, positive_weights):
    x = np.asarray(x, dtype=float)
    w = np.asarray(positive_weights, dtype=float)
    total = np.sum(w)
    if total <= 0:
        return np.nan, np.nan
    w = w / total
    centroid = float(np.sum(w * x))
    rms = float(np.sqrt(np.sum(w * (x - centroid) ** 2)))
    return centroid, rms


def cumulative_charge_at_window(times_ns, current, window_ns):
    """
    Integrate piecewise-constant current over [0, window_ns].

    The waveform samples are bin centers. We infer a uniform bin width and
    handle a partially covered last bin exactly under the piecewise-constant
    assumption.
    """
    times_ns = np.asarray(times_ns, dtype=float)
    current = np.asarray(current, dtype=float)

    if len(times_ns) < 2:
        raise RuntimeError("Need at least two waveform bins.")

    dt = float(np.median(np.diff(times_ns)))
    left_edges = times_ns - 0.5 * dt
    right_edges = times_ns + 0.5 * dt

    overlap = np.minimum(right_edges, window_ns) - np.maximum(left_edges, 0.0)
    overlap = np.clip(overlap, 0.0, None)

    return float(np.sum(current * overlap))


def main():
    args = parse_args()

    prefix = Path(args.prefix)
    wave_file = Path(str(prefix) + "_waveforms.csv")
    summary_file = Path(str(prefix) + "_summary.csv")

    wave = pd.read_csv(wave_file)
    summary = pd.read_csv(summary_file).sort_values("center_w_mm")

    times_ns = wave["time_ns"].to_numpy(dtype=float)
    if len(times_ns) < 2:
        raise RuntimeError("Waveform file has too few time bins.")

    dt_ns = float(np.median(np.diff(times_ns)))
    full_window_ns = float(times_ns[-1] + 0.5 * dt_ns)

    finite_windows = sorted(
        {
            float(x.strip())
            for x in args.windows_ns.split(",")
            if x.strip()
        }
    )
    if any(x <= 0 for x in finite_windows):
        raise ValueError("All integration windows must be positive.")

    windows = [x for x in finite_windows if x < full_window_ns]
    windows.append(full_window_ns)

    strip_ids = summary["strip"].to_numpy(dtype=int)
    centers = summary["center_w_mm"].to_numpy(dtype=float)

    current_matrix = []
    peak_abs = []
    peak_time = []

    for i in strip_ids:
        col = strip_label(int(i)) + "_fC_per_ns"
        if col not in wave.columns:
            raise KeyError(f"Missing waveform column {col!r}")
        y = wave[col].to_numpy(dtype=float)
        current_matrix.append(y)
        j = int(np.argmax(np.abs(y)))
        peak_abs.append(abs(float(y[j])))
        peak_time.append(float(times_ns[j]))

    current_matrix = np.asarray(current_matrix)
    peak_abs = np.asarray(peak_abs, dtype=float)
    peak_time = np.asarray(peak_time, dtype=float)

    peak_centroid, peak_rms = weighted_centroid_rms(centers, peak_abs)
    peak_frac = (
        peak_abs / np.sum(peak_abs)
        if np.sum(peak_abs) > 0
        else np.zeros_like(peak_abs)
    )

    rows = []
    window_profiles = {}

    for T in windows:
        q = np.array(
            [
                cumulative_charge_at_window(times_ns, current_matrix[j], T)
                for j in range(len(strip_ids))
            ],
            dtype=float,
        )
        absq = np.abs(q)
        total_absq = float(np.sum(absq))
        frac = absq / total_absq if total_absq > 0 else np.zeros_like(absq)
        centroid, rms = weighted_centroid_rms(centers, absq)

        window_profiles[T] = {
            "q": q,
            "abs_fraction": frac,
            "centroid_mm": centroid,
            "rms_mm": rms,
        }

        for j, i in enumerate(strip_ids):
            rows.append(
                {
                    "window_ns": T,
                    "is_full_window": bool(np.isclose(T, full_window_ns)),
                    "strip": int(i),
                    "center_w_mm": float(centers[j]),
                    "Q_fC": float(q[j]),
                    "abs_Q_fraction": float(frac[j]),
                    "peak_abs_current_fC_per_ns": float(peak_abs[j]),
                    "peak_time_ns": float(peak_time[j]),
                    "window_centroid_mm": centroid,
                    "window_rms_mm": rms,
                    "peak_profile_centroid_mm": peak_centroid,
                    "peak_profile_rms_mm": peak_rms,
                }
            )

    out = (
        Path(args.output_prefix)
        if args.output_prefix
        else Path(str(prefix) + "_time_windows")
    )

    table = pd.DataFrame(rows)
    table.to_csv(str(out) + "_by_strip.csv", index=False)

    # ------------------------------------------------------------
    # 1) Peak-current spatial profile.
    # ------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(9, 5.5))
    ax.bar(centers, peak_abs)
    ax.set_xlabel("strip center w [mm]")
    ax.set_ylabel("max |I_k(t)| [fC/ns]")
    ax.set_title(
        "Stage B: peak induced-current profile across strips\n"
        f"centroid = {peak_centroid:.3f} mm, RMS = {peak_rms:.3f} mm"
    )
    fig.tight_layout()
    fig.savefig(str(out) + "_peak_current_profile.png", dpi=200)
    plt.close(fig)

    # ------------------------------------------------------------
    # 2) Signed finite-window charge profiles.
    # ------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(10, 6))
    for T in windows:
        prof = window_profiles[T]
        label = (
            "full"
            if np.isclose(T, full_window_ns)
            else (f"{T / 1000.0:g} us" if T >= 1000 else f"{T:g} ns")
        )
        ax.plot(
            centers,
            prof["q"],
            marker="o",
            label=label,
        )
    ax.axhline(0.0, linewidth=0.8)
    ax.set_xlabel("strip center w [mm]")
    ax.set_ylabel("Q_k(T) [fC]")
    ax.set_title("Stage B: signed induced charge vs integration window")
    ax.legend()
    fig.tight_layout()
    fig.savefig(str(out) + "_signed_Q_profiles.png", dpi=200)
    plt.close(fig)

    # ------------------------------------------------------------
    # 3) Positive |Q|-sharing profiles.
    # ------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(10, 6))
    for T in windows:
        prof = window_profiles[T]
        label = (
            "full"
            if np.isclose(T, full_window_ns)
            else (f"{T / 1000.0:g} us" if T >= 1000 else f"{T:g} ns")
        )
        ax.plot(
            centers,
            100.0 * prof["abs_fraction"],
            marker="o",
            label=label,
        )
    ax.set_xlabel("strip center w [mm]")
    ax.set_ylabel("|Q_k(T)| / sum |Q_j(T)| [%]")
    ax.set_title("Stage B: time evolution of charge sharing")
    ax.legend()
    fig.tight_layout()
    fig.savefig(str(out) + "_abs_Q_sharing.png", dpi=200)
    plt.close(fig)

    # ------------------------------------------------------------
    # 4) Cluster RMS vs observation time.
    # ------------------------------------------------------------
    window_times = np.array(windows, dtype=float)
    rms_values = np.array(
        [window_profiles[T]["rms_mm"] for T in windows],
        dtype=float,
    )
    centroid_values = np.array(
        [window_profiles[T]["centroid_mm"] for T in windows],
        dtype=float,
    )

    fig, ax = plt.subplots(figsize=(9, 5.5))
    ax.semilogx(window_times, rms_values, marker="o")
    ax.axhline(
        peak_rms,
        linestyle="--",
        linewidth=1.0,
        label=f"peak-current RMS = {peak_rms:.3f} mm",
    )
    ax.set_xlabel("integration window T [ns]")
    ax.set_ylabel("|Q|-weighted cluster RMS [mm]")
    ax.set_title("Stage B: apparent cluster width vs observation time")
    ax.legend()
    fig.tight_layout()
    fig.savefig(str(out) + "_cluster_rms_vs_time.png", dpi=200)
    plt.close(fig)

    # ------------------------------------------------------------
    # 5) Compact per-window summary.
    # ------------------------------------------------------------
    compact_rows = []
    for T in windows:
        prof = window_profiles[T]
        compact_rows.append(
            {
                "window_ns": T,
                "is_full_window": bool(np.isclose(T, full_window_ns)),
                "centroid_mm": prof["centroid_mm"],
                "rms_mm": prof["rms_mm"],
                "central_abs_Q_fraction": float(
                    prof["abs_fraction"][strip_ids == 0][0]
                ),
            }
        )
    compact = pd.DataFrame(compact_rows)
    compact["peak_profile_centroid_mm"] = peak_centroid
    compact["peak_profile_rms_mm"] = peak_rms
    compact.to_csv(str(out) + "_summary.csv", index=False)

    print("\n=== STAGE B: ION-CLOUD TIME-WINDOW ANALYSIS ===")
    print(f"input prefix               : {prefix}")
    print(f"waveform dt                : {dt_ns:.6g} ns")
    print(f"full available window      : {full_window_ns:.6g} ns")
    print(f"peak-current centroid      : {peak_centroid:.6f} mm")
    print(f"peak-current RMS           : {peak_rms:.6f} mm")
    print("\nFinite-window charge widths:")
    print(
        compact[
            [
                "window_ns",
                "is_full_window",
                "centroid_mm",
                "rms_mm",
                "central_abs_Q_fraction",
            ]
        ].to_string(index=False)
    )

    print("\nWrote:")
    for suffix in [
        "_by_strip.csv",
        "_summary.csv",
        "_peak_current_profile.png",
        "_signed_Q_profiles.png",
        "_abs_Q_sharing.png",
        "_cluster_rms_vs_time.png",
    ]:
        print(" ", str(out) + suffix)


if __name__ == "__main__":
    main()
