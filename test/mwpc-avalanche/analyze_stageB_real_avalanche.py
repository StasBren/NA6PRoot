#!/usr/bin/env python3
"""
Compact, interpretable Stage-B analysis for one real microscopic Garfield avalanche.

The goal is NOT to produce many derived plots. We keep only three questions:

1) What do the raw electron and positive-ion currents actually look like?
2) How much accumulated induced charge comes from electrons vs ions as the
   observation window grows?
3) What spatial strip-sharing pattern do electrons, ions, and their sum produce
   in the early 25--200 ns window relevant for subsequent reconstruction work?

No electronics shaping, threshold, ADC, or VMM transfer function is applied.
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
        default="stageB_real_avalanche",
        help="Prefix produced by mwpc_phase_b_real_avalanche",
    )
    p.add_argument(
        "--windows-ns",
        default="25,50,100,200,500,1000,5000",
        help="Comma-separated observation windows in ns.",
    )
    p.add_argument(
        "--sharing-windows-ns",
        default="25,200",
        help="Two windows used for the spatial-sharing comparison.",
    )
    p.add_argument(
        "--output-prefix",
        default=None,
        help="Default: <prefix>_analysis",
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


def component_matrix(df, strip_ids, suffix):
    rows = []
    for i in strip_ids:
        col = strip_label(int(i)) + suffix
        if col not in df.columns:
            raise KeyError(f"Missing waveform column {col!r}")
        rows.append(df[col].to_numpy(dtype=float))
    return np.asarray(rows)


def spatial_fraction(q):
    w = np.abs(np.asarray(q, dtype=float))
    s = np.sum(w)
    return w / s if s > 0 else np.zeros_like(w)


def main():
    args = parse_args()
    prefix = Path(args.prefix)

    e = pd.read_csv(str(prefix) + "_electron_waveforms.csv")
    ion = pd.read_csv(str(prefix) + "_ion_waveforms.csv")
    strip_summary = pd.read_csv(
        str(prefix) + "_strip_summary.csv"
    ).sort_values("center_w_mm")
    event_summary = pd.read_csv(str(prefix) + "_event_summary.csv")

    e_times = e["time_ns"].to_numpy(dtype=float)
    i_times = ion["time_ns"].to_numpy(dtype=float)
    e_dt = infer_dt(e_times)
    i_dt = infer_dt(i_times)
    e_full = float(e_times[-1] + 0.5 * e_dt)
    i_full = float(i_times[-1] + 0.5 * i_dt)

    strip_ids = strip_summary["strip"].to_numpy(dtype=int)
    centers = strip_summary["center_w_mm"].to_numpy(dtype=float)

    e_current = component_matrix(
        e, strip_ids, "_electron_fC_per_ns"
    )
    i_current = component_matrix(
        ion, strip_ids, "_ion_fC_per_ns"
    )

    windows = sorted(
        {
            float(x.strip())
            for x in args.windows_ns.split(",")
            if x.strip()
        }
    )
    if any(T <= 0 for T in windows):
        raise ValueError("All observation windows must be positive.")
    windows = [T for T in windows if T <= i_full]

    sharing_windows = [
        float(x.strip())
        for x in args.sharing_windows_ns.split(",")
        if x.strip()
    ]
    if len(sharing_windows) != 2:
        raise ValueError("--sharing-windows-ns must contain exactly two values.")

    rows = []
    profiles = {}

    for T in windows:
        qe = np.array(
            [
                integrate_to(e_times, e_current[j], min(T, e_full))
                for j in range(len(strip_ids))
            ],
            dtype=float,
        )
        qi = np.array(
            [
                integrate_to(i_times, i_current[j], T)
                for j in range(len(strip_ids))
            ],
            dtype=float,
        )
        qt = qe + qi

        ae = float(np.sum(np.abs(qe)))
        ai = float(np.sum(np.abs(qi)))
        at = float(np.sum(np.abs(qt)))
        comp_denom = ae + ai

        ft = spatial_fraction(qt)

        def frac_for(strip):
            idx = np.where(strip_ids == strip)[0]
            if len(idx) == 0:
                return np.nan
            return float(ft[int(idx[0])])

        rows.append(
            {
                "window_ns": T,
                "sum_abs_Qe_fC": ae,
                "sum_abs_Qi_fC": ai,
                "sum_abs_Qtotal_fC": at,
                "electron_component_fraction": (
                    ae / comp_denom if comp_denom > 0 else np.nan
                ),
                "ion_component_fraction": (
                    ai / comp_denom if comp_denom > 0 else np.nan
                ),
                "total_spatial_fraction_strip_m2": frac_for(-2),
                "total_spatial_fraction_strip_m1": frac_for(-1),
                "total_spatial_fraction_strip_0": frac_for(0),
                "total_spatial_fraction_strip_p1": frac_for(1),
                "total_spatial_fraction_strip_p2": frac_for(2),
            }
        )

        profiles[T] = {
            "qe": qe,
            "qi": qi,
            "qt": qt,
            "fe": spatial_fraction(qe),
            "fi": spatial_fraction(qi),
            "ft": ft,
        }

    result = pd.DataFrame(rows)

    out = (
        Path(args.output_prefix)
        if args.output_prefix
        else Path(str(prefix) + "_analysis")
    )
    result.to_csv(str(out) + "_summary.csv", index=False)

    # ------------------------------------------------------------------
    # 1) Raw currents on their natural time scales.
    # Electron panel is AUTO-ZOOMED around the actual avalanche burst.
    # ------------------------------------------------------------------
    fig, axes = plt.subplots(2, 1, figsize=(9.5, 8.2))
    ax_e, ax_i = axes

    key_indices = []
    for sid in (-1, 0, 1):
        idx = np.where(strip_ids == sid)[0]
        if len(idx):
            key_indices.append((sid, int(idx[0])))

    global_e_peak = 0.0
    for _, j in key_indices:
        global_e_peak = max(
            global_e_peak,
            float(np.max(np.abs(e_current[j])))
        )

    if global_e_peak > 0:
        combined = np.max(
            np.abs(e_current[[j for _, j in key_indices], :]),
            axis=0,
        )
        active = np.flatnonzero(combined > 0.01 * global_e_peak)
        if len(active):
            i0 = max(0, int(active[0]) - 8)
            i1 = min(len(e_times) - 1, int(active[-1]) + 8)
            electron_xlim = (e_times[i0], e_times[i1])
        else:
            electron_xlim = (0.0, e_full)
    else:
        electron_xlim = (0.0, e_full)

    line_styles = {-1: "--", 0: "-", 1: ":"}

    for sid, j in key_indices:
        ax_e.plot(
            e_times,
            e_current[j],
            linestyle=line_styles[sid],
            label=f"strip {sid}",
        )
    ax_e.set_xlim(*electron_xlim)
    ax_e.set_xlabel("time [ns]")
    ax_e.set_ylabel("electron current [fC/ns]")
    ax_e.set_title(
        "Microscopic electron avalanche current (zoomed on the burst)"
    )
    ax_e.legend()

    for sid, j in key_indices:
        ax_i.plot(
            i_times / 1000.0,
            i_current[j],
            linestyle=line_styles[sid],
            label=f"strip {sid}",
        )
    ax_i.set_xlabel("time [us]")
    ax_i.set_ylabel("positive-ion current [fC/ns]")
    ax_i.set_title("Positive-ion current from the same avalanche")
    ax_i.legend()

    fig.suptitle(
        "Stage B: raw induced-current components from one real avalanche"
    )
    fig.tight_layout()
    fig.savefig(str(out) + "_current_components.png", dpi=200)
    plt.close(fig)

    # ------------------------------------------------------------------
    # 2) Signal composition versus observation time.
    # Top: absolute accumulated component magnitudes.
    # Bottom: the same information as a clean electron-vs-ion fraction.
    # ------------------------------------------------------------------
    fig, axes = plt.subplots(2, 1, figsize=(9.5, 8.0), sharex=True)
    ax0, ax1 = axes

    ax0.loglog(
        result["window_ns"],
        result["sum_abs_Qe_fC"],
        marker="o",
        label="electron component",
    )
    ax0.loglog(
        result["window_ns"],
        result["sum_abs_Qi_fC"],
        marker="o",
        label="positive-ion component",
    )
    ax0.loglog(
        result["window_ns"],
        result["sum_abs_Qtotal_fC"],
        marker="o",
        label="combined e + ion signal",
    )
    ax0.set_ylabel("sum over strips |Q_k(T)| [fC]")
    ax0.set_title("Accumulated induced signal: electrons, ions, and total")
    ax0.legend()

    ax1.semilogx(
        result["window_ns"],
        100.0 * result["electron_component_fraction"],
        marker="o",
        label="electron fraction",
    )
    ax1.semilogx(
        result["window_ns"],
        100.0 * result["ion_component_fraction"],
        marker="o",
        label="ion fraction",
    )
    ax1.set_ylim(0.0, 102.0)
    ax1.set_xlabel("observation window T [ns]")
    ax1.set_ylabel("component fraction [%]")
    ax1.set_title("Who supplies the accumulated induced charge?")
    ax1.legend()

    fig.suptitle("Stage B: electron + ion signal composition")
    fig.tight_layout()
    fig.savefig(str(out) + "_signal_composition_vs_time.png", dpi=200)
    plt.close(fig)

    # ------------------------------------------------------------------
    # 3) Spatial sharing at two early windows.
    # This is the resolution-relevant result: do e and i encode the same
    # strip pattern, and what does the total 3-strip cluster look like?
    # ------------------------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.2), sharey=True)

    for ax, T in zip(axes, sharing_windows):
        if T not in profiles:
            qe = np.array(
                [
                    integrate_to(e_times, e_current[j], min(T, e_full))
                    for j in range(len(strip_ids))
                ],
                dtype=float,
            )
            qi = np.array(
                [
                    integrate_to(i_times, i_current[j], T)
                    for j in range(len(strip_ids))
                ],
                dtype=float,
            )
            qt = qe + qi
            prof = {
                "fe": spatial_fraction(qe),
                "fi": spatial_fraction(qi),
                "ft": spatial_fraction(qt),
            }
        else:
            prof = profiles[T]

        ax.plot(
            centers,
            100.0 * prof["fe"],
            marker="o",
            linestyle="--",
            label="electron",
        )
        ax.plot(
            centers,
            100.0 * prof["fi"],
            marker="o",
            linestyle=":",
            label="positive ion",
        )
        ax.plot(
            centers,
            100.0 * prof["ft"],
            marker="o",
            linewidth=2.0,
            label="total e + ion",
        )
        ax.set_xlabel("strip center w [mm]")
        ax.set_title(f"T = {T:g} ns")

    axes[0].set_ylabel("|Q_k(T)| / sum_j |Q_j(T)| [%]")
    axes[0].legend()
    fig.suptitle(
        "Stage B: spatial strip sharing in the early signal window"
    )
    fig.tight_layout()
    fig.savefig(str(out) + "_early_spatial_sharing.png", dpi=200)
    plt.close(fig)

    ev = event_summary.iloc[0]

    print("\n=== STAGE B: REAL AVALANCHE - COMPACT INTERPRETATION ===")
    print(
        "Garfield avalanche e/ions : "
        f"{int(ev['garfield_electrons'])} / "
        f"{int(ev['garfield_ions'])}"
    )
    print(
        "ion-birth sigma_w         : "
        f"{1000.0 * ev['ion_birth_sigma_w_mm']:.3f} um"
    )
    print(
        "electron avalanche time   : "
        f"{ev['electron_start_min_ns']:.6g} .. "
        f"{ev['electron_end_max_ns']:.6g} ns"
    )
    print("\nSignal composition and total spatial sharing:")
    display_cols = [
        "window_ns",
        "electron_component_fraction",
        "ion_component_fraction",
        "total_spatial_fraction_strip_m1",
        "total_spatial_fraction_strip_0",
        "total_spatial_fraction_strip_p1",
    ]
    display = result[display_cols].copy()
    for col in display_cols[1:]:
        display[col] *= 100.0
    print(display.to_string(index=False))

    print("\nWrote three physics plots:")
    print(" ", str(out) + "_current_components.png")
    print(" ", str(out) + "_signal_composition_vs_time.png")
    print(" ", str(out) + "_early_spatial_sharing.png")
    print(" ", str(out) + "_summary.csv")


if __name__ == "__main__":
    main()
