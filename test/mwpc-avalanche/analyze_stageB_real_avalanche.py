#!/usr/bin/env python3
"""
Compact Stage-B analysis for one real microscopic Garfield avalanche.

Inputs:
  <prefix>_electron_waveforms.csv
  <prefix>_ion_waveforms.csv
  <prefix>_strip_summary.csv
  <prefix>_event_summary.csv

Outputs only three physics plots to avoid overloading the interpretation:
  1) current_components.png
     Real microscopic electron current and explicit positive-ion current
     on strip -1, 0, +1, shown on their natural time scales.

  2) cluster_width_vs_time.png
     |Q|-weighted cluster RMS for ions alone and electron+ion total.

  3) electron_effect_vs_time.png
     Top: fraction of accumulated absolute induced charge coming from
          electrons.
     Bottom: change in cluster RMS caused by adding electrons.

No electronics shaping is applied.
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
        help=(
            "Comma-separated observation windows in ns. "
            "The full ion window is added automatically."
        ),
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


def weighted_rms(x, weights):
    x = np.asarray(x, dtype=float)
    w = np.asarray(weights, dtype=float)
    total = np.sum(w)
    if total <= 0:
        return np.nan, np.nan
    w = w / total
    mean = float(np.sum(w * x))
    rms = float(np.sqrt(np.sum(w * (x - mean) ** 2)))
    return mean, rms


def component_matrix(df, strip_ids, suffix):
    rows = []
    for i in strip_ids:
        col = strip_label(int(i)) + suffix
        if col not in df.columns:
            raise KeyError(f"Missing waveform column {col!r}")
        rows.append(df[col].to_numpy(dtype=float))
    return np.asarray(rows)


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
    windows = [T for T in windows if T < i_full]
    windows.append(i_full)

    q_e_full = np.array(
        [integrate_to(e_times, e_current[j], e_full)
         for j in range(len(strip_ids))],
        dtype=float,
    )

    rows = []
    for T in windows:
        q_e = np.array(
            [
                integrate_to(
                    e_times, e_current[j], min(T, e_full)
                )
                for j in range(len(strip_ids))
            ],
            dtype=float,
        )
        q_i = np.array(
            [
                integrate_to(i_times, i_current[j], T)
                for j in range(len(strip_ids))
            ],
            dtype=float,
        )
        q_total = q_e + q_i

        _, rms_i = weighted_rms(centers, np.abs(q_i))
        _, rms_total = weighted_rms(centers, np.abs(q_total))

        denom = np.sum(np.abs(q_e)) + np.sum(np.abs(q_i))
        e_fraction = (
            float(np.sum(np.abs(q_e)) / denom)
            if denom > 0 else np.nan
        )

        rows.append(
            {
                "window_ns": T,
                "ion_rms_mm": rms_i,
                "electron_plus_ion_rms_mm": rms_total,
                "delta_rms_mm": rms_total - rms_i,
                "electron_abs_charge_fraction": e_fraction,
            }
        )

    result = pd.DataFrame(rows)

    out = (
        Path(args.output_prefix)
        if args.output_prefix
        else Path(str(prefix) + "_analysis")
    )
    result.to_csv(str(out) + "_summary.csv", index=False)

    # ------------------------------------------------------------
    # 1) The actual detector currents, on their natural time scales.
    # ------------------------------------------------------------
    fig, axes = plt.subplots(2, 1, figsize=(9.5, 8.0))
    ax_e, ax_i = axes

    for sid in (-1, 0, 1):
        mask = np.where(strip_ids == sid)[0]
        if len(mask) == 0:
            continue
        j = int(mask[0])
        ax_e.plot(
            e_times,
            e_current[j],
            label=f"strip {sid}",
        )
    ax_e.set_xlabel("time [ns]")
    ax_e.set_ylabel("electron current [fC/ns]")
    ax_e.set_title("Real microscopic avalanche: electron-induced current")
    ax_e.legend()

    for sid in (-1, 0, 1):
        mask = np.where(strip_ids == sid)[0]
        if len(mask) == 0:
            continue
        j = int(mask[0])
        ax_i.plot(
            i_times / 1000.0,
            i_current[j],
            label=f"strip {sid}",
        )
    ax_i.set_xlabel("time [us]")
    ax_i.set_ylabel("positive-ion current [fC/ns]")
    ax_i.set_title("Same avalanche: positive-ion induced current")
    ax_i.legend()

    fig.suptitle(
        "Stage B: one real Garfield avalanche - raw signal components"
    )
    fig.tight_layout()
    fig.savefig(str(out) + "_current_components.png", dpi=200)
    plt.close(fig)

    # ------------------------------------------------------------
    # 2) The quantity tied most directly to cluster formation.
    # ------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(9.5, 5.8))
    ax.semilogx(
        result["window_ns"],
        result["ion_rms_mm"],
        marker="o",
        label="positive ions only",
    )
    ax.semilogx(
        result["window_ns"],
        result["electron_plus_ion_rms_mm"],
        marker="o",
        label="electron + ion",
    )
    ax.set_xlabel("observation window T [ns]")
    ax.set_ylabel("|Q|-weighted cluster RMS [mm]")
    ax.set_title(
        "Stage B: cluster width from one real microscopic avalanche"
    )
    ax.legend()
    fig.tight_layout()
    fig.savefig(str(out) + "_cluster_width_vs_time.png", dpi=200)
    plt.close(fig)

    # ------------------------------------------------------------
    # 3) Direct answer to: what did electrons change?
    # ------------------------------------------------------------
    fig, axes = plt.subplots(2, 1, figsize=(9.5, 8.0), sharex=True)
    ax0, ax1 = axes

    ax0.semilogx(
        result["window_ns"],
        100.0 * result["electron_abs_charge_fraction"],
        marker="o",
    )
    ax0.set_ylabel("electron share of |Q| [%]")
    ax0.set_title(
        "How much of the accumulated induced charge comes from electrons?"
    )

    ax1.semilogx(
        result["window_ns"],
        result["delta_rms_mm"],
        marker="o",
    )
    ax1.axhline(0.0, linewidth=0.8)
    ax1.set_xlabel("observation window T [ns]")
    ax1.set_ylabel("Delta cluster RMS [mm]")
    ax1.set_title(
        "How much do electrons change the reconstructed cluster width?"
    )

    fig.suptitle("Stage B: electron effect in a real Garfield avalanche")
    fig.tight_layout()
    fig.savefig(str(out) + "_electron_effect_vs_time.png", dpi=200)
    plt.close(fig)

    ev = event_summary.iloc[0]

    print("\n=== STAGE B: REAL AVALANCHE ANALYSIS ===")
    print(f"input prefix             : {prefix}")
    print(
        "Garfield avalanche e/ions: "
        f"{int(ev['garfield_electrons'])} / "
        f"{int(ev['garfield_ions'])}"
    )
    print(
        "recorded ion births       : "
        f"{int(ev['recorded_ion_births'])}"
    )
    print(
        "electron time span        : "
        f"{ev['electron_start_min_ns']:.6g} .. "
        f"{ev['electron_end_max_ns']:.6g} ns"
    )
    print(
        "ion-birth sigma_w         : "
        f"{ev['ion_birth_sigma_w_mm']:.6g} mm"
    )
    print(
        "max ion end time          : "
        f"{ev['ion_end_max_ns'] / 1000.0:.6g} us"
    )

    print("\nWhat changes when electrons are included?")
    print(
        result[
            [
                "window_ns",
                "ion_rms_mm",
                "electron_plus_ion_rms_mm",
                "delta_rms_mm",
                "electron_abs_charge_fraction",
            ]
        ].to_string(index=False)
    )

    print("\nWrote only three physics plots:")
    print(" ", str(out) + "_current_components.png")
    print(" ", str(out) + "_cluster_width_vs_time.png")
    print(" ", str(out) + "_electron_effect_vs_time.png")
    print(" ", str(out) + "_summary.csv")


if __name__ == "__main__":
    main()
