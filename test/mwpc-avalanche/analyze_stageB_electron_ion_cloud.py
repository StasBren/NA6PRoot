#!/usr/bin/env python3
"""
Stage-B analysis of the synthetic electron + positive-ion cloud.

The generator stores electrons and ions on separate time grids:
  * prompt electrons: ps-scale fine grid;
  * slow positive ions: ns-scale grid covering the full drift.

For a finite observation window T, the detector-level induced charge is

    Q_total,k(T) = Q_e,k(T) + Q_i,k(T).

For the windows used in this study (50 ns and above), the electron pulse has
already ended, so Q_e,k(T) is simply the full prompt-electron induced charge.

This script deliberately applies NO frontend electronics, shaping, threshold,
ADC, or VMM model.

Important limitation of the SYNTHETIC electron cloud:
all electrons are created at t=0 with the same user-chosen initial energy and a
radially inward initial direction.  Therefore the absolute ps-scale electron
peak current is NOT a physical prediction of an avalanche.  It is extremely
sensitive to this artificial initial condition and to the electron time binning.

Robust observables in this synthetic e+i step are instead:
  * the full prompt-electron induced charge (Shockley-Ramo endpoint quantity);
  * finite-window charge sharing Q_e(T) + Q_i(T) for T well above the electron
    drift time;
  * cluster RMS versus observation time;
  * comparison to the ions-only baseline.

The raw electron waveform is retained only as a diagnostic.  A quantitative
electron peak-current prediction requires a real microscopic avalanche with
physical creation times and secondary-electron kinematics.
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
        default="stageB_ei_cloud_2mm",
        help="Prefix from mwpc_phase_b_electron_ion_cloud",
    )
    p.add_argument(
        "--ions-only-prefix",
        default="stageB_cloud_2mm",
        help=(
            "Optional ions-only baseline prefix. Use 'none' to disable "
            "(default: stageB_cloud_2mm)."
        ),
    )
    p.add_argument(
        "--windows-ns",
        default="50,100,200,500,1000,5000",
        help=(
            "Comma-separated finite integration windows in ns. "
            "The full ion waveform window is added automatically."
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


def weighted_centroid_rms(x, positive_weights):
    x = np.asarray(x, dtype=float)
    w = np.asarray(positive_weights, dtype=float)
    total = float(np.sum(w))
    if total <= 0:
        return np.nan, np.nan
    w = w / total
    centroid = float(np.sum(w * x))
    rms = float(np.sqrt(np.sum(w * (x - centroid) ** 2)))
    return centroid, rms


def infer_dt_ns(times):
    if len(times) < 2:
        raise RuntimeError("Need at least two waveform bins.")
    return float(np.median(np.diff(times)))


def integrate_to(times_ns, current, window_ns):
    """
    Integrate a uniform, piecewise-constant waveform over [0, window_ns].
    Samples are interpreted as bin centres.
    """
    times_ns = np.asarray(times_ns, dtype=float)
    current = np.asarray(current, dtype=float)
    dt = infer_dt_ns(times_ns)

    left = times_ns - 0.5 * dt
    right = times_ns + 0.5 * dt
    overlap = np.minimum(right, window_ns) - np.maximum(left, 0.0)
    overlap = np.clip(overlap, 0.0, None)
    return float(np.sum(current * overlap))


def load_component_matrix(df, strip_ids, suffix):
    rows = []
    for i in strip_ids:
        col = strip_label(int(i)) + suffix
        if col not in df.columns:
            raise KeyError(f"Missing column {col!r}")
        rows.append(df[col].to_numpy(dtype=float))
    return np.asarray(rows)


def load_ions_only_matrix(prefix, strip_ids):
    wave_file = Path(str(prefix) + "_waveforms.csv")
    if not wave_file.exists():
        return None, None
    wave = pd.read_csv(wave_file)
    times = wave["time_ns"].to_numpy(dtype=float)
    rows = []
    for i in strip_ids:
        col = strip_label(int(i)) + "_fC_per_ns"
        if col not in wave.columns:
            raise KeyError(f"Missing ions-only column {col!r}")
        rows.append(wave[col].to_numpy(dtype=float))
    return times, np.asarray(rows)


def window_label(T, full_window):
    if np.isclose(T, full_window):
        return "full"
    if T >= 1000:
        return f"{T / 1000.0:g} us"
    return f"{T:g} ns"


def main():
    args = parse_args()

    prefix = Path(args.prefix)
    electron_file = Path(str(prefix) + "_electron_waveforms.csv")
    ion_file = Path(str(prefix) + "_ion_waveforms.csv")
    summary_file = Path(str(prefix) + "_summary.csv")

    e = pd.read_csv(electron_file)
    ion = pd.read_csv(ion_file)
    summary = pd.read_csv(summary_file).sort_values("center_w_mm")

    e_times = e["time_ns"].to_numpy(dtype=float)
    i_times = ion["time_ns"].to_numpy(dtype=float)
    e_dt = infer_dt_ns(e_times)
    i_dt = infer_dt_ns(i_times)

    electron_window_ns = float(e_times[-1] + 0.5 * e_dt)
    ion_full_window_ns = float(i_times[-1] + 0.5 * i_dt)

    strip_ids = summary["strip"].to_numpy(dtype=int)
    centers = summary["center_w_mm"].to_numpy(dtype=float)

    e_current = load_component_matrix(
        e, strip_ids, "_electron_fC_per_ns"
    )
    i_current = load_component_matrix(
        ion, strip_ids, "_ion_fC_per_ns"
    )

    finite_windows = sorted(
        {
            float(x.strip())
            for x in args.windows_ns.split(",")
            if x.strip()
        }
    )
    if any(T <= 0 for T in finite_windows):
        raise ValueError("All windows must be positive.")

    windows = [T for T in finite_windows if T < ion_full_window_ns]
    windows.append(ion_full_window_ns)

    # ------------------------------------------------------------
    # Raw current profiles.
    # ------------------------------------------------------------
    peak_e = np.max(np.abs(e_current), axis=1)
    peak_i = np.max(np.abs(i_current), axis=1)

    # Controlled prompt-total approximation: add the essentially constant
    # initial ion current to the ps-scale electron waveform.
    ion_initial = i_current[:, 0]
    prompt_total = e_current + ion_initial[:, None]
    peak_prompt_total = np.max(np.abs(prompt_total), axis=1)

    e_peak_centroid, e_peak_rms = weighted_centroid_rms(centers, peak_e)
    i_peak_centroid, i_peak_rms = weighted_centroid_rms(centers, peak_i)
    t_peak_centroid, t_peak_rms = weighted_centroid_rms(
        centers, peak_prompt_total
    )

    # ------------------------------------------------------------
    # Finite-window charges.
    # ------------------------------------------------------------
    profiles = {}
    comparison_rows = []

    ions_only_prefix = args.ions_only_prefix.strip()
    baseline_times = None
    baseline_current = None
    if ions_only_prefix.lower() != "none":
        baseline_times, baseline_current = load_ions_only_matrix(
            Path(ions_only_prefix), strip_ids
        )

    # Full prompt-electron charge by strip.
    q_e_full = np.array(
        [
            integrate_to(e_times, e_current[j], electron_window_ns)
            for j in range(len(strip_ids))
        ],
        dtype=float,
    )

    for T in windows:
        q_e = np.array(
            [
                integrate_to(
                    e_times,
                    e_current[j],
                    min(T, electron_window_ns),
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

        abs_total = np.abs(q_total)
        frac_total = (
            abs_total / np.sum(abs_total)
            if np.sum(abs_total) > 0
            else np.zeros_like(abs_total)
        )
        total_centroid, total_rms = weighted_centroid_rms(
            centers, abs_total
        )

        ion_centroid, ion_rms = weighted_centroid_rms(
            centers, np.abs(q_i)
        )

        baseline_rms = np.nan
        baseline_centroid = np.nan
        baseline_central = np.nan
        max_ion_baseline_rel_diff = np.nan

        if baseline_times is not None:
            q_base = np.array(
                [
                    integrate_to(
                        baseline_times, baseline_current[j], T
                    )
                    for j in range(len(strip_ids))
                ],
                dtype=float,
            )
            baseline_centroid, baseline_rms = weighted_centroid_rms(
                centers, np.abs(q_base)
            )
            base_abs = np.abs(q_base)
            if np.sum(base_abs) > 0:
                base_frac = base_abs / np.sum(base_abs)
                baseline_central = float(
                    base_frac[strip_ids == 0][0]
                )

            denom = np.maximum(np.abs(q_base), 1.0e-30)
            max_ion_baseline_rel_diff = float(
                np.max(np.abs(q_i - q_base) / denom)
            )

        profiles[T] = {
            "qe": q_e,
            "qi": q_i,
            "qtotal": q_total,
            "fraction": frac_total,
            "centroid": total_centroid,
            "rms": total_rms,
            "ion_rms": ion_rms,
        }

        comparison_rows.append(
            {
                "window_ns": T,
                "is_full_window": bool(
                    np.isclose(T, ion_full_window_ns)
                ),
                "ei_centroid_mm": total_centroid,
                "ei_rms_mm": total_rms,
                "ei_central_abs_Q_fraction": float(
                    frac_total[strip_ids == 0][0]
                ),
                "ion_component_centroid_mm": ion_centroid,
                "ion_component_rms_mm": ion_rms,
                "ions_only_centroid_mm": baseline_centroid,
                "ions_only_rms_mm": baseline_rms,
                "ions_only_central_abs_Q_fraction": baseline_central,
                "delta_rms_ei_minus_ions_only_mm": (
                    total_rms - baseline_rms
                    if np.isfinite(baseline_rms)
                    else np.nan
                ),
                "max_rel_diff_ion_component_vs_baseline": (
                    max_ion_baseline_rel_diff
                ),
            }
        )

    result = pd.DataFrame(comparison_rows)

    out = (
        Path(args.output_prefix)
        if args.output_prefix
        else Path(str(prefix) + "_analysis")
    )
    result.to_csv(str(out) + "_summary.csv", index=False)

    # ------------------------------------------------------------
    # 1) Stage-B component spatial diagnostics.
    #
    # Do NOT present the synthetic electron raw peak as a physical avalanche
    # observable.  Instead show:
    #   - ion peak-current profile (well-defined in the current ion model);
    #   - full prompt-electron induced-charge profile (endpoint-robust).
    # ------------------------------------------------------------
    qe_centroid, qe_rms = weighted_centroid_rms(
        centers, np.abs(q_e_full)
    )

    fig, axes = plt.subplots(2, 1, figsize=(9.5, 8.5), sharex=True)
    ax0, ax1 = axes

    ax0.plot(
        centers,
        peak_i,
        marker="o",
        label="ion peak |I|",
    )
    ax0.set_ylabel("ion peak induced current [fC/ns]")
    ax0.set_title(
        "Ion peak-current spatial profile "
        f"(RMS = {i_peak_rms:.3f} mm)"
    )
    ax0.legend()

    ax1.plot(
        centers,
        q_e_full,
        marker="o",
        label="full prompt-electron Q",
    )
    ax1.axhline(0.0, linewidth=0.8)
    ax1.set_xlabel("strip center w [mm]")
    ax1.set_ylabel("electron induced charge [fC]")
    ax1.set_title(
        "Prompt-electron induced-charge spatial profile "
        f"(abs-weight RMS = {qe_rms:.3f} mm)"
    )
    ax1.legend()

    fig.suptitle(
        "Stage B: reliable component spatial diagnostics "
        "(synthetic e+i cloud)"
    )
    fig.tight_layout()
    fig.savefig(
        str(out) + "_component_spatial_diagnostics.png",
        dpi=200,
    )
    plt.close(fig)

    # ------------------------------------------------------------
    # 2) Signed e+i charge profiles versus observation time.
    # ------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(10, 6))
    for T in windows:
        ax.plot(
            centers,
            profiles[T]["qtotal"],
            marker="o",
            label=window_label(T, ion_full_window_ns),
        )
    ax.axhline(0.0, linewidth=0.8)
    ax.set_xlabel("strip center w [mm]")
    ax.set_ylabel("Q_e,k(T) + Q_i,k(T) [fC]")
    ax.set_title(
        "Stage B: electron + ion signed induced charge vs integration window"
    )
    ax.legend()
    fig.tight_layout()
    fig.savefig(str(out) + "_signed_Q_profiles.png", dpi=200)
    plt.close(fig)

    # ------------------------------------------------------------
    # 3) e+i charge sharing versus observation time.
    # ------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(10, 6))
    for T in windows:
        ax.plot(
            centers,
            100.0 * profiles[T]["fraction"],
            marker="o",
            label=window_label(T, ion_full_window_ns),
        )
    ax.set_xlabel("strip center w [mm]")
    ax.set_ylabel("|Q_total,k(T)| / sum |Q_total,j(T)| [%]")
    ax.set_title(
        "Stage B: electron + ion time evolution of charge sharing"
    )
    ax.legend()
    fig.tight_layout()
    fig.savefig(str(out) + "_abs_Q_sharing.png", dpi=200)
    plt.close(fig)

    # ------------------------------------------------------------
    # 4) Main comparison: cluster RMS vs observation time.
    # ------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(9.5, 5.8))
    ax.semilogx(
        result["window_ns"],
        result["ei_rms_mm"],
        marker="o",
        label="electron + ion",
    )
    ax.semilogx(
        result["window_ns"],
        result["ion_component_rms_mm"],
        marker="o",
        label="ion component in same run",
    )
    if np.any(np.isfinite(result["ions_only_rms_mm"])):
        ax.semilogx(
            result["window_ns"],
            result["ions_only_rms_mm"],
            marker="o",
            linestyle="--",
            label="previous ions-only baseline",
        )
    ax.axhline(
        i_peak_rms,
        linestyle=":",
        linewidth=1.2,
        label=f"ion peak-current RMS = {i_peak_rms:.3f} mm",
    )
    ax.set_xlabel("integration window T [ns]")
    ax.set_ylabel("|Q|-weighted cluster RMS [mm]")
    ax.set_title(
        "Stage B: does adding prompt electrons change the apparent cluster width?"
    )
    ax.legend()
    fig.tight_layout()
    fig.savefig(str(out) + "_cluster_rms_comparison.png", dpi=200)
    plt.close(fig)

    # ------------------------------------------------------------
    # 5) How important is the prompt-electron charge as the observation
    #    window grows?
    # ------------------------------------------------------------
    electron_abs_fraction = []
    for T in windows:
        qe = profiles[T]["qe"]
        qi = profiles[T]["qi"]
        denom = np.sum(np.abs(qe)) + np.sum(np.abs(qi))
        electron_abs_fraction.append(
            np.sum(np.abs(qe)) / denom if denom > 0 else np.nan
        )

    fig, ax = plt.subplots(figsize=(9.5, 5.8))
    ax.semilogx(
        windows,
        electron_abs_fraction,
        marker="o",
    )
    ax.set_xlabel("integration window T [ns]")
    ax.set_ylabel(
        "sum |Q_e,k| / (sum |Q_e,k| + sum |Q_i,k|)"
    )
    ax.set_title(
        "Stage B: relative prompt-electron charge contribution vs time"
    )
    fig.tight_layout()
    fig.savefig(
        str(out) + "_electron_charge_fraction_vs_time.png",
        dpi=200,
    )
    plt.close(fig)

    # ------------------------------------------------------------
    # 6) Full-charge component profile.
    # ------------------------------------------------------------
    Tfull = ion_full_window_ns
    q_i_full = profiles[Tfull]["qi"]
    q_total_full = profiles[Tfull]["qtotal"]

    fig, ax = plt.subplots(figsize=(9.5, 5.8))
    ax.plot(centers, q_e_full, marker="o", label="electron Q")
    ax.plot(centers, q_i_full, marker="o", label="ion Q")
    ax.plot(centers, q_total_full, marker="o", label="electron + ion Q")
    ax.axhline(0.0, linewidth=0.8)
    ax.set_xlabel("strip center w [mm]")
    ax.set_ylabel("full integrated induced charge [fC]")
    ax.set_title("Stage B: full-charge component decomposition")
    ax.legend()
    fig.tight_layout()
    fig.savefig(str(out) + "_full_charge_components.png", dpi=200)
    plt.close(fig)

    print("\n=== STAGE B: ELECTRON + ION CLOUD ANALYSIS ===")
    print(f"input prefix                  : {prefix}")
    print(f"electron waveform dt          : {e_dt:.6g} ns")
    print(f"electron waveform window      : {electron_window_ns:.6g} ns")
    print(f"ion waveform dt               : {i_dt:.6g} ns")
    print(f"ion waveform full window      : {ion_full_window_ns:.6g} ns")
    print(f"ion peak-profile RMS          : {i_peak_rms:.6f} mm")
    print(f"electron charge-profile RMS   : {qe_rms:.6f} mm")
    print(
        "electron raw peak current     : diagnostic only; not interpreted "
        "physically in the synthetic cloud"
    )
    if baseline_times is not None:
        print(f"ions-only baseline            : {ions_only_prefix}")

    print("\nFinite-window cluster comparison:")
    cols = [
        "window_ns",
        "ei_rms_mm",
        "ion_component_rms_mm",
        "ions_only_rms_mm",
        "delta_rms_ei_minus_ions_only_mm",
        "ei_central_abs_Q_fraction",
    ]
    print(result[cols].to_string(index=False))

    print("\nWrote:")
    for suffix in [
        "_summary.csv",
        "_component_spatial_diagnostics.png",
        "_signed_Q_profiles.png",
        "_abs_Q_sharing.png",
        "_cluster_rms_comparison.png",
        "_electron_charge_fraction_vs_time.png",
        "_full_charge_components.png",
    ]:
        print(" ", str(out) + suffix)


if __name__ == "__main__":
    main()
