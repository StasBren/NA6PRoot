#!/usr/bin/env python3

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


ELEMENTARY_CHARGE_FC = 1.602176634e-4


def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Plot Stage-B single-ion current, position, speed, cumulative charge "
            "and the two Shockley-Ramo factors v(t) and E_w,eff(t)."
        )
    )
    p.add_argument(
        "--prefix",
        default="phaseB_single_ion_2mm",
        help="Input file prefix (default: phaseB_single_ion_2mm)",
    )
    p.add_argument(
        "--strip",
        default="strip_0",
        help="Strip label to plot (default: strip_0)",
    )
    return p.parse_args()


def main():
    args = parse_args()
    prefix = Path(args.prefix)

    wave_file = Path(str(prefix) + "_waveforms.csv")
    drift_file = Path(str(prefix) + "_driftline.csv")

    wave = pd.read_csv(wave_file)
    drift = pd.read_csv(drift_file)

    current_col = f"{args.strip}_fC_per_ns"
    if current_col not in wave.columns:
        raise KeyError(
            f"Column {current_col!r} not found in {wave_file}. "
            f"Available columns: {list(wave.columns)}"
        )

    t_wave_ns = wave["time_ns"].to_numpy(dtype=float)
    current = wave[current_col].to_numpy(dtype=float)

    if len(t_wave_ns) < 2:
        raise RuntimeError("Need at least two waveform bins.")

    dt_ns = float(np.median(np.diff(t_wave_ns)))
    cumulative_q = np.cumsum(current) * dt_ns

    t_drift_ns = drift["time_ns"].to_numpy(dtype=float)
    v_mm = drift["v_mm"].to_numpy(dtype=float)

    # Local total drift speed from each RKF segment, used in the original
    # diagnostics panel.
    step_mm = drift["step_mm"].to_numpy(dtype=float)
    dt_seg_ns = drift["dt_ns"].to_numpy(dtype=float)
    speed_mm_per_us = np.full_like(step_mm, np.nan, dtype=float)
    good = dt_seg_ns > 0
    speed_mm_per_us[good] = 1000.0 * step_mm[good] / dt_seg_ns[good]

    valid_idx = np.flatnonzero(np.isfinite(speed_mm_per_us))
    if valid_idx.size:
        speed_mm_per_us[: valid_idx[0]] = speed_mm_per_us[valid_idx[0]]

    # Restrict waveform plots to the physical drift interval.
    t_end_ns = float(t_drift_ns[-1])
    mask = t_wave_ns <= t_end_ns + 0.5 * dt_ns
    t_wave_ns_phys = t_wave_ns[mask]
    t_wave_us = t_wave_ns_phys / 1000.0
    current_plot = current[mask]
    q_plot = cumulative_q[mask]

    t_drift_us = t_drift_ns / 1000.0

    # ------------------------------------------------------------------
    # Figure 1: original diagnostics.
    # ------------------------------------------------------------------
    fig, axes = plt.subplots(4, 1, figsize=(9, 10), sharex=True)

    axes[0].plot(t_wave_us, current_plot)
    axes[0].set_ylabel("I [fC/ns]")
    axes[0].set_title(f"Single positive ion: {args.strip}")
    axes[0].grid(alpha=0.3)

    axes[1].plot(t_drift_us, v_mm)
    axes[1].set_ylabel("v position [mm]")
    axes[1].grid(alpha=0.3)

    axes[2].plot(t_drift_us, speed_mm_per_us)
    axes[2].set_ylabel("speed [mm/us]")
    axes[2].grid(alpha=0.3)

    axes[3].plot(t_wave_us, q_plot)
    axes[3].set_ylabel("Q cumulative [fC]")
    axes[3].set_xlabel("time [us]")
    axes[3].grid(alpha=0.3)

    fig.tight_layout()

    out_png = Path(str(prefix) + f"_{args.strip}_diagnostics.png")
    fig.savefig(out_png, dpi=200)
    plt.close(fig)

    # ------------------------------------------------------------------
    # Figure 2: separate the two multiplicative Shockley-Ramo factors.
    #
    # On the central symmetry line (u=w=0), the motion is essentially along
    # local v and the relevant weighting-field component is E_w,v. We infer
    # an effective projected weighting field from
    #
    #   |E_w,eff| = |I| / (e |v_v|)
    #
    # using the Garfield waveform and the local v-velocity from the RKF
    # driftline. This is a diagnostic reconstruction, not an independent
    # field calculation; finite signal bins and interpolation limit the
    # precision close to the endpoints.
    # ------------------------------------------------------------------
    dt_pair_ns = np.diff(t_drift_ns)
    dv_pair_mm = np.diff(v_mm)
    pair_good = dt_pair_ns > 0

    t_mid_ns = 0.5 * (t_drift_ns[1:] + t_drift_ns[:-1])
    vv_mm_per_us = np.full_like(dt_pair_ns, np.nan, dtype=float)
    vv_mm_per_us[pair_good] = (
        1000.0 * np.abs(dv_pair_mm[pair_good]) / dt_pair_ns[pair_good]
    )

    finite_pair = np.isfinite(vv_mm_per_us) & (vv_mm_per_us > 0)
    if np.count_nonzero(finite_pair) < 2:
        raise RuntimeError("Not enough valid drift-line segments to reconstruct v(t).")

    vv_interp_mm_per_us = np.interp(
        t_wave_ns_phys,
        t_mid_ns[finite_pair],
        vv_mm_per_us[finite_pair],
    )
    vv_interp_mm_per_ns = vv_interp_mm_per_us / 1000.0

    ew_eff_per_mm = np.full_like(current_plot, np.nan, dtype=float)
    speed_ok = vv_interp_mm_per_ns > 0
    ew_eff_per_mm[speed_ok] = (
        np.abs(current_plot[speed_ok])
        / (ELEMENTARY_CHARGE_FC * vv_interp_mm_per_ns[speed_ok])
    )

    reconstructed_current = (
        ELEMENTARY_CHARGE_FC
        * vv_interp_mm_per_ns
        * ew_eff_per_mm
    )

    ramo_df = pd.DataFrame(
        {
            "time_us": t_wave_us,
            "current_fC_per_ns": current_plot,
            "v_speed_mm_per_us": vv_interp_mm_per_us,
            "Ew_eff_per_mm": ew_eff_per_mm,
            "reconstructed_abs_current_fC_per_ns": reconstructed_current,
        }
    )
    ramo_csv = Path(str(prefix) + f"_{args.strip}_ramo_factors.csv")
    ramo_df.to_csv(ramo_csv, index=False)

    fig2, axes2 = plt.subplots(3, 1, figsize=(9, 8), sharex=True)

    axes2[0].plot(t_wave_us, vv_interp_mm_per_us)
    axes2[0].set_ylabel("|v_v| [mm/us]")
    axes2[0].set_title(
        f"Shockley-Ramo factors along the ion path: {args.strip}"
    )
    axes2[0].grid(alpha=0.3)

    axes2[1].plot(t_wave_us, ew_eff_per_mm)
    axes2[1].set_ylabel("|Ew,eff| [1/mm]")
    axes2[1].grid(alpha=0.3)

    axes2[2].plot(t_wave_us, np.abs(current_plot), label="|Garfield current|")
    axes2[2].plot(
        t_wave_us,
        reconstructed_current,
        linestyle="--",
        label="e |v_v| |Ew,eff|",
    )
    axes2[2].set_ylabel("|I| [fC/ns]")
    axes2[2].set_xlabel("time [us]")
    axes2[2].grid(alpha=0.3)
    axes2[2].legend()

    fig2.tight_layout()

    ramo_png = Path(str(prefix) + f"_{args.strip}_ramo_factors.png")
    fig2.savefig(ramo_png, dpi=200)
    plt.close(fig2)

    final_q = float(q_plot[-1]) if len(q_plot) else float("nan")
    peak_idx = int(np.argmax(np.abs(current_plot))) if len(current_plot) else 0

    finite_ew = ew_eff_per_mm[np.isfinite(ew_eff_per_mm)]
    ew_min = float(np.min(finite_ew)) if finite_ew.size else float("nan")
    ew_max = float(np.max(finite_ew)) if finite_ew.size else float("nan")

    print(f"Input waveform : {wave_file}")
    print(f"Input driftline: {drift_file}")
    print(f"Strip          : {args.strip}")
    print(f"Drift time     : {t_end_ns/1000.0:.6f} us")
    print(
        f"Peak |current| : {abs(current_plot[peak_idx]):.9e} fC/ns "
        f"at {t_wave_us[peak_idx]:.6f} us"
    )
    print(f"Final Q        : {final_q:.9e} fC")
    print(f"Ew,eff range   : {ew_min:.6g} .. {ew_max:.6g} 1/mm")
    print("Note: Ew,eff is reconstructed from I/(e*v_v); it is not an")
    print("      independent Garfield weighting-field evaluation.")
    print(f"Wrote          : {out_png}")
    print(f"Wrote          : {ramo_png}")
    print(f"Wrote          : {ramo_csv}")


if __name__ == "__main__":
    main()
