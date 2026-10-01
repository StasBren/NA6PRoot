#!/usr/bin/env python3

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Plot the direct Stage-B Shockley-Ramo verification using Garfield's "
            "real field and weighting field along the ion drift line."
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
        help="Strip label used in the waveform (default: strip_0)",
    )
    return p.parse_args()


def main():
    args = parse_args()
    prefix = Path(args.prefix)

    check_file = Path(str(prefix) + "_ramo_check.csv")
    wave_file = Path(str(prefix) + "_waveforms.csv")

    check = pd.read_csv(check_file)
    wave = pd.read_csv(wave_file)

    wave_col = f"{args.strip}_fC_per_ns"
    if wave_col not in wave.columns:
        raise KeyError(
            f"Column {wave_col!r} not found in {wave_file}. "
            f"Available columns: {list(wave.columns)}"
        )

    t_us = check["t_mid_ns"].to_numpy(dtype=float) / 1000.0
    speed = check["speed_mm_per_us"].to_numpy(dtype=float)
    ereal = check["Ereal_mag_V_per_cm"].to_numpy(dtype=float)
    ew_mag_per_mm = check["Ew_mag_per_cm"].to_numpy(dtype=float) / 10.0

    i_field = check["i_ramo_field_fC_per_ns"].to_numpy(dtype=float)
    i_dphi = check["i_dphi_fC_per_ns"].to_numpy(dtype=float)

    t_wave_us = wave["time_ns"].to_numpy(dtype=float) / 1000.0
    i_garfield = wave[wave_col].to_numpy(dtype=float)

    # Agreement between the local direct-field expression and finite-step
    # weighting-potential difference. Exclude tiny values where relative
    # errors are not meaningful.
    scale = np.maximum(np.abs(i_field), np.abs(i_dphi))
    threshold = max(1e-30, 1e-4 * np.nanmax(scale))
    good = scale > threshold
    rel = np.full_like(i_field, np.nan, dtype=float)
    rel[good] = (i_field[good] - i_dphi[good]) / scale[good]

    fig, axes = plt.subplots(4, 1, figsize=(9, 11), sharex=True)

    axes[0].plot(t_us, ereal)
    axes[0].set_ylabel("|E real| [V/cm]")
    axes[0].set_title("Direct Garfield field check along one-ion trajectory")
    axes[0].grid(alpha=0.3)

    axes[1].plot(t_us, speed)
    axes[1].set_ylabel("speed [mm/us]")
    axes[1].grid(alpha=0.3)

    axes[2].plot(t_us, ew_mag_per_mm)
    axes[2].set_ylabel("|Ew| [1/mm]")
    axes[2].grid(alpha=0.3)

    axes[3].plot(
        t_us,
        i_field,
        label=r"$-e\,\mathbf{v}\cdot\mathbf{E}_w$ (direct field)",
    )
    axes[3].plot(
        t_us,
        i_dphi,
        linestyle="--",
        label=r"$e\,\Delta\phi_w/\Delta t$ (segment)",
    )

    # The Garfield waveform is time-binned (50 ns by default), whereas the
    # two curves above live on adaptive RKF trajectory segments. Plot it as
    # points so the different sampling is visually explicit.
    physical = t_wave_us <= (np.nanmax(t_us) + 0.1)
    axes[3].plot(
        t_wave_us[physical],
        i_garfield[physical],
        marker=".",
        linestyle="none",
        markersize=3,
        alpha=0.6,
        label="Garfield binned signal",
    )
    axes[3].set_ylabel("I [fC/ns]")
    axes[3].set_xlabel("time [us]")
    axes[3].grid(alpha=0.3)
    axes[3].legend()

    fig.tight_layout()
    out_png = Path(str(prefix) + f"_{args.strip}_direct_ramo_check.png")
    fig.savefig(out_png, dpi=200)
    plt.close(fig)

    fig2, ax = plt.subplots(figsize=(9, 4))
    ax.plot(t_us, rel)
    ax.axhline(0.0, linewidth=1)
    ax.set_xlabel("time [us]")
    ax.set_ylabel("(I_field - I_dphi) / max(|I|)")
    ax.set_title("Local consistency: weighting field vs weighting-potential change")
    ax.grid(alpha=0.3)
    fig2.tight_layout()

    out_rel = Path(str(prefix) + f"_{args.strip}_direct_ramo_relative_difference.png")
    fig2.savefig(out_rel, dpi=200)
    plt.close(fig2)

    finite_rel = rel[np.isfinite(rel)]
    if finite_rel.size:
        median_abs = float(np.median(np.abs(finite_rel)))
        p95_abs = float(np.percentile(np.abs(finite_rel), 95))
        max_abs = float(np.max(np.abs(finite_rel)))
    else:
        median_abs = p95_abs = max_abs = float("nan")

    print(f"Input direct check : {check_file}")
    print(f"Input waveform     : {wave_file}")
    print(f"Strip              : {args.strip}")
    print(f"Segments checked   : {len(check)}")
    print(f"Median |relative difference| : {median_abs:.6e}")
    print(f"95% |relative difference|    : {p95_abs:.6e}")
    print(f"Max |relative difference|    : {max_abs:.6e}")
    print("")
    print("The direct-field and dphi curves should agree closely when each")
    print("RKF segment is small enough. The Garfield waveform uses finite time")
    print("bins, so point-by-point agreement with the adaptive-segment curves")
    print("is not expected near sharp features.")
    print(f"Wrote: {out_png}")
    print(f"Wrote: {out_rel}")


if __name__ == "__main__":
    main()
