#!/usr/bin/env python3

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def parse_args():
    p = argparse.ArgumentParser(
        description="Plot Stage-B single-ion current, position, speed and cumulative charge."
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

    # Local drift speed from each RKF segment.
    step_mm = drift["step_mm"].to_numpy(dtype=float)
    dt_seg_ns = drift["dt_ns"].to_numpy(dtype=float)
    speed_mm_per_us = np.full_like(step_mm, np.nan, dtype=float)
    good = dt_seg_ns > 0
    speed_mm_per_us[good] = 1000.0 * step_mm[good] / dt_seg_ns[good]

    # The first point has no preceding segment. Fill it from the first valid speed
    # so the plotted curve starts at the ion initial time.
    valid_idx = np.flatnonzero(np.isfinite(speed_mm_per_us))
    if valid_idx.size:
        speed_mm_per_us[: valid_idx[0]] = speed_mm_per_us[valid_idx[0]]

    # Restrict waveform plots to the physical drift interval.
    t_end_ns = float(t_drift_ns[-1])
    mask = t_wave_ns <= t_end_ns + 0.5 * dt_ns
    t_wave_us = t_wave_ns[mask] / 1000.0
    current_plot = current[mask]
    q_plot = cumulative_q[mask]

    t_drift_us = t_drift_ns / 1000.0

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

    final_q = float(q_plot[-1]) if len(q_plot) else float("nan")
    peak_idx = int(np.argmax(np.abs(current_plot))) if len(current_plot) else 0

    print(f"Input waveform : {wave_file}")
    print(f"Input driftline: {drift_file}")
    print(f"Strip          : {args.strip}")
    print(f"Drift time     : {t_end_ns/1000.0:.6f} us")
    print(
        f"Peak |current| : {abs(current_plot[peak_idx]):.9e} fC/ns "
        f"at {t_wave_us[peak_idx]:.6f} us"
    )
    print(f"Final Q        : {final_q:.9e} fC")
    print(f"Wrote          : {out_png}")


if __name__ == "__main__":
    main()
