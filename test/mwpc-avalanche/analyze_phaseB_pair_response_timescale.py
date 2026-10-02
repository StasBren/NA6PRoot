#!/usr/bin/env python3

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Scan a generic electronics response time for the one-electron + "
            "one-ion signal. The response obeys dS/dt = I - S/tau, so small "
            "tau is current-like and large tau approaches charge integration."
        )
    )
    p.add_argument(
        "--prefix",
        default="phaseB_pair_2mm",
        help="Input prefix used by mwpc_phase_b_electron_ion_pair",
    )
    p.add_argument(
        "--electrode",
        default="strip_0",
        help=(
            "Electrode to analyse: strip_0, strip_p1, strip_m1, ... "
            "or readout_plane (default: strip_0)"
        ),
    )
    p.add_argument(
        "--taus-ns",
        default=(
            "0.005,0.01,0.02,0.05,0.1,0.2,0.5,1,2,5,10,20,50,100,"
            "200,500,1000,2000,5000,10000,20000,50000,100000"
        ),
        help="Comma-separated response time constants in ns",
    )
    p.add_argument(
        "--output-prefix",
        default=None,
        help="Output prefix (default: <prefix>_<electrode>_timescale)",
    )
    return p.parse_args()


def infer_dt_ns(times):
    if len(times) < 2:
        raise RuntimeError("Need at least two waveform bins.")
    return float(np.median(np.diff(times)))


def sample_piecewise_constant(t_ns, times_ns, values, dt_ns):
    t0 = float(times_ns[0] - 0.5 * dt_ns)
    t1 = float(times_ns[-1] + 0.5 * dt_ns)
    if t_ns < t0 or t_ns >= t1:
        return 0.0
    idx = int(np.floor((t_ns - t0) / dt_ns))
    idx = max(0, min(idx, len(values) - 1))
    return float(values[idx])


def response_peaks(e_times, e_current, i_times, i_current, tau_ns):
    dt_e = infer_dt_ns(e_times)
    dt_i = infer_dt_ns(i_times)

    e_edges = np.concatenate(
        ([e_times[0] - 0.5 * dt_e], e_times + 0.5 * dt_e)
    )
    i_edges = np.concatenate(
        ([i_times[0] - 0.5 * dt_i], i_times + 0.5 * dt_i)
    )

    boundaries = np.unique(
        np.concatenate(
            (
                np.maximum(e_edges, 0.0),
                np.maximum(i_edges, 0.0),
            )
        )
    )
    boundaries.sort()

    se = 0.0
    si = 0.0
    peak_e = 0.0
    peak_i = 0.0
    peak_total = 0.0
    t_peak_e = 0.0
    t_peak_i = 0.0
    t_peak_total = 0.0

    for a, b in zip(boundaries[:-1], boundaries[1:]):
        dt = float(b - a)
        if dt <= 0.0:
            continue
        tm = 0.5 * (a + b)

        ie = sample_piecewise_constant(
            tm, e_times, e_current, dt_e
        )
        ii = sample_piecewise_constant(
            tm, i_times, i_current, dt_i
        )

        decay = np.exp(-dt / tau_ns)
        gain = tau_ns * (1.0 - decay)

        se = se * decay + ie * gain
        si = si * decay + ii * gain
        st = se + si

        if abs(se) > peak_e:
            peak_e = abs(se)
            t_peak_e = b
        if abs(si) > peak_i:
            peak_i = abs(si)
            t_peak_i = b
        if abs(st) > peak_total:
            peak_total = abs(st)
            t_peak_total = b

    return {
        "tau_ns": tau_ns,
        "peak_e_fC": peak_e,
        "peak_i_fC": peak_i,
        "peak_total_fC": peak_total,
        "t_peak_e_ns": t_peak_e,
        "t_peak_i_ns": t_peak_i,
        "t_peak_total_ns": t_peak_total,
        "electron_to_ion_peak_ratio": (
            peak_e / peak_i if peak_i > 0.0 else np.inf
        ),
    }


def main():
    args = parse_args()
    prefix = Path(args.prefix)

    efile = Path(str(prefix) + "_electron_waveforms.csv")
    ifile = Path(str(prefix) + "_ion_waveforms.csv")

    e = pd.read_csv(efile)
    i = pd.read_csv(ifile)

    if args.electrode == "readout_plane":
        e_col = "readout_plane_electron_fC_per_ns"
        i_col = "readout_plane_ion_fC_per_ns"
    else:
        e_col = f"{args.electrode}_electron_fC_per_ns"
        i_col = f"{args.electrode}_ion_fC_per_ns"

    if e_col not in e.columns:
        raise KeyError(f"Missing {e_col!r} in {efile}")
    if i_col not in i.columns:
        raise KeyError(f"Missing {i_col!r} in {ifile}")

    e_times = e["time_ns"].to_numpy(dtype=float)
    e_current = e[e_col].to_numpy(dtype=float)
    i_times = i["time_ns"].to_numpy(dtype=float)
    i_current = i[i_col].to_numpy(dtype=float)

    taus = np.array(
        [float(x.strip()) for x in args.taus_ns.split(",") if x.strip()],
        dtype=float,
    )
    if len(taus) == 0 or np.any(taus <= 0.0):
        raise ValueError("All response times must be positive.")

    rows = [
        response_peaks(
            e_times, e_current, i_times, i_current, float(tau)
        )
        for tau in taus
    ]
    df = pd.DataFrame(rows)

    if args.output_prefix:
        out_prefix = Path(args.output_prefix)
    else:
        out_prefix = Path(
            str(prefix) + f"_{args.electrode}_timescale"
        )

    csv_out = Path(str(out_prefix) + ".csv")
    png_out = Path(str(out_prefix) + ".png")
    df.to_csv(csv_out, index=False)

    fig, axes = plt.subplots(2, 1, figsize=(8.5, 8.0), sharex=True)
    ax0, ax1 = axes

    ax0.loglog(
        df["tau_ns"],
        df["peak_e_fC"],
        marker="o",
        label="electron component",
    )
    ax0.loglog(
        df["tau_ns"],
        df["peak_i_fC"],
        marker="o",
        label="ion component",
    )
    ax0.loglog(
        df["tau_ns"],
        df["peak_total_fC"],
        marker="o",
        label="total response",
    )
    ax0.set_ylabel("peak response S [fC]")
    ax0.set_title(
        "Generic leaky-integrator response: dS/dt = I - S/tau"
    )
    ax0.grid(alpha=0.25, which="both")
    ax0.legend()

    ratio = df["electron_to_ion_peak_ratio"].replace(
        [np.inf, -np.inf], np.nan
    )
    ax1.semilogx(
        df["tau_ns"],
        ratio,
        marker="o",
    )
    ax1.axhline(1.0, linestyle="--", linewidth=1.0)
    ax1.set_xlabel("response time tau [ns]")
    ax1.set_ylabel("electron peak / ion peak")
    ax1.set_title(
        "Small tau = current-like; large tau = increasingly integrating"
    )
    ax1.grid(alpha=0.25, which="both")

    fig.suptitle(
        f"One-pair signal timescale scan: {args.electrode}",
        fontsize=14,
    )
    fig.tight_layout()
    fig.savefig(png_out, dpi=180)
    plt.close(fig)

    finite = df[np.isfinite(df["electron_to_ion_peak_ratio"])]
    crossover = None
    if len(finite) >= 2:
        x = finite["tau_ns"].to_numpy(dtype=float)
        y = finite["electron_to_ion_peak_ratio"].to_numpy(dtype=float) - 1.0
        crossings = np.flatnonzero(y[:-1] * y[1:] <= 0.0)
        if len(crossings):
            j = int(crossings[0])
            crossover = np.sqrt(x[j] * x[j + 1])

    print(f"Electrode      : {args.electrode}")
    print(f"Electron input : {efile}")
    print(f"Ion input      : {ifile}")
    print(f"Wrote          : {csv_out}")
    print(f"Wrote          : {png_out}")
    if crossover is not None:
        print(
            "Approx. electron/ion peak crossover tau "
            f": {crossover:.6g} ns"
        )
    else:
        print(
            "No electron/ion peak crossover found inside the scanned tau range."
        )


if __name__ == "__main__":
    main()
