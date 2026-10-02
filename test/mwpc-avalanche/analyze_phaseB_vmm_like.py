#!/usr/bin/env python3
"""
Apply a documented-VMM-timescale, third-order semi-Gaussian proxy shaper to the
Phase-B one-electron + one-ion current waveforms.

What is taken from VMM3a documentation:
  * third-order shaper;
  * selectable peaking times: 25, 50, 100, 200 ns;
  * selectable gain (handled here as a simple mV/fC scale factor).

What is NOT claimed:
  * this is NOT a transistor-level or exact DDF transfer function;
  * VMM3a mild tail cancellation / bipolar mode are NOT modeled here because
    the public user guide exposes the configuration bits but not enough circuit
    constants to reconstruct those transfer functions uniquely.

The proxy response to an instantaneous input charge Q is

    V(t) = G * Q * g(t),

with

    g(t) = (e^2 / 4) * (t/tau)^2 * exp(-t/tau),
    tau  = t_peak / 2,

so that g(t_peak) = 1.  Therefore an ideal delta-like charge pulse Q peaks at
G*Q, while a long detector current is convolved with the same shaping kernel.

The input current is charge-conservingly rebinned before convolution.  This is
important because the electron signal can be orders of magnitude shorter than
the VMM peaking time.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def parse_args():
    p = argparse.ArgumentParser(
        description="VMM-like third-order shaping of Phase-B pair currents."
    )
    p.add_argument(
        "--prefix",
        default="phaseB_pair_2mm",
        help="Input prefix used by mwpc_phase_b_electron_ion_pair.",
    )
    p.add_argument(
        "--electrode",
        default="strip_0",
        help=(
            "strip_0, strip_p1, strip_m1, ... or readout_plane "
            "(default: strip_0)"
        ),
    )
    p.add_argument(
        "--peaking-times-ns",
        default="25,50,100,200",
        help="Comma-separated peaking times in ns.",
    )
    p.add_argument(
        "--gain-mv-per-fc",
        type=float,
        default=1.0,
        help=(
            "VMM gain used only to convert equivalent input charge to mV. "
            "Default: 1 mV/fC."
        ),
    )
    p.add_argument(
        "--dt-ns",
        type=float,
        default=0.5,
        help="Common shaping grid step in ns (default: 0.5 ns).",
    )
    p.add_argument(
        "--analysis-tmax-us",
        type=float,
        default=None,
        help=(
            "Optional shaping-window end in us. By default it is chosen from "
            "the non-zero ion waveform plus 10 times the longest peaking time."
        ),
    )
    p.add_argument(
        "--early-window-us",
        type=float,
        default=1.0,
        help="Early-time plot window in us (default: 1 us).",
    )
    p.add_argument(
        "--output-prefix",
        default=None,
        help="Output prefix (default: <prefix>_<electrode>_vmm_like).",
    )
    return p.parse_args()


def infer_dt_ns(times: np.ndarray) -> float:
    if len(times) < 2:
        raise RuntimeError("Need at least two waveform bins.")
    return float(np.median(np.diff(times)))


def cumulative_charge_edges(
    times_ns: np.ndarray, current_fc_per_ns: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    dt = infer_dt_ns(times_ns)
    edges = np.concatenate(
        ([times_ns[0] - 0.5 * dt], times_ns + 0.5 * dt)
    )
    qcum = np.concatenate(
        ([0.0], np.cumsum(current_fc_per_ns * dt))
    )
    return edges, qcum


def rebin_to_charge(
    times_ns: np.ndarray,
    current_fc_per_ns: np.ndarray,
    dest_edges_ns: np.ndarray,
) -> np.ndarray:
    """Charge-conserving rebin of piecewise-constant current."""
    src_edges, src_qcum = cumulative_charge_edges(
        times_ns, current_fc_per_ns
    )
    q_at_dest = np.interp(
        dest_edges_ns,
        src_edges,
        src_qcum,
        left=0.0,
        right=float(src_qcum[-1]),
    )
    return np.diff(q_at_dest)


def kernel_third_order_semigaussian(
    peaking_time_ns: float,
    dt_ns: float,
    cutoff_peaking_times: float = 12.0,
) -> np.ndarray:
    """
    Unit-peak third-order semi-Gaussian proxy.

    g(t) = (e^2/4) (t/tau)^2 exp(-t/tau), tau = t_peak/2.
    Hence max(g)=1 at t=t_peak.
    """
    tau = 0.5 * peaking_time_ns
    tmax = cutoff_peaking_times * peaking_time_ns
    n = int(np.ceil(tmax / dt_ns)) + 1
    t = np.arange(n, dtype=float) * dt_ns
    x = t / tau
    return (np.e**2 / 4.0) * x**2 * np.exp(-x)


def fft_convolve_same_length(
    q_bins_fc: np.ndarray, kernel: np.ndarray
) -> np.ndarray:
    """Causal linear convolution, truncated to the input-grid length."""
    n_full = len(q_bins_fc) + len(kernel) - 1
    n_fft = 1 << int(np.ceil(np.log2(max(n_full, 2))))
    qf = np.fft.rfft(q_bins_fc, n=n_fft)
    kf = np.fft.rfft(kernel, n=n_fft)
    full = np.fft.irfft(qf * kf, n=n_fft)[:n_full]
    return full[: len(q_bins_fc)]


def choose_tmax_ns(
    ion_times_ns: np.ndarray,
    ion_current: np.ndarray,
    max_tp_ns: float,
    user_tmax_us: float | None,
) -> float:
    if user_tmax_us is not None:
        return 1000.0 * user_tmax_us

    peak = float(np.max(np.abs(ion_current))) if len(ion_current) else 0.0
    if peak > 0.0:
        active = np.flatnonzero(np.abs(ion_current) > peak * 1.0e-8)
        if len(active):
            last = float(ion_times_ns[active[-1]])
        else:
            last = float(ion_times_ns[-1])
    else:
        last = float(ion_times_ns[-1])

    return last + 10.0 * max_tp_ns


def component_columns(electrode: str) -> tuple[str, str]:
    if electrode == "readout_plane":
        return (
            "readout_plane_electron_fC_per_ns",
            "readout_plane_ion_fC_per_ns",
        )
    return (
        f"{electrode}_electron_fC_per_ns",
        f"{electrode}_ion_fC_per_ns",
    )


def abs_peak(y: np.ndarray) -> tuple[int, float]:
    idx = int(np.argmax(np.abs(y)))
    return idx, float(y[idx])


def main():
    args = parse_args()

    if args.dt_ns <= 0.0:
        raise ValueError("--dt-ns must be positive.")
    if args.gain_mv_per_fc <= 0.0:
        raise ValueError("--gain-mv-per-fc must be positive.")

    tps = np.array(
        [
            float(x.strip())
            for x in args.peaking_times_ns.split(",")
            if x.strip()
        ],
        dtype=float,
    )
    if len(tps) == 0 or np.any(tps <= 0.0):
        raise ValueError("All peaking times must be positive.")

    prefix = Path(args.prefix)
    efile = Path(str(prefix) + "_electron_waveforms.csv")
    ifile = Path(str(prefix) + "_ion_waveforms.csv")

    e = pd.read_csv(efile)
    i = pd.read_csv(ifile)

    e_col, i_col = component_columns(args.electrode)
    if e_col not in e.columns:
        raise KeyError(
            f"Missing {e_col!r} in {efile}. "
            "If using readout_plane, rerun the updated pair simulation first."
        )
    if i_col not in i.columns:
        raise KeyError(
            f"Missing {i_col!r} in {ifile}. "
            "If using readout_plane, rerun the updated pair simulation first."
        )

    e_times = e["time_ns"].to_numpy(dtype=float)
    e_current = e[e_col].to_numpy(dtype=float)
    i_times = i["time_ns"].to_numpy(dtype=float)
    i_current = i[i_col].to_numpy(dtype=float)

    electron_source_dt = infer_dt_ns(e_times)
    ion_source_dt = infer_dt_ns(i_times)

    tmax_ns = choose_tmax_ns(
        i_times,
        i_current,
        float(np.max(tps)),
        args.analysis_tmax_us,
    )

    edges = np.arange(
        0.0, tmax_ns + args.dt_ns, args.dt_ns, dtype=float
    )
    if edges[-1] < tmax_ns:
        edges = np.append(edges, tmax_ns)
    times = edges[:-1] + 0.5 * np.diff(edges)

    # Rebin detector current to charge per common time bin [fC].
    q_e = rebin_to_charge(e_times, e_current, edges)
    q_i = rebin_to_charge(i_times, i_current, edges)

    q_e_input = float(np.sum(q_e))
    q_i_input = float(np.sum(q_i))

    rows = []
    shaped = {}

    for tp in tps:
        g = kernel_third_order_semigaussian(
            float(tp), args.dt_ns
        )

        # Equivalent input-charge response [fC].
        s_e_fc = fft_convolve_same_length(q_e, g)
        s_i_fc = fft_convolve_same_length(q_i, g)
        s_t_fc = s_e_fc + s_i_fc

        # Convert to voltage using selectable VMM gain.
        v_e = args.gain_mv_per_fc * s_e_fc
        v_i = args.gain_mv_per_fc * s_i_fc
        v_t = args.gain_mv_per_fc * s_t_fc

        shaped[float(tp)] = (v_e, v_i, v_t)

        ie, pe = abs_peak(v_e)
        ii, pi = abs_peak(v_i)
        it, pt = abs_peak(v_t)

        e_at_total = float(v_e[it])
        i_at_total = float(v_i[it])
        total_at_peak = float(v_t[it])

        denom = abs(e_at_total) + abs(i_at_total)
        if denom > 0.0:
            e_fraction = abs(e_at_total) / denom
            i_fraction = abs(i_at_total) / denom
        else:
            e_fraction = np.nan
            i_fraction = np.nan

        rows.append(
            {
                "peaking_time_ns": float(tp),
                "gain_mV_per_fC": args.gain_mv_per_fc,
                "electron_input_charge_fC": q_e_input,
                "ion_input_charge_fC": q_i_input,
                "electron_component_peak_mV": pe,
                "electron_component_peak_time_ns": float(times[ie]),
                "ion_component_peak_mV": pi,
                "ion_component_peak_time_ns": float(times[ii]),
                "total_peak_mV": pt,
                "total_peak_time_ns": float(times[it]),
                "electron_at_total_peak_mV": e_at_total,
                "ion_at_total_peak_mV": i_at_total,
                "electron_abs_fraction_at_total_peak": e_fraction,
                "ion_abs_fraction_at_total_peak": i_fraction,
            }
        )

    result = pd.DataFrame(rows)

    if args.output_prefix:
        out = Path(args.output_prefix)
    else:
        out = Path(str(prefix) + f"_{args.electrode}_vmm_like")

    csv_out = Path(str(out) + "_summary.csv")
    early_png = Path(str(out) + "_early.png")
    full_png = Path(str(out) + "_full.png")
    summary_png = Path(str(out) + "_summary.png")
    result.to_csv(csv_out, index=False)

    # ---------------------------------------------------------------
    # Early-time figure: where the prompt electron contribution lives.
    # ---------------------------------------------------------------
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), sharex=True)
    axes = axes.ravel()
    early_ns = 1000.0 * args.early_window_us
    early_mask = times <= early_ns

    for ax, tp in zip(axes, tps):
        ve, vi, vt = shaped[float(tp)]
        ax.plot(
            times[early_mask] / 1000.0,
            ve[early_mask],
            label="electron",
        )
        ax.plot(
            times[early_mask] / 1000.0,
            vi[early_mask],
            label="ion",
        )
        ax.plot(
            times[early_mask] / 1000.0,
            vt[early_mask],
            label="total",
        )
        ax.set_title(f"t_peak = {tp:g} ns")
        ax.set_xlabel("time [us]")
        ax.set_ylabel("shaped output [mV]")
        ax.grid(alpha=0.25)

    axes[0].legend()
    fig.suptitle(
        f"VMM-like 3rd-order shaping: early response, {args.electrode}"
    )
    fig.tight_layout()
    fig.savefig(early_png, dpi=180)
    plt.close(fig)

    # ---------------------------------------------------------------
    # Full-time figure: shows how the slow ion contribution evolves.
    # ---------------------------------------------------------------
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), sharex=True)
    axes = axes.ravel()

    for ax, tp in zip(axes, tps):
        ve, vi, vt = shaped[float(tp)]
        ax.plot(times / 1000.0, ve, label="electron")
        ax.plot(times / 1000.0, vi, label="ion")
        ax.plot(times / 1000.0, vt, label="total")
        ax.set_title(f"t_peak = {tp:g} ns")
        ax.set_xlabel("time [us]")
        ax.set_ylabel("shaped output [mV]")
        ax.grid(alpha=0.25)

    axes[0].legend()
    fig.suptitle(
        f"VMM-like 3rd-order shaping: full response, {args.electrode}"
    )
    fig.tight_layout()
    fig.savefig(full_png, dpi=180)
    plt.close(fig)

    # ---------------------------------------------------------------
    # Summary figure: component peak sizes and composition of total peak.
    # ---------------------------------------------------------------
    fig = plt.figure(figsize=(9, 8))
    ax0 = fig.add_subplot(2, 1, 1)
    ax1 = fig.add_subplot(2, 1, 2)

    ax0.plot(
        result["peaking_time_ns"],
        np.abs(result["electron_component_peak_mV"]),
        marker="o",
        label="electron component peak",
    )
    ax0.plot(
        result["peaking_time_ns"],
        np.abs(result["ion_component_peak_mV"]),
        marker="o",
        label="ion component peak",
    )
    ax0.plot(
        result["peaking_time_ns"],
        np.abs(result["total_peak_mV"]),
        marker="o",
        label="total peak",
    )
    ax0.set_ylabel("|peak shaped output| [mV]")
    ax0.set_title("Component peak amplitudes")
    ax0.grid(alpha=0.25)
    ax0.legend()

    ax1.plot(
        result["peaking_time_ns"],
        result["electron_abs_fraction_at_total_peak"],
        marker="o",
        label="electron fraction",
    )
    ax1.plot(
        result["peaking_time_ns"],
        result["ion_abs_fraction_at_total_peak"],
        marker="o",
        label="ion fraction",
    )
    ax1.set_xlabel("VMM peaking time [ns]")
    ax1.set_ylabel("absolute fraction at total peak")
    ax1.set_ylim(-0.03, 1.03)
    ax1.set_title("Who contributes to the amplitude actually measured at the peak?")
    ax1.grid(alpha=0.25)
    ax1.legend()

    fig.suptitle(
        f"VMM-like shaping summary: {args.electrode}"
    )
    fig.tight_layout()
    fig.savefig(summary_png, dpi=180)
    plt.close(fig)

    print("\n=== VMM-LIKE THIRD-ORDER SHAPER ===")
    print(f"electrode              : {args.electrode}")
    print(f"gain                   : {args.gain_mv_per_fc:g} mV/fC")
    print(f"electron source dt     : {electron_source_dt:g} ns")
    print(f"ion source dt          : {ion_source_dt:g} ns")
    print(f"common shaping dt      : {args.dt_ns:g} ns")
    print(f"analysis window        : 0 .. {tmax_ns / 1000.0:.6g} us")
    print(f"electron input charge  : {q_e_input:.9e} fC")
    print(f"ion input charge       : {q_i_input:.9e} fC")
    if ion_source_dt > float(np.min(tps)) / 5.0:
        print(
            "\nWARNING: ion waveform is too coarsely sampled for the "
            "shortest requested peaking time. For 25 ns shaping, rerun the "
            "pair simulation with --ion-dt-ns 2 (or at most about 5 ns)."
        )

    print(
        "\nNOTE: third-order semi-Gaussian proxy matched to documented "
        "VMM peaking times; not exact DDF circuit simulation."
    )
    print(
        "NOTE: VMM mild tail cancellation / bipolar mode are not included "
        "because their exact transfer functions are not specified in the "
        "public user guide."
    )
    print("\n" + result.to_string(index=False))
    print("\nWrote:")
    print(f"  {csv_out}")
    print(f"  {early_png}")
    print(f"  {full_png}")
    print(f"  {summary_png}")


if __name__ == "__main__":
    main()
