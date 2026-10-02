#!/usr/bin/env python3

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Plot the separate electron and ion signals, their sum, and the "
            "cumulative induced charge for the Phase-B one-pair test."
        )
    )
    p.add_argument(
        "--prefix",
        default="phaseB_pair_2mm",
        help="Input prefix used by mwpc_phase_b_electron_ion_pair",
    )
    p.add_argument(
        "--strip",
        default="strip_0",
        help="Strip label to plot (default: strip_0)",
    )
    p.add_argument(
        "--output",
        default=None,
        help="Output PNG (default: <prefix>_<strip>_signals.png)",
    )
    return p.parse_args()


def main():
    args = parse_args()
    prefix = Path(args.prefix)

    efile = Path(str(prefix) + "_electron_waveforms.csv")
    ifile = Path(str(prefix) + "_ion_waveforms.csv")
    cfile = Path(str(prefix) + "_combined_waveforms.csv")
    sfile = Path(str(prefix) + "_summary.csv")

    e = pd.read_csv(efile)
    i = pd.read_csv(ifile)
    c = pd.read_csv(cfile)
    s = pd.read_csv(sfile)

    e_col = f"{args.strip}_electron_fC_per_ns"
    i_col = f"{args.strip}_ion_fC_per_ns"
    ce_col = f"{args.strip}_electron_fC_per_ns"
    ci_col = f"{args.strip}_ion_fC_per_ns"
    ct_col = f"{args.strip}_total_fC_per_ns"

    for df, col, name in [
        (e, e_col, efile),
        (i, i_col, ifile),
        (c, ce_col, cfile),
        (c, ci_col, cfile),
        (c, ct_col, cfile),
    ]:
        if col not in df.columns:
            raise KeyError(f"Missing column {col!r} in {name}")

    # Find the requested strip in the summary.
    if args.strip == "strip_0":
        strip_index = 0
    elif args.strip.startswith("strip_p"):
        strip_index = int(args.strip.replace("strip_p", ""))
    elif args.strip.startswith("strip_m"):
        strip_index = -int(args.strip.replace("strip_m", ""))
    else:
        raise ValueError(
            "--strip must look like strip_0, strip_p1, strip_m1, ..."
        )

    row = s.loc[s["strip"] == strip_index]
    if len(row) != 1:
        raise RuntimeError(
            f"Could not identify strip {strip_index} in {sfile}"
        )
    row = row.iloc[0]

    te_ns = e["time_ns"].to_numpy(dtype=float)
    ie = e[e_col].to_numpy(dtype=float)

    ti_us = i["time_ns"].to_numpy(dtype=float) / 1000.0
    ii = i[i_col].to_numpy(dtype=float)

    tc_ns = c["time_ns"].to_numpy(dtype=float)
    tc_us = tc_ns / 1000.0
    ice = c[ce_col].to_numpy(dtype=float)
    ici = c[ci_col].to_numpy(dtype=float)
    ict = c[ct_col].to_numpy(dtype=float)

    dt_c_ns = float(np.median(np.diff(tc_ns)))
    qe_cum = np.cumsum(ice) * dt_c_ns
    qi_cum = np.cumsum(ici) * dt_c_ns
    qt_cum = np.cumsum(ict) * dt_c_ns

    q_e_expected = float(row["q_electron_expected_fC"])
    q_i_expected = float(row["q_ion_expected_fC"])
    q_t_expected = float(row["q_total_expected_fC"])

    residual = ict - ice - ici
    max_residual = float(np.max(np.abs(residual)))

    fig, axes = plt.subplots(2, 2, figsize=(12, 7))
    ax_e, ax_i, ax_t, ax_q = axes.flat

    # Prompt electron on its own fine time scale.
    ax_e.plot(te_ns, ie, linewidth=1.7)
    ax_e.set_xlabel("time [ns]")
    ax_e.set_ylabel("electron current [fC/ns]")
    ax_e.set_title(f"{args.strip}: prompt electron")
    ax_e.grid(alpha=0.25)

    # Slow positive-ion component.
    ax_i.plot(ti_us, ii, linewidth=1.7)
    ax_i.set_xlabel("time [us]")
    ax_i.set_ylabel("ion current [fC/ns]")
    ax_i.set_title(f"{args.strip}: positive-ion tail")
    ax_i.grid(alpha=0.25)

    # Both contributions on the same coarse physical-time grid.
    ax_t.plot(tc_us, ice, label="electron", linewidth=1.4)
    ax_t.plot(tc_us, ici, label="ion", linewidth=1.4)
    ax_t.plot(tc_us, ict, label="total", linewidth=1.8)
    ax_t.set_xlabel("time [us]")
    ax_t.set_ylabel("induced current [fC/ns]")
    ax_t.set_title("Same electrode: separate components and sum")
    ax_t.legend()
    ax_t.grid(alpha=0.25)

    # Integrated signal charge makes the endpoint Ramo relation visible.
    ax_q.plot(tc_us, qe_cum, label="electron", linewidth=1.4)
    ax_q.plot(tc_us, qi_cum, label="ion", linewidth=1.4)
    ax_q.plot(tc_us, qt_cum, label="total", linewidth=1.8)
    ax_q.axhline(q_t_expected, linestyle="--", linewidth=1.0,
                 label="endpoint Ramo total")
    ax_q.set_xlabel("time [us]")
    ax_q.set_ylabel("cumulative signal charge [fC]")
    ax_q.set_title("Integrated induced charge")
    ax_q.legend(fontsize=8)
    ax_q.grid(alpha=0.25)

    fig.suptitle(
        "One electron + one positive ion: prompt and slow signal components",
        fontsize=14,
    )
    fig.tight_layout()

    if args.output:
        out = Path(args.output)
    else:
        out = Path(str(prefix) + f"_{args.strip}_signals.png")

    fig.savefig(out, dpi=180)
    plt.close(fig)

    print(f"Strip                  : {args.strip}")
    print(f"phi_start              : {float(row['phi_start']):.9g}")
    print(f"Qe expected            : {q_e_expected:.9e} fC")
    print(f"Qi expected            : {q_i_expected:.9e} fC")
    print(f"Qtotal expected        : {q_t_expected:.9e} fC")
    print(f"Qtotal Garfield coarse : {float(row['q_total_coarse_fC']):.9e} fC")
    print(f"max |Itot-Ie-Ii|       : {max_residual:.3e} fC/ns")
    print(f"Wrote                  : {out}")


if __name__ == "__main__":
    main()
