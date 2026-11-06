#!/usr/bin/env python3
"""Render a reproducible MS0/DET21 parent-event A comparison from saved LUT500 summary.

No Geant4, Heed, Garfield, or LUT reruns. Counts and moments are read as saved.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

DEFAULT_SUMMARY = Path(
    "test_runs/mwpc_stage_c/pi_lut500_parent_asymmetry_root_v1/summary.json"
)
DEFAULT_OUTPUT = Path("ms0_det21_parent_A_p3_LUT500.png")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path, default=DEFAULT_SUMMARY)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--station", type=int, default=0)
    parser.add_argument("--det", type=int, default=21)
    args = parser.parse_args()

    d = json.loads(args.summary.read_text(encoding="utf-8"))
    if d.get("partial_exposure"):
        raise ValueError("Refusing a partial-exposure/pilot summary")
    if not np.isclose(float(d["pitch_mm"]), 3.0):
        raise ValueError("This figure is specified for p=3 mm only")

    record = d["stations"][str(args.station)]["chambers"][str(args.det)]
    cfg = d["configs"]
    if len(cfg) != 2 or [c["field_wire_diameter_um"] for c in cfg] != [0, 75]:
        raise ValueError("Expected the matched noFW / grounded FW75 configurations")
    edges = np.asarray(d["asymmetry"]["hist_edges"], dtype=float)
    if len(edges) != 42 or not np.allclose(edges, np.linspace(-1, 1, 42)):
        raise ValueError("Expected the original 41-bin histogram on [-1, 1]")
    results = [record["configurations"][c["key"]] for c in cfg]

    fig, ax = plt.subplots(figsize=(9.2, 6.4))
    fig.subplots_adjust(left=0.115, right=0.97, top=0.77, bottom=0.19)
    for r, name, color, style in zip(
        results,
        ["No field wire", "Grounded field wire 75 μm"],
        ["#2457AE", "#D32D43"], ["-", "--"],
    ):
        counts = np.asarray(r["hist"], dtype=int)
        valid = int(r["valid_A"])
        if counts.size != 41 or counts.sum() != valid or valid <= 0:
            raise ValueError("Invalid saved histogram counts")
        ax.stairs(counts / valid, edges, color=color, linewidth=2,
                  linestyle=style,
                  label=f"{name}  ·  RMS = {r['rms_A']:.4f}")

    nofw, fw = results
    reduction = 100 * (1 - float(nofw["rms_A"]) / float(fw["rms_A"]))
    nslots = int(d["parents_generated_in_input_chunks"])
    fig.text(0.115, 0.953, f"MS{args.station} / DET {args.det} — parent-event P(A)",
             weight="bold", size=16)
    fig.text(0.115, 0.916,
             "NA60+ Simulation  ·  Work in progress  ·  p = 3 mm, LUT500",
             size=10)
    fig.text(0.115, 0.878,
             f"{nslots:,} source-parent slots  ·  noFW RMS is {reduction:.1f}% lower",
             size=10)
    ax.set(xlim=(-1, 1), xlabel="Electron-count asymmetry A (dimensionless)",
           ylabel="Fraction of valid parent events / 0.049")
    ax.set_ylim(0, max(float(np.max(r["hist"]) / r["valid_A"])
                       for r in results) * 1.19)
    ax.grid(axis="y", color="0.90", linewidth=0.8)
    ax.tick_params(top=True, right=True, direction="in")
    ax.legend(loc="upper right", frameon=False, fontsize=9)
    caption = (f"Valid A: {nofw['valid_A']:,} (noFW), {fw['valid_A']:,} (FW75); "
               f"mean A: {nofw['mean_A']:+.4f}, {fw['mean_A']:+.4f}. "
               "Parent-aggregated A, not cathode charge or spatial resolution.")
    fig.text(0.115, 0.080,
             f"Valid A: {nofw['valid_A']:,} (noFW), {fw['valid_A']:,} (FW75)"
             f"  ·  mean A: {nofw['mean_A']:+.4f}, {fw['mean_A']:+.4f}",
             fontsize=9, color="0.32")
    fig.text(0.115, 0.050,
             "Parent-aggregated electron-count asymmetry; not cathode charge or spatial resolution.",
             fontsize=8.4, color="0.32")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=210, bbox_inches="tight", pad_inches=0.20)
    print(f"Saved: {args.output}")
    print(caption)
    plt.close(fig)


if __name__ == "__main__":
    main()
