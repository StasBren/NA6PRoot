#!/usr/bin/env python3
"""Plot complete MS0 station-level per-J/psi-parent asymmetry from frozen LUT500 data.

One statistical entry is one J/psi parent with >=1 gas crossing in MS0.
For each parent, up/down electron endpoints are pooled over ALL MS0
chambers and crossings before A=(N_up-N_down)/(N_up+N_down) is formed.
The source is the immutable dashboard summary; no physics reruns.
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
DEFAULT_OUTPUT = Path(
    "test/mwpc-stage-c/figures/ms0_station_parent_A_p3_LUT500.png"
)


def make_plot(summary: Path, output: Path) -> None:
    data = json.loads(summary.read_text(encoding="utf-8"))
    if data.get("partial_exposure") or data.get("number_of_source_chunks") != 2:
        raise ValueError("Expected complete two-chunk, 100k-parent exposure")
    if data.get("parents_generated_in_input_chunks") != 100000:
        raise ValueError("Expected 100,000 generated J/psi parents")
    if abs(float(data["pitch_mm"]) - 3.) > 1e-9:
        raise ValueError("Expected p=3 mm")

    # Crucial: station_parent, not chambers[DET] and not crossing-level data.
    station = data["stations"]["0"]["station_parent"]
    configs = data["configs"]
    if sorted((int(c["field_wire_diameter_um"]) for c in configs)) != [0,75]:
        raise ValueError("Expected noFW and 75 micron grounded FW")
    nofw_cfg = next(c["key"] for c in configs if int(c["field_wire_diameter_um"]) == 0)
    fw_cfg = next(c["key"] for c in configs if int(c["field_wire_diameter_um"]) == 75)
    nofw = station["configurations"][nofw_cfg]
    fw = station["configurations"][fw_cfg]

    edges = np.asarray(data["asymmetry"]["hist_edges"], dtype=float)
    if edges.size != 42 or not np.allclose(edges, np.linspace(-1, 1, 42)):
        raise ValueError("Expected original 41-bin histogram")
    for name, r in (("noFW", nofw), ("FW75", fw)):
        hist = np.asarray(r["hist"], dtype=int)
        if hist.size != 41 or hist.sum() != int(r["valid_A"]):
            raise ValueError(f"Bad {name} bin normalization")
        if int(r["parents_with_gas_crossing"]) != int(r["valid_A"]) + int(r["undefined_A"]):
            raise ValueError(f"Bad {name} exposed parent accounting")

    fig, ax = plt.subplots(figsize=(9.5, 6.3))
    fig.subplots_adjust(left=0.11, right=0.97, top=0.79, bottom=0.21)
    for name, r, color, style in (
        ("No field wire", nofw, "#2557AD", "-"),
        ("Grounded field wire 75 μm", fw, "#D53047", "--"),
    ):
        y = np.asarray(r["hist"], dtype=float) / int(r["valid_A"])
        ax.stairs(y, edges, label=f"{name} (RMS {r['rms_A']:.4f})",
                  color=color, linestyle=style, linewidth=2.1)

    reduction = 100. * (1. - float(nofw["rms_A"]) / float(fw["rms_A"]))
    fig.text(0.11, 0.955, "MS0 — entire station: parent-event P(A)",
             fontsize=16, weight="bold")
    fig.text(0.11, 0.915, "NA60+ Simulation · p = 3 mm · LUT500 · Work in progress",
             fontsize=10)
    fig.text(0.11, 0.879,
             f"100,000 generated J/ψ parents · RMS reduction without field wire: {reduction:.1f}%",
             fontsize=10)

    ax.set_xlabel("Electron-count asymmetry A (dimensionless)", fontsize=11)
    ax.set_ylabel("Fraction of valid parent events / 0.049", fontsize=11)
    ax.set_xlim(-1, 1)
    ax.grid(axis="y", color="0.9")
    ax.tick_params(direction="in", top=True, right=True)
    ax.legend(frameon=False, loc="upper right", fontsize=9)
    ax.set_ylim(0, max(np.max(np.asarray(r["hist"],float) / r["valid_A"])
                       for r in (nofw, fw)) * 1.18)
    fig.text(0.11, 0.100,
             "Parents with MS0 crossing: "
             f"{nofw['parents_with_gas_crossing']:,}; "
             f"valid A: {nofw['valid_A']:,} / {fw['valid_A']:,} (noFW / FW75); "
             f"undefined: {nofw['undefined_A']} / {fw['undefined_A']}",
             fontsize=9)
    fig.text(0.11, 0.068,
             f"Mean A: {nofw['mean_A']:+.4f} / {fw['mean_A']:+.4f}  ·  "
             f"classified electrons: "
             f"{int(nofw['total_up'])+int(nofw['total_down']):,} / "
             f"{int(fw['total_up'])+int(fw['total_down']):,}",
             fontsize=9)
    fig.text(0.11, 0.038,
             "Primary-electron endpoint counts; not cathode signal or spatial resolution.",
             fontsize=9, color="0.38")
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=210, bbox_inches="tight", pad_inches=0.18)
    plt.close(fig)
    print(f"Saved station-level figure: {output}")
    print(f"Population RMS noFW/FW75: {nofw['rms_A']:.6f}/{fw['rms_A']:.6f}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path, default=DEFAULT_SUMMARY)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    make_plot(args.summary, args.output)


if __name__ == "__main__":
    main()
