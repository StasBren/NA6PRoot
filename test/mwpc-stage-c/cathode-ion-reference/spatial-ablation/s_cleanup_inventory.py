#!/usr/bin/env python3
"""Inventory of an S-ladder run for a FUTURE cleanup decision. Deletes nothing; prints paths, sizes and what would
be needed to reconstruct each item.

  s_cleanup_inventory.py RUN [--out INVENTORY.md]
"""
import argparse
from pathlib import Path


def size(p):
    p = Path(p)
    return p.stat().st_size if p.is_file() else sum(f.stat().st_size for f in p.rglob("*") if f.is_file())


RULES = [  # (glob relative to RUN, category, keep?, reconstruction)
    ("green_paths_*", "transport (RKF grid paths)", "KEEP", "p3_green_grid_paths + frozen exporter source + mobility; reused for any new strip geometry"),
    ("green_*/G.npy", "Green table (paths x weighting LUT)", "KEEP while kernels may be rebuilt", "s_green_table.py from green_paths_* + LUT (minutes)"),
    ("kernel_*", "harmonic kernels (production input)", "KEEP", "s_kernel.py build from green_*"),
    ("muon_*/frozen_export", "path export metadata (vertices already removed)", "candidate", "frozen exporter on ion_births.root (RKF, minutes per muon)"),
    ("muon_*/T1_export", "T1 path export incl. vertices.f64", "candidate after archiving the validation", "frozen exporter + t_time_collapse_export.py --mode primary"),
    ("muon_*/S0", "frozen replay of T1 (reference waveforms)", "KEEP for validated muons", "p3_replay_full_spectral_muon.py on T1_export (16-275 s per muon measured)"),
    ("muon_*/[GSK]*_*", "S-ladder variant waveforms", "candidate", "s_synthesize.py / s_kernel.py apply from green_* / kernel_* (seconds)"),
    ("muon_*/births_*", "birth-path validation", "KEEP (frozen evidence)", "p3_birth_moments + s_kernel.py apply"),
    ("muon_*/transport_validation_*", "grid-tool vs exporter check", "KEEP (small)", "run_s_ladder.sh validate"),
    ("report_*", "comparison reports", "KEEP", "run_s_ladder.sh report"),
    ("closure_*", "closure comparisons", "KEEP", "run_s_ladder.sh closure-compare"),
]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run", type=Path); ap.add_argument("--out", type=Path)
    a = ap.parse_args()
    L = [f"# Cleanup inventory of `{a.run}` (nothing deleted)", "",
         "Never candidates: the archives (`ion_births.root`, census CSVs), the sealed 20-muon sample, frozen sources.", "",
         "| path | category | size [MiB] | status | reconstruction |", "|---|---|---|---|---|"]
    tot = {}
    for g, cat, keep, rec in RULES:
        for p in sorted(a.run.glob(g)):
            if not g.endswith(".npy") and not p.is_dir():
                continue
            s = size(p); tot[keep] = tot.get(keep, 0) + s
            L.append(f"| {p.relative_to(a.run)} | {cat} | {s / 2**20:.1f} | {keep} | {rec} |")
    L += ["", "Totals by status: " + ", ".join(f"{k}: {v / 2**30:.2f} GiB" for k, v in tot.items())]
    txt = "\n".join(L) + "\n"
    if a.out:
        a.out.write_text(txt)
    print(txt)


if __name__ == "__main__":
    main()
