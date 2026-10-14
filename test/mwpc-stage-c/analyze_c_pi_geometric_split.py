#!/usr/bin/env python3
"""
Analyze the PI geometric Poisson/binomial wire-sharing output.

Inputs are produced by c_pi_geant_wire_split.py:
  PREFIX_crossings.csv
  PREFIX_wires.csv

The analysis intentionally stays at the simple PI-model level.  It checks:
  * how often a GEANT track crosses more than one geometric wire cell;
  * how often stochastic ionisation actually populates more than one wire;
  * mean path length / Poisson cluster count / Heed electron multiplicity;
  * agreement between geometric path fractions and sampled cluster/electron
    fractions for two-cell crossings;
  * dependence of sharing on the fractional path in the second wire.

No Garfield drift or electronics response enters here.
"""

from __future__ import annotations

import argparse
import csv
import math
import statistics as stats
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np


def read_csv(path: Path):
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def fmean(xs):
    xs = [float(x) for x in xs if math.isfinite(float(x))]
    return stats.fmean(xs) if xs else math.nan


def rms(xs):
    xs = np.asarray([float(x) for x in xs if math.isfinite(float(x))], dtype=float)
    return float(np.sqrt(np.mean(xs * xs))) if len(xs) else math.nan


def quantile(xs, q):
    xs = np.asarray([float(x) for x in xs if math.isfinite(float(x))], dtype=float)
    return float(np.quantile(xs, q)) if len(xs) else math.nan


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--prefix", required=True, type=Path)
    ap.add_argument("--output-dir", type=Path)
    args = ap.parse_args()

    prefix = args.prefix
    crossings_path = Path(str(prefix) + "_crossings.csv")
    wires_path = Path(str(prefix) + "_wires.csv")
    out = args.output_dir or prefix.parent
    out.mkdir(parents=True, exist_ok=True)

    crossings = read_csv(crossings_path)
    wires = read_csv(wires_path)
    if not crossings or not wires:
        raise SystemExit("Missing/empty geometric split inputs.")

    by_crossing = defaultdict(list)
    for r in wires:
        by_crossing[int(r["crossing_id"])].append(r)

    n = len(crossings)
    geom_multi = [r for r in crossings if int(r["wire_cells_crossed"]) > 1]
    stoch_multi = [r for r in crossings if int(r["multiwire_cluster_event"]) > 0]

    path_lengths = np.asarray([float(r["path_length_mm"]) for r in crossings])
    nclusters = np.asarray([int(r["n_clusters_total"]) for r in crossings])
    nelectrons = np.asarray([int(r["n_primary_electrons_total"]) for r in crossings])

    cell_counts = Counter(int(r["wire_cells_crossed"]) for r in crossings)

    two = []
    for r in geom_multi:
        cid = int(r["crossing_id"])
        wr = sorted(by_crossing[cid], key=lambda x: float(x["segment_length_mm"]), reverse=True)
        if len(wr) != 2:
            continue

        # Define "minor" as the geometrically shorter side.  This makes the
        # observable independent of which physical wire is labelled first.
        minor = wr[1]
        p = float(minor["path_fraction"])
        cfrac = float(minor["cluster_fraction"])
        efrac = float(minor["electron_fraction"])
        ncl = int(r["n_clusters_total"])
        ne = int(r["n_primary_electrons_total"])
        two.append((p, cfrac, efrac, ncl, ne))

    p = np.asarray([x[0] for x in two], dtype=float)
    cf = np.asarray([x[1] for x in two], dtype=float)
    ef = np.asarray([x[2] for x in two], dtype=float)
    finite_c = np.isfinite(cf)
    finite_e = np.isfinite(ef)

    # Binomial standardised residual for the cluster count.  For each event
    # sigma_f = sqrt(p(1-p)/N).  A correct multinomial implementation should
    # give z with width near one (away from p=0/1 and N=0).
    z_cluster = []
    for pp, cc, _, nn, _ in two:
        if nn <= 0 or not math.isfinite(cc) or pp <= 0 or pp >= 1:
            continue
        sig = math.sqrt(pp * (1.0 - pp) / nn)
        if sig > 0:
            z_cluster.append((cc - pp) / sig)

    # Bin p_minor to show the expected mean sharing curve.
    bins = np.linspace(0.0, 0.5, 11)
    binned = []
    for lo, hi in zip(bins[:-1], bins[1:]):
        mask = (p >= lo) & (p < hi if hi < 0.5 else p <= hi)
        if not np.any(mask):
            continue
        cvals = cf[mask & finite_c]
        evals = ef[mask & finite_e]
        binned.append({
            "p_lo": lo,
            "p_hi": hi,
            "p_mean": float(np.mean(p[mask])),
            "events": int(np.sum(mask)),
            "cluster_fraction_mean": float(np.mean(cvals)) if len(cvals) else math.nan,
            "cluster_fraction_std": float(np.std(cvals, ddof=1)) if len(cvals) > 1 else math.nan,
            "electron_fraction_mean": float(np.mean(evals)) if len(evals) else math.nan,
            "electron_fraction_std": float(np.std(evals, ddof=1)) if len(evals) > 1 else math.nan,
        })

    summary = {
        "crossings": n,
        "geometric_multi_cell": len(geom_multi),
        "geometric_multi_cell_fraction": len(geom_multi) / n,
        "stochastic_multi_wire": len(stoch_multi),
        "stochastic_multi_wire_fraction": len(stoch_multi) / n,
        "stochastic_given_geometric": len(stoch_multi) / len(geom_multi) if geom_multi else math.nan,
        "mean_path_length_mm": float(np.mean(path_lengths)),
        "median_path_length_mm": float(np.median(path_lengths)),
        "q05_path_length_mm": float(np.quantile(path_lengths, 0.05)),
        "q95_path_length_mm": float(np.quantile(path_lengths, 0.95)),
        "mean_clusters_per_crossing": float(np.mean(nclusters)),
        "variance_clusters_per_crossing": float(np.var(nclusters, ddof=1)),
        "mean_primary_electrons_per_crossing": float(np.mean(nelectrons)),
        "mean_primary_electrons_per_cluster": float(np.sum(nelectrons) / np.sum(nclusters)),
        "two_cell_crossings": len(two),
        "rms_cluster_fraction_minus_path_fraction": rms(cf[finite_c] - p[finite_c]),
        "rms_electron_fraction_minus_path_fraction": rms(ef[finite_e] - p[finite_e]),
        "mean_cluster_fraction_minus_path_fraction": fmean(cf[finite_c] - p[finite_c]),
        "mean_electron_fraction_minus_path_fraction": fmean(ef[finite_e] - p[finite_e]),
        "cluster_binomial_z_mean": fmean(z_cluster),
        "cluster_binomial_z_std": stats.stdev(z_cluster) if len(z_cluster) > 1 else math.nan,
        "wire_cells_crossed_counts": dict(sorted(cell_counts.items())),
    }

    summary_path = out / "c_pi_geometric_split_summary.txt"
    with summary_path.open("w", encoding="utf-8") as f:
        for k, v in summary.items():
            f.write(f"{k}={v}\n")

    binned_path = out / "c_pi_geometric_split_binned.csv"
    with binned_path.open("w", newline="", encoding="utf-8") as f:
        fields = list(binned[0]) if binned else []
        w = csv.DictWriter(f, fieldnames=fields)
        if fields:
            w.writeheader()
            w.writerows(binned)

    print("PI geometric split analysis")
    print(f"  crossings                         : {n}")
    print(f"  wire cells crossed                : {dict(sorted(cell_counts.items()))}")
    print(f"  geometric multi-cell fraction     : {len(geom_multi)/n:.6f}")
    print(f"  stochastic multi-wire fraction    : {len(stoch_multi)/n:.6f}")
    print(f"  P(stochastic multi | geom multi)  : {len(stoch_multi)/len(geom_multi):.6f}")
    print(f"  mean path length [mm]             : {np.mean(path_lengths):.6f}")
    print(f"  mean clusters / crossing          : {np.mean(nclusters):.6f}")
    print(f"  mean primary e- / crossing        : {np.mean(nelectrons):.6f}")
    print(f"  mean primary e- / cluster         : {np.sum(nelectrons)/np.sum(nclusters):.6f}")
    print(f"  two-cell crossings                : {len(two)}")
    print(f"  RMS(cluster fraction - path frac) : {summary['rms_cluster_fraction_minus_path_fraction']:.6f}")
    print(f"  RMS(electron fraction - path frac): {summary['rms_electron_fraction_minus_path_fraction']:.6f}")
    print(f"  binomial z mean/std               : {summary['cluster_binomial_z_mean']:.6f} / {summary['cluster_binomial_z_std']:.6f}")
    print(f"  summary                           : {summary_path}")
    print(f"  binned                            : {binned_path}")

    try:
        import matplotlib.pyplot as plt

        # Plot 1: sampled sharing vs geometric path fraction.
        fig, ax = plt.subplots(figsize=(7.2, 5.4))
        if len(two):
            take = np.arange(len(two))
            if len(take) > 12000:
                rng = np.random.default_rng(20261014)
                take = rng.choice(take, 12000, replace=False)
            ax.scatter(p[take], ef[take], s=7, alpha=0.16, label="event electron fraction")
            ax.plot([0, 0.5], [0, 0.5], linewidth=2, label="geometric expectation")
            if binned:
                ax.plot(
                    [r["p_mean"] for r in binned],
                    [r["electron_fraction_mean"] for r in binned],
                    marker="o",
                    linewidth=2,
                    label="binned mean",
                )
        ax.set_xlabel("Minor-wire geometric path fraction")
        ax.set_ylabel("Minor-wire primary-electron fraction")
        ax.set_xlim(0, 0.5)
        ax.set_ylim(0, 1)
        ax.legend()
        ax.set_title("PI geometric wire sharing: path length → sampled electrons")
        fig.tight_layout()
        p1 = out / "c_pi_geometric_sharing.png"
        fig.savefig(p1, dpi=180)
        plt.close(fig)

        # Plot 2: path-length distribution.
        fig, ax = plt.subplots(figsize=(7.2, 5.0))
        ax.hist(path_lengths, bins=80)
        ax.axvline(6.0, linestyle="--", linewidth=2, label="normal-incidence 6 mm")
        ax.set_xlabel("Muon path length in nominal gas segment [mm]")
        ax.set_ylabel("Crossings")
        ax.set_title("Realistic GEANT angles lengthen the 6-mm gas path")
        ax.legend()
        fig.tight_layout()
        p2 = out / "c_pi_geometric_path_lengths.png"
        fig.savefig(p2, dpi=180)
        plt.close(fig)

        print(f"  plot                              : {p1}")
        print(f"  plot                              : {p2}")
    except ImportError:
        print("  matplotlib not available: numerical outputs were still written.")


if __name__ == "__main__":
    main()
