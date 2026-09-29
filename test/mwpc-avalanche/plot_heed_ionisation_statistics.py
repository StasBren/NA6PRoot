#!/usr/bin/env python3
"""
Plot Heed primary-ionisation statistics.

Expected inputs are produced by mwpc_phase_a_heed_stats:
  <prefix>_tracks.csv
  <prefix>_clusters.csv

Outputs:
  <output-prefix>_total_electrons_per_track.png
  <output-prefix>_cluster_size_distribution.png
  <output-prefix>_track_dE_vs_electrons.png
  <output-prefix>_cluster_dE_vs_size.png
  <output-prefix>_summary.txt
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def finite(a):
    a = np.asarray(a, dtype=float)
    return a[np.isfinite(a)]


def percentile(a, q):
    a = finite(a)
    return float(np.percentile(a, q)) if len(a) else float("nan")


def save_total_electrons(tracks, out):
    x = tracks["electrons"].to_numpy(dtype=float)
    n = len(x)

    xmax = max(1.0, percentile(x, 99.0))
    xmax = float(np.ceil(xmax))
    overflow = int(np.sum(x > xmax))

    nbins = min(60, max(15, int(np.sqrt(max(n, 1)))))
    edges = np.linspace(0.0, xmax, nbins + 1)
    weights = np.ones_like(x) / max(n, 1)

    mean = float(np.mean(x))
    median = float(np.median(x))
    q10, q90 = np.percentile(x, [10, 90])

    fig, ax = plt.subplots(figsize=(10.0, 6.2))
    ax.hist(x, bins=edges, weights=weights, alpha=0.82)
    ax.axvline(mean, linestyle="--", linewidth=1.8, label=f"mean = {mean:.1f}")
    ax.axvline(median, linestyle=":", linewidth=2.0, label=f"median = {median:.1f}")

    ax.set_xlabel("Total Heed conduction electrons per muon track")
    ax.set_ylabel("Fraction of tracks per bin")
    ax.set_title("Primary-ionisation electrons per muon track")
    ax.grid(alpha=0.18)
    ax.legend(loc="upper right")

    txt = (
        f"N tracks = {n}\n"
        f"10–90% interval = {q10:.0f}–{q90:.0f} e⁻\n"
        f"{overflow} tracks above plotted 99th-percentile range\n"
        f"maximum = {int(np.max(x))} e⁻"
    )
    ax.text(
        0.98, 0.72, txt,
        transform=ax.transAxes,
        ha="right", va="top",
        bbox=dict(boxstyle="round,pad=0.35", facecolor="white", alpha=0.93),
    )

    fig.tight_layout()
    fig.savefig(out, dpi=180)
    plt.close(fig)


def save_cluster_size(clusters, out):
    x = clusters["electrons"].to_numpy(dtype=int)
    n = len(x)

    q99 = max(1, int(np.ceil(np.percentile(x, 99.0))))
    cap = min(max(q99, 6), 30)

    values = np.arange(1, cap + 1)
    counts = np.array([np.sum(x == k) for k in values], dtype=float)
    fractions = counts / max(n, 1)
    overflow = int(np.sum(x > cap))

    fig, ax = plt.subplots(figsize=(10.0, 6.2))
    ax.bar(values, fractions, width=0.82)

    ax.set_xlabel("Conduction electrons in one Heed ionisation cluster")
    ax.set_ylabel("Fraction of all clusters")
    ax.set_title("Heed ionisation-cluster size distribution")
    ax.set_xticks(values)
    ax.grid(axis="y", alpha=0.18)

    txt = (
        f"N clusters = {n}\n"
        f"mean size = {np.mean(x):.2f} e⁻/cluster\n"
        f"median size = {np.median(x):.0f} e⁻/cluster\n"
        f"{overflow} clusters with > {cap} e⁻\n"
        f"maximum = {int(np.max(x))} e⁻"
    )
    ax.text(
        0.98, 0.96, txt,
        transform=ax.transAxes,
        ha="right", va="top",
        bbox=dict(boxstyle="round,pad=0.35", facecolor="white", alpha=0.93),
    )

    fig.tight_layout()
    fig.savefig(out, dpi=180)
    plt.close(fig)


def save_track_correlation(tracks, out):
    ne = tracks["electrons"].to_numpy(dtype=float)
    de = tracks["total_dE_eV"].to_numpy(dtype=float)

    ratio = np.sum(de) / np.sum(ne) if np.sum(ne) > 0 else float("nan")

    fig, ax = plt.subplots(figsize=(9.5, 6.4))
    ax.scatter(ne, de / 1000.0, s=13, alpha=0.35)

    if np.isfinite(ratio):
        xx = np.linspace(0.0, max(np.max(ne), 1.0), 200)
        ax.plot(
            xx,
            ratio * xx / 1000.0,
            linestyle="--",
            linewidth=1.6,
            label=f"global ΣΔE/ΣNₑ = {ratio:.1f} eV/e⁻",
        )

    corr = np.corrcoef(ne, de)[0, 1] if len(ne) > 1 else float("nan")

    ax.set_xlabel("Total conduction electrons per track")
    ax.set_ylabel("Total Heed energy transfer per track [keV]")
    ax.set_title("Track-by-track ionisation: energy transfer vs electron yield")
    ax.grid(alpha=0.18)
    ax.legend(loc="upper left")
    ax.text(
        0.98, 0.05,
        f"Pearson r = {corr:.3f}",
        transform=ax.transAxes,
        ha="right", va="bottom",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.9),
    )

    fig.tight_layout()
    fig.savefig(out, dpi=180)
    plt.close(fig)


def save_cluster_correlation(clusters, out):
    ne = clusters["electrons"].to_numpy(dtype=int)
    de = clusters["energy_transfer_eV"].to_numpy(dtype=float)

    # Large ensembles can have many tens of thousands of clusters.
    # A deterministic evenly spaced sample keeps the scatter readable.
    max_points = 30000
    if len(ne) > max_points:
        idx = np.linspace(0, len(ne) - 1, max_points).astype(int)
        ne_scatter = ne[idx]
        de_scatter = de[idx]
    else:
        ne_scatter = ne
        de_scatter = de

    fig, ax = plt.subplots(figsize=(9.5, 6.4))
    ax.scatter(ne_scatter, de_scatter, s=9, alpha=0.18)

    max_size_for_summary = min(20, int(np.percentile(ne, 99.5)))
    ks, med, lo, hi = [], [], [], []
    for k in range(1, max_size_for_summary + 1):
        vals = de[ne == k]
        if len(vals) < 5:
            continue
        ks.append(k)
        med.append(np.median(vals))
        lo.append(np.percentile(vals, 10))
        hi.append(np.percentile(vals, 90))

    if ks:
        ks = np.asarray(ks)
        med = np.asarray(med)
        lo = np.asarray(lo)
        hi = np.asarray(hi)
        ax.plot(ks, med, marker="o", linewidth=1.8, label="median ΔE at fixed cluster size")
        ax.fill_between(ks, lo, hi, alpha=0.16, label="10–90% interval")

    positive = de[de > 0]
    if len(positive) and np.max(positive) / np.min(positive) > 50:
        ax.set_yscale("log")

    ax.set_xlabel("Conduction electrons in cluster")
    ax.set_ylabel("Heed energy transfer to cluster [eV]")
    ax.set_title("Local energy transfer vs ionisation-cluster size")
    ax.grid(alpha=0.18)
    ax.legend(loc="upper left")

    fig.tight_layout()
    fig.savefig(out, dpi=180)
    plt.close(fig)


def write_summary(tracks, clusters, out):
    ne_track = tracks["electrons"].to_numpy(dtype=float)
    de_track = tracks["total_dE_eV"].to_numpy(dtype=float)
    ne_cluster = clusters["electrons"].to_numpy(dtype=float)
    de_cluster = clusters["energy_transfer_eV"].to_numpy(dtype=float)

    with open(out, "w", encoding="utf-8") as f:
        f.write("HEED PRIMARY-IONISATION STATISTICS\n")
        f.write(f"tracks = {len(tracks)}\n")
        f.write(f"clusters = {len(clusters)}\n")
        f.write(f"mean electrons/track = {np.mean(ne_track):.6g}\n")
        f.write(f"median electrons/track = {np.median(ne_track):.6g}\n")
        f.write(f"mean clusters/track = {np.mean(tracks['clusters']):.6g}\n")
        f.write(f"mean electrons/cluster = {np.mean(ne_cluster):.6g}\n")
        f.write(f"median electrons/cluster = {np.median(ne_cluster):.6g}\n")
        f.write(f"max electrons/cluster = {np.max(ne_cluster):.6g}\n")
        f.write(f"mean dE/track eV = {np.mean(de_track):.6g}\n")
        f.write(f"global dE/electron eV = {np.sum(de_track)/np.sum(ne_track):.6g}\n")
        f.write(f"track corr(dE, Ne) = {np.corrcoef(de_track, ne_track)[0,1]:.6g}\n")
        f.write(f"cluster corr(dE, Ne) = {np.corrcoef(de_cluster, ne_cluster)[0,1]:.6g}\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tracks-csv", required=True)
    ap.add_argument("--clusters-csv", required=True)
    ap.add_argument("--output-prefix", default="heed_ionisation")
    args = ap.parse_args()

    tracks = pd.read_csv(args.tracks_csv)
    clusters = pd.read_csv(args.clusters_csv)

    required_tracks = {"clusters", "electrons", "total_dE_eV"}
    required_clusters = {"electrons", "energy_transfer_eV"}

    missing_t = required_tracks - set(tracks.columns)
    missing_c = required_clusters - set(clusters.columns)
    if missing_t:
        raise RuntimeError(f"Missing track columns: {sorted(missing_t)}")
    if missing_c:
        raise RuntimeError(f"Missing cluster columns: {sorted(missing_c)}")

    prefix = Path(args.output_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)

    save_total_electrons(
        tracks, f"{prefix}_total_electrons_per_track.png"
    )
    save_cluster_size(
        clusters, f"{prefix}_cluster_size_distribution.png"
    )
    save_track_correlation(
        tracks, f"{prefix}_track_dE_vs_electrons.png"
    )
    save_cluster_correlation(
        clusters, f"{prefix}_cluster_dE_vs_size.png"
    )
    write_summary(
        tracks, clusters, f"{prefix}_summary.txt"
    )

    print("Wrote:")
    print(f"  {prefix}_total_electrons_per_track.png")
    print(f"  {prefix}_cluster_size_distribution.png")
    print(f"  {prefix}_track_dE_vs_electrons.png")
    print(f"  {prefix}_cluster_dE_vs_size.png")
    print(f"  {prefix}_summary.txt")


if __name__ == "__main__":
    main()
