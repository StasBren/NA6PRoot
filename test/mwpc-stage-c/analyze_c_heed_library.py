#!/usr/bin/env python3
"""Summarise a Heed ionisation library and quantify the Poisson approximation."""

from __future__ import annotations

import argparse
import csv
import math
import statistics
from collections import Counter
from pathlib import Path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tracks", required=True, type=Path)
    ap.add_argument("--clusters", required=True, type=Path)
    ap.add_argument("--reference-length-mm", required=True, type=float)
    ap.add_argument("--output", required=True, type=Path)
    args = ap.parse_args()

    with args.tracks.open(newline="", encoding="utf-8") as f:
        tracks = list(csv.DictReader(f))
    with args.clusters.open(newline="", encoding="utf-8") as f:
        clusters = list(csv.DictReader(f))

    ncl = [int(float(r["clusters"])) for r in tracks]
    ne_track = [int(float(r["electrons"])) for r in tracks]
    ne_cluster = [int(float(r["electrons"])) for r in clusters]
    v = [float(r["v_mm"]) for r in clusters]
    u = [float(r["u_mm"]) for r in clusters]

    mean_ncl = statistics.fmean(ncl)
    var_ncl = statistics.variance(ncl) if len(ncl) > 1 else 0.0
    fano_like = var_ncl / mean_ncl if mean_ncl > 0 else math.nan
    lam = mean_ncl / args.reference_length_mm

    mean_ne_track = statistics.fmean(ne_track)
    mean_ne_cluster = statistics.fmean(ne_cluster)

    # Pearson chi2 diagnostic against Poisson(lambda=sample mean), grouping the
    # far tail until each expected bin is reasonably populated.
    obs = Counter(ncl)
    n = len(ncl)
    kmax = max(ncl)
    expected = {}
    p = math.exp(-mean_ncl)
    expected[0] = n * p
    for k in range(1, kmax + 1):
        p *= mean_ncl / k
        expected[k] = n * p

    bins = []
    run_o = 0.0
    run_e = 0.0
    for k in range(0, kmax + 1):
        run_o += obs.get(k, 0)
        run_e += expected.get(k, 0.0)
        if run_e >= 5.0:
            bins.append((run_o, run_e))
            run_o = 0.0
            run_e = 0.0
    if run_e > 0:
        if bins:
            o0, e0 = bins[-1]
            bins[-1] = (o0 + run_o, e0 + run_e)
        else:
            bins.append((run_o, run_e))

    chi2 = sum((o - e) ** 2 / e for o, e in bins if e > 0)
    ndof = max(len(bins) - 2, 0)  # one fitted mean + normalization

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as f:
        f.write("HEED IONISATION LIBRARY SUMMARY\n")
        f.write(f"tracks = {len(tracks)}\n")
        f.write(f"clusters = {len(clusters)}\n")
        f.write(f"reference_length_mm = {args.reference_length_mm:.9g}\n")
        f.write(f"mean_clusters_per_track = {mean_ncl:.9g}\n")
        f.write(f"variance_clusters_per_track = {var_ncl:.9g}\n")
        f.write(f"variance_over_mean = {fano_like:.9g}\n")
        f.write(f"lambda_clusters_per_mm = {lam:.9g}\n")
        f.write(f"mean_electrons_per_track = {mean_ne_track:.9g}\n")
        f.write(f"mean_electrons_per_cluster = {mean_ne_cluster:.9g}\n")
        f.write(f"v_min_mm = {min(v):.9g}\n")
        f.write(f"v_max_mm = {max(v):.9g}\n")
        f.write(f"u_min_mm = {min(u):.9g}\n")
        f.write(f"u_max_mm = {max(u):.9g}\n")
        f.write(f"poisson_chi2 = {chi2:.9g}\n")
        f.write(f"poisson_ndof = {ndof}\n")
        if ndof > 0:
            f.write(f"poisson_chi2_per_ndof = {chi2 / ndof:.9g}\n")

    print(args.output.read_text())


if __name__ == "__main__":
    main()
