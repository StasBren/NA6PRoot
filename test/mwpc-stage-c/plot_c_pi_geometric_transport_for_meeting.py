#!/usr/bin/env python3
"""
Create PI-facing plots for the geometric wire-sharing + Garfield-box study.

Required inputs:
  --segments
      c_pi_nominal_segments.csv
  --geometric-crossings
      c_pi_wire_split_crossings.csv
  --geometric-wires
      c_pi_wire_split_wires.csv
  --garfield-crossings
      c_pi_wire_transport_crossings.csv
  --garfield-clusters
      c_pi_wire_transport_clusters.csv
  --garfield-wires
      c_pi_wire_transport_wires.csv

The figures are deliberately separated: each plot has one main message.
"""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt


def read(path):
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def ff(r, k):
    try:
        return float(r[k])
    except (KeyError, TypeError, ValueError):
        return math.nan


def ii(r, k):
    return int(float(r[k]))


def save(fig, path):
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)
    print("  wrote:", path)


def binned_ratio(x, num, den, edges):
    x = np.asarray(x, dtype=float)
    num = np.asarray(num, dtype=float)
    den = np.asarray(den, dtype=float)
    rows = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        mask = (x >= lo) & (x < hi if hi < edges[-1] else x <= hi)
        if not np.any(mask):
            continue
        d = np.sum(den[mask])
        rows.append((
            float(np.mean(x[mask])),
            float(np.sum(num[mask]) / d) if d > 0 else math.nan,
            int(np.sum(mask)),
        ))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--segments", required=True, type=Path)
    ap.add_argument("--geometric-crossings", required=True, type=Path)
    ap.add_argument("--geometric-wires", required=True, type=Path)
    ap.add_argument("--garfield-crossings", required=True, type=Path)
    ap.add_argument("--garfield-clusters", required=True, type=Path)
    ap.add_argument("--garfield-wires", required=True, type=Path)
    ap.add_argument("--output-dir", required=True, type=Path)
    args = ap.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    seg = read(args.segments)
    gc = read(args.geometric_crossings)
    gw = read(args.geometric_wires)
    tc = read(args.garfield_crossings)
    tcl = read(args.garfield_clusters)
    tw = read(args.garfield_wires)

    # ------------------------------------------------------------------
    # 1-2. Real GEANT angular distributions.
    # ------------------------------------------------------------------
    theta_u = np.degrees(np.arctan2(
        [ff(r, "py_in_GeV") for r in seg],
        [ff(r, "pz_in_GeV") for r in seg],
    ))
    theta_w = np.degrees(np.arctan2(
        [ff(r, "px_in_GeV") for r in seg],
        [ff(r, "pz_in_GeV") for r in seg],
    ))

    def angle_plot(values, symbol, title, name):
        q05, q95 = np.quantile(values, [0.05, 0.95])
        fig, ax = plt.subplots(figsize=(7.2, 5.0))
        ax.hist(values, bins=120)
        ax.axvline(q05, linestyle="--", linewidth=1.5)
        ax.axvline(q95, linestyle="--", linewidth=1.5)
        ax.set_xlabel(f"{symbol} [deg]")
        ax.set_ylabel("Tracks")
        ax.set_title(title)
        ax.text(
            0.98, 0.96,
            f"mean = {np.mean(values):.2f} deg\n"
            f"std = {np.std(values, ddof=1):.2f} deg\n"
            f"5-95% = [{q05:.2f}, {q95:.2f}] deg",
            transform=ax.transAxes,
            ha="right", va="top",
        )
        save(fig, args.output_dir / name)

    angle_plot(
        theta_u,
        r"$\theta_u = \arctan(p_y/p_z)$",
        "Real GEANT track angles across the anode wires",
        "pi_01_theta_u_distribution.png",
    )
    angle_plot(
        theta_w,
        r"$\theta_w = \arctan(p_x/p_z)$",
        "Real GEANT track angles along the anode wires",
        "pi_02_theta_w_distribution.png",
    )

    # ------------------------------------------------------------------
    # 3. Gas path length: zoom where the population actually is.
    # ------------------------------------------------------------------
    path = np.asarray([ff(r, "path_length_mm") for r in seg])
    q99 = float(np.quantile(path, 0.99))
    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    shown = path[path <= max(6.8, q99)]
    ax.hist(shown, bins=80)
    ax.axvline(6.0, linestyle="--", linewidth=1.5, label="normal incidence")
    ax.set_xlabel("Muon path length in 6-mm gas [mm]")
    ax.set_ylabel("Tracks")
    ax.set_title("GEANT angles slightly lengthen the gas path")
    ax.legend()
    ax.text(
        0.98, 0.78,
        f"mean = {np.mean(path):.3f} mm\n"
        f"median = {np.median(path):.3f} mm\n"
        f"99% < {q99:.3f} mm\n"
        f"max = {np.max(path):.2f} mm",
        transform=ax.transAxes,
        ha="right", va="top",
    )
    save(fig, args.output_dir / "pi_03_path_length_zoom.png")

    # Index geometric/transport wire rows by crossing.
    gw_by = defaultdict(list)
    for r in gw:
        gw_by[ii(r, "crossing_id")].append(r)
    tw_by = defaultdict(dict)
    for r in tw:
        tw_by[ii(r, "crossing_id")][ii(r, "wire_index")] = r
    tc_by = {ii(r, "crossing_id"): r for r in tc}
    gc_by = {ii(r, "crossing_id"): r for r in gc}

    # Two-cell events and the geometrically shorter ("minor") side.
    minor_rows = []
    for cid, rows in gw_by.items():
        if len(rows) != 2 or cid not in tc_by:
            continue
        rows = sorted(rows, key=lambda r: ff(r, "segment_length_mm"))
        minor = rows[0]
        wire = ii(minor, "wire_index")
        p = ff(minor, "path_fraction")
        efrac = ff(minor, "electron_fraction")
        after = tw_by.get(cid, {}).get(wire)
        gfrac = ff(after, "garfield_collected_fraction") if after else 0.0
        minor_rows.append((cid, p, efrac, gfrac))

    pm = np.asarray([x[1] for x in minor_rows])
    em = np.asarray([x[2] for x in minor_rows])
    gm = np.asarray([x[3] for x in minor_rows])

    # ------------------------------------------------------------------
    # 4. Geometric expectation -> Heed electrons -> Garfield collected.
    # ------------------------------------------------------------------
    edges = np.linspace(0.0, 0.5, 11)
    centers, ep, gp, counts = [], [], [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        mask = (pm >= lo) & (pm < hi if hi < 0.5 else pm <= hi)
        if not np.any(mask):
            continue
        centers.append(float(np.mean(pm[mask])))
        ep.append(float(np.nanmean(em[mask])))
        gp.append(float(np.nanmean(gm[mask])))
        counts.append(int(np.sum(mask)))

    fig, ax = plt.subplots(figsize=(7.2, 5.2))
    ax.plot([0, 0.5], [0, 0.5], linestyle="--", linewidth=2,
            label="geometric path fraction")
    ax.plot(centers, ep, marker="o", linewidth=2,
            label="Poisson + Heed electrons")
    ax.plot(centers, gp, marker="s", linewidth=2,
            label="after Garfield transport")
    ax.set_xlim(0, 0.5)
    ax.set_ylim(0, 0.6)
    ax.set_xlabel("Minor-wire geometric path fraction")
    ax.set_ylabel("Mean fraction assigned/collected on that wire")
    ax.set_title("From midpoint geometry to microscopic electron transport")
    ax.legend()
    save(fig, args.output_dir / "pi_04_geometric_vs_garfield_sharing.png")

    # ------------------------------------------------------------------
    # 5. Probability that both wires are populated, before and after drift.
    # ------------------------------------------------------------------
    pvals, pre_both, post_both = [], [], []
    for cid, p, _, _ in minor_rows:
        pvals.append(p)
        pre_both.append(int(ii(gc_by[cid], "multiwire_cluster_event") > 0))
        post_both.append(int(ii(tc_by[cid], "garfield_active_wires") > 1))

    pvals = np.asarray(pvals)
    pre_both = np.asarray(pre_both)
    post_both = np.asarray(post_both)

    xbin, premean, postmean = [], [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        mask = (pvals >= lo) & (pvals < hi if hi < 0.5 else pvals <= hi)
        if not np.any(mask):
            continue
        xbin.append(float(np.mean(pvals[mask])))
        premean.append(float(np.mean(pre_both[mask])))
        postmean.append(float(np.mean(post_both[mask])))

    fig, ax = plt.subplots(figsize=(7.2, 5.2))
    ax.plot(xbin, premean, marker="o", linewidth=2,
            label="Poisson + Heed, midpoint assignment")
    ax.plot(xbin, postmean, marker="s", linewidth=2,
            label="after Garfield transport")
    ax.set_xlabel("Minor-wire geometric path fraction")
    ax.set_ylabel("Probability of >1 active wire")
    ax.set_ylim(0, 1.05)
    ax.set_title("When does a boundary-crossing muon actually populate two wires?")
    ax.legend()
    save(fig, args.output_dir / "pi_05_two_wire_probability.png")

    # ------------------------------------------------------------------
    # 6. This is where the small (u,v) boxes enter: migration vs local u.
    # ------------------------------------------------------------------
    uc = np.asarray([ff(r, "u_cell_mm") for r in tcl])
    coll = np.asarray([ii(r, "collected_electrons") for r in tcl], dtype=float)
    migr = np.asarray([ii(r, "migrated_electrons") for r in tcl], dtype=float)

    uedges = np.linspace(-2.0, 2.0, 21)
    migration_bins = binned_ratio(uc, migr, coll, uedges)

    fig, ax = plt.subplots(figsize=(7.2, 5.2))
    ax.plot(
        [r[0] for r in migration_bins],
        [r[1] for r in migration_bins],
        marker="o", linewidth=2,
    )
    ax.set_xlabel("Cluster local u within geometrically assigned wire cell [mm]")
    ax.set_ylabel("Fraction of collected electrons ending on another wire")
    ax.set_ylim(bottom=0)
    ax.set_title("Garfield correction to the hard midpoint rule")
    save(fig, args.output_dir / "pi_06_garfield_migration_vs_u.png")

    # ------------------------------------------------------------------
    # 7. Timing information from the same small-box lookup.
    # ------------------------------------------------------------------
    absu = np.abs(uc)
    b100 = np.asarray(
        [ii(r, "collected_electrons") * ff(r, "fraction_collected_before_100ns")
         if math.isfinite(ff(r, "fraction_collected_before_100ns")) else 0.0
         for r in tcl],
        dtype=float,
    )
    b200 = np.asarray(
        [ii(r, "collected_electrons") * ff(r, "fraction_collected_before_200ns")
         if math.isfinite(ff(r, "fraction_collected_before_200ns")) else 0.0
         for r in tcl],
        dtype=float,
    )
    aedges = np.linspace(0.0, 2.0, 11)
    t100 = binned_ratio(absu, b100, coll, aedges)
    t200 = binned_ratio(absu, b200, coll, aedges)

    fig, ax = plt.subplots(figsize=(7.2, 5.2))
    ax.plot(
        [r[0] for r in t100], [r[1] for r in t100],
        marker="o", linewidth=2, label="arrived by 100 ns",
    )
    ax.plot(
        [r[0] for r in t200], [r[1] for r in t200],
        marker="s", linewidth=2, label="arrived by 200 ns",
    )
    ax.set_xlabel("|local u| within wire cell [mm]")
    ax.set_ylabel("Fraction of collected electrons")
    ax.set_ylim(0, 1.05)
    ax.set_title("Drift-time information retained by the Garfield box map")
    ax.legend()
    save(fig, args.output_dir / "pi_07_timing_100_200_vs_abs_u.png")

    # ------------------------------------------------------------------
    # 8. Event-level effect of Garfield on wire multiplicity.
    # ------------------------------------------------------------------
    geom_counts = np.asarray([ii(r, "geometric_active_wires") for r in tc])
    garf_counts = np.asarray([ii(r, "garfield_active_wires") for r in tc])
    categories = ["1 wire", "2+ wires"]
    geom_pct = [
        100 * np.mean(geom_counts == 1),
        100 * np.mean(geom_counts > 1),
    ]
    garf_pct = [
        100 * np.mean(garf_counts == 1),
        100 * np.mean(garf_counts > 1),
    ]
    x = np.arange(len(categories))
    width = 0.36
    fig, ax = plt.subplots(figsize=(7.0, 5.0))
    ax.bar(x - width / 2, geom_pct, width, label="midpoint + Heed")
    ax.bar(x + width / 2, garf_pct, width, label="after Garfield")
    ax.set_xticks(x, categories)
    ax.set_ylabel("Events [%]")
    ax.set_title("Microscopic wire occupancy before electronics threshold")
    ax.legend()
    save(fig, args.output_dir / "pi_08_event_wire_multiplicity.png")


    # ------------------------------------------------------------------
    # 9. Threshold sensitivity of the microscopic multi-wire fraction.
    # ------------------------------------------------------------------
    # "Active wire" above means >=1 collected primary electron.  A real VMM
    # channel has avalanche gain and a discriminator threshold, so the measured
    # multi-wire fraction cannot be identified with that microscopic occupancy.
    # Scan a simple proxy threshold in collected primary electrons per wire.
    thresholds = np.asarray([1, 2, 3, 5, 10, 15, 20], dtype=int)
    by_crossing_wire = defaultdict(list)
    for r in tw:
        by_crossing_wire[ii(r, "crossing_id")].append(
            ii(r, "garfield_collected_electrons")
        )

    threshold_multi = []
    for threshold in thresholds:
        active_multi = 0
        for r in tc:
            cid = ii(r, "crossing_id")
            n_active = sum(
                n >= threshold for n in by_crossing_wire.get(cid, [])
            )
            active_multi += int(n_active >= 2)
        threshold_multi.append(active_multi / len(tc))

    fig, ax = plt.subplots(figsize=(7.2, 5.2))
    ax.plot(
        thresholds,
        100.0 * np.asarray(threshold_multi),
        marker="o",
        linewidth=2,
    )
    ax.set_xlabel("Proxy threshold per wire [collected primary electrons]")
    ax.set_ylabel("Events with >=2 wires [%]")
    ax.set_ylim(bottom=0)
    ax.set_title("Multi-wire fraction is threshold dependent")
    ax.text(
        0.98, 0.96,
        "No avalanche gain or VMM discriminator model yet",
        transform=ax.transAxes,
        ha="right", va="top",
    )
    save(fig, args.output_dir / "pi_09_multiwire_threshold_scan.png")

    # Event-topology transition counts for presentation notes.
    transition_1_to_multi = int(np.sum((geom_counts == 1) & (garf_counts > 1)))
    transition_multi_to_1 = int(np.sum((geom_counts > 1) & (garf_counts == 1)))
    dominant_changed = int(np.sum(
        np.asarray([ii(r, "dominant_wire_changed") for r in tc]) > 0
    ))

    # Numerical summary for presentation notes.
    out_summary = args.output_dir / "pi_plot_summary.txt"
    with out_summary.open("w", encoding="utf-8") as f:
        f.write(f"tracks={len(seg)}\n")
        for name, values in (("theta_u", theta_u), ("theta_w", theta_w)):
            f.write(f"{name}_mean_deg={np.mean(values):.9g}\n")
            f.write(f"{name}_std_deg={np.std(values, ddof=1):.9g}\n")
            f.write(f"{name}_mean_abs_deg={np.mean(np.abs(values)):.9g}\n")
            f.write(f"{name}_q05_deg={np.quantile(values, .05):.9g}\n")
            f.write(f"{name}_q95_deg={np.quantile(values, .95):.9g}\n")
            f.write(f"{name}_abs_q95_deg={np.quantile(np.abs(values), .95):.9g}\n")
            f.write(f"{name}_abs_q99_deg={np.quantile(np.abs(values), .99):.9g}\n")
            f.write(f"{name}_min_deg={np.min(values):.9g}\n")
            f.write(f"{name}_max_deg={np.max(values):.9g}\n")
        f.write(f"path_mean_mm={np.mean(path):.9g}\n")
        f.write(f"path_q99_mm={q99:.9g}\n")
        f.write(f"path_max_mm={np.max(path):.9g}\n")
        f.write(f"two_cell_events={len(minor_rows)}\n")
        f.write(f"garfield_multiwire_fraction={np.mean(garf_counts > 1):.9g}\n")
        f.write(f"geometric_multiwire_fraction={np.mean(geom_counts > 1):.9g}\n")
        f.write(f"transition_1wire_to_multi={transition_1_to_multi}\n")
        f.write(f"transition_multi_to_1wire={transition_multi_to_1}\n")
        f.write(f"dominant_wire_changed={dominant_changed}\n")
        f.write(f"mean_event_collection_fraction={np.mean([ff(r, 'collection_fraction') for r in tc]):.9g}\n")
        f.write(f"mean_event_F100={np.mean([ff(r, 'fraction_collected_before_100ns') for r in tc]):.9g}\n")
        f.write(f"mean_event_F200={np.mean([ff(r, 'fraction_collected_before_200ns') for r in tc]):.9g}\n")
        for threshold, frac in zip(thresholds, threshold_multi):
            f.write(f"multiwire_fraction_threshold_ge_{threshold}e={frac:.9g}\n")
    print("  wrote:", out_summary)


if __name__ == "__main__":
    main()
