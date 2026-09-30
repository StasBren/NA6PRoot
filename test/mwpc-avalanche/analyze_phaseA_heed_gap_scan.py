#!/usr/bin/env python3
"""Compare Heed conduction-electron statistics for 5 mm and 6 mm gas gaps.

Input files come from mwpc_phase_a_heed_stats (one CSV row per muon).
No detector efficiency is inferred here: probabilities below a primary-electron
cut are only ionisation-level proxies, NOT electronic-threshold efficiencies.
"""
import argparse
import bisect
import csv
import math
import statistics
from pathlib import Path


def quantile(sorted_values, fraction):
    if not sorted_values:
        raise ValueError("Quantile requested for an empty sample")
    p = fraction * (len(sorted_values) - 1)
    low, high = math.floor(p), math.ceil(p)
    return sorted_values[low] + (p - low) * (
        sorted_values[high] - sorted_values[low]
    )


def load_tracks(path):
    path = Path(path)
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        expected = {"clusters", "electrons", "total_dE_eV"}
        if not reader.fieldnames or not expected.issubset(reader.fieldnames):
            raise ValueError("Missing required Heed columns in " + str(path))
        rows = []
        for row in reader:
            n_clusters = int(row["clusters"])
            n_electrons = int(row["electrons"])
            energy = float(row["total_dE_eV"])
            if n_clusters < 0 or n_electrons < 0 or not math.isfinite(energy):
                raise ValueError("Invalid Heed event in " + str(path))
            rows.append((n_clusters, n_electrons, energy))
    if len(rows) < 2:
        raise ValueError("At least two tracks required: " + str(path))
    return rows


def summarise(gap_mm, rows):
    clusters = [row[0] for row in rows]
    electrons = sorted(row[1] for row in rows)
    energies = [row[2] for row in rows]
    mean_n = statistics.fmean(electrons)
    var_n = statistics.variance(electrons)
    return {
        "gap_mm": gap_mm,
        "tracks": len(rows),
        "mean_clusters": statistics.fmean(clusters),
        "mean_electrons": mean_n,
        "variance_electrons": var_n,
        "sd_electrons": math.sqrt(var_n),
        "cv_electrons": math.sqrt(var_n) / mean_n if mean_n else math.nan,
        "median_electrons": statistics.median(electrons),
        "p95_electrons": quantile(electrons, 0.95),
        "p99_electrons": quantile(electrons, 0.99),
        "max_electrons": max(electrons),
        "mean_dE_eV": statistics.fmean(energies),
        "electrons_sorted": electrons,
    }


def probability_below(sorted_values, threshold):
    return bisect.bisect_left(sorted_values, threshold) / len(sorted_values)


def write_csv(path, fieldnames, rows):
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row[key] for key in fieldnames})


def make_plots(a, b, outdir):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib unavailable; CSV summaries still written. Install matplotlib for plots.")
        return

    samples = [a, b]
    cap = max(30, int(math.ceil(max(s["p99_electrons"] for s in samples))))
    fig, ax = plt.subplots(figsize=(8, 5))
    for s in samples:
        x = [min(v, cap) for v in s["electrons_sorted"]]
        ax.hist(x, bins=50, range=(0, cap), histtype="step", density=True,
                label=str(s["gap_mm"]) + " mm")
    ax.set_xlabel("Conduction electrons per muon track (99th percentile shown)")
    ax.set_ylabel("Density (overflow in rightmost bin)")
    ax.set_title("Heed primary-electron count")
    ax.legend()
    fig.tight_layout()
    fig.savefig(outdir / "electron_count_histogram.png", dpi=170)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 5))
    xmax = max(100, int(math.ceil(max(s["p95_electrons"] for s in samples))))
    thresholds = list(range(xmax + 1))
    for s in samples:
        curve = [probability_below(s["electrons_sorted"], t) for t in thresholds]
        ax.plot(thresholds, curve, label=str(s["gap_mm"]) + " mm")
    ax.set_xlabel("Primary-electron threshold (N < threshold)")
    ax.set_ylabel("Fraction of muons below threshold")
    ax.set_ylim(0, 1)
    ax.set_title("Lower tail of ionisation distribution")
    ax.legend()
    fig.tight_layout()
    fig.savefig(outdir / "low_electron_tail.png", dpi=170)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot([5, 6], [a["cv_electrons"], b["cv_electrons"]], marker="o",
            label="Heed simulation")
    predicted_cv_6 = a["cv_electrons"] * math.sqrt(4.9 / 5.9)
    ax.plot([5, 6], [a["cv_electrons"], predicted_cv_6], linestyle="--",
            marker="x", label="1/sqrt(effective path length)")
    ax.set_xticks([5, 6])
    ax.set_xlabel("Nominal gas gap (mm)")
    ax.set_ylabel("SD(N) / mean(N)")
    ax.set_title("Relative primary-ionisation fluctuations")
    ax.legend()
    fig.tight_layout()
    fig.savefig(outdir / "relative_fluctuations.png", dpi=170)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--five", required=True, help="5 mm *_tracks.csv")
    parser.add_argument("--six", required=True, help="6 mm *_tracks.csv")
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()
    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    a = summarise(5, load_tracks(args.five))
    b = summarise(6, load_tracks(args.six))
    keys = [
        "gap_mm", "tracks", "mean_clusters", "mean_electrons",
        "variance_electrons", "sd_electrons", "cv_electrons",
        "median_electrons", "p95_electrons", "p99_electrons",
        "max_electrons", "mean_dE_eV"
    ]
    write_csv(outdir / "summary.csv", keys, [a, b])

    thresholds = [1, 5, 10, 20, 30, 40, 50, 75, 100, 150, 200]
    tail_rows = []
    for threshold in thresholds:
        tail_rows.append({
            "electron_threshold_exclusive": threshold,
            "p_5mm_N_below": probability_below(a["electrons_sorted"], threshold),
            "p_6mm_N_below": probability_below(b["electrons_sorted"], threshold),
        })
    write_csv(outdir / "low_tail.csv",
              ["electron_threshold_exclusive", "p_5mm_N_below", "p_6mm_N_below"],
              tail_rows)

    ratio_mean = b["mean_electrons"] / a["mean_electrons"]
    ratio_cv = b["cv_electrons"] / a["cv_electrons"]
    expected_mean_ratio = 5.9 / 4.9
    expected_cv_ratio = math.sqrt(4.9 / 5.9)

    report = (
        "HEED-ONLY MWPC GAS-GAP SCAN\n"
        "===========================\n"
        "10 GeV/c mu-, Ar/CO2 70:30, normal incidence, symmetric cathode gaps.\n"
        "Fixed starting offset: approximately 0.1 mm inside +v cathode.\n"
        "Nominal gas gaps: 5 and 6 mm. Effective straight paths: about 4.9 and 5.9 mm.\n\n"
        "                                 5 mm                6 mm\n"
        f"Tracks                       {a['tracks']:12d}        {b['tracks']:12d}\n"
        f"Mean electrons               {a['mean_electrons']:12.4f}        {b['mean_electrons']:12.4f}\n"
        f"SD electrons                 {a['sd_electrons']:12.4f}        {b['sd_electrons']:12.4f}\n"
        f"Relative SD (CV)             {a['cv_electrons']:12.4f}        {b['cv_electrons']:12.4f}\n"
        f"Mean clusters                {a['mean_clusters']:12.4f}        {b['mean_clusters']:12.4f}\n"
        f"99th percentile electrons    {a['p99_electrons']:12.1f}        {b['p99_electrons']:12.1f}\n\n"
        f"Observed mean N ratio 6/5: {ratio_mean:.4f}; path-length model: {expected_mean_ratio:.4f}\n"
        f"Observed CV ratio 6/5:     {ratio_cv:.4f}; 1/sqrt(L) model: {expected_cv_ratio:.4f}\n\n"
        "IMPORTANT:\n"
        "- This is primary ionisation only: no Garfield drift, avalanche, strips or electronics.\n"
        "- P(N < n) is NOT a detector inefficiency or a physical charge threshold.\n"
        "- Large delta-electron events make sample variances sensitive to statistics.\n"
        "- Changing drift length, electrode separation and ionisation length have distinct effects.\n"
        "- This executable uses its own idealised gas cell, not Geant4 input crossings.\n"
    )
    (outdir / "comparison.txt").write_text(report)
    print(report)
    make_plots(a, b, outdir)


if __name__ == "__main__":
    main()
