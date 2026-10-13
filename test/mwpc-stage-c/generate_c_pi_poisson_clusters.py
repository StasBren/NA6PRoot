#!/usr/bin/env python3
"""
Generate a simple PI-style ionisation-cluster population for realistic MWPC
track crossings.

This is intentionally a first-order model:
  * number of clusters is Poisson(lambda_cl * track_length);
  * clusters are uniform in path length through the asymmetric 2+4 mm gas gap;
  * the wire-cell phase u0 at v=0 is sampled uniformly within one 4-mm pitch,
    unless --wire-phase-from-y is explicitly requested;
  * cluster electron multiplicity can be fixed, approximated with a shifted
    Poisson, or sampled directly from an empirical Heed cluster library.

The preferred cluster density is derived from a Heed library.  Note that
phaseA_heed_statistics starts 0.1 mm inside the +v cathode, so a nominal
2+4 mm chamber samples 5.9 mm of track.  The derived density is then applied
to the full physical path length in the 6-mm chamber.

The preferred multiplicity mode for the PI study is heed-empirical, using the
actual Heed 'electrons' distribution from --heed-clusters-csv.

This file does NOT use the Garfield response map yet.  It prepares the cluster
population that will later be folded with P(wire|u,v), drift time and arrival
angle from c_transport_box_summary.csv / raw electron samples.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import numpy as np


def fold_cell(u_mm: float, pitch_mm: float) -> tuple[int, float]:
    # nearest-wire convention -> local coordinate in [-pitch/2, +pitch/2)
    wire = math.floor(u_mm / pitch_mm + 0.5)
    local = u_mm - wire * pitch_mm
    # protect exact floating endpoint
    if local >= 0.5 * pitch_mm:
        wire += 1
        local -= pitch_mm
    if local < -0.5 * pitch_mm:
        wire -= 1
        local += pitch_mm
    return wire, local


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tracks", required=True, type=Path)
    ap.add_argument("--output", required=True, type=Path)
    ap.add_argument("--seed", type=int, default=180101)
    ap.add_argument("--gap-minus-mm", type=float, default=2.0)
    ap.add_argument("--gap-plus-mm", type=float, default=4.0)
    ap.add_argument("--wire-pitch-mm", type=float, default=4.0)
    ap.add_argument(
        "--lambda-clusters-per-mm",
        type=float,
        default=21.77 / 6.0,
        help=(
            "Poisson primary-cluster density along track length. "
            "Ignored when --heed-tracks-csv is supplied."
        ),
    )
    ap.add_argument(
        "--cluster-electron-model",
        choices=("one", "shifted-poisson", "heed-empirical"),
        default="heed-empirical",
    )
    ap.add_argument(
        "--mean-electrons-per-cluster",
        type=float,
        default=2.3704,
    )
    ap.add_argument(
        "--heed-tracks-csv",
        type=Path,
        help=(
            "Optional phaseA_heed_statistics *_tracks.csv. "
            "Mean clusters/reference length defines lambda_cl."
        ),
    )
    ap.add_argument(
        "--heed-clusters-csv",
        type=Path,
        help=(
            "Optional phaseA_heed_statistics *_clusters.csv used by "
            "--cluster-electron-model heed-empirical."
        ),
    )
    ap.add_argument(
        "--heed-reference-length-mm",
        type=float,
        default=5.9,
        help=(
            "Effective path length sampled by phaseA_heed_statistics. "
            "For the 2+4 mm library the track starts 0.1 mm inside the +v "
            "cathode, so the Heed sample spans 5.9 mm."
        ),
    )
    ap.add_argument(
        "--wire-phase-from-y",
        action="store_true",
        help=(
            "Fold global y into a 4-mm wire cell assuming a global wire at "
            "y=0. Default is safer: random uniform wire phase."
        ),
    )
    ap.add_argument(
        "--max-tracks",
        type=int,
        default=0,
        help="0 = all rows in input sample",
    )
    args = ap.parse_args()

    if (
        args.gap_minus_mm <= 0
        or args.gap_plus_mm <= 0
        or args.wire_pitch_mm <= 0
        or args.lambda_clusters_per_mm <= 0
        or args.mean_electrons_per_cluster < 1
    ):
        raise SystemExit("Invalid geometry/cluster parameters.")

    rng = np.random.default_rng(args.seed)

    lambda_clusters_per_mm = args.lambda_clusters_per_mm
    if args.heed_tracks_csv:
        with args.heed_tracks_csv.open(newline="", encoding="utf-8") as f:
            heed_tracks = list(csv.DictReader(f))
        if not heed_tracks or "clusters" not in heed_tracks[0]:
            raise SystemExit("Invalid Heed tracks CSV: expected 'clusters' column.")
        if args.heed_reference_length_mm <= 0:
            raise SystemExit("--heed-reference-length-mm must be positive.")
        mean_heed_clusters = float(
            np.mean([float(r["clusters"]) for r in heed_tracks])
        )
        lambda_clusters_per_mm = (
            mean_heed_clusters / args.heed_reference_length_mm
        )

    empirical_cluster_sizes = None
    if args.cluster_electron_model == "heed-empirical":
        if not args.heed_clusters_csv:
            raise SystemExit(
                "heed-empirical mode requires --heed-clusters-csv."
            )
        with args.heed_clusters_csv.open(newline="", encoding="utf-8") as f:
            heed_clusters = list(csv.DictReader(f))
        if not heed_clusters or "electrons" not in heed_clusters[0]:
            raise SystemExit(
                "Invalid Heed clusters CSV: expected 'electrons' column."
            )
        empirical_cluster_sizes = np.asarray(
            [int(float(r["electrons"])) for r in heed_clusters],
            dtype=int,
        )
        empirical_cluster_sizes = empirical_cluster_sizes[
            empirical_cluster_sizes >= 1
        ]
        if len(empirical_cluster_sizes) == 0:
            raise SystemExit("No positive cluster sizes in Heed library.")

    with args.tracks.open(newline="", encoding="utf-8") as f:
        tracks = list(csv.DictReader(f))
    if args.max_tracks > 0:
        tracks = tracks[: args.max_tracks]
    if not tracks:
        raise SystemExit("No track rows found.")

    args.output.parent.mkdir(parents=True, exist_ok=True)

    fields = [
        "synthetic_track_id",
        "channel", "job", "event", "track", "pdg", "station",
        "p_GeV", "theta_u_deg", "theta_w_deg",
        "u0_cell_mm", "wire_phase_mode",
        "track_length_mm", "mean_clusters", "n_clusters_event",
        "cluster_index", "n_primary_electrons",
        "v_mm", "u_unfolded_mm", "wire_offset", "u_cell_mm",
        "w_relative_mm",
    ]

    n_clusters_total = 0
    n_e_total = 0
    n_zero_cluster_tracks = 0

    with args.output.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()

        for it, row in enumerate(tracks):
            tu = float(row["tan_theta_u"])
            tw = float(row["tan_theta_w"])
            norm = math.sqrt(1.0 + tu * tu + tw * tw)
            gap = args.gap_minus_mm + args.gap_plus_mm
            length_mm = gap * norm
            mean_clusters = lambda_clusters_per_mm * length_mm
            ncl = int(rng.poisson(mean_clusters))
            if ncl == 0:
                n_zero_cluster_tracks += 1
                continue

            if args.wire_phase_from_y:
                _, u0_cell = fold_cell(
                    10.0 * float(row["y_cm"]),
                    args.wire_pitch_mm,
                )
                phase_mode = "folded_global_y_assuming_phase0"
            else:
                u0_cell = float(
                    rng.uniform(
                        -0.5 * args.wire_pitch_mm,
                        +0.5 * args.wire_pitch_mm,
                    )
                )
                phase_mode = "uniform_random"

            # Uniform in v is equivalent to uniform in path length for a
            # straight track because ds/dv is constant.
            v_values = rng.uniform(
                -args.gap_minus_mm,
                +args.gap_plus_mm,
                size=ncl,
            )

            for ic, v_mm in enumerate(v_values):
                if args.cluster_electron_model == "one":
                    ne = 1
                elif args.cluster_electron_model == "shifted-poisson":
                    # At least one primary electron belongs to each cluster.
                    extra_mean = args.mean_electrons_per_cluster - 1.0
                    ne = 1 + int(rng.poisson(extra_mean))
                else:
                    ne = int(rng.choice(empirical_cluster_sizes))

                u_unfolded = u0_cell + v_mm * tu
                wire_offset, u_cell = fold_cell(
                    u_unfolded, args.wire_pitch_mm
                )
                w_rel = v_mm * tw

                w.writerow({
                    "synthetic_track_id": it,
                    "channel": row["channel"],
                    "job": row["job"],
                    "event": row["event"],
                    "track": row["track"],
                    "pdg": row["pdg"],
                    "station": row["station"],
                    "p_GeV": row["p_GeV"],
                    "theta_u_deg": row["theta_u_deg"],
                    "theta_w_deg": row["theta_w_deg"],
                    "u0_cell_mm": u0_cell,
                    "wire_phase_mode": phase_mode,
                    "track_length_mm": length_mm,
                    "mean_clusters": mean_clusters,
                    "n_clusters_event": ncl,
                    "cluster_index": ic,
                    "n_primary_electrons": ne,
                    "v_mm": v_mm,
                    "u_unfolded_mm": u_unfolded,
                    "wire_offset": wire_offset,
                    "u_cell_mm": u_cell,
                    "w_relative_mm": w_rel,
                })
                n_clusters_total += 1
                n_e_total += ne

    print("PI Poisson cluster generator")
    print(f"  tracks read                 : {len(tracks)}")
    print(
        f"  cluster density [/mm]       : "
        f"{lambda_clusters_per_mm:.6f}"
    )
    print(
        f"  electron multiplicity model : "
        f"{args.cluster_electron_model}"
    )
    print(f"  tracks with zero clusters   : {n_zero_cluster_tracks}")
    print(f"  clusters written            : {n_clusters_total}")
    print(f"  primary electrons (model)   : {n_e_total}")
    print(
        f"  mean clusters / input track : "
        f"{n_clusters_total / len(tracks):.4f}"
    )
    print(
        f"  mean electrons / cluster    : "
        f"{n_e_total / n_clusters_total:.4f}"
        if n_clusters_total
        else "  mean electrons / cluster    : nan"
    )
    print(f"  wrote                       : {args.output}")
    print()
    if not args.wire_phase_from_y:
        print(
            "Wire phase was sampled uniformly in one pitch. "
            "This is intentional until the old bore scan is tied to an exact "
            "per-chamber wire-origin phase."
        )


if __name__ == "__main__":
    main()
