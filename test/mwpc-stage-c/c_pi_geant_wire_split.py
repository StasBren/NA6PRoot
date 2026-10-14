#!/usr/bin/env python3
"""
PI-style geometric wire sharing from real GEANT MWPC gas crossings.

For each straight GEANT gas crossing:
  1) transform/use the saved local entry and exit positions (u,v,w);
  2) split the path exactly at mid-planes between 4-mm anode wires;
  3) generate a Poisson number of ionisation clusters along the full path;
  4) conditionally distribute those clusters between crossed wire cells with
     a binomial/multinomial law proportional to path length;
  5) assign each cluster a real empirical Heed mark:
       - number of primary electrons;
       - optional Heed energy transfer;
  6) sum cluster/electron/energy content per wire.

This is intentionally the simple PI model.  It does NOT drift electrons in
Garfield and does NOT use the Stage-C box transport map.

Mathematical equivalence:
  N_cl ~ Poisson(lambda L)
  (N_1,...,N_k) | N_cl ~ Multinomial(N_cl; L_1/L,...,L_k/L)

is equivalent to independent
  N_i ~ Poisson(lambda L_i).

GEANT Edep is kept as an independent diagnostic.  By default it is NOT used to
renormalise the Heed cluster population, avoiding double counting of energy-loss
fluctuations.
"""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path

import numpy as np


def nearest_wire(u_mm: float, pitch_mm: float) -> int:
    return math.floor(u_mm / pitch_mm + 0.5)


def split_segment_by_wire(
    u0: float,
    u1: float,
    total_length_mm: float,
    pitch_mm: float,
):
    """Return ordered (wire_index, s0, s1, length_mm) pieces."""
    if total_length_mm <= 0:
        return []

    cuts = [0.0, 1.0]
    du = u1 - u0
    if abs(du) > 1e-15:
        lo = min(u0, u1)
        hi = max(u0, u1)

        # Boundaries lie at (k + 1/2) * pitch.
        kmin = math.floor(lo / pitch_mm - 0.5) - 1
        kmax = math.ceil(hi / pitch_mm - 0.5) + 1
        for k in range(kmin, kmax + 1):
            boundary = (k + 0.5) * pitch_mm
            if lo < boundary < hi:
                s = (boundary - u0) / du
                if 0.0 < s < 1.0:
                    cuts.append(float(s))

    cuts = sorted(set(round(x, 15) for x in cuts))
    pieces = []
    for s0, s1 in zip(cuts, cuts[1:]):
        if s1 <= s0:
            continue
        smid = 0.5 * (s0 + s1)
        umid = u0 + smid * (u1 - u0)
        wire = nearest_wire(umid, pitch_mm)
        pieces.append(
            {
                "wire_index": wire,
                "s0": s0,
                "s1": s1,
                "length_mm": total_length_mm * (s1 - s0),
            }
        )

    # Merge adjacent pieces with the same wire, if floating-point boundaries
    # ever generated a redundant cut.
    merged = []
    for p in pieces:
        if merged and merged[-1]["wire_index"] == p["wire_index"]:
            merged[-1]["s1"] = p["s1"]
            merged[-1]["length_mm"] += p["length_mm"]
        else:
            merged.append(dict(p))
    return merged


def read_heed_lambda(path: Path, reference_length_mm: float) -> float:
    with path.open(newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    if not rows or "clusters" not in rows[0]:
        raise SystemExit("Invalid Heed tracks CSV: missing 'clusters'.")
    mean_n = float(np.mean([float(r["clusters"]) for r in rows]))
    return mean_n / reference_length_mm


def read_heed_marks(path: Path):
    with path.open(newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    if not rows or "electrons" not in rows[0]:
        raise SystemExit("Invalid Heed clusters CSV: missing 'electrons'.")

    electrons = np.asarray(
        [int(float(r["electrons"])) for r in rows],
        dtype=np.int64,
    )
    keep = electrons >= 1
    electrons = electrons[keep]

    energy = None
    if "energy_transfer_eV" in rows[0]:
        energy = np.asarray(
            [float(r["energy_transfer_eV"]) for r in rows],
            dtype=float,
        )[keep]
    return electrons, energy


def load_selection(path: Path | None, require_abs_pdg: int):
    if path is None:
        return None
    with path.open(newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    if not rows:
        raise SystemExit("Selection CSV is empty.")

    required = {"event", "track", "station"}
    if not required.issubset(rows[0]):
        raise SystemExit(
            "Selection CSV must contain event, track, station columns."
        )

    selected = set()
    for r in rows:
        if "pdg" in r and require_abs_pdg > 0:
            if abs(int(float(r["pdg"]))) != require_abs_pdg:
                continue
        selected.add(
            (
                int(float(r["event"])),
                int(float(r["track"])),
                int(float(r["station"])),
            )
        )
    return selected


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--crossings", required=True, type=Path)
    ap.add_argument("--heed-tracks-csv", required=True, type=Path)
    ap.add_argument("--heed-clusters-csv", required=True, type=Path)
    ap.add_argument("--output-prefix", required=True, type=Path)
    ap.add_argument("--selection-csv", type=Path)
    ap.add_argument(
        "--require-abs-pdg",
        type=int,
        default=13,
        help="Applied only when selection CSV has a pdg column.",
    )
    ap.add_argument("--wire-pitch-mm", type=float, default=4.0)
    ap.add_argument("--heed-reference-length-mm", type=float, default=5.9)
    ap.add_argument("--seed", type=int, default=200101)
    ap.add_argument("--max-crossings", type=int, default=0)
    args = ap.parse_args()

    if args.wire_pitch_mm <= 0 or args.heed_reference_length_mm <= 0:
        raise SystemExit("Pitch and Heed reference length must be positive.")

    rng = np.random.default_rng(args.seed)
    lam = read_heed_lambda(
        args.heed_tracks_csv,
        args.heed_reference_length_mm,
    )
    heed_ne, heed_de = read_heed_marks(args.heed_clusters_csv)
    selection = load_selection(args.selection_csv, args.require_abs_pdg)

    with args.crossings.open(newline="", encoding="utf-8") as f:
        crossings = list(csv.DictReader(f))
    if not crossings:
        raise SystemExit("Crossing CSV is empty.")

    needed = {
        "event", "track_id", "det_id", "station",
        "u_in_mm", "v_in_mm", "w_in_mm",
        "u_out_mm", "v_out_mm", "w_out_mm",
        "path_length_mm", "edep_eV",
    }
    missing = needed - set(crossings[0])
    if missing:
        raise SystemExit(
            "Crossing CSV missing columns: " + ", ".join(sorted(missing))
        )

    if selection is not None:
        crossings = [
            r for r in crossings
            if (
                int(float(r["event"])),
                int(float(r["track_id"])),
                int(float(r["station"])),
            ) in selection
        ]

    if args.max_crossings > 0:
        crossings = crossings[: args.max_crossings]
    if not crossings:
        raise SystemExit("No crossings remain after selection.")

    prefix = args.output_prefix
    prefix.parent.mkdir(parents=True, exist_ok=True)
    event_path = Path(str(prefix) + "_crossings.csv")
    wire_path = Path(str(prefix) + "_wires.csv")

    event_fields = [
        "crossing_id", "event", "track_id", "det_id", "station",
        "u_in_mm", "v_in_mm", "w_in_mm",
        "u_out_mm", "v_out_mm", "w_out_mm",
        "path_length_mm", "geant_edep_eV",
        "wire_cells_crossed",
        "first_wire_index", "last_wire_index",
        "mean_clusters_poisson", "n_clusters_total",
        "n_primary_electrons_total",
        "synthetic_heed_transfer_eV",
        "synthetic_over_geant_edep",
        "dominant_wire_index",
        "dominant_wire_electron_fraction",
        "multiwire_cluster_event",
    ]
    wire_fields = [
        "crossing_id", "event", "track_id", "det_id", "station",
        "wire_index", "wire_u_mm",
        "segment_length_mm", "path_fraction",
        "n_clusters", "cluster_fraction",
        "n_primary_electrons", "electron_fraction",
        "synthetic_heed_transfer_eV", "transfer_fraction",
    ]

    total_clusters = 0
    total_electrons = 0
    n_multi_geom = 0
    n_multi_cluster = 0

    with (
        event_path.open("w", newline="", encoding="utf-8") as fe,
        wire_path.open("w", newline="", encoding="utf-8") as fw,
    ):
        ew = csv.DictWriter(fe, fieldnames=event_fields)
        ww = csv.DictWriter(fw, fieldnames=wire_fields)
        ew.writeheader()
        ww.writeheader()

        for cid, r in enumerate(crossings):
            u0 = float(r["u_in_mm"])
            v0 = float(r["v_in_mm"])
            w0 = float(r["w_in_mm"])
            u1 = float(r["u_out_mm"])
            v1 = float(r["v_out_mm"])
            w1 = float(r["w_out_mm"])
            length = float(r["path_length_mm"])
            geant_edep = float(r["edep_eV"])

            if length <= 0:
                continue

            pieces = split_segment_by_wire(
                u0,
                u1,
                length,
                args.wire_pitch_mm,
            )
            if not pieces:
                continue
            if len(pieces) > 1:
                n_multi_geom += 1

            probs = np.asarray(
                [p["length_mm"] / length for p in pieces],
                dtype=float,
            )
            probs /= probs.sum()

            mean_n = lam * length
            ncl_total = int(rng.poisson(mean_n))
            if ncl_total > 0:
                ncl_by_piece = rng.multinomial(ncl_total, probs)
                mark_idx = rng.integers(
                    0,
                    len(heed_ne),
                    size=ncl_total,
                )
                marks_ne = heed_ne[mark_idx]
                marks_de = (
                    heed_de[mark_idx]
                    if heed_de is not None
                    else np.zeros(ncl_total, dtype=float)
                )
            else:
                ncl_by_piece = np.zeros(len(pieces), dtype=int)
                marks_ne = np.empty(0, dtype=int)
                marks_de = np.empty(0, dtype=float)

            # Randomise cluster ordering before partitioning it into multinomial
            # categories.  This makes the empirical Heed marks independent of
            # wire assignment.
            if ncl_total > 1:
                order = rng.permutation(ncl_total)
                marks_ne = marks_ne[order]
                marks_de = marks_de[order]

            offset = 0
            wire_rows = []
            for ip, (p, ncl) in enumerate(zip(pieces, ncl_by_piece)):
                ncl = int(ncl)
                sl = slice(offset, offset + ncl)
                ne = int(np.sum(marks_ne[sl])) if ncl else 0
                de = float(np.sum(marks_de[sl])) if ncl else 0.0
                offset += ncl
                wire_rows.append(
                    {
                        **p,
                        "n_clusters": ncl,
                        "n_primary_electrons": ne,
                        "synthetic_heed_transfer_eV": de,
                    }
                )

            ne_total = sum(x["n_primary_electrons"] for x in wire_rows)
            de_total = sum(x["synthetic_heed_transfer_eV"] for x in wire_rows)

            nonzero_wires = [
                x for x in wire_rows
                if x["n_primary_electrons"] > 0
            ]
            if len(nonzero_wires) > 1:
                n_multi_cluster += 1

            if nonzero_wires:
                dominant = max(
                    nonzero_wires,
                    key=lambda x: x["n_primary_electrons"],
                )
                dominant_wire = dominant["wire_index"]
                dominant_fraction = (
                    dominant["n_primary_electrons"] / ne_total
                    if ne_total > 0
                    else math.nan
                )
            else:
                dominant_wire = 999999
                dominant_fraction = math.nan

            ew.writerow(
                {
                    "crossing_id": cid,
                    "event": r["event"],
                    "track_id": r["track_id"],
                    "det_id": r["det_id"],
                    "station": r["station"],
                    "u_in_mm": u0,
                    "v_in_mm": v0,
                    "w_in_mm": w0,
                    "u_out_mm": u1,
                    "v_out_mm": v1,
                    "w_out_mm": w1,
                    "path_length_mm": length,
                    "geant_edep_eV": geant_edep,
                    "wire_cells_crossed": len(pieces),
                    "first_wire_index": pieces[0]["wire_index"],
                    "last_wire_index": pieces[-1]["wire_index"],
                    "mean_clusters_poisson": mean_n,
                    "n_clusters_total": ncl_total,
                    "n_primary_electrons_total": ne_total,
                    "synthetic_heed_transfer_eV": de_total,
                    "synthetic_over_geant_edep": (
                        de_total / geant_edep
                        if geant_edep > 0
                        else math.nan
                    ),
                    "dominant_wire_index": dominant_wire,
                    "dominant_wire_electron_fraction": dominant_fraction,
                    "multiwire_cluster_event": int(len(nonzero_wires) > 1),
                }
            )

            for x in wire_rows:
                ww.writerow(
                    {
                        "crossing_id": cid,
                        "event": r["event"],
                        "track_id": r["track_id"],
                        "det_id": r["det_id"],
                        "station": r["station"],
                        "wire_index": x["wire_index"],
                        "wire_u_mm": x["wire_index"] * args.wire_pitch_mm,
                        "segment_length_mm": x["length_mm"],
                        "path_fraction": x["length_mm"] / length,
                        "n_clusters": x["n_clusters"],
                        "cluster_fraction": (
                            x["n_clusters"] / ncl_total
                            if ncl_total > 0
                            else math.nan
                        ),
                        "n_primary_electrons": x["n_primary_electrons"],
                        "electron_fraction": (
                            x["n_primary_electrons"] / ne_total
                            if ne_total > 0
                            else math.nan
                        ),
                        "synthetic_heed_transfer_eV": (
                            x["synthetic_heed_transfer_eV"]
                        ),
                        "transfer_fraction": (
                            x["synthetic_heed_transfer_eV"] / de_total
                            if de_total > 0
                            else math.nan
                        ),
                    }
                )

            total_clusters += ncl_total
            total_electrons += ne_total

    print("PI geometric wire-sharing model COMPLETE")
    print(f"  GEANT crossings used          : {len(crossings)}")
    print(f"  Heed cluster density [/mm]    : {lam:.6f}")
    print(f"  geometric multi-cell tracks   : {n_multi_geom}")
    print(f"  stochastic multi-wire events  : {n_multi_cluster}")
    print(f"  clusters sampled              : {total_clusters}")
    print(f"  primary electrons sampled     : {total_electrons}")
    print(f"  crossing summary              : {event_path}")
    print(f"  wire summary                  : {wire_path}")
    print()
    print(
        "GEANT Edep is diagnostic only; the Poisson/Heed population is not "
        "renormalised to it."
    )


if __name__ == "__main__":
    main()
