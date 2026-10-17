#!/usr/bin/env python3
"""Resumable per-CHAMBER MWPC wire collection study for PI.

Input is the *historical GEANT state -> checkerboard gas-box* projection,
NOT actual Geant4 gas hits in the new 298-module geometry.

For EACH projected gas crossing (including double-module overlap tracks):
  Poisson(lambda_Heed * actual gas-path length) primary clusters
  -> empirical Heed primary electrons per cluster
  -> uniformly distributed cluster positions on that gas segment
  -> one sampled microscopic Garfield collection/wire/time result per electron.

Results are keyed by (source, station, det_id), not averaged across a station.
"Active wire" here means >=1 COLLECTED PRIMARY electron; no avalanche gain,
electronics shaping, noise, or hardware threshold is simulated.

Zero-cluster/zero-collected crossings remain in the denominator.
Wire-grid origin is unmeasured, so the position within a pitch is RANDOMIZED
per gas crossing, not invented as a known mechanical origin. No edge-wire
suppression is modeled. The Garfield 2+4 mm asymmetry is mapped onto the
6-mm GEANT gas box via a +1 mm normal-coordinate shift, documented below.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import time
from collections import defaultdict
from pathlib import Path

import numpy as np

SOURCES = ("jpsi40-bore", "phi40-bore", "omega40-bore")
PITCH_MM = 4.
FIELDS = (
    "source", "job", "event", "track", "station", "det_id",
    "gas_path_mm", "theta_u_deg", "theta_w_deg",
    "geom_wire_cells", "heed_clusters", "primary_electrons",
    "collected_electrons", "wires_ge1e", "wires_ge2e", "wires_ge5e",
    "wires_ge10e", "wires_ge1e_100ns", "wires_ge1e_200ns",
)
SUMMARY_FIELDS = (
    "source", "station", "det_id", "processed_chamber_crossings",
    "expected_all_crossings", "completed_source",
    "P_nwire_0", "P_nwire_1", "P_nwire_2",
    "P_nwire_3", "P_nwire_4", "P_nwire_5plus",
    "P_nwire_ge2", "P_nwire_ge3",
    "P_ge2w_threshold_2e", "P_ge2w_threshold_5e",
    "P_ge2w_threshold_10e", "P_ge2w_100ns",
    "P_ge2w_200ns", "mean_collected_primary_electrons",
    "mean_geometric_cells", "mean_heed_clusters",
    "mean_path_mm",
)
NEEDED = {
    "source", "job", "event", "track", "station", "det_id",
    "gas_path_mm", "gas_center_z_cm", "theta_u_deg", "theta_w_deg",
    "y_in_cm", "y_out_cm", "z_in_cm", "z_out_cm",
}


def read_heed(tracks_path, clusters_path):
    with tracks_path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        if "clusters" not in (reader.fieldnames or []):
            raise ValueError("Heed track CSV missing 'clusters' column")
        clusters = np.fromiter(
            (float(row["clusters"]) for row in reader), dtype=np.float64
        )
    if clusters.size == 0:
        raise ValueError("Empty Heed track library")
    with clusters_path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        if "electrons" not in (reader.fieldnames or []):
            raise ValueError("Heed clusters CSV missing 'electrons' column")
        marks = np.fromiter(
            (int(float(row["electrons"])) for row in reader), dtype=np.int32
        )
    marks = marks[marks >= 1]
    if marks.size == 0:
        raise ValueError("No empirical positive Heed electron multiplicities")
    return float(np.mean(clusters) / 5.9), marks


def read_garfield(path):
    """Pack variable-number per-box responses into one indexed sample pool."""
    slots = [[] for _ in range(20 * 30)]
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        required = {
            "iu", "iv", "collected", "wire_index", "drift_time_ns",
            "u_center_mm", "v_center_mm",
        }
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Garfield CSV missing {sorted(missing)}")
        for row in reader:
            iu, iv = int(float(row["iu"])), int(float(row["iv"]))
            if not (0 <= iu < 20 and 0 <= iv < 30):
                raise ValueError(f"Wrong Garfield indices: iu={iu}, iv={iv}")
            uc, vc = float(row["u_center_mm"]), float(row["v_center_mm"])
            if abs(uc - (-1.9 + .2 * iu)) > 1.e-4 or abs(
                    vc - (-1.9 + .2 * iv)) > 1.e-4:
                raise ValueError("Map is not the required [-2,2] x [-2,4] mm grid")
            collected = int(float(row["collected"])) == 1
            wire = int(float(row["wire_index"]))
            try:
                t = float(row["drift_time_ns"])
            except (TypeError, ValueError):
                t = math.nan
            slots[iu * 30 + iv].append((collected, wire, t))
    if any(not s for s in slots):
        raise ValueError(
            "Incomplete transport map: " +
            str(sum(not s for s in slots)) + " of 600 boxes empty"
        )
    offsets = np.zeros(600, dtype=np.int64)
    sizes = np.asarray([len(s) for s in slots], dtype=np.int64)
    offsets[1:] = np.cumsum(sizes[:-1])
    responses = np.asarray([t for s in slots for t in s], dtype=np.float64)
    if not np.all(np.isin(responses[:, 0], (0., 1.))):
        raise ValueError("Nonbinary collected flag")
    print(
        f"Garfield map: 600 boxes, {len(responses):,} sampled electron drifts "
        f"({sizes.min()}..{sizes.max()} per box)", flush=True
    )
    return offsets, sizes, responses


def simulate(row, rng, he_lambda, heed_marks, lookup):
    offsets, sizes, responses = lookup
    path = float(row["gas_path_mm"])
    if not (0 < path <= 1000. and math.isfinite(path)):
        raise ValueError(f"Unphysical gas path length {path}")
    zin = float(row["z_in_cm"])
    zout = float(row["z_out_cm"])
    yin = float(row["y_in_cm"])
    yout = float(row["y_out_cm"])
    zc = float(row["gas_center_z_cm"])
    if not all(math.isfinite(x) for x in (zin, zout, yin, yout, zc)):
        raise ValueError("Nonfinite gas coordinates")
    if zout < zin - 1.e-8:
        raise ValueError("Gas crossing z_out precedes z_in")

    # Wire phase is unknown for these historical projections.
    # Use one coherent random phase per *crossing*, then the actual y slope.
    phase_at_segment_midpoint = rng.uniform(-PITCH_MM / 2, PITCH_MM / 2)
    umid0 = phase_at_segment_midpoint + 10 * (yin - yout) / 2
    umid1 = phase_at_segment_midpoint + 10 * (yout - yin) / 2
    ngeom = 1 + abs(
        math.floor(umid1 / PITCH_MM + .5) -
        math.floor(umid0 / PITCH_MM + .5)
    )

    ncl = int(rng.poisson(he_lambda * path))
    result = dict(
        geom_wire_cells=ngeom, heed_clusters=ncl, primary_electrons=0,
        collected_electrons=0, wires_ge1e=0, wires_ge2e=0,
        wires_ge5e=0, wires_ge10e=0,
        wires_ge1e_100ns=0, wires_ge1e_200ns=0,
    )
    if ncl == 0:
        return result

    s = rng.random(ncl)
    marks = heed_marks[rng.integers(0, len(heed_marks), size=ncl)]
    nprim = int(marks.sum())
    result["primary_electrons"] = nprim
    frac = np.repeat(s, marks)

    u = phase_at_segment_midpoint + 10 * (
        yin + frac * (yout - yin) - .5 * (yin + yout)
    )
    # Transport map is asymmetric (-2..+4 mm); centered GEANT gas box
    # extends -3..+3 mm. The +1 mm translation aligns the two six-mm gaps.
    v = 10 * (zin + frac * (zout - zin) - zc) + 1.0
    if v.min() < -2.001 or v.max() > 4.001:
        raise ValueError(
            f"Projected ionization outside Garfield map: "
            f"v={v.min():.4f}..{v.max():.4f} mm"
        )

    wire0 = np.floor(u / PITCH_MM + .5).astype(np.int64)
    local_u = u - PITCH_MM * wire0
    iu = np.clip(np.floor((local_u + 2.) / .2).astype(np.int64), 0, 19)
    iv = np.clip(np.floor((v + 2.) / .2).astype(np.int64), 0, 29)
    ibox = 30 * iu + iv
    which = offsets[ibox] + np.floor(
        rng.random(nprim) * sizes[ibox]
    ).astype(np.int64)
    trial = responses[which]
    good = trial[:, 0] > .5
    ngood = int(good.sum())
    result["collected_electrons"] = ngood
    if ngood == 0:
        return result

    target_wires = wire0[good] + trial[good, 1].astype(np.int64)
    wires, counts = np.unique(target_wires, return_counts=True)
    result["wires_ge1e"] = len(wires)
    result["wires_ge2e"] = int(np.count_nonzero(counts >= 2))
    result["wires_ge5e"] = int(np.count_nonzero(counts >= 5))
    result["wires_ge10e"] = int(np.count_nonzero(counts >= 10))
    tgood = trial[good, 2]
    for bound in (100, 200):
        selected = np.isfinite(tgood) & (tgood <= bound)
        result[f"wires_ge1e_{bound}ns"] = (
            int(np.unique(target_wires[selected]).size)
        )
    return result


def checkpoint(out_csv, state_path, processed, rng, key):
    """Durably store CSV position BEFORE the state that points to it."""
    out_csv.flush()
    os.fsync(out_csv.fileno())
    info = {
        "source": key, "processed": processed,
        "csv_size_bytes": out_csv.tell(),
        "rng_state": rng.bit_generator.state,
    }
    temp = state_path.with_suffix(".tmp")
    temp.write_text(json.dumps(info) + "\n", encoding="utf-8")
    os.replace(temp, state_path)


def run_source(source, args, deadline, he_lambda, heed_marks, lookup):
    inp = args.crossings / f"crossings_{source}_ms2_530.csv"
    out_file = args.out / f"wire_crossings_{source}.csv"
    state_file = args.out / f"checkpoint_{source}.json"
    if not inp.is_file():
        raise FileNotFoundError(f"Missing chamber-projection CSV: {inp}")
    rng = np.random.default_rng(
        args.seed + SOURCES.index(source) * 10000019
    )
    processed = 0
    csv_size = 0
    if state_file.is_file():
        info = json.loads(state_file.read_text(encoding="utf-8"))
        if info["source"] != source:
            raise ValueError("Checkpoint source mismatch")
        processed = int(info["processed"])
        csv_size = int(info["csv_size_bytes"])
        rng.bit_generator.state = info["rng_state"]
        if not out_file.is_file() or out_file.stat().st_size < csv_size:
            raise ValueError(
                f"CSV lost bytes referenced by checkpoint: {out_file}"
            )
        print(f"Resuming {source} after {processed:,} crossings", flush=True)
    elif out_file.exists():
        raise ValueError(
            f"Output exists without checkpoint; refusing overwrite: {out_file}"
        )

    mode = "r+" if out_file.exists() else "w+"
    with inp.open(newline="", encoding="utf-8") as stream, out_file.open(
        mode, newline="", encoding="utf-8"
    ) as output:
        reader = csv.DictReader(stream)
        missing = NEEDED - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Crossing CSV missing columns {sorted(missing)}")
        writer = csv.DictWriter(output, fieldnames=FIELDS)
        if processed:
            output.seek(csv_size)
            output.truncate()
        else:
            writer.writeheader()
        for _ in range(processed):
            if next(reader, None) is None:
                raise ValueError("Checkpoint exceeds number of crossings")
        for row in reader:
            if row["source"] != source:
                raise ValueError("Crossing source mismatch")
            result = simulate(row, rng, he_lambda, heed_marks, lookup)
            writer.writerow({
                **{k: row[k] for k in FIELDS if k in row},
                **result,
            })
            processed += 1
            if processed % args.checkpoint_every == 0:
                checkpoint(output, state_file, processed, rng, source)
                print(
                    f"  {source}: {processed:,} chamber crossings, "
                    f"elapsed {((time.monotonic() - args.started)/60):.1f} min",
                    flush=True,
                )
                if args.max_per_source and processed >= args.max_per_source:
                    return False
                if time.monotonic() >= deadline:
                    return False
        checkpoint(output, state_file, processed, rng, source)

    # End-of-file validation: each crossed detector is counted separately.
    expected = {}
    with (args.crossings / "modules.csv").open(
        newline="", encoding="utf-8"
    ) as inp_modules:
        for row in csv.DictReader(inp_modules):
            if row["source"] == source and float(
                row["scenario_ms2_z_cm"]
            ) == 530.:
                expected[(int(row["station"]), int(row["det_id"]))] = int(
                    row["muon_crossings"]
                )
    if len(expected) != 298:
        raise ValueError(f"Expected 298 module reference counts for {source}")
    actual = defaultdict(int)
    with out_file.open(newline="", encoding="utf-8") as s:
        for row in csv.DictReader(s):
            actual[(int(row["station"]), int(row["det_id"]))] += 1
    bad = [(k, v, actual[k]) for k, v in expected.items() if actual[k] != v]
    if bad:
        raise ValueError(
            f"Per-chamber crossing count mismatch: {bad[:5]}"
        )
    done = args.out / f"complete_{source}.ok"
    done.write_text(
        f"{processed} crossings checked across 298 chambers\n",
        encoding="utf-8",
    )
    print(
        f"COMPLETE {source}: {processed:,} crossings, all 298 module counts agree",
        flush=True,
    )
    return True


def summaries(args):
    expected = {}
    with (args.crossings / "modules.csv").open(
        newline="", encoding="utf-8"
    ) as f:
        for row in csv.DictReader(f):
            if float(row["scenario_ms2_z_cm"]) == 530.:
                expected[(
                    row["source"], int(row["station"]), int(row["det_id"])
                )] = int(row["muon_crossings"])

    state = defaultdict(lambda: {
        "n": 0, "bins": np.zeros(6, dtype=np.int64),
        "ge2": 0, "ge3": 0, "t2": 0, "t5": 0, "t10": 0,
        "win100": 0, "win200": 0,
        "collected": 0., "geom": 0., "clusters": 0., "path": 0.,
    })
    for source in SOURCES:
        file = args.out / f"wire_crossings_{source}.csv"
        if not file.exists():
            continue
        with file.open(newline="", encoding="utf-8") as f:
            for row in csv.DictReader(f):
                k = (source, int(row["station"]), int(row["det_id"]))
                st = state[k]
                nw = int(row["wires_ge1e"])
                st["n"] += 1
                st["bins"][min(nw, 5)] += 1
                st["ge2"] += (nw >= 2)
                st["ge3"] += (nw >= 3)
                st["t2"] += (int(row["wires_ge2e"]) >= 2)
                st["t5"] += (int(row["wires_ge5e"]) >= 2)
                st["t10"] += (int(row["wires_ge10e"]) >= 2)
                st["win100"] += (int(row["wires_ge1e_100ns"]) >= 2)
                st["win200"] += (int(row["wires_ge1e_200ns"]) >= 2)
                st["collected"] += int(row["collected_electrons"])
                st["geom"] += int(row["geom_wire_cells"])
                st["clusters"] += int(row["heed_clusters"])
                st["path"] += float(row["gas_path_mm"])

    output = args.out / "pi_wire_multiplicity_per_chamber.csv"
    with output.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=SUMMARY_FIELDS)
        writer.writeheader()
        for key in sorted(expected):
            source, station, det_id = key
            st = state[key]
            n = st["n"]
            prob = lambda x: x / n if n else math.nan
            out = {
                "source": source, "station": station, "det_id": det_id,
                "processed_chamber_crossings": n,
                "expected_all_crossings": expected[key],
                "completed_source": int(
                    (args.out / f"complete_{source}.ok").is_file()
                ),
            }
            out.update({
                f"P_nwire_{tag}": prob(st["bins"][i])
                for i, tag in enumerate(("0", "1", "2", "3", "4", "5plus"))
            })
            out.update({
                "P_nwire_ge2": prob(st["ge2"]),
                "P_nwire_ge3": prob(st["ge3"]),
                "P_ge2w_threshold_2e": prob(st["t2"]),
                "P_ge2w_threshold_5e": prob(st["t5"]),
                "P_ge2w_threshold_10e": prob(st["t10"]),
                "P_ge2w_100ns": prob(st["win100"]),
                "P_ge2w_200ns": prob(st["win200"]),
                "mean_collected_primary_electrons": prob(st["collected"]),
                "mean_geometric_cells": prob(st["geom"]),
                "mean_heed_clusters": prob(st["clusters"]),
                "mean_path_mm": prob(st["path"]),
            })
            writer.writerow(out)
    print(f"Per-chamber summary saved: {output}", flush=True)
    return output


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--crossings", type=Path,
        default=Path("test_runs/mwpc_stage_c/checkerboard_full"),
    )
    ap.add_argument(
        "--heed-tracks", type=Path,
        default=Path("test_runs/mwpc_stage_c/c_heed_library_2p4/heed_2p4_tracks.csv"),
    )
    ap.add_argument(
        "--heed-clusters", type=Path,
        default=Path("test_runs/mwpc_stage_c/c_heed_library_2p4/heed_2p4_clusters.csv"),
    )
    ap.add_argument(
        "--garfield", type=Path, default=None,
        help="Override Garfield transport-box file. By default use the "
             "hybrid high-stat boundary map if it exists, otherwise the "
             "verified full 50-trials-per-box baseline.",
    )
    ap.add_argument(
        "--out", type=Path,
        default=Path("test_runs/mwpc_stage_c/pi_per_chamber_response"),
    )
    ap.add_argument("--hours", type=float, default=2.5)
    ap.add_argument("--checkpoint-every", type=int, default=10000)
    ap.add_argument("--seed", type=int, default=271828)
    ap.add_argument("--max-per-source", type=int, default=0,
                    help="Smoke/limited processing only; 0 = no limit")
    args = ap.parse_args()
    if args.hours <= 0 or args.checkpoint_every <= 0:
        ap.error("hours and checkpoint-every must be positive")
    if args.max_per_source < 0:
        ap.error("--max-per-source must be nonnegative")
    if args.garfield is None:
        hybrid = Path(
            "test_runs/mwpc_stage_c/c_transport_box_map_hybrid_boundary_2p4/"
            "c_transport_box_electrons.csv"
        )
        baseline = Path(
            "test_runs/mwpc_stage_c/c_transport_box_map_gm2p0_gp4p0/"
            "c_transport_box_electrons.csv"
        )
        args.garfield = hybrid if hybrid.is_file() else baseline
        print(f"Automatic Garfield map selection: {args.garfield}", flush=True)
    inputs = [args.heed_tracks, args.heed_clusters, args.garfield,
              args.crossings / "modules.csv"]
    missing = [str(p) for p in inputs if not p.is_file()]
    if missing:
        ap.error("Required existing files missing:\n" + "\n".join(missing))
    for source in SOURCES:
        expected_crossings = args.crossings / f"crossings_{source}_ms2_530.csv"
        if not expected_crossings.is_file():
            ap.error(f"Missing checkerboard crossing input: {expected_crossings}")
    args.out.mkdir(parents=True, exist_ok=True)
    descriptor = args.out / "settings.json"
    values = {
        "crossings": str(args.crossings.resolve()),
        "heed_tracks": str(args.heed_tracks.resolve()),
        "heed_clusters": str(args.heed_clusters.resolve()),
        "garfield": str(args.garfield.resolve()),
        "seed": args.seed,
        "fixed_ms2_cm": 530,
        "chamber_count": 298,
        "wire_pitch_mm": PITCH_MM,
        "wire_phase": "uniform random per gas crossing",
        "model": "Poisson(Heed density*actual path) + empirical Heed cluster marks"
                 " + sampled single-electron Garfield destination/time",
        "normal_coordinate_alignment": "v_map_mm=10*(z-z_gas_center)+1,"
            " because Garfield gap=[-2,+4]mm but projected gas=[-3,+3]mm",
        "limitations": [
            "NOT actual Geant4 hits/energy loss in modern geometry",
            "no avalanche gain, electronics threshold, ion signal or strip signal",
            "wires_ge5e uses primary collected electrons, not charge",
            "no cathode strips or wire-edge losses",
            "uniform randomized wire phase is not known actual wire origin",
            "Garfield B=0 Tesla reference map and 2+4 mm asymmetric gap",
            "different overlapping modules get separate simulated ionization",
            "timing refers to electron drift, not final front-end firing time",
        ],
    }
    if descriptor.exists():
        saved = json.loads(descriptor.read_text(encoding="utf-8"))
        if saved != values:
            ap.error("Input/settings changed since checkpoint. "
                     "Choose a fresh --out directory.")
    else:
        descriptor.write_text(
            json.dumps(values, indent=2) + "\n", encoding="utf-8"
        )
    lam, marks = read_heed(args.heed_tracks, args.heed_clusters)
    response = read_garfield(args.garfield)
    print(
        f"Heed: {lam:.6f} clusters/mm, {len(marks):,} observed cluster marks",
        flush=True,
    )
    args.started = time.monotonic()
    deadline = args.started + args.hours * 3600
    stopped = False
    try:
        for source in SOURCES:
            if (args.out / f"complete_{source}.ok").exists():
                print(f"SKIP completed source: {source}", flush=True)
                continue
            if time.monotonic() >= deadline:
                stopped = True
                break
            if not run_source(source, args, deadline, lam, marks, response):
                stopped = True
                break
    finally:
        summaries(args)
    if stopped:
        print(
            "PAUSED after checkpoint: wall-clock budget or per-source limit."
            " Re-run same command to continue without losing completed work.",
            flush=True,
        )
    else:
        print("ALL THREE CHANNELS COMPLETED and per-module counts validated.",
              flush=True)


if __name__ == "__main__":
    main()
