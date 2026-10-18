#!/usr/bin/env python3
"""Build one offline interactive NA60+ MWPC station dashboard.

No server, CDN, ROOT or Garfield required. Input is the ALREADY completed:
 - checkerboard_full/{modules.csv,stations.csv}
 - pi_fixed_wire_maps/pi_fixed_wire_counts.csv
 - pi_fixed_wire_maps/physical_validation/pi_geom_wire_counts.csv
 - pi_fixed_wire_maps/physical_validation/pi_crossing_density_xy.npz

Output: a single self-contained HTML; opens with file:// in Chrome/Firefox.
Numbers in the browser are real existing simulation outputs, not synthetic
mock data. Input consistency is checked for every source and every anode.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np

from mwpc_checkerboard_projection import SOURCES, build_chambers
from pi_fixed_wire_occupancy import n_wires

GEOMETRY = build_chambers()
CHAMBERS = {ch.det_id: ch for station in GEOMETRY for ch in station}
assert len(CHAMBERS) == 298
EXPECTED_WIRES = sum(n_wires(c.station) for c in CHAMBERS.values())
assert EXPECTED_WIRES == 36512


def geometry_payload():
    return [
        [{
            "id": c.det_id, "x": round(c.x, 4), "y": round(c.y, 4),
            "layer": c.layer, "row": c.row, "col": c.col,
        } for c in station]
        for station in GEOMETRY
    ]


def load_module_counts(path):
    counts = {(source, ch.det_id): None
              for source in SOURCES for ch in CHAMBERS.values()}
    with path.open(newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            if abs(float(r["scenario_ms2_z_cm"]) - 530.) > 1e-6:
                continue
            source, det = r["source"], int(r["det_id"])
            if (source, det) not in counts:
                raise ValueError(f"Unexpected module: {source}, {det}")
            if counts[source, det] is not None:
                raise ValueError("Duplicate reference module row")
            if int(r["station"]) != CHAMBERS[det].station:
                raise ValueError("Module/station mismatch")
            counts[source, det] = int(r["muon_crossings"])
    if any(n is None for n in counts.values()):
        raise ValueError("Reference module table missing chambers")
    return counts


def load_wire_tables(wire_file, geo_file, module_counts):
    payload = {source: {
        str(det): {"n": int(module_counts[source, det]),
                   "geom": [0]*n_wires(ch.station),
                   "e1": [0]*n_wires(ch.station),
                   "e5": [0]*n_wires(ch.station),
                   "t100": [0]*n_wires(ch.station)}
        for det, ch in CHAMBERS.items()
    } for source in SOURCES}
    seen = set()
    with wire_file.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        needed = {
            "source", "station", "det_id", "wire_index", "chamber_crossings",
            "n_ge1e", "n_ge5e", "n_le100ns",
        }
        if not needed.issubset(reader.fieldnames or []):
            raise ValueError(f"Wire CSV missing fields {needed-set(reader.fieldnames or [])}")
        for r in reader:
            source, det, j = r["source"], int(r["det_id"]), int(r["wire_index"])
            if source not in payload or det not in CHAMBERS:
                raise ValueError("Unrecognized wire source or module")
            c = CHAMBERS[det]
            if int(r["station"]) != c.station:
                raise ValueError("Wire station/det mismatch")
            if not 0 <= j < n_wires(c.station):
                raise ValueError("Out-of-range physical wire")
            key = source, det, j
            if key in seen:
                raise ValueError(f"Duplicate physical wire {key}")
            seen.add(key)
            chamber = payload[source][str(det)]
            n = chamber["n"]
            if int(r["chamber_crossings"]) != n:
                raise ValueError(f"Wire denominator mismatch for {key}")
            e1, e5, t100 = (int(r[k]) for k in (
                "n_ge1e", "n_ge5e", "n_le100ns"
            ))
            if not (0 <= e5 <= e1 <= n and 0 <= t100 <= e1):
                raise ValueError(f"Physically inconsistent wire counts {key}")
            chamber["e1"][j] = e1
            chamber["e5"][j] = e5
            chamber["t100"][j] = t100

    expected = len(SOURCES) * EXPECTED_WIRES
    if len(seen) != expected:
        raise ValueError(
            f"Expected {expected} fixed physical anodes, got {len(seen)}"
        )

    seen.clear()
    with geo_file.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        needed = {
            "source", "station", "det_id", "wire_index", "N_crossings",
            "n_geom",
        }
        if not needed.issubset(reader.fieldnames or []):
            raise ValueError(f"Geometry CSV missing fields {needed-set(reader.fieldnames or [])}")
        for r in reader:
            source, det, j = r["source"], int(r["det_id"]), int(r["wire_index"])
            key = source, det, j
            if key in seen:
                raise ValueError(f"Duplicate geometric wire: {key}")
            seen.add(key)
            if source not in payload or det not in CHAMBERS:
                raise ValueError("Unrecognized geometric wire source or module")
            ch = payload[source][str(det)]
            if not 0 <= j < len(ch["geom"]):
                raise ValueError("Out-of-range geometric wire")
            n = ch["n"]
            count = int(r["n_geom"])
            if int(r["N_crossings"]) != n or not 0 <= count <= n:
                raise ValueError(f"Geometric normalization mismatch: {key}")
            if int(r["station"]) != CHAMBERS[det].station:
                raise ValueError("Geometric station/det mismatch")
            ch["geom"][j] = count
    if len(seen) != expected:
        raise ValueError(
            f"Expected {expected} geometric wire identifiers, got {len(seen)}"
        )
    return payload


def stations_payload(path, xy_file, module_counts):
    payload = {source: [None]*6 for source in SOURCES}
    with path.open(newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            if abs(float(row["scenario_ms2_z_cm"]) - 530.) > 1e-6:
                continue
            src, st = row["source"], int(row["station"])
            if src not in payload or not 0 <= st < 6:
                raise ValueError("Unknown source or station")
            if payload[src][st] is not None:
                raise ValueError("Duplicate station summary")
            refs = int(row["reference_states"])
            covered = int(row["projected_muons_with_gas"])
            crossings = int(row["chamber_crossings_total"])
            expected = sum(module_counts[src, ch.det_id] for ch in GEOMETRY[st])
            if not 0 <= covered <= refs or crossings < covered or crossings != expected:
                raise ValueError(f"Station count conservation failed: {src} MS{st}")
            payload[src][st] = {
                "incident": refs, "withGas": covered, "crossings": crossings,
            }

    with np.load(xy_file, allow_pickle=False) as arrays:
        for source in SOURCES:
            for station in range(6):
                entry = payload[source][station]
                if entry is None:
                    raise ValueError(f"Missing station summary {source} MS{station}")
                key = f"{source}_MS{station}"
                for suffix in ("", "_xedges", "_yedges"):
                    if key + suffix not in arrays:
                        raise ValueError(f"Missing XY histogram {key+suffix}")
                grid = arrays[key]
                xedges, yedges = arrays[key+"_xedges"], arrays[key+"_yedges"]
                if grid.ndim != 2 or grid.shape != (
                    len(yedges)-1, len(xedges)-1
                ):
                    raise ValueError(f"Bad XY histogram shape for {key}")
                if int(grid.sum()) != entry["crossings"]:
                    raise ValueError(
                        f"XY histogram sum does not match module counts: {key}"
                    )
                entry["xy"] = {
                    "counts": grid.astype(int).tolist(),
                    "xedges": xedges.astype(float).round(5).tolist(),
                    "yedges": yedges.astype(float).round(5).tolist(),
                }
    return payload


def generate(args):
    files = [
        args.wire_counts, args.geom_counts, args.crossings/"modules.csv",
        args.crossings/"stations.csv", args.xy, args.template,
    ]
    missing = [str(p) for p in files if not p.is_file()]
    if missing:
        raise FileNotFoundError(
            "Dashboard needs these existing Stage C files:\n" +
            "\n".join(missing)
        )
    counts = load_module_counts(args.crossings/"modules.csv")
    chambers = load_wire_tables(args.wire_counts, args.geom_counts, counts)
    stations = stations_payload(args.crossings/"stations.csv", args.xy, counts)
    data = {
        "description": "200 000 событий/канал в историческом наборе; MS2=530 см",
        "sources": list(SOURCES),
        "geometry": geometry_payload(),
        "stations": stations,
        "chambers": chambers,
        "limitations": [
            "Historical Geant4 station reference states, not native gas hits",
            "Station and module gas crossings can double-count tracks at overlaps",
            "77/128 nominal physical wire positions with unverified wire origin",
            "Electron collection, not electronics discriminator firing",
            "Garfield zero-B reference map, asymmetric gap approximates 6-mm gas",
        ],
    }
    embedded = json.dumps(data, ensure_ascii=False, separators=(",", ":"))
    embedded = embedded.replace("</", "<\\/")
    template = args.template.read_text(encoding="utf-8")
    marker = "__PI_STATION_DATA__"
    if template.count(marker) != 1:
        raise ValueError("Template must contain exactly one JSON marker")
    html = template.replace(marker, embedded)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(html, encoding="utf-8")
    n_rows = 3*EXPECTED_WIRES
    print(f"Validated: {n_rows:,} anode rows, 298 modules, 18 station XY histograms")
    for source in SOURCES:
        print(
            f"  {source}: {sum(stations[source][s]['crossings'] for s in range(6)):,}"
            " gas-module crossings"
        )
    print(
        f"Offline interactive dashboard: {args.out.resolve()} "
        f"({args.out.stat().st_size/1048576:.1f} MiB)"
    )
    print("Open it as a normal local .html file in a modern web browser.")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    folder = Path("test_runs/mwpc_stage_c/pi_fixed_wire_maps")
    stage_c = Path("test/mwpc-stage-c")
    ap.add_argument("--crossings", type=Path,
                    default=Path("test_runs/mwpc_stage_c/checkerboard_full"))
    ap.add_argument("--wire-counts", type=Path,
                    default=folder/"pi_fixed_wire_counts.csv")
    ap.add_argument("--geom-counts", type=Path,
                    default=folder/"physical_validation/pi_geom_wire_counts.csv")
    ap.add_argument("--xy", type=Path,
                    default=folder/"physical_validation/pi_crossing_density_xy.npz")
    ap.add_argument("--template", type=Path,
                    default=stage_c/"pi_station_dashboard_template.html")
    ap.add_argument("--out", type=Path,
                    default=folder/"pi_station_interactive.html")
    args = ap.parse_args()
    generate(args)


if __name__ == "__main__":
    main()
