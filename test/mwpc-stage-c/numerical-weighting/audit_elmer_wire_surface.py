#!/usr/bin/env python3
"""Read-only geometric audit of actual Elmer second-order wire-surface nodes.

Input is the already-generated mesh.nodes + mesh.boundary; does not modify maps.
Elmer mesh.boundary: id, physical_boundary, parent1, parent2, elem_type, node_ids.
Elmer mesh.nodes: id, tag, u_mm, v_mm, w_mm.
Output distinguishes the 3 surface triangle vertices from 3 mid-edge nodes.
This measures actual nodes, NOT field accuracy or mesh-cell Jacobians.
"""
from __future__ import annotations

import argparse
import math
import statistics
import sys
from pathlib import Path


def percentile(values, q):
    seq = sorted(values)
    z = (len(seq) - 1) * q
    lo = int(z)
    hi = min(lo + 1, len(seq) - 1)
    return seq[lo] + (z - lo) * (seq[hi] - seq[lo])


def summarize(label, deviations_um, n_wires):
    if not deviations_um:
        raise RuntimeError(f"No {label} nodes were found")
    arr = deviations_um
    print(f"{label}: n={len(arr):,}; wire_ids={sorted(n_wires)}")
    print("  r-R [um]: min={:.6f}, p01={:.6f}, p50={:.6f}, "
          "p99={:.6f}, max={:.6f}, RMS={:.6f}".format(
              min(arr), percentile(arr, .01), statistics.median(arr),
              percentile(arr, .99), max(arr),
              math.sqrt(sum(v*v for v in arr) / len(arr))))
    print("  max inward departure [um]={:.6f}".format(-min(arr)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--map-dir", type=Path, required=True,
                        help="Existing FEM run directory containing chamber/mesh.*")
    parser.add_argument("--wire-pitch-mm", type=float, default=3.)
    parser.add_argument("--wire-radius-mm", type=float, default=.015)
    parser.add_argument("--half-wires", type=int, default=4)
    parser.add_argument("--wire-boundary-id", type=int, default=4)
    args = parser.parse_args()
    chamber = args.map_dir / "chamber"
    bound = chamber / "mesh.boundary"
    nodes = chamber / "mesh.nodes"
    if not bound.is_file() or not nodes.is_file():
        parser.error(f"Missing {bound} or {nodes}")
    corners, mids = set(), set()
    faces = 0
    element_types = set()
    with bound.open(encoding="utf-8") as stream:
        for line in stream:
            p = line.split()
            if len(p) < 5 or int(p[1]) != args.wire_boundary_id:
                continue
            etype = int(p[4])
            element_types.add(etype)
            if etype != 306 or len(p) != 11:
                raise RuntimeError(
                    f"Expected Elmer type-306 (6-node triangle) on wire surface, got {etype}: {line[:150]}")
            corners.update(map(int, p[5:8]))
            mids.update(map(int, p[8:11]))
            faces += 1
    if not faces:
        raise RuntimeError(
            f"No wire faces with boundary id {args.wire_boundary_id}; "
            "check map/manifest.json boundary_ids_after_elmergrid_autoclean")
    overlap = corners & mids
    if overlap:
        raise RuntimeError(f"Unexpected corner/mid-node ID overlap: {len(overlap)}")
    print(f"Wire surface: boundary_id={args.wire_boundary_id}, "
          f"6-node triangles={faces:,}, unique corner IDs={len(corners):,}, "
          f"mid-edge IDs={len(mids):,}")
    deviations = {"corners": [], "mid_edges": []}
    seen = {"corners": set(), "mid_edges": set()}
    wire_seen = {"corners": set(), "mid_edges": set()}
    with nodes.open(encoding="utf-8") as stream:
        for line in stream:
            p = line.split()
            if len(p) < 5:
                continue
            node_id = int(p[0])
            if node_id not in corners and node_id not in mids:
                continue
            u, v = float(p[2]), float(p[3])
            iwire = round(u / args.wire_pitch_mm)
            if abs(iwire) > args.half_wires:
                raise RuntimeError(
                    f"Boundary node {node_id} u={u} mm outside expected wire bank")
            radial = math.hypot(u - iwire * args.wire_pitch_mm, v)
            delta_um = (radial - args.wire_radius_mm) * 1000.
            key = "corners" if node_id in corners else "mid_edges"
            deviations[key].append(delta_um)
            seen[key].add(node_id)
            wire_seen[key].add(iwire)
    for key, wanted in (("corners", corners), ("mid_edges", mids)):
        if seen[key] != wanted:
            raise RuntimeError(f"Unresolved {key} node IDs: {len(wanted - seen[key])}")
        summarize(key, deviations[key], wire_seen[key])
    print("Physical wire radius [um] =", args.wire_radius_mm * 1000.)
    print("First near-wire test shell is r=1.1R; clearance [um] =",
          args.wire_radius_mm * 100.)
    print("NOTE: radial deficits quantify polygonization of the mesh nodes; "
          "this alone does not validate gradients, strip BCs or currents.")


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, RuntimeError) as exc:
        print(f"WIRE_SURFACE_AUDIT_ERROR: {exc}", file=sys.stderr)
        sys.exit(2)
