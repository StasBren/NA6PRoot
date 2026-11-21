#!/usr/bin/env python3
"""Gmsh 4.x geometry-only validity gate for a curved 10-node-tetrahedral mesh.

Checks exact-adaptive Jacobian determinant minDetJac and sampled minSJ.
Does not imply Elmer/Garfield map interpolation, physical mesh convergence,
or validity of the full nine-wire chamber model. Designed for a small trial.
"""
from __future__ import annotations

import argparse
import math
import sys


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("mesh", help="Fresh Gmsh .msh file, not an Elmer mesh")
    parser.add_argument("--chunk", type=int, default=8000)
    args = parser.parse_args()
    if args.chunk <= 0:
        parser.error("--chunk must be positive")
    try:
        import gmsh
    except ImportError as err:
        print("FEM MESH NOT VALIDATED: install/activate matching Gmsh Python API:", err)
        return 4
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.open(args.mesh)
        element_types, groups, _ = gmsh.model.mesh.getElements(3)
        if not element_types:
            raise RuntimeError("No 3D elements")
        if any(int(t) != 11 for t in element_types):
            raise RuntimeError(
                "Non-Tet10 volume element type(s): "
                + ",".join(map(str, map(int, element_types)))
            )
        n_total = n_bad_j = n_bad_sj = 0
        min_det = min_sj = math.inf
        first_bad = []
        for tags in groups:
            for k in range(0, len(tags), args.chunk):
                batch = tags[k : k + args.chunk]
                q_det = gmsh.model.mesh.getElementQualities(
                    batch, qualityName="minDetJac"
                )
                q_sj = gmsh.model.mesh.getElementQualities(
                    batch, qualityName="minSJ"
                )
                if len(q_det) != len(batch) or len(q_sj) != len(batch):
                    raise RuntimeError("Incomplete element-quality output")
                for tag, d, s in zip(batch, q_det, q_sj):
                    d, s = float(d), float(s)
                    if not (math.isfinite(d) and math.isfinite(s)):
                        raise RuntimeError(f"Nonfinite mesh metric: {int(tag)}")
                    n_total += 1
                    min_det = min(min_det, d)
                    min_sj = min(min_sj, s)
                    if d <= 0:
                        n_bad_j += 1
                        if len(first_bad) < 12:
                            first_bad.append(int(tag))
                    if s <= 0:
                        n_bad_sj += 1
        print(f"TET10 count={n_total}")
        print(f"minDetJac={min_det:.12g}   nonpositive_count={n_bad_j}")
        print(f"minSJ={min_sj:.12g}   nonpositive_count={n_bad_sj}")
        print(f"first_bad_tags={first_bad}")
        if n_total < 1 or n_bad_j or n_bad_sj:
            print("FAIL: negative/zero P2 tetrahedral Jacobian")
            return 3
        print("PASS geometry-only: all checked curved Tet10 Jacobians positive")
        print("NOT a full-detector FEM or Garfield field-accuracy certificate")
        return 0
    except Exception as err:
        print(f"FEM MESH NOT VALIDATED: {err}", file=sys.stderr)
        return 2
    finally:
        gmsh.finalize()


if __name__ == "__main__":
    raise SystemExit(main())
