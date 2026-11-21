#!/usr/bin/env python3
"""Read-only location report for failed curved-P2 Tet10 elements.

Run with the isolated gmsh 4.12.1 Python environment. Geometry in mm.
Uses the audited list of bad element tags; never changes the mesh.
"""
from __future__ import annotations
import argparse
import ast
import json
import math
from pathlib import Path

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mesh", type=Path, help="Existing chamber.msh")
    ap.add_argument("audit", type=Path, help="Prior Jacobian audit log")
    args = ap.parse_args()
    import gmsh
    gfile = args.mesh.resolve()
    cfgfile = gfile.parent.parent / "input_config.json"
    cfg = json.loads(cfgfile.read_text())
    geo = cfg["geometry"]
    pitch = float(geo["wire_pitch_mm"])
    half = int(geo["half_wires"])
    radius = 0.0005 * float(geo["wire_diameter_um"])
    wire_u = [i * pitch for i in range(-half, half + 1)]
    line = next((s for s in args.audit.read_text().splitlines()
                 if s.startswith("first_bad_tags=")), None)
    if not line:
        raise RuntimeError("Missing first_bad_tags in audit report")
    tags = ast.literal_eval(line.partition("=")[2])
    if not isinstance(tags, list) or not tags:
        raise RuntimeError("No listed bad elements; no localization necessary")
    print(f"Mesh: {gfile}")
    print(f"Anode axis u (mm): {wire_u}, v=0; radius={radius*1000:.3f} um")
    print(f"Diagnosing {len(tags)} bad-element IDs from previous audit")
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.open(str(gfile))
        print("Tag       minSJ     minDetJac       u_mm      v_mm      w_mm"
              "     wire_u_mm  center_r_um  min_node_r_um  max_node_r_um  bbox_u_mm  bbox_v_mm  bbox_w_mm")
        for tag in tags:
            etype, node_ids, _, _ = gmsh.model.mesh.getElement(int(tag))
            if int(etype) != 11 or len(node_ids) != 10:
                raise RuntimeError(f"Element {tag} is not Tet10 (type={etype})")
            points = [gmsh.model.mesh.getNode(int(n))[0] for n in node_ids]
            core = points[:4]  # Tet10 first four nodes are vertices
            c = [sum(p[k] for p in core)/4 for k in range(3)]
            uw = min(wire_u, key=lambda u: math.hypot(c[0]-u,c[1]))
            radii = [math.hypot(p[0]-uw, p[1]) for p in points]
            ranges = [max(p[k] for p in points)-min(p[k] for p in points)
                      for k in range(3)]
            sj = float(gmsh.model.mesh.getElementQualities([int(tag)], "minSJ")[0])
            det = float(gmsh.model.mesh.getElementQualities([int(tag)], "minDetJac")[0])
            r = math.hypot(c[0]-uw,c[1])
            print(f"{tag:8} {sj:9.5f} {det:13.5g}"
                  f" {c[0]:10.5f} {c[1]:10.5f} {c[2]:10.5f}"
                  f" {uw:12.5f} {r*1000:12.3f} {min(radii)*1000:14.3f}"
                  f" {max(radii)*1000:14.3f}"
                  f" {ranges[0]:10.5f} {ranges[1]:10.5f} {ranges[2]:10.5f}")
        print("Note: r values are distances to nearest wire axis; radius of metal wire = "
              f"{radius*1000:.3f} um.")
        print("These geometric locations are diagnostic, not an electrostatic certification.")
    finally:
        gmsh.finalize()

if __name__ == "__main__":
    main()
