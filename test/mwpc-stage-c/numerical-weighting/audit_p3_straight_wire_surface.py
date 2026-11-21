#!/usr/bin/env python3
"""Read-only as-built cylinder-boundary audit on existing Gmsh .msh.

For a straight-sided Tet10 mesh, the physical wire boundary consists of
planar Tri6 faces. Report observed radii of all 3 corner/3 mid-edge nodes
and sampled centroids for the physical group "anode_wires".
These are geometrical diagnostics, NOT field accuracy certification.
"""
from __future__ import annotations
import argparse
import json
import math
from pathlib import Path

def outstat(label, vals):
    if not vals:
        raise RuntimeError("No samples for " + label)
    print(
        f"{label}: count={len(vals)} min={min(vals):+.8f} um "
        f"max={max(vals):+.8f} um "
        f"max_abs={max(map(abs, vals)):.8f} um"
    )

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mesh", type=Path, help="Existing chamber.msh")
    args = ap.parse_args()
    msh = args.mesh.resolve()
    config_file = msh.parent.parent / "input_config.json"
    cfg = json.loads(config_file.read_text(encoding="utf-8"))
    geo = cfg["geometry"]
    nominal = float(geo["wire_diameter_um"]) * .0005
    pitch = float(geo["wire_pitch_mm"])
    half_wires = int(geo["half_wires"])
    centres = [i * pitch for i in range(-half_wires,half_wires+1)]
    import gmsh
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.open(str(msh))
        groups = {gmsh.model.getPhysicalName(d,t):t
                  for d,t in gmsh.model.getPhysicalGroups(dim=2)}
        if "anode_wires" not in groups:
            raise RuntimeError("No physical anode_wires group. Available: "+str(groups))
        surfaces = gmsh.model.getEntitiesForPhysicalGroup(2, groups["anode_wires"])
        corner_tags=set()
        mid_tags=set()
        all_triangles=[]
        for surf in surfaces:
            types, _, conn_sets = gmsh.model.mesh.getElements(dim=2, tag=int(surf))
            for element_type, conn_flat in zip(types, conn_sets):
                if int(element_type) != 9:
                    raise RuntimeError(f"Non-Tri6 element on wire surface: type {element_type}")
                if len(conn_flat) % 6:
                    raise RuntimeError("Invalid Tri6 connectivity")
                for k in range(0,len(conn_flat),6):
                    nodes = tuple(map(int,conn_flat[k:k+6]))
                    corner_tags.update(nodes[:3])
                    mid_tags.update(nodes[3:])
                    all_triangles.append(nodes)
        if not all_triangles:
            raise RuntimeError("No wire-surface triangles")
        wanted=corner_tags|mid_tags
        coords={}
        for surf in surfaces:
            tags,xyz,_=gmsh.model.mesh.getNodes(dim=2,tag=int(surf),
                                                includeBoundary=True)
            for k, tag in enumerate(tags):
                key=int(tag)
                if key in wanted:
                    coords[key]=(float(xyz[3*k]),float(xyz[3*k+1]),
                                 float(xyz[3*k+2]))
        for t in wanted-coords.keys():
            xyz,_,_,_=gmsh.model.mesh.getNode(t)
            coords[t]=(float(xyz[0]),float(xyz[1]),float(xyz[2]))
        def radius_error(xyz):
            u,v,_=xyz
            uw=min(centres,key=lambda w:abs(u-w))
            return (math.hypot(u-uw,v)-nominal)*1000.0
        cerrors=[radius_error(coords[t]) for t in corner_tags]
        merrors=[radius_error(coords[t]) for t in mid_tags]
        centerrors=[]
        max_dphi=0.
        for nodes in all_triangles:
            p=[coords[x] for x in nodes[:3]]
            c=tuple(sum(q[i] for q in p)/3.0 for i in range(3))
            centerrors.append(radius_error(c))
            for i,j in ((0,1),(1,2),(2,0)):
                q1,q2=p[i],p[j]
                umean=(q1[0]+q2[0])/2
                uw=min(centres,key=lambda w:abs(umean-w))
                a=math.atan2(q1[1],q1[0]-uw)
                b=math.atan2(q2[1],q2[0]-uw)
                d=abs(math.atan2(math.sin(a-b), math.cos(a-b)))
                max_dphi=max(max_dphi,d)
        print(f"Mesh: {msh}")
        print(f"Wires: {len(centres)}, nominal radius={nominal*1000:.4f} um")
        print(f"Physical wire Tri6 faces: {len(all_triangles)}")
        print(f"Surface corner nodes: {len(corner_tags)}, midside nodes: {len(mid_tags)}")
        outstat("Corner radial error",cerrors)
        outstat("Midside radial error",merrors)
        outstat("Triangle-centroid radial error (sample)",centerrors)
        print(f"Maximum observed triangle-edge azimuthal span: "
              f"{max_dphi*180/math.pi:.5f} degrees")
        print("This reports sampled boundary positions; not a rigorous maximum over "
              "all points of each face, nor field/Ramo convergence.")
    finally:
        gmsh.finalize()
if __name__=="__main__":
    main()
