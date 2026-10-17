#!/usr/bin/env python3
"""Physical geometry-only control for PI wire maps (no Garfield/Heed rerun).

Revisit all stored GEANT->checkerboard gas crossings to construct:
  - exact finite-wire CELL intersections per physical anode, as a control
    for the separate previously-computed collected-electron wire occupancy
  - XY muon gas-crossing density on a fixed 5-cm grid, per channel/station.

These are reference-plane extrapolations, not native Geant4 volume hits.
This control is deterministic, counts an entire crossed cell at most once,
and keeps overlapping module crossings separate in each det_id.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import time
from pathlib import Path

import numpy as np

from mwpc_checkerboard_projection import (
    SOURCES, GAS_X_CM, GAS_Y_CM, MS0_GAS_Y_CM, build_chambers,
)
from pi_fixed_wire_occupancy import LOOKUP, n_wires, nominal_wire_index

BIN_CM = 5.0


def bounds(st):
    chamber = build_chambers()[st]
    hx = GAS_X_CM/2.
    hy = (MS0_GAS_Y_CM if st == 0 else GAS_Y_CM)/2.
    xmin = math.floor((min(c.x for c in chamber)-hx-4.)/BIN_CM)*BIN_CM
    xmax = math.ceil((max(c.x for c in chamber)+hx+4.)/BIN_CM)*BIN_CM
    ymin = math.floor((min(c.y for c in chamber)-hy-4.)/BIN_CM)*BIN_CM
    ymax = math.ceil((max(c.y for c in chamber)+hy+4.)/BIN_CM)*BIN_CM
    return xmin, xmax, ymin, ymax


def cell_indices(row, station):
    """Finite wire cells intersected by projected segment in chamber local Y."""
    cy = float(row["chamber_center_y_cm"])
    y0, y1 = float(row["y_in_cm"]), float(row["y_out_cm"])
    a = int(nominal_wire_index(10*(y0-cy), station, 0.0))
    b = int(nominal_wire_index(10*(y1-cy), station, 0.0))
    lo = max(0, min(a, b))
    hi = min(n_wires(station)-1, max(a, b))
    return range(lo, hi+1) if hi >= lo else range(0)


def load_expected(path):
    expected = {(source, d): None for source in SOURCES for d in LOOKUP}
    with path.open(newline="", encoding="utf-8") as stream:
        rd = csv.DictReader(stream)
        for row in rd:
            if abs(float(row["scenario_ms2_z_cm"])-530) > 1e-5:
                continue
            src, d = row["source"], int(row["det_id"])
            if (src, d) not in expected:
                raise ValueError(f"Unknown chamber {src} {d}")
            if expected[src, d] is not None:
                raise ValueError(f"Duplicate module row: {src} {d}")
            expected[src, d] = int(row["muon_crossings"])
    if any(v is None for v in expected.values()):
        raise ValueError("Incomplete expected module table")
    return expected


def produce(crossings, out, max_rows=0):
    out.mkdir(parents=True, exist_ok=True)
    expected = load_expected(crossings/"modules.csv")
    counts = np.zeros((3, 298, 128), dtype=np.int64)
    exposures = np.zeros((3, 298), dtype=np.int64)
    grids = {}
    edges = {}
    for isrc in range(3):
        for st in range(6):
            xmin,xmax,ymin,ymax = bounds(st)
            xedges=np.arange(xmin,xmax+BIN_CM*.5,BIN_CM)
            yedges=np.arange(ymin,ymax+BIN_CM*.5,BIN_CM)
            edges[isrc,st]=(xedges,yedges)
            grids[isrc,st]=np.zeros((len(yedges)-1,len(xedges)-1),dtype=np.int64)
    start=time.monotonic()
    for isrc,source in enumerate(SOURCES):
        input_path=crossings/f"crossings_{source}_ms2_530.csv"
        if not input_path.is_file():raise FileNotFoundError(input_path)
        with input_path.open(newline="",encoding="utf-8") as f:
            rd=csv.DictReader(f)
            required={"source","station","det_id","chamber_center_y_cm",
                      "y_in_cm","y_out_cm","x_in_cm","x_out_cm"}
            if not required.issubset(rd.fieldnames or []):
                raise ValueError(f"Missing crossing fields: {required-set(rd.fieldnames or [])}")
            for i,row in enumerate(rd,1):
                det=int(row["det_id"]);st=int(row["station"])
                if row["source"]!=source or det not in LOOKUP or LOOKUP[det].station!=st:
                    raise ValueError(f"Source/station/detector mismatch: {row}")
                if abs(float(row["chamber_center_y_cm"])-LOOKUP[det].y)>1e-4:
                    raise ValueError(f"Wrong chamber Y centre det={det}")
                exposures[isrc,det]+=1
                for w in cell_indices(row,st):
                    counts[isrc,det,w]+=1
                x=.5*(float(row["x_in_cm"])+float(row["x_out_cm"]))
                y=.5*(float(row["y_in_cm"])+float(row["y_out_cm"]))
                xe,ye=edges[isrc,st]
                ix=int(math.floor((x-xe[0])/BIN_CM))
                iy=int(math.floor((y-ye[0])/BIN_CM))
                if not (0<=ix<grids[isrc,st].shape[1] and
                        0<=iy<grids[isrc,st].shape[0]):
                    raise ValueError(f"Out of station XY grid source={source} MS{st}: {x},{y}")
                grids[isrc,st][iy,ix]+=1
                if max_rows and i>=max_rows: break
                if i%250000==0:
                    print(f"  {source}: {i:,} gas crossings after {(time.monotonic()-start)/60:.1f} min",flush=True)
        if not max_rows:
            for det in LOOKUP:
                if exposures[isrc,det]!=expected[source,det]:
                    raise ValueError(f"Bad chamber count {source} DET{det}: "
                                     f"{exposures[isrc,det]} vs {expected[source,det]}")
            print(f"PASS: {source}, all 298 module counts agree",flush=True)
    for isrc in range(3):
        for st in range(6):
            expected_sum=int(sum(exposures[isrc,c.det_id] for c in build_chambers()[st]))
            if int(grids[isrc,st].sum())!=expected_sum:
                raise ValueError(f"XY histogram did not conserve crossings {(isrc,st)}")

    with (out/"pi_geom_wire_counts.csv").open("w",newline="",encoding="utf-8") as f:
        fields=("source","station","det_id","wire_index","N_crossings",
                "n_geom","P_geom")
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader()
        for isrc,src in enumerate(SOURCES):
            for det,ch in sorted(LOOKUP.items()):
                n=int(exposures[isrc,det])
                for j in range(n_wires(ch.station)):
                    v=int(counts[isrc,det,j])
                    if v>n:raise AssertionError("One geometric cell counted >1x per crossing")
                    w.writerow({"source":src,"station":ch.station,"det_id":det,
                                "wire_index":j,"N_crossings":n,"n_geom":v,
                                "P_geom":v/n if n else math.nan})
    arr={}
    for isrc,src in enumerate(SOURCES):
        for st in range(6):
            key=f"{src}_MS{st}"
            xe,ye=edges[isrc,st]
            arr[key]=grids[isrc,st]
            arr[key+"_xedges"]=xe
            arr[key+"_yedges"]=ye
    np.savez_compressed(out/"pi_crossing_density_xy.npz",**arr)
    (out/"pi_geom_wire_settings.json").write_text(json.dumps({
        "crossings":str(crossings.resolve()),
        "total_crossings_per_source":[int(exposures[i].sum()) for i in range(3)],
        "grid_bin_cm":BIN_CM, "nominal_wire_grid_phase_mm":0,
        "mode":"full" if max_rows==0 else "smoke",
        "scope":"GEANT reference-plane straight-line chamber gas projections",
        "geom_definition":"wire cell intersects finite y segment, no Heed or Garfield",
        "xy_definition":"midpoint of projected chamber gas segment, one entry per det_id",
        "rate_denominator":"sum of all chamber crossings in the station and source",
    },indent=2)+"\n",encoding="utf-8")
    print(f"Geometric wire control and XY crossing density saved in {out}",flush=True)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--crossings",type=Path,default=Path("test_runs/mwpc_stage_c/checkerboard_full"))
    ap.add_argument("--out",type=Path,default=Path("test_runs/mwpc_stage_c/pi_fixed_wire_maps/physical_validation"))
    ap.add_argument("--max-rows",type=int,default=0)
    args=ap.parse_args()
    produce(args.crossings,args.out,args.max_rows)


if __name__=="__main__":
    main()
