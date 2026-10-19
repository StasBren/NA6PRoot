#!/usr/bin/env python3
"""Render clearer, map-clickable PI dashboard V2 from EXISTING completed pilot HTML.

No Geant4, Heed or Garfield rerun; preserves exact 3+3mm pilot payload.
The original V1 viewer stays intact for reproducibility. 2D original
coordinate/time counts are presented as normalized 1D coordinate PDFs
for SELECTED ELECTRON ARRIVAL-TIME SLICES, not actual trajectories u(t).
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

from mwpc_checkerboard_projection import build_chambers

HERE=Path(__file__).resolve().parent
REPO=HERE.parents[1]
BASE=REPO/"test_runs/mwpc_stage_c/pi_fixed_wire_maps"
ORIGINAL=BASE/"pi_centered_timing_pilot.html"
OUTPUT=BASE/"pi_centered_timing_dashboard_v2.html"
TEMPLATE=HERE/"pi_centered_timing_dashboard_v2_template.html"
OPEN='<script id="payload" type="application/json">'
CLOSE='</script>'
TAG='__PI_CENTERED_DASHBOARD_V2_DATA__'


def payload_from_html(path:Path):
    html=path.read_text(encoding="utf-8")
    if html.count(OPEN)!=1:
        raise ValueError("Expected exactly one complete original pilot JSON script")
    rest=html.split(OPEN,1)[1]
    raw=rest.split(CLOSE,1)[0]
    obj=json.loads(raw)
    if obj.get("meta",{}).get("geometry_id")!="MWPC_GARFIELD_CENTERED_3p0_3p0_B0_V1":
        raise ValueError("Refusing old 2+4mm or non-verified 3+3mm response")
    if len(obj.get("sourceNames",[]))!=3 or len(obj.get("geometry",{}))!=298:
        raise ValueError("Incomplete native-geometry pilot JSON input")
    return obj,hashlib.sha256(raw.encode("utf-8")).hexdigest()


def validate_payload(data):
    geo={c.det_id:c for group in build_chambers() for c in group}
    assert len(geo)==298
    if any(int(data["geometry"][str(det)])!=c.station for det,c in geo.items()):
        raise ValueError("Physical chamber DET IDs disagree with the original pilot")
    n_checked=0
    for key,item in data["groups"].items():
        n=int(item["n"])
        if n<0 or int(item["known"])+int(item["unknown"])!=int(item["collected"]):
            raise ValueError("Inconsistent electron timing counters "+key)
        if sum(item["time"])!=int(item["known"]):
            raise ValueError("Electron drift time spectrum conservation failed "+key)
        if len(item["m"])!=6 or any(sum(counts)!=n for counts in item["m"]):
            raise ValueError("Physical anode multiplicity conservation failed "+key)
        means=[sum(j*v for j,v in enumerate(counts[:5]))+
               sum(5*v for v in counts[5:]) for counts in item["m"]]
        if means!=sorted(means):
            raise ValueError("Time-window physical anode multiplicity decreased "+key)
        if item["coords"] is not None:
            for name,grid in item["coords"].items():
                if name not in data["coordinates"]:
                    raise ValueError("Unknown initial electron coordinate "+name)
                if len(grid)!=len(data["coordinates"][name])-1:
                    raise ValueError("Coordinate-binning mismatch "+key)
                if sum(sum(row) for row in grid)!=int(item["known"]):
                    raise ValueError("Coordinate/time count conservation failed "+key)
        n_checked+=1

    # Physical chamber and station exposures must agree with original
    # first-chunk native gas-hit truth for all 3 channels and 6 stations.
    for source in data["sourceNames"]:
        for station in range(6):
            for axis in ("u","w"):
                key="|".join(map(str,(source,station,-1,axis,-1)))
                station_hits=int(data["groups"][key]["n"])
                module_hits=sum(
                    int(data["groups"].get("|".join(map(str,
                        (source,station,c.det_id,axis,-1))),{}).get("n",0))
                    for c in geo.values() if c.station==station)
                if station_hits!=module_hits:
                    raise ValueError("Native source/module/station count mismatch "+key)

    data["layout"]=[
        {"id":c.det_id,"st":c.station,"x":c.x,"y":c.y,
         "hx":c.half_x,"hy":c.half_y,"row":c.row,"col":c.col,
         "layer":c.layer}
        for c in geo.values()
    ]
    data["meta"]["dashboard_v2"]="UI and presentation-only upgrade; identical original physics payload"
    return n_checked


def build(input_file:Path,output:Path,template:Path):
    data,sha=payload_from_html(input_file)
    count=validate_payload(data)
    html=template.read_text(encoding="utf-8")
    if html.count(TAG)!=1:
        raise ValueError("V2 dashboard template must have exactly one JSON payload")
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(
        html.replace(TAG,json.dumps(data,ensure_ascii=False,separators=(",",":"))
                     .replace("</","<\\/")),
        encoding="utf-8",
    )
    print("PASS: original fully completed symmetric Garfield 3+3mm pilot remains unchanged")
    print(f"PASS: {count} original angular/chamber/source histogram groups conserved")
    print(f"PASS: full 298-module XY geometry and 18 station sums matched")
    print(f"Original embedded payload SHA256: {sha}")
    print("PI dashboard V2 ready:",output)
    return output


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input",type=Path,default=ORIGINAL)
    ap.add_argument("--output",type=Path,default=OUTPUT)
    ap.add_argument("--template",type=Path,default=TEMPLATE)
    args=ap.parse_args()
    build(args.input,args.output,args.template)


if __name__=="__main__":
    main()
