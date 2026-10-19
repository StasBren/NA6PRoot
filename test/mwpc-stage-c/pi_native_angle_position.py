#!/usr/bin/env python3
"""Interactive full Geant4 native angular PDFs vs each MWPC module's XY location.

Uses completed read-only native PI QA CSV; NO ROOT/Garfield/Heed and no
transport rerun. Maps signed-theta probability fractions at EXACT saved
histogram-bin boundaries. No fabricated means or sub-bin angle resolution.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from mwpc_checkerboard_projection import SOURCES, build_chambers
from pi_native_full_campaign_qa import ANGLE

HERE=Path(__file__).resolve().parent
REPO=HERE.parents[1]
DEFAULT_IN=REPO/"test_runs/mwpc_stage_c/pi_native_full_quality"
DEFAULT_OUT=REPO/"test_runs/mwpc_stage_c/pi_fixed_wire_maps/pi_native_angle_position.html"
GEOMETRY=build_chambers()
FIELDS=("theta_u_in_deg","theta_w_in_deg",
        "theta_u_out_deg","theta_w_out_deg")
NBIN=len(ANGLE)-1


def build(input_dir,out,template):
    qa=json.loads((input_dir/"pi_native_all_QA_report.json").read_text(encoding="utf-8"))
    if not qa["full_campaign_analyzed"] or qa["processed_chunks_per_source"]!=67:
        raise ValueError("Native angular spatial dashboard requires full Geant4 QA")
    exposure={(source,c.det_id):None for source in SOURCES for group in GEOMETRY for c in group}
    with (input_dir/"pi_native_all_chamber_quality.csv").open(newline="") as stream:
        for row in csv.DictReader(stream):
            key=row["source"],int(row["det_id"])
            if key not in exposure or exposure[key] is not None:
                raise ValueError("Invalid/duplicate native module exposure "+str(key))
            exposure[key]=int(row["gas_hit_rows"])
    if any(v is None for v in exposure.values()):
        raise ValueError("Missing 894 native module exposures")
    station_exposure={}
    with (input_dir/"pi_native_all_station_quality.csv").open(newline="") as stream:
        for row in csv.DictReader(stream):
            key=row["source"],int(row["station"])
            if key in station_exposure or key[0] not in SOURCES or key[1] not in range(6):
                raise ValueError("Invalid/duplicate source×station native QA row")
            station_exposure[key]=int(row["gas_hit_rows"])
    if len(station_exposure)!=len(SOURCES)*6:
        raise ValueError("Missing full 18 source×station QA exposures")
    data={src:{str(ch.det_id):{
        "st":ch.station,"n":exposure[src,ch.det_id],
        "angles":{name:[None]*NBIN for name in FIELDS}
    } for station in GEOMETRY for ch in station} for src in SOURCES}
    with (input_dir/"pi_native_all_histograms.csv").open(newline="") as stream:
        for row in csv.DictReader(stream):
            field=row["variable"]
            if field not in FIELDS:
                continue
            src=row["source"];det=int(row["det_id"])
            if src not in data or str(det) not in data[src]:
                raise ValueError(f"Unrecognized native angle module: {src} {det}")
            bucket=data[src][str(det)]
            n=int(row["gas_hit_exposure"])
            if n!=bucket["n"]:
                raise ValueError("Native angle histogram exposure mismatch")
            lo=float(row["bin_lower"]);hi=float(row["bin_upper"])
            try:
                index=next(i for i in range(NBIN)
                           if abs(ANGLE[i]-lo)<1e-8 and abs(ANGLE[i+1]-hi)<1e-8)
            except StopIteration:
                raise ValueError(f"Invalid native ANGLE bin {lo}..{hi}") from None
            if bucket["angles"][field][index] is not None:
                raise ValueError("Repeated module/angle histogram bin")
            bucket["angles"][field][index]=int(row["hit_count"])
    for src in SOURCES:
        for group in GEOMETRY:
            for c in group:
                item=data[src][str(c.det_id)]
                for name in FIELDS:
                    counts=item["angles"][name]
                    if any(x is None for x in counts):
                        raise ValueError(f"Missing angular PDF {src} DET{c.det_id} {name}")
                    if sum(counts)!=item["n"]:
                        raise ValueError(f"Angular PDF exposure mismatch {src} DET{c.det_id} {name}")
                if item["n"]<0:
                    raise ValueError("Negative module native gas-hit exposure")
        for st in range(6):
            n=sum(data[src][str(c.det_id)]["n"] for c in GEOMETRY[st])
            expected=station_exposure[src,st]
            if n!=expected:
                raise ValueError(
                    f"Actual native angle module count != station QA "
                    f"{src} MS{st}: modules={n}, station={expected}"
                )
    payload={
        "sources":SOURCES,
        "edges":ANGLE,
        "fields":FIELDS,
        "geometry":[{
            "id":c.det_id,"st":c.station,
            "x":c.x,"y":c.y,"row":c.row,"col":c.col,"layer":c.layer,
            "hx":c.half_x,"hy":c.half_y,
        } for station in GEOMETRY for c in station],
        "module":data,
        "meta":{
            "scope":"FULL 201-chunk native Geant4, 200k events/channel",
            "denominator":"Only muons with TRUE Geant4 gas-hit segment in selected module",
            "angles":"Signed detector TGeo-local angles, theta_u=atan2(pu,pv), theta_w=atan2(pw,pv)",
            "statistics":"Exact bin-integrated probabilities at 0, +-2, +-10 deg; original coarsely binned per-module angle histograms, NOT exact mean or median",
            "not":"NOT Heed/Garfield drift, NOT primary-electron time, NOT the number of incident station muons",
        },
    }
    html=template.read_text(encoding="utf-8")
    marker="__ANGLE_POSITION_DATA__"
    if html.count(marker)!=1:
        raise ValueError("Exactly one angle payload marker is required")
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(html.replace(marker,json.dumps(payload,ensure_ascii=False,
           separators=(",",":")).replace("</","<\\/")),encoding="utf-8")
    print(f"PASS: native Geant4 full 201-chunk angular maps: "
          f"{len(SOURCES)*sum(map(len,GEOMETRY))} module/source records, "
          f"4 validated signed-angle histograms each")
    print("Interactive full-run angle-vs-chamber XY map:",out)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input",type=Path,default=DEFAULT_IN)
    ap.add_argument("--out",type=Path,default=DEFAULT_OUT)
    ap.add_argument("--template",type=Path,default=HERE/"pi_native_angle_position_template.html")
    args=ap.parse_args()
    build(args.input,args.out,args.template)


if __name__=="__main__":
    main()
