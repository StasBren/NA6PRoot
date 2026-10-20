#!/usr/bin/env python3
"""PI 1D TRUE time-sliced Garfield electron cloud (u(t),v(t),w(t)).

Uses REAL microscopic Garfield drift line points for electron origins
sampled along actual native Geant4 muon tracks by signed geometry angle.
Positions at fixed snapshot times are interpolated ONLY BETWEEN SUCCESSIVE
RECORDED MICROSCOPIC DRIFT POINTS (collision stride, not birth->endpoint).
At t>final transport time, electron is absent from surviving cloud.
Outputs stratified P(position | alive at t, native muon angle).
This is a limited angle-stratified pilot, NOT a complete event electron
response or a full Geant4/Heed/Garfield response across 200k parent events.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from collections import defaultdict
import numpy as np

from pi_centered_native_drift_births import OUT, SOURCES, CATEGORIES
from pi_centered_electron_sampler import GEOMETRY_ID

HERE=Path(__file__).resolve().parent
TEMPLATE=HERE/"pi_centered_true_1d_cloud_template.html"
TIMES=(0.,5.,10.,20.,40.,60.,80.,100.,150.)
EDGES={
 "u_relative_wire_mm":np.arange(-16.,16.0001,.5),
 "v_mm":np.linspace(-3.1,3.1,32),
 "w_displacement_mm":np.arange(-15.,15.0001,.5),
 "w_along_wire_mm":np.arange(-450.,450.0001,15.),
}


def load_micro(path):
    if not path.is_file():raise FileNotFoundError(path)
    lines=defaultdict(list);expected_n={};statuses={}
    with path.open(newline='',encoding='utf-8') as f:
        reader=csv.DictReader(f)
        required={"electron_id","point_index","t_rel_birth_ns",
                  "u_cell_mm","v_mm","w_mm","end_status","number_points"}
        if not required.issubset(reader.fieldnames or ()):
            raise ValueError("Incomplete REAL Garfield drift-line point schema")
        for row in reader:
            key=int(row["electron_id"]);ip=int(row["point_index"])
            n=int(row["number_points"])
            if n<1 or ip<0 or key<0:
                raise ValueError("Invalid native microscopic drift-line identifiers")
            if key in expected_n and (expected_n[key]!=n or
                   statuses[key]!=int(row["end_status"])):
                raise ValueError("Inconsistent Garfield end status/point count")
            expected_n[key]=n;statuses[key]=int(row["end_status"])
            values=tuple(float(row[k]) for k in (
                "t_rel_birth_ns","u_cell_mm","v_mm","w_mm"))
            if not all(math.isfinite(x) for x in values):
                raise ValueError("Nonfinite microscopic Garfield drift-line point")
            lines[key].append((ip,*values))
    for key,points in lines.items():
        if len(points)!=expected_n[key]:
            raise ValueError(f"Missing microscopic drift points, e={key}")
        points.sort()
        if any(points[i][0]!=i for i in range(len(points))):
            raise ValueError(f"Repeated/skipped Garfield collision point, e={key}")
        if abs(points[0][1])>1.e-5:
            raise ValueError("Garfield drift line does not start at electron birth t=0")
        if any(points[i][1]<points[i-1][1]-1.e-8
               for i in range(1,len(points))):
            raise ValueError("Garfield trace time not monotonically increasing")
    if not lines:
        raise ValueError("No microscopic Garfield trajectories")
    return lines,statuses


def position_at(points,t):
    """Only interpolate two neighboring real microscopic Garfield points."""
    if t>points[-1][1]+1.e-10:
        return None
    if t<=points[0][1]:
        return points[0][2:5]
    times=[p[1] for p in points]
    j=int(np.searchsorted(times,t,side='right'))
    if j>=len(points):return points[-1][2:5]
    a=points[j-1];b=points[j]
    dt=b[1]-a[1]
    q=(t-a[1])/dt if dt>1.e-12 else 0.
    return tuple(a[k]+(b[k]-a[k])*q for k in (2,3,4))


def analyze(metadata,tracks,manifest):
    groups=defaultdict(lambda:{
      "born":0,"alive":[0]*len(TIMES),
      "hist":{k:np.zeros((len(TIMES),len(v)-1),dtype=np.int64)
              for k,v in EDGES.items()},
    })
    miss=0;outside=defaultdict(int)
    for row in metadata:
        i=int(row["electron_id"])
        if i not in tracks:
            miss+=1;continue
        src=row["source"];st=int(row["station"]);g=int(row["angle_group"])
        if src not in SOURCES or not 0<=st<6 or not 0<=g<len(CATEGORIES):
            raise ValueError("Unknown native gas-hit angle/source/station label")
        k=f"{src}|{st}|{g}"
        bucket=groups[k];bucket["born"]+=1
        wire_center=float(row["wire_center_u_mm"])
        w_birth=float(row["w_physical_birth_mm"])
        points=tracks[i]
        for j,t in enumerate(TIMES):
            p=position_at(points,t)
            if p is None:continue
            u,v,w=p
            bucket["alive"][j]+=1
            coords={
              "u_relative_wire_mm":u,
              "v_mm":v,
              "w_displacement_mm":w,
              "w_along_wire_mm":w_birth+w,
            }
            for field,value in coords.items():
                edges=EDGES[field]
                if not edges[0]<=value<=edges[-1]:
                    outside[field]+=1
                    continue
                hist=bucket["hist"][field]
                bin_index=min(len(edges)-2,int(np.searchsorted(
                    edges,value,side='right')-1))
                hist[j,bin_index]+=1
    if miss or len(tracks)!=len(metadata):
        raise ValueError(f"Microscopic drift metadata not complete "
                         f"(missing={miss}, drift={len(tracks)}, meta={len(metadata)})")
    if outside:
        raise ValueError("Microscopic drift outside prepared physically "
                         f"declared axes: {dict(outside)}. Widen, never clamp.")
    out={}
    for k,bucket in groups.items():
        if any(int(arr[i].sum())!=bucket["alive"][i]
               for arr in bucket["hist"].values()
               for i in range(len(TIMES))):
            raise AssertionError("True microscopic 1D histogram conservation failed")
        if any(bucket["alive"][i]>bucket["alive"][i-1]
               for i in range(1,len(TIMES))):
            raise AssertionError("An electron came back after being collected/stopped")
        out[k]={
          "born":bucket["born"],"alive":bucket["alive"],
          "hist":{name:arr.tolist() for name,arr in bucket["hist"].items()},
        }
    return out


def build(args):
    paths=[args.metadata,args.trace,args.starts_manifest,
           args.map_manifest,args.template]
    for p in paths:
        if not p.is_file():raise FileNotFoundError(p)
    starts=json.loads(args.starts_manifest.read_text(encoding='utf-8'))
    model=json.loads(args.map_manifest.read_text(encoding='utf-8'))
    if starts.get("geometry_id")!=GEOMETRY_ID or (
        model.get("geometry_id")!=GEOMETRY_ID
    ) or model.get('provenance_log_checked') is not True:
        raise ValueError("Refusing non-centred 3+3mm or unverified Garfield map")
    if hashlib.sha256(args.metadata.read_bytes()).hexdigest()!=starts['metadata_sha256']:
        raise ValueError("Original native Heed start provenance SHA mismatch")
    with args.metadata.open(newline='',encoding='utf-8') as f:
        metadata=list(csv.DictReader(f))
    lines,statuses=load_micro(args.trace)
    if len(metadata)!=starts["generated_native_Heed_electron_starts"]:
        raise ValueError("Native electron birth sample count changed")
    groups=analyze(metadata,lines,starts)
    data={
      "sources":list(SOURCES),
      "angleCategories":[c[2] for c in CATEGORIES],
      "times_ns":list(TIMES),
      "edges":{k:v.tolist() for k,v in EDGES.items()},
      "groups":groups,
      "meta":{
        "geometry_id":GEOMETRY_ID,
        "scope":"Angle-stratified subset of true native Geant4 first-chunk "
                "muon gas segments, with newly re-simulated REAL Garfield "
                "microscopic electron trajectories, not full campaign",
        "method":"Linear interpolation only between consecutively STORED "
                 "Garfield microscopic collision subsample points; NO "
                 "birth-to-endpoint interpolation",
        "normalization":"P(position_bin | electron still drifting at "
                         "snapshot time, selected muon entrance angle)",
        "disappearing":"Collected/lost electrons removed from surviving "
                       "cloud for times after their final drift point",
        "magnetic_field_T":0.,
        "no_electronics_or_cathodes":True,
      },
    }
    template=args.template.read_text(encoding='utf-8')
    marker="__PI_TRUE_DRIFT_1D_DATA__"
    if template.count(marker)!=1:
        raise ValueError("True 1D Garfield template marker mismatch")
    args.html.parent.mkdir(parents=True,exist_ok=True)
    args.html.write_text(template.replace(marker,json.dumps(
        data,separators=(',',':'),ensure_ascii=False).replace('</','<\\/')),
        encoding='utf-8')
    report={
      "geometry_id":GEOMETRY_ID,
      "number_traced_electrons":len(lines),
      "source_native_Heed_metadata_sha256":starts['metadata_sha256'],
      "true_microscopic_Garfield_drift_points_sha256":
          hashlib.sha256(args.trace.read_bytes()).hexdigest(),
      "selected_angular_groups":len(groups),
      "snapshot_times_ns":list(TIMES),
      "no_original_endpoints_relabelled":True,
      "interpolation":"Only between adjacent microscopic recorded points",
      "not_full_200k_Geant4":True,
    }
    (args.html.parent/"pi_centered_true_1d_cloud_manifest.json").write_text(
        json.dumps(report,indent=2,ensure_ascii=False)+"\n",encoding='utf-8')
    print(f"PASS: {len(lines):,} new microscopic Garfield electron trajectories"
          f" grouped into {len(groups)} native muon angle×station selections")
    print("True time-resolved PI 1D u(t),v(t),w(t) pilot:",args.html)
    return report


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--metadata',type=Path,
                   default=OUT/"native_heed_start_metadata.csv")
    p.add_argument('--trace',type=Path,
                   default=OUT/"centered_native_microscopic_drift_points.csv")
    p.add_argument('--starts-manifest',type=Path,
                   default=OUT/"native_heed_trace_manifest.json")
    p.add_argument('--map-manifest',type=Path,
                   default=OUT.parent/"c_transport_box_map_gm3p0_gp3p0"/
                           "pi_centered_transport_manifest.json")
    p.add_argument('--html',type=Path,
                   default=OUT.parent/"pi_fixed_wire_maps"/
                           "pi_centered_true_1d_cloud.html")
    p.add_argument('--template',type=Path,default=TEMPLATE)
    args=p.parse_args()
    build(args)


if __name__=='__main__':
    main()
