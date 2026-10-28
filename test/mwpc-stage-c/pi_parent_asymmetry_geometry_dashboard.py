#!/usr/bin/env python3
"""Extend the ORIGINAL parent-event ROOT-style P(A) dashboard with p3/p4 LUT20.

Strictly keep parent-level definition across ALL selectable geometries:
sum N_up/N_down over all physical muon crossings of one (chunk,event)
in one MS station or one DET, THEN calculate A for that parent. Never
merge per-crossing histograms or compare mismatched statistical units.

The original p=3 LUT500 parent dashboard is read, not modified.
For p=3 and p=4 LUT20, reaggregate the complete existing 100k-parent
'full34/crossings.csv.gz' without rerunning Garfield, Geant4 or Heed.
Original sources/manifest SHA256 remain unchanged. Output is a NEW HTML.
"""
from __future__ import annotations
import argparse
import csv
import gzip
import hashlib
import json
import os
from collections import Counter,defaultdict
from pathlib import Path

from pi_lut500_event_asymmetry_root_dashboard import Distribution,Paired
from pi_pitch_fieldwire_pilot import tag

HERE=Path(__file__).resolve().parent
REPO=HERE.parents[1]
ROOT=REPO/"test_runs/mwpc_stage_c"
OLD=ROOT/"pi_lut500_parent_asymmetry_root_v1"
CROSS=ROOT/"pi_lut20_p4_muon_crossing_asymmetry_v1/full34"
TEMPLATE=HERE/"pi_lut500_event_asymmetry_root_template.html"
OUTPUT=OLD/"dashboard_geometry_lut_compare.html"
PITCHES=(3.,4.)
GEOMS=[(p,fw,tag(p,fw)) for p in PITCHES for fw in (0.,75.)]
SCHEMA="PI_JPSI_SIX_STATIONS_PER_DET_PARENT_ASYM_LUT500_AND_LUT20_V1"
RAW_FIELDS=("source","chunk","event","hit_index","track_id","station","det_id",
            "theta_u_deg","theta_w_deg","u_in_mm","u_out_mm","pitch_mm",
            "lut_trials_per_cell","field_wire_um","n_primary","n_raw_anode",
            "n_anode_finite","n_edge_lost","up","down","hemisphere_boundary",
            "top","bottom","side","A_hemisphere","A_top_bottom","f_top_bottom")

def digest(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for data in iter(lambda:f.read(1<<20),b""):h.update(data)
    return h.hexdigest()

def read_sources(old_dir,cross_dir):
    old_dir=Path(old_dir);cross_dir=Path(cross_dir)
    old=json.loads((old_dir/"summary.json").read_text())
    om=json.loads((old_dir/"manifest.json").read_text())
    source=json.loads((cross_dir/"summary.json").read_text())
    cm=json.loads((cross_dir/"manifest.json").read_text())
    if (old.get("schema")!="PI_JPSI_SIX_STATIONS_PER_DET_PARENT_ASYM_P3_LUT500_ROOT_V1"
        or old.get("pitch_mm")!=3. or old.get("parents_generated_in_input_chunks")!=100000
        or old.get("computed_stations")!=list(range(6))
        or old.get("partial_exposure") is not False
        or om.get("schema")!=old["schema"] or om.get("complete") is not True):
        raise ValueError("Old parent-level p3 LUT500 must be the verified 100k full source")
    hashes=om.get("files_sha256",{})
    if hashes.get("summary.json")!=digest(old_dir/"summary.json"):
        raise ValueError("Original LUT500 parent summary SHA256 mismatch")
    if (source.get("schema")!="PI_JPSI_NATIVE_MUON_GAS_CROSSING_ASYM_P3_P4_LUT20_V1"
        or source.get("source")!="jpsi40-bore"
        or source.get("lut_trials_per_cell")!=20
        or source.get("pitches_mm")!=[3.,4.]
        or source.get("stations_selected")!=list(range(6))
        or source.get("generated_parent_slots")!=100000
        or source.get("partial_exposure") is not False
        or cm.get("schema")!=source["schema"]
        or cm.get("complete_six_station") is not True):
        raise ValueError("Expected immutable complete 100k-parent p3/p4 crossing LUT20")
    ch=cm.get("sha256_outputs",{})
    for name in ("summary.json","crossings.csv.gz"):
        if ch.get(name)!=digest(cross_dir/name):
            raise ValueError("Frozen LUT20 full34 SHA256 mismatch: "+name)
    if cm.get("garfield_lut_sha256")!=source.get("LUT20"):
        raise ValueError("LUT20 frozen Garfield manifest mismatch")
    for st in range(6):
        saved_det=set(old["stations"][str(st)]["chambers"])
        source_det=set(source["stations"][str(st)]["det"])
        if not source_det.issubset(saved_det):
            raise ValueError(f"MS{st} contains a DET absent from p3 layout")
    return old,source

def aggregate_lut20(path,original,source):
    """O(Ncrossings) streaming reducer, O(one-parent) temporary event storage."""
    st_dist={p:{st:{g[2]:Distribution() for g in GEOMS if g[0]==p}
                for st in range(6)} for p in PITCHES}
    st_pair={p:{st:Paired() for st in range(6)} for p in PITCHES}
    det_dist={p:{st:{det:{g[2]:Distribution() for g in GEOMS if g[0]==p}
                     for det in original["stations"][str(st)]["chambers"]}
                for st in range(6)} for p in PITCHES}
    det_pair={p:{st:{det:Paired()
                    for det in original["stations"][str(st)]["chambers"]}
                for st in range(6)} for p in PITCHES}
    crossings=Counter()
    primaries=Counter()
    summed=defaultdict(lambda:[0,0])
    parents_seen=Counter()
    current_parent=None
    parent_station={}
    parent_det={}
    crossing_key=None
    seen=[]
    nrows=0

    def empty_counts():
        return [[0,0] for _ in GEOMS]

    def finish_parent():
        if current_parent is None:return
        for st,counts in parent_station.items():
            parents_seen[st]+=1
            for p in PITCHES:
                j=int(p==4.)*2
                st_dist[p][st][GEOMS[j][2]].add(*counts[j])
                st_dist[p][st][GEOMS[j+1][2]].add(*counts[j+1])
                st_pair[p][st].add([*counts[j],*counts[j+1]])
        for (st,det),counts in parent_det.items():
            for p in PITCHES:
                j=int(p==4.)*2
                det_dist[p][st][str(det)][GEOMS[j][2]].add(*counts[j])
                det_dist[p][st][str(det)][GEOMS[j+1][2]].add(*counts[j+1])
                det_pair[p][st][str(det)].add([*counts[j],*counts[j+1]])

    def finish_crossing():
        if crossing_key is not None and seen!=list(range(len(GEOMS))):
            raise ValueError(f"Missing/duplicated LUT20 geometry for gas crossing {crossing_key}: {seen}")

    with gzip.open(path,"rt",encoding="utf-8",newline="") as f:
        reader=csv.DictReader(f)
        if tuple(reader.fieldnames or ())!=RAW_FIELDS:
            raise ValueError("Unexpected full34 crossing CSV column schema")
        for row in reader:
            nrows+=1
            parent=(int(row["chunk"]),int(row["event"]))
            if current_parent!=parent:
                finish_crossing()
                if current_parent is not None and parent<=current_parent:
                    raise ValueError("Native source parent ids are not monotonically grouped")
                finish_parent()
                parent_station={};parent_det={}
                current_parent=parent
                crossing_key=None
                seen=[]
            st=int(row["station"]);det=int(row["det_id"])
            if st not in range(6) or str(det) not in det_dist[3.][st]:
                raise ValueError("Unexpected physical station/DET")
            key=(parent,int(row["hit_index"]),int(row["track_id"]),st,det)
            if crossing_key!=key:
                finish_crossing()
                crossing_key=key
                seen=[]
            if len(seen)>=len(GEOMS):
                raise ValueError("Duplicated original muon gas crossing")
            geom=(float(row["pitch_mm"]),float(row["field_wire_um"]),
                  tag(float(row["pitch_mm"]),float(row["field_wire_um"])))
            if geom!=GEOMS[len(seen)] or int(row["lut_trials_per_cell"])!=20:
                raise ValueError("Wrong LUT20 configuration ordering or trial count")
            j=len(seen)
            seen.append(j)
            u=int(row["up"]);d=int(row["down"])
            b=int(row["hemisphere_boundary"]);fin=int(row["n_anode_finite"])
            if min(u,d,b,fin)<0 or u+d+b!=fin or row["source"]!="jpsi40-bore":
                raise ValueError("Invalid parent or electron conservation")
            sc=parent_station.setdefault(st,empty_counts())
            dc=parent_det.setdefault((st,det),empty_counts())
            sc[j][0]+=u;sc[j][1]+=d
            dc[j][0]+=u;dc[j][1]+=d
            summed[(st,geom[2])][0]+=u
            summed[(st,geom[2])][1]+=d
            if j==0:
                crossings[st]+=1
                primaries[st]+=int(row["n_primary"])
    finish_crossing()
    finish_parent()
    expected=sum(int(source["stations"][str(st)]["gas_crossings"]) for st in range(6))
    if nrows!=expected*len(GEOMS):
        raise ValueError("Unexpected number of raw crossing/config records")
    for st in range(6):
        oldst=original["stations"][str(st)]
        crossed=int(source["stations"][str(st)]["gas_crossings"])
        if crossings[st]!=crossed or crossings[st]!=oldst["gas_crossings"]:
            raise ValueError(f"MS{st} gas crossing mismatch")
        if primaries[st]!=int(oldst["heed_primaries"]):
            raise ValueError(f"MS{st} Heed birth count mismatch")
        for _,_,cfg in GEOMS:
            s=source["stations"][str(st)]["configs"][cfg]
            if summed[(st,cfg)]!=[s["n_upper"],s["n_lower"]]:
                raise ValueError(f"MS{st}/{cfg} primary electron conservation failure")
        if parents_seen[st]!=oldst["station_parent"]["configurations"][
             original["configs"][0]["key"]]["parents_with_gas_crossing"]:
            raise ValueError(f"MS{st} parent-with-crossing count mismatch")
    result={}
    for p in PITCHES:
        per_station={}
        cfgs=[g[2] for g in GEOMS if g[0]==p]
        for st in range(6):
            base=original["stations"][str(st)]
            vals={cfg:st_dist[p][st][cfg].result() for cfg in cfgs}
            dets={}
            for det in base["chambers"]:
                dets[det]={
                    "configurations":{cfg:det_dist[p][st][det][cfg].result() for cfg in cfgs},
                    "paired":det_pair[p][st][det].result()}
            for cfg in cfgs:
                for field in ("total_up","total_down"):
                    if vals[cfg][field]!=sum(x["configurations"][cfg][field] for x in dets.values()):
                        raise ValueError("Physical DET-to-station electron mismatch")
            per_station[str(st)]={
                "gas_crossings":crossings[st],
                "heed_primaries":primaries[st],
                "layout":base["layout"],
                "station_parent":{"configurations":vals,"paired":st_pair[p][st].result()},
                "chambers":dets,
                "n_physical_chambers":base["n_physical_chambers"]}
        result[p]={
            "schema":SCHEMA,
            "source":"jpsi40-bore",
            "pitch_mm":p,
            "lut_trials_per_cell":20,
            "B_t":0,
            "number_of_source_chunks":2,
            "parents_generated_in_input_chunks":100000,
            "partial_exposure":False,
            "computed_stations":list(range(6)),
            "configs":[{"key":cfgs[0],"label":"No field wire",
                       "field_wire_diameter_um":0},
                       {"key":cfgs[1],"label":"Field wire 75 um",
                       "field_wire_diameter_um":75}],
            "asymmetry":original["asymmetry"],
            "stations":per_station}
    print(f"PASS LUT20 full34 parent regroup: {expected:,} muon crossings, "
          f"{sum(parents_seen.values()):,} station-level parent observations")
    return result

def build(old_dir=OLD,cross_dir=CROSS,out=OUTPUT,force=False):
    old_dir=Path(old_dir);cross_dir=Path(cross_dir);out=Path(out)
    original,source=read_sources(old_dir,cross_dir)
    originals={
        "parent_summary":digest(old_dir/"summary.json"),
        "parent_manifest":digest(old_dir/"manifest.json"),
        "cross_summary":digest(cross_dir/"summary.json"),
        "cross_manifest":digest(cross_dir/"manifest.json"),
        "cross_raw":digest(cross_dir/"crossings.csv.gz")}
    combined=aggregate_lut20(cross_dir/"crossings.csv.gz",original,source)
    original=dict(original)
    original["lut_trials_per_cell"]=500
    datasets=[{"label":"3.0 mm · LUT500 (original)", "data":original}]
    datasets.extend({"label":f"{p:.1f} mm · LUT20 (100k parents)",
                     "data":combined[p]} for p in PITCHES)
    payload={"schema":SCHEMA,"datasets":datasets,
             "note":"All curves are J/psi PARENT-event A. Separate per-muon crossing dashboards exist.",
             "source_sha256":originals}
    html=TEMPLATE.read_text(encoding="utf-8")
    if html.count("__PAYLOAD__")!=1 or "id=\"geometry\"" not in html:
        raise ValueError("Updated ROOT-style geometry template required")
    html=html.replace("__PAYLOAD__",
        json.dumps(payload,separators=(",",":"),ensure_ascii=False,allow_nan=False)
        .replace("</","<\\/"))
    if out.exists() and not force:
        if out.read_text(encoding="utf-8")==html:
            print("PASS existing enhanced dashboard matches, no changes:",out)
            return out
        raise FileExistsError(f"Refusing to overwrite existing comparison HTML: {out}")
    out.parent.mkdir(parents=True,exist_ok=True)
    stage=out.with_name(out.name+".tmp")
    if stage.exists():raise FileExistsError(stage)
    stage.write_text(html,encoding="utf-8")
    os.replace(stage,out)
    (out.with_suffix(".provenance.json")).write_text(json.dumps({
        "schema":SCHEMA,"dashboard_sha256":digest(out),
        "sources_sha256":originals,"unit":"one J/psi parent per station or DET",
        "geometry_labels":[item["label"] for item in datasets],
        "electron_count_definition":"total_up+total_down, classified primary-electron endpoints on physical anode wires; excludes boundaries"
    },indent=2)+"\n")
    print("PASS: "+str(out))
    print("Physical and original dashboard summaries and manifests unchanged.")
    return out

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--old",type=Path,default=OLD)
    p.add_argument("--full34",type=Path,default=CROSS)
    p.add_argument("--out",type=Path,default=OUTPUT)
    p.add_argument("--replace-html",action="store_true")
    a=p.parse_args()
    build(a.old,a.full34,a.out,a.replace_html)

if __name__=="__main__":main()
