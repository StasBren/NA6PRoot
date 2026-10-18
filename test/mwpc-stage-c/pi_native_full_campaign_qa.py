#!/usr/bin/env python3
"""Read-only PI quality analysis of ALL native Geant4 MWPC gas-hit chunks.

MCKine-remapped actual muon hits only. NO detector digitizer, Garfield,
cathodes or historical extrapolated gas intersections are mixed into this
report. Uses one temporary CSV at a time, independently validated chunk
checkpoints and streaming histogram aggregators (not full-run event arrays).

Usage after Stage-C library activation:
  python3 test/mwpc-stage-c/pi_native_full_campaign_qa.py

Finished per-chunk summaries are checkpointed so interruption is resumable.
This report does NOT equate muon-gas-crossing counts with incoming muons at
reference station planes. Count station incident muons separately before
calling a quantity "per incoming muon".
"""
from __future__ import annotations

import argparse
import bisect
from collections import Counter, defaultdict
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import time

from audit_pi_geant4_campaign import audit, DEFAULT as CAMPAIGN, CHANNELS
from mwpc_checkerboard_projection import build_chambers

HERE=Path(__file__).resolve().parent
REPO=HERE.parents[1]
EXPECTED={c.det_id:c.station for station in build_chambers() for c in station}
SOURCES={"Jpsi":"jpsi40-bore","Phi":"phi40-bore","Omega":"omega40-bore"}
EDEP=(0.,100.,200.,300.,400.,500.,750.,1000.,1250.,1500.,2000.,
      2500.,3000.,4000.,5000.,7500.,10000.,15000.,25000.,50000.,
      100000.,250000.,1.e20)
ANGLE=(-180.,-90.,-60.,-45.,-30.,-20.,-10.,-5.,-2.,-1.,0.,1.,
       2.,5.,10.,20.,30.,45.,60.,90.,180.)
DELTA=(-180.,-45.,-20.,-10.,-5.,-2.,-1.,-.25,-.05,0.,.05,.25,
       1.,2.,5.,10.,20.,45.,180.)
PATH=(0.,1.,2.,3.,4.,5.,5.5,6.,6.1,6.25,6.5,7.,8.,10.,
      15.,25.,50.,100.,1.e20)
VARIABLES={
    "edep_eV":EDEP, "path_mm":PATH,
    "theta_u_in_deg":ANGLE, "theta_u_out_deg":ANGLE,
    "theta_w_in_deg":ANGLE, "theta_w_out_deg":ANGLE,
    "delta_theta_u_deg":DELTA, "delta_theta_w_deg":DELTA,
}


def bin_of(value,edges):
    if not math.isfinite(value):
        return -1
    return max(0,min(len(edges)-2,bisect.bisect_right(edges,value)-1))


def fresh():
    return {
        "n":0,"edep_sum":0.,"edep_squared":0.,
        "path_sum":0.,"path_squared":0.,"path_edep":0.,
        "positive_edep":0,
        "h":{v:[0]*(len(edges)-1) for v,edges in VARIABLES.items()},
        "anomalies":Counter(),
    }


def update(acc,row):
    edep=float(row["edep_eV"])
    path=float(row["path_length_mm"])
    upIn=float(row["theta_u_in_deg"])
    upOut=float(row["theta_u_out_deg"])
    wpIn=float(row["theta_w_in_deg"])
    wpOut=float(row["theta_w_out_deg"])
    # Wrapped, signed angle differences are safer than raw +/- 180 jumps.
    du=((upOut-upIn+180.)%360.)-180.
    dw=((wpOut-wpIn+180.)%360.)-180.
    vals={
        "edep_eV":edep,"path_mm":path,
        "theta_u_in_deg":upIn,"theta_u_out_deg":upOut,
        "theta_w_in_deg":wpIn,"theta_w_out_deg":wpOut,
        "delta_theta_u_deg":du,"delta_theta_w_deg":dw,
    }
    bad=[k for k,v in vals.items() if not math.isfinite(v)]
    if bad:
        for v in bad:
            acc["anomalies"]["nonfinite_"+v]+=1
        return False
    if path<=0 or edep<0:
        acc["anomalies"]["nonphysical_edep_or_path"]+=1
        return False
    acc["n"]+=1
    acc["edep_sum"]+=edep
    acc["edep_squared"]+=edep*edep
    acc["path_sum"]+=path
    acc["path_squared"]+=path*path
    acc["path_edep"]+=path*edep
    acc["positive_edep"]+=edep>0
    if path>50:
        acc["anomalies"]["long_gas_segment_over_50mm"]+=1
    if abs(du)>5 or abs(dw)>5:
        acc["anomalies"]["delta_angle_over_5deg"]+=1
    for name,value in vals.items():
        index=bin_of(value,VARIABLES[name])
        if index>=0:
            acc["h"][name][index]+=1
    return True


def collapse(acc):
    out={k:v for k,v in acc.items() if k!="anomalies"}
    out["anomalies"]=dict(acc["anomalies"])
    return out


def absorb(dest,src):
    for key in ("n","edep_sum","edep_squared","path_sum","path_squared",
                "path_edep","positive_edep"):
        dest[key]+=src[key]
    for variable,values in src["h"].items():
        for index,value in enumerate(values):
            dest["h"][variable][index]+=value
    dest["anomalies"].update(src["anomalies"])


def process_csv(path):
    summary={"total_detector_hits":0, "muon_hits":0,"nonmuon_hits":0,
             "invalid_muon_hits":0, "unique_tracks_with_gas":0,
             "tracks_crossing_gas_in_all_six":0,
             "duplicate_hit_keys":0, "station":{}, "module":{},
             "station_unique_gas_tracks":{}, "station_overlap_tracks":{},
             "station_extra_chamber_hits":{}}
    stats_station=defaultdict(fresh)
    stats_module=defaultdict(fresh)
    seen_hits=set()
    tracks=defaultdict(int)
    per_station=defaultdict(set)
    with path.open(newline="",encoding="utf-8") as stream:
        rd=csv.DictReader(stream)
        required={"event","hit_index","track_id","det_id","station","pdg",
                  "edep_eV","path_length_mm",
                  "theta_u_in_deg","theta_u_out_deg",
                  "theta_w_in_deg","theta_w_out_deg"}
        missing=required-set(rd.fieldnames or [])
        if missing:
            raise ValueError(f"Native CSV missing columns {missing}")
        for row in rd:
            summary["total_detector_hits"]+=1
            ev=int(row["event"])
            tid=int(row["track_id"])
            st=int(row["station"])
            det=int(row["det_id"])
            pdg=int(row["pdg"])
            ident=(ev,int(row["hit_index"]))
            if ident in seen_hits:
                summary["duplicate_hit_keys"]+=1
            else:
                seen_hits.add(ident)
            if det not in EXPECTED or EXPECTED[det]!=st:
                raise ValueError(f"Wrong physical detector location: MS{st} DET{det}")
            if abs(pdg)!=13:
                summary["nonmuon_hits"]+=1
                continue
            summary["muon_hits"]+=1
            tracks[(ev,tid)]|=1<<st
            per_station[(st,ev,tid)].add(det)
            if not update(stats_station[st],row):
                summary["invalid_muon_hits"]+=1
            if not update(stats_module[det],row):
                raise ValueError("Station/module acceptance disagree")
    summary["unique_tracks_with_gas"]=len(tracks)
    summary["tracks_crossing_gas_in_all_six"]=sum(m==63 for m in tracks.values())
    for st in range(6):
        groups=[dets for (station,_,_),dets in per_station.items()
                if station==st]
        summary["station_unique_gas_tracks"][str(st)]=len(groups)
        summary["station_overlap_tracks"][str(st)]=sum(len(g)>1 for g in groups)
        summary["station_extra_chamber_hits"][str(st)]=(
            sum(len(g)-1 for g in groups)
        )
    summary["station"]={str(k):collapse(v)
                        for k,v in stats_station.items()}
    summary["module"]={str(k):collapse(v)
                       for k,v in stats_module.items()}
    if summary["total_detector_hits"]!=(summary["muon_hits"]+
                                         summary["nonmuon_hits"]):
        raise AssertionError("Detector hit count not conserved")
    if sum(v["n"] for v in stats_station.values())!=(
            summary["muon_hits"]-summary["invalid_muon_hits"]):
        raise AssertionError("Station accepted hit count not conserved")
    if sum(v["n"] for v in stats_module.values())!=(
            summary["muon_hits"]-summary["invalid_muon_hits"]):
        raise AssertionError("Module accepted hit count not conserved")
    for st in range(6):
        name=str(st)
        unique=summary["station_unique_gas_tracks"][name]
        hits=sum(1 for (station,_,_) in per_station if station==st)
        assert hits==unique
        if stats_station[st]["n"]<unique:
            raise ValueError("Fewer valid hits than unique station muon tracks")
    if summary["duplicate_hit_keys"]:
        raise ValueError("Non-unique (event,hit_index) in native ROOT export")
    return summary


def save_atomic(path,data):
    tmp=path.with_suffix(path.suffix+".tmp")
    with tmp.open("w",encoding="utf-8") as stream:
        json.dump(data,stream,ensure_ascii=False,indent=2,allow_nan=False)
        stream.write("\n")
    os.replace(tmp,path)


def file_id(path):
    stat=path.stat()
    return {"name":path.name,"bytes":stat.st_size,"mtime_ns":stat.st_mtime_ns}


def export_chunk(root,chunk,out,log):
    mac=HERE/"export_c_geant_muon_crossings.C"
    cmd=f'{mac}+("{chunk}/HitsMuonSpec.root","{chunk}/geometry.root",'+\
        f'"{out}","{chunk}/MCKine.root")'
    with log.open("w",encoding="utf-8") as stream:
        p=subprocess.run([root,"-l","-b","-q",cmd],stdout=stream,
                         stderr=subprocess.STDOUT,check=False,cwd=REPO)
    if p.returncode or not out.is_file() or not out.stat().st_size or (
       "GEANT MWPC crossing export COMPLETE" not in log.read_text(errors="replace")
    ):
        raise RuntimeError(f"Native ROOT export failed: {log}")
    return mac


def all_chunks(args):
    status=audit(args.campaign,False)
    if status["status"]!="COMPLETE":
        raise RuntimeError("Full Geant4 campaign not complete")
    frozen=json.loads((args.campaign/"campaign.json").read_text())["frozen_inputs"]
    chunks=math.ceil(int(frozen["events_per_source"])/int(frozen["chunk_size"]))
    selected=chunks if not args.first_chunks else min(chunks,args.first_chunks)
    args.out.mkdir(parents=True,exist_ok=True)
    checkpoints=args.out/"checkpoints"
    checkpoints.mkdir(exist_ok=True)
    temp=args.out/"_export_in_progress.csv"
    logfile=args.out/"_export_in_progress.log"
    codepath=HERE/"export_c_geant_muon_crossings.C"
    # Epoch freezes content-derived exporter revision. Results generated
    # using different macro versions CANNOT share one output directory.
    epoch=hashlib.sha256(codepath.read_bytes()).hexdigest()
    root=os.environ.get("PI_ROOT_EXECUTABLE","root")
    combined=defaultdict(fresh)
    module=defaultdict(fresh)
    global_metrics=defaultdict(Counter)
    print("PI QA VERSION:",epoch[:12],flush=True)
    for channel in CHANNELS:
        source=SOURCES[channel]
        for i in range(selected):
            chunk=args.campaign/channel.lower()/f"chunk{i+1:04d}"
            target=checkpoints/f"{channel.lower()}_chunk{i+1:04d}.json"
            inputs={name:file_id(chunk/name)
                    for name in ("HitsMuonSpec.root","MCKine.root","geometry.root")}
            if target.is_file():
                saved=json.loads(target.read_text(encoding="utf-8"))
                if saved.get("exporter_sha256")!=epoch or saved.get("input")!=inputs:
                    raise ValueError(f"Checkpoint source changed: {target}; "
                                     "use a NEW --out directory")
                this=saved["result"]
                label="SKIP"
            else:
                if temp.exists():
                    temp.unlink() # never interpret a prior crash temp as a checkpoint
                print(f"RUN {channel} chunk{i+1:04d}/{selected}: ROOT export...",flush=True)
                export_chunk(root,chunk,temp,logfile)
                this=process_csv(temp)
                # One step-level hit may have been stopped at the chamber
                # boundary and have zero geometric path. Such anomalies
                # are visible in invalid_muon_hits and not silently replaced.
                save_atomic(target,{
                    "exporter_sha256":epoch,
                    "input":inputs,"result":this,
                    "channel":channel,"chunk":i+1
                })
                temp.unlink()
                label="PASS"
            for k,v in this["station"].items():
                absorb(combined[(source,int(k))],v)
            for k,v in this["module"].items():
                absorb(module[(source,int(k))],v)
            for k in ("total_detector_hits","muon_hits","nonmuon_hits",
                      "invalid_muon_hits","unique_tracks_with_gas",
                      "tracks_crossing_gas_in_all_six"):
                global_metrics[source][k]+=int(this[k])
            for st in range(6):
                key=(source,st)
                for metric in ("station_unique_gas_tracks","station_overlap_tracks",
                               "station_extra_chamber_hits"):
                    global_metrics[key][metric]+=int(this[metric][str(st)])
            print(f"{label} {channel} chunk{i+1:04d}: "
                  f"muon gas hits {this['muon_hits']:,}, "
                  f"nonmuon {this['nonmuon_hits']:,}",flush=True)
    if temp.exists():
        temp.unlink()
    def corr(r):
        n=r["n"]
        if n<=1:
            return None
        dx=r["path_squared"]/n-(r["path_sum"]/n)**2
        dy=r["edep_squared"]/n-(r["edep_sum"]/n)**2
        if dx<=0 or dy<=0:
            return None
        return ((r["path_edep"]/n-(r["path_sum"]/n)*(r["edep_sum"]/n))/
                math.sqrt(dx*dy))

    with (args.out/"pi_native_all_chamber_quality.csv").open(
       "w",newline="",encoding="utf-8") as stream:
        fields=("source","station","det_id","gas_hit_rows",
                "mean_edep_eV","mean_path_mm","positive_edep_fraction",
                "edep_path_pearson","nonphysical_count",
                "long_path_over_50mm","large_direction_change")
        w=csv.DictWriter(stream,fieldnames=fields);w.writeheader()
        for source in SOURCES.values():
            for det,st in sorted(EXPECTED.items()):
                r=module[(source,det)]
                n=r["n"]
                w.writerow({
                    "source":source,"station":st,"det_id":det,
                    "gas_hit_rows":n,
                    "mean_edep_eV":r["edep_sum"]/n if n else "",
                    "mean_path_mm":r["path_sum"]/n if n else "",
                    "positive_edep_fraction":r["positive_edep"]/n if n else "",
                    "edep_path_pearson":corr(r) if corr(r) is not None else "",
                    "nonphysical_count":r["anomalies"].get(
                        "nonphysical_edep_or_path",0),
                    "long_path_over_50mm":r["anomalies"].get(
                        "long_gas_segment_over_50mm",0),
                    "large_direction_change":r["anomalies"].get(
                        "delta_angle_over_5deg",0),
                })
    with (args.out/"pi_native_all_station_quality.csv").open(
      "w",newline="",encoding="utf-8") as stream:
        fields=("source","station","gas_hit_rows","unique_gas_crossing_muon_tracks",
                "overlap_muon_tracks","extra_chamber_hits_per_overlaps",
                "mean_edep_eV","mean_path_mm","edep_path_pearson")
        w=csv.DictWriter(stream,fieldnames=fields);w.writeheader()
        for source in SOURCES.values():
            for st in range(6):
                r=combined[(source,st)]
                m=global_metrics[(source,st)]
                n=r["n"]
                w.writerow({
                    "source":source,"station":st,"gas_hit_rows":n,
                    "unique_gas_crossing_muon_tracks":
                       m["station_unique_gas_tracks"],
                    "overlap_muon_tracks":m["station_overlap_tracks"],
                    "extra_chamber_hits_per_overlaps":
                       m["station_extra_chamber_hits"],
                    "mean_edep_eV":r["edep_sum"]/n if n else "",
                    "mean_path_mm":r["path_sum"]/n if n else "",
                    "edep_path_pearson":corr(r) if corr(r) is not None else "",
                })
    with (args.out/"pi_native_all_histograms.csv").open(
      "w",newline="",encoding="utf-8") as stream:
        fields=("source","station","det_id","variable","bin_lower","bin_upper",
                "hit_count","gas_hit_exposure","conditional_bin_fraction")
        w=csv.DictWriter(stream,fieldnames=fields);w.writeheader()
        for source in SOURCES.values():
            for det,st in sorted(EXPECTED.items()):
                r=module[(source,det)]
                for variable,edges in VARIABLES.items():
                    for b,count in enumerate(r["h"][variable]):
                        w.writerow({
                            "source":source,"station":st,"det_id":det,
                            "variable":variable,"bin_lower":edges[b],
                            "bin_upper":edges[b+1],"hit_count":count,
                            "gas_hit_exposure":r["n"],
                            "conditional_bin_fraction":
                                 count/r["n"] if r["n"] else "",
                        })
    report={
        "analysis":"True Geant4-MCKine joined gas-hit quality control, "
                   "no simulated anode firing",
        "source_campaign":str(args.campaign),
        "exporter_sha256":epoch,
        "processed_chunks_per_source":selected,
        "total_chunks_per_source":chunks,
        "full_campaign_analyzed":selected==chunks,
        "generated_events_per_source":(
            int(frozen["events_per_source"]) if selected==chunks
            else min(int(frozen["events_per_source"]),
                     selected*int(frozen["chunk_size"]))
        ),
        "per_source":{src:dict(global_metrics[src])
                      for src in SOURCES.values()},
        "limitations":[
            "Reference-plane incident muon denominator has not yet been "
            "joined into this report. Gas-hit tracks are conditional on "
            "touching active gas, not incoming station reference muons.",
            "First 298 chambers use current fixed Stage-C layout; no "
            "alternative-layout Geant4 Edep can be inferred from this.",
            "A Geant4 gas hit may occur multiple times for one track "
            "in a chamber. Duplicate chamber-hit segments require "
            "per-event deduplication before defining 0/1/2/3 wire PDF.",
            "No electrons, drift times, avalanche gain, anode thresholds "
            "or time-dependent wire multiplicity are measured here.",
            "Histograms include every valid native muon-hit segment; "
            "variable bins with zero entries are preserved.",
        ]
    }
    save_atomic(args.out/"pi_native_all_QA_report.json",report)
    print("PASS: PI FULL native gas Edep, in/out U/W angle, path, "
          "delta-angle, station/chamber overlap QA")
    print("Wrote:",args.out)
    if not report["full_campaign_analyzed"]:
        print("NOTE: subset only, full campaign NOT yet summarized",flush=True)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--campaign",type=Path,default=CAMPAIGN)
    ap.add_argument("--out",type=Path,default=REPO/
         "test_runs/mwpc_stage_c/pi_native_full_quality")
    ap.add_argument("--first-chunks",type=int,default=0,
         help="Smoke (e.g. 1) OR 0 for all 67 chunks/source")
    args=ap.parse_args()
    if args.first_chunks<0:
        ap.error("--first-chunks cannot be negative")
    try:
        all_chunks(args)
    except Exception as exc:
        print("PI_NATIVE_QA_ERROR:",exc,file=sys.stderr,flush=True)
        raise SystemExit(1)


if __name__=="__main__":
    main()
