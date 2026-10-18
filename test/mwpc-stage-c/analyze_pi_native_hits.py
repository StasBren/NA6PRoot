#!/usr/bin/env python3
"""Summarize actual Geant4 gas-hit Edep, incidence angle and path length.

Input is the small first-chunk ONLY native Geant4+MCKine-tagged input to
the anode event browser; not historical geometric extrapolations.
No ROOT or detector simulation is rerun.
"""
from __future__ import annotations
import argparse
import csv
import json
import math
from collections import Counter,defaultdict
from pathlib import Path

from mwpc_checkerboard_projection import SOURCES

REPO=Path(__file__).resolve().parents[2]
DEFAULT=REPO/"test_runs/mwpc_stage_c/pi_native_anode_first_chunk/crossings"


def mean(xs):
    return sum(xs)/len(xs) if xs else math.nan


def sample_quantile(xs,q):
    if not xs:
        return math.nan
    x=sorted(xs)
    v=q*(len(x)-1)
    lo=int(math.floor(v));hi=int(math.ceil(v))
    return x[lo]*(hi-v)+x[hi]*(v-lo) if hi!=lo else x[lo]


def analyze(root,out):
    out.mkdir(parents=True,exist_ok=True)
    group=defaultdict(lambda:{"energy":[],"length":[],"angleU":[],"angleW":[],
                              "tracks":set(),"corr":[]})
    modules=defaultdict(lambda:[0,0.,0.])
    station_tracks=defaultdict(set)
    nrows=0
    for src in SOURCES:
        path=root/f"crossings_{src}_ms2_530.csv"
        if not path.is_file():
            raise FileNotFoundError(path)
        with path.open(newline="",encoding="utf-8") as stream:
            for row in csv.DictReader(stream):
                if row.get("native_hit")!="1":
                    raise ValueError("Refusing historical projected crossings")
                station=int(row["station"]);det=int(row["det_id"])
                if not 0<=station<6 or not 0<=det<298:
                    raise ValueError("Invalid station/det")
                en=float(row["edep_eV"]);length=float(row["gas_path_mm"])
                if not all(map(math.isfinite,(en,length))) or en<0 or length<=0:
                    raise ValueError("Nonphysical gas hit")
                key=(src,station)
                label=(row["job"],row["event"],row["track"])
                d=group[key]
                d["energy"].append(en);d["length"].append(length)
                d["angleU"].append(abs(float(row["theta_u_deg"])))
                d["angleW"].append(abs(float(row["theta_w_deg"])))
                d["tracks"].add(label)
                d["corr"].append((length,en))
                station_tracks[(src,*label)].add(station)
                acc=modules[(src,station,det)]
                acc[0]+=1;acc[1]+=en;acc[2]+=length
                nrows+=1
    station_out=out/"pi_native_station_edep_angles.csv"
    fields=("source","station","gas_hit_rows","unique_muon_tracks",
            "mean_edep_eV","median_edep_eV","p90_edep_eV",
            "mean_path_mm","median_path_mm",
            "mean_abs_theta_u_deg","mean_abs_theta_w_deg",
            "corr_edep_vs_path")
    with station_out.open("w",newline="",encoding="utf-8") as stream:
        writer=csv.DictWriter(stream,fieldnames=fields);writer.writeheader()
        for key,d in sorted(group.items()):
            n=len(d["energy"]);ls=d["length"];es=d["energy"]
            el=mean([v*w for v,w in zip(es,ls)])
            de=mean([x*x for x in es])-mean(es)**2
            dl=mean([x*x for x in ls])-mean(ls)**2
            corr=((el-mean(es)*mean(ls))/math.sqrt(de*dl)
                  if de>0 and dl>0 else math.nan)
            writer.writerow(dict(source=key[0],station=key[1],
                gas_hit_rows=n,unique_muon_tracks=len(d["tracks"]),
                mean_edep_eV=mean(es),median_edep_eV=sample_quantile(es,.5),
                p90_edep_eV=sample_quantile(es,.9),
                mean_path_mm=mean(ls),median_path_mm=sample_quantile(ls,.5),
                mean_abs_theta_u_deg=mean(d["angleU"]),
                mean_abs_theta_w_deg=mean(d["angleW"]),
                corr_edep_vs_path=corr))
            print(f"{key[0]} MS{key[1]}: "
                  f"{n:,} actual gas hits / {len(d['tracks']):,} muons; "
                  f"Edep median {sample_quantile(es,.5):.1f} eV; "
                  f"L median {sample_quantile(ls,.5):.2f} mm",
                  flush=True)
    modfile=out/"pi_native_chamber_edep.csv"
    with modfile.open("w",newline="",encoding="utf-8") as stream:
        fields=("source","station","det_id","gas_hit_rows",
                "total_edep_eV","mean_edep_eV","mean_path_mm")
        writer=csv.DictWriter(stream,fieldnames=fields);writer.writeheader()
        for key,v in sorted(modules.items()):
            writer.writerow(dict(source=key[0],station=key[1],det_id=key[2],
                gas_hit_rows=v[0],total_edep_eV=v[1],
                mean_edep_eV=v[1]/v[0],mean_path_mm=v[2]/v[0]))
    complete={}
    for src in SOURCES:
        tracks=[stations for (source,*_),stations in station_tracks.items()
                if source==src]
        complete[src]=sum(len(x)==6 for x in tracks)
        print(f"{src}: {len(tracks):,} distinct muon track IDs with ≥1 gas hit; "
              f"{complete[src]} have hits in all 6 stations",flush=True)
    (out/"pi_native_first_chunk_summary.json").write_text(
        json.dumps({
            "input":"TRUE current-geometry Geant4 MuonSpec gas hits, first chunk/channel",
            "source_filter":"abs(PDG)=13 after valid remapped MCKine join",
            "total_muon_hit_rows":nrows,
            "unique_tracks_crossing_gas_all_six_stations":complete,
            "limitations":[
                "Per-chamber gas hits include overlaps, so they double-count some muons",
                "Unique track identity is (source,job,event,remapped track id)",
                "No inference of all station-plane incident counts from gas hits alone",
                "First chunk per source, not 200000/channel campaign statistics",
                "Heed primary electron realization and analog signal are separate from Geant4 Edep",
            ],
        },ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print("PASS: true Geant4 gas-hit Edep/path/angle first-chunk summaries:",out)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--crossings",type=Path,default=DEFAULT)
    ap.add_argument("--out",type=Path,default=DEFAULT.parent/"summaries")
    args=ap.parse_args()
    analyze(args.crossings,args.out)


if __name__=="__main__":
    main()
