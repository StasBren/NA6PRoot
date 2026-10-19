#!/usr/bin/env python3
"""Anode-only primary-electron DRIFT TIME PDFs on TRUE native Geant4 gas hits.

All first 3000 generated events per Jpsi/Phi/Omega (not cherry-picked
Muon Track Explorer events). Native Geant4 gas segment geometry and
Heed+already-saved Garfield electron timing trials. NO new Geant4,
Garfield microscopic simulation, cathode strips or analog thresholds.

This is a PILOT for the PI timing question. Results are CONDITIONAL on
a native muon gas-hit SEGMENT. This is NOT the completed 200k/channel
timing-fold, not incident-station-muon normalized and not hardware
threshold crossing time.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import sys

import numpy as np

from build_pi_muon_track_explorer import stable_u64
from mwpc_checkerboard_projection import SOURCES, build_chambers
from pi_fixed_wire_occupancy import (
    cluster_native_gas_segment, collected_electron_endpoints, n_wires,
)
from run_pi_per_chamber_response import read_heed, read_garfield

HERE=Path(__file__).resolve().parent
REPO=HERE.parents[1]
DEFAULT_NATIVE=REPO/"test_runs/mwpc_stage_c/pi_native_anode_first_chunk/crossings"
DEFAULT_RESULTS=REPO/"test_runs/mwpc_stage_c/pi_native_drift_time_first_chunk"
DEFAULT_HTML=REPO/"test_runs/mwpc_stage_c/pi_fixed_wire_maps/pi_native_drift_times.html"
DEFAULT_HEED=REPO/"test_runs/mwpc_stage_c/c_heed_library_2p4"
GEOMETRY={c.det_id:c for st in build_chambers() for c in st}
# 100 bins of width 20 ns, plus one unbounded overflow bin.
EDGES=np.array(list(range(0,2001,20))+[np.inf],dtype=float)
TIMED_WINDOWS=(50,100,200,400,800)
WINDOWS=tuple(map(str,TIMED_WINDOWS))+("all",)
REQUIRED={"source","native_hit","job","event","track","pdg",
          "station","det_id","hit_index","gas_path_mm",
          "u_in_mm","u_out_mm","v_in_mm","v_out_mm"}


def blank():
    return {
        "n":0,"collected":0,"validTimes":0,"unknownTimes":0,
        "edgeLoss":0,"zeroCollectedHits":0,
        "hist":np.zeros(len(EDGES)-1,dtype=np.int64),
        "multiplicity":np.zeros((len(WINDOWS),6),dtype=np.int64),
    }


def add_hit(bucket,wires,times,lost):
    """One independently modeled native gas-hit segment (including N=0)."""
    wires=np.asarray(wires,dtype=np.int64)
    times=np.asarray(times,dtype=float)
    if wires.shape!=times.shape:
        raise ValueError("Electron wire/time vectors differ in length")
    if lost<0:
        raise ValueError("Negative finite-anode edge electron loss")
    bucket["n"]+=1
    bucket["collected"]+=wires.size
    bucket["edgeLoss"]+=lost
    bucket["zeroCollectedHits"]+=wires.size==0
    valid=np.isfinite(times)&(times>=0)
    good=times[valid]
    bucket["validTimes"]+=good.size
    bucket["unknownTimes"]+=wires.size-good.size
    bucket["hist"]+=np.histogram(good,bins=EDGES)[0]
    for iw,window in enumerate(TIMED_WINDOWS):
        multiplicity=np.unique(wires[valid & (times<=window)]).size
        bucket["multiplicity"][iw,min(multiplicity,5)]+=1
    all_multiplicity=np.unique(wires).size
    bucket["multiplicity"][-1,min(all_multiplicity,5)]+=1
    if not np.all(bucket["multiplicity"].sum(axis=1)==bucket["n"]):
        raise AssertionError("Multiplicity denominators must count every gas hit")
    if int(bucket["hist"].sum())!=bucket["validTimes"]:
        raise AssertionError("Time histogram must account for all valid drift times")


def encode(bucket):
    return {
        "n":bucket["n"],"collected":int(bucket["collected"]),
        "validTimes":int(bucket["validTimes"]),
        "unknownTimes":int(bucket["unknownTimes"]),
        "edgeLoss":int(bucket["edgeLoss"]),
        "zeroCollectedHits":int(bucket["zeroCollectedHits"]),
        "hist":bucket["hist"].tolist(),
        "multiplicity":bucket["multiplicity"].tolist(),
    }


def build(args):
    for p in (args.native,args.heed_tracks,args.heed_clusters,
              args.garfield,args.template):
        if not p.exists():
            raise FileNotFoundError(
                f"Required offline file missing: {p}. "
                "Build first-chunk native inputs with "
                "run_pi_native_anode_first_chunk.sh if needed.")
    lam, marks=read_heed(args.heed_tracks,args.heed_clusters)
    transport=read_garfield(args.garfield)
    stations={source:[blank() for _ in range(6)] for source in SOURCES}
    modules={source:{} for source in SOURCES}
    readout=[]
    for source in SOURCES:
        path=args.native/f"crossings_{source}_ms2_530.csv"
        if not path.is_file():
            raise FileNotFoundError(path)
        n=0
        with path.open(newline="",encoding="utf-8") as stream:
            reader=csv.DictReader(stream)
            if not REQUIRED.issubset(reader.fieldnames or []):
                raise ValueError(f"{path}: missing required fields")
            for row in reader:
                if row["native_hit"]!="1" or row["source"]!=source or abs(int(row["pdg"]))!=13:
                    raise ValueError("Timing PDF accepts ONLY native MCKine-joined muon gas hits")
                st=int(row["station"]);det=int(row["det_id"])
                if det not in GEOMETRY or GEOMETRY[det].station!=st:
                    raise ValueError(f"Native station/physical chamber mismatch {st}:{det}")
                seed=stable_u64("PI_NATIVE_TIMING_FIRST_CHUNK_V1",
                                source,row["job"],row["event"],row["track"],
                                st,det,row["hit_index"])
                rng=np.random.default_rng(seed)
                u,v,_,_=cluster_native_gas_segment(row,rng,lam,marks)
                wires,times,lost=collected_electron_endpoints(
                    u,v,st,0.,transport,rng.random(u.size)
                )
                add_hit(stations[source][st],wires,times,lost)
                mod=modules[source].setdefault(str(det),blank())
                add_hit(mod,wires,times,lost)
                n+=1
                if n%5000==0:
                    print(f"{source}: {n:,} native muon gas hits...",flush=True)
        if n<=0:
            raise ValueError(f"No native gas-hit rows for {source}")
        if sum(x["n"] for x in stations[source])!=n or sum(
            x["n"] for x in modules[source].values()
        )!=n:
            raise AssertionError("Station/module/source hit-count disagreement")
        for st in range(6):
            nmods=sum(v["n"] for det,v in modules[source].items()
                      if GEOMETRY[int(det)].station==st)
            if nmods!=stations[source][st]["n"]:
                raise AssertionError(f"Station{st} != sum of physical modules")
        readout.append((source,n,sum(x["validTimes"] for x in stations[source])))
        print(f"{source}: PASS {n:,} unbiased first-chunk muon gas segments / "
              f"{readout[-1][2]:,} drift-timed collected electrons",
              flush=True)
    payload={
        "sources":list(SOURCES),
        "edges":EDGES.tolist()[:-1]+["overflow"],
        "windows":list(WINDOWS),
        "station":{src:[encode(x) for x in stations[src]] for src in SOURCES},
        "modules":{src:{det:encode(x) for det,x in modules[src].items()} for src in SOURCES},
        "geometry":{str(det):c.station for det,c in GEOMETRY.items()},
        "meta":{
            "scope":"ALL NATIVE MUON GAS HIT SEGMENTS, FIRST 3000 GENERATED EVENTS PER SOURCE ONLY, not full 200000",
            "notFullRun":True,
            "usesActualGeantGasGeometry":True,
            "model":"LEGACY ASYMMETRIC 2+4mm Garfield electric field, NOT centered 3+3mm: Heed sampled cluster ionization, B=0 saved Garfield single-electron transit samples, fixed 77/128 anodes and +1mm 6mm-gap approximation",
            "meaning":"Electron drift time from ionization to anode, measured in ns; NOT Geant4 particle flight time or electronics discriminator time",
            "threshold":"≥1 collected primary electron per anode. No avalanche, analog shaping, noise, or ASIC threshold",
            "exposure":"Per REAL Geant4 gas HIT SEGMENT, not incoming station muon",
            "unknownTimes":"Collected Garfield outcomes without valid nonnegative finite drift times remain in all-collected metric but excluded from timed windows",
        }
    }
    args.out.mkdir(parents=True,exist_ok=True)
    with (args.out/"pi_native_drift_time_histograms.csv").open(
        "w",newline="",encoding="utf-8"
    ) as stream:
        fields=("source","station","det_id","t_min_ns","t_max_ns",
                "collected_electron_count","valid_collected_electrons","muon_gas_hit_rows")
        writer=csv.DictWriter(stream,fieldnames=fields);writer.writeheader()
        for src in SOURCES:
            for det in range(298):
                st=GEOMETRY[det].station
                b=modules[src].get(str(det),blank())
                for i,n in enumerate(b["hist"]):
                    writer.writerow(dict(
                        source=src,station=st,det_id=det,
                        t_min_ns=float(EDGES[i]),
                        t_max_ns=float(EDGES[i+1]) if np.isfinite(EDGES[i+1]) else "inf",
                        collected_electron_count=int(n),
                        valid_collected_electrons=b["validTimes"],
                        muon_gas_hit_rows=b["n"],
                    ))
    with (args.out/"pi_native_time_window_multiplicity.csv").open(
        "w",newline="",encoding="utf-8"
    ) as stream:
        fields=("source","station","det_id","cutoff_ns","muon_gas_hit_rows",
                "N0","N1","N2","N3","N4","N5plus","P_ge2",
                "collected_primary_electrons","valid_drift_times",
                "unknown_drift_times","edge_lost_electrons")
        wr=csv.DictWriter(stream,fieldnames=fields);wr.writeheader()
        for src in SOURCES:
            for det in range(298):
                st=GEOMETRY[det].station
                b=modules[src].get(str(det),blank())
                for k,w in enumerate(WINDOWS):
                    counts=b["multiplicity"][k]
                    wr.writerow(dict(
                        source=src,station=st,det_id=det,
                        cutoff_ns=w,
                        muon_gas_hit_rows=b["n"],
                        **{key:int(v) for key,v in zip(
                            ("N0","N1","N2","N3","N4","N5plus"),counts)},
                        P_ge2=float(counts[2:].sum())/b["n"] if b["n"] else "",
                        collected_primary_electrons=b["collected"],
                        valid_drift_times=b["validTimes"],
                        unknown_drift_times=b["unknownTimes"],
                        edge_lost_electrons=b["edgeLoss"],
                    ))
    summary={"model":payload["meta"],"samples":readout,
             "garfield_sha256":hashlib.sha256(args.garfield.read_bytes()).hexdigest(),
             "heed_library_tracks_sha256":hashlib.sha256(args.heed_tracks.read_bytes()).hexdigest(),
             "heed_library_clusters_sha256":hashlib.sha256(args.heed_clusters.read_bytes()).hexdigest()}
    (args.out/"pi_native_drift_time_provenance.json").write_text(
        json.dumps(summary,indent=2,ensure_ascii=False)+"\n",encoding="utf-8"
    )
    marker="__TIMING_PAYLOAD__"
    html=args.template.read_text(encoding="utf-8")
    if html.count(marker)!=1:
        raise ValueError("Expected exactly one HTML payload marker")
    args.html.parent.mkdir(parents=True,exist_ok=True)
    args.html.write_text(html.replace(marker,json.dumps(
        payload,ensure_ascii=False,separators=(",",":")).replace("</","<\\/")
    ),encoding="utf-8")
    print("PASS: 894 modules, six time windows, unselected full first-chunk"
          " native muon hit sample; no Geant4/Garfield rerun")
    print("Interactive drift-time graph:",args.html)
    print("Station and module raw timing tables:",args.out)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--native",type=Path,default=DEFAULT_NATIVE)
    ap.add_argument("--out",type=Path,default=DEFAULT_RESULTS)
    ap.add_argument("--html",type=Path,default=DEFAULT_HTML)
    ap.add_argument("--template",type=Path,default=HERE/"pi_native_drift_times_template.html")
    ap.add_argument("--heed-tracks",type=Path,default=DEFAULT_HEED/"heed_2p4_tracks.csv")
    ap.add_argument("--heed-clusters",type=Path,default=DEFAULT_HEED/"heed_2p4_clusters.csv")
    ap.add_argument("--garfield",type=Path,default=None)
    args=ap.parse_args()
    if args.garfield is None:
        a=REPO/"test_runs/mwpc_stage_c/c_transport_box_map_hybrid_boundary_2p4/c_transport_box_electrons.csv"
        b=REPO/"test_runs/mwpc_stage_c/c_transport_box_map_gm2p0_gp4p0/c_transport_box_electrons.csv"
        args.garfield=a if a.is_file() else b
    try:
        build(args)
    except Exception as exc:
        print("PI_NATIVE_DRIFT_TIMES_ERROR:",exc,file=sys.stderr,flush=True)
        raise SystemExit(1) from exc


if __name__=="__main__":
    main()
