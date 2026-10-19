#!/usr/bin/env python3
"""Centered 3+3mm MWPC anode timing-angle PILOT on TRUE native Geant4 gas hits.

Heed primary ionization along actual 3D TGeo-local (u,v,w) gas segments.
Garfield sample drawn ONLY from the newly computed and SHA-validated 3+3mm
response map. Data are grouped by signed native muon entrance angle.

Outputs: 1D electron arrival-time histograms and birth-position/time density
for u, v, w plus angle-specific P(N_anodes=0..5+) at time gates.

NOT a drift trajectory position u(t),v(t),w(t): snapshots require additional
microscopic step recording. NO Geant4 rerun, cathode strips or electronics.
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
from pi_fixed_wire_occupancy import nominal_wire_index, n_wires
from pi_validate_centered_box_map import validate_map
from run_pi_per_chamber_response import read_heed, read_garfield

HERE=Path(__file__).resolve().parent
REPO=HERE.parents[1]
NATIVE=REPO/"test_runs/mwpc_stage_c/pi_native_anode_first_chunk/crossings"
TRANSPORT=REPO/"test_runs/mwpc_stage_c/c_transport_box_map_gm3p0_gp3p0"
OUT=REPO/"test_runs/mwpc_stage_c/pi_centered_timing_native_first_chunk"
ANGLE_EDGES=(-180.,-30.,-10.,-5.,-2.,0.,2.,5.,10.,30.,180.)
TIMES=(25,50,100,200,400,800)
GATES=tuple(map(str,TIMES))+("all",)
# 25 x 10ns plus overflow [250,inf): all measured values remain counted.
T_EDGES=np.array(list(range(0,251,10))+[np.inf],dtype=float)
COORD_EDGES={
    "u":np.linspace(-2,2,21), # u relative to nearest wire, pitch 4mm
    "v":np.linspace(-3,3,31), # normal, centered 3+3 mm field
    "w":np.linspace(-340,340,21), # along wire, TGeo-local; 66.56cm gas length
}
GEOMETRY={c.det_id:c for group in build_chambers() for c in group}
NEEDED={"source","native_hit","job","event","track","hit_index",
        "pdg","station","det_id","gas_path_mm","u_in_mm","u_out_mm",
        "v_in_mm","v_out_mm","w_in_mm","w_out_mm",
        "theta_u_in_deg","theta_w_in_deg"}


def blank():
    return {
        "gas_hits":0,"clusters":0,"primary_electrons":0,
        "collected_electrons":0,"valid_drift_times":0,
        "unknown_drift_times":0,"edge_lost_electrons":0,
        "t":np.zeros(len(T_EDGES)-1,dtype=np.int64),
        "coord":{k:np.zeros((len(e)-1,len(T_EDGES)-1),dtype=np.int64)
                 for k,e in COORD_EDGES.items()},
        "multi":np.zeros((len(GATES),6),dtype=np.int64),
    }


def angle_index(theta):
    if not math.isfinite(theta) or not -180 <= theta <= 180:
        raise ValueError(f"Invalid signed TGeo entrance angle {theta}")
    return max(0,min(len(ANGLE_EDGES)-2,
                     int(np.searchsorted(ANGLE_EDGES,theta,side="right")-1)))


def born_on_segment(row,rng,lam,marks):
    """Return 3D true local birth points and clusters, WITHOUT +1mm v shift."""
    path=float(row["gas_path_mm"])
    if not(math.isfinite(path) and 0<path<1000):
        raise ValueError(f"Invalid gas path: {path}")
    ncl=int(rng.poisson(lam*path))
    if ncl==0:
        return (np.empty(0),np.empty(0),np.empty(0),0)
    n_per_cluster=marks[rng.integers(len(marks),size=ncl)]
    loc=np.repeat(rng.random(ncl),n_per_cluster)
    coords=[]
    for name in ("u","v","w"):
        a=float(row[name+"_in_mm"])
        b=float(row[name+"_out_mm"])
        coords.append(a+loc*(b-a))
    u,v,w=coords
    if len(v) and (v.min() < -3.05 or v.max() > 3.05):
        raise ValueError(f"Native symmetric 3+3mm gas v outside [-3,+3]: "
                         f"{v.min()}..{v.max()}; bad TGeo frame?")
    return u,v,w,ncl


def centered_endpoints(u,v,station,transport,rng):
    """Response-map electron endpoints for v in -3..+3, NEVER legacy +1mm."""
    if not len(u):
        return (np.empty(0,dtype=np.int64),np.empty(0),np.empty(0),0,
                np.empty(0,dtype=bool))
    offsets,sizes,responses=transport
    n=n_wires(station)
    closest=nominal_wire_index(u,station,0.)
    urel=u-(closest-(n-1.)/2.)*4.
    if np.any(urel < -2.00001) or np.any(urel > 2.00001):
        raise ValueError("Nearest-wire u phase inconsistent")
    if np.any(v < -3.05) or np.any(v > 3.05):
        raise ValueError("Incorrect symmetric Garfield normal v coordinate")
    iu=np.clip(np.floor((urel+2.)/.2).astype(np.int64),0,19)
    iv=np.clip(np.floor((v+3.)/.2).astype(np.int64),0,29)
    box=30*iu+iv
    if np.any(sizes[box]<=0):
        raise ValueError("Incomplete centered Garfield response map")
    trial=offsets[box]+np.minimum((rng.random(len(u))*sizes[box]).astype(
        np.int64),sizes[box]-1)
    choice=responses[trial]
    collected=choice[:,0]>.5
    wire=closest+choice[:,1].astype(np.int64)
    good=collected&(wire>=0)&(wire<n)
    lost=int(np.count_nonzero(collected&~good))
    return wire[good],choice[good,2],urel[good],lost,good


def consume(a,wire,times,urel,v,w,ncl,nprimary,lost):
    """Single native gas segment, including ionization-free (N0) segments."""
    a["gas_hits"]+=1
    a["clusters"]+=ncl
    a["primary_electrons"]+=nprimary
    a["collected_electrons"]+=len(wire)
    a["edge_lost_electrons"]+=lost
    valid=np.isfinite(times)&(times>=0)
    a["valid_drift_times"]+=int(valid.sum())
    a["unknown_drift_times"]+=int((~valid).sum())
    tv=times[valid]
    a["t"]+=np.histogram(tv,bins=T_EDGES)[0]
    for coord,vals in (("u",urel),("v",v),("w",w)):
        # Exactly one born-position/time cell per known-t collected electron.
        hist,_,_=np.histogram2d(vals[valid],tv,bins=(
            COORD_EDGES[coord],T_EDGES))
        a["coord"][coord]+=hist.astype(np.int64)
    for i,t in enumerate(TIMES):
        multiplicity=np.unique(wire[valid & (times<=t)]).size
        a["multi"][i,min(multiplicity,5)]+=1
    a["multi"][-1,min(np.unique(wire).size,5)]+=1
    if np.any(a["multi"].sum(axis=1)!=a["gas_hits"]):
        raise AssertionError("Every timing gate needs the same gas-hit denominator")


def serialize(a):
    if int(a["t"].sum())!=a["valid_drift_times"]:
        raise AssertionError("Unaccounted known electron arrival times")
    for variable,h in a["coord"].items():
        if int(h.sum())!=a["valid_drift_times"]:
            raise AssertionError(f"Out-of-range electron birth {variable} coordinate")
    return dict((k,int(v)) for k,v in a.items() if isinstance(v,(int,np.integer))) | {
        "arrival_time_counts":a["t"].tolist(),
        "birth_vs_time":{k:v.tolist() for k,v in a["coord"].items()},
        "multiplicity_counts":a["multi"].tolist(),
    }


def build(args):
    with args.manifest.open() as stream:
        manifest=json.load(stream)
    if manifest.get("geometry_id")!="MWPC_GARFIELD_CENTERED_3p0_3p0_B0_V1":
        raise ValueError("Expected validated NEW symmetric 3+3mm Garfield map")
    if args.garfield.resolve()!=args.manifest.parent.joinpath(
            "c_transport_box_electrons.csv").resolve():
        raise ValueError("Map and manifest must belong to the same geometry folder")
    sha=hashlib.sha256(args.garfield.read_bytes()).hexdigest()
    if sha!=manifest.get("response_csv_sha256"):
        raise ValueError("Garfield centered response CSV differs from validated manifest")
    check=validate_map(args.garfield,int(manifest["electrons_per_box"]))
    if check["response_csv_sha256"]!=sha:
        raise ValueError("Centered map validation changed or failed")
    lam,marks=read_heed(args.heed_tracks,args.heed_clusters)
    transport=read_garfield(args.garfield)
    aggregated={}
    counts={}
    unique=[]
    for source in SOURCES:
        nvalid=0
        nrows=0
        skipped=0
        fname=args.native/f"crossings_{source}_ms2_530.csv"
        seen=set()
        with fname.open(newline="",encoding="utf-8") as file:
            reader=csv.DictReader(file)
            missing=NEEDED-set(reader.fieldnames or [])
            if missing:
                raise ValueError(f"{fname}: missing {sorted(missing)}")
            for row in reader:
                if (row["native_hit"]!="1" or row["source"]!=source or
                    abs(int(row["pdg"]))!=13):
                    raise ValueError("Expected Geant4 true muon hits joined to MCKine")
                key=(row["event"],row["hit_index"])
                if key in seen:
                    raise ValueError(f"Duplicate real gas hit {source} {key}")
                seen.add(key)
                if args.max_per_source and nrows>=args.max_per_source:
                    break
                nrows+=1
                st=int(row["station"]);det=int(row["det_id"])
                if det not in GEOMETRY or GEOMETRY[det].station!=st:
                    raise ValueError(f"Unexpected physical chamber {source} {st=} {det=}")
                if float(row["gas_path_mm"])<=0:
                    skipped+=1
                    continue
                rng=np.random.default_rng(stable_u64(
                    "NATIVE_CENTERED_PI_BIRTH_TIME_V1",source,"chunk0001",
                    row["event"],row["track"],row["hit_index"],st,det))
                u,v,w,ncl=born_on_segment(row,rng,lam,marks)
                nprimary=len(u)
                if nprimary:
                    wire,times,urel,lost,good=centered_endpoints(
                        u,v,st,transport,rng)
                    vgood,wgood=v[good],w[good]
                else:
                    wire=np.empty(0,dtype=np.int64)
                    times=np.empty(0)
                    urel=np.empty(0)
                    vgood=wgood=np.empty(0)
                    lost=0
                nvalid+=1
                for axis,theta_key in (("theta_u","theta_u_in_deg"),
                                       ("theta_w","theta_w_in_deg")):
                    abin=angle_index(float(row[theta_key]))
                    for chosen_bin in (-1,abin):
                        key=(source,st,axis,chosen_bin)
                        bucket=aggregated.setdefault(key,blank())
                        consume(bucket,wire,times,urel,vgood,wgood,
                                ncl,nprimary,lost)
                if nvalid%5000==0:
                    print(f"{source}: {nvalid:,} centered true gas-hit samples",
                          flush=True)
        counts[source]={"gas_hits_processed":nvalid,
                        "invalid_gas_paths":skipped,
                        "native_rows_read":nrows}
        for axis in ("theta_u","theta_w"):
            tally=sum(a["gas_hits"] for (src,_,ax,abin),a in aggregated.items()
                      if src==source and ax==axis and abin==-1)
            if tally!=nvalid:
                raise AssertionError("Native Geant4 gas-hit exposure missing from angle bins")
        print(f"PASS {source}: {nvalid:,} first-chunk true native muon gas hits "
              f"× 2 signed angular analyses; invalid {skipped}",flush=True)
    payload={
        "meta":{
            "geometry_id":manifest["geometry_id"],
            "b_tesla":0.,
            "response_sha256":sha,
            "scope":"FIRST 3,000 GENERATED EVENTS PER CHANNEL; not full 200k",
            "coordinate_meaning":"u_birth is relative nearest physical wire; v_birth is centered chamber normal; w_birth is along anode",
            "physics":"Geant4 REAL native gas hit + independently sampled Heed ionisation + NEW symmetric Garfield electron endpoints",
            "limits":"This is birth-position versus electron-ARRIVAL time, not real u(t),v(t),w(t) drift trajectories; no avalanche/electronics, no cathode strips, no incoming reference muon normalization",
        },
        "sources":SOURCES,
        "edges_t":list(map(float,T_EDGES[:-1]))+["infinity"],
        "edges_angle":list(ANGLE_EDGES),
        "edges_coordinate":{k:list(map(float,v)) for k,v in COORD_EDGES.items()},
        "time_gates":GATES,
        "counts":counts,
        "groups":{"|".join(map(str,key)):serialize(v) for key,v in aggregated.items()},
    }
    args.out.mkdir(parents=True,exist_ok=True)
    out_json=args.out/"pi_centered_native_first_chunk_angle_timing.json"
    out_json.write_text(json.dumps(payload,separators=(",",":"))+"\n",encoding="utf-8")
    prefix="__CENTERED_TIMING_DATA__"
    tpl=args.template.read_text(encoding="utf-8")
    if tpl.count(prefix)!=1:
        raise ValueError("Centered timing HTML missing data marker")
    args.html.parent.mkdir(parents=True,exist_ok=True)
    args.html.write_text(tpl.replace(prefix,json.dumps(payload,separators=(",",":"))),
                         encoding="utf-8")
    print("PASS: centered 3+3mm Garfield & real Geant4 first-chunk"
          " origin-coordinate/arrival-time angular pilot")
    print("Raw aggregated data:",out_json)
    print("Interactive graph:",args.html)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--native",type=Path,default=NATIVE)
    parser.add_argument("--manifest",type=Path,default=TRANSPORT/"pi_centered_transport_manifest.json")
    parser.add_argument("--garfield",type=Path,default=TRANSPORT/"c_transport_box_electrons.csv")
    parser.add_argument("--heed-tracks",type=Path,
        default=REPO/"test_runs/mwpc_stage_c/c_heed_library_2p4/heed_2p4_tracks.csv")
    parser.add_argument("--heed-clusters",type=Path,
        default=REPO/"test_runs/mwpc_stage_c/c_heed_library_2p4/heed_2p4_clusters.csv")
    parser.add_argument("--out",type=Path,default=OUT)
    parser.add_argument("--html",type=Path,
        default=REPO/"test_runs/mwpc_stage_c/pi_fixed_wire_maps/pi_centered_native_timing.html")
    parser.add_argument("--template",type=Path,
        default=HERE/"pi_centered_native_timing_template.html")
    parser.add_argument("--max-per-source",type=int,default=0)
    args=parser.parse_args()
    if args.max_per_source<0:
        parser.error("--max-per-source must be nonnegative")
    try:
        build(args)
    except Exception as exc:
        print("CENTERED_TIMING_ERROR:",str(exc),file=sys.stderr,flush=True)
        raise SystemExit(1) from exc


if __name__=="__main__":
    main()
