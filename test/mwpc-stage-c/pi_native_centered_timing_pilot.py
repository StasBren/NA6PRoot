#!/usr/bin/env python3
"""Native first-chunk 3+3mm centered anode primary-electron timing pilot.

True first 3000 generated events/channel: raw MCKine-joined Geant4 gas-hit
entry/exit (including local w and signed entrance muon angles).
Ionisation: independently sampled Heed clusters. Transport: verified NEW
Garfield centered 3+3mm map, NOT LEGACY 2+4mm.
Real no-electron gas segments remain in N0. Histogram PDFs are conditional
on collected primary electrons with known drift time. This is not analog
electronics, cannot yield real intermediate electron positions u(t),v(t),w(t).
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

from pi_centered_electron_sampler import (
    GEOMETRY_ID, REQUIRED_HIT, stable_seed, load_centered_map,
    ionize_and_drift, multiplicities,
)
from run_pi_per_chamber_response import read_heed
from mwpc_checkerboard_projection import build_chambers

HERE=Path(__file__).resolve().parent
REPO=HERE.parents[1]
BASE=REPO/"test_runs/mwpc_stage_c"
EXPORT=BASE/"pi_native_anode_first_chunk/native_exports"
CENTERED=BASE/"c_transport_box_map_gm3p0_gp3p0"
OUT=BASE/"pi_centered_anode_timing_pilot"
HTML=BASE/"pi_fixed_wire_maps/pi_centered_timing_pilot.html"
SOURCES={"jpsi40-bore":"jpsi","phi40-bore":"phi","omega40-bore":"omega"}
FIRST_CHUNK_EXPECTED={"jpsi40-bore":29488,"phi40-bore":23301,
                      "omega40-bore":21533}
GEOMETRY={c.det_id:c for group in build_chambers() for c in group}
ANGLES=np.array([-180,-30,-10,-5,-2,0,2,5,10,30,180],dtype=float)
CUTOFFS=(50,100,200,400,800)
WINDOW_NAMES=("50","100","200","400","800","all")
TIME_EDGES=np.array(list(range(0,301,10))+[np.inf],dtype=float)
# Birth coordinate u relative to nearest infinite wire (0 at anode),
# chamber-local v (normal), w (along wires), and full local u.
COORDINATES={
    "u_cell_mm":np.linspace(-2,2,9),
    "v_birth_mm":np.linspace(-3,3,13),
    "w_birth_mm":np.linspace(-350,350,29),
    "u_birth_mm":np.linspace(-300,300,31),
}


def angle_bin(theta):
    if not math.isfinite(theta) or not -180<=theta<=180:
        raise ValueError(f"Invalid signed native local entrance angle {theta}")
    return min(len(ANGLES)-2,max(0,int(np.searchsorted(
        ANGLES,theta,side="right")-1)))


def fresh(module=-1):
    return {
        "module":module,"n":0,"nprim":0,"collected":0,"known":0,
        "unknown":0,"edge":0,"n0":0,
        "m":np.zeros((len(WINDOW_NAMES),6),dtype=np.int64),
        "time":np.zeros(len(TIME_EDGES)-1,dtype=np.int64),
        "coords":({name:np.zeros((len(edge)-1,len(TIME_EDGES)-1),
                  dtype=np.int64) for name,edge in COORDINATES.items()}
                  if module==-1 else None),
    }


def update(bucket,hit,nw,tbins,valid,coordinate_bins):
    """One physical gas segment, NOT a station reference-plane track."""
    bucket["n"]+=1
    bucket["nprim"]+=hit["nprimary"]
    bucket["collected"]+=hit["ncollected"]
    bucket["edge"]+=hit["nedge"]
    bucket["unknown"]+=int((~valid).sum())
    bucket["known"]+=int(valid.sum())
    bucket["n0"]+=hit["ncollected"]==0
    for i,k in enumerate(nw):
        bucket["m"][i,min(k,5)]+=1
    if np.any(valid):
        bucket["time"]+=np.bincount(tbins[valid],minlength=len(TIME_EDGES)-1)
        if bucket["coords"] is not None:
            for name,positions in coordinate_bins.items():
                arr=bucket["coords"][name]
                ids=positions[valid]
                hist=np.bincount(ids*(len(TIME_EDGES)-1)+tbins[valid],
                                 minlength=arr.size).reshape(arr.shape)
                arr+=hist


def finish(bucket):
    n=bucket["n"]
    if not np.all(bucket["m"].sum(axis=1)==n):
        raise AssertionError("Zero-electron inclusive gas-hit multiplicity mismatch")
    if np.any(bucket["m"]<0) or int(bucket["time"].sum())!=bucket["known"]:
        raise AssertionError("Electron drift-time histogram conservation failed")
    if bucket["known"]+bucket["unknown"]!=bucket["collected"]:
        raise AssertionError("Collected endpoint/missing-time disagreement")
    if bucket["coords"] is not None:
        for name,arr in bucket["coords"].items():
            if int(arr.sum())!=bucket["known"]:
                raise AssertionError("Coordinate-time histogram lost electrons: "+name)
    return dict(
        n=n,nprim=bucket["nprim"],collected=bucket["collected"],
        known=bucket["known"],unknown=bucket["unknown"],
        edge=bucket["edge"],n0=bucket["n0"],
        m=bucket["m"].tolist(),time=bucket["time"].tolist(),
        coords={k:v.tolist() for k,v in bucket["coords"].items()}
                if bucket["coords"] is not None else None,
    )


def run(args):
    for p in (args.export_dir,args.garfield,args.manifest,
              args.heed_tracks,args.heed_clusters,args.template):
        if not p.exists():
            raise FileNotFoundError(f"Missing pipeline input: {p}")
    pool,field=load_centered_map(args.garfield,args.manifest)
    lam,marks=read_heed(args.heed_tracks,args.heed_clusters)
    containers={}
    hit_totals={}
    args.out.mkdir(parents=True,exist_ok=True)
    detailed=args.out/"pi_centered_native_hit_multiplicities.csv"
    columns=("source","chunk","event","hit_index","track_id",
             "station","det_id","theta_u_in_deg","theta_w_in_deg",
             "angle_bin_u","angle_bin_w","gas_path_mm","edep_eV",
             "nclusters","nprimary","ncollected","n_uncollected",
             "n_lost_finite_anode","n_missing_time",
             "N_t50","N_t100","N_t200","N_t400","N_t800","N_all")
    with detailed.open("w",newline="",encoding="utf-8") as stream:
        writer=csv.DictWriter(stream,fieldnames=columns)
        writer.writeheader()
        for source,channel in SOURCES.items():
            raw=args.export_dir/f"{channel}_geant_gas_crossings.csv"
            if not raw.is_file():
                raise FileNotFoundError("Missing first-chunk native ROOT export: "+str(raw))
            n=0
            with raw.open(newline="",encoding="utf-8") as inp:
                reader=csv.DictReader(inp)
                if not REQUIRED_HIT.issubset(reader.fieldnames or ()):
                    raise ValueError(f"{raw}: missing true TGeo-local fields: "
                                     f"{REQUIRED_HIT-set(reader.fieldnames or ())}")
                for row in reader:
                    pdg=int(row["pdg"])
                    if abs(pdg)!=13:
                        continue
                    st=int(row["station"]);det=int(row["det_id"])
                    if det not in GEOMETRY or GEOMETRY[det].station!=st:
                        raise ValueError("Native 298-module detector/station map mismatch")
                    event=int(row["event"])
                    index=int(row["hit_index"])
                    track=int(row["track_id"])
                    if not (0<=event<3000 and index>=0 and track>=0):
                        raise ValueError("Incorrect native first-chunk event/hit index")
                    seed=stable_seed("centered_native_pilot_v1",source,
                                     "chunk0001",event,index,track)
                    rng=np.random.default_rng(seed)
                    hit=ionize_and_drift(row,rng,lam,marks,pool)
                    mult=multiplicities(hit,CUTOFFS)
                    th_u=float(row["theta_u_in_deg"])
                    th_w=float(row["theta_w_in_deg"])
                    iu=angle_bin(th_u);iw=angle_bin(th_w)
                    valid=hit["valid_times"]
                    # Individual primary electrons only; t=0 is not used
                    # for a missing drift-time trial.
                    tbins=np.clip(np.searchsorted(
                        TIME_EDGES,hit["times"],side="right")-1,
                        0,len(TIME_EDGES)-2).astype(np.int64)
                    cbins={
                        "u_cell_mm":hit["u_cell"],
                        "v_birth_mm":hit["v"],
                        "w_birth_mm":hit["w"],
                        "u_birth_mm":hit["u"],
                    }
                    cbins={
                        name:np.clip(np.searchsorted(edges,v,side="right")-1,
                           0,len(edges)-2).astype(np.int64)
                        for name,(edges,v) in (
                            (name,(COORDINATES[name],values))
                            for name,values in cbins.items()
                        )
                    }
                    # Two independent signed angular axes, always with an
                    # angle-integrated bin (-1) and one signed narrow bin.
                    for axis,abin in (("u",iu),("w",iw)):
                        for module in (-1,det):
                            for idx in (-1,abin):
                                key=(source,st,module,axis,idx)
                                bucket=containers.setdefault(key,fresh(module))
                                update(bucket,hit,mult,tbins,valid,cbins)
                    writer.writerow(dict(
                        source=source,chunk="chunk0001",event=event,hit_index=index,
                        track_id=track,station=st,det_id=det,
                        theta_u_in_deg=th_u,theta_w_in_deg=th_w,
                        angle_bin_u=iu,angle_bin_w=iw,
                        gas_path_mm=row["path_length_mm"],edep_eV=row["edep_eV"],
                        nclusters=hit["nclusters"],nprimary=hit["nprimary"],
                        ncollected=hit["ncollected"],
                        n_uncollected=hit["nuncollected"],
                        n_lost_finite_anode=hit["nedge"],
                        n_missing_time=int((~valid).sum()),
                        **{f"N_t{cutoff}":k for cutoff,k in zip(CUTOFFS,mult[:-1])},
                        N_all=mult[-1],
                    ))
                    n+=1
                    if n%5000==0:
                        print(f"{source}: {n:,} native gas-hit segments processed",flush=True)
            if n!=FIRST_CHUNK_EXPECTED[source]:
                raise ValueError(f"{source}: found {n:,} native muon gas-hit segments,"
                                 f" frozen first-chunk truth is "
                                 f"{FIRST_CHUNK_EXPECTED[source]:,}")
            hit_totals[source]=n
            print(f"PASS {source}: {n:,} MCKine-confirmed true native muon gas segments",flush=True)
    payload_groups={}
    sumrows=[]
    for (src,st,mod,axis,bin_index),bucket in sorted(containers.items()):
        key="|".join(map(str,(src,st,mod,axis,bin_index)))
        encoded=finish(bucket)
        payload_groups[key]=encoded
        n=encoded["n"]
        for i,window in enumerate(WINDOW_NAMES):
            counts=encoded["m"][i]
            sumrows.append(dict(source=src,station=st,det_id=mod,
                angle_axis=axis,angle_bin=bin_index,window_ns=window,
                gas_hit_segments=n,N0=counts[0],N1=counts[1],
                N2=counts[2],N3=counts[3],N4=counts[4],
                N5plus=counts[5],P_ge2=sum(counts[2:])/n if n else "",
                known_electron_drift_times=encoded["known"],
                unknown_electron_drift_times=encoded["unknown"],
                collected_electrons=encoded["collected"]))
    with (args.out/"pi_centered_time_angle_multiplicity.csv").open(
        "w",newline="",encoding="utf-8"
    ) as stream:
        writer=csv.DictWriter(stream,fieldnames=list(sumrows[0]))
        writer.writeheader();writer.writerows(sumrows)
    # Four ways of updating the same native hit (two angle axes ×
    # all-inclusive bin + one signed bin). The station & module integrated
    # bins across axes must independently reproduce frozen hit counts.
    for src in SOURCES:
        for st in range(6):
            for axis in ("u","w"):
                station=containers[(src,st,-1,axis,-1)]["n"]
                allmods=sum(containers[key]["n"] for key in containers
                        if key[0]==src and key[1]==st and key[2]>=0
                        and key[3]==axis and key[4]==-1)
                if station!=allmods:
                    raise AssertionError("Physical module-to-station hit conservation")
        summed=sum(containers[(src,st,-1,"u",-1)]["n"] for st in range(6))
        if summed!=hit_totals[src]:
            raise AssertionError("Channel true gas-hit conservation failed")
    payload={
        "meta":{
            "scope":"FIRST 3000 generated events/channel, ALL TRUE native muon gas hits",
            "geometry_id":GEOMETRY_ID,"map_sha256":field["response_csv_sha256"],
            "model":"Fresh symmetric 3+3mm Garfield B0; Heed independent of Geant4 Edep",
            "cutoffs_ns":list(CUTOFFS)+["all"],
            "sample":"MCKine-joined Geant4 native gas HIT SEGMENTS, not station incoming muons",
            "signal":"anode with >=1 collected PRIMARY electron, not analog threshold",
            "times":"t_drift relative to electron birth, NOT Geant4 hit time or ASIC trigger time",
            "coords":"Birth coordinates only; these are NOT true u(t),v(t),w(t) trajectories",
        },
        "sourceNames":list(SOURCES),
        "angles":ANGLES.tolist(),"tEdges":TIME_EDGES[:-1].tolist()+["overflow"],
        "coordinates":{key:edges.tolist() for key,edges in COORDINATES.items()},
        "geometry":{str(det):c.station for det,c in GEOMETRY.items()},
        "groups":payload_groups,
    }
    marker="__CENTERED_PILOT_PAYLOAD__"
    template=args.template.read_text(encoding="utf-8")
    if template.count(marker)!=1:
        raise ValueError("Incorrect centered pilot HTML marker")
    args.html.parent.mkdir(parents=True,exist_ok=True)
    args.html.write_text(template.replace(marker,json.dumps(payload,
        ensure_ascii=False,separators=(",",":")).replace("</","<\\/")),
        encoding="utf-8")
    manifest={
        "geometry_id":GEOMETRY_ID,
        "scope":"NATIVE FIRST CHUNK ONLY, 3000 events/source",
        "native_muon_hits":hit_totals,
        "sha256_centered_garfield":field["response_csv_sha256"],
        "sha256_heed_tracks":hashlib.sha256(args.heed_tracks.read_bytes()).hexdigest(),
        "sha256_heed_clusters":hashlib.sha256(args.heed_clusters.read_bytes()).hexdigest(),
        "sha256_code":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "model":"Heed PRIMARY electron + verified symmetric Garfield drift; no analog electronics",
        "no_cathode_strips":True,
    }
    (args.out/"pi_centered_timing_manifest.json").write_text(
        json.dumps(manifest,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print("PASS: NEW centered 3+3mm Geant4/Heed/Garfield first-chunk pilot")
    print("Native true muon gas hits:",hit_totals)
    print("Open:",args.html)
    print("PI time×angle×module CSV:",args.out)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--export-dir",type=Path,default=EXPORT)
    ap.add_argument("--garfield",type=Path,
                    default=CENTERED/"c_transport_box_electrons.csv")
    ap.add_argument("--manifest",type=Path,
                    default=CENTERED/"pi_centered_transport_manifest.json")
    ap.add_argument("--heed-tracks",type=Path,
                    default=BASE/"c_heed_library_2p4/heed_2p4_tracks.csv")
    ap.add_argument("--heed-clusters",type=Path,
                    default=BASE/"c_heed_library_2p4/heed_2p4_clusters.csv")
    ap.add_argument("--out",type=Path,default=OUT)
    ap.add_argument("--html",type=Path,default=HTML)
    ap.add_argument("--template",type=Path,
                    default=HERE/"pi_centered_timing_template.html")
    args=ap.parse_args()
    run(args)


if __name__=="__main__":
    main()
