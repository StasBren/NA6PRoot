#!/usr/bin/env python3
"""PI Stage-C: NATIVE GEANT4 + shared Heed births + six Garfield LUTs.

Physically: true first-chunk (3000 parent events per Jpsi/Phi/Omega source)
TGeo-local muon gas-hit segments from existing native_exports; fresh
independent Heed draws along EACH actual segment, reused identically
across all six p={3,3.5,4}mm x field={absent,75um} Garfield responses.
No new Garfield, ROOT, or Geant4 jobs. t is ELECTRON drift time, not
muon time or front-end discriminator time.

Includes full 10ns timing PDFs at source/station/chamber/physical-anode
level; signed muon entrance-angle bins for station/chamber. Parent-event
bootstrap uncertainties are NOT estimated in this FIRST visualization pass.
Finite physical gas widths remain exactly 307.2mm (MS0), 512mm (MS1..5);
number of centered physical anode positions varies with pitch.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
import os
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
from pi_centered_electron_sampler import REQUIRED_HIT, stable_seed
from pi_native_centered_timing_pilot import (
    FIRST_CHUNK_EXPECTED, SOURCES, GEOMETRY, ANGLES
)
from pi_pitch_fieldwire_pilot import PITCHES, tag
from run_pi_per_chamber_response import read_heed

HERE=Path(__file__).resolve().parent
REPO=HERE.parents[1]
BASE=REPO/"test_runs/mwpc_stage_c"
EXPORT=BASE/"pi_native_anode_first_chunk/native_exports"
MAPS=BASE/"pi_pitch_fieldwire_pilot_v1/pilot"
OUT=BASE/"pi_six_geometry_native_timing_v1"
CUTOFF=300
DT=10
NBINS=CUTOFF//DT+1  # index 30 is t>=300ns, not silently dropped
WIDTHS=(307.2,512.,512.,512.,512.,512.)
HEED_TRACKS=BASE/"c_heed_library_2p4/heed_2p4_tracks.csv"
HEED_CLUSTERS=BASE/"c_heed_library_2p4/heed_2p4_clusters.csv"
MODEL_VERSION="PI_SIX_GEOM_NATIVE_FIRST_CHUNK_V1"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def n_anodes(station,pitch):
    """Symmetric centered lattice, nearly full original physical aperture.

    N=round(gas_u_width/pitch) reproduces 77/128 at p=4mm.
    Edge gap is (W-(N-1)*p)/2, and must be 0<gap<=p.
    """
    if station not in range(6) or pitch not in PITCHES:
        raise ValueError("Unknown station or pitch")
    n=int(math.floor(WIDTHS[station]/pitch+0.5))
    margin=(WIDTHS[station]-(n-1)*pitch)/2
    if n<2 or not (0 < margin <= pitch+1.e-10):
        raise ValueError("Invalid finite physical anode aperture")
    return n


def signed_bin(theta):
    if not math.isfinite(theta) or not -180<=theta<=180:
        raise ValueError(f"Invalid signed entrance angle: {theta}")
    return min(len(ANGLES)-2,
               max(0,int(np.searchsorted(ANGLES,theta,side="right")-1)))


def load_map(pitch,diameter,root=MAPS):
    directory=root/tag(pitch,diameter)
    config=json.loads((directory/"frozen_config.json").read_text())
    if (config.get("mode")!="pilot" or
        config.get("pitch_mm")!=pitch or
        config.get("field_wire_diam_um")!=diameter or
        config.get("gap_minus_mm")!=3. or
        config.get("gap_plus_mm")!=3. or
        config.get("anode_voltage_v")!=1800. or
        config.get("anode_diameter_um")!=30. or
        config.get("magnetic_field_t")!=0. or
        config.get("gas")!="Ar/CO2 70:30, 293.15 K, 760 Torr" or
        config.get("du_mm")!=.1 or config.get("dv_mm")!=.2 or
        config.get("n_u")!=round(pitch/.1) or config.get("n_v")!=30):
        raise ValueError(f"Unexpected Garfield frozen physics: {directory}")
    trials=int(config["samples_per_box"])
    if trials<2:
        raise ValueError("Pilot LUT must contain >=2 trials per box")
    nboxes=int(config["n_u"])*30
    path=directory/"garfield_electrons.csv"
    # Exactly one row per box/trial; use an immutable, isolated per-geometry LUT.
    pool=np.full((nboxes,trials,5),np.nan,dtype=np.float64)
    found=np.zeros((nboxes,trials),dtype=bool)
    with path.open(newline="",encoding="utf-8") as stream:
        reader=csv.DictReader(stream)
        needed={"box_id","sample","iu","iv","collected","wire_index",
                "drift_time_ns","incidence_normal_deg","phi_hit_deg"}
        if not needed.issubset(reader.fieldnames or []):
            raise ValueError(f"Missing Garfield LUT columns: {path}")
        for row in reader:
            bid=int(row["box_id"]);j=int(row["sample"])
            if (not 0<=bid<nboxes or not 0<=j<trials or
                int(row["iu"])*30+int(row["iv"])!=bid or found[bid,j]):
                raise ValueError(f"Duplicate/bad Garfield (box,sample): {bid,j}")
            found[bid,j]=True
            yes=int(row["collected"])
            if yes not in (0,1):
                raise ValueError("Invalid Garfield collection")
            pool[bid,j,0]=yes
            pool[bid,j,1]=int(row["wire_index"]) if yes else 0
            if yes:
                t=float(row["drift_time_ns"])
                pool[bid,j,2]=t if math.isfinite(t) and t>=0 else np.nan
                pool[bid,j,3]=float(row["incidence_normal_deg"])
                pool[bid,j,4]=float(row["phi_hit_deg"])
    if not found.all():
        raise ValueError(f"Partial Garfield map: {int(found.sum())}/{found.size}")
    return pool,{"csv_sha256":sha(path),"frozen_sha256":sha(
        directory/"frozen_config.json"),"samples_per_box":trials,
        "n_u":config["n_u"],"pitch_mm":pitch,"field_wire_diam_um":diameter}


class Count:
    __slots__=("segments","primaries","collected","known","unknown","edge",
               "hist","summed","squared","incidence_sum","incidence_squared",
               "incidence_n","phi_sum","phi_squared","phi_n")
    def __init__(self):
        self.segments=0
        self.primaries=0
        self.collected=0
        self.known=0
        self.unknown=0
        self.edge=0
        self.hist=np.zeros(NBINS,dtype=np.int64)
        self.summed=self.squared=0.
        self.incidence_sum=self.incidence_squared=0.
        self.incidence_n=0
        self.phi_sum=self.phi_squared=0.
        self.phi_n=0

    def counts(self,segments=0,primaries=0,collected=0,known=0,
               unknown=0,edge=0):
        self.segments+=segments
        self.primaries+=primaries
        self.collected+=collected
        self.known+=known
        self.unknown+=unknown
        self.edge+=edge

    def add_electrons(self,times,angle,phi):
        valid=np.isfinite(times)&(times>=0)
        t=times[valid]
        self.hist+=np.bincount(np.minimum((t//DT).astype(np.int64),
                                           NBINS-1),minlength=NBINS)
        self.summed+=float(t.sum())
        self.squared+=float(np.square(t).sum())
        a=angle[np.isfinite(angle)&(angle>=0)&(angle<=180)]
        self.incidence_sum+=float(a.sum())
        self.incidence_squared+=float(np.square(a).sum())
        self.incidence_n+=len(a)
        # phi is signed in [-180,+180]; arithmetic moments are not
        # circular statistics; preserve mean sin/cos in later study.
        ph=phi[np.isfinite(phi)&(phi>=-180)&(phi<=180)]
        self.phi_sum+=float(ph.sum())
        self.phi_squared+=float(np.square(ph).sum())
        self.phi_n+=len(ph)

    def done(self):
        if self.collected != self.known+self.unknown or \
                int(self.hist.sum())!=self.known:
            raise AssertionError("Collected/time known/overflow conservation")
        def moments(s,q,n):
            if n==0:return None,None
            m=s/n
            return m,math.sqrt(max(0,q/n-m*m))
        mean,rms=moments(self.summed,self.squared,self.known)
        am,ar=moments(self.incidence_sum,self.incidence_squared,
                      self.incidence_n)
        return dict(segments=self.segments,primaries=self.primaries,
            collected=self.collected,known=self.known,unknown=self.unknown,
            edge=self.edge,hist=self.hist.tolist(),mean_ns=mean,rms_ns=rms,
            angle_n=self.incidence_n,mean_angle_deg=am,rms_angle_deg=ar)


def first_indices(u,pitch,n):
    nearest=np.floor(u/pitch+(n-1)/2+.5).astype(np.int64)
    cell=u-(nearest-(n-1)/2)*pitch
    if len(cell) and np.any(abs(cell)>pitch/2+1.e-7):
        raise AssertionError("Unphysical nearest-wire local coordinate")
    return nearest,cell


def response_for_births(u,v,q,pitch,station,pool):
    """Sample one LUT response per shared Heed electron in finite chamber."""
    n=n_anodes(station,pitch)
    nearest,cell=first_indices(u,pitch,n)
    n_u=round(pitch/.1)
    iu=np.clip(np.floor((cell+pitch/2)/.1).astype(np.int64),0,n_u-1)
    iv=np.clip(np.floor((v+3.)/.2).astype(np.int64),0,29)
    trial=np.minimum((q*pool.shape[1]).astype(np.int64),pool.shape[1]-1)
    samples=pool[iu*30+iv,trial]
    yes=samples[:,0]==1
    wire=nearest+samples[:,1].astype(np.int64)
    finite=yes&(wire>=0)&(wire<n)
    return finite,yes,wire,samples


def birth_sample(row, lam, marks, seed):
    rng=np.random.default_rng(seed)
    length=float(row["path_length_mm"])
    if not 0<length<1000 or not math.isfinite(length):
        raise ValueError("Invalid native gas-path length")
    coords=np.array([[float(row[f"{ax}_in_mm"]) for ax in "uvw"],
                     [float(row[f"{ax}_out_mm"]) for ax in "uvw"]])
    if not np.isfinite(coords).all() or abs(coords[:,1]).max()>3.05:
        raise ValueError("Invalid native gas endpoints")
    nclusters=int(rng.poisson(lam*length))
    multiplicity=marks[rng.integers(0,len(marks),size=nclusters)] \
        if nclusters else np.empty(0,dtype=int)
    t=np.repeat(rng.random(nclusters),multiplicity)
    uvw=coords[0]+t[:,None]*(coords[1]-coords[0]) \
        if len(t) else np.empty((0,3))
    quantiles=rng.random(len(uvw))
    return uvw,quantiles


def aggregate(records,root,lam,marks,limit=0):
    """Consume CSV dict rows, no change to native Geant4 or Garfield data.

    Stable source parent/event/segment seed. Group key:
     (source, config, station, det or -1, anode or -1, axis, angle-bin)
    'all' source duplicates one source histogram for overview only.
    """
    summaries=defaultdict(Count)
    processed=defaultdict(int)
    for source,row in records:
        if abs(int(row["pdg"]))!=13:
            continue
        if limit and processed[source]>=limit:
            continue
        st=int(row["station"]);det=int(row["det_id"])
        if det not in GEOMETRY or GEOMETRY[det].station!=st:
            raise ValueError("Native module geometry identity failure")
        th_u=float(row["theta_u_in_deg"])
        th_w=float(row["theta_w_in_deg"])
        bu,bw=signed_bin(th_u),signed_bin(th_w)
        seed=stable_seed("native-six-pitch-shared-Heed-v1",source,
            "chunk0001",row["event"],row["hit_index"],row["track_id"])
        positions,q=birth_sample(row,lam,marks,seed)
        nprim=len(q)
        processed[source]+=1
        u,v=positions[:,0],positions[:,1]
        for config,(pitch,diam,pool) in root.items():
            finite,rawcol,wires,samples=response_for_births(
                u,v,q,pitch,st,pool)
            idx=np.flatnonzero(finite)
            wire_good=wires[idx]
            times=samples[idx,2]
            angles=samples[idx,3]
            phis=samples[idx,4]
            known=np.isfinite(times)&(times>=0)
            edge=int(np.count_nonzero(rawcol&~finite))
            for src in (source,"all"):
                # ALL station, chamber, station-angle and chamber-angle.
                groups=[
                    (src,config,st,-1,-1,"all",-1),
                    (src,config,st,det,-1,"all",-1),
                    (src,config,st,-1,-1,"u",bu),
                    (src,config,st,-1,-1,"w",bw),
                    (src,config,st,det,-1,"u",bu),
                    (src,config,st,det,-1,"w",bw),
                ]
                for key in groups:
                    bucket=summaries[key]
                    bucket.counts(segments=1,primaries=nprim,
                        collected=len(idx),known=int(known.sum()),
                        unknown=len(idx)-int(known.sum()),edge=edge)
                    bucket.add_electrons(times,angles,phis)
                # Actual ANODE number in this pitch-dependent chamber.
                for wire in np.unique(wire_good):
                    mask=wire_good==wire
                    tt=times[mask]
                    key=(src,config,st,det,int(wire),"all",-1)
                    bucket=summaries[key]
                    bucket.counts(segments=1,collected=int(mask.sum()),
                        known=int(np.count_nonzero(np.isfinite(tt)&(tt>=0))),
                        unknown=int(mask.sum())-int(np.count_nonzero(
                            np.isfinite(tt)&(tt>=0))))
                    bucket.add_electrons(tt,angles[mask],phis[mask])
    return summaries,dict(processed)


def validate_conservation(groups):
    """Every physical-anode electron counted once; every chamber in station.

    This is a strict global sanity gate on native finite-anode assignments,
    NOT a check of incoming reference-plane muon efficiency.
    """
    per_chamber=defaultdict(lambda:[0,0,np.zeros(NBINS,dtype=np.int64)])
    per_station=defaultdict(lambda:[0,0,np.zeros(NBINS,dtype=np.int64)])
    for key,b in groups.items():
        src,config,st,det,anode,axis,anglebin=key
        if axis!="all" or anglebin!=-1:
            continue
        if anode>=0:
            k=(src,config,st,det)
            slot=per_chamber[k]
            slot[0]+=b.collected
            slot[1]+=b.known
            slot[2]+=b.hist
        elif det>=0:
            k=(src,config,st)
            slot=per_station[k]
            slot[0]+=b.collected
            slot[1]+=b.known
            slot[2]+=b.hist
    for key,b in groups.items():
        src,config,st,det,anode,axis,anglebin=key
        if axis!="all" or anglebin!=-1:
            continue
        if det>=0 and anode==-1:
            observed=per_chamber[(src,config,st,det)]
        elif det==-1 and anode==-1:
            observed=per_station[(src,config,st)]
        else:
            continue
        if (b.collected!=observed[0] or b.known!=observed[1]
                or not np.array_equal(b.hist,observed[2])):
            raise AssertionError(
                f"Station/chamber/anode finite-wire electron conservation "
                f"failed: {key}, counted {observed[0]} vs {b.collected}")
    return True


def write_report(out,groups,processed,meta):
    validate_conservation(groups)
    out.mkdir(parents=True,exist_ok=True)
    response=[]
    for key,bucket in sorted(groups.items()):
        source,config,station,det,anode,axis,anglebin=key
        row=dict(source=source,config=config,station=station,det_id=det,
                 anode=anode,axis=axis,anglebin=anglebin,**bucket.done())
        response.append(row)
    # Save analysis as machine-readable JSON (browser embeds subset).
    payload=dict(meta=meta,group_count=len(response),groups=response)
    target=out/"pi_six_geometry_native_timing.json"
    temp=out/"pi_six_geometry_native_timing.json.tmp"
    temp.write_text(json.dumps(payload,ensure_ascii=False,
                  separators=(",",":"))+"\n",encoding="utf-8")
    os.replace(temp,target)
    print("PASS: first-chunk shared-Heed six-geometry histograms")
    print("Native source gas segments:",processed)
    print("Groups:",len(response),"output:",target,flush=True)
    return payload


def run(args):
    t0=time.perf_counter()
    for f in (args.export_dir,args.maps,args.heed_tracks,args.heed_clusters):
        if not f.exists():
            raise FileNotFoundError(f)
    maps={}
    provenance={}
    for pitch in PITCHES:
        for diameter in (0.,75.):
            name=tag(pitch,diameter)
            pool,meta=load_map(pitch,diameter,args.maps)
            maps[name]=(pitch,diameter,pool)
            provenance[name]=meta
    lam,marks=read_heed(args.heed_tracks,args.heed_clusters)
    def native_rows():
        for source,ch in SOURCES.items():
            path=args.export_dir/f"{ch}_geant_gas_crossings.csv"
            with path.open(newline="",encoding="utf-8") as file:
                rdr=csv.DictReader(file)
                if not REQUIRED_HIT.issubset(rdr.fieldnames or ()):
                    raise ValueError("Missing native Geant4 gas-hit columns")
                for row in rdr:
                    yield source,row
    groups,counts=aggregate(native_rows(),maps,lam,marks,
                            limit=args.limit_per_source)
    if args.limit_per_source == 0:
        if set(counts)!=set(SOURCES):
            raise ValueError("Full native first-chunk missing a source")
        for source,n in counts.items():
            if n!=FIRST_CHUNK_EXPECTED[source]:
                raise ValueError(f"Native first-chunk count mismatch: {source} {n}")
    meta={
        "model_version":MODEL_VERSION,
        "data_scope":"NATIVE FIRST CHUNK (3000 generated parent events per source)"
                     if args.limit_per_source==0 else
                     f"EXPLORATORY first {args.limit_per_source} muon gas "
                     "segments per source, NOT full native first chunk",
        "gas_hit_counts":counts,
        "native_geant4_export_sha256":{
            source:sha(args.export_dir/f"{channel}_geant_gas_crossings.csv")
            for source,channel in SOURCES.items()},
        "heed_clusters_per_mm":lam,
        "heed_tracks_sha":sha(args.heed_tracks),
        "heed_clusters_sha":sha(args.heed_clusters),
        "garfield_luts":provenance,
        "physical_gas_width_mm":list(WIDTHS),
        "finite_anodes_by_pitch":{str(p):[n_anodes(st,p)
                                        for st in range(6)] for p in PITCHES},
        "wire_grid_rule":"centered finite (N-1) wire positions, "
                         "N=floor(W/p+0.5), physical gas aperture W fixed",
        "time_edges_ns":list(range(0,CUTOFF+1,DT))+["overflow"],
        "time_bin_width_ns":DT,
        "angle_bins_deg":ANGLES.tolist(),
        "angle_semantics":"signed muon entrance theta_u,theta_w Geant4; "
                 "electron incidence_normal_deg from LAST microscopic "
                 "Garfield collision step; phi_hit not circular-averaged",
        "input_model":"Heed cluster density + empirical Heed marks; "
                "same birth coordinates/quantiles all six Garfield maps. "
                "Single-electron LUT sampling (20 trial/cell), not "
                "direct Garfield drift for each native Heed electron.",
        "simulation_limitations":"B=0 ideal analytic field, no avalanche "
                  "gain/ion current/ASIC threshold/cathode strips; NOT "
                  "a proof of full detector timing efficiency",
        "statistical_errors":"NOT YET REPORTED: parent-event cluster bootstrap "
                  "required before PI-facing error claims. Avoid treating "
                  "individual electrons as independent source events.",
        "wall_seconds":time.perf_counter()-t0,
    }
    return write_report(args.out,groups,counts,meta)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--export-dir",type=Path,default=EXPORT)
    parser.add_argument("--maps",type=Path,default=MAPS)
    parser.add_argument("--heed-tracks",type=Path,default=HEED_TRACKS)
    parser.add_argument("--heed-clusters",type=Path,default=HEED_CLUSTERS)
    parser.add_argument("--out",type=Path,default=OUT)
    parser.add_argument("--limit-per-source",type=int,default=0,
                        help="Smoke with 200; full native first chunk with 0")
    args=parser.parse_args()
    if args.limit_per_source < 0:
        parser.error("limit-per-source must be >= 0")
    run(args)


if __name__=="__main__":
    main()
