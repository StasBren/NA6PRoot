#!/usr/bin/env python3
"""Read-only 1000-muon Penning-ON timing/avalanche correlation cache.

IMPORTANT: t10 is the 10th percentile of ion BIRTH times saved from
<=128 uniformly sampled ionisations in one Garfield avalanche, relative
to the original MUON ENTRY time (Geant4/Heed t=0). It is NOT an independently
measured electron ARRIVAL time, drift time, or incoming angle.

This script NEVER runs Garfield, modifies sealed source data, or accesses
the ongoing 30k muon campaign. The new auxiliary NPZ is independent of
the existing 12-plot Streamlit dashboard cache.
"""
from __future__ import annotations
import csv
import hashlib
import json
import math
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

REPO=Path(__file__).resolve().parents[2]
ROOT=REPO/"test_runs/mwpc_stage_c/pi_p3_noFW_jpsi1000_penning057_v1"
DIR=ROOT/"interactive_dashboard_v1"
ORIGINAL=DIR/"dashboard_data.npz"
ORIGMETA=DIR/"provenance.json"
CACHE=DIR/"temporal_correlations_v1.npz"
META=DIR/"temporal_correlations_v1.json"
PITCH=3.
PHI_EDGES=np.linspace(-180.,180.,73)

KINDS={
  "t10_gain":"Mean gain vs early ion-birth time",
  "t10_asym":"Mean |ion asymmetry| vs early ion-birth time",
  "t10_angles":"Ion-birth angular distributions: early vs late",
  "t10_start":"Early ion-birth time vs initial |v|",
  "t10_gain_position":"Gain conditioned on early time and initial |v|",
  "t10_timeline":"One real muon's cluster/avalanche timeline",
}
DEFAULTS={
  "t10_gain":("Avalanche gain vs early ion-birth time",
              "Early ion-birth time since muon entry [ns]",
              "Mean collected electrons per starting electron",
              "Mean collected gain"),
  "t10_asym":("Ion-birth asymmetry vs early ion-birth time",
              "Early ion-birth time since muon entry [ns]",
              "Mean |A_ion,hemi|",
              "Mean |A|"),
  "t10_angles":("Ion-birth azimuth in early and late avalanches",
                "Birth azimuth around nearest anode [degrees]",
                "Ion-weighted probability per angular bin",
                "Ion births"),
  "t10_start":("Early ion birth time versus initial |v|",
               "Initial electron |v| [mm]",
               "Early ion-birth time since muon entry [ns]",
               "Heed electron starts"),
  "t10_gain_position":("Gain vs early ion birth time and starting |v|",
                       "Early ion-birth time since muon entry [ns]",
                       "Initial electron |v| [mm]",
                       "Mean collected gain"),
  "t10_timeline":("One muon: cluster and avalanche chronology",
                  "Time since muon gas entry [ns]",
                  "Original TrackHeed cluster index",
                  "Electron collection and exact ion asymmetry"),
}

def must(x,msg):
    if not x:raise RuntimeError(msg)
def sha(path):
    d=hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda:f.read(1<<20),b""):
            d.update(block)
    return d.hexdigest()
def rows(path,required):
    with Path(path).open(newline="") as f:
        r=csv.DictReader(f)
        must(r.fieldnames is not None and
             len(set(r.fieldnames))==len(r.fieldnames) and
             set(required).issubset(r.fieldnames),
             "Incorrect raw native CSV schema "+str(path))
        yield from r
def seed_key(row):
    return (row["crossing_id"],int(row["primary_index"]))
def wire_offset(u):
    n=math.floor(u/PITCH+.5) if u>=0 else math.ceil(u/PITCH-.5)
    return u-PITCH*n
def phi_bin(u,v):
    theta=math.degrees(math.atan2(v,wire_offset(u)))
    return max(0,min(71,int((theta+180.)/5.)))
def sample_quantiles(tt):
    t=np.asarray(tt,dtype=float)
    must(t.size>0 and np.isfinite(t).all(),"Invalid sampled physical birth time")
    return np.quantile(t,[.10,.50,.90])
def validate_original_cache():
    must(ORIGINAL.is_file() and ORIGMETA.is_file(),
         "Run the original dashboard 'prepare' once before timing prepare")
    meta=json.loads(ORIGMETA.read_text())
    must(meta.get("schema")=="P3_NOFW_PENNING_ON_INTERACTIVE_DASHBOARD_V1"
         and meta.get("n_independent_Geant4_Jpsi_parents")==1000
         and sha(ORIGINAL)==meta.get("cache_sha256"),
         "Verified original Penning-ON dashboard cache invalid")
    return meta

def prepare():
    source=validate_original_cache()
    must(not CACHE.exists() and not META.exists(),
         "Refusing to overwrite auxiliary timing-cache data: "+str(CACHE))
    variables=("parent","primary","cluster","gain","ions","A",
               "u0","v0","t_heed","t10","t50","t90","span","sample_count")
    out={field:[] for field in variables}
    angular=[]
    muons=[]
    ions_all=0
    points_all=0
    n_bad_time=0
    for block in range(50):
        d=ROOT/"native_chunks"/f"chunk_{block:04d}"
        sealed=d/"complete.json"
        must(sealed.is_file(),"Missing sealed Garfield production checkpoint "+str(d))
        rec=json.loads(sealed.read_text())
        must(rec["chunk"]==block and rec["n_parents"]==20 and
             rec["stats"]["n_capped_seeds"]==0 and
             rec["stats"]["same_Heed_electron_positions_and_Garfield_seeds"] is True,
             "Garfield physics seal/paired-Heed validation failed")
        for file in ("muon_reference_crossing.csv","muon_reference_heed_primaries.csv",
                     "muon_reference_heed_clusters.csv","muon_reference_seed.csv",
                     "muon_reference_ion_birth_sample.csv"):
            must(rec["file_sha256"].get(file)==sha(d/file),
                 "SHA mismatch in immutable Penning-ON "+str(d/file))
        crossings=list(rows(d/"muon_reference_crossing.csv",(
            "source|chunk|event|hit_index|track_id|station|det_id",
            "heed_primaries","heed_clusters",
            "theta_u_in_deg","theta_w_in_deg","path_mm",
            "collected_avalanche_electrons","ion_births")))
        must(len(crossings)==20,"Penning-on chunk does not contain 20 muons")
        index={}
        for r in crossings:
            name=r["source|chunk|event|hit_index|track_id|station|det_id"]
            must(name not in index,"Duplicate real Geant4 parent")
            index[name]=len(muons)
            muons.append({
                "id":name,
                "n_Heed":int(r["heed_primaries"]),
                "n_clusters":int(r["heed_clusters"]),
                "theta_u_deg":float(r["theta_u_in_deg"]),
                "theta_w_deg":float(r["theta_w_in_deg"]),
                "path_mm":float(r["path_mm"]),
                "collected_e":int(r["collected_avalanche_electrons"]),
                "ion_births":int(r["ion_births"]),
            })
        heed={}
        n_seen=defaultdict(int)
        for r in rows(d/"muon_reference_heed_primaries.csv",(
            "crossing_id","primary_index","cluster_index","t_ns")):
            k=seed_key(r)
            must(k not in heed,"Duplicate original Heed primary electron")
            heed[k]=(int(r["cluster_index"]),float(r["t_ns"]))
            n_seen[k[0]]+=1
        cluster_count=defaultdict(int)
        for r in rows(d/"muon_reference_heed_clusters.csv",(
            "crossing_id","cluster_index","n_accepted","t_ns")):
            cluster_count[r["crossing_id"]]+=int(r["n_accepted"])
        for r in crossings:
            name=r["source|chunk|event|hit_index|track_id|station|det_id"]
            must(n_seen.get(name,0)==int(r["heed_primaries"])==
                 cluster_count.get(name,0),
                 "Original Heed clusters and primaries inconsistent")
        seeds=list(rows(d/"muon_reference_seed.csv",(
            "crossing_id","primary_index","seed_u_mm","seed_v_mm",
            "collected_e","ion_births","ion_up","ion_down","capped")))
        lookup={}
        must(len(seeds)==len(heed),
             "Garfield seed count and original Heed population differ")
        for j,r in enumerate(seeds):
            k=seed_key(r)
            must(k in heed and k[0] in index and k not in lookup,
                 "Garfield primary identity invalid")
            lookup[k]=j
        # The microscopic ion-birth CSV can be interspersed with ionless
        # avalanches. Group by exact (crossing_id,primary_index); never
        # assume every seed has a reservoir or that rows are dense.
        samples=defaultdict(list)
        angle_counts=defaultdict(lambda:np.zeros(72,dtype=np.float32))
        weights=defaultdict(float)
        for r in rows(d/"muon_reference_ion_birth_sample.csv",(
            "crossing_id","primary_index","u_mm","v_mm",
            "t_ns","reservoir_weight","seed_ion_births_total")):
            k=seed_key(r)
            must(k in lookup,"Unmatched sampled Garfield ion")
            num=int(seeds[lookup[k]]["ion_births"])
            m=min(128,num)
            must(num>0 and int(r["seed_ion_births_total"])==num,
                 "Sampled ion birth has invalid avalanche ion census")
            w=float(r["reservoir_weight"])
            must(math.isfinite(w) and abs(w-num/m)<1.e-7*max(1.,num/m),
                 "Incorrect reservoir weight or sample size")
            t=float(r["t_ns"])
            must(math.isfinite(t),"Nonfinite microscopic physical time")
            samples[k].append(t)
            angle_counts[k][phi_bin(float(r["u_mm"]),float(r["v_mm"]))]+=1.
            weights[k]+=w
        must(sum(int(r["ion_births"]) for r in seeds)==
             sum(int(r["ion_births"]) for r in crossings),
             "Exact native seed and real muon ion census mismatch")
        must(sum(int(r["collected_e"]) for r in seeds)==
             sum(int(r["collected_avalanche_electrons"]) for r in crossings),
             "Exact Garfield collected charges disagree with muon sums")
        for r in seeds:
            k=seed_key(r)
            p=index[k[0]]
            cl,th=heed[k]
            n=int(r["ion_births"])
            gain=int(r["collected_e"])
            up,down=int(r["ion_up"]),int(r["ion_down"])
            must(min(n,gain,up,down)>=0 and int(r["capped"])==0,
                 "Invalid/capped microscopic avalanche")
            m=min(128,n)
            t=samples.get(k,[])
            must(len(t)==m,
                 "Sampled ion birth reservoir does not match exact ion census")
            must(abs(weights.get(k,0.)-n)<1.e-6*max(1.,n),
                 "Reservoir total does not reproduce exact ion birth count")
            if n:
                q=sample_quantiles(t)
                if q[0]<th-1.e-6:n_bad_time+=1
                onset,p50,p90=map(float,q)
                span=float(max(t)-min(t))
                fraction=angle_counts[k]/float(m)
                must(abs(float(fraction.sum())-1.)<1e-5,
                     "Angular ion-birth samples do not conserve unit weight")
                ions_all+=n
                points_all+=m
            else:
                onset=p50=p90=span=float("nan")
                fraction=np.zeros(72,dtype=np.float32)
            A=(up-down)/(up+down) if up+down else float("nan")
            vals=dict(parent=p,primary=k[1],cluster=cl,gain=gain,ions=n,
                A=A,u0=wire_offset(float(r["seed_u_mm"])),
                v0=float(r["seed_v_mm"]),t_heed=th,
                t10=onset,t50=p50,t90=p90,
                span=span,sample_count=m)
            for field in variables:out[field].append(vals[field])
            angular.append(fraction)
        del samples,angle_counts,weights
        print(f"PASS sealed Penning ON {block+1}/50 muon chunks; "
              f"{len(seeds)} original Garfield avalanches",flush=True)
    must(len(muons)==1000 and len(out["gain"])>20_000,
         "1000 real muons/Heed electron population missing")
    must(n_bad_time==0,
         f"{n_bad_time} sampled early birth times predate the initiating Heed electron")
    arrays={}
    for field,values in out.items():
        dtype=(np.int32 if field in ("parent","primary","cluster","gain",
                                    "ions","sample_count") else np.float64)
        arrays[field]=np.asarray(values,dtype=dtype)
    arrays["angle_fraction"]=np.asarray(angular,dtype=np.float32)
    n=len(arrays["gain"])
    must(arrays["angle_fraction"].shape==(n,72),
         "Incorrect per-avalanche angular sample dimensions")
    with np.load(ORIGINAL,allow_pickle=False) as source_data:
        must(np.array_equal(source_data["gain"],arrays["gain"]) and
             np.array_equal(source_data["births"],arrays["ions"]),
             "New timing data do not match existing verified Penning gain/birth arrays")
        must(len(source_data["muon_heed"])==len(muons) and
             np.array_equal(source_data["muon_heed"],
                            np.asarray([x["n_Heed"] for x in muons])),
             "New timing muon identity order differs from existing dashboard")
    must(sum(x["ion_births"] for x in muons)==ions_all,
         "Full 1000-muon ion birth total disagrees with 128-point reservoirs")
    must(ions_all==int(source["n_true_ion_births"]),
         "New temporal ion census disagrees with old verified dashboard")
    np.savez_compressed(CACHE,**arrays)
    metadata={
        "schema":"P3_NOFW_PENNING_ON_TEMPORAL_CORRELATIONS_V1",
        "muons":muons,
        "n_independent_true_Geant4_muons":1000,
        "n_exact_Garfield_avalanches":n,
        "n_true_ion_births":ions_all,
        "n_sampled_ion_birth_positions":points_all,
        "recommended_sample_count_threshold":128,
        "t10_definition":"10th percentile of <=128 uniformly sampled ion-birth timestamps after original muon entry; NOT electron drift or arrival time",
        "physical_setup":"Penning ON r=0.57, p3_noFW, pitch=3mm, +/-3mm cathodes, 1800V",
        "original_dashboard_cache_sha256":sha(ORIGINAL),
        "auxiliary_cache_sha256":sha(CACHE),
    }
    META.write_text(json.dumps(metadata,indent=2)+"\n")
    print(f"PASS new immutable temporal cache: 1000 muons, "
          f"{n:,} avalanches, {ions_all:,} exact ion births, "
          f"{points_all:,} sampled positions")
    print("NEW CACHE",CACHE,flush=True)

def selftest():
    must(phi_bin(0.,1.)==54 and phi_bin(0.,-1.)==18
         and abs(wire_offset(3.2)-.2)<1e-12,
         "Time study nearest-wire angular geometry incorrect")
    q=sample_quantiles([1.,2.,3.,4.])
    must(np.allclose(q,[1.3,2.5,3.7]),
         "Physical 10/50/90 percentile estimator incorrect")
    # A zero-ion avalanche has no birth samples and no defined onset.
    ions=0; onset=float("nan") if ions==0 else 0.
    must(np.isnan(onset),"Zero-ion onset must remain undefined")
    print("PASS nearest-anode, sampled onset, zero-ion selftests; "
          "NO Garfield/raw data writes.")

def load():
    must(CACHE.is_file() and META.is_file(),
         "Temporal cache missing: run timing-correlations script 'prepare' first")
    m=json.loads(META.read_text())
    must(m.get("schema")=="P3_NOFW_PENNING_ON_TEMPORAL_CORRELATIONS_V1"
         and sha(CACHE)==m.get("auxiliary_cache_sha256"),
         "Temporal cache checksum invalid")
    with np.load(CACHE,allow_pickle=False) as z:
        return {k:z[k] for k in z.files},m

if __name__=="__main__":
    mode=sys.argv[1] if len(sys.argv)>1 else "selftest"
    if mode=="selftest":selftest()
    elif mode=="prepare":prepare()
    else:raise SystemExit("Usage: python SCRIPT selftest | prepare")
