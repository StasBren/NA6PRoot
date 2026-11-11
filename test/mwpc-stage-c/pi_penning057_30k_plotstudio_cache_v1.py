#!/usr/bin/env python3
"""Prepare reusable 10k/channel + pooled 30k Plot Studio caches (read-only input).

This is the same 3mm, no-field-wire, 1800V, Penning r=0.57 TrackHeed and
full Garfield microscopic model as the legacy 1000 J/psi Plot Studio.
Preparation uses only frozen original-2M Geant4 crossings and their SHA-sealed
native outputs. It never runs or modifies physics, and never overwrites a cache.

Run one PBS preparation per channel, then a single combine task:
  python this.py prepare --channel Jpsi
  python this.py prepare --channel Phi
  python this.py prepare --channel Omega
  python this.py combine
Produced: /storage/cgrp/stanis/mwpc_penning057/plotstudio_30k_v1/<scope>/
where scope is jpsi/phi/omega/all and contains dashboard_data.npz,
provenance.json, temporal_correlations_v1.npz, temporal_correlations_v1.json.

Exact integer per-seed gains and ion sector counts are distinguished from
sampled ion-birth space/time (<128 points per avalanche) weighted by true ion
multiplicity. Per-muon gain is NaN if its Heed-primary count is zero.
The pooled view is a descriptive equal-10k-muon-per-channel mixture;
per-avalanche pooling weights channels by their actual Heed multiplicities.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path
import numpy as np

from pi_penning057_original2m_cluster10k_v1 import (
    OUTPUTS, SCHEMA, check_input, digest, metadata, require, readcsv,
    crossing_id,
)
from pi_penning057_original2m_full_recovery_audit_v1 import sealed_files
from pi_p3_noFW_penning057_interactive_dashboard_v1 import (
    PHI_EDGES, U_EDGES, V_EDGES, T_EDGES, T_PHI_EDGES,
    nearest_wire_offset, phi_idx, binned, ahemi,
)

DEFAULT_BASE = Path("/storage/cgrp/stanis/mwpc_penning057")
CHANNELS = ("Jpsi", "Phi", "Omega")
CROSS_ID = "source|chunk|event|hit_index|track_id|station|det_id"
DASHBOARD_SCHEMA = "P3_NOFW_PENNING_ON_30K_PLOTSTUDIO_V1"
TEMPORAL_SCHEMA = "P3_NOFW_PENNING_ON_TEMPORAL_CORRELATIONS_V1"
SEED_FIELDS = (
    "parent","primary","cluster","gain","ions","A","u0","v0",
    "t_heed","t10","t50","t90","span","sample_count",
)
INT_SEED_FIELDS = {"parent","primary","cluster","gain","ions","sample_count"}
MAIN_SEED_FIELDS = (
    "gain","raw_ne","births","seed_u","seed_v","seed_A",
    "seed_absA","seed_span_ns","seed_mean_birth_dt_ns",
)
MAIN_MUON_FIELDS = ("muon_gain","muon_heed","muon_births","muon_A","muon_absA")
SUM_KEYS = ("angle_equal","angle_ion","xy_ion","start_sum","start_n",
            "t_weights","t_phi")
EDGES = ("phi_edges","u_edges","v_edges","t_edges","t_phi_edges")


def rows(path, required):
    names, values = readcsv(path)
    require(set(required).issubset(names),
            f"Missing required native columns: {path}")
    return values


def stream_rows(path, required):
    """Never materialise a potentially multi-GB ion birth CSV in RAM."""
    with Path(path).open(newline="") as file:
        reader=csv.DictReader(file)
        require(reader.fieldnames is not None and
                len(reader.fieldnames)==len(set(reader.fieldnames)) and
                set(required).issubset(reader.fieldnames),
                f"Missing/duplicated native CSV columns: {path}")
        yield from reader


def successful(d):
    return (d.is_dir() and (d/"exit_status.txt").is_file()
            and (d/"exit_status.txt").read_text().strip()=="0"
            and (d/"verified_report.json").is_file())


def chunk_dir(camp,index):
    four=f"{index:04d}"
    original=camp/"results"/f"chunk_{four}"
    recovered=camp/"recovery_runtime_v1"/f"chunk_{four}"
    oo,rr=successful(original),successful(recovered)
    require(not (oo and rr),f"Duplicate successful original/recovery: {camp}:{four}")
    if recovered.exists():
        require(rr and not oo, f"Unfinished or invalid recovery: {camp}:{four}")
        return recovered,True
    require(oo, f"Missing or failed chunk without recovery: {camp}:{four}")
    return original,False


def write_cache(path,arrays):
    require(not path.exists(),f"Refusing to overwrite cache: {path}")
    np.savez_compressed(path,**arrays)
    require(path.is_file() and path.stat().st_size>1000,
            f"Cache write did not complete: {path}")


def prepare(base, channel, report_every, output_root, skip_binary_check):
    camp=base/f"original2m_{channel.lower()}_10000_v1"
    meta=metadata(camp)
    require(meta["channel"]==channel and meta["n_chunks"]==500
            and meta["schema"]==SCHEMA,"Wrong channel/physics campaign")
    frozen=camp/"frozen_original2m_10000_input.csv"
    require(digest(frozen)==meta["frozen_source_sha256"],
            "Original 10k Geant4 input changed")
    require(digest(camp/"selection_manifest.json")==
            meta["original_2m_selection_manifest_sha256"],
            "Original 2M selection manifest changed")
    if not skip_binary_check:
        require(digest(Path(meta["binary_path"]))==meta["binary_sha256"],
                "Native Penning-on executable changed")
    else:
        # Off-cluster read-only analysis: absolute binary path exists only on
        # remote host. Every SHA-sealed chunk report still must carry the
        # exact frozen binary_sha256 from the physics campaign manifest.
        print("NOTE: skipping local executable-file SHA; checking frozen "
              "manifest binary SHA against EVERY signed per-chunk report.",
              flush=True)

    out=output_root/channel.lower()
    require(not out.exists(), f"Refusing overwrite existing plot cache: {out}")

    m={k:[] for k in MAIN_SEED_FIELDS+MAIN_MUON_FIELDS}
    t={k:[] for k in SEED_FIELDS}
    angular_chunks=[]
    angle_equal=np.zeros(72,dtype=np.float64)
    angle_ion=np.zeros(72,dtype=np.float64)
    xy_ion=np.zeros((len(V_EDGES)-1,len(U_EDGES)-1),dtype=np.float64)
    start_sum=np.zeros_like(xy_ion)
    start_n=np.zeros_like(xy_ion)
    time_ion=np.zeros(len(T_EDGES)-1,dtype=np.float64)
    time_phi=np.zeros((len(T_PHI_EDGES)-1,72),dtype=np.float64)
    muons=[]
    collected=births=primaries=sampled=ion_avalanches=0
    zero_collected=zero_heed=recoveries=0

    for block in range(500):
        outchunk,recovered=chunk_dir(camp,block)
        recoveries+=int(recovered)
        _,src=check_input(camp,block)
        sealed_files(outchunk)  # Full SHA of seven CSVs + verified report
        report=json.loads((outchunk/"verified_report.json").read_text())
        require(report["chunk_index"]==block and
                report["n_muons"]==20 and
                report["schema"]==SCHEMA and
                report["binary_sha256"]==meta["binary_sha256"] and
                report["input_sha256"]==meta["chunks"][f"{block:04d}"]["sha256"] and
                report["n_capped_avalanches"]==0,
                f"Invalid frozen chunk scientific provenance: {channel}:{block}")

        cr=rows(outchunk/"penning_crossing.csv",
                (CROSS_ID,"heed_primaries","heed_clusters","ion_births",
                 "collected_avalanche_electrons","ion_up","ion_down",
                 "theta_u_in_deg","theta_w_in_deg","path_mm"))
        require(len(cr)==20 and [r[CROSS_ID] for r in cr]==
                [crossing_id(r) for r in src],
                "Geant4 crossing order/identity changed")
        imap={}
        for j,r in enumerate(cr):
            cid=r[CROSS_ID]; parent=block*20+j
            require(cid not in imap,"Duplicate crossing ID")
            imap[cid]=parent
            hp=int(r["heed_primaries"])
            nc=int(r["heed_clusters"])
            col=int(r["collected_avalanche_electrons"])
            ions=int(r["ion_births"])
            up,down=int(r["ion_up"]),int(r["ion_down"])
            require(min(hp,nc,col,ions,up,down)>=0,"Negative per-muon count")
            if hp==0:
                require(col==0 and ions==0 and up==0 and down==0,
                        "Zero-Heed crossing has nonzero Garfield products")
                zero_heed+=1
            A=ahemi(up,down)
            m["muon_heed"].append(hp)
            m["muon_gain"].append(col/hp if hp else float("nan"))
            m["muon_births"].append(ions)
            m["muon_A"].append(A)
            m["muon_absA"].append(abs(A) if math.isfinite(A) else float("nan"))
            muons.append({
                "id":cid,"channel":channel,"n_Heed":hp,"n_clusters":nc,
                "theta_u_deg":float(r["theta_u_in_deg"]),
                "theta_w_deg":float(r["theta_w_in_deg"]),
                "path_mm":float(r["path_mm"]),
                "collected_e":col,"ion_births":ions,
            })

        prim=rows(outchunk/"penning_heed_primaries.csv",
                  ("crossing_id","primary_index","cluster_index","t_ns"))
        heeds={}
        for r in prim:
            key=(r["crossing_id"],int(r["primary_index"]))
            require(key not in heeds and key[0] in imap,
                    "Duplicated/unmatched Heed primary")
            heeds[key]=(int(r["cluster_index"]),float(r["t_ns"]))
        cluster_count=defaultdict(int)
        for r in rows(outchunk/"penning_heed_clusters.csv",
                      ("crossing_id","cluster_index","n_accepted")):
            cluster_count[r["crossing_id"]]+=int(r["n_accepted"])
        for r in cr:
            cid=r[CROSS_ID]
            require(sum(1 for k in heeds if k[0]==cid)==
                    cluster_count[cid]==int(r["heed_primaries"]),
                    "Heed cluster/primary discrepancy")

        seeds=rows(outchunk/"penning_seed.csv",
                   ("crossing_id","primary_index","collected_e","ion_births",
                    "garfield_ne","seed_u_mm","seed_v_mm",
                    "ion_up","ion_down","capped"))
        require(len(seeds)==len(heeds)==report["n_true_heed_primaries"],
                f"Heed/Garfield/verified report mismatch: {channel}:{block}")
        keys={}
        for j,r in enumerate(seeds):
            k=(r["crossing_id"],int(r["primary_index"]))
            require(k in heeds and k not in keys,"Invalid seed identity")
            keys[k]=j

        nseed=len(seeds)
        times=[[] for _ in range(nseed)]
        frac=np.zeros((nseed,72),dtype=np.float32)
        sample_counts=np.zeros(nseed,dtype=np.int32)
        delay_min=np.full(nseed,np.inf,dtype=float)
        delay_max=np.full(nseed,-np.inf,dtype=float)
        delay_sum=np.zeros(nseed,dtype=np.float64)
        weighted_sum=0.
        for r in stream_rows(outchunk/"penning_ion_birth_sample.csv",
                      ("crossing_id","primary_index","u_mm","v_mm","t_ns",
                       "seed_ion_births_total","reservoir_weight")):
            k=(r["crossing_id"],int(r["primary_index"]))
            require(k in keys,"Unmatched ion-birth sample")
            j=keys[k]
            b=int(seeds[j]["ion_births"])
            n=min(128,b)
            require(b>0 and int(r["seed_ion_births_total"])==b,
                    "Invalid sampled ion avalanche multiplicity")
            w=float(r["reservoir_weight"])
            require(math.isfinite(w) and abs(w-b/n)<=max(1.,b/n)*1e-7,
                    "Ion-birth sample weight inconsistent")
            birth=float(r["t_ns"])
            delay=birth-heeds[k][1]
            require(math.isfinite(delay) and delay>=-1e-6,
                    "Ion birth predates initial Heed electron")
            dt=max(0.,delay)
            u,v=float(r["u_mm"]),float(r["v_mm"])
            phi=phi_idx(u,v)
            frac[j,phi]+=1.
            angle_equal[phi]+=1./n
            angle_ion[phi]+=w
            xi=binned(nearest_wire_offset(u),U_EDGES)
            yi=binned(v,V_EDGES)
            if xi>=0 and yi>=0:xy_ion[yi,xi]+=w
            ti=binned(dt,T_EDGES)
            qi=binned(dt,T_PHI_EDGES)
            require(ti>=0 and qi>=0,
                    "Physical ion-birth delay outside original histogram edges")
            time_ion[ti]+=w
            time_phi[qi,phi]+=w
            sample_counts[j]+=1
            weighted_sum+=w
            times[j].append(birth)
            delay_min[j]=min(delay_min[j],dt)
            delay_max[j]=max(delay_max[j],dt)
            delay_sum[j]+=dt
            sampled+=1

        block_col=block_birth=block_ionsamp=0
        for j,r in enumerate(seeds):
            k=(r["crossing_id"],int(r["primary_index"]))
            b=int(r["ion_births"])
            col=int(r["collected_e"])
            ne=int(r["garfield_ne"])
            up,down=int(r["ion_up"]),int(r["ion_down"])
            require(min(b,col,ne,up,down)>=0 and int(r["capped"])==0,
                    "Capped/invalid microscopic avalanche")
            n=min(128,b)
            require(sample_counts[j]==n,"Missing reservoir samples in avalanche")
            if n:
                q=np.quantile(times[j],[.1,.5,.9])
                span=float(delay_max[j]-delay_min[j])
                mean_dt=float(delay_sum[j]/n)
                frac[j]/=n
                ion_avalanches+=1
            else:
                q=np.array([np.nan,np.nan,np.nan])
                span=mean_dt=float("nan")
            A=ahemi(up,down)
            u0=nearest_wire_offset(float(r["seed_u_mm"]))
            v0=float(r["seed_v_mm"])
            xi=binned(u0,U_EDGES)
            yi=binned(v0,V_EDGES)
            if xi>=0 and yi>=0:
                start_n[yi,xi]+=1.
                start_sum[yi,xi]+=col
            vals={
                "gain":col,"raw_ne":ne,"births":b,"seed_u":u0,
                "seed_v":v0,"seed_A":A,
                "seed_absA":abs(A) if math.isfinite(A) else float("nan"),
                "seed_span_ns":span,"seed_mean_birth_dt_ns":mean_dt,
            }
            for name in MAIN_SEED_FIELDS:m[name].append(vals[name])
            cluster,th=heeds[k]
            tv={
                "parent":imap[k[0]],"primary":k[1],"cluster":cluster,
                "gain":col,"ions":b,"A":A,"u0":u0,"v0":v0,
                "t_heed":th,"t10":q[0],"t50":q[1],"t90":q[2],
                "span":span,"sample_count":n,
            }
            for name in SEED_FIELDS:t[name].append(tv[name])
            collected+=col;births+=b;primaries+=1
            zero_collected+=int(col==0)
            block_col+=col;block_birth+=b;block_ionsamp+=n
        angular_chunks.append(frac)
        require(block_col==report["collected_electrons"] and
                block_birth==report["ion_births"] and
                abs(weighted_sum-block_birth)<1e-5*max(1,block_birth) and
                block_col==sum(int(r["collected_avalanche_electrons"]) for r in cr)
                and block_birth==sum(int(r["ion_births"]) for r in cr),
                "Raw/verified physical charge or ion census mismatch")
        if (block+1)%report_every==0:
            print(f"{channel}: verified/processed {block+1}/500 chunks "
                  f"({primaries:,} Garfield avalanches)",flush=True)

    require(len(muons)==10000 and primaries==len(t["gain"])==
            len(m["gain"]) and int(sum(m["muon_heed"]))==primaries
            and len(m["muon_heed"])==10000 and
            abs(angle_ion.sum()-births)<1e-5*max(1,births) and
            abs(time_ion.sum()-births)<1e-5*max(1,births) and
            abs(angle_equal.sum()-ion_avalanches)<1e-5*max(1,ion_avalanches),
            "Full source/ion/Heed census mismatch")

    main={}
    for name,values in m.items():
        dtype=np.float32 if name.startswith("seed_") or name.startswith("muon_") else np.int32
        if name in ("muon_heed","muon_births"):dtype=np.int32
        main[name]=np.asarray(values,dtype=dtype)
    main.update(angle_equal=angle_equal,angle_ion=angle_ion,xy_ion=xy_ion,
                start_sum=start_sum,start_n=start_n,t_weights=time_ion,
                t_phi=time_phi,phi_edges=PHI_EDGES,u_edges=U_EDGES,
                v_edges=V_EDGES,t_edges=T_EDGES,t_phi_edges=T_PHI_EDGES)
    temporal={}
    for name,values in t.items():
        temporal[name]=np.asarray(
            values,dtype=np.int32 if name in INT_SEED_FIELDS else np.float64)
    temporal["angle_fraction"]=np.concatenate(angular_chunks,axis=0)
    require(temporal["angle_fraction"].shape==(primaries,72),
            "Loss of per-avalanche angular samples")

    out.mkdir(parents=True,exist_ok=False)
    n0=int(sum(v==0 for v in m["muon_heed"]))
    dashboard_meta={
        "schema":DASHBOARD_SCHEMA,"dataset":channel.lower(),
        "dataset_label":channel,"n_muons":10000,
        "n_independent_Geant4_muons":10000,
        "n_muons_with_zero_retained_Heed_primaries":n0,
        "n_muons_with_retained_Heed_primaries":10000-n0,
        "n_microscopic_TrackHeed_primary_Garfield_avalanches":primaries,
        "n_true_ion_births":births,
        "n_ion_producing_avalanches":ion_avalanches,
        "n_saved_ion_birth_coordinate_samples":sampled,
        "n_zero_collected_avalanche_seeds":zero_collected,
        "n_capped_avalanches":0,
        "sum_collected_electrons":collected,
        "mean_collected_electrons_per_aval":collected/primaries,
        "recovery_chunks":recoveries,
        "physics":"p3_noFW 3mm pitch 1800V +/-3mm cathodes Ar/CO2 70:30 Penning r=0.57",
        "source":"10,000 authentic independent parent-event Geant4 gas crossings selected from original 2M/channel",
        "frozen_input_sha256":meta["frozen_source_sha256"],
        "frozen_binary_sha256":meta["binary_sha256"],
        "source_manifest_sha256":digest(camp/"selection_manifest.json"),
        "normalization":"Per-muon denominators exclude undefined 0/0; avalanche pooling weights actual retained Heed seeds",
        "cathode_signal":"NOT simulated; ion births and collected electron counts only",
    }
    write_cache(out/"dashboard_data.npz",main)
    dashboard_meta["cache_sha256"]=digest(out/"dashboard_data.npz")
    (out/"provenance.json").write_text(json.dumps(dashboard_meta,indent=2)+"\n")
    temporal_meta={
        "schema":TEMPORAL_SCHEMA,"dataset":channel.lower(),
        "muons":muons,"n_independent_true_Geant4_muons":10000,
        "n_exact_Garfield_avalanches":primaries,
        "n_true_ion_births":births,
        "n_sampled_ion_birth_positions":sampled,
        "recommended_sample_count_threshold":128,
        "t10_definition":"10th percentile of <=128 uniform sampled ion birth timestamps relative to original muon gas-entry time, NOT electron arrival",
        "original_dashboard_cache_sha256":dashboard_meta["cache_sha256"],
    }
    write_cache(out/"temporal_correlations_v1.npz",temporal)
    temporal_meta["auxiliary_cache_sha256"]=digest(out/"temporal_correlations_v1.npz")
    (out/"temporal_correlations_v1.json").write_text(
        json.dumps(temporal_meta,indent=2)+"\n")
    print(f"PASS {channel}: 10,000 independent Geant4 muons, "
          f"{primaries:,} full Garfield avalanches, {births:,} ion births, "
          f"{recoveries} separate recovered chunks; full SHA and charge conservation",flush=True)
    print("CACHE:",out,flush=True)


def combine(output_root):
    folder=output_root
    target=folder/"all"
    require(not target.exists(),"Refusing overwrite of pooled 30k cache")
    metas=[];tmetas=[];data=[];temps=[]
    for channel in CHANNELS:
        source=folder/channel.lower()
        pm=source/"provenance.json"
        tm=source/"temporal_correlations_v1.json"
        require(pm.is_file() and tm.is_file(),"Prepare all three 10k caches before pooling")
        m=json.loads(pm.read_text());t=json.loads(tm.read_text())
        require(m["dataset"]==channel.lower() and m["n_muons"]==10000
                and t["n_independent_true_Geant4_muons"]==10000,
                "Incorrect per-channel cache provenance")
        pf=source/"dashboard_data.npz"
        tf=source/"temporal_correlations_v1.npz"
        require(digest(pf)==m["cache_sha256"] and
                digest(tf)==t["auxiliary_cache_sha256"],
                "Input cache SHA changed prior to pooling")
        with np.load(pf,allow_pickle=False) as f:
            data.append({k:f[k] for k in f.files})
        with np.load(tf,allow_pickle=False) as f:
            temps.append({k:f[k] for k in f.files})
        metas.append(m);tmetas.append(t)
    require(all(set(d)==set(data[0]) for d in data) and
            all(set(t)==set(temps[0]) for t in temps),
            "Incompatible 10k dashboard schemas")

    merged={}
    for k in data[0]:
        if k in EDGES:
            require(all(np.array_equal(d[k],data[0][k]) for d in data[1:]),
                    "Different bin edges, cannot pool histograms")
            merged[k]=data[0][k]
        elif k in SUM_KEYS:merged[k]=sum((d[k] for d in data[1:]),data[0][k].copy())
        else:merged[k]=np.concatenate([d[k] for d in data])
    combined={}
    for k in temps[0]:
        if k=="parent":
            combined[k]=np.concatenate([t[k]+n*10000
                                        for n,t in enumerate(temps)])
        else:combined[k]=np.concatenate([t[k] for t in temps],axis=0)

    muons=[entry for t in tmetas for entry in t["muons"]]
    nseed=sum(m["n_microscopic_TrackHeed_primary_Garfield_avalanches"] for m in metas)
    nion=sum(m["n_true_ion_births"] for m in metas)
    require(len(muons)==30000 and len(merged["muon_heed"])==30000
            and len(merged["gain"])==nseed and
            len(combined["gain"])==nseed and
            int(sum(merged["muon_heed"]))==nseed and
            int(sum(merged["births"]))==nion,
            "Pooled data do not conserve primary/ion counts")
    target.mkdir(parents=True,exist_ok=False)
    write_cache(target/"dashboard_data.npz",merged)
    meta={
        "schema":DASHBOARD_SCHEMA,"dataset":"all",
        "dataset_label":"All (J/ψ + φ + ω)",
        "n_muons":30000,"n_independent_Geant4_muons":30000,
        "n_muons_with_zero_retained_Heed_primaries":
            sum(m["n_muons_with_zero_retained_Heed_primaries"] for m in metas),
        "n_microscopic_TrackHeed_primary_Garfield_avalanches":nseed,
        "n_true_ion_births":nion,
        "n_ion_producing_avalanches":sum(m["n_ion_producing_avalanches"] for m in metas),
        "n_saved_ion_birth_coordinate_samples":
            sum(m["n_saved_ion_birth_coordinate_samples"] for m in metas),
        "n_zero_collected_avalanche_seeds":
            sum(m["n_zero_collected_avalanche_seeds"] for m in metas),
        "n_capped_avalanches":0,
        "sum_collected_electrons":sum(m["sum_collected_electrons"] for m in metas),
        "mean_collected_electrons_per_aval":
            sum(m["sum_collected_electrons"] for m in metas)/nseed,
        "recovery_chunks":sum(m["recovery_chunks"] for m in metas),
        "source":"10,000 genuine original-Geant4 parents/channel, descriptive 1:1:1 mixture",
        "normalization":"Per muon: 30k independent crossing parents; per avalanche: all 1.55M actual Heed seeds (not 3 equal-weight channel means)",
        "physics":"p3_noFW 3mm pitch 1800V +/-3mm cathodes Ar/CO2 70:30 Penning r=0.57",
        "cathode_signal":"NOT simulated; ion births and collected electrons only",
        "source_cache_shas":[m["cache_sha256"] for m in metas],
        "cache_sha256":digest(target/"dashboard_data.npz"),
    }
    (target/"provenance.json").write_text(json.dumps(meta,indent=2)+"\n")
    write_cache(target/"temporal_correlations_v1.npz",combined)
    tmeta={
        "schema":TEMPORAL_SCHEMA,"dataset":"all","muons":muons,
        "n_independent_true_Geant4_muons":30000,
        "n_exact_Garfield_avalanches":nseed,
        "n_true_ion_births":nion,
        "n_sampled_ion_birth_positions":sum(t["n_sampled_ion_birth_positions"] for t in tmetas),
        "recommended_sample_count_threshold":128,
        "t10_definition":"10th percentile of <=128 sampled ion-birth timestamps after original muon entry, NOT arrival time",
        "original_dashboard_cache_sha256":meta["cache_sha256"],
        "auxiliary_cache_sha256":digest(target/"temporal_correlations_v1.npz"),
    }
    (target/"temporal_correlations_v1.json").write_text(json.dumps(tmeta,indent=2)+"\n")
    print(f"PASS ALL 30,000 muons, {nseed:,} full avalanches, {nion:,} ion births; "
          f"4 independently selectable Plot Studio datasets",flush=True)
    print("CACHE:",folder,flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("command",choices=("prepare","combine"))
    p.add_argument("--base",type=Path,default=DEFAULT_BASE,
                   help="Read-only original native Garfield campaign archive")
    p.add_argument("--output-root",type=Path,
                   help="Separate WRITEABLE output location for derived caches; "
                        "defaults to BASE/plotstudio_30k_v1 on cluster")
    p.add_argument("--skip-binary-check",action="store_true",
                   help="Off-cluster only: native executable is unavailable; "
                        "cross-check its SHA from frozen manifest in every chunk report")
    p.add_argument("--channel",choices=CHANNELS)
    p.add_argument("--report-every",type=int,default=25)
    a=p.parse_args()
    require(a.report_every>0,"Invalid reporting interval")
    output_root=a.output_root or (a.base/"plotstudio_30k_v1")
    if a.command=="prepare":
        require(a.channel is not None,"--channel required for prepare")
        prepare(a.base,a.channel,a.report_every,output_root,a.skip_binary_check)
    else:
        combine(output_root)


if __name__=="__main__":
    main()
