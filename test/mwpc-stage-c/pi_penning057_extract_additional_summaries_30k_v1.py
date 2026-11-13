#!/usr/bin/env python3
"""Read-only 30k Penning-ON physics-summary extraction via SSHFS.

Produces per-avalanche longitudinal ion BIRTH position summaries, authentic
Geant4 3D crossing geometry, ALL TrackHeed cluster multiplicities (including
zero-accepted clusters), original per-avalanche/per-muon WALL times.

NO new Garfield calculation. NO raw archive changes. Resumable one 20-muon
chunk at a time; atomic local compressed checkpoints. Only standard-library
Python; original SHA-256 seals are verified *while streaming each CSV once*.

Run:
    python this.py extract --root ~/na6p/cluster_mwpc --out ... --channel All
    python this.py status  --out ...
    python this.py combine --out ...
The resulting *.csv.gz under OUT/summary_csv_gz are shareable summaries.

Warning: 128 uniformly sampled ion-birth positions are archived per avalanche
at most. For N_ions > 128, longitudinal mean/sigma/quantiles are SAMPLE
estimates, NOT exact full ion birth width. Times are measured elapsed WALL
seconds (std::chrono::steady_clock), not process CPU seconds. No matched
pre-multiplication arrival state was recorded by this campaign.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import math
import os
import re
import shutil
import sys
import time
from collections import defaultdict
from pathlib import Path
from zipfile import ZipFile, ZIP_DEFLATED

from pi_penning057_original2m_cluster10k_v1 import (
    OUTPUTS, SCHEMA as PHYSICS_SCHEMA, check_input, crossing_id, digest,
    metadata, require,
)


def _success(d):
    return (d.is_dir() and (d/"exit_status.txt").is_file()
            and (d/"exit_status.txt").read_text().strip()=="0"
            and (d/"verified_report.json").is_file())


def chunk_dir(campaign,index):
    token=f"{index:04d}"
    original=campaign/"results"/f"chunk_{token}"
    repaired=campaign/"recovery_runtime_v1"/f"chunk_{token}"
    oo,rr=_success(original),_success(repaired)
    require(not (oo and rr),"Duplicate completed original and recovery")
    if repaired.exists():
        require(rr and not oo,"Unfinished or invalid recovery chunk")
        return repaired,True
    require(oo,"Incomplete/missing native chunk without recovery")
    return original,False


SCHEMA="NA60_PENNING057_ORIGINAL2M_30K_ADDITIONAL_SUMMARIES_V1"
CHANNELS=("Jpsi","Phi","Omega")
NCHUNKS=500
OUTPUT_TABLES={
 "avalanche_longitudinal.csv":(
    "channel","parent_index","chunk_index","crossing_id","primary_index",
    "cluster_index","w0_mm","n_ion_births_exact","n_collected_electrons_exact",
    "n_retained_births","sampling_fraction","longitudinal_statistics_exact",
    "w_birth_mean_mm","w_birth_sigma_population_mm",
    "w_birth_p05_mm","w_birth_p95_mm"),
 "crossings_3d.csv":(
    "channel","parent_index","crossing_id","source","chunk","event",
    "hit_index","track_id","station","det_id","pdg","u_in_mm","v_in_mm",
    "w_in_mm","u_out_mm","v_out_mm","w_out_mm",
    "theta_u_in_deg","theta_w_in_deg","path_length_mm","edep_eV",
    "momentum_GeV","pu_in_GeV","pv_in_GeV","pw_in_GeV",
    "n_Heed_clusters","n_Heed_primaries","n_zero_accepted_clusters"),
 "heed_clusters.csv":(
    "channel","parent_index","crossing_id","cluster_index",
    "u_rel_cm","v_abs_cm","w_rel_cm",
    "u_global_mm","v_global_mm","w_global_mm",
    "t_ns","n_original","n_accepted"),
 "avalanche_wall_times.csv":(
    "channel","parent_index","crossing_id","primary_index",
    "cluster_index","wall_seconds","garfield_ne","garfield_ni",
    "n_ion_births_exact","n_collected_electrons_exact",
    "success","capped"),
 "muon_wall_times.csv":(
    "channel","parent_index","crossing_id","station",
    "n_Heed_clusters","n_Heed_primaries",
    "heed_wall_seconds","garfield_and_output_wall_seconds",
    "muon_total_wall_seconds"),
}

REQUIRED_NATIVE={
 "penning_crossing.csv":(
    "source|chunk|event|hit_index|track_id|station|det_id","heed_clusters",
    "heed_primaries","ion_births","collected_avalanche_electrons"),
 "penning_seed.csv":(
    "crossing_id","primary_index","seed_w_mm","ion_births","collected_e",
    "garfield_ne","garfield_ni","successful","capped"),
 "penning_ion_birth_sample.csv":(
    "crossing_id","primary_index","w_mm","seed_ion_births_total",
    "reservoir_weight"),
 "penning_heed_primaries.csv":(
    "crossing_id","primary_index","cluster_index"),
 "penning_heed_clusters.csv":(
    "crossing_id","cluster_index","u_rel_cm","v_abs_cm","w_rel_cm",
    "t_ns","n_original","n_accepted"),
 "penning_timing.csv":(
    "crossing_id","primary_index","wall_seconds",
    "avalanche_electrons","avalanche_ions","success","capped"),
 "penning_muon_timing.csv":(
    "crossing_id","station","n_Heed_clusters","n_Heed_primaries",
    "heed_wall_seconds","garfield_and_output_wall_seconds",
    "muon_total_wall_seconds"),
}
SOURCE_COLUMNS=(
 "source","chunk","event","hit_index","track_id","station","det_id",
 "pdg","u_in_mm","v_in_mm","w_in_mm","u_out_mm","v_out_mm","w_out_mm",
 "theta_u_in_deg","theta_w_in_deg","path_length_mm","edep_eV",
 "momentum_GeV","pu_in_GeV","pv_in_GeV","pw_in_GeV")


def _signed_entries(folder):
    seal=folder/"outputs.sha256"
    require(seal.is_file(),"Missing original SHA seal: "+str(seal))
    hashes={}
    for line in seal.read_text().splitlines():
        require(len(line)>=67 and line[64:66]=="  " and
                re.fullmatch(r"[a-f0-9]{64}",line[:64]) is not None,
                "Malformed original SHA line: "+str(seal))
        name=line[66:]
        require(name not in hashes and name in set(OUTPUTS)|{"verified_report.json"},
                "Unknown/duplicate original sealed file: "+name)
        hashes[name]=line[:64]
    require(set(hashes)==set(OUTPUTS)|{"verified_report.json"},
            "Original seal does not contain all native products")
    return hashes


def _checked_rows(file,sha,columns):
    """Verify original SHA in the SAME single pass as decoding huge native CSVs."""
    digest_obj=hashlib.sha256()
    with Path(file).open("rb") as stream:
        def decoded_lines():
            for line in stream:
                digest_obj.update(line)
                yield line.decode("utf-8")
        reader=csv.DictReader(decoded_lines())
        require(reader.fieldnames is not None and
                len(reader.fieldnames)==len(set(reader.fieldnames)) and
                set(columns).issubset(reader.fieldnames),
                "Original native CSV schema has changed: "+str(file))
        for row in reader:
            yield row
    require(digest_obj.hexdigest()==sha,
            "SHA mismatch in frozen native Garfield file: "+str(file))


def _native(folder,hashes,name):
    return _checked_rows(folder/name,hashes[name],REQUIRED_NATIVE[name])


def _quantile(sorted_values,q):
    if not sorted_values:
        return math.nan
    p=(len(sorted_values)-1)*q
    lo=math.floor(p)
    hi=math.ceil(p)
    return sorted_values[lo]*(hi-p)+sorted_values[hi]*(p-lo) if hi!=lo \
           else sorted_values[lo]


def _csv_bytes(name,rows):
    stream=io.StringIO(newline="")
    writer=csv.DictWriter(stream,fieldnames=OUTPUT_TABLES[name],
                          extrasaction="ignore",lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return stream.getvalue().encode("utf-8")


def _rows_for_chunk(campaign,channel,index):
    meta,source=check_input(campaign,index)
    require(meta["channel"]==channel and meta["schema"]==PHYSICS_SCHEMA
            and meta["n_chunks"]==NCHUNKS,
            "Not a frozen original-2M Penning-on 10k source")
    folder,recovered=chunk_dir(campaign,index)
    signed=_signed_entries(folder)
    require(digest(folder/"verified_report.json")==signed["verified_report.json"],
            "Original verified physics report SHA has changed")
    rep=json.loads((folder/"verified_report.json").read_text())
    require(rep["schema"]==PHYSICS_SCHEMA and
            rep["chunk_index"]==index and rep["n_muons"]==20 and
            rep["input_sha256"]==meta["chunks"][f"{index:04d}"]["sha256"] and
            rep["binary_sha256"]==meta["binary_sha256"] and
            rep["n_capped_avalanches"]==0,
            "Invalid signed input/native binary/physics provenance")

    ids=[crossing_id(row) for row in source]
    parent={cid:index*20+j for j,cid in enumerate(ids)}
    require(len(ids)==len(set(ids))==20,"Not 20 independent Geant4 crossings")
    crossings={}
    for row in _native(folder,signed,"penning_crossing.csv"):
        cid=row["source|chunk|event|hit_index|track_id|station|det_id"]
        require(cid in parent and cid not in crossings,
                "Unmatched/duplicate Garfield crossing")
        crossings[cid]=row
    require(len(crossings)==20,"Native run has missing Geant4 parents")

    clusters=defaultdict(list)
    by_cluster={}
    for row in _native(folder,signed,"penning_heed_clusters.csv"):
        cid=row["crossing_id"]
        require(cid in parent,"Cluster with foreign crossing identity")
        j=int(row["cluster_index"])
        key=(cid,j)
        require(key not in by_cluster,"Duplicate original TrackHeed cluster")
        n_original=int(row["n_original"])
        n_accepted=int(row["n_accepted"])
        require(n_original>=n_accepted>=0,
                "Physically invalid TrackHeed cluster census")
        by_cluster[key]=row
        clusters[cid].append(row)

    primaries={}
    for row in _native(folder,signed,"penning_heed_primaries.csv"):
        cid=row["crossing_id"]
        j=int(row["primary_index"])
        key=(cid,j)
        require(cid in parent and key not in primaries,
                "Duplicate/unmatched original TrackHeed electron")
        cl=int(row["cluster_index"])
        require((cid,cl) in by_cluster,
                "Heed primary belongs to nonexistent TrackHeed cluster")
        primaries[key]=cl

    seeds={}
    for row in _native(folder,signed,"penning_seed.csv"):
        cid=row["crossing_id"]
        primary_index=int(row["primary_index"])
        key=(cid,primary_index)
        require(key in primaries and key not in seeds,
                "Duplicated/missing actual Garfield avalanche")
        require(int(row["capped"])==0,
                "Capped microscopic avalanche in frozen original campaign")
        require(int(row["ion_births"])>=0
                and int(row["collected_e"])>=0,
                "Negative exact avalanche census")
        seeds[key]=row
    require(len(seeds)==len(primaries)==rep["n_true_heed_primaries"],
            "True retained Heed electron/avalanche count mismatch")

    heeds_by_muon=defaultdict(int)
    for cid,j in seeds:
        heeds_by_muon[cid]+=1
    for cid in ids:
        cc=clusters[cid]
        require(sorted(int(row["cluster_index"]) for row in cc)==
                list(range(len(cc))),
                "Original zero-accepted clusters dropped or reordered")
        require(int(crossings[cid]["heed_clusters"])==len(cc)
                and int(crossings[cid]["heed_primaries"])==
                heeds_by_muon[cid]==sum(int(r["n_accepted"]) for r in cc),
                "Geant4 / TrackHeed / Garfield per-parent mismatch")
    require(sum(int(s["ion_births"]) for s in seeds.values())==
            rep["ion_births"],
            "Original Garfield ion birth total mismatch")

    # Work on only ONE 20-muon chunk at a time, not millions of full-archive rows.
    # There are at most 128 stored birth coordinates per native avalanche.
    observed_w=defaultdict(list)
    seen_births=0
    for row in _native(folder,signed,"penning_ion_birth_sample.csv"):
        key=(row["crossing_id"],int(row["primary_index"]))
        require(key in seeds,"Ion birth row without real Garfield seed")
        n=int(seeds[key]["ion_births"])
        require(n>0 and int(row["seed_ion_births_total"])==n,
                "Mismatched original ion reservoir multiplicity")
        w=float(row["w_mm"])
        require(math.isfinite(w),"Non-finite longitudinal ion birth position")
        observed_w[key].append(w)
        require(len(observed_w[key])<=min(128,n),
                "Unexpected reservoir size; not original 128-point sampling")
        seen_births+=1
    for key,seed in seeds.items():
        require(len(observed_w[key])==min(128,int(seed["ion_births"])),
                "Missing sampled longitudinal birth coordinates")

    timings={}
    not_successful=0
    for row in _native(folder,signed,"penning_timing.csv"):
        key=(row["crossing_id"],int(row["primary_index"]))
        require(key in seeds and key not in timings,
                f"Duplicate/missing per-avalanche WALL timing: {channel} "
                f"chunk={index:04d} seed={key}")
        elapsed=float(row["wall_seconds"])
        success,capped=int(row["success"]),int(row["capped"])
        # 'success' is the return value of Garfield's AvalancheElectron().
        # The ORIGINAL 30k campaign validator required zero capped seeds,
        # but did NOT require that this return value be 1 for every seed.
        # Retain these flagged seed records rather than silently censoring
        # their physical populations or discarding their wall timings.
        require(success in (0,1) and capped==0 and
                math.isfinite(elapsed) and elapsed>=0.,
                f"Invalid native avalanche timing: channel={channel} "
                f"chunk={index:04d} crossing={key[0]} "
                f"primary={key[1]} success={success} capped={capped} "
                f"wall_seconds={row['wall_seconds']!r}")
        seed=seeds[key]
        require(int(row["avalanche_electrons"])==int(seed["garfield_ne"])
                and int(row["avalanche_ions"])==int(seed["garfield_ni"])
                and int(row["ion_births"])==int(seed["ion_births"])
                and int(row["collected_e"])==int(seed["collected_e"])
                and success==int(seed["successful"]),
                f"Original Garfield physics/timing census mismatch: "
                f"{channel} chunk={index:04d} primary={key[1]}")
        not_successful+=int(success==0)
        timings[key]=row
    require(len(timings)==len(seeds),"Native per-seed timing sample missing")

    muon_times={}
    for row in _native(folder,signed,"penning_muon_timing.csv"):
        cid=row["crossing_id"]
        require(cid in parent and cid not in muon_times,
                "Missing/duplicate full Geant4 parent timing")
        for field in ("heed_wall_seconds","garfield_and_output_wall_seconds",
                      "muon_total_wall_seconds"):
            require(math.isfinite(float(row[field])) and float(row[field])>=0.,
                    "Invalid native steady-clock wall time")
        muon_times[cid]=row
    require(len(muon_times)==20,"Not all 20 muon WALL timers recorded")

    out={name:[] for name in OUTPUT_TABLES}
    for key,seed in seeds.items():
        cid,j=key
        w=sorted(observed_w[key])
        n=int(seed["ion_births"])
        w0=float(seed["seed_w_mm"])
        require(math.isfinite(w0),"Non-finite Heed starting w")
        if w:
            mean=sum(w)/len(w)
            sigma=math.sqrt(sum((v-mean)**2 for v in w)/len(w))
            p05,p95=_quantile(w,.05),_quantile(w,.95)
        else:
            mean=sigma=p05=p95=""
        t=timings[key]
        require(int(t["ion_births"])==n and
                int(t["collected_e"])==int(seed["collected_e"]),
                "Native seed and independent timer physics counters differ")
        out["avalanche_longitudinal.csv"].append({
            "channel":channel,"parent_index":parent[cid],
            "chunk_index":index,"crossing_id":cid,"primary_index":j,
            "cluster_index":primaries[key],"w0_mm":w0,
            "n_ion_births_exact":n,
            "n_collected_electrons_exact":int(seed["collected_e"]),
            "n_retained_births":len(w),
            "sampling_fraction":(len(w)/n if n else ""),
            "longitudinal_statistics_exact":int(n<=128),
            "w_birth_mean_mm":mean,
            "w_birth_sigma_population_mm":sigma,
            "w_birth_p05_mm":p05,"w_birth_p95_mm":p95,
        })
        out["avalanche_wall_times.csv"].append({
            "channel":channel,"parent_index":parent[cid],
            "crossing_id":cid,"primary_index":j,
            "cluster_index":primaries[key],
            "wall_seconds":t["wall_seconds"],
            "garfield_ne":seed["garfield_ne"],
            "garfield_ni":seed["garfield_ni"],
            "n_ion_births_exact":n,
            "n_collected_electrons_exact":int(seed["collected_e"]),
            "success":t["success"],"capped":t["capped"],
        })

    for j,src in enumerate(source):
        cid=ids[j]
        nzero=sum(int(row["n_accepted"])==0 for row in clusters[cid])
        out["crossings_3d.csv"].append({
            "channel":channel,"parent_index":parent[cid],
            "crossing_id":cid,
            **{field:src[field] for field in SOURCE_COLUMNS},
            "n_Heed_clusters":len(clusters[cid]),
            "n_Heed_primaries":heeds_by_muon[cid],
            "n_zero_accepted_clusters":nzero,
        })
        tt=muon_times[cid]
        require(int(tt["n_Heed_clusters"])==len(clusters[cid])
                and int(tt["n_Heed_primaries"])==heeds_by_muon[cid],
                "Original muon WALL timer differs from cluster/seed census")
        out["muon_wall_times.csv"].append({
            "channel":channel,"parent_index":parent[cid],
            "crossing_id":cid,
            "station":tt["station"],
            "n_Heed_clusters":len(clusters[cid]),
            "n_Heed_primaries":heeds_by_muon[cid],
            "heed_wall_seconds":tt["heed_wall_seconds"],
            "garfield_and_output_wall_seconds":
                tt["garfield_and_output_wall_seconds"],
            "muon_total_wall_seconds":tt["muon_total_wall_seconds"],
        })
        for row in clusters[cid]:
            out["heed_clusters.csv"].append({
                "channel":channel,"parent_index":parent[cid],
                "crossing_id":cid,"cluster_index":row["cluster_index"],
                **{fld:row[fld] for fld in
                   ("u_rel_cm","v_abs_cm","w_rel_cm","t_ns",
                    "n_original","n_accepted")},
                "u_global_mm":float(src["u_in_mm"])+
                              10.*float(row["u_rel_cm"]),
                "v_global_mm":10.*float(row["v_abs_cm"]),
                "w_global_mm":float(src["w_in_mm"])+
                              10.*float(row["w_rel_cm"]),
            })
    totals={name:len(vals) for name,vals in out.items()}
    require(totals["avalanche_longitudinal.csv"]==len(seeds)
            and totals["avalanche_wall_times.csv"]==len(seeds)
            and totals["crossings_3d.csv"]==20
            and totals["muon_wall_times.csv"]==20,
            "Incomplete output after frozen native join")
    manifest={
        "schema":SCHEMA,"channel":channel,"chunk_index":index,
        "physics_schema":PHYSICS_SCHEMA,"recovered":recovered,
        "frozen_input_sha256":meta["chunks"][f"{index:04d}"]["sha256"],
        "physics_binary_sha256":meta["binary_sha256"],
        "original_outputs_seal_sha256":digest(folder/"outputs.sha256"),
        "counts":totals,
        "n_sampled_ion_birth_points":seen_births,
        "n_ion_births_exact":sum(int(s["ion_births"]) for s in seeds.values()),
        "n_garfield_avalancheelectron_false":not_successful,
        "interpretation":{
            "w0":"Exact original TrackHeed electron start w, mm, not arrival w",
            "w_statistics":"Uniform reservoir sample, max 128 per avalanche",
            "w_sigma":"Population standard deviation of retained sample",
            "wall_time":"std::chrono::steady_clock elapsed wall seconds; not CPU time",
            "arrival_surface":"Not recorded; requires a new calibration",
        }
    }
    return out,manifest


def _chunk_path(out,channel,index):
    return Path(out)/"chunk_checkpoints"/channel.lower()/f"chunk_{index:04d}.zip"


def _existing_valid(file,channel,index):
    try:
        with ZipFile(file) as zf:
            info=json.loads(zf.read("manifest.json"))
            require(info["schema"]==SCHEMA and
                    info["channel"]==channel and
                    info["chunk_index"]==index and
                    set(zf.namelist())==set(OUTPUT_TABLES)|{"manifest.json"} and
                    zf.testzip() is None,
                    "Incomplete or corrupted local chunk checkpoint")
        return True
    except Exception as exc:
        raise RuntimeError(f"Existing local checkpoint is invalid, "
                           f"refusing to overwrite: {file}. {exc}") from exc


def _save_chunk(file,channel,index,tables,manifest):
    file.parent.mkdir(parents=True,exist_ok=True)
    require(not file.exists(),"Refusing to overwrite local summary checkpoint")
    temporary=file.with_name(file.name+f".{os.getpid()}.partial")
    require(not temporary.exists(),"Stale temporary checkpoint exists: "+str(temporary))
    try:
        with ZipFile(temporary,"w",compression=ZIP_DEFLATED,
                     compresslevel=5,allowZip64=True) as archive:
            for name,rows in tables.items():
                archive.writestr(name,_csv_bytes(name,rows))
            archive.writestr("manifest.json",
                             json.dumps(manifest,indent=2)+"\n")
        with ZipFile(temporary) as archive:
            require(archive.testzip() is None,
                    "Output checkpoint failed ZIP CRC test")
        require(not file.exists(),
                "Another process already created this checkpoint")
        os.replace(temporary,file)
    finally:
        if temporary.exists():
            temporary.unlink()


def extract(root,out,channels,max_chunks,report_every):
    out=Path(out)
    require(out.resolve()!=Path(root).resolve(),
            "Output root MUST NOT be the read-only original archive")
    start=time.monotonic()
    for channel in channels:
        campaign=Path(root)/f"original2m_{channel.lower()}_10000_v1"
        m=metadata(campaign)
        require(m["channel"]==channel and m["n_chunks"]==500,
                "Invalid 10k physics channel")
        completed=0
        for index in range(max_chunks):
            path=_chunk_path(out,channel,index)
            if path.exists():
                _existing_valid(path,channel,index)
            else:
                tables,manifest=_rows_for_chunk(campaign,channel,index)
                _save_chunk(path,channel,index,tables,manifest)
            completed+=1
            if completed%report_every==0 or completed==max_chunks:
                print(f"{channel}: {completed}/{max_chunks} read-only "
                      f"verified 20-muon chunks · "
                      f"{time.monotonic()-start:.1f} s elapsed",flush=True)
    print("PASS extracted. These are LOCAL checkpoints; run combine after "
          "500/500 chunks per channel.",flush=True)


def status(out):
    for channel in CHANNELS:
        indices=[i for i in range(NCHUNKS)
                 if _chunk_path(out,channel,i).is_file()]
        absent=[i for i in range(NCHUNKS) if i not in set(indices)]
        print(f"{channel}: {len(indices)}/{NCHUNKS} chunk checkpoints; "
              f"first missing: {absent[0] if absent else 'none'}")
    folder=Path(out)/"summary_csv_gz"
    if folder.is_dir():
        for path in sorted(folder.glob("*.csv.gz")):
            print(f"OUTPUT {path.name}: {path.stat().st_size/1024/1024:.1f} MiB")


def combine(out):
    out=Path(out)
    checkpoints=[]
    for channel in CHANNELS:
        for index in range(NCHUNKS):
            p=_chunk_path(out,channel,index)
            require(p.is_file(),f"Missing checkpoint: {p} "
                    "(rerun extract; it is resumable)")
            _existing_valid(p,channel,index)
            checkpoints.append(p)
    folder=out/"summary_csv_gz"
    folder.mkdir(parents=True,exist_ok=True)
    require(not any((folder/(name+".gz")).exists() for name in OUTPUT_TABLES),
            "Final summary already exists. Refusing to overwrite outputs.")
    overall={k:0 for k in OUTPUT_TABLES}
    births=sampled=recovered=garfield_false=0
    legacy_status_manifests=0
    for k,name in enumerate(OUTPUT_TABLES):
        destination=folder/(name+".gz")
        temp=destination.with_name(destination.name+f".{os.getpid()}.partial")
        require(not temp.exists(),"Temporary combined output already exists")
        try:
            with gzip.open(temp,"wb",compresslevel=5) as target:
                header=None
                for file in checkpoints:
                    with ZipFile(file) as zf:
                        with zf.open(name) as stream:
                            first=stream.readline()
                            if header is None:
                                header=first
                                target.write(header)
                            else:
                                require(first==header,
                                        f"CSV schema inconsistent: {file} {name}")
                            shutil.copyfileobj(stream,target,1024*1024)
                        if k==0:
                            m=json.loads(zf.read("manifest.json"))
                            for table,num in m["counts"].items():
                                overall[table]+=num
                            births+=m["n_ion_births_exact"]
                            sampled+=m["n_sampled_ion_birth_points"]
                            recovered+=int(m["recovered"])
                            if "n_garfield_avalancheelectron_false" in m:
                                garfield_false+=int(
                                    m["n_garfield_avalancheelectron_false"])
                            else:
                                # Early 30k checkpoints predate the status
                                # counter, but still contain the original
                                # per-avalanche success column. Recover it
                                # exactly from the LOCAL checkpoint; do not
                                # re-read SSHFS and do not treat missing
                                # metadata as zero unsuccessful avalanches.
                                legacy_status_manifests+=1
                                with zf.open("avalanche_wall_times.csv") as raw:
                                    with io.TextIOWrapper(
                                            raw,encoding="utf-8",
                                            newline="") as decoded:
                                        rows=csv.DictReader(decoded)
                                        require(rows.fieldnames is not None
                                                and "success" in rows.fieldnames,
                                                "Old checkpoint lacks saved "
                                                f"Garfield success flags: {file}")
                                        n_rows=n_false=0
                                        for row in rows:
                                            flag=int(row["success"])
                                            require(flag in (0,1),
                                                    "Non-binary Garfield "
                                                    f"success in {file}")
                                            n_false+=int(flag==0)
                                            n_rows+=1
                                require(n_rows==m["counts"][
                                            "avalanche_wall_times.csv"],
                                        "Checkpoint manifest count differs "
                                        f"from timing rows: {file}")
                                garfield_false+=n_false
            require(not destination.exists(),"Refusing overwrite while combining")
            os.replace(temp,destination)
        finally:
            if temp.exists():
                temp.unlink()
        print(f"COMBINED {name} · {destination.stat().st_size/1024/1024:.1f} MiB",
              flush=True)
    print(f"Recovered original Garfield success flags from "
          f"{legacy_status_manifests:,} earlier local checkpoints.",
          flush=True)
    require(overall["crossings_3d.csv"]==30000
            and overall["muon_wall_times.csv"]==30000
            and overall["avalanche_longitudinal.csv"]==
                overall["avalanche_wall_times.csv"]==1554258
            and births==5072644729 and recovered==3,
            "Combined census differs from frozen 30k physics campaign")
    manifest={
        "schema":SCHEMA,
        "channels":CHANNELS,
        "total_original_real_Geant4_muons":30000,
        "total_microscopic_Garfield_avalanches":1554258,
        "total_exact_positive_ion_births":births,
        "sampled_ion_birth_positions":sampled,
        "recovered_20_muon_chunks":recovered,
        "n_garfield_avalancheelectron_false":garfield_false,
        "tables":overall,
        "limitations":[
            "128 uniformly sampled spatial births at most per avalanche",
            "For >128 births, w_mean/sigma/p05/p95 are sample-derived, not exact",
            "For zero-ion avalanches longitudinal statistics are undefined",
            "Original recorded times are wall-clock elapsed, not CPU times",
            "Fast model speedup requires matched separate fast-model benchmark",
            "Matched electron pre-multiplication arrival states were not recorded",
        ],
        "sha256":{
            name+".gz":digest(folder/(name+".gz")) for name in OUTPUT_TABLES
        }
    }
    (folder/"provenance.json").write_text(json.dumps(manifest,indent=2)+"\n")
    print("PASS 30,000 real parent muons; 1,554,258 exact native avalanches; "
          f"{births:,} ion births; {recovered} repaired chunks; "
          f"{garfield_false:,} Garfield AvalancheElectron false returns; "
          "five independent compact tables + provenance",flush=True)
    print("SUMMARY DIRECTORY:",folder.resolve(),flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command",choices=("extract","status","combine"))
    parser.add_argument("--root",type=Path,
                        default=Path.home()/"na6p/cluster_mwpc")
    parser.add_argument("--out",type=Path,required=True)
    parser.add_argument("--channel",choices=(*CHANNELS,"All"),default="All")
    parser.add_argument("--max-chunks",type=int,default=500,
                        help="Pilot: --max-chunks 1; full: 500")
    parser.add_argument("--report-every",type=int,default=1)
    options=parser.parse_args()
    require(1<=options.max_chunks<=500 and options.report_every>0,
            "Invalid chunk count/progress interval")
    channels=CHANNELS if options.channel=="All" else (options.channel,)
    if options.command=="extract":
        extract(options.root,options.out,channels,options.max_chunks,
                options.report_every)
    elif options.command=="status":
        status(options.out)
    else:
        combine(options.out)


if __name__=="__main__":
    try:
        main()
    except (RuntimeError,ValueError,OSError,KeyError,AssertionError) as error:
        print("FAILED additional-summary extraction: "+str(error),
              file=sys.stderr,flush=True)
        raise SystemExit(1)
