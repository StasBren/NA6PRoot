#!/usr/bin/env python3
"""Original-2M 3x10k native Geant4/TrackHeed/Garfield wire-truth production.

Operations: prepare, lock-binary, runtime-check, verify, audit, gate, status,
report, export. All native physics runs are performed by the PBS shell worker.
Original physics archives are never written. No implicit PBS submissions.
"""
from __future__ import annotations
import argparse
import csv
import gzip
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

from pi_penning057_original2m_cluster10k_v1 import (
    SCHEMA as OLD_SCHEMA, check_input, crossing_id, digest, metadata, require,
)

CHANNELS=("Jpsi","Phi","Omega")
ROOT_DEFAULT=Path("/storage/cgrp/stanis/mwpc_penning057_wiretruth30k_v1")
OLD_DEFAULT=Path("/storage/cgrp/stanis/mwpc_penning057")
REPO=Path(__file__).resolve().parents[2]
SOURCE=REPO/"test/mwpc-avalanche/phaseB_native_muon_p3_noFW_penning057_wiretruth_v1.cxx"
BINARY=REPO/"test/mwpc-avalanche/build-wiretruth30k-cvmfs/mwpc_p3_noFW_penning057_wiretruth_v1"
SCHEMA="NA60_WIRETRUTH30K_PENNING057_PRODUCTION_V1"
PARTS=("crossing","seed","ion_birth_sample","heed_primaries","heed_clusters",
       "timing","muon_timing","wire_seed","wire_muon")
RESERVOIR=128
PITCH_MM=3.
OUTPUT_COLUMNS={
 "wire_seed":(
    "crossing_id","station","theta_u_deg","theta_w_deg","primary_index",
    "cluster_index","wire_index_global","wire_u_global_mm","collected_e",
    "ion_births","ion_up","ion_down","ion_boundary","ion_top","ion_bottom",
    "ion_side","A_ion_hemi","A_ion_TB","mean_w_birth_global_mm",
    "sigma_w_birth_mm","mean_t_birth_ns","sigma_t_birth_ns",
    "mean_t_collected_ns","sigma_t_collected_ns"),
 "seed":("crossing_id","primary_index","garfield_seed","ion_births","collected_e",
         "garfield_ne","garfield_ni","attached","other","successful","capped"),
 "crossing":("source|chunk|event|hit_index|track_id|station|det_id",
             "heed_clusters","heed_primaries","collected_avalanche_electrons",
             "ion_births","zero_seeds","capped_seeds","failed_seeds"),
 "heed_primaries":("crossing_id","primary_index","cluster_index",
                   "u_rel_cm","v_abs_cm","w_rel_cm","t_ns"),
 "heed_clusters":("crossing_id","cluster_index","n_original","n_accepted",
                  "u_rel_cm","v_abs_cm","w_rel_cm","t_ns"),
 "timing":("crossing_id","primary_index","wall_seconds","garfield_seed",
           "ion_births","collected_e","success","capped"),
 "muon_timing":("crossing_id","n_Heed_clusters","n_Heed_primaries",
                "heed_wall_seconds","garfield_and_output_wall_seconds",
                "muon_total_wall_seconds"),
 "ion_birth_sample":("crossing_id","primary_index","u_mm","v_mm","w_mm",
                      "t_ns","seed_ion_births_total","reservoir_weight"),
}


def _json(path,obj):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    require(not path.exists(),"Refusing to overwrite immutable file: "+str(path))
    temporary=path.with_name(path.name+f".{os.getpid()}.part")
    with temporary.open("x") as f:
        json.dump(obj,f,sort_keys=True,indent=2)
        f.write("\n");f.flush();os.fsync(f.fileno())
    os.rename(temporary,path)


def _load(path):
    return json.loads(Path(path).read_text())


def _run(argv):
    return subprocess.check_output(argv,text=True,stderr=subprocess.STDOUT).strip()


def _git():
    return _run(["git","-C",str(REPO),"rev-parse","HEAD"])


def _path(root,channel,index):
    return Path(root)/"results"/channel.lower()/f"chunk_{index:04d}"


def _input(root,channel,index):
    return Path(root)/"inputs"/channel.lower()/f"chunk_{index:04d}.csv"


def _native(work,name):
    return Path(work)/f"microscopic_{name}.csv"


def _read(path,fields=()):
    with Path(path).open(newline="") as f:
        rd=csv.DictReader(f)
        require(rd.fieldnames is not None and
                len(rd.fieldnames)==len(set(rd.fieldnames)) and
                set(fields).issubset(rd.fieldnames),
                f"Invalid native CSV header: {path}")
        yield from rd


def _original_dir(old,channel,index):
    camp=Path(old)/f"original2m_{channel.lower()}_10000_v1"
    chunk=f"chunk_{index:04d}"
    original=camp/"results"/chunk
    recovery=camp/"recovery_runtime_v1"/chunk
    def valid(d):
        return (d/"exit_status.txt").is_file() and (
            d/"exit_status.txt").read_text().strip()=="0" and (
            d/"verified_report.json").is_file()
    o,r=valid(original),valid(recovery)
    require(not(o and r),"Both old originals and recoveries succeeded")
    require(o or r,f"Old reference {channel} chunk {index} not complete")
    return (recovery if r else original),bool(r)


def _round_cpp(x):
    return math.floor(x+.5) if x>=0 else math.ceil(x-.5)


def _wire_centre(index,station):
    require(0<=station<=5,"Invalid MWPC station")
    n=math.floor((307.2 if station==0 else 512.)/PITCH_MM+.5)
    phase=PITCH_MM/2 if n%2==0 else 0.
    return phase+index*PITCH_MM


def _wire_base(u,station):
    phase=_wire_centre(0,station)
    return _round_cpp((u-phase)/PITCH_MM)


def _metadata_campaign(root):
    return Path(root)/"locked"/"selection.json"


def prepare(root,old):
    root=Path(root);old=Path(old)
    require(root.resolve()!=old.resolve()
            and not str(root.resolve()).startswith(str(old.resolve())+"/"),
            "New wire-truth output must be outside original 30k campaign")
    root.mkdir(parents=True,exist_ok=True)
    (root/"logs").mkdir(exist_ok=True)
    for channel in CHANNELS:
        (root/"inputs"/channel.lower()).mkdir(parents=True,exist_ok=True)
        (root/"results"/channel.lower()).mkdir(parents=True,exist_ok=True)
        (root/"attempts"/channel.lower()).mkdir(parents=True,exist_ok=True)
        (root/"claims"/channel.lower()).mkdir(parents=True,exist_ok=True)
    manifest={
        "schema":SCHEMA,"git_commit":_git(),
        "cxx_sha256":digest(SOURCE),
        "physics":{
            "geometry":"p3_noFW 9-wire analytic local patch",
            "pitch_mm":3.,"anode_diameter_um":30.,"cathode_v_mm":[-3.,3.],
            "anode_voltage_V":1800.,"gas":{"Ar":70,"CO2":30},
            "penning_r":0.57,"penning_lambda_cm":0.,
            "magnetic_field_T":0.,"track":"TrackHeed delta transport",
            "microscopic":"Garfield::AvalancheMicroscopic",
            "max_primaries_per_muon":50000,"cap_per_seed":200000,
            "thermal_start_eV":0.1,"birth_reservoir_per_seed":RESERVOIR,
            "ion_wire":"nearest local anode to ion birth u",
            "electron_wire":"actual Garfield electron endpoint nearest wire",
            "wire_global_index":"station phase + global integer*k*3mm",
            "w_reference":"global w in mm, ion birth w offset from Geant4 entry",
            "time_reference":"ion birth and collected electron times in ns",
        },
        "source_campaigns":{}, "n_muons":30000,"chunks_per_channel":500,
        "input_sha256":{},
    }
    for channel in CHANNELS:
        oldc=old/f"original2m_{channel.lower()}_10000_v1"
        m=metadata(oldc)
        require(m["schema"]==OLD_SCHEMA and m["channel"]==channel
                and m["n_muons"]==10000 and m["n_chunks"]==500,
                f"Old frozen 10k selection invalid: {channel}")
        manifest["source_campaigns"][channel]={
            "source_manifest_sha256":digest(oldc/"manifest.json"),
            "source_selection_manifest_sha256":digest(oldc/"selection_manifest.json"),
            "source_frozen_input_sha256":m["frozen_source_sha256"],
            "previous_physics_binary_sha256":m["binary_sha256"],
        }
        for index in range(500):
            _,rows=check_input(oldc,index)
            require(len(rows)==20,
                    "Old original Geant4 input does not contain 20 tracks")
            src=oldc/"chunks"/f"chunk_{index:04d}.csv"
            dst=_input(root,channel,index)
            sha=digest(src)
            if dst.exists():
                require(digest(dst)==sha,
                        f"Existing staged original Geant4 crossing changed: {dst}")
            else:
                # Original archive remains read-only. Atomic local copy.
                tmp=dst.with_name(dst.name+f".{os.getpid()}.part")
                shutil.copyfile(src,tmp)
                require(digest(tmp)==sha,"Copied Geant4 crossing SHA differs")
                os.rename(tmp,dst)
            manifest["input_sha256"][f"{channel}:{index:04d}"]=sha
    lock=_metadata_campaign(root)
    if lock.exists():
        require(_load(lock)==manifest,
                "Selection/source/git SHA lock mismatch. New campaign name needed.")
    else:
        _json(lock,manifest)
    print("PASS prepared 1,500 old-identical Geant4 input chunks "
          "(10,000 Jpsi + 10,000 Phi + 10,000 Omega), locked to",
          manifest["git_commit"],flush=True)


def _env_data():
    x={
        "gpp_version":_run(["g++","-dumpfullversion"]),
        "root_version":_run(["root-config","--version"]),
        "cmake_version":_run(["cmake","--version"]).splitlines()[0],
        "garfield_install":os.environ.get("GARFIELD_INSTALL",""),
    }
    require(x["garfield_install"],
            "Source the same CVMFS Garfield++-env.sh used in old PBS jobs")
    x["garfield_install_realpath"]=str(Path(x["garfield_install"]).resolve())
    # Preserve exact shared-library fingerprints, not just HEAD symlink names.
    ld=_run(["ldd",str(BINARY)])
    require("not found" not in ld,"Missing Garfield/ROOT shared libraries")
    libs={}
    for line in ld.splitlines():
        if "Garfield" in line:
            rhs=line.split("=>")
            require(len(rhs)==2 and rhs[1].strip().startswith("/"),
                    "Garfield++ shared-library location unknown: "+line)
            file=Path(rhs[1].strip().split(" (",1)[0])
            require(file.is_file(),"Missing Garfield runtime shared library")
            libs[file.name]={"path":str(file.resolve()),"sha256":digest(file)}
    x["garfield_shared_libraries"]=libs
    x["ldd"]=ld
    return x


def lock_binary(root):
    selection=_load(_metadata_campaign(root))
    require(selection["git_commit"]==_git()
            and selection["cxx_sha256"]==digest(SOURCE),
            "Production source differs from frozen Geant4 selection")
    require(BINARY.is_file(),"Build the existing wiretruth CMake target first")
    stamp={
        "schema":SCHEMA,"git_commit":_git(),"cxx_sha256":digest(SOURCE),
        "binary_path":str(BINARY),"binary_sha256":digest(BINARY),
        "toolchain":_env_data(),"utc_approval":"BUILD_FINISHED_GATE_REQUIRED",
    }
    dest=Path(root)/"locked"/"binary.json"
    if dest.exists():
        require(_load(dest)==stamp,
                "Binary/runtime differs from original locked production build")
    else:
        _json(dest,stamp)
    print("PASS locked full microscopic wiretruth binary SHA",
          stamp["binary_sha256"],flush=True)


def runtime_check(root):
    lock=_load(Path(root)/"locked"/"binary.json")
    sel=_load(_metadata_campaign(root))
    require(lock["git_commit"]==sel["git_commit"]==_git()
            and lock["cxx_sha256"]==sel["cxx_sha256"]==digest(SOURCE),
            "Production C++/Git version changed")
    require(lock["binary_sha256"]==digest(BINARY),
            "CVMFS native executable differs from locked binary")
    env=_env_data()
    for key in ("gpp_version","root_version","garfield_install_realpath"):
        require(env[key]==lock["toolchain"][key],
                "Runtime changed since build: "+key)
    require(env["garfield_shared_libraries"]==
            lock["toolchain"]["garfield_shared_libraries"],
            "Garfield library SHA/path differs; cannot mix physical versions")
    print("PASS locked commit, binary SHA, ROOT, compiler, Garfield runtime")


def _source_rows(root,channel,index):
    dst=_input(root,channel,index)
    sel=_load(_metadata_campaign(root))
    require(digest(dst)==sel["input_sha256"][f"{channel}:{index:04d}"],
            "Original Geant4 input chunk has changed")
    rs=list(_read(dst))
    ids=[crossing_id(r) for r in rs]
    require(len(rs)==len(set(ids))==20,
            "Source selection lost independent Geant4 parents")
    return {cid:r for cid,r in zip(ids,rs)}


def _asym(row,a,b,key):
    A=row[key]
    x,y=int(row[a]),int(row[b])
    if x+y==0:
        require(A=="","An empty-sector asymmetry should be undefined")
    else:
        require(A!="" and math.isclose(
                    float(A),(x-y)/(x+y),rel_tol=1e-8,abs_tol=1e-9),
                "Native signed ion-birth asymmetry definition changed")


def _verify_wire_row(r,source):
    cid=r["crossing_id"]
    require(cid in source,"Wire row not in selected Geant4 crossings")
    h=source[cid]
    st=int(h["station"])
    k=int(r["wire_index_global"])
    b=_wire_base(float(h["u_in_mm"]),st)
    require(abs(k-b)<=4,"Anode is outside physically simulated 9-wire patch")
    require(int(r["station"])==st and
            math.isclose(float(r["wire_u_global_mm"]),
                         _wire_centre(k,st),abs_tol=1e-8),
            "Global physical collecting/anode-region wire index mismatch")
    q,n=int(r["collected_e"]),int(r["ion_births"])
    require(q>=0 and n>=0 and (q or n),
            "Unphysical or duplicate empty per-wire row")
    ups=int(r["ion_up"])+int(r["ion_down"])+int(r["ion_boundary"])
    top=int(r["ion_top"])+int(r["ion_bottom"])+int(r["ion_side"])
    require(ups==top==n,"Anode ion-sector census differs from full birth count")
    _asym(r,"ion_up","ion_down","A_ion_hemi")
    _asym(r,"ion_top","ion_bottom","A_ion_TB")
    for prefix,num in (
          (("mean_w_birth_global_mm","sigma_w_birth_mm"),n),
          (("mean_t_birth_ns","sigma_t_birth_ns"),n),
          (("mean_t_collected_ns","sigma_t_collected_ns"),q)):
        for field in prefix:
            if num:
                v=float(r[field])
                require(math.isfinite(v) and
                        ("sigma_" not in field or v>=0),
                        "Invalid full-population anode moment: "+field)
            else:
                require(r[field]=="",
                        "Empty charge/ion population has nonempty moment")
    return q,n


def _original_heed_replay(old,channel,index,newrows,seedrows):
    """Compare new TrackHeed against previous *same-crossing* archived primaries.

    Matching starts, times and cluster indices are evidence of reproducibility;
    if not, mark the fresh Heed sample explicitly as independent.
    """
    directory,recovered=_original_dir(old,channel,index)
    file=directory/"penning_heed_primaries.csv"
    oldrows={(r["crossing_id"],int(r["primary_index"])):r
             for r in _read(file,OUTPUT_COLUMNS["heed_primaries"])}
    mismatch=0;max_abs_coord_cm=0.;max_abs_t_ns=0.
    if len(oldrows)!=len(newrows):
        mismatch+=abs(len(oldrows)-len(newrows))
    for key,r in newrows.items():
        o=oldrows.get(key)
        if o is None:
            mismatch+=1;continue
        if int(o["cluster_index"])!=int(r["cluster_index"]):
            mismatch+=1;continue
        diffs=[abs(float(o[f])-float(r[f])) for f in
               ("u_rel_cm","v_abs_cm","w_rel_cm")]
        dt=abs(float(o["t_ns"])-float(r["t_ns"]))
        max_abs_coord_cm=max(max_abs_coord_cm,*diffs)
        max_abs_t_ns=max(max_abs_t_ns,dt)
        if max(diffs)>1e-9 or dt>1e-7:
            mismatch+=1
    return {"matching":"matched" if mismatch==0 else "independent_new_Heed",
            "n_mismatching_primaries":mismatch,
            "old_recovery_chunk":recovered,
            "max_abs_coordinate_difference_cm":max_abs_coord_cm,
            "max_abs_time_difference_ns":max_abs_t_ns}


def verify(work,root,old,channel,index,pbs_id="unknown"):
    """Native CSV audit (stream huge birth reservoir); return scientific receipt."""
    require(channel in CHANNELS and 0<=index<500,"Channel/chunk index invalid")
    root=Path(root)
    source=_source_rows(root,channel,index)
    prefix=Path(work)/"microscopic"
    paths={p:_native(work,p) for p in PARTS}
    for p,file in paths.items():
        require(file.is_file(),"Missing native full Garfield output: "+str(file))
    bymu={}
    for r in _read(paths["crossing"],OUTPUT_COLUMNS["crossing"]):
        cid=r["source|chunk|event|hit_index|track_id|station|det_id"]
        require(cid in source and cid not in bymu,
                "Missing/duplicated original Geant4 crossing")
        bymu[cid]=r
    require(set(bymu)==set(source),"Not all 20 frozen Geant4 muons completed")
    clusters=defaultdict(list); primary={}
    for r in _read(paths["heed_clusters"],OUTPUT_COLUMNS["heed_clusters"]):
        cid=r["crossing_id"];j=int(r["cluster_index"])
        require(cid in source and j>=0 and int(r["n_original"])>=
                int(r["n_accepted"])>=0,"Invalid native Heed cluster")
        clusters[cid].append(r)
    for cid in source:
        cc=clusters[cid]
        require(sorted(int(r["cluster_index"]) for r in cc)==
                list(range(len(cc))) and len(cc)==int(bymu[cid]["heed_clusters"]),
                "Original Heed cluster identity/count mismatch")
    for r in _read(paths["heed_primaries"],OUTPUT_COLUMNS["heed_primaries"]):
        key=(r["crossing_id"],int(r["primary_index"]))
        require(key[0] in source and key not in primary,"Duplicate Heed electron")
        j=int(r["cluster_index"])
        require(0<=j<len(clusters[key[0]]),"Heed primary with missing cluster")
        primary[key]=r
    for cid in source:
        ps=[j for c,j in primary if c==cid]
        require(sorted(ps)==list(range(len(ps))) and
                len(ps)==int(bymu[cid]["heed_primaries"])==
                sum(int(r["n_accepted"]) for r in clusters[cid]),
                "Original Heed count/primary index mismatch")
    seeds={}
    for r in _read(paths["seed"],OUTPUT_COLUMNS["seed"]):
        key=(r["crossing_id"],int(r["primary_index"]))
        require(key in primary and key not in seeds,"Unknown or repeated Garfield seed")
        i=key[1]
        h=source[key[0]]
        rng=1+((int(h["heed_seed"])+(i+1)*104729)%2000000000)
        require(int(r["garfield_seed"])==rng,
                "Microscopic avalanche RNG seed changed")
        require(int(r["capped"])==0,"CAPPED Garfield event: reject campaign chunk")
        seeds[key]=r
    require(len(seeds)==len(primary),"Not all Heed electrons propagated")
    # Validate exact independent collecting-endpoint wire and nearest-birth
    # wire totals, including populations where q=0 or n=0.
    seedwire={}; seedq=defaultdict(int);seedn=defaultdict(int)
    seedsectors=defaultdict(lambda:Counter())
    for r in _read(paths["wire_seed"],OUTPUT_COLUMNS["wire_seed"]):
        key=(r["crossing_id"],int(r["primary_index"]))
        require(key in seeds and
                int(r["cluster_index"])==int(primary[key]["cluster_index"]),
                "Per-wire seed/cluster identity mismatch")
        wk=key+(int(r["wire_index_global"]),)
        require(wk not in seedwire,"Duplicate microscopic collecting/ion-region wire")
        q,n=_verify_wire_row(r,source)
        seedwire[wk]=r;seedq[key]+=q;seedn[key]+=n
        for suffix in ("ion_up","ion_down","ion_top","ion_bottom"):
            seedsectors[key][suffix]+=int(r[suffix])
    for key,r in seeds.items():
        require(seedq[key]==int(r["collected_e"]) and
                seedn[key]==int(r["ion_births"]),
                "Per-seed full Garfield endpoint/birth wire conservation fails")
        for field in ("ion_up","ion_down","ion_top","ion_bottom"):
            require(seedsectors[key][field]==int(r[field]),
                    "Per-seed signed ion asymmetry changed: "+field)
    muwire={}; muq=defaultdict(int);mun=defaultdict(int)
    for r in _read(paths["wire_muon"],OUTPUT_COLUMNS["wire_seed"]):
        cid=r["crossing_id"]
        require(int(r["primary_index"])==-1 and
                int(r["cluster_index"])==-1,
                "Muon-wire rows must be aggregated across seeds")
        wk=(cid,int(r["wire_index_global"]))
        require(wk not in muwire,"Duplicate physical wire for one muon")
        q,n=_verify_wire_row(r,source)
        muwire[wk]=r;muq[cid]+=q;mun[cid]+=n
        # Verify the muon wire aggregates exactly equal the independent
        # original Heed seed totals on the SAME physical global anode.
        subtotal=[s for key,s in seedwire.items()
                  if key[0]==cid and key[2]==wk[1]]
        require(q==sum(int(s["collected_e"]) for s in subtotal)
                and n==sum(int(s["ion_births"]) for s in subtotal),
                "Per-muon and per-seed wire census differ")
        # Reconstruct full-population per-wire moments independently from
        # each microscopic seed's exact moments (no ion sampling involved).
        for nfield,meanfield,sigmafield in (
            ("ion_births","mean_w_birth_global_mm","sigma_w_birth_mm"),
            ("ion_births","mean_t_birth_ns","sigma_t_birth_ns"),
            ("collected_e","mean_t_collected_ns","sigma_t_collected_ns")):
            population=int(r[nfield])
            if not population:
                continue
            mean=sum(int(s[nfield])*float(s[meanfield])
                     for s in subtotal if int(s[nfield]))/population
            # Centre on the aggregate mean to avoid cancellation for
            # a narrow ~micrometre cloud at large global w coordinate.
            variance=sum(
                int(s[nfield])*(float(s[sigmafield])**2+
                                (float(s[meanfield])-mean)**2)
                for s in subtotal if int(s[nfield]))/population
            sigma=math.sqrt(max(0.,variance))
            require(math.isclose(float(r[meanfield]),mean,
                                 rel_tol=1e-7,abs_tol=1e-6) and
                    math.isclose(float(r[sigmafield]),sigma,
                                 rel_tol=1e-5,abs_tol=1e-5),
                    "Original exact wire moments differ after per-muon "
                    "weighted seed aggregation: "+meanfield)
    for cid,r in bymu.items():
        require(muq[cid]==int(r["collected_avalanche_electrons"])==
                sum(int(s["collected_e"]) for (c,j),s in seeds.items() if c==cid)
                and mun[cid]==int(r["ion_births"])==
                sum(int(s["ion_births"]) for (c,j),s in seeds.items() if c==cid),
                "Per-muon original Garfield conservation fails")
        zero=sum(int(s["collected_e"])==0 for (c,j),s in seeds.items() if c==cid)
        fail=sum(int(s["successful"])!=1 for (c,j),s in seeds.items() if c==cid)
        require(zero==int(r["zero_seeds"]) and fail==int(r["failed_seeds"])
                and int(r["capped_seeds"])==0,
                "Zero-gain/failed/capped physics status mismatch")
    # Streaming ion reservoir audit. Its weights are N_exact/n_sample for
    # a UNIFORM sample across the entire individual avalanche; not per-wire.
    sample=defaultdict(int)
    for r in _read(paths["ion_birth_sample"],OUTPUT_COLUMNS["ion_birth_sample"]):
        key=(r["crossing_id"],int(r["primary_index"]))
        require(key in seeds,"Reservoir birth without Garfield primary")
        n=int(seeds[key]["ion_births"])
        require(n>0 and int(r["seed_ion_births_total"])==n,
                "Reservoir physical ion birth count mismatch")
        require(all(math.isfinite(float(r[f])) for f in
                    ("u_mm","v_mm","w_mm","t_ns","reservoir_weight")),
                "Non-finite uniform birth sample")
        expected_weight=n/min(n,RESERVOIR)
        require(math.isclose(float(r["reservoir_weight"]),
                             expected_weight,rel_tol=1e-7,abs_tol=1e-7),
                "Reservoir statistical weight is not N_ion/N_retained")
        sample[key]+=1
        require(sample[key]<=min(n,RESERVOIR),
                "More birth positions than saved unbiased sample capacity")
    for key,seed in seeds.items():
        require(sample[key]==min(int(seed["ion_births"]),RESERVOIR),
                "Incomplete uniform sample of ion births")
    # Timing counters remain independently joined to the original seeds.
    times={}
    for r in _read(paths["timing"],OUTPUT_COLUMNS["timing"]):
        key=(r["crossing_id"],int(r["primary_index"]))
        require(key in seeds and key not in times,
                "Original per-primary timing duplicated or lost")
        s=seeds[key]
        require(int(r["garfield_seed"])==int(s["garfield_seed"])
                and int(r["ion_births"])==int(s["ion_births"])
                and int(r["collected_e"])==int(s["collected_e"])
                and int(r["capped"])==0
                and int(r["success"])==int(s["successful"])
                and math.isfinite(float(r["wall_seconds"]))
                and float(r["wall_seconds"])>=0.,
                "Original full microscopic timing physics mismatch")
        times[key]=r
    require(len(times)==len(seeds),"Per-avalanche wall time rows incomplete")
    mu_times={}
    for r in _read(paths["muon_timing"],OUTPUT_COLUMNS["muon_timing"]):
        cid=r["crossing_id"]
        require(cid in source and cid not in mu_times,
                "Original muon wall-clock time duplicated")
        for f in ("heed_wall_seconds","garfield_and_output_wall_seconds",
                  "muon_total_wall_seconds"):
            require(math.isfinite(float(r[f])) and float(r[f])>=0.,
                    "Invalid original muon wall-clock duration")
        mu_times[cid]=r
    require(len(mu_times)==20,"Missing per-muon full microphysics times")
    for cid,r in mu_times.items():
        require(int(r["n_Heed_clusters"])==len(clusters[cid])
                and int(r["n_Heed_primaries"])==
                sum(c==cid for c,j in seeds),
                "Muon original Heed timing/count mismatch")
    replay=_original_heed_replay(old,channel,index,primary,seeds)
    # A zero-charge muon remains part of the independent 20 Geant4 muon sample.
    one=two=many=zero=0
    q2rows=[]
    for cid,src in source.items():
        q=sorted((int(r["collected_e"]) for (c,k),r in muwire.items()
                 if c==cid and int(r["collected_e"])>0),reverse=True)
        if len(q)==0:zero+=1
        elif len(q)==1:one+=1
        elif len(q)==2:two+=1
        else:many+=1
        q2rows.append({
            "crossing_id":cid, "channel":channel,
            "theta_u_deg":float(src["theta_u_in_deg"]),
            "n_active_wires":len(q),"Qtotal":sum(q),
            "Q2_over_Qtotal":q[1]/sum(q) if len(q)>=2 else 0.,
        })
    require(one+two+many+zero==20,"Inconsistent multi-anode census")
    counts={
       "muons":20,"heed_primaries":len(seeds),
       "heed_clusters":sum(map(len,clusters.values())),
       "collected_e":sum(int(r["collected_e"]) for r in seeds.values()),
       "ion_births":sum(int(r["ion_births"]) for r in seeds.values()),
       "zero_gain":sum(int(r["collected_e"])==0 for r in seeds.values()),
       "failed_return":sum(int(r["successful"])!=1 for r in seeds.values()),
       "capped":0,"no_charge_muons":zero,"single_wire_muons":one,
       "two_wire_muons":two,"three_or_more_wire_muons":many,
       "reservoir_points":sum(sample.values()),
       "positive_wire_seed_rows":len(seedwire),
       "positive_wire_muon_rows":len(muwire),
    }
    paths_sha={p:{"sha256":digest(f),"bytes":f.stat().st_size}
               for p,f in paths.items()}
    return {
        "schema":SCHEMA,"channel":channel,"chunk":index,"counts":counts,
        "input_sha256":digest(_input(root,channel,index)),
        "binary_sha256":_load(Path(root)/"locked"/"binary.json")["binary_sha256"],
        "git_commit":_git(),"cxx_sha256":digest(SOURCE),
        "source_old_heed_replay":replay,
        "files":paths_sha,"q2_event_summary":q2rows,
        "pbs_job_id":pbs_id,
    }


def verify_command(opts):
    x=verify(opts.work,opts.root,opts.old_root,opts.channel,opts.index,
             os.environ.get("PBS_JOBID","local"))
    receipt=Path(opts.work)/"verified_report.json"
    _json(receipt,x)
    print("PASS verified 20 native muons; conservation of real endpoints, "
          "ion births, station wire IDs, seed→muon sums, reservoir weights; "
          f"new Heed replay={x['source_old_heed_replay']['matching']}; "
          f"{x['counts']['heed_primaries']} microscopic avalanches",flush=True)


def audit_work(root,channel,index,dest,hash_files=True):
    dest=Path(dest)
    receipt=dest/"verified_report.json"
    require(receipt.is_file(),f"Missing completed receipt: {receipt}")
    x=_load(receipt)
    require(x["schema"]==SCHEMA and x["channel"]==channel
            and x["chunk"]==index,
            "Wrong completed native production identity")
    sel=_load(_metadata_campaign(root))
    binary=_load(Path(root)/"locked"/"binary.json")
    require(x["input_sha256"]==sel["input_sha256"][f"{channel}:{index:04d}"]
            and x["binary_sha256"]==binary["binary_sha256"]
            and x["git_commit"]==binary["git_commit"]
            and x["cxx_sha256"]==binary["cxx_sha256"],
            "Chunk was processed with different physics source/binary")
    if hash_files:
        for p,item in x["files"].items():
            file=_native(dest,p)
            require(file.is_file() and file.stat().st_size==item["bytes"]
                    and digest(file)==item["sha256"],
                    "Native full Garfield result checksum changed: "+str(file))
    return x


def audit_one(root,channel,index,hash_files=True):
    return audit_work(root,channel,index,
                      _path(root,channel,index),hash_files)


def status(root):
    for ch in CHANNELS:
        good=[];partial=[]
        for idx in range(500):
            p=_path(root,ch,idx)
            if (p/"verified_report.json").is_file():
                good.append(idx)
            elif p.exists():partial.append(idx)
        missing=[i for i in range(500) if i not in set(good)]
        print(f"{ch}: {len(good):3d}/500 receipts, "
              f"first_missing={missing[0] if missing else 'none'}; "
              f"incomplete={partial[:20]}")


def gate(root,gate_name):
    """Only creates a reviewable pilot gate; it does NOT submit jobs."""
    root=Path(root)
    if gate_name=="20":
        requested={"Jpsi":[0]}
    elif gate_name=="1000":
        requested={"Jpsi":list(range(20)),"Phi":list(range(15)),
                   "Omega":list(range(15))}
    elif gate_name=="full":
        requested={c:list(range(500)) for c in CHANNELS}
    else:
        raise RuntimeError("Invalid scientific gate stage")
    records=[]
    for ch,indices in requested.items():
        for i in indices:
            records.append(audit_one(root,ch,i,hash_files=True))
    sums=Counter()
    for r in records:sums.update(r["counts"])
    require(sums["muons"]==sum(map(len,requested.values()))*20 and
            sums["capped"]==0,"Pilot incomplete or has capped physics")
    for k in ("single_wire_muons","two_wire_muons",
              "three_or_more_wire_muons","no_charge_muons"):
        require(sums[k]>=0,"Negative active-wire census")
    # HPC constraints are actual limits from the cluster's PBS script:
    # 4 GB/job, 8 h/job. Measuring 1000 muons is a resource prerequisite,
    # never infer full production memory/time solely from average physics.
    resources=[]
    for ch,indices in requested.items():
        for i in indices:
            file=_path(root,ch,i)/"resource_time.txt"
            require(file.is_file(),
                    "Missing /usr/bin/time -v resource measurement: "+str(file))
            content=file.read_text()
            match=re.search(r"Maximum resident set size \(kbytes\):\s*(\d+)",
                            content)
            require(match is not None,
                    "Cannot read pilot peak memory from "+str(file))
            # Prefer measured PBS-native worker's resource file. elapsed may
            # also appear as h:mm:ss or m:ss; use external wall time if needed.
            elapsed=re.search(r"Elapsed \(wall clock\) time[^\n]*?:\s*"
                              r"((?:\d+:){1,2}\d+(?:\.\d+)?)\s*$",
                              content,re.MULTILINE)
            raw=elapsed.group(1).strip() if elapsed else ""
            seconds=0.
            if raw:
                nums=[float(n) for n in raw.split(":")]
                for number in nums:seconds=seconds*60+number
            require(seconds>0.,
                    "Cannot parse true Garfield elapsed wall time: "+str(file))
            row=next(t for t in records
                     if t["channel"]==ch and t["chunk"]==i)
            resources.append({
                "channel":ch,"chunk":i,"max_rss_kib":int(match.group(1)),
                "garfield_wall_seconds":seconds,
                "native_csv_bytes":sum(
                    v["bytes"] for v in row["files"].values()),
            })
    require(max(r["max_rss_kib"] for r in resources)<3800000
            and max(r["garfield_wall_seconds"] for r in resources)<8*3600,
            "Pilot exceeded 4 GB/8h PBS limits with safety margin; "
            "adjust job request only after cluster resource review")
    import statistics
    avg=statistics.mean(r["garfield_wall_seconds"] for r in resources)
    avg_bytes=statistics.mean(r["native_csv_bytes"] for r in resources)
    gate_report={
        "schema":SCHEMA,"gate":gate_name,"chunks":len(records),
        "measured_resources":{
            "max_rss_kib":max(r["max_rss_kib"] for r in resources),
            "max_chunk_garfield_wall_seconds":max(
                r["garfield_wall_seconds"] for r in resources),
            "mean_chunk_garfield_wall_seconds":avg,
            "projected_1500_chunk_cpu_wall_hours":avg*1500/3600,
            "projected_1500_chunk_raw_GiB":avg_bytes*1500/(1024**3),
            "based_on_exactly_n_tested_chunks":len(resources),
            "caution":"Projection only; high-Heed tail/time and disk distribution "
                      "must be monitored during full PBS campaign",
        },
        "counts":dict(sums),
        "source_git_commit":_git(),
        "binary_sha256":_load(root/"locked"/"binary.json")["binary_sha256"],
        "heed_replay_match":{
           "matched":sum(r["source_old_heed_replay"]["matching"]=="matched"
                         for r in records),
           "independent_new_Heed":sum(
                 r["source_old_heed_replay"]["matching"]!="matched"
                 for r in records),
        },
    }
    dest=root/"gates"/f"gate_{gate_name}.json"
    if dest.exists():
        require(_load(dest)==gate_report,
                "Existing pilot gate differs: inspect, never overwrite")
    else:
        _json(dest,gate_report)
    print(json.dumps(gate_report,indent=2),flush=True)
    print("PASS scientific audit gate",gate_name,"recorded at",dest,
          "(No PBS jobs submitted.)",flush=True)


def report(root,old_root):
    root=Path(root)
    sums={c:Counter() for c in CHANNELS}
    bins=[0,5,10,15,20,30,60,180]
    hist=defaultdict(lambda:{"events":0,"charge_events":0,
                            "second_gt_5pct":0,"q2_sum":0.,
                            "n0":0,"n1":0,"n2":0,"n3plus":0})
    run_list=[]
    for c in CHANNELS:
        for i in range(500):
            r=audit_one(root,c,i,hash_files=False)
            sums[c].update(r["counts"])
            run_list.append({
                "channel":c,"chunk":i,"pbs_job_id":r["pbs_job_id"],
                "git_commit":r["git_commit"],
                "binary_sha256":r["binary_sha256"],
                "files":r["files"],"replay":r["source_old_heed_replay"],
            })
            for ev in r["q2_event_summary"]:
                a=abs(ev["theta_u_deg"])
                angle=next((f"{bins[j]}-{bins[j+1]}"
                            for j in range(len(bins)-1)
                            if bins[j]<=a<bins[j+1]),">=180")
                for group in (c,"All"):
                    v=hist[(group,angle)]
                    v["events"]+=1
                    v["charge_events"]+=int(ev["Qtotal"]>0)
                    v["second_gt_5pct"]+=int(ev["Q2_over_Qtotal"]>0.05)
                    v["q2_sum"]+=ev["Q2_over_Qtotal"]
                    active=int(ev["n_active_wires"])
                    v["n0"]+=int(active==0)
                    v["n1"]+=int(active==1)
                    v["n2"]+=int(active==2)
                    v["n3plus"]+=int(active>=3)
    require(sum(s["muons"] for s in sums.values())==30000,
            "Cannot report incomplete campaign")
    old_sums={ch:Counter() for ch in CHANNELS}
    for ch in CHANNELS:
        for index in range(500):
            olddir,repaired=_original_dir(old_root,ch,index)
            cross=olddir/"penning_crossing.csv"
            seal=olddir/"outputs.sha256"
            frozen={}
            for line in seal.read_text().splitlines():
                if len(line)>=67 and line[64:66]=="  ":
                    frozen[line[66:]]=line[:64]
            require(frozen.get("penning_crossing.csv")==digest(cross),
                    "Original 30k crossing SHA mismatch: "+str(cross))
            original=list(_read(cross,OUTPUT_COLUMNS["crossing"]))
            require(len(original)==20,
                    "Original 30k chunk no longer contains 20 real muons")
            old_sums[ch]["muons"]+=len(original)
            old_sums[ch]["heed_primaries"]+=sum(
                int(x["heed_primaries"]) for x in original)
            old_sums[ch]["ion_births"]+=sum(
                int(x["ion_births"]) for x in original)
            old_sums[ch]["collected_e"]+=sum(
                int(x["collected_avalanche_electrons"]) for x in original)
    folder=root/"compact"
    folder.mkdir(parents=True,exist_ok=True)
    with (folder/"previous_30k_sanity_comparison.csv").open(
            "x",newline="") as f:
        writer=csv.writer(f)
        writer.writerow(("channel","metric","old_30k","new_wiretruth_30k",
                         "difference","ratio_new_over_old"))
        for ch in CHANNELS:
            for field in ("muons","heed_primaries","collected_e","ion_births"):
                oldn=old_sums[ch][field];newn=sums[ch][field]
                writer.writerow((
                    ch,field,oldn,newn,newn-oldn,
                    newn/oldn if oldn else "",
                ))
    with (folder/"wire_sharing_by_angle.csv").open("x",newline="") as f:
        w=csv.writer(f);w.writerow([
            "channel","abs_theta_u_bin_deg","muons","muons_with_charge",
            "n0_active_wires","n1_active_wire","n2_active_wires",
            "n3plus_active_wires","n_Q2_fraction_gt_0p05",
            "mean_Q2_fraction_inclusive","muon_fraction_Q2_gt_0p05"])
        for (channel,angle),v in sorted(hist.items()):
            require(v["n0"]+v["n1"]+v["n2"]+v["n3plus"]==
                    v["events"],"Active-wire angle-bin census mismatch")
            w.writerow([channel,angle,v["events"],v["charge_events"],
                       v["n0"],v["n1"],v["n2"],v["n3plus"],
                       v["second_gt_5pct"],v["q2_sum"]/v["events"],
                       v["second_gt_5pct"]/v["events"]])
    _json(folder/"run_manifest.json",{
        "schema":SCHEMA,"source_selection":_load(_metadata_campaign(root)),
        "compiled_binary":_load(root/"locked"/"binary.json"),
        "channels":{c:dict(s) for c,s in sums.items()},
        "previous_30k_archived_counts":{c:dict(old_sums[c]) for c in CHANNELS},
        "completed_pbs_jobs_and_file_checksums":run_list,
        "n_jobs":len(run_list),"comparison_rule":
          "Jpsi calibration only; Phi/Omega holdout independent validation",
        "physics_note":
          "Wire-collected electron charge is endpoint truth; ion-region "
          "nearest wire is a separate ion-birth classification. "
          "Ion birth reservoir max 128 per seed; convergence for "
          "cathode signals NOT established.",
    })
    print("PASS completed 30k production",json.dumps(
          {c:dict(s) for c,s in sums.items()},indent=2),flush=True)


def export(root):
    """Compress physics tables only, retaining original native outputs on HPC."""
    root=Path(root)
    folder=root/"compact"
    require((folder/"run_manifest.json").exists(),
            "Run and validate report after full 30k scientific gate")
    names=("crossing","seed","wire_seed","wire_muon","heed_primaries",
           "heed_clusters","timing","muon_timing")
    for ch in CHANNELS:
        for name in names:
            dest=folder/f"{ch.lower()}_{name}.csv.gz"
            require(not dest.exists(),
                    "Refuse to overwrite existing compact export: "+str(dest))
            temp=dest.with_name(dest.name+f".{os.getpid()}.part")
            with gzip.open(temp,"wt",newline="",compresslevel=5) as f:
                first=True;header=None
                for index in range(500):
                    r=audit_one(root,ch,index,hash_files=False)
                    with _native(_path(root,ch,index),name).open() as src:
                        h=src.readline()
                        if first:
                            header=h;f.write(h);first=False
                        else:
                            require(h==header,
                                    "Native CSV column order drift across chunks")
                        shutil.copyfileobj(src,f)
            os.rename(temp,dest)
            print("EXPORTED",dest,dest.stat().st_size,flush=True)
    # Full authentic Geant4 entrance/exit geometry and native Heed seeds:
    # include the frozen original selection inputs, not only C++ aggregate
    # crossing outputs (which have less spatial information).
    for ch in CHANNELS:
        dst=folder/f"{ch.lower()}_geant4_crossings.csv.gz"
        require(not dst.exists(),"Refusing to overwrite original G4 selection export")
        with gzip.open(dst,"wt",newline="",compresslevel=5) as f:
            first=True;header=None
            for index in range(500):
                src=_input(root,ch,index)
                with src.open() as stream:
                    h=stream.readline()
                    if first:
                        header=h;f.write(h);first=False
                    else:
                        require(h==header,
                                "Original input column order changed")
                    shutil.copyfileobj(stream,f)
        print("EXPORTED frozen 3D Geant4 crossing rows",dst,flush=True)
    manifest=_load(folder/"run_manifest.json")
    table_hashes={p.name:digest(p) for p in folder.glob("*.csv.gz")}
    _json(folder/"export_sha256.json",table_hashes)
    print("PASS 27 compact gzip physics tables, including authentic Geant4 "
          "inputs; excluding full ion reservoirs")


def selftest():
    # Physical global anode map: station 0 has even N_anodes and a +1.5
    # mm phase; all other stations have odd N_anodes and zero phase.
    assert math.isclose(_wire_centre(0,0),1.5)
    assert math.isclose(_wire_centre(0,1),0.)
    assert math.isclose(_wire_centre(-1,0),-1.5)
    assert _wire_base(1.0,0)==0
    assert _wire_base(-1.6,0)==-1
    assert _wire_base(1.0,1)==0
    assert _round_cpp(-.5)==-1 and _round_cpp(.5)==1
    crossing="synthetic|1|2|3|4|0|1"
    src={crossing:{"station":"0","u_in_mm":"1.0"}}
    row={
        "crossing_id":crossing,"station":"0",
        "wire_index_global":"0","wire_u_global_mm":"1.5",
        "collected_e":"5","ion_births":"7",
        "ion_up":"3","ion_down":"4","ion_boundary":"0",
        "ion_top":"3","ion_bottom":"2","ion_side":"2",
        "A_ion_hemi":str(-1/7),"A_ion_TB":str(1/5),
        "mean_w_birth_global_mm":"2.","sigma_w_birth_mm":".5",
        "mean_t_birth_ns":"1.5","sigma_t_birth_ns":".3",
        "mean_t_collected_ns":"3.","sigma_t_collected_ns":".4",
    }
    assert _verify_wire_row(row,src)==(5,7)
    bad={**row,"wire_index_global":"1"}
    try:
        _verify_wire_row(bad,src)
    except RuntimeError:
        pass
    else:
        raise RuntimeError("Global physical wire phase guard did not reject "
                           "a mismapped anode")
    for text,seconds in (
        ("Elapsed (wall clock) time (h:mm:ss or m:ss): 0:10:00.32",600.32),
        ("Elapsed (wall clock) time (h:mm:ss or m:ss): 12:34.5",754.5),
    ):
        match=re.search(r"Elapsed \(wall clock\) time[^\n]*?:\s*"
                        r"((?:\d+:){1,2}\d+(?:\.\d+)?)\s*$",
                        text,re.MULTILINE)
        require(match is not None,"GNU time -v elapsed parse self-test")
        parsed=0.
        for number in [float(n) for n in match.group(1).split(":")]:
            parsed=parsed*60+number
        require(math.isclose(parsed,seconds,abs_tol=1e-6),
                "GNU elapsed resource units self-test failed")
    print("PASS physical station-0 half-pitch versus other-station anode "
          "phase, global wire IDs, charge/ion census and elapsed resource "
          "units, without starting Garfield++")


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("command",choices=("selftest","prepare","lock-binary","runtime-check",
                                       "verify","audit-one","audit-work",
                                       "status","gate","report","export"))
    p.add_argument("--root",type=Path,default=ROOT_DEFAULT)
    p.add_argument("--old-root",type=Path,default=OLD_DEFAULT)
    p.add_argument("--channel",choices=CHANNELS)
    p.add_argument("--index",type=int)
    p.add_argument("--work",type=Path)
    p.add_argument("--stage",choices=("20","1000","full"))
    x=p.parse_args()
    require(x.command in ("selftest","prepare","lock-binary",
                           "runtime-check","status","report","export")
            or x.channel or x.stage,
            "Missing --channel or --stage")
    if x.command=="selftest":selftest()
    elif x.command=="prepare":prepare(x.root,x.old_root)
    elif x.command=="lock-binary":lock_binary(x.root)
    elif x.command=="runtime-check":runtime_check(x.root)
    elif x.command=="verify":
        require(x.work is not None and x.index is not None and x.channel,
                "verify requires --work --index --channel")
        verify_command(x)
    elif x.command in ("audit-one","audit-work"):
        require(x.index is not None and x.channel,
                "audit requires --index and --channel")
        if x.command=="audit-work":
            require(x.work is not None,"audit-work needs --work")
            r=audit_work(x.root,x.channel,x.index,x.work,True)
        else:
            r=audit_one(x.root,x.channel,x.index,True)
        print("PASS sealed chunk SHA/physics identity audit",
              r["channel"],r["chunk"])
    elif x.command=="status":status(x.root)
    elif x.command=="gate":
        require(x.stage,"gate requires --stage 20/1000/full")
        gate(x.root,x.stage)
    elif x.command=="report":report(x.root,x.old_root)
    elif x.command=="export":export(x.root)


if __name__=="__main__":
    try: main()
    except (OSError,RuntimeError,ValueError,KeyError,AssertionError,
            subprocess.CalledProcessError) as exc:
        print("WIRETRUTH30K FAIL:",exc,file=sys.stderr,flush=True)
        sys.exit(1)
