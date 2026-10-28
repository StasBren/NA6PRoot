#!/usr/bin/env python3
"""Bounded overnight native Garfield AvalancheElectron run on frozen LUT500 starts.

Physics: ORIGINAL exact saved uniformly sampled electron starts in every one
of 4200 pitch/field-wire cells, 20 distinct LUT500 trial indices per cell.
Each avalanche is a NEW independent full-start microscopic simulation; it is
NOT a continuation of the old drift and NOT a true Geant4-muon event.

Safety: finite 8.75h default wall budget; pass-by-pass random interleaving;
400-trial immutable SHA-stamped chunks; resume by verifying each completed
chunk; 5 ion-birth reservoir positions per full microscopic avalanche;
per-geometry full CSV outputs and per-cell statistics for Google Drive.
No original LUT/Geant4 or earlier B1 run files are modified.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import random
import shutil
import statistics
import subprocess
import sys
import time
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

from pi_audit_avalanche_lut_reuse import BASE, CONFIGS, EXPORT_FIELDS, locations, sha
from pi_run_avalanche_lut_sidecar import (
    COLUMNS, BIRTH_COLUMNS, TIME_COLUMNS, load_index, validate_outputs,
)
from pi_run_native_muon_b1 import CXX, EXE, check_executable, commit
from pi_garfield_runtime_env import garfield_runtime_env

STUDY="PI_LUT500_OVERNIGHT_20_PER_CELL_V1"
ROOT=BASE/"pi_lut500_avalanche_overnight_v1"/"run_20percell_firstpass_v1"
PER_CELL=20
CHUNK_SIZE=400
RESERVOIR=5
STRIDE=137  # coprime to 500, giving distinct samples for every pass
ORDER_SEED=20261029
GEOM_NAMES=tuple(x[0] for x in CONFIGS)


def canonical(obj):
    return json.dumps(obj,indent=2,sort_keys=True)+"\n"


def offset_for(geom,box_id):
    payload=f"NA60-B1-LUT500-OVERNIGHT-V1|{geom}|{box_id}".encode()
    return int.from_bytes(hashlib.sha256(payload).digest()[:8],"big")%500


def selected_j(geom,box_id,pass_id):
    if not 0<=pass_id<PER_CELL:
        raise ValueError("Bad avalanche pass")
    return (offset_for(geom,box_id)+STRIDE*pass_id)%500


def source_maps():
    # The already certified LUT20/500 sources and 432-index SHA provenance
    # must still be byte-for-byte identical on the user's WSL.
    frozen,already_indexed=load_index()
    meta={x["geometry"]:x for x in frozen["input_luts"] if x["trials"]==500}
    if set(meta)!=set(GEOM_NAMES):
        raise RuntimeError("No four SHA-pinned old LUT500 source grids")
    return meta, frozen["index_csv_sha256"]


def generate_selected(meta):
    """Stream all original LUT500 trials; retain ONLY 20 per cell.

    All four geometries and every old sample have equal opportunity.
    Each retained electron is identified by old LUT SHA, box_id, sample.
    No old 'collected' flag is used as a gate.
    """
    passes=[[] for _ in range(PER_CELL)]
    total=0
    for geom,pitch,fw in CONFIGS:
        conf,raw=locations(500,pitch,fw)
        if sha(raw)!=meta[geom]["lut_csv_sha256"]:
            raise RuntimeError("Frozen LUT500 bytes changed during scan: "+str(raw))
        nbox=round(pitch/.1)*30
        mapping=[{
            selected_j(geom,b,p):p for p in range(PER_CELL)
        } for b in range(nbox)]
        for b in range(nbox):
            if len(mapping[b])!=PER_CELL:
                raise RuntimeError("LUT500 selection index collision")
        found=bytearray(nbox*PER_CELL)
        local=0
        with raw.open(newline="",encoding="utf-8") as stream:
            rd=csv.DictReader(stream)
            for key in ("box_id","sample","iu","iv","random_seed",
                        "u_start_mm","v_start_mm","w_start_mm",
                        "drift_ok","collected","field_wire_collected",
                        "wire_index","phi_hit_deg","drift_time_ns",
                        "last_step_energy_eV"):
                if key not in (rd.fieldnames or ()):
                    raise ValueError("Missing ORIGINAL LUT500 column "+key)
            for row in rd:
                b=int(row["box_id"]);s=int(row["sample"])
                if b<0 or b>=nbox or s<0 or s>=500:
                    raise ValueError("Out-of-range LUT500 trial")
                p=mapping[b].get(s)
                if p is None:continue
                at=b*PER_CELL+p
                if found[at]:raise ValueError("Duplicate old LUT trial in one pass")
                found[at]=1
                iu=int(row["iu"]);iv=int(row["iv"])
                if b!=iu*30+iv:
                    raise ValueError("Old LUT500 box coordinate mismatch")
                u=float(row["u_start_mm"]);v=float(row["v_start_mm"])
                if (not math.isfinite(u) or not math.isfinite(v)
                    or not (-pitch/2+iu*.1-1.e-6<=u<=-pitch/2+(iu+1)*.1+1.e-6)
                    or not (-3+iv*.2-1.e-6<=v<=-3+(iv+1)*.2+1.e-6)):
                    raise ValueError("Old LUT500 exact electron start outside original box")
                rec={
                    "lut_trials":500,"geometry":geom,
                    "original_lut_sha256":meta[geom]["lut_csv_sha256"],
                    "box_id":b,"iu":iu,"iv":iv,"sample":s,
                    "random_seed":int(row["random_seed"]),
                    "u_start_mm":row["u_start_mm"],
                    "v_start_mm":row["v_start_mm"],
                    "w_start_mm":row["w_start_mm"],
                    "old_drift_ok":row["drift_ok"],
                    "old_collected":row["collected"],
                    "old_field_wire_collected":row["field_wire_collected"],
                    "old_wire_index":row["wire_index"],
                    "old_phi_hit_deg":row["phi_hit_deg"],
                    "old_drift_time_ns":row["drift_time_ns"],
                    "old_last_step_energy_eV":row["last_step_energy_eV"],
                }
                passes[p].append(rec)
                total+=1;local+=1
        if not all(found) or local!=nbox*PER_CELL:
            raise RuntimeError(f"Incomplete old LUT500 stratified selection {geom}: {local}")
        print(f"PASS {geom}: selected {local} independent exact original "
              "LUT500 trial starts, no old collected gating",flush=True)
    if total!=84000 or any(len(x)!=4200 for x in passes):
        raise RuntimeError("Not exactly 20 trials of 4200 original cells")
    # Random mixing INSIDE each full-cell pass. If time expires, all
    # prior passes are complete; a partial current pass is unbiased
    # in position/geometry order, unlike a sequential geometry scan.
    for p,items in enumerate(passes):
        random.Random(ORDER_SEED+p).shuffle(items)
    return [r for items in passes for r in items]


def write_csv(path,header,rows):
    with path.open("w",newline="") as f:
        wr=csv.DictWriter(f,fieldnames=header)
        wr.writeheader();wr.writerows(rows)


def batches(data):
    return [data[i:i+CHUNK_SIZE] for i in range(0,len(data),CHUNK_SIZE)]


def ensure_plan(run_dir,meta,index_sha,selection,binary_sha):
    expected={
       "schema":STUDY,
       "n_original_lut500_trials":500,
       "n_cells":4200,"per_cell":PER_CELL,"total_selected":len(selection),
       "n_chunks":len(batches(selection)),
       "chunk_size":CHUNK_SIZE,"reservoir_per_avalanche":RESERVOIR,
       "selection":"SHA256 deterministic offset per geometry+box; every cell uses 20 distinct sample indices (offset+137*pass)%500",
       "order":"20 complete passes, deterministically shuffled each pass across 4200 four-geometry cells",
       "original_lut500_csv_sha256":{k:v["lut_csv_sha256"] for k,v in meta.items()},
       "frozen_432_index_csv_sha256":index_sha,
       "compiled_cpp_sha256":sha(CXX),
       "compiled_executable_sha256":binary_sha,
       "physics":"One new independent full AvalancheElectron from every exact original LUT500 start (0.1eV isotropic, ±3mm cathodes, p3/p4, FW75/noFW, +1800V, Ar:CO2 70:30, B=0). OLD drift collected flag only provenance, never gating.",
       "no_full_sensor_signal":"Only primary-electron avalanche and full ion-birth angular census; weighting fields/readout not modeled",
       "created_in_repo_commit":commit(),
    }
    frozen=run_dir/"run_config.json"
    if frozen.exists():
        old=json.loads(frozen.read_text())
        if {k:v for k,v in old.items() if k!="created_in_repo_commit"}!={k:v for k,v in expected.items() if k!="created_in_repo_commit"}:
            raise ValueError("Existing immutable overnight physics configuration differs")
    else:
        frozen.write_text(canonical(expected))
    return expected


def run_chunk(run_dir,i,selected,env):
    chunk=run_dir/"chunks"/f"chunk_{i:04d}"
    manifest_path=chunk/"complete.json"
    if chunk.exists():
        if not manifest_path.is_file():
            raise RuntimeError("Incomplete old physics chunk retained; DO NOT overwrite: "+str(chunk))
        saved=json.loads(manifest_path.read_text())
        if saved.get("status")!="complete" or saved.get("chunk")!=i:
            raise RuntimeError("Prior chunk has invalid status: "+str(chunk))
        original=chunk/"original_exact_starts.csv"
        if saved.get("selected_sha256")!=sha(original):
            raise RuntimeError("Prior chunk exact-start input changed: "+str(original))
        if (saved.get("n_original_starts")!=len(selected)
            or saved.get("source_compiled_sha256")!=sha(EXE)):
            raise RuntimeError("Prior complete chunk incompatible with current binary/plan")
        for name,digest in saved["outputs_sha256"].items():
            if sha(chunk/name)!=digest:
                raise RuntimeError("Prior completed physics chunk output altered: "+name)
        stats=validate_outputs(chunk,selected,reservoir_limit=RESERVOIR)
        # A size-capped tail is retained as CENSORED data with an explicit
        # warning, rather than abandoning an otherwise complete night.
        return stats,False

    chunk.mkdir(parents=True)
    index=chunk/"original_exact_starts.csv"
    write_csv(index,EXPORT_FIELDS,selected)
    # C++ parser reads all old physical start keys BEFORE its gas calculation;
    # pilot chunk must never include more than the source's 432 cap.
    if len(selected)>432:raise RuntimeError("C++ row cap violation")
    process=[str(EXE),"--lut-sidecar-compact5",
             str(index),str(chunk/"lut_avalanche")]
    started=time.perf_counter()
    with (chunk/"physics.log").open("w") as log:
        code=subprocess.run(process,cwd=chunk,stdout=log,
                            stderr=subprocess.STDOUT,env=env,check=False).returncode
    if code:
        raise RuntimeError("LUT500 overnight Garfield C++ FAILED in "+str(chunk)
          +"\nLast log lines:\n"
          +"\n".join((chunk/"physics.log").read_text(errors="replace").splitlines()[-15:]))
    stats=validate_outputs(chunk,selected,reservoir_limit=RESERVOIR)
    if stats["n_cap"] or stats["n_fail"]:
        print("WARNING: chunk",i,"contains",stats["n_cap"],
              "capped and",stats["n_fail"],
              "unsuccessful avalanches (retained as censored/diagnostic)",
              flush=True)
    paths=[chunk/"original_exact_starts.csv",chunk/"physics.log"]
    paths += [chunk/("lut_avalanche_"+suffix+".csv")
             for suffix in ("trials","birth_reservoir","timing")]
    state={
        "schema":"PI_LUT500_COMPACT5_SINGLE_CHUNK_V1",
        "status":"complete","chunk":i,"n_original_starts":len(selected),
        "selected_sha256":sha(index),"source_compiled_sha256":sha(EXE),
        "chunk_wall_seconds":time.perf_counter()-started,
        "stats":{k:v for k,v in stats.items() if k!="files_sha256"},
        "outputs_sha256":{p.name:sha(p) for p in paths},
    }
    manifest_path.write_text(canonical(state))
    return stats,True


def recover_existing_smoke_chunk(target, selected):
    """Explicit forensic salvage of 8 completed Garfield trials after Python
    only validation failed. Never run Garfield again or change raw data.
    Requires unique exact source identity, complete output, and C++ PASS.
    """
    chunk=target/"chunks"/"chunk_0000"
    if not chunk.is_dir() or (chunk/"complete.json").exists():
        raise RuntimeError("Expected existing, unsealed smoke chunk; no files modified")
    original=chunk/"original_exact_starts.csv"
    with original.open(newline="") as f:
        rd=csv.DictReader(f)
        if tuple(rd.fieldnames or ())!=EXPORT_FIELDS:
            raise ValueError("Failed smoke source schema mismatch")
        saved=list(rd)
    if len(saved)!=len(selected):
        raise ValueError("Recovered smoke has wrong number of original starts")
    for source,expected in zip(saved,selected):
        for key in EXPORT_FIELDS:
            if source[key]!=str(expected[key]):
                raise ValueError("Recovered smoke source start changed: "+key)
    log=chunk/"physics.log"
    if "PASS: independent full-start microscopic avalanche sidecar" not in log.read_text(errors="replace"):
        raise RuntimeError("Original Garfield C++ did not report complete sidecar run")
    stats=validate_outputs(chunk,selected,reservoir_limit=RESERVOIR)
    if stats["n_cap"] or stats["n_fail"]:
        raise RuntimeError("Refusing to seal capped or failed smoke")
    outputs=[original,log]+[chunk/("lut_avalanche_"+suffix+".csv")
                            for suffix in ("trials","birth_reservoir","timing")]
    state={
        "schema":"PI_LUT500_COMPACT5_SINGLE_CHUNK_V1",
        "status":"complete","chunk":0,"n_original_starts":len(selected),
        "selected_sha256":sha(original),
        "source_compiled_sha256":sha(EXE),
        "chunk_wall_seconds":None,
        "recovery_note":"Forensic seal after Python int-vs-string identity comparison bug; EXACT original raw Garfield outputs preserved",
        "stats":{k:v for k,v in stats.items() if k!="files_sha256"},
        "outputs_sha256":{p.name:sha(p) for p in outputs},
    }
    (chunk/"complete.json").write_text(canonical(state))
    print("PASS recovered",len(selected),
          "ORIGINAL microscopic Garfield avalanche rows "
          "without rerun/overwrite; SHA-protected smoke checkpoint",flush=True)
    return stats


def aggregate(run_dir,plan,n_completed,elapsed):
    """Publishable, per-geometry full responses + per-cell moments.

    All raw checkpoint chunks are kept in run_dir/chunks. No approximation
    changes the FULL ion births per avalanche; only 5 randomly selected
    birth coordinates per avalanche are stored.
    """
    geoms={g:{"n":0,"ion":0,"zero":0,"time":0.0,
              "up":0,"down":0,"top":0,"bottom":0,
              "capped":0,"failed":0} for g in GEOM_NAMES}
    cells=defaultdict(lambda:{"births":[],"up":0,"down":0,"top":0,
                              "bottom":0,"zero":0,"time":0.0})
    record_path={}
    chosen_headers={}
    try:
        for geom in GEOM_NAMES:
            for category,cols in (("trials",COLUMNS),("timing",TIME_COLUMNS)):
                f=(run_dir/f"{category}_{geom}.csv").open("w",newline="")
                writer=csv.DictWriter(f,fieldnames=cols)
                writer.writeheader();record_path[(category,geom)]=(f,writer)
            for part in range(2):
                f=(run_dir/f"birth_reservoir_{geom}_part{part}.csv").open("w",newline="")
                writer=csv.DictWriter(f,fieldnames=BIRTH_COLUMNS)
                writer.writeheader()
                record_path[("births",geom,part)]=(f,writer)
        ids=set()
        for i in range(n_completed):
            base=run_dir/"chunks"/f"chunk_{i:04d}"
            with (base/"lut_avalanche_trials.csv").open(newline="") as f:
                rd=csv.DictReader(f)
                if tuple(rd.fieldnames or ())!=COLUMNS:raise ValueError("Wrong aggregated trial schema")
                for row in rd:
                    geom=row["configuration"]
                    if geom not in geoms:raise ValueError("Unexpected geometry")
                    identity=(geom,int(row["box_id"]),int(row["sample"]))
                    if identity in ids:raise ValueError("Duplicate original LUT500 trial")
                    ids.add(identity)
                    record_path[("trials",geom)][1].writerow(row)
                    ion=int(row["ion_births"]);s=cells[(geom,int(row["box_id"]))]
                    if ion<0:raise ValueError("Negative ion birth")
                    s["births"].append(ion)
                    s["up"]+=int(row["ion_up"]);s["down"]+=int(row["ion_down"])
                    s["top"]+=int(row["ion_top"]);s["bottom"]+=int(row["ion_bottom"])
                    s["zero"]+=int(ion==0)
                    g=geoms[geom];g["n"]+=1;g["ion"]+=ion;g["zero"]+=int(ion==0)
                    g["capped"]+=int(row["capped"])
                    g["failed"]+=int(row["successful"])!=1
                    g["up"]+=int(row["ion_up"]);g["down"]+=int(row["ion_down"])
                    g["top"]+=int(row["ion_top"]);g["bottom"]+=int(row["ion_bottom"])
            with (base/"lut_avalanche_timing.csv").open(newline="") as f:
                rd=csv.DictReader(f)
                if tuple(rd.fieldnames or ())!=TIME_COLUMNS:raise ValueError("Wrong timing schema")
                for row in rd:
                    geom=row["configuration"]
                    dur=float(row["wall_seconds"])
                    if not math.isfinite(dur) or dur<0:raise ValueError("Bad timing")
                    record_path[("timing",geom)][1].writerow(row)
                    geoms[geom]["time"]+=dur
                    box=int(row["crossing_id"].split("|")[2])
                    cells[(geom,box)]["time"]+=dur
            with (base/"lut_avalanche_birth_reservoir.csv").open(newline="") as f:
                rd=csv.DictReader(f)
                if tuple(rd.fieldnames or ())!=BIRTH_COLUMNS:raise ValueError("Wrong birth schema")
                for row in rd:
                    geom=row["configuration"]
                    box=int(row["crossing_id"].split("|")[2])
                    record_path[("births",geom,box%2)][1].writerow(row)
    finally:
        for f,_ in record_path.values():f.close()
    if sum(g["n"] for g in geoms.values())!=n_completed*CHUNK_SIZE and n_completed!=plan["n_chunks"]:
        raise RuntimeError("Completed chunk source count mismatch")
    if sum(g["n"] for g in geoms.values())>84000:
        raise RuntimeError("More than 84k output avalanches")
    fields=("geometry","box_id","iu","iv","n","zero","p_zero",
        "mean_ion_births","sd_ion_births","median_ion_births",
        "mean_timing_seconds","ion_weighted_A_hemi","ion_weighted_A_TB")
    with (run_dir/"cell_statistics.csv").open("w",newline="") as f:
        wr=csv.DictWriter(f,fieldnames=fields);wr.writeheader()
        for geom,pitch,_ in CONFIGS:
            for b in range(round(pitch/.1)*30):
                s=cells.get((geom,b))
                n=len(s["births"]) if s else 0
                hemi=s["up"]+s["down"] if s else 0
                tb=s["top"]+s["bottom"] if s else 0
                wr.writerow({
                    "geometry":geom,"box_id":b,"iu":b//30,"iv":b%30,"n":n,
                    "zero":s["zero"] if s else 0,
                    "p_zero":(s["zero"]/n if n else ""),
                    "mean_ion_births":(statistics.mean(s["births"]) if n else ""),
                    "sd_ion_births":(statistics.pstdev(s["births"]) if n else ""),
                    "median_ion_births":(statistics.median(s["births"]) if n else ""),
                    "mean_timing_seconds":(s["time"]/n if n else ""),
                    "ion_weighted_A_hemi":((s["up"]-s["down"])/hemi if hemi else ""),
                    "ion_weighted_A_TB":((s["top"]-s["bottom"])/tb if tb else ""),
                })
    completed=sum(g["n"] for g in geoms.values())
    summary={
        "schema":STUDY,"status":"complete" if completed==84000 else "time_budget_partial",
        "n_original_frozen_lut500_cells":4200,"original_trials_available_per_cell":500,
        "planned_fresh_avalanches_per_cell":PER_CELL,
        "n_completed_avalanches":completed,"n_planned_avalanches":84000,
        "n_completed_chunks":n_completed,"n_planned_chunks":plan["n_chunks"],
        "elapsed_seconds_this_invocation":elapsed,
        "full_ion_births":sum(x["ion"] for x in geoms.values()),
        "total_capped":sum(x["capped"] for x in geoms.values()),
        "total_failed":sum(x["failed"] for x in geoms.values()),
        "birth_reservoir_limit_per_avalanche":RESERVOIR,
        "geometries":{},
        "note":"NO native Geant4 muon physics and NO induced cathode signals. Every trial is an independent full-start microscopic Garfield avalanche from the exact frozen LUT500 initial electron position; old collected is never used as survival gate.",
        "source_cpp_sha256":plan["compiled_cpp_sha256"],
        "source_executable_sha256":plan["compiled_executable_sha256"],
        "old_lut500_csv_sha256":plan["original_lut500_csv_sha256"],
    }
    if summary["total_capped"] or summary["total_failed"]:
        summary["physics_warning"]="Capped trials are right-censored in gain; failed trials are diagnostics. Do not infer unbiased gain tails from them."
    for geom,g in geoms.items():
        hemi=g["up"]+g["down"]
        tb=g["top"]+g["bottom"]
        summary["geometries"][geom]={
            "n":g["n"],"zero":g["zero"],"capped":g["capped"],"failed":g["failed"],
            "p_zero":(g["zero"]/g["n"] if g["n"] else None),
            "total_ion_births":g["ion"],
            "mean_ion_births_per_start":g["ion"]/g["n"] if g["n"] else None,
            "mean_avalanche_seconds":g["time"]/g["n"] if g["n"] else None,
            "ion_weighted_A_hemi":(g["up"]-g["down"])/hemi if hemi else None,
            "ion_weighted_A_TB":(g["top"]-g["bottom"])/tb if tb else None,
        }
    (run_dir/"summary.json").write_text(canonical(summary))
    outputs={p.name:sha(p) for p in sorted(run_dir.iterdir())
             if p.is_file() and p.name not in ("manifest.json",)}
    summary["files_sha256"]=outputs
    (run_dir/"manifest.json").write_text(canonical(summary))
    # Guard publication against silently skipping oversize individual files
    # or an oversize total, per the EXISTING publication script's policy.
    from pi_publish_run_to_drive import MAX_FILE_BYTES, MAX_TOTAL_BYTES
    allowed=[x for x in run_dir.iterdir() if x.is_file()]
    files_for_upload=[p for p in allowed if p.name!="manifest.json"]
    over=[(p.name,p.stat().st_size) for p in files_for_upload if p.stat().st_size>MAX_FILE_BYTES]
    total=sum(p.stat().st_size for p in files_for_upload)
    if over or total>MAX_TOTAL_BYTES:
        raise RuntimeError(
            f"Run finished, raw outputs preserved, but Drive bundle limits exceeded: "
            f"oversized={over}, total={total} / {MAX_TOTAL_BYTES}. "
            "Do not discard raw data; publish reduced/split copy later.")
    return summary


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("mode",choices=("smoke","recover-smoke","overnight"))
    p.add_argument("--hours",type=float,default=8.75,
                   help="Wall-time budget incl setup/aggregation, default 8.75h")
    p.add_argument("--resume",action="store_true",
                   help="Explicitly continue a previously interrupted or time-budget-partial overnight directory")
    p.add_argument("--upload",action="store_true",
                   help="After final validation publish immutable full per-geometry CSV to existing na60drive: (rclone)")
    a=p.parse_args()
    if not math.isfinite(a.hours) or not .1<=a.hours<=10:
        p.error("--hours must be between 0.1 and 10")
    started=time.monotonic()
    binary=check_executable()
    env=garfield_runtime_env()
    print(f"PASS compiled B1 Garfield source SHA and Heed database {env['HEED_DATABASE']}",flush=True)
    meta,ixsha=source_maps()
    print("PASS existing 8 original LUT20/500 immutable source SHA checks",flush=True)
    starts=generate_selected(meta)
    # Small isolated first test of the compact 5-ion reservoir format.
    if a.mode in ("smoke","recover-smoke"):
        target=ROOT.parent/"smoke8_compact5_v1"
        if a.mode=="smoke":
            if target.exists():
                raise FileExistsError(
                    "Preserving prior smoke; if Python validation failed AFTER "
                    "Garfield completed, use 'recover-smoke' to verify it "
                    "without repeating the physics: "+str(target))
            target.mkdir(parents=True)
        elif not target.is_dir():
            raise FileNotFoundError("No failed smoke folder to recover: "+str(target))
        subset=[]
        # TWO distinct real LUT500 starts per geometry, all four fields.
        # No Heed is run in this Garfield-only frozen-grid test.
        for pass_id in (0,1):
            segment=starts[pass_id*4200:(pass_id+1)*4200]
            for geom in GEOM_NAMES:
                candidate=next((x for x in segment if x["geometry"]==geom),
                               None)
                if candidate is None:
                    raise RuntimeError("Smoke input missing geometry "+geom)
                subset.append(candidate)
        if a.mode=="recover-smoke":
            stats=recover_existing_smoke_chunk(target,subset)
        else:
            stats,_=run_chunk(target,0,subset,env)
        # Crucial physics-equivalence preflight: the compact output samples
        # only 5 instead of 100 ion birth positions. It MUST reproduce
        # every FULL microscopic avalanche's seed, gain and ion census.
        # Reservoir sampling has its own RNG and must not alter Garfield.
        reference=target/"reference_full100"
        if reference.exists():
            raise FileExistsError(
                "Reference run already exists; preserving original reference files: "
                +str(reference))
        reference.mkdir()
        input_csv=target/"chunks"/"chunk_0000"/"original_exact_starts.csv"
        shutil.copy2(input_csv,reference/"original_exact_starts.csv")
        with (reference/"physics.log").open("w") as log:
            status=subprocess.run([str(EXE),"--lut-sidecar",
                    str(reference/"original_exact_starts.csv"),
                    str(reference/"lut_avalanche")],
                    cwd=reference,stdout=log,stderr=subprocess.STDOUT,
                    check=False,env=env).returncode
        if status:
            raise RuntimeError("Compact-output equivalence full100 reference failed: "
                               +str(reference/"physics.log"))
        ref_stats=validate_outputs(reference,subset,reservoir_limit=100)
        def trial_rows(file):
            with file.open(newline="") as f:
                rd=csv.DictReader(f)
                if tuple(rd.fieldnames or ())!=COLUMNS:
                    raise ValueError("Wrong reference trial schema")
                return list(rd)
        compact_rows=trial_rows(target/"chunks"/"chunk_0000"/"lut_avalanche_trials.csv")
        full_rows=trial_rows(reference/"lut_avalanche_trials.csv")
        if compact_rows!=full_rows:
            raise RuntimeError("CRITICAL: compact5 changed microscopic avalanche results! "
                               "Long run BLOCKED.")
        verified={
            "schema":"PI_OVERNIGHT_COMPACT_ION_SAMPLING_EQUIVALENCE_V1",
            "status":"complete","n_trials":len(subset),
            "full100_trials_sha256":sha(reference/"lut_avalanche_trials.csv"),
            "compact5_trials_sha256":sha(target/"chunks"/"chunk_0000"/"lut_avalanche_trials.csv"),
            "full_ion_births":stats["total_ion_births"],
            "reference_ion_births":ref_stats["total_ion_births"],
            "physics_identical":True,
        }
        (target/"smoke_verified.json").write_text(canonical(verified))
        print("PASS compact5 smoke AND exact 100-vs-5 ion reservoir "
              "microscopic physics identity:",stats["n_trials"],
              "original LUT500 starts; full ion births:",
              stats["total_ion_births"],flush=True)
        return
    smoke_root=ROOT.parent/"smoke8_compact5_v1"
    smoke=smoke_root/"smoke_verified.json"
    if not smoke.is_file():
        raise RuntimeError("First PASS full100-vs-compact5 smoke; no long job started")
    sm=json.loads(smoke.read_text())
    if (sm.get("physics_identical") is not True or
        sm.get("status")!="complete" or
        sm.get("compact5_trials_sha256")!=sha(smoke_root/"chunks"/"chunk_0000"/"lut_avalanche_trials.csv") or
        sm.get("full100_trials_sha256")!=sha(smoke_root/"reference_full100"/"lut_avalanche_trials.csv")):
        raise RuntimeError("Full100-vs-compact5 physics preflight provenance changed")
    if shutil.disk_usage(ROOT.parent).free<1_000_000_000:
        raise RuntimeError("Less than 1GB available for immutable checkpoint outputs")
    if ROOT.exists() and not a.resume:
        raise FileExistsError(
            "Overnight output directory already exists; preserving old chunks. "
            "Pass --resume only after reviewing its current manifest: "+str(ROOT))
    if (ROOT/"manifest.json").is_file():
        current=json.loads((ROOT/"manifest.json").read_text())
        if current.get("status")=="complete":
            raise FileExistsError("Overnight 84k run already complete; never recalculate: "+str(ROOT))
    ROOT.mkdir(parents=True,exist_ok=True)
    plan=ensure_plan(ROOT,meta,ixsha,starts,binary)
    all_chunks=batches(starts)
    print(f"OVERNIGHT: {len(starts)} new Garfield avalanches, {len(all_chunks)}"
          f" immutable chunks of {CHUNK_SIZE}, deadline {a.hours:.2f}h",flush=True)
    completed=0
    for i,subset in enumerate(all_chunks):
        # Reserve 5 min for aggregation, checksum audits and optional upload.
        if time.monotonic()-started >= a.hours*3600-300:
            print("SAFE STOP: wall-time budget reached, preserving all finished chunks",flush=True)
            break
        stats,was_new=run_chunk(ROOT,i,subset,env)
        completed+=1
        elapsed=time.monotonic()-started
        if (completed%5==0 or completed==len(all_chunks)):
            print(f"CHECKPOINT {completed}/{len(all_chunks)}, {completed*CHUNK_SIZE}"
                  f" real original LUT500 starts, elapsed={elapsed/3600:.2f}h,"
                  f" current chunk birth total={stats['total_ion_births']},"
                  f" {'new' if was_new else 'verified existing'}",
                  flush=True)
    if not completed:
        raise RuntimeError("No complete chunk to aggregate; source/case preserved")
    result=aggregate(ROOT,plan,completed,time.monotonic()-started)
    print("FINAL",result["status"],result["n_completed_avalanches"],
          "/84000 microscopic avalanche trials; raw chunks preserved",flush=True)
    for geom,s in result["geometries"].items():
        print(geom,s,flush=True)
    if a.upload:
        from pi_publish_run_to_drive import prepare,publish,EXPECTED_REMOTE
        run_id=f"lut500_overnight_20261029_c{completed:03d}_v1"
        bundle=prepare(ROOT,run_id)
        publish(bundle,EXPECTED_REMOTE,upload=True)
        print("PASS uploaded and independently checked full output on Google Drive",flush=True)


if __name__=="__main__":
    main()
