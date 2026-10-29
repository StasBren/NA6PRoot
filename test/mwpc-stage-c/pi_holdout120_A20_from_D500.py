#!/usr/bin/env python3
"""INDEPENDENT Geant4 chunk-2 muon holdout to validate A20@D500.

A20@D500 = 20 NEW independent full microscopic Garfield avalanches per
0.1mm x 0.2mm cell, sampled from 500 OLD original DRIFT starting points.
It is NOT an avalanche LUT with 500 trials per cell.

120 distinct actual J/psi Geant4 muon tracks are SHA-ranked/stratified
from frozen ORIGINAL chunk2, disjoint from the 20-muon chunk1 pilot.
Each true TrackHeed source is generated ONCE per muon, and full
Garfield is run from its exact Heed primary electrons in all four fields.
The resulting exact SAME Heed electron inputs are subsequently sampled
through A20@D500 200 times. No Geant4 particle/angle is synthesized.

The microscopic job checkpoints 10 muons per chunk, is resumable without
overwriting existing data and has an optional between-chunk time budget.
All original LUT500, Garfield results and 20-muon pilot are read-only.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import random
import shutil
import statistics
import subprocess
import tempfile
import time
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path

from pi_2m_native_six_geometry_fold import crossing_provenance, SOURCE
from pi_prepare_muon_avalanche_b1 import FIELDS,validate,sha,seed
from pi_run_native_muon_b1 import (
    CXX, EXE, GEOMS, CROSSING_COLUMNS, check_executable,
    validate_native_csv_before_physics, commit,
)
from pi_garfield_runtime_env import garfield_runtime_env
from pi_run_muon_lut_validation import (
    N_U,PRIMARY_FIELDS,CLUSTER_FIELDS,SEED_FIELDS,
    load_lut,spatial_cell,stable_seed,summarize_samples,asym,
    validate_reference,read_csv,write_json,id_of,percent,round_away,
)
from pi_publish_run_to_drive import prepare as drive_prepare,publish,EXPECTED_REMOTE
from pi_audit_avalanche_lut_reuse import BASE

ROOT=BASE/"pi_muon_holdout120_A20_from_D500_v1"
INPUT=ROOT/"input_chunk2_120"
REF=ROOT/"full_Garfield_120"
COMP=ROOT/"comparison_A20_from_D500_120"
N_MUONS=120
GROUP=10
N_CHUNKS=N_MUONS//GROUP
REPLICATES=200
SOURCE_CHUNK=2
FIELDS_TO_MERGE={
    "muon_reference_crossing.csv":CROSSING_COLUMNS,
    "muon_reference_seed.csv":SEED_FIELDS,
    "muon_reference_heed_primaries.csv":PRIMARY_FIELDS,
    "muon_reference_heed_clusters.csv":CLUSTER_FIELDS,
}
# Birth samples retained in chunk originals and aggregated for archival
# (could be large but current birth reservoir is at most five per electron).
BIRTH_FIELDS=("crossing_id","configuration","primary_index","u_mm",
              "v_mm","w_mm","t_ns","seed_ion_births_total","reservoir_weight")
FIELDS_TO_MERGE["muon_reference_ion_birth_sample.csv"]=BIRTH_FIELDS
KEY=CROSSING_COLUMNS[0]

def time_now():
    return datetime.now(timezone.utc).isoformat()

def write_csv(path,fields,rows):
    with path.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=fields)
        w.writeheader();w.writerows(rows)

def payload_rows():
    path=INPUT/"native_muon_crossings_input.csv"
    metaf=INPUT/"input_manifest.json"
    if not path.is_file() or not metaf.is_file():
        raise FileNotFoundError("First run prepare to select 120 original chunk-2 muons")
    m=json.loads(metaf.read_text())
    if (m.get("schema")!="PI_INDEPENDENT_GEANT4_MUON_HOLDOUT120_V1"
        or m.get("status")!="complete"
        or m.get("n_muons")!=120 or m.get("source_chunk")!=SOURCE_CHUNK
        or m.get("sha256_input")!=sha(path)):
        raise ValueError("Original 120-muon input manifest/SHA mismatch")
    rows=read_csv(path,FIELDS)
    if len(rows)!=120 or len({id_of(r) for r in rows})!=120:
        raise ValueError("Wrong number of source gas crossings")
    unique_tracks={(r["source"],r["chunk"],r["event"],r["track_id"])
                   for r in rows}
    if len(unique_tracks)!=120:
        raise ValueError("One physical Geant4 track selected more than once")
    for r in rows:
        if r["source"]!=SOURCE["Jpsi"] or str(r["chunk"])!=str(SOURCE_CHUNK):
            raise ValueError("Non-holdout Geant4 source in fixed holdout input")
        if validate(r) is None:raise ValueError("Geant4 original hit now invalid")
    return rows,m

def choose_holdout_from_pool(pool,per_stratum=10):
    """Pick original distinct GEANT4 physical muons, balanced across
    all 12 station/signed-angle strata, preserving original SHA ranking.

    PURE selection helper; unit tested without Geant4/Garfield.
    """
    expected={(st,side) for st in range(6) for side in (-1,1)}
    if not expected.issubset(pool):
        raise ValueError("Missing signed-angle/station Geant4 strata")
    selected=[];tracks=set();counts=Counter()
    for k in range(per_stratum):
        for mirror in (False,True):
            for st in range(6):
                sign=(1 if st%2==0 else -1)*(-1 if mirror else 1)
                stratum=(st,sign)
                candidate=None
                for _,_,rec in pool[stratum]:
                    physical=(str(rec["source"]),str(rec["chunk"]),
                              str(rec["event"]),str(rec["track_id"]))
                    if physical not in tracks:
                        candidate=(physical,rec)
                        break
                if candidate is None:
                    raise ValueError("No independent original Geant4 track for "+str(stratum))
                physical,rec=candidate
                tracks.add(physical);selected.append(rec)
                counts[stratum]+=1
    if len(selected)!=len(expected)*per_stratum or len(tracks)!=len(selected) or any(
            counts[stratum]!=per_stratum for stratum in expected):
        raise AssertionError("Unbalanced/duplicated physical muons")
    return selected,counts


def prepare():
    if ROOT.exists():
        raise FileExistsError("Preserve prior holdout directory; use status, not overwrite: "+str(ROOT))
    path,provenance=crossing_provenance("Jpsi",SOURCE_CHUNK)
    source=SOURCE["Jpsi"]
    # Keep 48 best original gas crossings per signed-angle/station stratum.
    # Bounded memory regardless of input file's number of rows.
    pool=defaultdict(list);counts=Counter()
    with path.open(newline="") as f:
        reader=csv.DictReader(f)
        if not set(FIELDS)-{"source","chunk","momentum_GeV","heed_seed"}<=set(reader.fieldnames or ()):
            raise ValueError("Frozen Geant4 chunk2 missing original native input fields")
        for r in reader:
            data=validate(r)
            if data is None:continue
            moment,_=data
            st=int(r["station"])
            sign=-1 if float(r["theta_u_in_deg"])<0 else 1
            key=(st,sign)
            counts[key]+=1
            srcid=(source,SOURCE_CHUNK,r["event"],r["hit_index"],
                   r["track_id"],st,r["det_id"])
            score=seed("A20D500-120-independent-holdout-v1",*srcid)
            rec={k:r[k] for k in FIELDS if k not in
                 ("source","chunk","momentum_GeV","heed_seed")}
            rec["source"]=source;rec["chunk"]=SOURCE_CHUNK
            rec["momentum_GeV"]=f"{moment:.12g}"
            rec["heed_seed"]=seed("real-Heed-once-v1",*srcid)
            q=pool[key]
            q.append((score,tuple(map(str,srcid)),rec))
            q.sort(key=lambda v:(v[0],v[1]))
            if len(q)>48:q.pop()
    selected,by_stratum=choose_holdout_from_pool(pool,per_stratum=10)
    expected={(st,side) for st in range(6) for side in (-1,1)}
    if len(selected)!=N_MUONS:
        raise AssertionError("Not all 120 independent physical Geant4 muons selected")
    INPUT.mkdir(parents=True)
    outfile=INPUT/"native_muon_crossings_input.csv"
    write_csv(outfile,FIELDS,selected)
    manifest={
        "schema":"PI_INDEPENDENT_GEANT4_MUON_HOLDOUT120_V1",
        "status":"complete","created_utc":time_now(),
        "source_channel":"Jpsi","source_chunk":SOURCE_CHUNK,
        "source_id":source,"n_muons":120,"n_independent_physical_tracks":120,
        "per_station_signed_angle_stratum":10,
        "selection":"SHA-ranked 48 best per signed-angle×station stratum, 10 balanced passes; no duplicate physical Geant4 track, physical incidence-angle PDF NOT represented",
        "sha256_input":sha(outfile),
        "frozen_source_csv":str(path),
        "frozen_source_sha256":provenance["csv_sha256"],
        "original_root_sha256":provenance["immutable_input"]["input_root_sha256"],
        "original_parent_count":provenance["immutable_input"]["n_generated_parent_events"],
        "original_lut_model":"A20@D500: 20 independent full Garfield avalanches drawn from 500 original DRIFT starts per cell, 4200 cells/84000 new avalanches",
        "selected_original_ids":[id_of(r) for r in selected],
        "stratum_counts":{f"MS{s}_{'negative' if sign<0 else 'nonnegative'}":v
                          for (s,sign),v in sorted(by_stratum.items())},
        "station_input_counts":{f"MS{s}_{'negative' if sign<0 else 'nonnegative'}":counts[s,sign]
                          for s,sign in sorted(expected)},
        "git_commit":commit(),
    }
    write_json(INPUT/"input_manifest.json",manifest)
    print("PASS prepared 120 DIFFERENT real Geant4 Jpsi muon tracks "
          "from UNSEEN ORIGINAL chunk2; 10 per signed-angle×station stratum",flush=True)
    print("INPUT",outfile,flush=True)

def groups(rows):
    if len(rows)!=120:raise ValueError("Requires complete fixed holdout")
    return [rows[i*GROUP:(i+1)*GROUP] for i in range(N_CHUNKS)]

def preflight():
    rows,manifest=payload_rows()
    exe=check_executable()
    print("PASS full B1 C++ binary/source SHA",exe,flush=True)
    # Validate every true Geant4 chord and Heed cluster/electron generation
    # in the SAME compiled binary which will later run the real avalanche.
    with tempfile.TemporaryDirectory() as tmp:
        for i,subset in enumerate(groups(rows)):
            p=Path(tmp)/f"holdout{i:02d}.csv"
            write_csv(p,FIELDS,subset)
            validate_native_csv_before_physics(p)
            print(f"PASS true Geant4 + TrackHeed preflight {i+1}/{N_CHUNKS}",flush=True)
    print("PASS all 120 native Geant4 muons produce correct Heed "
          "physics in compiled binary (NO full avalanches yet)",flush=True)

def chunk_dir(i):
    return REF/"chunks"/f"chunk_{i:02d}"

def read_seal(c,i,subset,binary):
    seal=c/"complete.json"
    if not seal.is_file():raise RuntimeError("Unsealed existing full Garfield chunk: "+str(c))
    obj=json.loads(seal.read_text())
    if (obj.get("status")!="complete" or obj.get("chunk")!=i or
        obj.get("cpp_executable_sha256")!=binary or
        obj.get("n_muons")!=len(subset)):
        raise ValueError("Full Garfield immutable checkpoint identity mismatch: "+str(c))
    if obj.get("input_sha256")!=sha(c/"native_muon_crossings_input.csv"):
        raise ValueError("Frozen Garfield chunk input SHA changed")
    if [id_of(x) for x in read_csv(c/"native_muon_crossings_input.csv",FIELDS)]!=[id_of(x) for x in subset]:
        raise ValueError("Different real Geant4 muons in saved checkpoint")
    for name,digest in obj["outputs_sha256"].items():
        if sha(c/name)!=digest:
            raise ValueError(f"Full Garfield original file changed: {c/name}")
    return obj

def run(hours, max_chunks):
    if hours<=0 or max_chunks<1 or max_chunks>N_CHUNKS:
        raise ValueError("Use --hours >0 and --max-chunks 1..12")
    rows,meta=payload_rows()
    binary=check_executable()
    if REF.exists() and not (REF/"run_config.json").is_file():
        raise RuntimeError("Unrecognized/incomplete full Garfield root; preserve it for forensics")
    if not REF.exists():
        REF.mkdir(parents=True)
        write_json(REF/"run_config.json",{
            "schema":"PI_INDEPENDENT_GEANT4_HEED_GARFIELD_HOLDOUT120_V1",
            "status":"running","created_utc":time_now(),
            "input_sha256":meta["sha256_input"],
            "source_geant4_sha256":meta["frozen_source_sha256"],
            "cpp_source_sha256":sha(CXX),
            "cpp_executable_sha256":binary,
            "description":"120 unseen physical Geant4 chunk2 muons, 12 chunks×10 true Heed/4 independent full Garfield geometries; A20@D500 validation; no LUT altered",
        })
    cfg=json.loads((REF/"run_config.json").read_text())
    if (cfg.get("input_sha256")!=meta["sha256_input"]
        or cfg.get("cpp_source_sha256")!=sha(CXX)
        or cfg.get("cpp_executable_sha256")!=binary):
        raise ValueError("Rebuild/source/input changed since original holdout Garfield start")
    # A subsequent safe rerun must only REVERIFY the previously sealed
    # calculation. Never attempt a second aggregate or overwrite it.
    if (REF/"manifest.json").is_file():
        verify_full()
        print("PASS existing full 120-muon Garfield reference already SHA-sealed; "
              "nothing recomputed or overwritten",flush=True)
        return
    start=time.monotonic()
    done=0
    for i,subset in enumerate(groups(rows)):
        c=chunk_dir(i)
        if c.exists():
            read_seal(c,i,subset,binary)
            done+=1
            print(f"PASS reused ORIGINAL SHA-valid Garfield checkpoint {i+1}/12",flush=True)
            continue
        if done>=max_chunks or time.monotonic()-start>=hours*3600:
            print(f"SAFE STOP: {done}/12 original 10-muon chunks complete; "
                  "restart same run mode to resume without recomputing",flush=True)
            return
        c.mkdir(parents=True)
        write_csv(c/"native_muon_crossings_input.csv",FIELDS,subset)
        prefix=c/"muon_reference"
        print(f"REAL FULL Garfield chunk {i+1}/12 started, "
              f"10 unseen physical muons, 4 geometries",flush=True)
        with (c/"physics.log").open("w") as out:
            cp=subprocess.run(
                [str(EXE),"--muon-cache-physics",
                 str(c/"native_muon_crossings_input.csv"),str(prefix)],
                stdout=out,stderr=subprocess.STDOUT,
                cwd=c,check=False,env=garfield_runtime_env())
        if cp.returncode or "PASS: immutable same-Heed primary/cluster cache" not in (
                c/"physics.log").read_text(errors="replace"):
            raise RuntimeError(f"True Garfield chunk FAILED. Raw files preserved: {c}; tail:\n"+
                "\n".join((c/"physics.log").read_text(errors="replace").splitlines()[-40:]))
        stat=validate_reference(c,subset)
        if stat["n_native_muons"]!=GROUP or stat["n_heed_primaries"]<=0:
            raise ValueError("Full Garfield chunk incomplete/not physical")
        raw=[c/"native_muon_crossings_input.csv",c/"physics.log"]+[
            c/name for name in FIELDS_TO_MERGE]
        write_json(c/"complete.json",{
            "schema":"PI_10_REAL_GEANT4_MUON_GARFIELD_SHA_CHECKPOINT_V1",
            "status":"complete","chunk":i,"n_muons":len(subset),
            "cpp_executable_sha256":binary,
            "input_sha256":sha(raw[0]),
            "stats":stat,
            "outputs_sha256":{p.name:sha(p) for p in raw},
        })
        done+=1
        print(f"PASS CHECKPOINT {done}/12: "
              f"native Heed clusters={stat['n_heed_clusters']}, "
              f"full Garfield primary×geometry={stat['n_full_avalanche_trials']}; "
              "zero cap/fail, SHA sealed",flush=True)
    finish()

def finish():
    rows,meta=payload_rows()
    binary=check_executable()
    if (REF/"manifest.json").exists():
        raise FileExistsError("Full 120-muon immutable reference already sealed")
    for i,subset in enumerate(groups(rows)):
        read_seal(chunk_dir(i),i,subset,binary)
    # NEVER delete or modify original chunk files. Build separate root
    # aggregate, only once, enforcing complete CSV row schemas.
    if any((REF/name).exists() for name in FIELDS_TO_MERGE):
        raise FileExistsError("Incomplete existing aggregate; inspect without overwriting")
    for name,fields in FIELDS_TO_MERGE.items():
        with (REF/name).open("w",newline="") as out:
            w=csv.DictWriter(out,fieldnames=fields)
            w.writeheader()
            for i in range(N_CHUNKS):
                c=chunk_dir(i)
                with (c/name).open(newline="") as f:
                    reader=csv.DictReader(f)
                    if tuple(reader.fieldnames or ())!=fields:
                        raise ValueError("Garfield source chunk CSV schema mismatch: "+str(c/name))
                    for r in reader:w.writerow(r)
    write_csv(REF/"native_muon_crossings_input.csv",FIELDS,rows)
    # This EXACT old B2 physical validator checks Heed electron and cluster
    # IDs and full Garfield muon ion birth conservation across 480 rows.
    stat=validate_reference(REF,rows)
    if (stat["n_native_muons"]!=120 or
        stat["n_full_avalanche_trials"]!=4*stat["n_heed_primaries"]):
        raise ValueError("Aggregated full Garfield/Heed sample count mismatch")
    summary={
        "schema":"PI_INDEPENDENT_GEANT4_HEED_GARFIELD_HOLDOUT120_V1",
        "status":"complete","n_physically_distinct_muons":120,
        "source_chunk":SOURCE_CHUNK,"source_geant4_sha256":meta["frozen_source_sha256"],
        "original_root_sha256":meta["original_root_sha256"],
        "cpp_source_sha256":sha(CXX),"cpp_executable_sha256":binary,
        "original_geant_input_sha256":meta["sha256_input"],
        "n_completed_chunks":12,
        "A20D500_definition":"20 NEW full avalanche realizations/cell from 500 ORIGINAL DRIFT starts/cell; 4200 cells, 84000 microscopic avalanche trials",
        "stat":stat,"gas":"Ar/CO2 70/30, 293.15K 760 Torr",
        "geometry":"3mm/4mm anode spacing; noFW/FW75; anode +1800V; cathodes ±3mm; B=0",
        "note":"Real Geant4->TrackHeed clusters and 120 real muons. ONE independent full Garfield avalanche for each Heed primary per geometry. No cathode-induced strip signal. Selection is balanced station/sign not true incidence PDF.",
    }
    write_json(REF/"summary.json",summary)
    files={p.name:sha(p) for p in REF.iterdir() if p.is_file() and p.name!="manifest.json"}
    write_json(REF/"manifest.json",{
        **summary,"outputs_sha256":files,
        "original_chunks_sha256":{
            f"chunk_{i:02d}":sha(chunk_dir(i)/"complete.json") for i in range(N_CHUNKS)
        },
    })
    print("PASS FULL 120 REAL INDEPENDENT Geant4 muons x four Garfield "
          "geometries complete and SHA-sealed. "
          f"Real Heed clusters={stat['n_heed_clusters']}, "
          f"primaries={stat['n_heed_primaries']}.",flush=True)

def verify_full():
    rows,meta=payload_rows()
    manifest=REF/"manifest.json"
    if not manifest.is_file():
        raise FileNotFoundError("Full 120-muon Garfield reference unfinished. Run reference.")
    m=json.loads(manifest.read_text())
    if m.get("status")!="complete" or m.get("n_physically_distinct_muons")!=120:
        raise ValueError("Not a complete 120-muon reference")
    for name,digest in m.get("outputs_sha256",{}).items():
        if sha(REF/name)!=digest:
            raise ValueError("Original Garfield aggregated SHA changed: "+name)
    for i,subset in enumerate(groups(rows)):
        read_seal(chunk_dir(i),i,subset,check_executable())
        if sha(chunk_dir(i)/"complete.json")!=m["original_chunks_sha256"][f"chunk_{i:02d}"]:
            raise ValueError("Original full Garfield chunk checkpoint modified")
    stat=validate_reference(REF,rows)
    return rows,stat,sha(manifest)

def evaluate(replicates):
    if not 20<=replicates<=1000:raise ValueError("Replicates 20..1000")
    rows,ref_stat,refsha=verify_full()
    lut,lutsha=load_lut() # verifies ALL frozen 84k microavalanches per cell and hashes
    if COMP.exists():raise FileExistsError("Preserve earlier comparison: "+str(COMP))
    COMP.mkdir(parents=True)
    primary=defaultdict(list)
    for r in read_csv(REF/"muon_reference_heed_primaries.csv",PRIMARY_FIELDS):
        primary[r["crossing_id"]].append(r)
    clusters=defaultdict(list)
    for r in read_csv(REF/"muon_reference_heed_clusters.csv",CLUSTER_FIELDS):
        clusters[r["crossing_id"]].append(r)
    full_seeds=defaultdict(dict)
    for r in read_csv(REF/"muon_reference_seed.csv",SEED_FIELDS):
        key=(r["crossing_id"],r["configuration"])
        i=int(r["primary_index"])
        if i in full_seeds[key]:raise ValueError("Duplicate true Garfield electron identity")
        full_seeds[key][i]=r
    full_hits={(r[KEY],r["configuration"]):r for r in
               read_csv(REF/"muon_reference_crossing.csv",CROSSING_COLUMNS)}
    muons=[];cluster_rows=[];replica_rows=[]
    # Hold-out physical response. Re-read common, exact Heed source from
    # full reference, NO second Heed call.
    for j,mu in enumerate(rows):
        id_=id_of(mu)
        p=sorted(primary[id_],key=lambda r:int(r["primary_index"]))
        cs=clusters[id_]
        if len(p)!=len(full_seeds[id_,"p3_noFW"]):
            raise ValueError("Muon Heed electron list differs between physical fields")
        for geom in GEOMS:
            phase_inputs=[];nearest_u=3. if geom.startswith("p3_") else 4.
            for e in p:
                idx=int(e["primary_index"]);cluster=int(e["cluster_index"])
                box,local_u,local_v=spatial_cell(
                    geom,int(mu["station"]),float(mu["u_in_mm"]),
                    float(e["u_rel_cm"]),float(e["v_abs_cm"]))
                full=full_seeds[id_,geom][idx]
                # C++ true-muon phase and Python LUT phase identity.
                true_u=float(full["seed_u_mm"])
                local_true=true_u-round_away(true_u/nearest_u)*nearest_u
                if abs(local_true-local_u)>1.e-3 and abs(
                    abs(local_true-local_u)-nearest_u)>1.e-3:
                    raise ValueError("Real Garfield electron and LUT spatial phase differ")
                if abs(float(full["seed_v_mm"])-local_v)>1.e-4:
                    raise ValueError("Real Garfield electron v differs from LUT sample")
                phase_inputs.append((idx,cluster,box))
            full_cluster=Counter()
            for idx,cl,box in phase_inputs:
                full_cluster[cl]+=int(full_seeds[id_,geom][idx]["ion_births"])
            rng=random.Random(stable_seed("HOLDOUT120-AVALANCHELUT20-FROM-DRIFT500-V1",id_,geom))
            ions=[];zeros=[];ah=[];at=[];cluster_draws=defaultdict(list)
            for rep in range(replicates):
                totals=[0,0,0,0,0];nz=0;bycl=Counter()
                for idx,cl,box in phase_inputs:
                    trials=lut[geom][box]
                    n,u,d,t,b=trials[rng.randrange(20)]
                    totals[0]+=n;totals[1]+=u;totals[2]+=d
                    totals[3]+=t;totals[4]+=b
                    nz+=int(n==0);bycl[cl]+=n
                for cl in cs:cluster_draws[int(cl["cluster_index"])].append(
                    bycl[int(cl["cluster_index"])])
                n,u,d,t,b=totals
                ions.append(n);zeros.append(nz)
                ah.append(asym(u,d));at.append(asym(t,b))
                replica_rows.append({
                    "crossing_id":id_,"configuration":geom,"replicate":rep,
                    "n_Heed_primaries":len(p),"ion_births":n,"zero_birth_Heed_seeds":nz,
                    "ion_up":u,"ion_down":d,"ion_top":t,"ion_bottom":b,
                    "A_hemi":ah[-1] if ah[-1] is not None else "",
                    "A_TB":at[-1] if at[-1] is not None else "",
                })
            r=full_hits[id_,geom]
            fion=int(r["ion_births"])
            fzero=sum(int(int(e["ion_births"])==0) for e in full_seeds[id_,geom].values())
            q=summarize_samples(ions);z=summarize_samples(zeros)
            h=summarize_samples([x for x in ah if x is not None]) if any(x is not None for x in ah) else None
            t=summarize_samples([x for x in at if x is not None]) if any(x is not None for x in at) else None
            fah=asym(int(r["ion_up"]),int(r["ion_down"]))
            fat=asym(int(r["ion_top"]),int(r["ion_bottom"]))
            muons.append({
                "crossing_id":id_,"configuration":geom,
                "station":mu["station"],"theta_u_deg":mu["theta_u_in_deg"],
                "theta_w_deg":mu["theta_w_in_deg"],
                "path_length_mm":mu["path_length_mm"],
                "n_Heed_clusters":len(cs),"n_Heed_primaries":len(p),
                "full_ion_births":fion,"lut_mean_ion_births":q["mean"],
                "lut_sd_ion_births":q["sd"],"lut_p05":q["p05"],"lut_p50":q["p50"],
                "lut_p95":q["p95"],"full_in_lut_90pct":int(q["p05"]<=fion<=q["p95"]),
                "full_zero_birth_Heed_seeds":fzero,"lut_mean_zero_birth_Heed_seeds":z["mean"],
                "full_A_hemi":fah if fah is not None else "",
                "lut_mean_A_hemi":h["mean"] if h else "",
                "lut_A_hemi_p05":h["p05"] if h else "",
                "lut_A_hemi_p95":h["p95"] if h else "",
                "full_A_TB":fat if fat is not None else "",
                "lut_mean_A_TB":t["mean"] if t else "",
                "lut_A_TB_p05":t["p05"] if t else "",
                "lut_A_TB_p95":t["p95"] if t else "",
                "A_hemi_in_lut_90pct":int(h["p05"]<=fah<=h["p95"]) if h and fah is not None else "",
                "A_TB_in_lut_90pct":int(t["p05"]<=fat<=t["p95"]) if t and fat is not None else "",
            })
            for cl in cs:
                k=int(cl["cluster_index"])
                if len(cluster_draws[k])!=replicates:
                    raise ValueError("Real Heed cluster response lost")
                vals=summarize_samples(cluster_draws[k]);fq=full_cluster[k]
                cluster_rows.append({
                    "crossing_id":id_,"configuration":geom,
                    "cluster_index":k,"cluster_t_ns":cl["t_ns"],
                    "cluster_u_rel_cm":cl["u_rel_cm"],"cluster_v_abs_cm":cl["v_abs_cm"],
                    "cluster_w_rel_cm":cl["w_rel_cm"],
                    "n_Heed_e_original":cl["n_original"],
                    "n_Heed_e_accepted":cl["n_accepted"],
                    "full_cluster_ion_births":fq,
                    "lut_cluster_mean_ion_births":vals["mean"],
                    "lut_cluster_sd_ion_births":vals["sd"],
                    "lut_cluster_p05":vals["p05"],"lut_cluster_p50":vals["p50"],
                    "lut_cluster_p95":vals["p95"],
                    "full_cluster_in_lut_90pct":int(vals["p05"]<=fq<=vals["p95"]),
                })
            if sum(full_cluster.values())!=fion:
                raise ValueError("Heed cluster sum != full Garfield muon ion births")
        if (j+1)%10==0:
            print(f"LUT HOLDOUT {j+1}/120 ORIGINAL Geant4 muons, "
                  f"all 4 geometries x {replicates} conditionally resampled",flush=True)
    if len(muons)!=480 or len(replica_rows)!=480*replicates:
        raise ValueError("Incomplete independent muon comparison")
    bygeom={};boots=random.Random(20261029)
    for geom in GEOMS:
        g=[r for r in muons if r["configuration"]==geom]
        full=[r["full_ion_births"] for r in g]
        means=[r["lut_mean_ion_births"] for r in g]
        # Total variance decomposition includes both between-muon Heed
        # ionization variance AND stochastic avalanche variance at fixed Heed.
        full_sd=statistics.stdev(full)
        lut_cond_var=statistics.mean(r["lut_sd_ion_births"]**2 for r in g)
        lut_total_sd=math.sqrt(statistics.pvariance(means)+lut_cond_var)
        ratios=[];differences_hemi=[];differences_tb=[]
        for k in range(1500):
            ids=[boots.randrange(len(g)) for _ in g]
            den=sum(full[i] for i in ids)
            ratios.append(sum(means[i] for i in ids)/den if den else float("nan"))
            d=[g[i]["lut_mean_A_hemi"]-g[i]["full_A_hemi"] for i in ids
               if g[i]["full_A_hemi"]!="" and g[i]["lut_mean_A_hemi"]!=""]
            differences_hemi.append(statistics.mean(d) if d else float("nan"))
            d=[g[i]["lut_mean_A_TB"]-g[i]["full_A_TB"] for i in ids
               if g[i]["full_A_TB"]!="" and g[i]["lut_mean_A_TB"]!=""]
            differences_tb.append(statistics.mean(d) if d else float("nan"))
        coverage=[r["full_in_lut_90pct"] for r in g]
        ah=[r["A_hemi_in_lut_90pct"] for r in g if r["A_hemi_in_lut_90pct"]!=""]
        at=[r["A_TB_in_lut_90pct"] for r in g if r["A_TB_in_lut_90pct"]!=""]
        bygeom[geom]={
            "n_distinct_true_Geant4_muons":len(g),
            "n_real_Heed_electrons":sum(int(r["n_Heed_primaries"]) for r in g),
            "mean_Garfield_ion_births_per_muon":statistics.mean(full),
            "mean_A20D500_ion_births_per_muon":statistics.mean(means),
            "A20D500_over_full_mean":sum(means)/sum(full),
            "ratio_mean_paired_bootstrap_95pct":[percent(sorted(ratios),.025),
                                                   percent(sorted(ratios),.975)],
            "Garfield_muon_gain_sd":full_sd,
            "A20D500_total_muon_gain_sd":lut_total_sd,
            "A20D500_over_full_sd":lut_total_sd/full_sd if full_sd else None,
            "n_full_in_A20D500_90pct":sum(coverage),
            "fraction_full_in_A20D500_90pct":statistics.mean(coverage),
            "n_full_A_hemi_in_A20D500_90pct":sum(ah),
            "n_defined_full_A_hemi":len(ah),
            "n_full_A_TB_in_A20D500_90pct":sum(at),
            "n_defined_full_A_TB":len(at),
            "bias_mean_A_hemi_A20_minus_full":
                statistics.mean(r["lut_mean_A_hemi"]-r["full_A_hemi"] for r in g
                                if r["full_A_hemi"]!=""),
            "bias_mean_A_TB_A20_minus_full":
                statistics.mean(r["lut_mean_A_TB"]-r["full_A_TB"] for r in g
                                if r["full_A_TB"]!=""),
            "bias_A_hemi_bootstrap_95pct":[percent(sorted(differences_hemi),.025),
                                          percent(sorted(differences_hemi),.975)],
            "bias_A_TB_bootstrap_95pct":[percent(sorted(differences_tb),.025),
                                        percent(sorted(differences_tb),.975)],
        }
    write_csv(COMP/"muon_comparison.csv",tuple(muons[0]),muons)
    write_csv(COMP/"heed_cluster_comparison.csv",tuple(cluster_rows[0]),cluster_rows)
    write_csv(COMP/"muon_A20D500_replicas.csv",tuple(replica_rows[0]),replica_rows)
    info={
        "schema":"PI_GEANT4_CHUNK2_HOLDOUT120_A20_FROM_D500_COMPARISON_V1",
        "status":"complete",
        "n_new_distinct_Geant4_muons":120,
        "n_geometries":4,"replicas_per_muon_geom":replicates,
        "n_LUT_conditional_muon_replicas":len(replica_rows),
        "n_real_Heed_cluster_geom_observations":len(cluster_rows),
        "full_Garfield_manifest_sha256":refsha,
        "A20D500_overnight_original_manifest_sha256":lutsha,
        "by_geometry":bygeom,
        "method":"Holdout 120 NEW independent Geant4 muons from source Jpsi chunk2; one exact original TrackHeed cluster/electron source per muon; same exact primaries go through full Garfield and through 200 conditionally independent 20-microscopic-avalanche empirical cell draws, four fields. 20 avalanche trials/cell were drawn from ORIGINAL 500 drift starts/cell, not 500 avalanche trials/cell.",
        "limits":"One microscopic Garfield realization per muon; bootstrap paired across 120 real independent source tracks. 200 LUT replicas reuse the same 20 empirical avalanches per spatial box, so they do NOT remove finite A20 training uncertainty. Stratified station/signed-angle pilot not the true flux angle PDF; cell bin is 0.1mm x 0.2mm (NO continuous interpolation); B=0; no strip weighting field/induced electronics or ion drift. Strong local/tail asymmetry conclusions need targeted higher statistics.",
        "git_commit":commit(),
    }
    write_json(COMP/"summary.json",info)
    write_json(COMP/"manifest.json",{
        **info,"outputs_sha256":{p.name:sha(p) for p in COMP.iterdir()
                                  if p.is_file() and p.name!="manifest.json"},
    })
    print("PASS HOLDOUT 120 UNSEEN real muons, "
          "four full Garfield fields versus A20@D500 "
          f"{replicates}× conditional replicates each (total {len(replica_rows)}); "
          "statistical summary and immutable SHA manifest",flush=True)
    for geom,s in bygeom.items():
        print(geom,json.dumps(s,sort_keys=True),flush=True)

def publish_all():
    rows,_,_=verify_full()
    mfile=COMP/"manifest.json"
    if not mfile.is_file():raise FileNotFoundError("First compare holdout, not published yet")
    m=json.loads(mfile.read_text())
    if m.get("status")!="complete":raise ValueError("Holdout comparison incomplete")
    for name,checksum in m.get("outputs_sha256",{}).items():
        if sha(COMP/name)!=checksum:raise ValueError("Holdout comparison SHA changed")
    for source,run_id in (
        (REF,"real_muon120_holdout_chunk2_full_Garfield_v1"),
        (COMP,"real_muon120_holdout_A20_from_D500_compare_v1"),
    ):
        pub=drive_prepare(source,run_id)
        publish(pub,EXPECTED_REMOTE,upload=True)
    print("PASS published two original holdout120 A20@D500 analysis bundles "
          "to Drive with immutable SHA verified rclone readback",flush=True)

def status():
    for d in (INPUT,REF,COMP):
        print("DIR",d,"present" if d.is_dir() else "missing")
    if (INPUT/"input_manifest.json").exists():
        print("INPUT",json.loads((INPUT/"input_manifest.json").read_text()).get("n_muons"))
    if REF.is_dir():
        done=[i for i in range(N_CHUNKS) if (chunk_dir(i)/"complete.json").is_file()]
        bad=[i for i in range(N_CHUNKS) if chunk_dir(i).exists()
             and not (chunk_dir(i)/"complete.json").is_file()]
        print("FULL_GARFIELD",len(done),"/12 completed physical checkpoints")
        if bad:print("WARNING incomplete chunk(s), inspect physics.log:",bad)
        print("REFERENCE_SEALED",(REF/"manifest.json").is_file())
    if COMP.is_dir():print("COMPARISON_SEALED",(COMP/"manifest.json").is_file())

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("mode",choices=("prepare","preflight","run","finish","compare","publish","pipeline","status"))
    ap.add_argument("--hours",type=float,default=7.75)
    ap.add_argument("--max-chunks",type=int,default=12)
    ap.add_argument("--replicates",type=int,default=REPLICATES)
    ap.add_argument("--upload",action="store_true",
                    help="Only pipeline: upload two immutable Drive archives after all 12 physical checkpoints pass")
    a=ap.parse_args()
    if a.upload and a.mode!="pipeline":
        ap.error("--upload is only valid with pipeline mode")
    if a.mode=="prepare":prepare()
    elif a.mode=="preflight":preflight()
    elif a.mode=="run":run(a.hours,a.max_chunks)
    elif a.mode=="finish":finish()
    elif a.mode=="compare":evaluate(a.replicates)
    elif a.mode=="publish":publish_all()
    elif a.mode=="pipeline":
        run(a.hours,a.max_chunks)
        if not (REF/"manifest.json").is_file():
            print("SAFE STOP pipeline partial: full Garfield not yet complete; "
                  "NO LUT validation or publication performed. Resume same command.",flush=True)
            return
        if not (COMP/"manifest.json").is_file():
            evaluate(a.replicates)
        else:
            print("PASS existing SHA-sealed holdout comparison; "
                  "will verify before optional publication",flush=True)
        if a.upload:
            publish_all()
        else:
            print("PASS pipeline local Garfield and A20@D500 comparison complete; "
                  "no upload requested (use --upload to publish)",flush=True)
    elif a.mode=="status":status()

if __name__=="__main__":
    main()
