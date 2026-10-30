#!/usr/bin/env python3
"""Targeted second microscopic avalanche sample for A20@D500.

20 original independent Garfield avalanches/cell remain IMMUTABLE.
Add 80 original D500 starting positions (passes 20..99) to ONLY the 50
largest expected finite-A20 source-response uncertainties in each of
four geometries (200 selected cells, 16,000 NEW avalanches total).

Rank uses ONLY frozen TrackHeed POSITIONS of 120 source muons, and old
A20 sample variances, NEVER full Garfield holdout outcomes/errors.
The 120-muon reference is unmodified and may later be used as a
position-selected validation, though another untouched source sample
is preferable for definitive generalization claims.

Modes: prepare | status | run | finish
Run is SHA-checkpointed in 400-row true Garfield chunks, resumable.
No old 84k A20 avalanches, frozen 500 drift points or 120 real muons
are rewritten. All in-metal nontransport outcomes remain explicit.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
import random
import subprocess
import tempfile
import time
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path

from pi_run_native_muon_b1 import EXE,CXX,check_executable,commit
from pi_publish_run_to_drive import (
    prepare as drive_prepare,publish as drive_publish,EXPECTED_REMOTE,
)
from pi_garfield_runtime_env import garfield_runtime_env
from pi_prepare_muon_avalanche_b1 import FIELDS,sha
from pi_audit_avalanche_lut_reuse import BASE,EXPORT_FIELDS,CONFIGS,locations
from pi_run_lut500_overnight import (
    ROOT as A20_ORIGINAL, source_maps,selected_j,
)
from pi_run_avalanche_lut_sidecar import (
    COLUMNS,BIRTH_COLUMNS,TIME_COLUMNS,validate_outputs,
)
from pi_run_muon_lut_validation import (
    spatial_cell,PRIMARY_FIELDS,read_csv,id_of,
    field_wire_blocked_zero_start,
)
from pi_holdout120_A20_from_D500 import (
    REF as FULL120,INPUT as HEED120,payload_rows,verify_full,
)
ROOT=BASE/"pi_A100_targeted_50cells_per_geometry_v1"
INPUT=ROOT/"selected_exact_old_D500_starts.csv"
RANK=ROOT/"cell_ranking.csv"
PLAN=ROOT/"target_cells.csv"
MANIFEST=ROOT/"plan_manifest.json"
CHUNK=400
TOP=50
ADD=80
TOTAL=TOP*4*ADD
N_CHUNKS=TOTAL//CHUNK
GEOMS=tuple(x[0] for x in CONFIGS)
RANK_FIELDS=("geometry","box_id","rank","n_Heed_electrons",
             "n_occupied_muons","A20_sd_ion_births",
             "finite_A20_mean_variance_score","fraction_of_geometry_score")
OUTFILE_FIELDS=("geometry","box_id","iu","iv","extra_n","old_n","new_n")

def write_csv(p,cols,records):
    with p.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=cols)
        w.writeheader();w.writerows(records)

def write_json(p,data):
    p.write_text(json.dumps(data,sort_keys=True,indent=2)+"\n")

def verify_plan():
    if not MANIFEST.is_file():
        raise FileNotFoundError("First run prepare; no original A100 targeted plan")
    m=json.loads(MANIFEST.read_text())
    if (m.get("schema")!="PI_A100_TARGETED_50CELLS_x4_V1"
        or m.get("n_new_full_avalanche_trials")!=TOTAL or
        m.get("n_selected_cells")!=TOP*4 or
        m.get("A20_old_original_manifest_sha256")!=sha(A20_ORIGINAL/"manifest.json")
        or m.get("orig_120_full_Garfield_manifest_sha256")!=sha(FULL120/"manifest.json")
        or m.get("CXX_sha256")!=sha(CXX)
        or m.get("compiled_executable_sha256")!=check_executable()):
        raise ValueError("Targeted campaign original source/model SHA or cardinality changed")
    for file,digest in m["plan_sha256"].items():
        if sha(ROOT/file)!=digest:
            raise ValueError("Targeted campaign original plan altered: "+file)
    rows=read_csv(INPUT,EXPORT_FIELDS)
    if len(rows)!=TOTAL:
        raise ValueError("Not all 16k original D500 source points in plan")
    all_ids={(r["geometry"],int(r["box_id"]),int(r["sample"])) for r in rows}
    if len(all_ids)!=TOTAL:
        raise ValueError("Duplicate physics trial in targeted plan")
    targets=read_csv(PLAN,OUTFILE_FIELDS)
    if len(targets)!=TOP*4:
        raise ValueError("Not exactly 50 selected cells in each of four fields")
    tkeys={(r["geometry"],int(r["box_id"])) for r in targets}
    for r in rows:
        g,b,s=r["geometry"],int(r["box_id"]),int(r["sample"])
        if (g,b) not in tkeys:
            raise ValueError("Targeted D500 start not in frozen priority cell list")
        if s in {selected_j(g,b,p) for p in range(20)}:
            raise ValueError("New avalanche start overlaps already run A20 starts")
        if s not in {selected_j(g,b,p) for p in range(20,100)}:
            raise ValueError("Extra trial not deterministic D500 pass 20..99")
    return m,rows,targets

def prepare():
    if ROOT.exists():
        raise FileExistsError("Existing targeted experiment MUST NOT be overwritten: "+str(ROOT))
    _,meta_original,src=payload_rows()
    _,full_stats,_=verify_full() # SHA verification of immutable full physics
    a20=json.loads((A20_ORIGINAL/"manifest.json").read_text())
    if (a20.get("status")!="complete" or
        a20.get("n_completed_avalanches")!=84000 or
        a20.get("total_capped")!=0 or a20.get("total_failed")!=0):
        raise RuntimeError("Original 84k A20 microscopic physics is not a completed exact baseline")
    cellfile=A20_ORIGINAL/"cell_statistics.csv"
    if sha(cellfile)!=a20["files_sha256"]["cell_statistics.csv"]:
        raise ValueError("Original A20 per-cell variances were modified")
    std={(r["geometry"],int(r["box_id"])):float(r["sd_ion_births"])
         for r in read_csv(cellfile,
           ("geometry","box_id","iu","iv","n","zero","p_zero",
            "mean_ion_births","sd_ion_births","median_ion_births",
            "mean_timing_seconds","ion_weighted_A_hemi","ion_weighted_A_TB"))}
    if len(std)!=4200:
        raise ValueError("Original A20 cell grid not 4200 full covered boxes")
    original={id_of(r):r for r in src}
    grouped=defaultdict(lambda:Counter())
    pri=read_csv(FULL120/"muon_reference_heed_primaries.csv",PRIMARY_FIELDS)
    for e in pri:
        mu=original[e["crossing_id"]]
        for g in GEOMS:
            b,_,_=spatial_cell(g,int(mu["station"]),float(mu["u_in_mm"]),
                               float(e["u_rel_cm"]),float(e["v_abs_cm"]))
            grouped[g,e["crossing_id"]][b]+=1
    score=Counter();counts=Counter();n_muons=Counter()
    for (g,mu),cnt in grouped.items():
        for box,n in cnt.items():
            # This is the finite A20 estimate of common cell-MEAN error.
            # Error repeated n times across one physical muon is
            # correlated: uncertainty grows as n^2, not as n.
            score[g,box]+=n*n*std[g,box]**2/20.
            counts[g,box]+=n
            n_muons[g,box]+=1
    all_rank=[];selected=[];groups={}
    for g in GEOMS:
        s=sorted(((b,score[g,b]) for (gg,b) in score if gg==g),
                 key=lambda item:(-item[1],item[0]))
        total=sum(v for b,v in s)
        if len(s)<TOP or total<=0:
            raise RuntimeError("Insufficient original Heed source coverage")
        for j,(b,val) in enumerate(s):
            all_rank.append({
                "geometry":g,"box_id":b,"rank":j+1,
                "n_Heed_electrons":counts[g,b],
                "n_occupied_muons":n_muons[g,b],
                "A20_sd_ion_births":std[g,b],
                "finite_A20_mean_variance_score":val,
                "fraction_of_geometry_score":val/total,
            })
        groups[g]=set(b for b,val in s[:TOP])
        for b,val in s[:TOP]:
            selected.append({"geometry":g,"box_id":b,"iu":b//30,
                             "iv":b%30,"extra_n":ADD,"old_n":20,
                             "new_n":100})
        print("TARGET",g,"top",TOP,"cover",
              f"{sum(x[1] for x in s[:TOP])/total:.3%}",
              "of estimated finite-A20 cell-mean uncertainty",
              flush=True)
    # Verify all 4 original frozen D500 CSVs and the source index SHA.
    orig_meta,index_sha=source_maps()
    chosen={}
    for g,pitch,fw in CONFIGS:
        _,raw=locations(500,pitch,fw)
        if sha(raw)!=orig_meta[g]["lut_csv_sha256"]:
            raise ValueError("Original immutable D500 input SHA changed: "+g)
        bysample={b:{selected_j(g,b,k):k for k in range(20,100)}
                  for b in groups[g]}
        seen=set()
        with raw.open(newline="") as f:
            rd=csv.DictReader(f)
            needed=("box_id","sample","iu","iv","random_seed",
                    "u_start_mm","v_start_mm","w_start_mm","drift_ok",
                    "collected","field_wire_collected","wire_index",
                    "phi_hit_deg","drift_time_ns",
                    "last_step_energy_eV")
            if not set(needed).issubset(rd.fieldnames or []):
                raise ValueError("Missing immutable original D500 seed coordinates")
            for r in rd:
                b=int(r["box_id"])
                if b not in bysample:continue
                sample=int(r["sample"])
                passage=bysample[b].get(sample)
                if passage is None:continue
                ident=(g,b,sample)
                if ident in seen:raise ValueError("Duplicate original D500 sample")
                seen.add(ident)
                if b!=int(r["iu"])*30+int(r["iv"]):
                    raise ValueError("Original D500 seed box identity malformed")
                if abs(float(r["w_start_mm"]))>1e-8:
                    raise ValueError("Original D500 seed out of 2D field")
                rec={
                    "lut_trials":500,"geometry":g,
                    "original_lut_sha256":orig_meta[g]["lut_csv_sha256"],
                    "box_id":b,"iu":r["iu"],"iv":r["iv"],
                    "sample":sample,"random_seed":r["random_seed"],
                    "u_start_mm":r["u_start_mm"],
                    "v_start_mm":r["v_start_mm"],
                    "w_start_mm":r["w_start_mm"],
                    "old_drift_ok":r["drift_ok"],
                    "old_collected":r["collected"],
                    "old_field_wire_collected":r["field_wire_collected"],
                    "old_wire_index":r["wire_index"],
                    "old_phi_hit_deg":r["phi_hit_deg"],
                    "old_drift_time_ns":r["drift_time_ns"],
                    "old_last_step_energy_eV":r["last_step_energy_eV"],
                }
                chosen[g,b,passage]=rec
        if len(seen)!=TOP*ADD:
            raise ValueError("Incomplete selected 80 original D500 seeds per important cell: "+g)
    # Pass-major order ensures every targeted cell gains balanced
    # independent physics as partial checkpoints accumulate.
    ordered=[]
    for passage in range(20,100):
        starts=[chosen[g,b,passage] for g in GEOMS for b in sorted(groups[g])]
        r=random.Random(20261030+passage)
        r.shuffle(starts)
        ordered.extend(starts)
    if len(ordered)!=TOTAL:
        raise AssertionError("Not exactly 16k additional physical starts")
    ROOT.mkdir(parents=True)
    write_csv(RANK,RANK_FIELDS,all_rank)
    write_csv(PLAN,OUTFILE_FIELDS,selected)
    write_csv(INPUT,EXPORT_FIELDS,ordered)
    write_json(MANIFEST,{
        "schema":"PI_A100_TARGETED_50CELLS_x4_V1",
        "created_utc":datetime.now(timezone.utc).isoformat(),
        "status":"prepared",
        "new_model":"A100@D500 for selected 50/geometry; A20@D500 in other original cells",
        "n_original_500_drift_starts_per_cell":500,
        "n_original_avalanche_trials_per_cell":20,
        "n_new_avalanche_trials_per_selected_cell":80,
        "n_selected_cells":TOP*4,
        "n_new_full_avalanche_trials":TOTAL,
        "n_original_full_avalanche_trials":84000,
        "n_chunks":N_CHUNKS,"chunk_size":CHUNK,
        "ranking":"Frozen 120-muon TrackHeed starting coordinates + frozen original 20-avalanche cell population SD only: sum over muons (n_primary_in_cell_per_muon**2 * A20_cell_sd**2 /20), do not use 120-muon full Garfield ion outcomes or model residuals",
        "new_trials":"Original D500 sample index (SHA_offset+137*pass)%500; pass=20..99; does not reuse the A20 passes 0..19",
        "orig_120_full_Garfield_manifest_sha256":sha(FULL120/"manifest.json"),
        "A20_old_original_manifest_sha256":sha(A20_ORIGINAL/"manifest.json"),
        "A20_cell_stats_sha256":sha(cellfile),
        "original_D500_by_geometry_sha256":{g:m["lut_csv_sha256"]
                                            for g,m in orig_meta.items()},
        "source_index_sha256":index_sha,
        "CXX_sha256":sha(CXX),
        "compiled_executable_sha256":check_executable(),
        "plan_sha256":{p.name:sha(p) for p in (RANK,PLAN,INPUT)},
        "git_commit":commit(),
        "scope":"One new microscopic Garfield avalanche per original saved D500 start. No native Geant4 or new Heed, B=0, no ion drift or induced strips. All old data immutable.",
    })
    print("PASS PREPARED immutable targeted 16,000 Garfield starts "
          "(50 cells x 80 new starts x four geometries),"
          "original A20/Geant4/Garfield reference unchanged; "
          "safe to run SHA-checkpointed 'run' mode.",flush=True)

def preflight():
    """Actual SHA-pinned C++ parser for ORIGINAL D500 source rows;
    NO Garfield/Heed physics or modifications to prepared source.
    """
    m,rows,targets=verify_plan()
    env=garfield_runtime_env()
    with tempfile.TemporaryDirectory() as td:
        td=Path(td)
        for i,batch in enumerate((rows[:CHUNK],rows[-CHUNK:])):
            p=td/f"source_preflight_{i}.csv"
            write_csv(p,EXPORT_FIELDS,batch)
            proc=subprocess.run(
                [str(EXE),"--validate-lut-input",str(p)],
                capture_output=True,text=True,check=False,env=env)
            marker=f"PASS parsed {CHUNK} ORIGINAL frozen LUT20/500 exact-start rows in C++"
            if proc.returncode!=0 or marker not in proc.stdout:
                raise RuntimeError(
                    "Actual SHA-pinned C++ D500 start parser refused "
                    f"segment {i}; NO physics performed:\n"+
                    proc.stdout+"\n"+proc.stderr)
    print("PASS compiled C++ D500 original exact-start preflight, "
          "first + last 400 rows parsed, NO Garfield/Heed physics",
          flush=True)


def chunk_dir(i):return ROOT/"chunks"/f"chunk_{i:03d}"

def validate_chunk(c,i,rows,binary):
    inp=c/"original_exact_starts.csv"
    if sha(inp)!=sha_csv_rows(EXPORT_FIELDS,rows):
        raise ValueError("Original chunk exact D500 starts differ")
    st=c/"complete.json"
    if not st.is_file():
        raise RuntimeError("Incomplete original Garfield chunk; inspect without overwriting: "+str(c))
    m=json.loads(st.read_text())
    if (m.get("status")!="complete" or m.get("chunk")!=i
        or m.get("n_rows")!=len(rows) or
        m.get("executable_sha256")!=binary):
        raise ValueError("Original campaign chunk model/binary changed")
    for fn,d in m["file_sha256"].items():
        if sha(c/fn)!=d:
            raise ValueError("SHA corruption in original targeted chunk: "+fn)
    validate_chunk_physics(c,rows)
    return m

def sha_csv_rows(cols,rows):
    import io
    b=io.StringIO(newline="")
    w=csv.DictWriter(b,fieldnames=cols)
    w.writeheader();w.writerows(rows)
    return hashlib.sha256(b.getvalue().encode()).hexdigest()

def validate_chunk_physics(c,rows):
    stat=validate_outputs(c,rows,reservoir_limit=5)
    if stat["n_cap"]:
        raise RuntimeError(f"Capped Garfield avalanche in {c}; preserve original files")
    tr=read_csv(c/"lut_avalanche_trials.csv",COLUMNS)
    blocked=[]
    for r in tr:
        if int(r["capped"]):raise RuntimeError("CAPPED single-electron physics")
        if int(r["successful"])!=1:
            if not field_wire_blocked_zero_start(r):
                raise RuntimeError(
                    "UNEXPLAINED unsuccessful original-start Garfield: "
                    f"{r['configuration']} cell={r['box_id']} "
                    f"sample={r['sample']} "
                    f"u={r['seed_u_mm']} v={r['seed_v_mm']}, "
                    f"ne={r['garfield_ne']} ni={r['garfield_ni']} "
                    f"births={r['ion_births']} other={r['other']}. "
                    "Original raw results preserved for investigation.")
            blocked.append((r["configuration"],int(r["box_id"]),int(r["sample"])))
    if stat["n_fail"]!=len(blocked):
        raise RuntimeError("Unexplained failed count mismatch")
    return stat,blocked

def run(hours,maxchunks):
    if not 0<hours<=168 or not 1<=maxchunks<=N_CHUNKS:
        raise ValueError("Unreasonable run budget/count")
    m,rows,targets=verify_plan()
    if (ROOT/"manifest.json").is_file():
        print("PASS original complete targeted A100 run already sealed; no Garfield repeated")
        return
    start=time.monotonic();n_done=0;n_new=0
    for i in range(N_CHUNKS):
        batch=rows[i*CHUNK:(i+1)*CHUNK]
        c=chunk_dir(i)
        if c.exists():
            validate_chunk(c,i,batch,m["compiled_executable_sha256"])
            n_done+=1
            print(f"PASS SHA-reused targeted Garfield chunk {i+1}/{N_CHUNKS}",flush=True)
            continue
        if n_new>=maxchunks or time.monotonic()-start>=hours*3600:
            print(f"SAFE STOP: {n_done}/{N_CHUNKS} original "
                  "40×400 Garfield chunks sealed; repeat identical command.",flush=True)
            return
        c.mkdir(parents=True)
        inp=c/"original_exact_starts.csv"
        write_csv(inp,EXPORT_FIELDS,batch)
        first=time.perf_counter()
        with (c/"physics.log").open("w") as log:
            cp=subprocess.run(
                [str(EXE),"--lut-sidecar-compact5",
                 str(inp),str(c/"lut_avalanche")],
                cwd=c,stdout=log,stderr=subprocess.STDOUT,
                env=garfield_runtime_env(),check=False)
        wall=time.perf_counter()-first
        if cp.returncode:
            raise RuntimeError("Garfield terminated in new targeted block "
                +str(c)+":\n"+
                "\n".join((c/"physics.log").read_text(errors="replace").splitlines()[-25:]))
        if "PASS: independent full-start microscopic avalanche sidecar" not in (
                c/"physics.log").read_text(errors="replace"):
            raise RuntimeError("Native C++ completion marker missing: "+str(c))
        stat,blocked=validate_chunk_physics(c,batch)
        files=[inp,c/"physics.log"]+[
            c/f"lut_avalanche_{tag}.csv" for tag in
            ("trials","birth_reservoir","timing")]
        write_json(c/"complete.json",{
            "schema":"PI_A100_TARGETED_GARFIELD_400_ORIGINAL_D500_STARTS_V1",
            "status":"complete","chunk":i,"n_rows":len(batch),
            "executable_sha256":m["compiled_executable_sha256"],
            "original_input_sha256":sha(inp),
            "full_physics_wall_seconds":wall,
            "physics_stats":{k:v for k,v in stat.items() if k!="files_sha256"},
            "n_verified_in_metal_no_transport":len(blocked),
            "verified_in_metal_original_ids":blocked,
            "file_sha256":{p.name:sha(p) for p in files},
        })
        n_done+=1;n_new+=1
        print(f"PASS TARGETED {n_done}/{N_CHUNKS}, "
              f"+{len(batch)} independent Garfield avalanche starts "
              f"({wall:.1f}s), physically embedded metal starts "
              f"{len(blocked)}; original SHA checkpoint sealed",flush=True)
    finish()

def finish():
    m,rows,targets=verify_plan()
    out=ROOT/"manifest.json"
    if out.exists():
        raise FileExistsError("Original targeted final reference already completed")
    for geom in GEOMS:
        if (ROOT/f"trials_extra_{geom}.csv").exists():
            raise FileExistsError(
                "Previous aggregate was interrupted; preserve original "
                "files and inspect before retrying, do not overwrite: "
                +str(ROOT/f"trials_extra_{geom}.csv"))
    if (ROOT/"summary.json").exists():
        raise FileExistsError("Preserve existing original A100 summary")
    bundles=[]
    for i in range(N_CHUNKS):
        batch=rows[i*CHUNK:(i+1)*CHUNK]
        bundles.append(validate_chunk(chunk_dir(i),i,batch,
                                      m["compiled_executable_sha256"]))
    paths={}
    count=Counter();blocked=Counter()
    try:
        for g in GEOMS:
            p=ROOT/f"trials_extra_{g}.csv"
            f=p.open("w",newline="")
            w=csv.DictWriter(f,fieldnames=COLUMNS)
            w.writeheader();paths[g]=(f,w)
        seen=set()
        for i in range(N_CHUNKS):
            for r in read_csv(chunk_dir(i)/"lut_avalanche_trials.csv",COLUMNS):
                g=r["configuration"];b=int(r["box_id"]);s=int(r["sample"])
                if (g,b,s) in seen or g not in paths:
                    raise ValueError("Duplicate/unknown new independent original D500 avalanche")
                seen.add((g,b,s))
                count[g,b]+=1
                if int(r["successful"])!=1:blocked[g]+=1
                paths[g][1].writerow(r)
    finally:
        for f,w in paths.values():f.close()
    expect={(r["geometry"],int(r["box_id"])) for r in targets}
    if (len(seen)!=TOTAL or set(count)!=expect or
        any(v!=ADD for v in count.values())):
        raise RuntimeError("New A100 per-cell coverage incomplete")
    result={
        "schema":"PI_A100_TARGETED_50CELLS_x4_COMPLETE_V1",
        "status":"complete",
        "new_full_Garfield_avalanches":TOTAL,
        "distinct_targeted_cells":TOP*4,
        "additional_Garfield_per_selected_cell":ADD,
        "total_original_A20plus_extra_avalanches":84000+TOTAL,
        "n_physical_Garfield_chunks":N_CHUNKS,
        "original_target_plan_sha256":sha(MANIFEST),
        "original_A20_manifest_sha256":m["A20_old_original_manifest_sha256"],
        "original_120_full_Garfield_manifest_sha256":
            m["orig_120_full_Garfield_manifest_sha256"],
        "genuine_new_physics_wall_seconds":
            sum(x["full_physics_wall_seconds"] for x in bundles),
        "metal_start_zero_births_by_geometry":dict(blocked),
        "no_censored_nonmetal_or_capped_trials":True,
        "A100_vs_A20":"Selected 50 cells per geometry now have 100 independent full Garfield starts, other 4,000 cells retain 20 original starts. The 16k original extra results are in four trials_extra_[geometry].csv files, separated from the immutable 84k.",
    }
    write_json(ROOT/"summary.json",result)
    result["files_sha256"]={p.name:sha(p) for p in ROOT.iterdir()
        if p.is_file() and p.name!="manifest.json"}
    result["chunk_complete_manifest_sha256"]={
        str(i):sha(chunk_dir(i)/"complete.json") for i in range(N_CHUNKS)}
    write_json(out,result)
    print("PASS TARGETED A100 full microscopic Garfield complete, "
          "16,000 new original D500 starts, SHA-sealed without "
          "modifying original 84k or 120 real muons",flush=True)

def publish_results():
    _,rows,targets=verify_plan()
    f=ROOT/"manifest.json"
    if not f.exists():
        raise FileNotFoundError("Targeted original Garfield 16k full run not finished")
    data=json.loads(f.read_text())
    if data.get("status")!="complete" or data.get("new_full_Garfield_avalanches")!=TOTAL:
        raise ValueError("Not a completed original targeted 16k Garfield source")
    for filename,digest in data["files_sha256"].items():
        if sha(ROOT/filename)!=digest:
            raise ValueError("Original A100 physics changed after SHA seal: "+filename)
    bundle=drive_prepare(ROOT,"mwpc_A100D500_targeted50x4_original_Garfield_v1")
    drive_publish(bundle,EXPECTED_REMOTE,upload=True)
    print("PASS published frozen 16k original targeted A100 Garfield + "
          "per-cell targets and source hashes (no source overwrite)",flush=True)


def status():
    print("TARGET ROOT",ROOT,"present" if ROOT.is_dir() else "absent")
    if MANIFEST.is_file():
        print("PLAN prepared",json.loads(MANIFEST.read_text()).get(
            "n_new_full_avalanche_trials"))
    complete=[i for i in range(N_CHUNKS)
              if (chunk_dir(i)/"complete.json").is_file()]
    partial=[i for i in range(N_CHUNKS)
             if chunk_dir(i).exists() and i not in complete]
    print(f"TARGETED GARFIELD {len(complete)}/{N_CHUNKS} blocks")
    if partial:print("WARNING incomplete original chunk(s):",partial)
    print("SEALED COMPLETE", (ROOT/"manifest.json").is_file())

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("mode",choices=("prepare","preflight","run","finish","publish","status"))
    p.add_argument("--hours",type=float,default=6.75)
    p.add_argument("--max-chunks",type=int,default=N_CHUNKS)
    a=p.parse_args()
    if a.mode=="prepare":prepare()
    elif a.mode=="preflight":preflight()
    elif a.mode=="run":run(a.hours,a.max_chunks)
    elif a.mode=="finish":finish()
    elif a.mode=="publish":publish_results()
    elif a.mode=="status":status()

if __name__=="__main__":main()
