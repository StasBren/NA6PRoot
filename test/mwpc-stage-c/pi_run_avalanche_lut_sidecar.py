#!/usr/bin/env python3
"""Avalanche sidecar for ORIGINAL LUT20 and LUT500 saved electron STARTS.

No Geant4 or TrackHeed here: this is a NEW independent microscopic
AvalancheElectron realization from the exact saved original start of each
frozen single-electron Garfield DRIFT trial. It is NOT a continuation
from an anode endpoint and NOT a native-muon physics measurement.

Modes: check | smoke8 | full432. The full432 phase is deliberately separate.
No source LUTs or prior runs are altered.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import shutil
import subprocess
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

from pi_audit_avalanche_lut_reuse import (
    BASE, CONFIGS, EXPORT_FIELDS, locations, sha,
)
from pi_run_native_muon_b1 import CXX, EXE, check_executable, commit
from pi_garfield_runtime_env import garfield_runtime_env

ROOT=BASE/"pi_avalanche_lut_reuse_v1"
INDEX=ROOT/"frozen_exact_starts_all"
COLUMNS=("lut_trials","original_geometry","original_lut_sha256","box_id",
    "iu","iv","sample","original_random_seed","original_u_start_mm",
    "original_v_start_mm","original_w_start_mm","old_collected",
    "old_field_wire_collected","old_wire_index","old_phi_hit_deg",
    "old_drift_time_ns","old_last_step_energy_eV",
    "crossing_id","configuration","primary_index","garfield_seed",
    "seed_u_mm","seed_v_mm","seed_w_mm","seed_energy_eV",
    "garfield_ne","garfield_ni","ion_births","collected_e",
    "attached","other","successful","capped","ion_up","ion_down",
    "ion_top","ion_bottom","A_ion_hemi","A_ion_TB")
BIRTH_COLUMNS=("crossing_id","configuration","primary_index","u_mm",
               "v_mm","w_mm","t_ns","seed_ion_births_total","reservoir_weight")
TIME_COLUMNS=("crossing_id","station","theta_u_deg","theta_w_deg",
              "configuration","primary_index","n_true_heed_primaries",
              "u_seed_mm","v_seed_mm","garfield_seed","wall_seconds",
              "avalanche_electrons","avalanche_ions","ion_births",
              "collected_e","success","capped")

def load_index():
    manifest_path=INDEX/"index_manifest.json"
    csv_path=INDEX/"lut_exact_start_index.csv"
    if not manifest_path.is_file() or not csv_path.is_file():
        raise FileNotFoundError("Run pi_audit_avalanche_lut_reuse.py --trials all --export first")
    m=json.loads(manifest_path.read_text())
    if (m.get("schema")!="NA60_PI_LUT_AVALANCHE_EXACT_START_AUDIT_V1"
        or m.get("n_index_rows")!=432
        or m.get("index_csv_sha256")!=sha(csv_path)):
        raise ValueError("Frozen exact-start index SHA/identity changed")
    source_infos=m.get("input_luts",[])
    if len(source_infos)!=8:
        raise ValueError("Expected exact 8 frozen LUT20/LUT500 source identities")
    known=set()
    for info in source_infos:
        kind=info["trials"]
        geom=info["geometry"]
        pitch_fw={name:(pitch,fw) for name,pitch,fw in CONFIGS}
        if kind not in (20,500) or geom not in pitch_fw or (kind,geom) in known:
            raise ValueError("Incorrect/duplicate original LUT identity")
        known.add((kind,geom))
        path_conf,path_csv=locations(kind,*pitch_fw[geom])
        if (not path_csv.is_file() or not path_conf.is_file()
            or sha(path_csv)!=info["lut_csv_sha256"]
            or sha(path_conf)!=info["frozen_config_sha256"]):
            raise ValueError("Immutable LUT source changed: "+str(path_csv))
    with csv_path.open(newline="") as f:
        rd=csv.DictReader(f)
        if tuple(rd.fieldnames or ())!=EXPORT_FIELDS:
            raise ValueError("Frozen sidecar index CSV columns changed")
        rows=list(rd)
    if len(rows)!=432:
        raise ValueError("Wrong exact-start index cardinality")
    ids=set()
    shas={(info["trials"],info["geometry"]):info["lut_csv_sha256"]
          for info in source_infos}
    count=Counter()
    for r in rows:
        kind=int(r["lut_trials"]);geom=r["geometry"]
        key=(kind,geom)
        if key not in shas or r["original_lut_sha256"]!=shas[key]:
            raise ValueError("Original SHA/geometry mismatch on exact-start row")
        id_=(kind,geom,int(r["box_id"]),int(r["sample"]))
        if id_ in ids:raise ValueError("Duplicate frozen electron start")
        ids.add(id_);count[key]+=1
        if int(r["old_collected"])+int(r["old_field_wire_collected"])>1:
            raise ValueError("Inconsistent original anode/FW classification")
    if set(count)!=set(shas) or any(n!=54 for n in count.values()):
        raise ValueError("Incorrect frozen 54 exact-start samples per LUT geometry")
    return m,rows

def choose_rows(rows, mode):
    if mode=="full432":return rows
    if mode!="smoke8":raise ValueError(mode)
    chosen=[]
    for kind in (20,500):
        for geom,pitch,fw in CONFIGS:
            target=(kind,geom,round(pitch/.1)//2,15,kind//2)
            row=next((r for r in rows if (
                int(r["lut_trials"]),r["geometry"],int(r["iu"]),
                int(r["iv"]),int(r["sample"]))==target),None)
            if row is None:raise ValueError("Missing central representative frozen seed "+str(target))
            chosen.append(row)
    return chosen

def validate_outputs(run_dir, chosen, *, reservoir_limit=100):
    if not isinstance(reservoir_limit,int) or reservoir_limit<0 or reservoir_limit>100:
        raise ValueError("Invalid ion birth reservoir sample limit")
    keys={(int(r["lut_trials"]),r["geometry"],int(r["box_id"]),
           int(r["sample"])):r for r in chosen}
    files=[run_dir/("lut_avalanche"+suffix) for suffix in (
        "_trials.csv","_birth_reservoir.csv","_timing.csv")]
    data=[]
    for name,expected in zip(files,(COLUMNS,BIRTH_COLUMNS,TIME_COLUMNS)):
        with name.open(newline="") as f:
            rd=csv.DictReader(f)
            if tuple(rd.fieldnames or ())!=expected:
                raise ValueError("Wrong sidecar CSV columns "+str(name))
            data.append(list(rd))
    trials,births,times=data
    if len(trials)!=len(chosen) or len(times)!=len(chosen):
        raise ValueError("Incomplete sidecar: one avalanche/timing row per original LUT start required")
    seen=set();totals={};caps=0;failed=0
    for r in trials:
        key=(int(r["lut_trials"]),r["original_geometry"],int(r["box_id"]),
             int(r["sample"]))
        if key not in keys or key in seen:
            raise ValueError("Wrong/duplicate physical LUT trial in sidecar")
        seen.add(key)
        orig=keys[key]
        if (r["original_lut_sha256"]!=orig["original_lut_sha256"] or
            r["configuration"]!=orig["geometry"] or
            int(r["original_random_seed"])!=int(orig["random_seed"]) or
            r["old_collected"]!=orig["old_collected"] or
            r["old_field_wire_collected"]!=orig["old_field_wire_collected"]):
            raise ValueError("LUT identities and endpoint flags no longer paired")
        for axis in "uvw":
            expected=float(orig[f"{axis}_start_mm"])
            observed=float(r[f"seed_{axis}_mm"])
            if not math.isclose(expected,observed,abs_tol=1.e-7):
                raise ValueError("Avalanche no longer starts from original saved exact (u,v,w)")
        ni=int(r["ion_births"])
        up,down,boundary=int(r["ion_up"]),int(r["ion_down"]),None
        if ni<0 or up+down>ni:
            raise ValueError("Inconsistent avalanche ion hemisphere census")
        top=int(r["ion_top"]);bottom=int(r["ion_bottom"])
        if top+bottom>ni:raise ValueError("Invalid ion polar sector census")
        if int(r["garfield_ni"])>0 and ni==0:
            raise ValueError("Garfield positive ions with no ionization callback")
        caps+=int(r["capped"])
        failed+=int(r["successful"])!=1
        totals[r["crossing_id"]]=ni
    if seen!=set(keys):
        raise ValueError("Missing LUT sidecar trials")
    time_ids={r["crossing_id"] for r in times}
    if time_ids!=set(totals) or len(time_ids)!=len(times):
        raise ValueError("Timing missing or duplicated unique single-electron labels")
    for row in times:
        if (not math.isfinite(float(row["wall_seconds"])) or
            float(row["wall_seconds"])<0):
            raise ValueError("Invalid microscopic wall-clock time")
    grouped=Counter()
    for r in births:
        k=r["crossing_id"]
        if k not in totals:raise ValueError("Reservoir contains unknown seed")
        if int(r["seed_ion_births_total"])!=totals[k]:
            raise ValueError("Reservoir full-ion denominator mismatch")
        grouped[k]+=1
    for key,n in totals.items():
        if grouped[key]!=min(n,reservoir_limit):
            raise ValueError("Ion-birth reservoir is incomplete/biased")
    return {
        "n_trials":len(trials),"n_timed":len(times),
        "n_reservoir_ion_positions":len(births),
        "n_cap":caps,"n_fail":failed,
        "by_geometry":dict(Counter(r["configuration"] for r in trials)),
        "by_original_lut":dict(Counter(str(r["lut_trials"]) for r in trials)),
        "total_ion_births":sum(totals.values()),
        "mean_seed_wall_seconds":sum(float(r["wall_seconds"]) for r in times)/len(times),
        "files_sha256":{p.name:sha(p) for p in files}
    }

def run(mode):
    index_manifest,rows=load_index()
    print("PASS 8 original frozen LUT20/500 SHA checks; 432 exact-start indexed electrons")
    exe_sha=check_executable()
    # Parse the ACTUAL SHA-verified input in the current C++ executable
    # before allocating a new immutable physics output directory.
    cpp=subprocess.run([str(EXE),"--validate-lut-input",
                        str(INDEX/"lut_exact_start_index.csv")],
                       capture_output=True,text=True,check=False,
                       env=garfield_runtime_env())
    if cpp.returncode or "PASS parsed 432 ORIGINAL frozen LUT20/500" not in cpp.stdout:
        raise RuntimeError("B1 C++ exact-start schema preflight failed:\n"+
                           cpp.stdout+"\n"+cpp.stderr)
    print(cpp.stdout.strip())
    if mode=="check":return
    selected=choose_rows(rows,mode)
    run_dir=ROOT/f"run_{mode}_independent_full_start_avalanche_v1"
    if run_dir.exists():
        raise FileExistsError("Preserving existing sidecar pilot (success or failure): "+str(run_dir))
    if mode=="full432":
        smoke=ROOT/"run_smoke8_independent_full_start_avalanche_v1"/"manifest.json"
        if not smoke.is_file() or json.loads(smoke.read_text()).get("status")!="complete":
            raise RuntimeError("First validate smoke8 before running all 432 independent avalanches")
    runtime_env=garfield_runtime_env()
    run_dir.mkdir(parents=True)
    subset=run_dir/"frozen_exact_start_subset.csv"
    with subset.open("w",newline="") as f:
        wr=csv.DictWriter(f,fieldnames=EXPORT_FIELDS)
        wr.writeheader();wr.writerows(selected)
    manifest={
      "schema":"PI_LUT20_LUT500_INDEPENDENT_FULL_START_MICROSCOPIC_AVALANCHE_V1",
      "status":"incomplete","created_utc":datetime.now(timezone.utc).isoformat(),
      "mode":mode,"n_original_exact_starts":len(selected),
      "source_index_csv_sha256":index_manifest["index_csv_sha256"],
      "source_index_manifest_sha256":sha(INDEX/"index_manifest.json"),
      "selected_exact_start_sha256":sha(subset),
      "original_lut_sha256_by_geometry":index_manifest["input_luts"],
      "compiled_cxx_sha256":sha(CXX),"compiled_executable_sha256":exe_sha,
      "git_commit":commit(),"heed_database":runtime_env["HEED_DATABASE"],
      "model":"One new full microscopic AvalancheElectron per ORIGINAL saved (geometry,box_id,sample) START, with independent random seed; NOT a continuation of the old drift endpoint",
      "limitations":"NOT native-muon or physical angle PDF; stratified artificial single-electron grid; never double gate on old collected status; no field signal or reconstruction",
    }
    (run_dir/"run_config.json").write_text(json.dumps(manifest,indent=2)+"\n")
    cmd=[str(EXE),"--lut-sidecar",str(subset),str(run_dir/"lut_avalanche")]
    with (run_dir/"physics.log").open("w") as log:
        status=subprocess.run(cmd,cwd=run_dir,stdout=log,
                              stderr=subprocess.STDOUT,env=runtime_env,check=False).returncode
    if status!=0:
        raise RuntimeError(f"LUT sidecar failed with code {status}; preserved {run_dir}. Last log lines:\n"+
                           "\n".join((run_dir/"physics.log").read_text(errors="replace").splitlines()[-25:]))
    stats=validate_outputs(run_dir,selected)
    manifest.update(stats)
    manifest["status"]="complete" if not stats["n_cap"] and not stats["n_fail"] else "complete_with_warning"
    (run_dir/"summary.json").write_text(json.dumps(manifest,indent=2)+"\n")
    manifest["output_sha256"]={
        p.name:sha(p) for p in sorted(run_dir.iterdir())
        if p.is_file() and p.name!="manifest.json"
    }
    (run_dir/"manifest.json").write_text(json.dumps(manifest,indent=2)+"\n")
    print("PASS independent original-LUT exact-start avalanche pilot",len(selected),"trials")
    print("  by geometry:",stats["by_geometry"])
    print("  capped:",stats["n_cap"],"failed:",stats["n_fail"],
          "ion births:",stats["total_ion_births"])
    print("  mean avalanche wall seconds:",round(stats["mean_seed_wall_seconds"],4))
    print("  output:",run_dir)
    if manifest["status"]!="complete":
        print("WARNING: caps/failures detected; do not treat this as a validated physics LUT")
    print("NOT muon-driven. Do not condition on old_collected and count losses twice.")

if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("mode",choices=("check","smoke8","full432"))
    run(p.parse_args().mode)
