#!/usr/bin/env python3
"""Read-only source verification and safe first true-muon B1 smoke execution.

No Geant4 production rerun; no old LUT/Stage-A/B overwrite.
One *physical* Geant4 muon crossing is one analyzed observation.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
import os
import shutil
import subprocess
import sys
from datetime import datetime,timezone
from pathlib import Path
from collections import Counter

from pi_prepare_muon_avalanche_b1 import OUT, FIELDS, sha, make, REPO

HERE=Path(__file__).resolve().parent
CXX=REPO/"test/mwpc-avalanche/phaseB_native_muon_3p3_avalanche.cxx"
EXE=REPO/"test/mwpc-avalanche/build/mwpc_phase_b_native_muon_3p3"
def run_dir(n, retry=0):
    # Failed runs remain byte-for-byte intact for forensics. A C++ gas
    # source/medium fix needs its OWN output directory, not overwrite the
    # earlier original/CRLF failure records.
    if retry==0:
        name=f"run_smoke{n}_chunk0001"
    elif retry==1:
        name=f"run_smoke{n}_crlf_fix1_chunk0001"
    elif retry==2:
        name=f"run_smoke{n}_gas_medium_fix2_chunk0001"
    else:
        raise ValueError("Unrecognized immutable B1 retry identifier")
    return OUT/name
def input_dir(n):return OUT/f"smoke{n}_chunk0001"
GEOMS=("p3_noFW","p3_FW75","p4_noFW","p4_FW75")
PREFIX="native_muon_b1"
CROSSING_COLUMNS=("source|chunk|event|hit_index|track_id|station|det_id",
  "configuration","station","det_id","theta_u_in_deg","theta_w_in_deg",
  "path_mm","momentum_GeV","heed_clusters","heed_primaries","heed_discarded",
  "heed_failed","collected_avalanche_electrons","ion_births","garfield_ion_count",
  "zero_seeds","capped_seeds","failed_seeds","ion_up","ion_down","ion_boundary",
  "ion_top","ion_bottom","ion_side","A_ion_hemi","A_ion_TB","e_up","e_down",
  "e_top","e_bottom","A_e_hemi","A_e_TB","G_eff_per_Heed_primary")

def commit():
    return subprocess.check_output(["git","-C",str(REPO),"rev-parse","HEAD"],text=True).strip()

def source_check():
    if not CXX.is_file():raise FileNotFoundError(CXX)
    s=CXX.read_text()
    for needed in ('SetUserHandleIonisation(OnIonisation)',
                   'heed.GetCluster','heed.GetElectron',
                   'p3_noFW','p3_FW75','p4_noFW','p4_FW75'):
        if needed not in s:raise RuntimeError("Native-muon B1 C++ source incomplete: "+needed)
    return sha(CXX)

def check_input(n):
    base=input_dir(n)
    selection=base/"native_muon_crossings_input.csv"
    meta=base/"input_manifest.json"
    if not selection.is_file() or not meta.is_file():
        raise FileNotFoundError("First run 'prepare' for 8 original native Geant4 muons")
    m=json.loads(meta.read_text())
    if m.get("schema")!="PI_TRUE_GEANT4_MUON_B1_INPUT_V1" or m.get("status")!="complete":
        raise ValueError("Wrong/incomplete Geant4 selection provenance")
    if sha(selection)!=m.get("input_sha256") or m.get("n_crossings")!=n:
        raise ValueError("Native Geant4 muon input changed or wrong n")
    with selection.open(newline="") as f:
        reader=csv.DictReader(f)
        if tuple(reader.fieldnames or ())!=FIELDS:raise ValueError("Native CSV changed schema")
        rows=list(reader)
    if len(rows)!=n or any(abs(int(r["pdg"]))!=13 for r in rows):
        raise ValueError("Unexpected native muon selection")
    return selection,m,rows

def check_executable():
    source_check()
    if not EXE.is_file():raise FileNotFoundError("Build Garfield C++ target first: "+str(EXE))
    if EXE.stat().st_mtime_ns < CXX.stat().st_mtime_ns:
        raise RuntimeError("C++ source newer than executable; rebuild before running")
    return sha(EXE)

def validate_native_csv_before_physics(selection):
    # Test the ACTUAL C++ CSV decoder on the frozen original Python output
    # before starting any Heed or Garfield physics work.
    res=subprocess.run([str(EXE),"--validate-input",str(selection)],
                       capture_output=True,text=True,check=False)
    if res.returncode!=0:
        raise RuntimeError("C++ pre-physics input validation failed:\n"+
                           res.stdout+"\n"+res.stderr)
    if ("PASS: parsed" not in res.stdout or
        "PASS: Heed gas source medium is ionisable" not in res.stdout):
        raise RuntimeError("B1 C++ missing HEED gas-medium preflight marker:\n"+
                           res.stdout)
    print(res.stdout.strip())

def check_crossings(path,selected):
    with path.open(newline="") as f:
        reader=csv.DictReader(f)
        if tuple(reader.fieldnames or ())!=CROSSING_COLUMNS:
            raise ValueError("Wrong B1 Geant4 crossing output CSV columns")
        records=list(reader)
    if len(records)!=len(selected)*4:
        raise ValueError("B1 did not finish every crossing in all four geometries")
    original={
        "|".join(str(r[k]) for k in ("source","chunk","event","hit_index","track_id",
                                     "station","det_id")):r
        for r in selected
    }
    counts=Counter()
    by_geom={name:{"n_crossings":0,"n_heed_primaries":0,"n_avalanche_ion_births":0,
                   "n_collected_electrons":0,"n_zero_seeds":0,
                   "n_capped_seeds":0,"n_failed_seeds":0,"defined_A_hemi":0,
                   "defined_A_TB":0,"mean_A_hemi":0.,"mean_A_TB":0.}
             for name in GEOMS}
    matched={}
    for row in records:
        id_=row[CROSSING_COLUMNS[0]]
        name=row["configuration"]
        if id_ not in original or name not in GEOMS:
            raise ValueError("Unknown true Geant4 crossing/configuration")
        key=id_,name
        counts[key]+=1
        if counts[key]!=1:raise ValueError("Duplicate muon/config")
        original_row=original[id_]
        st=int(row["station"]);det=int(row["det_id"])
        if st!=int(original_row["station"]) or det!=int(original_row["det_id"]):
            raise ValueError("Muon station/detID mismatch")
        for coordinate in ("theta_u_in_deg","theta_w_in_deg","path_mm"):
            original_coordinate = "path_length_mm" if coordinate=="path_mm" else coordinate
            if abs(float(row[coordinate])-float(original_row[original_coordinate]))>1.e-6:
                raise ValueError("Original muon angle or gas path length changed: "+coordinate)
        if abs(float(row["momentum_GeV"])-float(original_row["momentum_GeV"]))>1.e-6:
            raise ValueError("Original muon momentum changed")
        nprim=int(row["heed_primaries"])
        ncluster=int(row["heed_clusters"])
        if nprim<0 or nprim>1000 or ncluster<0 or ncluster>nprim+10:
            raise ValueError("Impossible Heed cluster/electron multiplicity")
        if id_ in matched and matched[id_]!=(nprim,ncluster):
            raise ValueError("Heed primary ionization differs across geometries")
        matched[id_]=(nprim,ncluster)
        up=int(row["ion_up"]);down=int(row["ion_down"]);bound=int(row["ion_boundary"])
        top=int(row["ion_top"]);bot=int(row["ion_bottom"]);side=int(row["ion_side"])
        born=int(row["ion_births"]);ne=int(row["collected_avalanche_electrons"])
        if up+down+bound!=born or top+bot+side!=born:
            raise ValueError("B1 angular ionisation-birth counts violate conservation")
        for a,aa,bb in (("A_ion_hemi",up,down),("A_ion_TB",top,bot)):
            s=row[a]
            if aa+bb==0:
                if s!="":raise ValueError("Undefined A stored as numeric")
            else:
                val=float(s)
                if not math.isfinite(val) or abs(val-(aa-bb)/(aa+bb))>5.e-9:
                    raise ValueError("Angular A inconsistent with ion counts")
        for a,upkey,downkey in (("A_e_hemi","e_up","e_down"),
                                 ("A_e_TB","e_top","e_bottom")):
            u=int(row[upkey]);d=int(row[downkey])
            if u+d>ne:raise ValueError("Endpoint angular count exceeds collected electrons")
            if (u+d==0)!=(row[a]==""):
                raise ValueError("Undefined amplified-electron A not marked empty")
        eff=row["G_eff_per_Heed_primary"]
        if nprim:
            if not eff or abs(float(eff)-ne/nprim)>1.e-7:
                raise ValueError("Wrong full-hit gain normalization")
        elif eff:
            raise ValueError("G undefined for zero-Heed-primary crossing")
        d=by_geom[name]
        d["n_crossings"]+=1;d["n_heed_primaries"]+=nprim
        d["n_avalanche_ion_births"]+=born;d["n_collected_electrons"]+=ne
        d["n_zero_seeds"]+=int(row["zero_seeds"])
        d["n_capped_seeds"]+=int(row["capped_seeds"])
        d["n_failed_seeds"]+=int(row["failed_seeds"])
        for a,label in (("A_ion_hemi","hemi"),("A_ion_TB","TB")):
            if row[a]!="":
                d["defined_A_"+label]+=1
                d["mean_A_"+label]+=float(row[a])
    if len(counts)!=len(selected)*4:raise ValueError("Missing geometry/crossing pair")
    for name,d in by_geom.items():
        for label in ("hemi","TB"):
            n=d["defined_A_"+label]
            d["mean_A_"+label]=d["mean_A_"+label]/n if n else None
        d["G_eff_total_collected_over_total_Heed_primary"]=(
            d["n_collected_electrons"]/d["n_heed_primaries"]
            if d["n_heed_primaries"] else None)
    return by_geom

def run(mode):
    if mode.startswith("prepare"):
        n=int(mode.replace("prepare","") or "2")
        if n not in (2,8):raise ValueError("Only prepare2 or prepare8 is allowed")
        source_check()
        make(max_events=n,chunk=1)
        return
    if mode=="status":
        for n in (2,8):
            for d in (input_dir(n),run_dir(n)):
                print(d, "exists" if d.exists() else "missing")
            if n==2:
                for retry in (1,2):
                    d=run_dir(n,retry=retry)
                    print(d,"exists" if d.exists() else "missing")
            if (run_dir(n)/"manifest.json").is_file():
                m=json.loads((run_dir(n)/"manifest.json").read_text())
                print("Completed:",m["status"],"crossings:",m["n_true_geant4_crossings"])
        return
    if mode.startswith("check"):
        n=int(mode.replace("check","") or "2")
        print("Source SHA:",source_check())
        print("Binary SHA:",check_executable())
        selection, input_meta, _ = check_input(n)
        print("Input SHA:",input_meta["input_sha256"])
        validate_native_csv_before_physics(selection)
        print(f"PASS static B1 executable and {n} immutable real native muon inputs")
        return
    if mode not in ("smoke2","smoke8","smoke2retry1","smoke2retry2"):
        raise ValueError("Use prepare2|prepare8|check2|check8|smoke2|smoke2retry1|smoke2retry2|smoke8|status")
    retry = 2 if mode=="smoke2retry2" else (1 if mode=="smoke2retry1" else 0)
    n = 2 if retry else int(mode.replace("smoke",""))
    binary_sha=check_executable()
    selection,meta,selected=check_input(n)
    validate_native_csv_before_physics(selection)
    run_dir_=run_dir(n, retry=retry)
    if run_dir_.exists():raise FileExistsError(
        "Preserving prior B1 output (including interrupted): "+str(run_dir_))
    run_dir_.mkdir(parents=True)
    shutil.copy2(selection,run_dir_/"native_muon_crossings_input.csv")
    shutil.copy2(input_dir(n)/"input_manifest.json",run_dir_/"input_manifest_provenance.json")
    out_prefix=str(run_dir_/PREFIX)
    args=[str(EXE),str(run_dir_/"native_muon_crossings_input.csv"),out_prefix,"500"]
    audit={
        "schema":"PI_MUON_DRIVEN_TRUE_HEED_AVALANCHE_B1_SMOKE_V1",
        "status":"incomplete","timestamp_utc":datetime.now(timezone.utc).isoformat(),
        "channel":"Jpsi","n_true_geant4_crossings":n,"n_geometries":4,
        "geometry":[{"pitch_mm":3.0,"field_wire_um":0.0},
                    {"pitch_mm":3.0,"field_wire_um":75.0},
                    {"pitch_mm":4.0,"field_wire_um":0.0},
                    {"pitch_mm":4.0,"field_wire_um":75.0}],
        "cathode_gap_minus_mm":3.,"cathode_gap_plus_mm":3.,
        "gas":"Ar/CO2 70:30; 293.15 K; 760 Torr",
        "anode_diameter_um":30.,"anode_voltage_V":1800.,"B_T":0.,
        "analysis_unit":"one true GEANT4 MUON gas crossing; shared true TrackHeed primaries across all four geometries",
        "sample_kind":"signed-angle, station-stratified small SMOKE (not physical angle PDF)",
        "input_selector_manifest":meta,
        "git_commit":commit(),
        "b1_CXX_sha256":sha(CXX),"b1_binary_sha256":binary_sha,
        "input_csv_sha256":sha(run_dir_/"native_muon_crossings_input.csv")
    }
    (run_dir_/"run_config.json").write_text(json.dumps(audit,indent=2)+"\n")
    print(f"Running true Geant4→TrackHeed→AvalancheMicroscopic smoke, {n} native muons × four geometries.",flush=True)
    with (run_dir_/"physics.log").open("w") as logfile:
        result=subprocess.run(args,stdout=logfile,stderr=subprocess.STDOUT,
                              cwd=run_dir_,check=False)
    if result.returncode:
        logfile=run_dir_/"physics.log"
        tail="\n".join(logfile.read_text(errors="replace").splitlines()[-50:])
        raise RuntimeError(
            f"B1 failed with code {result.returncode}; original log: {logfile}\n"
            f"--- LAST 50 PHYSICS LOG LINES ---\n{tail}\n"
            "--- END LOG, PRESERVED WITHOUT OVERWRITE ---")
    if "PASS: muon-driven B1 all crossing/geometry rows written" not in (
        run_dir_/"physics.log").read_text():
        raise RuntimeError("B1 did not write completion marker")
    outputs=[run_dir_/(PREFIX+x) for x in ("_crossing.csv","_seed.csv","_ion_birth_sample.csv")]
    if not all(p.is_file() for p in outputs):
        raise RuntimeError("Incomplete B1 original per-crossing and per-primary output")
    stats=check_crossings(outputs[0],selected)
    if all(d["n_avalanche_ion_births"]==0 for d in stats.values()):
        raise RuntimeError("No avalanche ionization recorded for any geometry: inspect callback, gas field and seed locations")
    audit["status"]="complete"
    audit["stats"]=stats
    (run_dir_/"run_config.json").write_text(json.dumps(audit,indent=2)+"\n")
    (run_dir_/"summary.json").write_text(json.dumps(audit,indent=2)+"\n")
    audit["sha256_outputs"]={p.name:sha(p) for p in sorted(run_dir_.iterdir())
                              if p.is_file() and p.name!="manifest.json"}
    (run_dir_/"manifest.json").write_text(json.dumps(audit,indent=2)+"\n")
    print(f"PASS B1 completed SHA-manifested {n} ORIGINAL muon gas crossings in four candidate 3+3 geometries")
    for name,d in stats.items():
        print(name, "Heed primaries",d["n_heed_primaries"],
              "ion births",d["n_avalanche_ion_births"],
              "gain",d["G_eff_total_collected_over_total_Heed_primary"],
              "caps",d["n_capped_seeds"])
    print("Run directory:",run_dir_)
    print("WARNING: this is a smoke only; its angle strata DO NOT reproduce true muon incidence PDF.")

if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("mode",choices=("prepare2","prepare8","check2","check8","smoke2","smoke2retry1","smoke2retry2","smoke8","status"))
    args=p.parse_args()
    run(args.mode)
