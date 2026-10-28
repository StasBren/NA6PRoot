#!/usr/bin/env python3
"""Safe 3+3 mm, four-geometry Phase-A Garfield avalanche gain comparison.

Reuses test/mwpc-avalanche/phaseA_gain_fluctuations.cxx.
NO claim of an induced cathode signal or full-muon hit.
Each row is a SINGLE incoming electron seed from a fixed location.
"""
from __future__ import annotations
import csv
import hashlib
import json
import math
import re
import statistics
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

HERE=Path(__file__).resolve().parent
REPO=HERE.parents[1]
SRC=REPO/"test/mwpc-avalanche/phaseA_gain_fluctuations.cxx"
EXE=REPO/"test/mwpc-avalanche/build/mwpc_phase_a_gain"
OUT=REPO/"test_runs/mwpc_stage_c/avalanche_3p3_gain_v1"
CONFIGS=(
    ("p3_noFW",3.0,0.0),
    ("p3_FW75",3.0,75.0),
    ("p4_noFW",4.0,0.0),
    ("p4_FW75",4.0,75.0),
)
COLS=("event","avalanche_electrons","avalanche_ions","collected_on_wire",
      "attached","other_endpoints","mean_collection_time_ns")
COMMON={"gap_minus_mm":3.0,"gap_plus_mm":3.0,"anode_diam_um":30.0,
        "anode_voltage_V":1800.0,"B_T":0.0,"gas":"Ar/CO2 70:30",
        "T_K":293.15,"P_Torr":760.0,
        "seed_u_mm":0.0,"seed_v_mm":1.0,
        "base_seed":280103,"avalanche_limit":200000}
def sha(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda:f.read(1<<20),b""): h.update(chunk)
    return h.hexdigest()
def git_sha():
    return subprocess.check_output(["git","-C",str(REPO),"rev-parse","HEAD"],
                                   text=True).strip()
def preflight():
    if not SRC.is_file():raise FileNotFoundError(SRC)
    content=SRC.read_text()
    for required in ('"--field-wire-diam-um"','"--base-seed"',
                     '"--gap-minus-mm"','"field_wire"'):
        if required not in content:raise RuntimeError("Missing C++ FW/seed support: "+required)
    if not EXE.is_file():raise FileNotFoundError(
        f"Build mwpc_phase_a_gain first: {EXE}")
    if EXE.stat().st_mtime_ns < SRC.stat().st_mtime_ns:
        raise RuntimeError("Existing executable is OLDER than updated C++ source: rebuild mwpc_phase_a_gain")
    if not EXE.stat().st_mode & 0o111:
        raise RuntimeError("Garfield executable is not executable")
    return {"source_sha256":sha(SRC),"binary_sha256":sha(EXE),"git_commit":git_sha()}

def analyze(path,n_expected):
    with path.open(newline="") as f:
        rd=csv.DictReader(f)
        if tuple(rd.fieldnames or ())!=COLS:raise ValueError("Unexpected gain CSV schema")
        rows=list(rd)
    if len(rows)!=n_expected:raise ValueError("Wrong number of single-seed gain events")
    numbers={col:[] for col in COLS[1:]}
    for i,row in enumerate(rows):
        if int(row["event"])!=i:raise ValueError("Unordered/repeated seed ID")
        for col in COLS[1:]:
            v=float(row[col])
            if not math.isfinite(v) or v<0:raise ValueError("Bad event value in "+col)
            if col!="mean_collection_time_ns" and not v.is_integer():
                raise ValueError("Noninteger count: "+col)
            numbers[col].append(v)
    gains=numbers["collected_on_wire"]
    positive=[x for x in gains if x>0]
    return {"single_seed_events":len(rows),
            "effective_collected_gain_mean":statistics.mean(gains),
            "effective_collected_gain_rms":statistics.pstdev(gains),
            "mean_gain_given_positive":statistics.mean(positive) if positive else None,
            "zero_collected_events":len(gains)-len(positive),
            "zero_collected_fraction":1-len(positive)/len(gains),
            "mean_avalanche_ne":statistics.mean(numbers["avalanche_electrons"]),
            "mean_avalanche_ni":statistics.mean(numbers["avalanche_ions"]),
            "mean_attached":statistics.mean(numbers["attached"]),
            "nonsignal_endpoint_mean":statistics.mean(numbers["other_endpoints"]),
            "gain_p50":statistics.median(gains),
            "gain_p90":sorted(gains)[math.ceil(.9*len(gains))-1]}

def run(mode):
    if mode=="status":
        for name in ("smoke5","pilot100"):
            d=OUT/name
            print(f"{name}: manifest={'YES' if (d/'manifest.json').is_file() else 'NO'}; "
                  f"directory={'YES' if d.exists() else 'NO'}")
        return
    info=preflight()
    if mode=="check":
        print("PASS: freshly built four-geometry Garfield executable, source/commit and supported parameters")
        print(json.dumps(info,indent=2))
        return
    if mode not in ("smoke5","pilot100"):raise ValueError("mode must be check/smoke5/pilot100/status")
    events=5 if mode=="smoke5" else 100
    if mode=="pilot100" and not (OUT/"smoke5/manifest.json").is_file():
        raise RuntimeError("Review the 5-seed smoke before 100-seed pilot")
    dest=OUT/mode
    if dest.exists():raise FileExistsError(f"Preserve existing run: {dest}")
    dest.mkdir(parents=True) # An interrupted directory is intentionally not resumed automatically.
    report={"schema":"NA60_MWPC_AVALANCHE_P3P4_3PLUS3_GAIN_V1",
            "status":"incomplete","run_id":mode,
            "created_utc":datetime.now(timezone.utc).isoformat(),
            "physical_scope":"single electron seed ensemble; NOT Jpsi crossing/cathode current",
            "baseline":"p3_noFW","common":COMMON,
            "n_events_per_geometry":events,"build":info,
            "geometries":{}}
    (dest/"run_config.json").write_text(json.dumps(report,indent=2)+"\n")
    for name,pitch,fw in CONFIGS:
        out=dest/(name+".csv")
        log=dest/(name+".log")
        args=[str(EXE),
              "--pitch-mm",str(pitch),"--wire-diam-um","30",
              "--field-wire-diam-um",str(fw),
              "--gap-minus-mm","3","--gap-plus-mm","3",
              "--hv","1800","--b","0","--events",str(events),
              "--x0-mm","0","--y0-mm","1",
              "--base-seed",str(COMMON["base_seed"]),
              "--avalanche-limit",str(COMMON["avalanche_limit"]),
              "--progress-every","0","--output",str(out)]
        print(f"Running {name}: {events} microscopic single-electron avalanches.",flush=True)
        with log.open("w") as logf:
            completed=subprocess.run(args,cwd=dest,stdout=logf,stderr=subprocess.STDOUT)
        if completed.returncode!=0:
            raise RuntimeError(f"{name} failed (exit {completed.returncode}); inspect {log}")
        stats=analyze(out,events)
        text=log.read_text()
        m=re.search(r"size-limit-like evt\s*:\s*(\d+)",text)
        if m is None:raise RuntimeError("Missing size-limit audit: "+str(log))
        stats["size_limit_like_events"]=int(m.group(1))
        stats["field_wire_diameter_um"]=fw
        stats["pitch_mm"]=pitch
        stats["csv_sha256"]=sha(out);stats["log_sha256"]=sha(log)
        report["geometries"][name]=stats
        (dest/"run_config.json").write_text(json.dumps(report,indent=2)+"\n")
    report["status"]="complete"
    (dest/"summary.json").write_text(json.dumps(report,indent=2)+"\n")
    report["sha256_outputs"]={p.name:sha(p) for p in sorted(dest.iterdir())
                              if p.is_file() and p.name not in ("manifest.json","run_config.json")}
    (dest/"manifest.json").write_text(json.dumps(report,indent=2)+"\n")
    print("PASS: "+str(dest))
    for name,s in report["geometries"].items():
        print(f"{name}: <G_eff>={s['effective_collected_gain_mean']:.3f}, "
              f"zero={s['zero_collected_fraction']:.3f}, "
              f"limit-like={s['size_limit_like_events']}/{events}")
    if any(s["size_limit_like_events"] for s in report["geometries"].values()):
        print("WARNING: avalanche size caps reached; do not interpret mean gain as uncensored.")
    print("No induced cathode charges and no full Geant4/Heed muon-crossing simulation.")

if __name__=="__main__":
    if len(sys.argv)!=2:
        raise SystemExit("Usage: python pi_avalanche_3p3_gain_campaign.py check|smoke5|pilot100|status")
    run(sys.argv[1])
