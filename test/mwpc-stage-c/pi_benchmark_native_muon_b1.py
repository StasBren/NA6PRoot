#!/usr/bin/env python3
"""Bounded timing of REAL TrackHeed electrons' microscopic Garfield avalanches.

BENCHMARK ONLY. Never use these partial crossing totals as a physical
measurement of the full-muon ion-birth asymmetry A_ion(mu).
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import platform
import resource
import statistics
import subprocess
import time
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

from pi_run_native_muon_b1 import (
    OUT, EXE, CXX, GEOMS, check_input, check_executable,
    validate_native_csv_before_physics, commit, sha,
)

def percent(values, pct):
    values=sorted(values)
    if not values:return None
    i=(len(values)-1)*pct
    lo=int(math.floor(i));hi=int(math.ceil(i))
    return values[lo]+(values[hi]-values[lo])*(i-lo)

def summarize(rows, original):
    by=defaultdict(list)
    count=defaultdict(set)
    per_muon=defaultdict(dict)
    primary_map={}
    for row in original:
        id_="|".join(str(row[k]) for k in (
            "source","chunk","event","hit_index","track_id","station","det_id"))
        primary_map[id_]=row
    for row in rows:
        geom=row["configuration"];id_=row["crossing_id"]
        if geom not in GEOMS or id_ not in primary_map:
            raise ValueError("Unexpected B1 timing configuration or muon crossing")
        index=int(row["primary_index"])
        key=id_,geom
        if index in count[key]:
            raise ValueError("Duplicate actual Heed electron index in timing results")
        count[key].add(index)
        sec=float(row["wall_seconds"])
        if not math.isfinite(sec) or sec<0:
            raise ValueError("Invalid per-seed time")
        for k in ("avalanche_electrons","avalanche_ions","ion_births","collected_e"):
            if int(row[k])<0:raise ValueError("Negative Garfield particle count")
        by[geom].append(sec)
        per_muon[geom][id_]=per_muon[geom].get(id_,0.)+sec
    if set(by)!=set(GEOMS):
        raise ValueError("All FOUR B1 geometries must have independent timing samples")
    n_per={k:len(v) for k,v in count.items()}
    if len(n_per)!=len(original)*4:
        raise ValueError("Missing genuine muon/geometry benchmark pair")
    expected={geom:{str(k):v for k,v in per_muon[geom].items()} for geom in GEOMS}
    summary={}
    for geom in GEOMS:
        samples=by[geom]
        s={
            "n_seed_benchmarks":len(samples),
            "mean_seconds_per_seed":statistics.mean(samples),
            "median_seconds_per_seed":statistics.median(samples),
            "p90_seconds_per_seed":percent(samples,.90),
            "min_seconds_per_seed":min(samples),
            "max_seconds_per_seed":max(samples),
            "std_seconds_per_seed":statistics.pstdev(samples),
            "sampled_muon_ids":sorted(per_muon[geom]),
            "seed_seconds_per_muon":expected[geom],
        }
        summary[geom]=s
    return summary

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--seeds",type=int,default=3,help="1..8 distinct real primary-electron indices per muon & geometry; default 3")
    a=p.parse_args()
    if not 1<=a.seeds<=8:p.error("--seeds must be between 1 and 8")
    # This is an independent, immutable diagnostic dataset. Never write
    # into a completed or failed B1 production directory.
    output=OUT/f"benchmark_true_Heed_k{a.seeds}_chunk0001_v1"
    if output.exists():raise FileExistsError("Preserving prior timing benchmark: "+str(output))
    binary_sha=check_executable()
    selection,original_meta,original=check_input(2)
    # This invokes the actual Garfield C++ parsing, ionisable-medium
    # check AND per-muon true TrackHeed generation (NO avalanche).
    validate_native_csv_before_physics(selection)
    output.mkdir(parents=True)
    csvout=output/"timings_per_primary.csv"
    log=output/"benchmark.log"
    manifest={
        "schema":"PI_TRUE_MUON_B1_MICROSCOPIC_BENCHMARK_V1",
        "status":"incomplete",
        "created_utc":datetime.now(timezone.utc).isoformat(),
        "geometry":list(GEOMS),"samples_per_muon_geometry":a.seeds,
        "n_physical_muon_crossings":len(original),
        "selection":"up to k quantile-spaced true Heed primaries along the REAL physical muon crossing",
        "physical_observable":"NONE: partial seed samples not valid for full-hit angular asymmetry",
        "compiled_cxx_sha256":sha(CXX),"compiled_binary_sha256":binary_sha,
        "original_muon_input_sha256":original_meta["input_sha256"],
        "original_geant4_manifest":original_meta,
        "git_commit":commit(),
        "host":{"system":platform.platform(),"python":platform.python_version(),
                "logical_cpu_count":os.cpu_count()},
        "units":"wall seconds per full microscopic single-primary-electron avalanche",
        "avalanche_cap":200000,
    }
    (output/"run_config.json").write_text(json.dumps(manifest,indent=2)+"\n")
    args=[str(EXE),"--benchmark",str(selection),str(csvout),str(a.seeds)]
    before=time.perf_counter()
    ru0=resource.getrusage(resource.RUSAGE_CHILDREN)
    with log.open("w") as stream:
        result=subprocess.run(args,cwd=output,stdout=stream,
                              stderr=subprocess.STDOUT,check=False)
    manifest["total_wall_seconds"]=time.perf_counter()-before
    ru1=resource.getrusage(resource.RUSAGE_CHILDREN)
    manifest["child_user_cpu_seconds"]=max(0.,ru1.ru_utime-ru0.ru_utime)
    manifest["child_system_cpu_seconds"]=max(0.,ru1.ru_stime-ru0.ru_stime)
    manifest["child_maxrss_kib"]=ru1.ru_maxrss
    manifest["return_code"]=result.returncode
    if result.returncode:
        (output/"run_config.json").write_text(json.dumps(manifest,indent=2)+"\n")
        raise RuntimeError("B1 timing failed; original files retained. Last log lines:\n"+
                           "\n".join(log.read_text(errors="replace").splitlines()[-30:]))
    with csvout.open(newline="") as f:
        reader=csv.DictReader(f)
        data=list(reader)
        expected_columns=("crossing_id","station","theta_u_deg","theta_w_deg",
            "configuration","primary_index","u_seed_mm","v_seed_mm",
            "garfield_seed","wall_seconds","avalanche_electrons",
            "avalanche_ions","ion_births","collected_e","success","capped")
        if tuple(reader.fieldnames or ())!=expected_columns:
            raise ValueError("Wrong benchmark timing CSV schema")
    summary=summarize(data,original)
    if any(int(x["capped"]) for x in data):
        manifest["warning"]="Size cap hit: measured times are censored; no unqualified extrapolation"
    if any(int(x["success"])!=1 for x in data):
        manifest["warning"]="At least one avalanche returned unsuccessful: inspect individual timing rows"
    manifest["status"]="complete"
    manifest["timing_by_geometry"]=summary
    manifest["n_timed_avalanches"]=len(data)
    # Stable files FIRST, hash them SECOND, and hash list in manifest LAST.
    # Never self-reference the manifest, or put stale hashes in summary.
    (output/"run_config.json").write_text(json.dumps(manifest,indent=2)+"\n")
    (output/"summary.json").write_text(json.dumps(manifest,indent=2)+"\n")
    manifest["sha256_outputs"]={x.name:sha(x) for x in sorted(output.iterdir())
                              if x.is_file() and x.name!="manifest.json"}
    (output/"manifest.json").write_text(json.dumps(manifest,indent=2)+"\n")
    print("PASS B1 real-muon avalanche BENCHMARK ONLY")
    print("  physical muon crossings:",len(original),"timed avalanches:",len(data))
    print("  entire process wall seconds:",f'{manifest["total_wall_seconds"]:.2f}')
    for g,s in summary.items():
        print(f"  {g:8s}: n={s['n_seed_benchmarks']}"
              f" mean={s['mean_seconds_per_seed']:.3f}s"
              f" median={s['median_seconds_per_seed']:.3f}s"
              f" p90={s['p90_seconds_per_seed']:.3f}s"
              f" max={s['max_seconds_per_seed']:.3f}s")
    print("  output:",output)
    print("IMPORTANT: k seed avalanches do NOT yield a physical per-muon A_ion; this is a timing diagnostic")

if __name__=="__main__":
    main()
