#!/usr/bin/env python3
"""ADD 10,000 disjoint genuine Geant4 J/psi muons to original 1000 p3_noFW.

No new Geant4 generation, no LUT surrogates: one authentic frozen gas chord
per unique generated J/psi parent -> real TrackHeed -> FULL microscopic
Garfield++ AvalancheElectron for ALL Heed electrons, same compiled p3_noFW
benchmark physics as the first 1000. Read-only original 1000 and J/psi source.

Uses the previously proven 20-parent chunk writer/strict native SHA auditor,
with a distinct immutable input/results root and source/parent exclusions.
No Garfield compilation or earlier source/binary modification required.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

from pi_2m_native_six_geometry_fold import SOURCE
from pi_prepare_muon_avalanche_b1 import FIELDS, seed, validate, sha
import pi_p3_noFW_jpsi1000_native_benchmark_v1 as base

# Change ONLY Python orchestration state in the known-validated benchmark.
# All actual source physics (TrackHeed/Garfield) remains in the SAME binary.
ROOT = base.REPO/"test_runs/mwpc_stage_c/pi_p3_noFW_jpsi_extra10000_native_v1"
OLD = base.REPO/"test_runs/mwpc_stage_c/pi_p3_noFW_jpsi1000_native_v1"
N = 10_000
GROUP = 20
TAG = "NA60-P3NOFW-JPSI-EXTRA10000-DISJOINT-PARENTS-20261105-V1"
BASE_SCHEMA = "NA60_JPSI_EXTRA10000_P3NOFW_BENCH_V1"
ROOT_FILES = ROOT/"input"
base.ROOT = ROOT
base.INPUT = ROOT_FILES
base.CHUNKS = ROOT/"native_chunks"
base.N = N
base.GROUP = GROUP
base.TAG = TAG


def ensure(condition, message):
    if not condition:
        raise ValueError(message)


def original_1000():
    p = OLD/"input/native_muon_crossings_input.csv"
    m = OLD/"input/input_manifest.json"
    ensure(p.is_file() and m.is_file(),
           "Previous complete 1000-muon immutable input missing: "+str(OLD))
    meta = json.loads(m.read_text())
    ensure(meta.get("n_parents")==1000 and sha(p)==meta.get("selected_input_sha256"),
           "Original 1000 source changed since successful SHA-sealed campaign")
    rows = list(base.csvread(p, FIELDS))
    ids = {base.parent_key(x) for x in rows}
    ensure(len(ids)==1000 and len(rows)==1000,
           "Previous 1000 source has duplicate parent identities")
    return ids, {"input_path":str(p.relative_to(base.REPO)),
                 "input_sha256":sha(p), "n_excluded_parent_events":1000}


def origin(chunk):
    """Only SHA-verified original 50k generated J/psi input chunks."""
    # Prior scripts supply this original provenance verification.
    return base.source(chunk)


def select_10000(start_chunk:int, end_chunk:int):
    exclude, excluded_proof = original_1000()
    ensure(3<=start_chunk<=end_chunk<=40,
           "Use frozen J/psi 50k-parent chunks in [3,40]; chunk2 is reserved")
    chosen = {}
    provenance=[]
    inspected=[]
    original_crossings=0
    for chunk in range(start_chunk,end_chunk+1):
        try:
            p,meta = origin(chunk)
        except FileNotFoundError:
            # Chunks that were not exported from the original 2M Geant4
            # simulation are NOT silently synthesized or substituted.
            inspected.append(chunk)
            continue
        provenance.append({
            "chunk":chunk,
            "csv_sha256":meta["csv_sha256"],
            "geant_root_sha256":meta["immutable_input"]["input_root_sha256"],
            "input_csv_relative_path":str(p.relative_to(base.REPO)),
        })
        neligible=0
        for r in base.csvread(p, set(FIELDS)-
                              {"source","chunk","momentum_GeV","heed_seed"}):
            ck=validate(r)
            if ck is None:
                continue
            momentum, vals=ck
            if abs(vals["u_out_mm"]-vals["u_in_mm"])>8. or \
               abs(vals["w_out_mm"]-vals["w_in_mm"])>40.:
                continue
            original_crossings += 1
            neligible += 1
            orig=(SOURCE["Jpsi"],chunk,r["event"],r["hit_index"],
                  r["track_id"],r["station"],r["det_id"])
            parent=(str(SOURCE["Jpsi"]),str(chunk),str(r["event"]))
            if parent in exclude:
                continue
            # One eligible gas crossing per original J/psi parent. Choose
            # independently of all subsequent Heed/Garfield outcomes.
            rank=hashlib.sha256(
                (TAG+"|crossing|"+"|".join(map(str,orig))).encode()
            ).hexdigest()
            old=chosen.get(parent)
            if old is None or rank<old[0]:
                record={k:r[k] for k in FIELDS
                        if k not in ("source","chunk","momentum_GeV","heed_seed")}
                record.update({
                    "source":SOURCE["Jpsi"],"chunk":chunk,
                    "momentum_GeV":f"{momentum:.12g}",
                    "heed_seed":seed(TAG,"heed",*orig),
                })
                chosen[parent]=(rank,record)
        print(f"Verified frozen J/psi chunk {chunk}: {neligible} eligible gas crossings; "
              f"{len(chosen)} new independent candidate parents",flush=True)
        # To avoid scanning millions of unnecessary original crossings,
        # freeze the exact set of chunks inspected once >= N are available.
        # The 10k selection is SHA-random among *that* fixed prefix.
        if len(chosen)>=N:
            break

    ensure(len(chosen)>=N,
           f"Only {len(chosen)} eligible NEW disjoint J/psi parents found "
           f"in available frozen chunks {start_chunk}..{end_chunk}; "
           "NO campaign was created. Export additional original Geant4 "
           "native crossing chunks (without new Geant4 simulation), then retry.")
    ranking=[
        (hashlib.sha256((TAG+"|parent|"+"|".join(parent)).encode()).hexdigest(),
         parent,record)
        for parent,(_,record) in chosen.items()
    ]
    ranking.sort()
    selected=[record for _,_,record in ranking[:N]]
    newids={base.parent_key(x) for x in selected}
    ensure(len(selected)==len(newids)==N and newids.isdisjoint(exclude),
           "Source has duplicate or original-1000-overlapping J/psi parents")
    stat={}
    for r in selected:
        k=str(r["station"])
        stat[k]=stat.get(k,0)+1
    selection_meta={
        "original_1000_exclusion":excluded_proof,
        "used_frozen_source_chunks":provenance,
        "skipped_missing_chunks":inspected,
        "n_eligible_crossings_scanned":original_crossings,
        "n_eligible_new_independent_parents":len(chosen),
        "selected_per_station":stat,
        "selection_rule":"SHA-random per independent Jpsi parent; one SHA-random "
                         "eligible gas crossing per parent; exclude all original "
                         "1000 parents; source chunks >=3 only",
        "limitations":"One crossing per eligible original parent, NOT unweighted "
                      "all-crossings Jpsi detector flux; no cathode simulation.",
    }
    return selected,selection_meta


def prepare(start_chunk:int,end_chunk:int):
    ensure(not ROOT.exists(),
           "Existing SHA-frozen extra10000 root. Never overwrite: "+str(ROOT))
    selected,selectmeta=select_10000(start_chunk,end_chunk)
    phys=base.check_binary()
    ROOT.mkdir(parents=True)
    base.INPUT.mkdir()
    output=base.INPUT/"native_muon_crossings_input.csv"
    base.csvwrite(output,selected)
    meta={
        "schema":BASE_SCHEMA,
        "n_parents":N,
        "group_size":GROUP,
        "geometry":"p3_noFW",
        "reservoir_per_avalanche":128,
        "selected_input_sha256":sha(output),
        "source":"Jpsi 2M original Geant4; frozen SHA-verified original chunks",
        "source_selection":selectmeta,
        "executable_provenance_prepared":phys,
        "campaign_program":str(Path(__file__).name),
    }
    base.jwrite(base.INPUT/"input_manifest.json",meta)
    check_input()
    print("PASS PREPARED 10000 NEW independent J/psi parent muons",flush=True)
    print("DISJOINT from the prior completed 1000; 500 x 20-parent SHA chunks")
    print("p3_noFW only; unchanged full microscopic Garfield native executable")
    print("NO Garfield physics in prepare")


def check_input():
    meta=base.jread(base.INPUT/"input_manifest.json")
    p=base.INPUT/"native_muon_crossings_input.csv"
    ensure(meta["schema"]==BASE_SCHEMA and meta["n_parents"]==N and
           meta["group_size"]==GROUP and meta["geometry"]=="p3_noFW" and
           sha(p)==meta["selected_input_sha256"],
           "Extra 10000 frozen selection/input SHA mismatch")
    oldids,oldproof=original_1000()
    ensure(oldproof==meta["source_selection"]["original_1000_exclusion"],
           "Original 1000 parent exclusion proof changed")
    for chunk in meta["source_selection"]["used_frozen_source_chunks"]:
        src,srcmeta=origin(chunk["chunk"])
        ensure(sha(src)==chunk["csv_sha256"]==srcmeta["csv_sha256"] and
               srcmeta["immutable_input"]["input_root_sha256"]==
                    chunk["geant_root_sha256"],
               "Original Geant4 input chunk changed: "+str(src))
    rr=list(base.csvread(p,FIELDS))
    identities={base.parent_key(x) for x in rr}
    ensure(len(rr)==len(identities)==N and identities.isdisjoint(oldids),
           "Wrong 10000 NEW independent parent identities")
    allowed={x["chunk"] for x in
             meta["source_selection"]["used_frozen_source_chunks"]}
    ensure(all(r["source"]==SOURCE["Jpsi"] and int(r["chunk"]) in allowed
               and validate(r) is not None for r in rr),
           "Invalid Geant4 parent source/chord/velocity")
    return rr,meta


# Reuse existing production auditor, runner, preflight and timing summaries.
# Their global config/lookups are monkeypatched only in this NEW Python module.
base.check_input = check_input


def status(audit=False):
    """Full 10k campaign status, no misleading 1000-muon inherited label."""
    results,pending=base.scanned_status(full_audit=audit)
    parents=sum(r["stats"]["n_parents"] for r in results)
    avalanches=sum(r["stats"]["n_full_garfield_avalanches"] for r in results)
    prim=sum(r["stats"]["n_trackheed_electrons"] for r in results)
    ion_births=sum(r["stats"]["n_ion_births"] for r in results)
    elapsed=sum(r["actual_native_wall_seconds"] for r in results)
    failed=sum(r["stats"]["n_unsuccessful_seeds"] for r in results)
    capped=sum(r["stats"]["n_capped_seeds"] for r in results)
    print("SHA-SEALED P3 NOFW EXTRA JPSI MUONS",parents,"/",N)
    print("SEALED CHUNKS",len(results),"/",N//GROUP)
    print("REAL MICROSCOPIC GARFIELD AVALANCHES",avalanches)
    print("TRACKHEED ELECTRONS",prim)
    print("FULL ION-BIRTH CENSUS",ion_births)
    print("ACTUAL NATIVE PHYSICS+IO WALL SECONDS",round(elapsed,2))
    print("SIZE-CAPPED CALLS",capped,"UNSUCCESSFUL CALLS",failed)
    if parents:
        mean=elapsed/parents
        print("MEASURED s/MUON",round(mean,4))
        print("MEASURED s/GARFIELD AVALANCHE",
              round(sum(r["stats"]["sum_seed_garfield_wall_seconds"]
                        for r in results)/max(1,avalanches),5))
        print("PROJECTED ADDITIONAL 10000 MUONS sequential HOURS",
              round(N*mean/3600,3))
        print("PROJECTED REMAINING sequential HOURS",
              round((N-parents)*mean/3600,3))
        muon_times=sorted(
            float(x["muon_total_wall_seconds"])
            for i in range(len(results))
            for x in base.csvread(base.folder(i)/"muon_reference_muon_timing.csv",
                                  ("muon_total_wall_seconds",))
        )
        ensure(len(muon_times)==parents,"Incomplete real-muon timing rows")
        def q(p):
            f=(parents-1)*p
            lo=math.floor(f)
            hi=math.ceil(f)
            return muon_times[lo]*(hi-f)+muon_times[hi]*(f-lo)
        print("TRUE MUON TIMES s p50/p90/p99/max",
              *(round(q(p),3) for p in (.5,.9,.99)),
              round(muon_times[-1],3))
    if pending is not None:
        print("UNSEALED LAST CHUNK",pending,
              "may be actively computing; do not overwrite/delete")
    elif len(results)==N//GROUP:
        print("PASS ALL 10000 NEW JPSI P3 NOFW MUONS SHA-SEALED")
    else:
        print("NEXT UNUSED CHUNK",len(results))
    if audit:
        print("PASS all completed native result file hashes independently re-audited")


def selftest():
    ensure(N==10000 and GROUP==20 and N//GROUP==500,
           "Wrong extra 10k target count/chunk plan")
    ensure(ROOT!=OLD and base.G=="p3_noFW" and base.RESERVOIR==128,
           "Prior full Garfield source or geometry changed")
    original,proof=original_1000()
    ensure(len(original)==1000,"Invalid old disjointness fence")
    sample={("Jpsi",str(3+i//1000),str(i)) for i in range(10000)}
    ensure(len(sample)==10000 and sample.isdisjoint(original),
           "Synthetic parent exclusion test failed")
    print("PASS disjoint selection, 10000 genuine J/psi parents, 500x20 plan")
    print("SAME pre-existing p3_noFW executable; no physics run in selftest")


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("mode",choices=("selftest","prepare","preflight","run",
                                   "status","audit"))
    p.add_argument("--start-chunk",type=int,default=3)
    p.add_argument("--end-chunk",type=int,default=40)
    p.add_argument("--max-chunks",type=int,default=1,
                   help="Number of NEW 20-muon chunks in this run")
    p.add_argument("--hours",type=float,default=20.,
                   help="Soft wall budget, checked BETWEEN completed 20-muon chunks")
    a=p.parse_args()
    ensure(a.max_chunks>0 and a.hours>0,"Invalid time/chunk budget")
    if a.mode=="selftest":
        selftest()
    elif a.mode=="prepare":
        prepare(a.start_chunk,a.end_chunk)
    elif a.mode=="preflight":
        base.preflight()
    elif a.mode=="run":
        base.run(a.max_chunks,a.hours)
    elif a.mode=="status":
        status()
    else:
        status(audit=True)


if __name__=="__main__":
    main()
