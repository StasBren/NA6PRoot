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
import fcntl
import re
import shutil
import subprocess
import time
from pathlib import Path

from pi_garfield_runtime_env import garfield_runtime_env

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

# Recover exceptional >1000-Heed-electron native events WITHOUT modifying the
# already SHA-sealed 426+ existing 20-parent chunks or the original binary.
# The failed attempt is PRESERVED BYTE-FOR-BYTE in chunk_XXXX/; the full
# scientific re-run lives in chunk_XXXX/high_heed_replay_v1/.
RECOVERY_TAG = "high_heed_replay_v1"
SPECIAL_CXX = base.REPO / (
    "test/mwpc-avalanche/"
    "phaseB_native_muon_p3_noFW_jpsi_highheed_recovery_v1.cxx"
)
SPECIAL_EXE = base.REPO / (
    "test/mwpc-avalanche/build/"
    "mwpc_phase_b_native_p3_noFW_jpsi_highheed_recovery_v1"
)
SPECIAL_MAX_PRIMARY = 50000
ORIGINAL_AUDIT_CHUNK = base.audit_chunk


def special_physics():
    ensure(SPECIAL_CXX.is_file() and SPECIAL_EXE.is_file(),
           "Build HIGH-Heed-only target before recover: " + str(SPECIAL_EXE))
    text_source = SPECIAL_CXX.read_text()
    ensure("highMultiplicityMaxPrimary=50000" in text_source
           and "RunGeometry(gas,GEOMS.front()" in text_source
           and "p3_noFW ONLY" in text_source,
           "Recovery C++ is not the reviewed p3_noFW real-muon code")
    return {"cpp_sha256":sha(SPECIAL_CXX),"exe_sha256":sha(SPECIAL_EXE)}


def recovery_chunk_audit(i, batch, phys, verify_files=True):
    """Accept normal frozen binary OR formally SHA-receipted rare-event replay."""
    d=base.folder(i)
    seal=d/"complete.json"
    if not seal.is_file():
        return ORIGINAL_AUDIT_CHUNK(i,batch,phys,verify_files=verify_files)
    m=base.jread(seal)
    if m.get("recovery_layout") != RECOVERY_TAG:
        return ORIGINAL_AUDIT_CHUNK(i,batch,phys,verify_files=verify_files)

    rec=d/RECOVERY_TAG
    ensure(rec.is_dir() and m["schema"]=="NA60_JPSI_20PARENT_P3NOFW_SHA_V1"
           and m["chunk"]==i and m["n_parents"]==len(batch)
           and m["base_physics"]==phys
           and m["physics"]==special_physics(),
           "High-Heed recovery physics provenance or 20-parent count changed")
    ensure(m["source_sha256"]==sha(d/base.FILES[0])
           and m["recovery_source_sha256"]==sha(rec/base.FILES[0]),
           "High-Heed recovered source identity SHA changed")
    ids=[base.id_of(x) for x in batch]
    for src in (d/base.FILES[0],rec/base.FILES[0]):
        ensure([base.id_of(x) for x in base.csvread(src,FIELDS)]==ids,
               "Replayed or original source rows differ from frozen J/psi input")
    ensure(m["original_partial_sha256"].keys() <=
           set(p.name for p in d.iterdir() if p.is_file()),
           "Part of original failed attempt has disappeared")
    if verify_files:
        ensure(set(m["file_sha256"])==set(base.FILES),
               "High-Heed recovered final physics file list invalid")
        for name,h in m["file_sha256"].items():
            ensure(sha(rec/name)==h,
                   "High-Heed recovered full physics SHA mismatch: "+str(rec/name))
        for name,h in m["original_partial_sha256"].items():
            ensure(sha(d/name)==h,
                   "Original failed chunk was modified: "+str(d/name))
    return m


# The old benchmark's own functions use globals of its module, so patching
# only audit_chunk is sufficient for BOTH its run() and scanned_status().
base.audit_chunk = recovery_chunk_audit


def diagnose_high_heed(i):
    """No avalanche calls and no output changes: exact TrackHeed counts."""
    rr,_=check_input()
    ensure(0<=i<N//GROUP,"Invalid 20-muon chunk index")
    d=base.folder(i)
    src=d/"native_muon_crossings_input.csv"
    ensure(src.is_file(),"No failed native chunk input: "+str(src))
    batch=rr[i*GROUP:(i+1)*GROUP]
    ensure([base.id_of(x) for x in base.csvread(src,FIELDS)]==
           [base.id_of(x) for x in batch],
           "Original failed native source does not match frozen selected 10000")
    env=garfield_runtime_env()
    cp=subprocess.run(
        [str(base.EXE),"--validate-heed",str(src)],
        capture_output=True,text=True,env=env)
    print(cp.stdout,flush=True)
    if cp.stderr:
        print("STDERR:",cp.stderr[-2000:],flush=True)
    ensure(cp.returncode==0 and
           "PASS: true TrackHeed source generated" in cp.stdout,
           "Heed-only diagnostic failed; NO Garfield replay")
    observed=[]
    for row in cp.stdout.splitlines():
        m=re.search(r"Heed native crossing event=(\d+).*?\belectrons=(\d+)",row)
        if m:
            observed.append((int(m.group(1)),int(m.group(2))))
    ensure(len(observed)==GROUP,"TrackHeed diagnostic did not cover all 20 parents")
    peak=max(n for _,n in observed)
    print(f"TRUE Heed high-multiplicity diagnostic: max={peak}, "
          f"n_over_1000={sum(n>1000 for _,n in observed)} "
          f"full microscopic avalanches needed={sum(n for _,n in observed)}",
          flush=True)
    return peak,observed


def recover_high_heed(i):
    """Complete one failed chunk via separate higher-budget C++ binary, safely."""
    rr,_=check_input()
    base_phys=base.checked_physics()
    special_phys=special_physics()
    ensure(0<=i<N//GROUP,"Invalid recovery 20-parent chunk index")
    base.CHUNKS.mkdir(parents=True,exist_ok=True)
    lock=base.CHUNKS/".writer.lock"
    with lock.open("a+") as fh:
        try:
            fcntl.flock(fh,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:
            raise RuntimeError("Normal or recovery Garfield writer is ACTIVE")
        batches=base.chunks(rr)
        for k in range(i):
            ensure((base.folder(k)/"complete.json").is_file(),
                   "Earlier chunk is not already fully SHA-sealed: "+str(k))
            recovery_chunk_audit(k,batches[k],base_phys,verify_files=False)
        d=base.folder(i)
        ensure(d.is_dir() and not (d/"complete.json").exists(),
               "Recovery requires an existing UNSEALED failed 20-parent chunk")
        log=d/"physics.log"
        ensure(log.is_file() and "Heed primary budget exceeded" in
               log.read_text(errors="replace"),
               "Recovery is permitted ONLY for a verified 1000-Heed budget failure")
        src=d/"native_muon_crossings_input.csv"
        ensure([base.id_of(x) for x in base.csvread(src,FIELDS)]==
               [base.id_of(x) for x in batches[i]],
               "Failed chunk native input not identical to frozen source")
        peak,observed=diagnose_high_heed(i)
        ensure(1000<peak<=SPECIAL_MAX_PRIMARY,
               "Heed primary count outside validated recovery range; "
               f"max={peak}. Do NOT truncate or skip this physical J/psi event.")
        recovered=d/RECOVERY_TAG
        ensure(not recovered.exists(),
               "A recovery attempt already exists: "+str(recovered)+
               ". Do not overwrite any potentially partial physics.")
        # Record every original partial byte so a future audit also proves
        # the interrupted 35-avalanche event was NOT overwritten/deleted.
        old_hash={
            p.name:sha(p) for p in sorted(d.iterdir()) if p.is_file()
            and p.name!="complete.json"
        }
        recovered.mkdir()
        r_input=recovered/"native_muon_crossings_input.csv"
        shutil.copyfile(src,r_input)
        ensure(sha(src)==sha(r_input),"Frozen source copy differs")
        env=garfield_runtime_env()
        print(f"REPLAY full 20 parents p3_noFW (highest Heed={peak}) "
              f"with physically identical special binary; old partial "
              f"remains at {d}",flush=True)
        t0=time.monotonic()
        with (recovered/"physics.log").open("x") as f:
            cp=subprocess.run(
                [str(SPECIAL_EXE),"--muon-cache-physics-p3-noFW",
                 str(r_input),str(recovered/"muon_reference")],
                stdout=f,stderr=subprocess.STDOUT,cwd=recovered,env=env)
        wall=time.monotonic()-t0
        stdout=(recovered/"physics.log").read_text(errors="replace")
        ensure(cp.returncode==0 and
               "PASS: immutable same-Heed primary/cluster cache and full "
               "Garfield microscopic reference emitted; p3_noFW ONLY" in stdout,
               "Recovery full Garfield incomplete: preserve "+str(recovered))
        stats=base.validate_native(recovered,batches[i])
        # Replay the original first completed native muon and all its 35
        # original electron avalanches exactly; no silent physics changes.
        partial_hit=list(base.csvread(d/"muon_reference_crossing.csv",
                                      CROSSING_COLUMNS))
        full_hit=list(base.csvread(recovered/"muon_reference_crossing.csv",
                                   CROSSING_COLUMNS))
        ensure(len(full_hit)==GROUP
               and full_hit[:len(partial_hit)]==partial_hit,
               "Recovered original completed muon counts differ: "
               "possible physics or RNG change; NOT SHA-sealing")
        partial_seeds=list(base.csvread(d/"muon_reference_seed.csv",
                                        ("crossing_id","primary_index")))
        replay_seeds=list(base.csvread(recovered/"muon_reference_seed.csv",
                                       ("crossing_id","primary_index")))
        ensure(replay_seeds[:len(partial_seeds)]==partial_seeds,
               "Recovered partial Garfield avalanche counts differ: "
               "possible RNG change; NOT SHA-sealing")
        source_hash=sha(src)
        receipt={
            "schema":"NA60_JPSI_20PARENT_P3NOFW_SHA_V1",
            "chunk":i,"n_parents":len(batches[i]),
            "recovery_layout":RECOVERY_TAG,
            "base_physics":base_phys,
            "physics":special_phys,
            "source_sha256":source_hash,
            "recovery_source_sha256":sha(r_input),
            "original_partial_sha256":old_hash,
            "original_first_completed_crossings_replayed_identically":len(partial_hit),
            "original_partial_seeds_replayed_identically":len(partial_seeds),
            "Heed_max_primaries_in_one_muon":peak,
            "Heed_primaries_in_all_20":sum(n for _,n in observed),
            "actual_native_wall_seconds":wall,
            "stats":stats,
            "file_sha256":{name:sha(recovered/name) for name in base.FILES},
        }
        # Publish the chunk as complete ONLY AFTER checking all real
        # Geant4/Heed/Garfield conservation and physical replay identity.
        tmp=d/"complete.json.tmp"
        base.jwrite(tmp,receipt)
        tmp.replace(d/"complete.json")
        recovery_chunk_audit(i,batches[i],base_phys,verify_files=True)
        print(f"PASS HIGH-HEED RECOVERY SHA chunk {i}: "
              f"{stats['n_parents']} true muons, "
              f"{stats['n_full_garfield_avalanches']} full avalanches, "
              f"{stats['n_ion_births']} ion births, "
              f"max Heed {peak}, wall {wall:.1f} s",flush=True)
        print("Original partial native physics remains unmodified; "
              "resume ordinary campaign with run, status, audit",flush=True)


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
    sample={(str(SOURCE["Jpsi"]),str(30+i//1000),str(i))
            for i in range(10000)}
    ensure(len(sample)==10000 and sample.isdisjoint(original),
           "Synthetic parent exclusion test failed")
    print("PASS disjoint selection, 10000 genuine J/psi parents, 500x20 plan")
    print("SAME pre-existing p3_noFW executable; no physics run in selftest")


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("mode",choices=("selftest","prepare","preflight","run",
                                   "status","audit","diagnose-high-heed",
                                   "recover-high-heed"))
    p.add_argument("--chunk",type=int,default=426,
                   help="Exact failed 20-muon chunk index for high-Heed recovery")
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
    elif a.mode=="diagnose-high-heed":
        diagnose_high_heed(a.chunk)
    elif a.mode=="recover-high-heed":
        recover_high_heed(a.chunk)
    else:
        status(audit=True)


if __name__=="__main__":
    main()
