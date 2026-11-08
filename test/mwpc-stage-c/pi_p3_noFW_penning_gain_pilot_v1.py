#!/usr/bin/env python3
"""Isolated Penning-ON vs frozen Penning-OFF p3_noFW gain control (real J/psi).

Does not modify/recompute any sealed 11k parent Garfield reference physics.
Reuses a FEW EXACT Geant4 crossings from sealed original-1000 chunk0000.
Runs native TrackHeed and full Garfield once in a SEPARATE Penning r=0.57
binary. Checks paired original Heed electrons/seeds and collects mean gas
gain before claiming any cause of the observed order-of-magnitude gap.

This tests ONLY the Penning effect at the original CURRENT 3-mm pitch,
symmetric +/-3-mm cathode planes and 1800 V anodes. It is NOT a
prototype-3 4-mm-pitch, 2+4-mm asymmetry test.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,subprocess,sys,datetime
from pathlib import Path
from collections import defaultdict
import numpy as np
from pi_garfield_runtime_env import garfield_runtime_env

REPO=Path(__file__).resolve().parents[2]
BASE=REPO/"test_runs/mwpc_stage_c"
OLD=BASE/"pi_p3_noFW_jpsi1000_native_v1/native_chunks/chunk_0000"
OUT=BASE/"pi_p3_noFW_jpsi_penning057_gain_pilot_v1"
ORIG_CPP=REPO/"test/mwpc-avalanche/phaseB_native_muon_p3_noFW_jpsi_runtime_v1.cxx"
PEN_CPP=REPO/"test/mwpc-avalanche/phaseB_native_muon_p3_noFW_jpsi_penning057_pilot_v1.cxx"
EXE=REPO/"test/mwpc-avalanche/build/mwpc_p3_noFW_jpsi_penning057_pilot_v1"
ID_FIELDS=("source","chunk","event","hit_index","track_id","station","det_id")
HEED_FIELDS=("crossing_id","primary_index","cluster_index","u_rel_cm","v_abs_cm",
             "w_rel_cm","t_ns")
SEED_FIELDS=("crossing_id","primary_index","garfield_seed","seed_u_mm","seed_v_mm",
             "garfield_ne","ion_births","collected_e","successful","capped")
INPUT="native_muon_crossings_input.csv"

def must(x,msg):
    if not x:raise RuntimeError(msg)
def sha(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda:f.read(1<<20),b""):h.update(b)
    return h.hexdigest()
def rows(p):
    with p.open(newline="") as f:
        r=csv.DictReader(f)
        must(r.fieldnames and len(set(r.fieldnames))==len(r.fieldnames),
             "Invalid/duplicate CSV columns "+str(p))
        return r.fieldnames,list(r)
def parentid(row):
    return "|".join([row[ID_FIELDS[0]]]+[str(int(row[k])) for k in ID_FIELDS[1:]])
def code_identity():
    base=ORIG_CPP.read_text()
    pen=PEN_CPP.read_text()
    needle="    gas.SetMaxElectronEnergy(200.);gas.Initialise(false);"
    must(base.count(needle)==1,"Frozen source changed; audit Penning-only diff")
    # The exact Penning-only pilot-source insertion must match the commit
    # created from the original one-geometry C++ source.
    added='''    gas.SetMaxElectronEnergy(200.);gas.Initialise(false);
    // CRITICAL PHYSICS TEST: Ar(70)/CO2(30) Penning excitation transfer.
    // Garfield++ Ar/CO2 standard r=0.57, lambda=0cm; unlike the frozen
    // reference 11k campaign, which intentionally remains UNCHANGED.
    gas.EnablePenningTransfer(0.57, 0., "ar");
    std::cout<<"PHYSICS_PENNING_AR_CO2_70_30_R=0.57 LAMBDA_CM=0\\n";'''
    hdr='''// PENNING-ON isolated physics gain PILOT: only difference from original
// 11k p3_noFW executable is enabled Ar Penning transfer r=0.57.
// Original SHA-sealed microscopic 11k campaign is NEVER overwritten.
'''
    must(pen==hdr+base.replace(needle,added),
         "Penning C++ differs from sealed reference beyond documented change")
    must(EXE.is_file(),"Separate Penning pilot binary is not built: "+str(EXE))
    return {"frozen_cpp_sha256":sha(ORIG_CPP),
            "penning_cpp_sha256":sha(PEN_CPP),
            "penning_executable_sha256":sha(EXE),
            "single_physics_change":"gas.EnablePenningTransfer(0.57,0.,\"ar\")"}
def selection(n):
    must(1<=n<=5,"Pilot must use 1..5 true muons for safe first test")
    seal_path=OLD/"complete.json"
    must(seal_path.is_file(),"Original frozen 1000 reference chunk not sealed")
    seal=json.loads(seal_path.read_text())
    must(seal["chunk"]==0 and seal["n_parents"]==20,
         "Original 20-muon baseline seal invalid")
    for name in (INPUT,"muon_reference_crossing.csv",
                 "muon_reference_seed.csv","muon_reference_heed_primaries.csv"):
        must(sha(OLD/name)==seal["file_sha256"][name],
             "Frozen Garfield baseline SHA mismatch: "+name)
    fields,inps=rows(OLD/INPUT)
    must(len(inps)==20,"Baseline input must have exactly 20 Geant4 parents")
    _,crossings=rows(OLD/"muon_reference_crossing.csv")
    data={r["source|chunk|event|hit_index|track_id|station|det_id"]:r for r in crossings}
    must(len(data)==20,"Original 20 real-muon crossing census incorrect")
    selected=[]
    # Protect pilot runtime: use muons with 25..60 true Heed electrons.
    for row in inps:
        key=parentid(row)
        s=data[key]
        count=int(s["heed_primaries"])
        if 25<=count<=60 and int(s["capped_seeds"])==0:
            selected.append(row)
            if len(selected)>=n:break
    must(len(selected)==n,"Not enough 25..60-electron real muons "
         "in sealed chunk0000; inspect physically, do NOT change inputs silently")
    return fields,selected,seal_path
def write_csv(p,fields,rr):
    with p.open("x",newline="") as f:
        wr=csv.DictWriter(f,fieldnames=fields,lineterminator="\n")
        wr.writeheader()
        wr.writerows(rr)
def means(vals):
    a=np.asarray(vals,dtype=np.float64)
    must(a.size>0,"Missing true Garfield electron counts")
    alive=a[a>0]
    return {"n":len(a),"n_zero":int(np.sum(a==0)),
            "p_zero":float(np.mean(a==0)),
            "mean_all":float(np.mean(a)),
            "mean_given_positive":float(np.mean(alive)) if alive.size else None,
            "RMS_all":float(np.std(a,ddof=0)),
            "max":int(np.max(a))}
def analyse(frozen,pen,n):
    oldbase=OLD/"muon_reference_seed.csv"
    newbase=OUT/"penning_seed.csv"
    _,oldrows=rows(oldbase)
    _,newrows=rows(newbase)
    target={parentid(r) for r in frozen}
    oldrows=[r for r in oldrows if r["crossing_id"] in target]
    newrows=[r for r in newrows if r["crossing_id"] in target]
    must(oldrows and newrows,"No original/replay real avalanche results")
    oldmap={(r["crossing_id"],r["primary_index"]):r for r in oldrows}
    newmap={(r["crossing_id"],r["primary_index"]):r for r in newrows}
    must(set(oldmap)==set(newmap),"Real Heed primary count/identity differs")
    must(len(oldmap)==sum(int(r["heed_primaries"]) for r in
             rows(OLD/"muon_reference_crossing.csv")[1]
             if r["source|chunk|event|hit_index|track_id|station|det_id"] in target),
         "Original full Heed electron count changed")
    for k,x in oldmap.items():
        y=newmap[k]
        for field in ("garfield_seed","seed_u_mm","seed_v_mm"):
            if field=="garfield_seed":
                same=int(x[field])==int(y[field])
            else:
                same=abs(float(x[field])-float(y[field]))<1e-9
            must(same,f"Not the same original starting electron {k}: {field}")
        must(int(y["capped"])==0,
             "Penning-on avalanche capped, cannot use pilot mean as unbiased")
    origHeed=rows(OLD/"muon_reference_heed_primaries.csv")[1]
    newHeed=rows(OUT/"penning_heed_primaries.csv")[1]
    origHeed=[r for r in origHeed if r["crossing_id"] in target]
    must(len(origHeed)==len(newHeed),
         "Penning changed TrackHeed primary census; not pairable")
    oldprim={(r["crossing_id"],r["primary_index"]):r for r in origHeed}
    newprim={(r["crossing_id"],r["primary_index"]):r for r in newHeed}
    must(set(oldprim)==set(newprim),
         "Penning changed TrackHeed electron indices; cannot pair")
    for k,x in oldprim.items():
        y=newprim[k]
        must(x["cluster_index"]==y["cluster_index"],
             "Original TrackHeed cluster index changed")
        for field in ("u_rel_cm","v_abs_cm","w_rel_cm","t_ns"):
            must(abs(float(x[field])-float(y[field]))<1e-9,
                 "Original TrackHeed starting coordinate/time changed")
    output={
        "schema":"P3_NOFW_PAIRED_REAL_JPSI_PENNING057_PILOT_V1",
        "n_genuine_muons":n,
        "n_exact_matched_TrackHeed_primaries":len(oldmap),
        "baseline":{"collected_e":means([int(x["collected_e"]) for x in oldmap.values()]),
                    "ion_births":means([int(x["ion_births"]) for x in oldmap.values()])},
        "Penning057":{"collected_e":means([int(y["collected_e"]) for y in newmap.values()]),
                      "ion_births":means([int(y["ion_births"]) for y in newmap.values()])},
        "paired_identical_primary_positions_and_seeds":True,
        "no_sealed_11k_reference_changed":True,
        "geometry":"3mm pitch; 30um wires; cathodes +/-3mm; HV=1800V; ArCO2 70:30; B=0",
        "caveat":"3-mm symmetric geometry is NOT prototype3 p4, 2+4mm chamber",
        "physical_interpretation":"PENNING parameter r=0.57 is a literature approximation, not validated by 2026 prototype gain curve",
    }
    g0=output["baseline"]["collected_e"]["mean_all"]
    g1=output["Penning057"]["collected_e"]["mean_all"]
    output["gain_ratio_penning_on_over_off"]=g1/g0 if g0>0 else None
    output["gain_uncapped"]=True
    return output

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("mode",choices=["preflight","run","status"])
    ap.add_argument("--muons",type=int,default=3)
    a=ap.parse_args()
    physics=code_identity()
    fields,chosen,seal_path=selection(a.muons)
    if a.mode=="preflight":
        print("PASS isolated Penning-only C++ source verified identical except r=0.57")
        print("Reference sealed chunk:",seal_path)
        print("selected true J/psi parents:",*[parentid(r) for r in chosen],sep="\n  ")
        print("NO source or physics output modified.")
        return
    if a.mode=="status":
        p=OUT/"comparison.json"
        must(p.is_file(),"Penning pilot not yet completed; preserve unsealed folder")
        print(p.read_text())
        return
    must(not OUT.exists(),"Refusing to overwrite previous Penning pilot: "+str(OUT))
    OUT.mkdir(parents=True)
    src=OUT/INPUT
    write_csv(src,fields,chosen)
    record={
        "source_sha256":sha(src),
        "old_chunk_complete_sha256":sha(seal_path),
        "physics_sha256":physics,
        "script_sha256":sha(Path(__file__)),
        "selected_parent_ids":[parentid(r) for r in chosen],
        "started_utc":datetime.datetime.now(datetime.timezone.utc).isoformat()
    }
    (OUT/"input_provenance.json").write_text(json.dumps(record,indent=2)+"\n")
    outprefix=OUT/"penning"
    with (OUT/"physics.log").open("x") as logfile:
        cp=subprocess.run(
            [str(EXE),"--muon-cache-physics-p3-noFW",str(src),str(outprefix)],
            cwd=OUT,env=garfield_runtime_env(),stdout=logfile,
            stderr=subprocess.STDOUT)
    raw=(OUT/"physics.log").read_text(errors="replace")
    must(cp.returncode==0 and
         "PHYSICS_PENNING_AR_CO2_70_30_R=0.57" in raw and
         "PASS: immutable same-Heed primary/cluster cache and full Garfield microscopic reference emitted; p3_noFW ONLY" in raw,
         "Full Penning Garfield pilot incomplete. Preserve all raw results: "+str(OUT))
    stats=analyse(chosen,OUT,a.muons)
    stats["input_provenance"]=record
    stats["pilot_file_sha256"]={
        p.name:sha(p) for p in OUT.iterdir()
        if p.is_file() and p.name!="comparison.json"}
    (OUT/"comparison.json").write_text(json.dumps(stats,indent=2)+"\n")
    print("PASS EXACT paired Penning ON/OFF initial J/psi Heed electrons")
    print(f"Penning OFF mean collected e/seed: {stats['baseline']['collected_e']['mean_all']:.3f}")
    print(f"Penning ON  mean collected e/seed: {stats['Penning057']['collected_e']['mean_all']:.3f}")
    print(f"Penning ON/OFF ratio: {stats['gain_ratio_penning_on_over_off']:.3f}")
    print(f"Penning ON zero fraction: {stats['Penning057']['collected_e']['p_zero']:.5f}")
    print("OUTPUT:",OUT/"comparison.json")
if __name__=="__main__":
    main()
