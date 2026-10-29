#!/usr/bin/env python3
"""Real Geant4 muons -> one cached TrackHeed population -> full Garfield
versus stochastic frozen-LUT500 avalanche response.

Modes: prepare20 | check20 | reference20 | compare20 | publish20 | status.
Reads original immutable Geant4 source and previously published 84k LUT
WITHOUT altering either. No random muon directions are fabricated.
Full Garfield reference has ONE stochastic realization per Heed primary.
Fast LUT has 200 INDEPENDENT draws per same frozen Heed population.
A pair is matched by primaries, NOT by microscopic avalanche RNG.
This first pass validates ion birth counts and angular sectors only;
neither ion drift/induced cathode strip signal nor timing emulation is claimed.
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
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path

from pi_prepare_muon_avalanche_b1 import FIELDS,OUT,make,sha
from pi_run_native_muon_b1 import (
    CXX,EXE,GEOMS,CROSSING_COLUMNS,check_executable,
    check_crossings,validate_native_csv_before_physics,commit,
)
from pi_garfield_runtime_env import garfield_runtime_env
from pi_audit_avalanche_lut_reuse import BASE,CONFIGS
from pi_publish_run_to_drive import prepare as drive_prepare,publish,EXPECTED_REMOTE

ROOT=BASE/"pi_muon_lut_validation_v1"
INPUT=ROOT/"input_20_native_geant4_v1"
REFERENCE=ROOT/"reference_20_muon_full_garfield_v1"
COMPARISON=ROOT/"compare_20_muon_lut500_v1"
OVERNIGHT=BASE/"pi_lut500_avalanche_overnight_v1"/"run_20percell_firstpass_v1"
N_MUONS=20
REPLICATES=200
CROSSING_KEY=CROSSING_COLUMNS[0]
PRIMARY_FIELDS=(
    "crossing_id","primary_index","cluster_index",
    "u_rel_cm","v_abs_cm","w_rel_cm","t_ns",
    "raw_heed_energy_eV","raw_dx","raw_dy","raw_dz")
CLUSTER_FIELDS=(
    "crossing_id","cluster_index","u_rel_cm","v_abs_cm",
    "w_rel_cm","t_ns","n_original","n_accepted")
SEED_FIELDS=("crossing_id","configuration","primary_index","garfield_seed",
  "seed_u_mm","seed_v_mm","seed_w_mm","seed_energy_eV",
  "garfield_ne","garfield_ni","ion_births","collected_e",
  "attached","other","successful","capped",
  "ion_up","ion_down","ion_top","ion_bottom","A_ion_hemi","A_ion_TB")
LUT_FIELDS=("lut_trials","original_geometry","original_lut_sha256","box_id",
    "iu","iv","sample","original_random_seed","original_u_start_mm",
    "original_v_start_mm","original_w_start_mm","old_collected",
    "old_field_wire_collected","old_wire_index","old_phi_hit_deg",
    "old_drift_time_ns","old_last_step_energy_eV",
    "crossing_id","configuration","primary_index","garfield_seed",
    "seed_u_mm","seed_v_mm","seed_w_mm","seed_energy_eV",
    "garfield_ne","garfield_ni","ion_births","collected_e",
    "attached","other","successful","capped",
    "ion_up","ion_down","ion_top","ion_bottom","A_ion_hemi","A_ion_TB")
N_U={"p3_noFW":30,"p3_FW75":30,"p4_noFW":40,"p4_FW75":40}

def stable_seed(*terms):
    payload="\x1f".join(map(str,terms)).encode()
    return int.from_bytes(hashlib.sha256(payload).digest()[:8],"big")

def write_json(path,obj):
    path.write_text(json.dumps(obj,sort_keys=True,indent=2)+"\n")

def read_csv(path,fields):
    with path.open(newline="") as f:
        rd=csv.DictReader(f)
        if tuple(rd.fieldnames or ())!=fields:
            raise ValueError(f"CSV schema mismatch: {path}")
        return list(rd)

def id_of(row):
    return "|".join(str(row[key]) for key in
        ("source","chunk","event","hit_index","track_id","station","det_id"))

def percent(sorted_values,q):
    if not sorted_values:return None
    x=(len(sorted_values)-1)*q
    lo=math.floor(x);hi=math.ceil(x)
    return sorted_values[lo]*(hi-x)+sorted_values[hi]*(x-lo) if lo!=hi else sorted_values[lo]

def round_away(v):
    return math.floor(v+.5) if v>=0 else math.ceil(v-.5)

def entry_phase_mm(entry_u_mm,station,pitch_mm):
    # EXACT legacy C++ NAnodes/EntrancePhase, including wire-grid parity.
    width=307.2 if station==0 else 512.
    n=math.floor(width/pitch_mm+.5)
    phase=pitch_mm*.5 if n%2==0 else 0.
    nearest=round_away((entry_u_mm-phase)/pitch_mm)*pitch_mm+phase
    return entry_u_mm-nearest

def spatial_cell(geom,station,entry_u_mm,heed_u_rel_cm,heed_v_abs_cm):
    pitch=3. if geom.startswith("p3_") else 4.
    phase=entry_phase_mm(entry_u_mm,station,pitch)
    u_mm=10.*heed_u_rel_cm+phase
    local_u=u_mm-round_away(u_mm/pitch)*pitch
    v_mm=10.*heed_v_abs_cm
    if not (math.isfinite(local_u) and math.isfinite(v_mm)
            and -3.<=v_mm<=3.):
        raise ValueError("Nonphysical real Heed primary for frozen LUT500 gas span")
    # Source LUT cells are [−pitch/2,+pitch/2] in 0.1mm u bins
    # and [−3,+3] in 0.2mm v bins, with nearest wire wrapping.
    nu=N_U[geom]
    iu=min(nu-1,max(0,int(math.floor((local_u+pitch/2)/.1))))
    iv=min(29,max(0,int(math.floor((v_mm+3.)/.2))))
    if not (0<=iu<nu and 0<=iv<30):
        raise ValueError("Invalid native Heed->LUT500 bin")
    return iu*30+iv,local_u,v_mm

def input_verified():
    f=INPUT/"native_muon_crossings_input.csv"
    m=INPUT/"input_manifest.json"
    if not f.is_file() or not m.is_file():
        raise FileNotFoundError("Run prepare20 (native Geant4 source selection) first")
    meta=json.loads(m.read_text())
    if (meta.get("schema")!="PI_TRUE_GEANT4_MUON_B1_INPUT_V1"
        or meta.get("status")!="complete"
        or meta.get("n_crossings")!=N_MUONS
        or meta.get("input_sha256")!=sha(f)):
        raise ValueError("Immutable 20 real Geant4 input selection changed")
    rows=read_csv(f,FIELDS)
    if len(rows)!=N_MUONS or len({id_of(x) for x in rows})!=N_MUONS:
        raise ValueError("Not twenty distinct native muon gas crossings")
    return f,meta,rows

def run_reference():
    binary=check_executable()
    inp,meta,selected=input_verified()
    validate_native_csv_before_physics(inp)
    if REFERENCE.exists():
        raise FileExistsError("Preserving original full muon Garfield run: "+str(REFERENCE))
    REFERENCE.mkdir(parents=True)
    shutil.copy2(inp,REFERENCE/inp.name)
    shutil.copy2(INPUT/"input_manifest.json",REFERENCE/"input_manifest_provenance.json")
    cfg={
        "schema":"NA60_MUON_HEED_COMMON_PRIMARIES_FULL_GARFIELD_REFERENCE_V1",
        "status":"incomplete","n_true_geant4_crossings":N_MUONS,
        "native_input_sha256":meta["input_sha256"],
        "source_geant4_raw_sha256":meta["original_native_sha256"],
        "input_root_sha256":meta["original_root_sha256"],
        "cpp_source_sha256":sha(CXX),"cpp_executable_sha256":binary,
        "git_commit":commit(),"n_geom":4,
        "gas":"Ar/CO2 70:30, 293.15 K, 760 Torr",
        "anode_voltage_V":1800,"anode_diameter_um":30,
        "cathode_gap_mm":3,"magnetic_field_T":0,
        "physics":"One native Geant4 crossing -> TrackHeed cluster cache ONCE -> every accepted Heed electron as exactly one full Garfield AvalancheElectron in all four fields. 0.1eV thermal isotropic initial direction; ion position reservoir 5; full per-seed ion sector census.",
    }
    write_json(REFERENCE/"run_config.json",cfg)
    prefix=REFERENCE/"muon_reference"
    with (REFERENCE/"physics.log").open("w") as log:
        result=subprocess.run(
            [str(EXE),"--muon-cache-physics",
             str(REFERENCE/inp.name),str(prefix)],
             cwd=REFERENCE,stdout=log,stderr=subprocess.STDOUT,
             env=garfield_runtime_env(),check=False)
    if result.returncode:
        raise RuntimeError(f"Full physical Garfield failed (raw original files untouched): "
            f"{REFERENCE/'physics.log'}\n" +
            "\n".join((REFERENCE/"physics.log").read_text(errors="replace").splitlines()[-30:]))
    if "PASS: immutable same-Heed primary/cluster cache" not in (
            REFERENCE/"physics.log").read_text():
        raise RuntimeError("Missing C++ real-Heed-cache completion marker")
    checks=validate_reference(REFERENCE,selected)
    cfg["status"]="complete";cfg["stats"]=checks
    write_json(REFERENCE/"summary.json",cfg)
    files={p.name:{"sha256":sha(p),"bytes":p.stat().st_size}
           for p in sorted(REFERENCE.iterdir()) if p.is_file()
           and p.name!="manifest.json"}
    write_json(REFERENCE/"manifest.json",{**cfg,"outputs":files})
    print("PASS real Geant4×Heed×full Garfield reference for 20 fixed "
          "muon crossings in four geometries")
    for geom,stat in checks["by_geometry"].items():
        print("  ",geom,stat,flush=True)

def validate_reference(ref,selected):
    src={id_of(x):x for x in selected}
    primary=read_csv(ref/"muon_reference_heed_primaries.csv",PRIMARY_FIELDS)
    clusters=read_csv(ref/"muon_reference_heed_clusters.csv",CLUSTER_FIELDS)
    seeds=read_csv(ref/"muon_reference_seed.csv",SEED_FIELDS)
    hits=read_csv(ref/"muon_reference_crossing.csv",CROSSING_COLUMNS)
    stats=check_crossings(ref/"muon_reference_crossing.csv",selected)
    nprimary=Counter();ncluster=Counter();nkeep=Counter();pn=set();cln=set()
    for r in primary:
        x=r["crossing_id"]
        j=int(r["primary_index"])
        ci=int(r["cluster_index"])
        if x not in src or j<0 or (x,j) in pn:
            raise ValueError("Cached primary has unknown/duplicate native muon ID")
        pn.add((x,j));nprimary[x]+=1
        if ci<0:raise ValueError("Heed electron lacks accepted cluster identity")
        nkeep[x,ci]+=1
        if not all(math.isfinite(float(r[k])) for k in
                   ("u_rel_cm","v_abs_cm","w_rel_cm","t_ns")):
            raise ValueError("Nonfinite Heed coordinate/time")
    for x,n in nprimary.items():
        if {i for key,i in pn if key==x}!=set(range(n)):
            raise ValueError("Discontinuous cached Heed electron index")
    for r in clusters:
        x=r["crossing_id"];i=int(r["cluster_index"])
        if x not in src or (x,i) in cln or i<0:
            raise ValueError("Cached cluster malformed/duplicate")
        cln.add((x,i));ncluster[x]+=1
        if int(r["n_accepted"])!=nkeep[x,i]:
            raise ValueError("Cached cluster-to-electron membership count mismatch")
        if int(r["n_original"])<int(r["n_accepted"]):
            raise ValueError("Heed cluster has more electrons than created")
    for x in src:
        if {i for z,i in cln if z==x}!=set(range(ncluster[x])):
            raise ValueError("Heed cluster index lost")
    hit_lookup={(r[CROSSING_KEY],r["configuration"]):r for r in hits}
    if len(hit_lookup)!=len(selected)*4:
        raise ValueError("Missing full physical muon geometry rows")
    for (x,g),r in hit_lookup.items():
        if (int(r["heed_primaries"])!=nprimary[x] or
            int(r["heed_clusters"])!=ncluster[x]):
            raise ValueError("Heed cluster/primary cache is not the full reference input")
    seen=set();sums=defaultdict(lambda:Counter())
    for row in seeds:
        x=row["crossing_id"];geom=row["configuration"]
        idx=int(row["primary_index"]);key=(x,geom,idx)
        if x not in src or geom not in GEOMS or key in seen or (x,idx) not in pn:
            raise ValueError("Unknown/duplicate native full Garfield primary seed")
        seen.add(key)
        if int(row["capped"]) or int(row["successful"])!=1:
            raise RuntimeError("Capped/failed full Garfield microscopic muon reference")
        ion=int(row["ion_births"])
        a,b=int(row["ion_up"]),int(row["ion_down"])
        top,bot=int(row["ion_top"]),int(row["ion_bottom"])
        if ion<0 or a+b>ion or top+bot>ion:
            raise ValueError("Full Garfield per-seed ion angular conservation")
        s=sums[(x,geom)]
        for key_,val in (("ions",ion),("up",a),("down",b),
            ("top",top),("bottom",bot),("nonzero",int(ion>0))):
            s[key_]+=val
        s["n"]+=1
    for (x,g),r in hit_lookup.items():
        if (sums[(x,g)]["n"]!=nprimary[x] or
            sums[(x,g)]["ions"]!=int(r["ion_births"])):
            raise ValueError("Per-Heed-primary Garfield births do not sum to muon")
    if len(seen)!=sum(nprimary.values())*4:
        raise ValueError("Missing any per-Heed-primary reference Garfield avalanche")
    if sum(nprimary.values())==0:
        raise ValueError("No genuine Heed primary ionization in sample")
    return {
        "n_native_muons":len(selected),
        "n_heed_clusters":sum(ncluster.values()),
        "n_heed_primaries":sum(nprimary.values()),
        "n_full_avalanche_trials":len(seeds),
        "by_geometry":stats,
    }

def verify_ref():
    inp,meta,selected=input_verified()
    if not (REFERENCE/"manifest.json").is_file():
        raise FileNotFoundError("First run reference20; no original full muon reference")
    manifest=json.loads((REFERENCE/"manifest.json").read_text())
    if manifest.get("status")!="complete":
        raise ValueError("Full Garfield reference incomplete")
    for name,record in manifest.get("outputs",{}).items():
        if sha(REFERENCE/name)!=record["sha256"]:
            raise ValueError("Original full Garfield reference SHA changed: "+name)
    return selected,validate_reference(REFERENCE,selected),sha(REFERENCE/"manifest.json")

def load_lut():
    """Verify all original overnight full 84k single-electron trials.

    No use of OLD drift 'collected' flag: ONLY independent new avalanches.
    Each native Heed location is coarse-binned to one 0.1×0.2mm box.
    """
    man_file=OVERNIGHT/"manifest.json"
    if not man_file.is_file():
        raise FileNotFoundError("Missing complete overnight original LUT500 physics")
    man=json.loads(man_file.read_text())
    if (man.get("schema")!="PI_LUT500_OVERNIGHT_20_PER_CELL_V1"
        or man.get("status")!="complete"
        or man.get("n_completed_avalanches")!=84000
        or man.get("total_capped")!=0 or man.get("total_failed")!=0):
        raise ValueError("Original 84k LUT500 not a complete error-free physics run")
    out={}
    for geom in GEOMS:
        name=f"trials_{geom}.csv"
        path=OVERNIGHT/name
        if sha(path)!=man["files_sha256"][name]:
            raise ValueError("Frozen independent Garfield LUT500 SHA changed: "+name)
        allrows=read_csv(path,LUT_FIELDS)
        expected=N_U[geom]*30*20
        if len(allrows)!=expected:
            raise ValueError(f"Wrong 20-per-cell LUT500 coverage in {geom}")
        cells=defaultdict(list);used=set()
        for r in allrows:
            box=int(r["box_id"]);s=int(r["sample"])
            key=(box,s)
            if key in used or not 0<=box<N_U[geom]*30 or not 0<=s<500:
                raise ValueError("Duplicate/out-of-range original LUT500 trial")
            used.add(key)
            if (r["original_geometry"]!=geom or r["configuration"]!=geom
                or int(r["lut_trials"])!=500
                or r["original_lut_sha256"]!=man["old_lut500_csv_sha256"][geom]
                or int(r["capped"])!=0 or int(r["successful"])!=1):
                raise ValueError("Unmatched original complete physical LUT500 provenance")
            ni=int(r["ion_births"]);a=int(r["ion_up"]);b=int(r["ion_down"])
            top=int(r["ion_top"]);bot=int(r["ion_bottom"])
            if ni<0 or a+b>ni or top+bot>ni:
                raise ValueError("Bad full microscopic ion census in LUT500")
            cells[box].append((ni,a,b,top,bot))
        if len(cells)!=N_U[geom]*30 or any(len(x)!=20 for x in cells.values()):
            raise ValueError("LUT500 has incomplete cell sampling")
        out[geom]=cells
        print("PASS frozen LUT500",geom,expected,
              "independent full-start Garfield trials, 20 per cell",flush=True)
    return out,sha(man_file)

def summarize_samples(vals):
    vals=sorted(vals)
    return {
        "mean":statistics.mean(vals),
        "sd":statistics.pstdev(vals),
        "p05":percent(vals,.05),
        "p50":percent(vals,.5),
        "p95":percent(vals,.95),
    }

def asym(a,b):
    return (a-b)/(a+b) if a+b else None

def compare(replicates):
    if not 10<=replicates<=1000:raise ValueError("Use 10..1000 replicas")
    selected,ref_stats,ref_hash=verify_ref()
    lut,lut_hash=load_lut()
    if COMPARISON.exists():
        raise FileExistsError("Immutable LUT-vs-true-muon comparison exists: "+str(COMPARISON))
    COMPARISON.mkdir(parents=True)
    original={id_of(x):x for x in selected}
    primaries=defaultdict(list)
    for p in read_csv(REFERENCE/"muon_reference_heed_primaries.csv",PRIMARY_FIELDS):
        primaries[p["crossing_id"]].append(p)
    refhits={(r[CROSSING_KEY],r["configuration"]):r
         for r in read_csv(REFERENCE/"muon_reference_crossing.csv",CROSSING_COLUMNS)}
    ref_seeds=defaultdict(list)
    for r in read_csv(REFERENCE/"muon_reference_seed.csv",SEED_FIELDS):
        ref_seeds[(r["crossing_id"],r["configuration"])].append(r)
    report=[]
    ensemble=[]
    for crossing in selected:
        x=id_of(crossing)
        for geom in GEOMS:
            pitch=3. if geom.startswith("p3_") else 4.
            positions=[]
            for p in sorted(primaries[x],key=lambda s:int(s["primary_index"])):
                box,du,v=spatial_cell(
                    geom,int(crossing["station"]),float(crossing["u_in_mm"]),
                    float(p["u_rel_cm"]),float(p["v_abs_cm"]))
                positions.append((int(p["primary_index"]),int(p["cluster_index"]),
                                  box,du,v))
            refseed=ref_seeds[(x,geom)]
            if len(refseed)!=len(positions):
                raise ValueError("Full Garfield has different Heed population")
            byidx={int(r["primary_index"]):r for r in refseed}
            for i,cluster,box,du,v in positions:
                row=byidx[i]
                # The FULL reference seed might be outside the nearest local
                # cell, so only its periodic relative wire phase is invariant.
                recorded=float(row["seed_u_mm"])
                local_recorded=recorded-round_away(recorded/pitch)*pitch
                if abs(local_recorded-du)>1.e-3:
                    raise ValueError("REFERENCE Garfield and LUT sampling have different local u phase")
                if abs(float(row["seed_v_mm"])-v)>1.e-3:
                    raise ValueError("REFERENCE Garfield and LUT sampling have different v coordinate")
            seed_rng=random.Random(stable_seed("NA60-MUON-LUT500-V1",x,geom))
            gain_values=[];a_hemi_values=[];a_tb_values=[];zero_counts=[]
            for rep in range(replicates):
                counts=[0,0,0,0,0]
                zero=0
                for i,cluster,box,du,v in positions:
                    source=lut[geom][box]
                    n,up,down,top,bot=source[seed_rng.randrange(len(source))]
                    counts[0]+=n;counts[1]+=up;counts[2]+=down
                    counts[3]+=top;counts[4]+=bot
                    zero+=int(n==0)
                n,up,down,top,bot=counts
                gain_values.append(n);zero_counts.append(zero)
                a_hemi_values.append(asym(up,down))
                a_tb_values.append(asym(top,bot))
                ensemble.append({
                    "crossing_id":x,"configuration":geom,"replicate":rep,
                    "n_Heed_primaries":len(positions),"n_ion_births":n,
                    "n_no_birth_Heed_seeds":zero,
                    "ion_up":up,"ion_down":down,"ion_top":top,"ion_bottom":bot,
                    "A_ion_hemi":a_hemi_values[-1] if a_hemi_values[-1] is not None else "",
                    "A_ion_TB":a_tb_values[-1] if a_tb_values[-1] is not None else "",
                })
            r=refhits[(x,geom)]
            full_ion=int(r["ion_births"])
            full_zero=sum(int(int(a["ion_births"])==0) for a in refseed)
            full_hemi=asym(int(r["ion_up"]),int(r["ion_down"]))
            full_tb=asym(int(r["ion_top"]),int(r["ion_bottom"]))
            score=summarize_samples(gain_values)
            zero_stats=summarize_samples(zero_counts)
            a_hemi_finite=[v for v in a_hemi_values if v is not None]
            a_tb_finite=[v for v in a_tb_values if v is not None]
            h=summarize_samples(a_hemi_finite) if a_hemi_finite else None
            t=summarize_samples(a_tb_finite) if a_tb_finite else None
            report.append({
                "crossing_id":x,"configuration":geom,
                "station":crossing["station"],
                "theta_u_deg":crossing["theta_u_in_deg"],
                "theta_w_deg":crossing["theta_w_in_deg"],
                "path_mm":crossing["path_length_mm"],
                "n_Heed_clusters":r["heed_clusters"],
                "n_Heed_primaries":len(positions),
                "full_ion_births":full_ion,
                "lut_mean_ion_births":score["mean"],
                "lut_sd_ion_births":score["sd"],
                "lut_p05_ion_births":score["p05"],
                "lut_p50_ion_births":score["p50"],
                "lut_p95_ion_births":score["p95"],
                "full_within_lut_90pct":int(score["p05"]<=full_ion<=score["p95"]),
                "full_no_birth_seeds":full_zero,
                "lut_mean_no_birth_seeds":zero_stats["mean"],
                "full_A_hemi":full_hemi if full_hemi is not None else "",
                "lut_mean_A_hemi":h["mean"] if h else "",
                "full_A_TB":full_tb if full_tb is not None else "",
                "lut_mean_A_TB":t["mean"] if t else "",
            })
        print(f"Muon {len({r['crossing_id'] for r in report})}/{N_MUONS}: "
              f"{x}, TrackHeed primary electrons={len(primaries[x])}",flush=True)
    if len(report)!=N_MUONS*len(GEOMS):
        raise ValueError("Incomplete true-muon × four-geometry LUT comparison")
    bygeom={}
    for geom in GEOMS:
        rs=[r for r in report if r["configuration"]==geom]
        ref_mean=statistics.mean(r["full_ion_births"] for r in rs)
        lut_mean=statistics.mean(r["lut_mean_ion_births"] for r in rs)
        bygeom[geom]={
            "n_independent_real_Geant4_muons":len(rs),
            "n_true_Heed_primary_electrons":sum(int(r["n_Heed_primaries"]) for r in rs),
            "full_reference_mean_ion_births_per_muon":ref_mean,
            "lut_mean_ion_births_per_muon":lut_mean,
            "lut_over_full_mean_ratio":lut_mean/ref_mean if ref_mean else None,
            "fraction_full_reference_within_5_95pct_LUT_interval":
                statistics.mean(r["full_within_lut_90pct"] for r in rs),
            "mean_full_reference_Heed_seeds_with_zero_births":
                statistics.mean(r["full_no_birth_seeds"] for r in rs),
            "mean_LUT_Heed_seeds_with_zero_births":
                statistics.mean(r["lut_mean_no_birth_seeds"] for r in rs),
        }
    fieldnames=tuple(report[0])
    with (COMPARISON/"muon_comparison.csv").open("w",newline="") as f:
        wr=csv.DictWriter(f,fieldnames=fieldnames);wr.writeheader();wr.writerows(report)
    with (COMPARISON/"muon_lut_replicas.csv").open("w",newline="") as f:
        wr=csv.DictWriter(f,fieldnames=tuple(ensemble[0]))
        wr.writeheader();wr.writerows(ensemble)
    info={
        "schema":"NA60_GEANT4_HEED_TO_GARFIELD_VS_EMPIRICAL_LUT500_MUON_V1",
        "status":"complete",
        "n_true_geant4_muon_crossings":len(selected),
        "n_candidate_geometry":4,"n_lut_replicates_per_muon_geometry":replicates,
        "n_synthetic_lut_muon_replicas":len(ensemble),
        "reference_manifest_sha256":ref_hash,
        "overnight_lut_manifest_sha256":lut_hash,
        "original_lut500_overnight_avalanches":84000,
        "source_cpp_sha256":sha(CXX),"git_commit":commit(),
        "reference_dir":str(REFERENCE),
        "empirical_lut_dir":str(OVERNIGHT),
        "method":"Common real Geant4->single TrackHeed cluster+primary-electron population; full Garfield ONE independent stochastic realization per electron; LUT200 samples of per-cell empirical 20 avalanche census (with replacement), independent RNG per muon and geometry. Finite bin 0.1mm u × 0.2mm v, anode phase by station, B=0, 1800 V. Do NOT compare electron-by-electron identical Garfield random streams; conditional distributions compared across muons.",
        "limitations":"NOT a physical incident-muon angle PDF (balanced signed-station strata); not a fully interpolated continuous LUT; only 20 avalanches per LUT spatial cell; ref has ONE Garfield realization per Heed primary; LUT replicas are conditioned on fixed Heed ionization, neglect all inter-avalanche correlations; ion birth angle sectors, NOT induced cathode signals or final electronics; no validated ion-birth time LUT yet.",
        "by_geometry":bygeom,
    }
    write_json(COMPARISON/"summary.json",info)
    outputs={p.name:sha(p) for p in sorted(COMPARISON.iterdir()) if p.is_file()}
    write_json(COMPARISON/"manifest.json",{**info,"outputs_sha256":outputs})
    print("PASS 20 real Geant4×Heed muons, 4 fields, 200 independent LUT "
          "muon-response replicas each; same exact Heed primary input.",flush=True)
    for geom,v in bygeom.items():print(geom,json.dumps(v,sort_keys=True),flush=True)

def publish_all():
    selected,_,_=verify_ref()
    if not (COMPARISON/"manifest.json").is_file():
        raise FileNotFoundError("First compare20; no completed LUT muon comparison")
    comp=json.loads((COMPARISON/"manifest.json").read_text())
    if comp.get("status")!="complete":
        raise ValueError("Muon-level LUT comparison not complete")
    for name,digest in comp.get("outputs_sha256",{}).items():
        if sha(COMPARISON/name)!=digest:
            raise ValueError("Muon comparison changed since manifest: "+name)
    for source,run_id in (
      (REFERENCE,"real_muon20_geant4_heed_full_garfield_v1_20261029"),
      (COMPARISON,"real_muon20_lut500_vs_full_garfield_v1_20261029")):
        bundle=drive_prepare(source,run_id)
        publish(bundle,EXPECTED_REMOTE,upload=True)
    print("PASS both immutable REAL-muon Garfield/reference and independent "
          "empirical LUT comparison published and rclone readback verified")

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode",choices=(
        "prepare20","check20","reference20","compare20","publish20","status"))
    parser.add_argument("--replicates",type=int,default=REPLICATES)
    args=parser.parse_args()
    if args.mode=="status":
        for p in (INPUT,REFERENCE,COMPARISON):
            print(p, "exists" if p.exists() else "missing")
        return
    if args.mode=="prepare20":
        make(max_events=20,chunk=1,output=INPUT)
        input_verified()
        print("PASS twenty ORIGINAL distinct Geant4 Jpsi muon gas crossings; "
              "no Garfield run")
    elif args.mode=="check20":
        inp,_,_=input_verified()
        check_executable()
        validate_native_csv_before_physics(inp)
        print("PASS REAL native Geant4 input and deterministic TrackHeed preflight")
    elif args.mode=="reference20":
        run_reference()
    elif args.mode=="compare20":
        compare(args.replicates)
    elif args.mode=="publish20":
        publish_all()

if __name__=="__main__":main()
