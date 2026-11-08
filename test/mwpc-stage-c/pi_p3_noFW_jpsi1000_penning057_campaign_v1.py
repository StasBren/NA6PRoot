#!/usr/bin/env python3
"""Exactly same 1,000 frozen Geant4 J/psi muon parents, Penning ON (r=0.57).

This is a NEW, SHA-sealed, read-only-to-original p3_noFW campaign.
It reuses the successfully tested independent Penning C++ executable and
the rigorously validated historical 50x20 native benchmark runner.

Physics: 3-mm pitch, 30-um wires, cathodes +/-3 mm, 1800 V,
Ar/CO2 (70:30), B=0. Penning Ar excitation-to-ionisation r=0.57.
NOT the later p4 2+4-mm asymmetric actual prototype.

Modes:
selftest / prepare / preflight / run / gate / status / audit / plots

IMPORTANT: First run ALWAYS at most ONE 20-parent chunk. Subsequent
run requires a passing paired Penning-on/off gain sanity gate.
Never erase or overwrite an unsealed native physics chunk.
"""
from __future__ import annotations
import argparse
import csv
import json
import math
import shutil
from collections import defaultdict
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import pi_p3_noFW_jpsi1000_native_benchmark_v1 as base
from pi_p3_noFW_penning_gain_pilot_v1 import code_identity

REPO=Path(__file__).resolve().parents[2]
ORIG=REPO/"test_runs/mwpc_stage_c/pi_p3_noFW_jpsi1000_native_v1"
ROOT=REPO/"test_runs/mwpc_stage_c/pi_p3_noFW_jpsi1000_penning057_v1"
PLOTS=ROOT/"plots_collected_gain_and_ion_birth_angles_v1"
PILOT=REPO/"test_runs/mwpc_stage_c/pi_p3_noFW_jpsi_penning057_gain_pilot_v1/comparison.json"
PEN_CXX=REPO/"test/mwpc-avalanche/phaseB_native_muon_p3_noFW_jpsi_penning057_pilot_v1.cxx"
PEN_EXE=REPO/"test/mwpc-avalanche/build/mwpc_p3_noFW_jpsi_penning057_pilot_v1"
ORIG_CPP=REPO/"test/mwpc-avalanche/phaseB_native_muon_p3_noFW_jpsi_runtime_v1.cxx"
N=1000
GROUP=20
BINS_PHI=np.linspace(-180.,180.,73)

# Preserve original base physics functions; reconfigure only THIS new module's
# invocation of the frozen runner. Its source/output formats are unchanged.
base.ROOT=ROOT
base.INPUT=ROOT/"input"
base.CHUNKS=ROOT/"native_chunks"
base.CXX=PEN_CXX
base.EXE=PEN_EXE
_original_validate=base.validate_native

def must(c,m):
    if not c:raise RuntimeError(m)

def original_sealed_chunk(i):
    d=ORIG/"native_chunks"/f"chunk_{i:04d}"
    p=d/"complete.json"
    must(p.is_file(),"Historical Penning-OFF chunk is missing seal: "+str(d))
    seal=base.jread(p)
    must(seal["schema"]=="NA60_JPSI_20PARENT_P3NOFW_SHA_V1"
         and seal["chunk"]==i and seal["n_parents"]==GROUP,
         "Historical original source chunk is not the expected 20 parents")
    must(seal["physics"]["cpp_sha256"]==base.sha(ORIG_CPP),
         "Original Penning-OFF source has changed")
    for f in ("native_muon_crossings_input.csv",
              "muon_reference_crossing.csv",
              "muon_reference_seed.csv",
              "muon_reference_heed_primaries.csv"):
        must(base.sha(d/f)==seal["file_sha256"][f],
             "Original source or original real Garfield physics SHA changed: "+str(d/f))
    return d,seal

def pilot_gate():
    # This verifies the previously user-executed 3-muon Penning test,
    # rather than trusting a physical setting without running it.
    must(PILOT.is_file(),"The paired Penning ON/OFF 3-muon pilot is not sealed: "+str(PILOT))
    p=base.jread(PILOT)
    must(p["schema"]=="P3_NOFW_PAIRED_REAL_JPSI_PENNING057_PILOT_V1"
         and p["paired_identical_primary_positions_and_seeds"] is True
         and p["n_exact_matched_TrackHeed_primaries"]>=60
         and 2.<p["gain_ratio_penning_on_over_off"]<30.
         and p["gain_uncapped"] is True,
         "The prior 3-muon physical gain pilot is missing/failed")
    code_identity()  # Enforces no non-Penning change versus frozen C++.
    print(f"PASS known paired Penning pilot ON/OFF ratio "
          f"{p['gain_ratio_penning_on_over_off']:.3f}; code identical except transfer",flush=True)

def prepare():
    must(not ROOT.exists(),"Refusing to overwrite Penning1000 campaign "+str(ROOT))
    pilot_gate()
    oldsrc=ORIG/"input/native_muon_crossings_input.csv"
    oldmeta_path=ORIG/"input/input_manifest.json"
    meta=base.jread(oldmeta_path)
    must(meta["schema"]=="NA60_JPSI_1000_P3NOFW_GARFIELD_BENCH_V1"
         and meta["n_parents"]==N and meta["group_size"]==GROUP
         and meta["geometry"]=="p3_noFW"
         and base.sha(oldsrc)==meta["selected_input_sha256"],
         "Historical native Geant4 input is not the frozen 1000-muon source")
    # Before any expensive computation, establish all baseline sources exist.
    for i in range(N//GROUP):
        original_sealed_chunk(i)
    orig=list(base.csvread(oldsrc,base.FIELDS))
    must(len(orig)==N and len({base.parent_key(r) for r in orig})==N,
         "Original 1000 Geant4 parents are not unique")
    ROOT.mkdir()
    base.INPUT.mkdir()
    new_src=base.INPUT/"native_muon_crossings_input.csv"
    shutil.copyfile(oldsrc,new_src)
    must(base.sha(new_src)==base.sha(oldsrc),"Real selected parent input copy changed")
    newmeta=dict(meta)
    newmeta.update({
        "executable_provenance_prepared":base.check_binary(),
        "schema":"NA60_JPSI_1000_P3NOFW_GARFIELD_BENCH_V1",
        "physical_variant":"Penning ON r=0.57, Ar 70 / CO2 30, original 1800V",
        "paired_baseline_manifest_sha256":base.sha(oldmeta_path),
        "paired_baseline_source_sha256":base.sha(oldsrc),
        "paired_baseline":"pi_p3_noFW_jpsi1000_native_v1",
        "penning_code_diff_audited":True
    })
    base.jwrite(base.INPUT/"input_manifest.json",newmeta)
    base.check_input()
    print("PASS 1000 IDENTICAL original Geant4 J/psi parents prepared",flush=True)
    print("PASS 50 old Penning-OFF blocks SHA-verified; no Garfield run",flush=True)
    print("NEW separate output",ROOT,flush=True)

def kseed(x):
    return x["crossing_id"],int(x["primary_index"])

def paired_validate(d,batch):
    """Called BEFORE SHA sealing EVERY 20-parent Penning ON block."""
    stats=_original_validate(d,batch) # exact complete Garfield per-seed and births
    i=int(d.name.split("_")[-1])
    old,oldseal=original_sealed_chunk(i)
    must(base.sha(old/"native_muon_crossings_input.csv")==
         base.sha(d/"native_muon_crossings_input.csv"),
         "ON and OFF used different Geant4 real-muon crossings")
    log=(d/"physics.log").read_text(errors="replace")
    must("PHYSICS_PENNING_AR_CO2_70_30_R=0.57" in log,
         "Production binary did not report Penning enabled")
    must(stats["n_capped_seeds"]==0,
         "Penning ON produced a SIZE-CAPPED avalanche; preserve raw chunk, "
         "DO NOT count as real gain or continue.")

    _,oldseed=readcsv(old/"muon_reference_seed.csv")
    _,newseed=readcsv(d/"muon_reference_seed.csv")
    must(len(oldseed)==len(newseed)==stats["n_full_garfield_avalanches"],
         "Not an identical count of true Heed initial electrons")
    for x,y in zip(oldseed,newseed):
        must(kseed(x)==kseed(y),"Original TrackHeed primary ordering differs")
        must(int(x["garfield_seed"])==int(y["garfield_seed"]),
             "Garfield random seed changed across paired physics")
        for field in ("seed_u_mm","seed_v_mm","seed_w_mm","seed_energy_eV"):
            must(abs(float(x[field])-float(y[field]))<1.e-9,
                 "Real initial electron state changed in field "+field)

    _,oldp=readcsv(old/"muon_reference_heed_primaries.csv")
    _,newp=readcsv(d/"muon_reference_heed_primaries.csv")
    must(len(oldp)==len(newp)==stats["n_trackheed_electrons"],
         "TrackHeed electron population differs with Penning ON")
    for x,y in zip(oldp,newp):
        must(kseed(x)==kseed(y) and x["cluster_index"]==y["cluster_index"],
             "Different Heed electron identity or cluster")
        for field in ("u_rel_cm","v_abs_cm","w_rel_cm","t_ns",
                      "raw_heed_energy_eV","raw_dx","raw_dy","raw_dz"):
            must(abs(float(x[field])-float(y[field]))<1.e-9,
                 "TrackHeed real electron differs in "+field)
    on=np.array([int(x["collected_e"]) for x in newseed],dtype=np.float64)
    off=np.array([int(x["collected_e"]) for x in oldseed],dtype=np.float64)
    must(len(on)>0 and float(np.mean(off))>0,"Invalid paired gain statistics")
    stats.update({
        "Penning_on_sum_collected_e":int(np.sum(on)),
        "Penning_off_paired_sum_collected_e":int(np.sum(off)),
        "Penning_on_mean_collected_e":float(np.mean(on)),
        "Penning_off_paired_mean_collected_e":float(np.mean(off)),
        "Penning_on_off_paired_gain_ratio":float(np.mean(on)/np.mean(off)),
        "Penning_on_zero_collected":int(np.sum(on==0)),
        "same_Heed_electron_positions_and_Garfield_seeds":True,
        "original_Penning_off_seal_sha256":base.sha(old/"complete.json"),
    })
    print(f"PASS paired real TrackHeed starts chunk {i:02d}: "
          f"ON gain {np.mean(on):.2f}, OFF gain {np.mean(off):.2f}, "
          f"ON/OFF {np.mean(on)/np.mean(off):.3f}, "
          f"zero ON {int(np.sum(on==0))}/{len(on)}",flush=True)
    return stats

def readcsv(path):
    with Path(path).open(newline="") as f:
        r=csv.DictReader(f)
        must(r.fieldnames is not None,"Empty real Garfield CSV: "+str(path))
        return r.fieldnames,list(r)

# The base writer calls validate_native before creating complete.json.
base.validate_native=paired_validate

def gate():
    rr,_=base.check_input()
    phys=base.checked_physics()
    d=base.folder(0)
    must((d/"complete.json").is_file(),
         "First run exactly ONE 20-muon chunk, then use gate")
    receipt=base.audit_chunk(0,rr[:GROUP],phys)
    s=receipt["stats"]
    must(s["same_Heed_electron_positions_and_Garfield_seeds"] is True
         and s["n_capped_seeds"]==0
         and s["n_full_garfield_avalanches"]>=200
         and 2.<s["Penning_on_off_paired_gain_ratio"]<30.,
         "First 20 matched Penning muons FAILED gain/cap sanity gate. "
         "Do NOT launch the remaining 980 muons.")
    print("PASS PHYSICS GATE: 20 exact paired parents, "
          f"Penning gain {s['Penning_on_mean_collected_e']:.2f}, "
          f"OFF {s['Penning_off_paired_mean_collected_e']:.2f}, "
          f"ratio {s['Penning_on_off_paired_gain_ratio']:.3f}, "
          "uncapped, SHA complete",flush=True)

def run(max_chunks,hours):
    pilot_gate()
    rr,_=base.check_input()
    phys=base.checked_physics()
    first=base.folder(0)/"complete.json"
    if first.exists():
        gate()
    else:
        must(max_chunks==1 and not base.folder(0).exists(),
             "First physics run MUST be --max-chunks 1 (20 real muons) "
             "with no existing unsealed chunk")
    base.run(max_chunks,hours)

def status(full=False):
    pilot_gate()
    print("PHYSICS VARIANT: P3 NO FIELD WIRES, PENNING ON, r=0.57")
    base.status(audit=full)
    rr,_=base.check_input()
    phy=base.checked_physics()
    totals=defaultdict(float)
    chunks_done=0
    for i,batch in enumerate(base.chunks(rr)):
        d=base.folder(i)
        if not (d/"complete.json").is_file():break
        seal=base.audit_chunk(i,batch,phy,verify_files=full)
        ss=seal["stats"]
        for k in ("Penning_on_sum_collected_e",
                  "Penning_off_paired_sum_collected_e",
                  "n_full_garfield_avalanches","Penning_on_zero_collected"):
            totals[k]+=ss[k]
        chunks_done+=1
    if chunks_done:
        n=totals["n_full_garfield_avalanches"]
        on=totals["Penning_on_sum_collected_e"]/n
        off=totals["Penning_off_paired_sum_collected_e"]/n
        print(f"PAIRED PINNED 1000-muon subset: {GROUP*chunks_done}/{N} complete")
        print(f"Mean collected gain Penning ON {on:.3f}, OFF {off:.3f}; "
              f"ON/OFF {on/off:.3f}, zero gain {totals['Penning_on_zero_collected']/n:.3%}")
        print("No Garfield source or archive was modified.")

def cpp_round(x):
    return math.floor(x+.5) if x>=0 else math.ceil(x-.5)

def angle_bin(u,v):
    wire=cpp_round(u/3.)*3.
    phi=math.degrees(math.atan2(v,u-wire))
    return max(0,min(71,int((phi+180.)*72/360.)))

def gain_stats(x):
    a=np.asarray(x,dtype=np.float64)
    active=a[a>0]
    must(len(active)>0,"No real positive collected gain")
    return dict(n=len(a),n_positive=len(active),
                zero=int(np.sum(a==0)),zero_frac=float(np.mean(a==0)),
                mean_all=float(np.mean(a)),rms_all=float(np.std(a,ddof=0)),
                mean_positive=float(np.mean(active)),
                rms_positive=float(np.std(active,ddof=0)),
                max=int(np.max(a)))

def gain_figures(collected,off):
    g=gain_stats(collected)
    old=gain_stats(off)
    width=50
    top=int(max(collected))
    e=np.r_[0.,np.arange(1,top+width+1,width,dtype=np.float64)]
    count=np.histogram(collected,bins=e)[0]
    must(sum(count)==len(collected) and count[0]==g["zero"],
         "N=0 and full collected gain histogram do not conserve counts")
    for exclude,filename in [
        (False,"01_penning_on_avalanche_gain_with_zero.png"),
        (True,"02_penning_on_avalanche_gain_positive_only.png")
    ]:
        hist=count[1:] if exclude else count
        bins=e[1:] if exclude else e
        denom=g["n_positive"] if exclude else g["n"]
        prob=hist/denom
        avg=g["mean_positive"] if exclude else g["mean_all"]
        sd=g["rms_positive"] if exclude else g["rms_all"]
        must(abs(float(np.sum(prob))-1.)<1e-10,"Histogram wrong normalization")
        legend=(f"Penning ON collected e; mean={avg:,.1f}; "
                f"RMS={sd:,.1f}; N={denom:,}")
        fig,axes=plt.subplots(1,2,figsize=(13.5,5.1),sharex=True)
        for ax in axes:
            val=np.r_[prob,prob[-1]]
            if ax is axes[1]:val=np.where(val>0,val,np.nan)
            ax.step(bins,val,where="post",lw=1.6,color="#e46b08",label=legend)
            if not exclude:
                ax.scatter([.5],[g["zero_frac"]],edgecolor="#e46b08",
                           facecolor="white",s=45,zorder=5)
            ax.grid(alpha=.15)
            ax.set_xlim(bins[0],bins[-1])
            ax.set_xlabel("Anode-collected electrons per Heed primary electron")
            ax.set_ylabel("Probability per bin")
            ax.legend(fontsize=8)
        axes[0].set_ylim(bottom=0);axes[0].set_title("Linear probability")
        axes[1].set_yscale("log");axes[1].set_title("Logarithmic probability")
        note=(f"N=0 bin [0,1); positive bins {width} electrons; "
              f"zero fraction {g['zero_frac']:.2%}" if not exclude else
              f"N=0 excluded; positive bins {width} electrons; conditional N>0")
        fig.suptitle("p3_noFW Penning ON (r=0.57), 1,000 original real J/psi muons\n"
                     +note)
        fig.tight_layout()
        fig.savefig(PLOTS/filename,dpi=185);plt.close(fig)
    return g,old

def angle_figures():
    ion_weighted=np.zeros(72)
    avalanche_weighted=np.zeros(72)
    nonzero=0
    npoints=0
    full_ion_births=0
    for i in range(N//GROUP):
        d=base.folder(i)
        seeds={}
        for r in base.csvread(d/"muon_reference_seed.csv",
                             ("crossing_id","primary_index","ion_births")):
            b=int(r["ion_births"])
            seeds[(r["crossing_id"],int(r["primary_index"]))]=b
            if b>0:
                nonzero+=1
                full_ion_births+=b
        points_per=defaultdict(int)
        for r in base.csvread(d/"muon_reference_ion_birth_sample.csv",
                             ("crossing_id","primary_index","u_mm","v_mm",
                              "seed_ion_births_total","reservoir_weight")):
            k=(r["crossing_id"],int(r["primary_index"]))
            must(k in seeds and seeds[k]>0,"Missing original one-electron Garfield seed")
            b=seeds[k]
            m=min(128,b)
            must(int(r["seed_ion_births_total"])==b,"Ion birth sample count mismatch")
            w=float(r["reservoir_weight"])
            must(abs(w-b/m)<1e-7*max(1.,b/m),"Ion birth reservoir biased")
            idx=angle_bin(float(r["u_mm"]),float(r["v_mm"]))
            avalanche_weighted[idx]+=1./m
            ion_weighted[idx]+=b/m
            points_per[k]+=1
            npoints+=1
        must(all(points_per[k]==min(128,b) for k,b in seeds.items() if b>0),
             "Not all 128 unbiased sampled ion births were saved")
        if (i+1)%10==0:
            print(f"PASS original Penning Garfield ion-angle reservoirs {i+1}/50",flush=True)
    must(abs(sum(avalanche_weighted)-nonzero)<1e-5*nonzero,
         "Equal-avalanche angular probability fails conservation")
    must(abs(sum(ion_weighted)-full_ion_births)<1e-5*full_ion_births,
         "Ion-weighted angle probability fails full-birth conservation")
    c=(BINS_PHI[:-1]+BINS_PHI[1:])/2.
    for p,label,y in [
        ("03_penning_on_ion_birth_angular_equal_avalanche.png",
         "Each ion-producing avalanche equally weighted",
         avalanche_weighted/nonzero),
        ("04_penning_on_ion_birth_angular_ion_weighted.png",
         "Each ion birth equally weighted (128-point unbiased sampling)",
         ion_weighted/full_ion_births)
    ]:
        fig,ax=plt.subplots(figsize=(11.5,5.6))
        ax.step(c,y,where="mid",lw=1.85,color="#e46b08",label=label)
        ax.set_xlim(-180,180)
        ax.set_xlabel("Birth azimuth phi around nearest 3-mm anode [degrees]")
        ax.set_ylabel("Probability per 5-degree angular bin")
        ax.set_title("p3_noFW, Penning ON r=0.57: true Garfield ion-birth directions\n"
                     "128 uniformly sampled birth positions per avalanche, not full cloud")
        ax.legend(fontsize=9);ax.grid(alpha=.17)
        fig.tight_layout();fig.savefig(PLOTS/p,dpi=185);plt.close(fig)
    with (PLOTS/"angular_distribution.csv").open("x",newline="") as f:
        w=csv.writer(f)
        w.writerow(("phi_bin_center_deg","equal_avalanche_probability",
                    "ion_weighted_probability"))
        for j,x in enumerate(c):
            w.writerow((x,avalanche_weighted[j]/nonzero,
                        ion_weighted[j]/full_ion_births))
    return {"n_ion_producing_avalanches":nonzero,
            "full_true_ion_births":full_ion_births,
            "n_reservoir_birth_coordinates":npoints,
            "angular_bin_width_deg":5,
            "weighting":"Nonzero ion-producing avalanche equal or exact-ion-count weighted",
            "statistics_note":"Only 128 sampled positions per avalanche, not exact individual ion-birth positions"}

def plots():
    rr,_=base.check_input()
    phy=base.checked_physics()
    must(not PLOTS.exists(),"Refusing to overwrite previous validated plots: "+str(PLOTS))
    # Prove every 20-parent source+final physics file is SHA sealed before use.
    for i,batch in enumerate(base.chunks(rr)):
        base.audit_chunk(i,batch,phy,verify_files=True)
    must((base.folder(49)/"complete.json").is_file(),
         "1000 muons incomplete; cannot produce final Penning figures")
    all_on=[];all_off=[]
    for i in range(N//GROUP):
        new=base.folder(i)
        old,_=original_sealed_chunk(i)
        seed_on=list(base.csvread(new/"muon_reference_seed.csv",
                                 ("crossing_id","primary_index","collected_e")))
        seed_off=list(base.csvread(old/"muon_reference_seed.csv",
                                  ("crossing_id","primary_index","collected_e")))
        must(len(seed_on)==len(seed_off) and
             all(kseed(a)==kseed(b) for a,b in zip(seed_on,seed_off)),
             "True initial Garfield seeds not matched for gain comparison")
        all_on.extend(int(x["collected_e"]) for x in seed_on)
        all_off.extend(int(x["collected_e"]) for x in seed_off)
    must(len(all_on)==len(all_off)>20000,"Missing true microscopic Garfield avalanches")
    PLOTS.mkdir(parents=True)
    g,on_off=gain_figures(all_on,all_off)
    ang=angle_figures()
    report={
        "schema":"NA60_P3_NOFW_JPSI1000_PENNING_ON_GAIN_ANGLE_V1",
        "n_true_real_Geant4_Jpsi_muon_parents":N,
        "n_full_Garfield_avalanches":len(all_on),
        "Penning_on_gain":g,
        "same_primaries_Penning_off_gain":on_off,
        "Penning_on_to_off_mean_gain_ratio":g["mean_all"]/on_off["mean_all"],
        "angles":ang,
        "physics_variant":"Ar/CO2 70/30, Penning r=0.57, p=3mm, cathodes +/-3mm, +1800V",
        "caveat":"NOT the real p4 2+4-mm prototype",
        "baseline_input_sha256":base.sha(ORIG/"input/native_muon_crossings_input.csv"),
        "penning_input_sha256":base.sha(ROOT/"input/native_muon_crossings_input.csv"),
        "penning_executable_sha256":base.sha(PEN_EXE),
        "script_sha256":base.sha(Path(__file__)),
        "figures_sha256":{p.name:base.sha(p) for p in PLOTS.glob("*.png")},
    }
    base.jwrite(PLOTS/"summary.json",report)
    print("PASS 1000 PAIRED true Geant4 J/psi parent Penning ON/OFF physics",flush=True)
    print(f"GAIN ON mean(all) {g['mean_all']:.3f}, "
          f"mean(>0) {g['mean_positive']:.3f}, zero {g['zero_frac']:.3%}")
    print(f"GAIN OFF mean(all) {on_off['mean_all']:.3f}; "
          f"paired ON/OFF ratio {report['Penning_on_to_off_mean_gain_ratio']:.3f}")
    print("Angular full ion births:",ang["full_true_ion_births"])
    print("OUTPUT",PLOTS,flush=True)

def selftest():
    base.selftest()
    assert angle_bin(0.,1.)==54 and angle_bin(0.,-1.)==18
    assert angle_bin(3.,1.)==54 and angle_bin(-3.,1.)==54
    assert base.GROUP==GROUP and base.N==N
    print("PASS 1000 SAME original parent and 3mm nearest-anode angle selftest")

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("mode",choices=("selftest","prepare","preflight",
                                   "run","gate","status","audit","plots"))
    p.add_argument("--max-chunks",type=int,default=1)
    p.add_argument("--hours",type=float,default=24.)
    a=p.parse_args()
    must(1<=a.max_chunks<=50 and 0<a.hours<=72,"Invalid safe chunk/runtime budget")
    if a.mode=="selftest":selftest()
    elif a.mode=="prepare":prepare()
    elif a.mode=="preflight":
        pilot_gate()
        base.preflight()
    elif a.mode=="run":run(a.max_chunks,a.hours)
    elif a.mode=="gate":gate()
    elif a.mode=="status":status()
    elif a.mode=="audit":status(full=True)
    else:plots()

if __name__=="__main__":
    main()
