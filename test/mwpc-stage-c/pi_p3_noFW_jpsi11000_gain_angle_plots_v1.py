#!/usr/bin/env python3
"""3 requested physics plots for SHA-frozen 11,000 real J/psi p3_noFW muons.

1. Real electron number per microscopic Garfield avalanche: GetAvalancheSize
   ne AND genuine anode-collected endpoints; zero events kept separately.
2. Real birth-angle azimuth around nearest 3 mm anode: <=128 uniform ion
   births per avalanche, weighted by TRUE ion count; separately show
   each nonzero avalanche with equal weight and circular-mean directions.
3. Mean effective anode gain = sum collected electron endpoints / N REAL
   Heed primary electrons. Conditional mean gain versus (u_start,v_start).

Never call Garfield, never modify SHA-frozen physics, never resimulate.
Handles SHA-sealed high-Heed-count recovery chunk_0426 via subdirectory.
Each observed Geant4 event has exactly one original parent; 11k disjoint.
Angular distributions are reservoir estimates, NOT full 333-million-birth
coordinates and NOT final ion motion/induced cathode currents.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from collections import defaultdict

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO=Path(__file__).resolve().parents[2]
DATA=REPO/"test_runs/mwpc_stage_c"
OLD=DATA/"pi_p3_noFW_jpsi1000_native_v1"
NEW=DATA/"pi_p3_noFW_jpsi_extra10000_native_v1"
OUT=DATA/"pi_p3_noFW_jpsi11000_gain_angle_linear_v2"
ELECTRONS_ONLY_OUT=DATA/"pi_p3_noFW_jpsi11000_electrons_linear_v2"
GNAME="p3_noFW"
NR=72
PHI_EDGES=np.linspace(-180.,180.,NR+1)
PHI_CENTRE=(PHI_EDGES[:-1]+PHI_EDGES[1:])/2
WIRE_PITCH_MM=3.
K_AVAL=("crossing_id","primary_index")
K_ION=("crossing_id","configuration","primary_index","u_mm","v_mm",
       "w_mm","t_ns","seed_ion_births_total","reservoir_weight")
K_SEED=("crossing_id","configuration","primary_index","garfield_ne",
        "garfield_ni","ion_births","collected_e","attached","other",
        "successful","capped","seed_u_mm","seed_v_mm")

def must(flag,text):
    if not flag:raise RuntimeError(text)

def sha(file):
    h=hashlib.sha256()
    with file.open("rb") as f:
        for part in iter(lambda:f.read(1<<20),b""):
            h.update(part)
    return h.hexdigest()

def csv_rows(path,required):
    with path.open(newline="") as f:
        rdr=csv.DictReader(f)
        must(set(required).issubset(rdr.fieldnames or []),
             "Wrong physics CSV schema "+str(path))
        yield from rdr

def parse_int(x):
    z=int(x)
    must(z>=0,"Negative microscopic Garfield count")
    return z

def key(r):
    return r["crossing_id"],int(r["primary_index"])

def circular_phi(u_mm,v_mm):
    # Match C++ std::round (ties away from zero). Standard round()
    # in Python uses ties-to-even; physics uses nearest wire.
    n=math.floor(u_mm/WIRE_PITCH_MM+.5) if u_mm>=0 else math.ceil(u_mm/WIRE_PITCH_MM-.5)
    return math.degrees(math.atan2(v_mm,u_mm-n*WIRE_PITCH_MM))

def physics_files():
    selected_parent_ids=set()
    manifest=[]
    for root,nparents in ((OLD,1000),(NEW,10000)):
        inp=root/"input/native_muon_crossings_input.csv"
        m=root/"input/input_manifest.json"
        must(inp.is_file() and m.is_file(),"Missing already-frozen p3-noFW campaign "+str(root))
        meta=json.loads(m.read_text())
        must(meta["n_parents"]==nparents and meta["geometry"]=="p3_noFW"
             and meta["group_size"]==20 and meta["selected_input_sha256"]==sha(inp),
             "Input-selection provenance failed "+str(root))
        # Strict Geant4 PARENT identity, not the individual TrackHeed electron.
        with inp.open(newline="") as f:
            r=csv.DictReader(f)
            n=0
            for x in r:
                parent=(x["source"],x["chunk"],x["event"])
                must(parent not in selected_parent_ids,
                     "Original vs extra J/psi source parent overlap")
                selected_parent_ids.add(parent)
                n+=1
        must(n==nparents,"Frozen Geant4 input event count changed")
        count=0
        for chunk in range(nparents//20):
            d=root/"native_chunks"/f"chunk_{chunk:04d}"
            rec_path=d/"complete.json"
            must(rec_path.is_file(),"Unsealed original full-Garfield chunk "+str(d))
            rec=json.loads(rec_path.read_text())
            must(rec["chunk"]==chunk and rec["n_parents"]==20,
                 "Wrong completed microscopic physics block")
            layout=rec.get("recovery_layout")
            must(layout in (None,"high_heed_replay_v1"),
                 "Unknown recovery layout in original Garfield SHA checkpoint")
            p=d/layout if layout else d
            # The long ion-birth sample file MUST be checked as well.
            for name in ("muon_reference_seed.csv",
                         "muon_reference_crossing.csv",
                         "muon_reference_ion_birth_sample.csv"):
                must(name in rec["file_sha256"] and
                     sha(p/name)==rec["file_sha256"][name],
                     "Physics file SHA differs from original sealed Garfield "+str(p/name))
            manifest.append({"source":root.name,"chunk":chunk,
                             "complete_sha256":sha(rec_path),
                             "storage":"high_heed_replay_v1" if layout else "normal"})
            count+=1
            yield root,chunk,p
        must(count==nparents//20,"Incorrect native checkpoint count")
    must(len(selected_parent_ids)==11000,
         "Not 11000 unique original J/psi parent events")


def plot_electron_counts_linear(raw_ne, coll, filename, num_muons=11000):
    """Full-range, linear X and Y; ROOT-style RMS is std around mean.

    Zeros are INCLUDED in both the linear-bin probability histogram and
    empirical CDF. Means and RMSs use ALL avalanche calls, never a clipped
    or positive-only subset. Values in x data are actual integer electrons.
    """
    raw_ne=np.asarray(raw_ne,dtype=np.float64)
    coll=np.asarray(coll,dtype=np.float64)
    must(len(raw_ne)==len(coll)>0,"Electron populations must match")
    must(np.all(np.isfinite(raw_ne)) and np.all(np.isfinite(coll))
         and np.min(raw_ne)>=0 and np.min(coll)>=0,
         "Bad original Garfield electron multiplicity")
    N=len(raw_ne)
    xmax=max(float(np.max(raw_ne)),float(np.max(coll)),1.)
    # Fixed linear bins across the ENTIRE physical range, including tails;
    # no percentile truncation and no hidden logarithmic normalization.
    edges=np.linspace(0.,xmax+1.,121)
    fig,axs=plt.subplots(1,2,figsize=(13.5,5.1))
    info=[
        (raw_ne,"Garfield ne","-"),
        (coll,"Anode-collected e","--"),
    ]
    statistics={}
    for x,label,linestyle in info:
        mu=float(np.mean(x))
        rms=float(np.std(x,ddof=0))
        zeros=int(np.sum(x==0))
        name="garfield_ne" if label=="Garfield ne" else "collected_e"
        statistics[name]={"mean":mu,"RMS":rms,"zeros":zeros,"n_avalanche_calls":N}
        legend=f"{label}: mean={mu:,.2f}, RMS={rms:,.2f} (zero={zeros:,})"
        # Probability per linear bin sums to 1, including zero-G avalanches.
        axs[0].hist(x,bins=edges,histtype="step",
                    weights=np.full(N,1./N),linewidth=1.7,
                    linestyle=linestyle,label=legend)
        values,counts=np.unique(x,return_counts=True)
        cdf=np.cumsum(counts)/N
        axs[1].step(values,cdf,where="post",linewidth=1.7,
                    linestyle=linestyle,label=f"{label}: mean={mu:,.2f}, RMS={rms:,.2f}")
    axs[0].set_xlim(0,xmax+1)
    axs[0].set_ylim(bottom=0)
    axs[0].set_xlabel("Electrons in a single Garfield avalanche")
    axs[0].set_ylabel("Probability per linear-width bin (zero included)")
    axs[0].set_title("Electron multiplicity — full linear range")
    axs[1].set_xlim(0,xmax+1)
    axs[1].set_ylim(0,1.02)
    axs[1].set_xlabel("Electrons in a single Garfield avalanche")
    axs[1].set_ylabel("Empirical cumulative probability")
    axs[1].set_title("CDF — full linear range")
    for ax in axs:
        ax.legend(fontsize=8)
        ax.grid(alpha=.13)
    fig.suptitle(
        f"Real Garfield p3_noFW: {N:,} avalanches from {num_muons:,} J/psi muons\n"
        "Both axes linear; mean and RMS include zero avalanches")
    fig.tight_layout()
    fig.savefig(filename,dpi=190)
    plt.close(fig)
    return statistics


def plot_electrons_only():
    """Fast: verify source SHA and read true seed counters, not 60M ion points."""
    must(not ELECTRONS_ONLY_OUT.exists(),
         "Refusing to overwrite already generated electron plot: "+str(ELECTRONS_ONLY_OUT))
    raw_ne=[];collected=[];n_muons=0
    for root,ch,p in physics_files():
        events=list(csv_rows(p/"muon_reference_crossing.csv",
                             ("configuration","heed_primaries")))
        must(len(events)==20 and all(e["configuration"]==GNAME for e in events),
             "Missing original 20 Geant4 muons / wrong geometry")
        n_muons+=len(events)
        n_seed=0
        for r in csv_rows(p/"muon_reference_seed.csv",K_SEED):
            must(r["configuration"]==GNAME and parse_int(r["capped"])==0,
                 "Unexpected geometry or capped Garfield avalanche")
            raw_ne.append(parse_int(r["garfield_ne"]))
            collected.append(parse_int(r["collected_e"]))
            n_seed+=1
        must(n_seed==sum(int(e["heed_primaries"]) for e in events),
             "Original real Heed seed count vs Garfield calls mismatch")
        if (ch+1)%50==0:
            print(f"PASS SHA seed census {root.name} chunk {ch+1}",flush=True)
    must(n_muons==11000 and len(raw_ne)==588560,
         "11k full Garfield source/electron sample changed unexpectedly")
    ELECTRONS_ONLY_OUT.mkdir(parents=True,exist_ok=False)
    png=ELECTRONS_ONLY_OUT/"01_electrons_per_avalanche_linear.png"
    stats=plot_electron_counts_linear(raw_ne,collected,png,n_muons)
    record={
        "schema":"P3_NOFW_JPSI11000_LINEAR_ELECTRON_MULTIPLICITY_V2",
        "physics":"Original true full Garfield + Heed; no new simulations",
        "n_true_Jpsi_muon_parents":n_muons,
        "n_full_Garfield_avalanche_calls":len(raw_ne),
        "RMS_definition":"np.std(ddof=0), ROOT-style population standard deviation",
        "range":"All nonnegative values (no clipping), zeros included",
        "statistics":stats,
        "plot_sha256":sha(png),
        "script_sha256":sha(Path(__file__))
    }
    with (ELECTRONS_ONLY_OUT/"statistics.json").open("x") as f:
        json.dump(record,f,indent=2)
        f.write("\n")
    print("PASS LINEAR electron avalanche figure, no microscopic Garfield calls")
    for name,z in stats.items():
        print(f"{name}: mean={z['mean']:.4f} RMS={z['RMS']:.4f} "
              f"zero={z['zeros']} / {len(raw_ne)}")
    print("OUTPUT",ELECTRONS_ONLY_OUT,flush=True)


def analyze():
    must(not OUT.exists(),"Refusing to overwrite previously generated plots "+str(OUT))
    collected=[];raw_ne=[];ion_births=[];raw_ni=[];us=[];vs=[]
    # Each original avalanche = 1 true TrackHeed conduction electron
    # and is an independent Garfield call (but NOT an independent parent).
    angles_ion=np.zeros(NR,float)
    angles_equal_avalanche=np.zeros(NR,float)
    centroid=np.zeros(NR,float)
    sample_total=0
    actual_births=0
    missing_direction=0
    saved_cloud=0
    chk=[]
    zero_collected=0
    unsuccessful=0
    caps=0
    passed_muons=0
    for root,ch,p in physics_files():
        expected={}
        totals=[]
        for r in csv_rows(p/"muon_reference_crossing.csv",
                          ("ion_births","heed_primaries","configuration")):
            must(r["configuration"]==GNAME,"Incorrect geometry in native muon summary")
            passed_muons+=1
        for r in csv_rows(p/"muon_reference_seed.csv",K_SEED):
            must(r["configuration"]==GNAME,"Incorrect true Garfield seed geometry")
            k=key(r)
            must(k not in expected,"Duplicated individual Garfield primary")
            e=parse_int(r["collected_e"])
            n=parse_int(r["garfield_ne"])
            b=parse_int(r["ion_births"])
            ni=parse_int(r["garfield_ni"])
            must(int(r["capped"])==0,
                 "Size-capped avalanche: must not treat it as complete physics")
            expected[k]=b
            collected.append(e);raw_ne.append(n);ion_births.append(b);raw_ni.append(ni)
            us.append(float(r["seed_u_mm"]));vs.append(float(r["seed_v_mm"]))
            caps+=int(r["capped"])
            unsuccessful+=1-int(r["successful"])
            zero_collected+=int(e==0)
            actual_births+=b
        # Sample coordinates were drawn by unbiased independent RNG from
        # every individual microscopic avalanche. Their interpretation:
        # 1 ion = one vote in "ion weighted" distribution, so EACH saved
        # point has weight N_true_births / min(N_true_births,128).
        # 1 avalanche = one vote in equal-avalanche distribution, so each
        # point has weight 1 / n_saved for positive-ion avalanches.
        counts=defaultdict(int)
        sumcos=defaultdict(float)
        sumsin=defaultdict(float)
        allbins=defaultdict(int)
        for r in csv_rows(p/"muon_reference_ion_birth_sample.csv",K_ION):
            must(r["configuration"]==GNAME,"Wrong geometry of saved Garfield ion birth")
            k=key(r)
            must(k in expected and expected[k]>0,
                 "Sampled ion birth has no corresponding real full Garfield avalanche")
            nb=expected[k]
            must(int(r["seed_ion_births_total"])==nb,"Sample changed true Garfield ion count")
            m=min(nb,128)
            weight=float(r["reservoir_weight"])
            must(math.isfinite(weight) and abs(weight-nb/m)<max(1,nb/m)*1e-7,
                 "Unbiased individual-avalanche ion weight corrupt")
            phi=circular_phi(float(r["u_mm"]),float(r["v_mm"]))
            b=max(0,min(NR-1,int((phi+180.)*NR/360.)))
            angles_ion[b]+=nb/m
            angles_equal_avalanche[b]+=1./m
            counts[k]+=1
            a=math.radians(phi)
            sumcos[k]+=math.cos(a)
            sumsin[k]+=math.sin(a)
            sample_total+=1
        for k,nb in expected.items():
            m=min(nb,128)
            must(counts[k]==m,"Incomplete per-avalanche Garfield birth reservoir")
            if nb==0:continue
            saved_cloud+=1
            resultant=math.hypot(sumcos[k],sumsin[k])/m
            if resultant<0.2:
                missing_direction+=1
                continue
            phi=math.degrees(math.atan2(sumsin[k],sumcos[k]))
            binidx=max(0,min(NR-1,int((phi+180.)*NR/360.)))
            centroid[binidx]+=1
        print(f"PASS SHA + birth census {root.name} chunk {ch:04d} "
              f"avalanche_seeds={len(expected)} saved_ion_points="
              f"{sum(counts.values())}",flush=True)
    must(passed_muons==11000,"Unexpected whole-muon count from exact Garfield")
    N=len(collected)
    must(N==588560,"Original 1000 + additional 10000 Garfield avalanche number changed")
    raw_ne=np.array(raw_ne,dtype=float)
    coll=np.array(collected,dtype=float)
    ion=np.array(ion_births,dtype=float)
    ni=np.array(raw_ni,dtype=float)
    u=np.array(us,dtype=float);v=np.array(vs,dtype=float)
    must(int(ion.sum())==actual_births,"Missing full Garfield positive-ion count")
    must(abs(angles_ion.sum()-actual_births)<1e-6*max(1,actual_births),
         "Normalized weighted full-ion azimuth does not sum to true full ion count")
    must(abs(angles_equal_avalanche.sum()-saved_cloud)<1e-5*max(1,saved_cloud),
         "Equal-avalanche angular density not normalized")
    must(np.all(np.isfinite(u)) and np.all(np.isfinite(v)),
         "Nonfinite original TrackHeed starts")
    OUT.mkdir(parents=True,exist_ok=False)

    # Figure 1: true full-range linear electron multiplicity, mean and RMS.
    plot_electron_counts_linear(raw_ne,coll,OUT/"01_electrons_per_avalanche.png",N and passed_muons)

    # Figure 2: correct per-ION vs per-avalanche angular comparisons.
    # Both angles are of ion creation around nearest physical anode.
    fig,axs=plt.subplots(1,2,figsize=(12,4.8))
    axs[0].step(PHI_CENTRE,angles_ion/angles_ion.sum(),where="mid",
                label="Each ion birth equally weighted")
    axs[0].step(PHI_CENTRE,angles_equal_avalanche/saved_cloud,where="mid",
                label="Each nonzero avalanche equally weighted",ls="--")
    axs[0].set_xlabel("Birth azimuth phi around nearest anode [degree]")
    axs[0].set_ylabel("Probability / 5-degree bin")
    axs[0].legend(fontsize=8);axs[0].grid(alpha=.13)
    axs[1].step(PHI_CENTRE,centroid/max(1,centroid.sum()),where="mid")
    axs[1].set_xlabel("Circular mean direction of ion births per avalanche [degree]")
    axs[1].set_ylabel("Probability / 5-degree bin")
    axs[1].set_title(f"Coherent avalanches: {int(centroid.sum()):,} / {saved_cloud:,}"
                     "\n(vector-resultant threshold R>=0.2)")
    axs[1].grid(alpha=.13)
    fig.suptitle("True full microscopic Garfield ion births: 128-point/avalanche "
                 "unbiased angular reservoir")
    fig.tight_layout();fig.savefig(OUT/"02_ion_birth_angles_around_anode.png",dpi=190);plt.close(fig)

    # Figure 3: REAL effective collected-electron gain (seed electron=1).
    # Conditional mean gain map at exact primary Heed coordinates.
    ulocal=u-3*np.where(u>=0,np.floor(u/3+.5),np.ceil(u/3-.5))
    Uedges=np.linspace(-1.5,1.5,31);Vedges=np.linspace(-3,3,31)
    count,xE,yE=np.histogram2d(ulocal,v,bins=(Uedges,Vedges))
    sumG,_,_=np.histogram2d(ulocal,v,bins=(Uedges,Vedges),weights=coll)
    with np.errstate(invalid="ignore",divide="ignore"):
        meanG=sumG/count
    meanG[count<40]=np.nan
    mean_overall=float(coll.mean())
    mean_ion=float(ion.mean())
    balance=float(coll.sum()/max(1.,ion.sum()))
    fig,axs=plt.subplots(1,2,figsize=(11.5,5))
    # Linear color normalization: NO logarithmic color scale.
    im=axs[0].pcolormesh(Uedges,Vedges,meanG.T,shading="auto",
                          vmin=float(np.nanpercentile(meanG,2)),
                          vmax=float(np.nanpercentile(meanG,98)),
                          cmap="viridis")
    fig.colorbar(im,ax=axs[0],label="Mean collected electrons / Heed seed")
    axs[0].set_xlabel("Primary electron u relative nearest anode [mm]")
    axs[0].set_ylabel("Primary electron v [mm]")
    axs[0].set_title("Conditional mean effective gas gain\n"
                     "cells with >=40 real Garfield avalanches")
    # 1D gain along perpendicular distance v, with true electron count.
    vBin=np.digitize(v,Vedges)-1
    centres=(Vedges[:-1]+Vedges[1:])/2
    mcurve=[];ncurve=[]
    for b in range(len(centres)):
        z=coll[vBin==b];mcurve.append(z.mean() if len(z)>40 else np.nan)
        ncurve.append(len(z))
    axs[1].plot(centres,mcurve,marker=".",lw=1.4,label="Mean collected_e")
    axs[1].axhline(mean_overall,ls="--",lw=1,label=f"Global mean gain {mean_overall:.2f}")
    axs[1].set_xlabel("Initial electron v (normal to anode plane) [mm]")
    axs[1].set_ylabel("Mean anode-collected electrons per primary")
    axs[1].legend(fontsize=8);axs[1].grid(alpha=.13)
    fig.suptitle("Mean effective gas gain vs exact TrackHeed start, not ne/N_ionization")
    fig.text(.5,.016,
        f"Mean collected/seed = {mean_overall:,.2f}    "
        f"Mean Garfield ne/seed = {raw_ne.mean():,.2f}    "
        f"Mean ion pairs/seed = {mean_ion:,.2f}    "
        f"collected / ion pairs = {balance:.4g} (NOT gas gain)",
        fontsize=9,ha="center")
    fig.tight_layout(rect=[0,.05,1,.96])
    fig.savefig(OUT/"03_mean_collected_gas_gain.png",dpi=190);plt.close(fig)

    report={
        "schema":"P3_NOFW_JPSI11000_ELECTRONS_ION_AZIMUTH_GAIN_V1",
        "n_real_Geant4_parents":passed_muons,
        "n_real_microscopic_Garfield_calls":N,
        "n_total_ion_births":actual_births,
        "n_saved_real_ion_birth_coordinates":sample_total,
        "n_nonzero_ion_birth_avalanches":saved_cloud,
        "n_zero_collected":zero_collected,
        "n_unsuccessful_avalanches":unsuccessful,
        "n_capped_avalanches":caps,
        "mean_Garfield_ne_per_seed":float(raw_ne.mean()),
        "rms_Garfield_ne_per_seed":float(raw_ne.std(ddof=0)),
        "rms_collected_e_per_Heed_seed":float(coll.std(ddof=0)),
        "RMS_definition":"Population standard deviation (ddof=0), including zeros",
        "electron_histogram_axes":"linear x and linear y; all counts including zero",
        "mean_collected_e_per_Heed_seed_effective_gas_gain":mean_overall,
        "mean_created_e_ion_pairs_per_Heed_seed":mean_ion,
        "total_collected_e_div_total_created_pairs_NOT_a_gain":balance,
        "gain_definition":"Total anode-collected electron endpoints / total TrackHeed initial electrons",
        "GetAvalancheSize_ne_is_not_necessarily_anode_collected_charge":True,
        "angles":"Ion BIRTH phi=atan2(v,u-3mm*std::round(u/3mm)) around nearest anode",
        "angular_census":"Unbiased 128-point-per-avalanche ion-birth reservoir; weights recover ion-count total, but local bins approximate",
        "circular_mean_R_cut":0.2,
        "omitted_near_isotropic_centroids":missing_direction,
        "notes":"No ion drift, no cathode signal. Parent J/psi is independent physical event, not 588560 independent muons.",
        "plots_sha256":{p.name:sha(p) for p in OUT.glob("*.png")},
        "script_sha256":sha(Path(__file__)),
    }
    with (OUT/"summary.json").open("x") as f:
        json.dump(report,f,indent=2);f.write("\n")
    # Save angle bins and 30x30 gain grid for future re-analysis.
    with (OUT/"angular_distribution.csv").open("x",newline="") as f:
        w=csv.writer(f);w.writerow(("phi_bin_center_deg","ion_weighted_probability",
                                   "equal_nonzero_avalanche_probability",
                                   "circular_mean_direction_probability"))
        for j,c in enumerate(PHI_CENTRE):
            w.writerow((c,angles_ion[j]/angles_ion.sum(),
                        angles_equal_avalanche[j]/saved_cloud,
                        centroid[j]/max(1,centroid.sum())))
    with (OUT/"mean_gain_map.csv").open("x",newline="") as f:
        w=csv.writer(f);w.writerow(("u_bin_center_mm","v_bin_center_mm",
                                   "n_primary_electrons",
                                   "mean_collected_e_gain"))
        for i in range(30):
            for j in range(30):
                w.writerow(((Uedges[i]+Uedges[i+1])/2,
                            (Vedges[j]+Vedges[j+1])/2,int(count[i,j]),
                            meanG[i,j] if np.isfinite(meanG[i,j]) else ""))
    print("PASS 11,000 original J/psi p3_noFW physical muons",flush=True)
    print("FULL GARFIELD AVALANCHES",N,"ION BIRTHS",actual_births)
    print(f"MEAN gas gain (anode collected / real 1-e Heed seed): {mean_overall:.6g}")
    print(f"MEAN GetAvalancheSize(ne) per seed: {raw_ne.mean():.6g}")
    print(f"MEAN CREATED electron-ion pairs per seed: {mean_ion:.6g}")
    print(f"COLLECTED ELECTRONS / CREATED PAIRS: {balance:.6g} (NOT gain)")
    print("SAVED birth positions used:",sample_total)
    print("OUTPUT",OUT,flush=True)

def selftest():
    must(circular_phi(0.,1.)==90 and circular_phi(0.,-1.)==-90,
         "Real 3mm wire angle convention broken")
    must(abs(circular_phi(3.,1.)-90.)<1e-12,"Nearest anode folding broken")
    must(abs(circular_phi(-3.,1.)-90.)<1e-12,"Negative wire folding broken")
    assert len(PHI_CENTRE)==72
    print("PASS exact C++ nearest-anode geometry and 5-degree angular bins")
    print("NO physics, file writing or Garfield calls")

if __name__=="__main__":
    a=argparse.ArgumentParser(description=__doc__)
    a.add_argument("mode",choices=("selftest","plot","electrons-linear"))
    z=a.parse_args()
    if z.mode=="selftest":selftest()
    elif z.mode=="electrons-linear":plot_electrons_only()
    else:analyze()
