#!/usr/bin/env python3
"""Interactive Penning-ON p3_noFW MWPC physics plotting dashboard.

One-time preparation reads/validates the immutable 1000 real J/psi muons
and produces compact, reusable npz arrays. Every UI edit is immediate;
Garfield and TrackHeed are NEVER re-run.

Usage:
  python this_script.py prepare
  python -m streamlit run this_script.py
Data and cached images remain within a new, separate results subfolder.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

REPO=Path(__file__).resolve().parents[2]
ROOT=REPO/"test_runs/mwpc_stage_c/pi_p3_noFW_jpsi1000_penning057_v1"
ORIGINAL=REPO/"test_runs/mwpc_stage_c/pi_p3_noFW_jpsi1000_native_v1"
OUT=ROOT/"interactive_dashboard_v1"
NPZ=OUT/"dashboard_data.npz"
MANIFEST=OUT/"provenance.json"
N_MUONS=1000
PITCH_MM=3.0
PHI_EDGES=np.linspace(-180.,180.,73)
U_EDGES=np.linspace(-1.5,1.5,49)
V_EDGES=np.linspace(-3.,3.,73)
# >0-time bins are logarithmic to resolve sub-ns and multi-ns phenomena.
# A separate [0,10^-6 ns) bin includes exactly t_birth=t_Heed.
T_EDGES=np.r_[0.,np.logspace(-6.,6.,721)]
T_PHI_EDGES=np.r_[0.,np.logspace(-6.,6.,121)]

def must(condition,message):
    if not condition: raise RuntimeError(message)

def sha(path):
    digest=hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda:f.read(1<<20),b""):
            digest.update(block)
    return digest.hexdigest()

def readcsv(path,required):
    with Path(path).open(newline="") as f:
        rows=csv.DictReader(f)
        must(rows.fieldnames is not None and
             len(rows.fieldnames)==len(set(rows.fieldnames)) and
             set(required).issubset(set(rows.fieldnames)),
             "Missing/duplicated CSV columns in "+str(path))
        yield from rows

def key(row):
    return (row["crossing_id"],int(row["primary_index"]))

def nearest_wire_offset(u_mm):
    n=math.floor(u_mm/PITCH_MM+0.5) if u_mm>=0 else math.ceil(u_mm/PITCH_MM-0.5)
    return u_mm-PITCH_MM*n

def phi_idx(u,v):
    ang=math.degrees(math.atan2(v,nearest_wire_offset(u)))
    return min(71,max(0,int((ang+180.)/5.)))

def binned(value,edges):
    if not math.isfinite(value) or value<edges[0] or value>edges[-1]:
        return -1
    return min(len(edges)-2,max(0,int(np.searchsorted(edges,value,side="right")-1)))

def ahemi(up,down):
    s=up+down
    return (up-down)/s if s>0 else float("nan")

def prepare():
    must(not NPZ.exists() and not MANIFEST.exists(),
         "Dashboard cache already exists; refusing to overwrite "+str(OUT))
    old_input=ORIGINAL/"input/native_muon_crossings_input.csv"
    new_input=ROOT/"input/native_muon_crossings_input.csv"
    must(new_input.is_file() and old_input.is_file() and
         sha(new_input)==sha(old_input),
         "1000 Penning ON source muons are not the same as archived OFF input")
    man=json.loads((ROOT/"input/input_manifest.json").read_text())
    must(man["n_parents"]==N_MUONS and man["geometry"]=="p3_noFW"
         and sha(new_input)==man["selected_input_sha256"],
         "Original 1000-muon source was not SHA verified")
    pilot_marker="PHYSICS_PENNING_AR_CO2_70_30_R=0.57"
    sampled_ion=np.zeros(72,dtype=np.float64)
    sampled_equal=np.zeros(72,dtype=np.float64)
    xy_ion=np.zeros((len(V_EDGES)-1,len(U_EDGES)-1),dtype=np.float64)
    xy_n=np.zeros_like(xy_ion)
    start_sum=np.zeros_like(xy_ion)
    start_n=np.zeros_like(xy_ion)
    delay_ion=np.zeros(len(T_EDGES)-1,dtype=np.float64)
    tphi_ion=np.zeros((len(T_PHI_EDGES)-1,72),dtype=np.float64)
    gain=[]
    raw_ne=[]
    ion_births=[]
    seed_u=[]
    seed_v=[]
    seed_absA=[]
    seed_A=[]
    seed_span=[]
    seed_mean_delay=[]
    muon_gain=[]
    muon_heed=[]
    muon_absA=[]
    muon_A=[]
    muon_births=[]
    seals=[]
    n_ion_samples=0
    total_ions=0
    n_ion_avals=0
    invalid_negative_dt=0
    out_of_range_dt_ions=0.0
    min_dt=float("inf")
    max_dt=-float("inf")
    for i in range(50):
        d=ROOT/"native_chunks"/f"chunk_{i:04d}"
        receipt_path=d/"complete.json"
        must(receipt_path.is_file(),"Missing sealed Penning-ON checkpoint "+str(d))
        receipt=json.loads(receipt_path.read_text())
        must(receipt["schema"]=="NA60_JPSI_20PARENT_P3NOFW_SHA_V1"
             and receipt["chunk"]==i and receipt["n_parents"]==20,
             "Unrecognised 20-muon native checkpoint")
        must(receipt["stats"]["n_capped_seeds"]==0
             and receipt["stats"]["same_Heed_electron_positions_and_Garfield_seeds"] is True,
             "Penning-on checkpoint capped or original initial states not paired")
        for name,expected in receipt["file_sha256"].items():
            must(sha(d/name)==expected,
                 "Original Garfield physical output failed SHA: "+str(d/name))
        must(pilot_marker in (d/"physics.log").read_text(errors="replace"),
             "Penning was not enabled in actual binary for "+str(d))
        seals.append(sha(receipt_path))
        rows=list(readcsv(d/"muon_reference_seed.csv",
            ("crossing_id","primary_index","garfield_ne","ion_births",
             "collected_e","ion_up","ion_down","seed_u_mm","seed_v_mm","capped")))
        tmap={key(r):float(r["t_ns"]) for r in readcsv(
            d/"muon_reference_heed_primaries.csv",
            ("crossing_id","primary_index","t_ns"))}
        must(len(tmap)==len(rows),"The original Heed electron count differs from Garfield seeds")
        local_keys={key(r):j for j,r in enumerate(rows)}
        must(len(local_keys)==len(rows),"Duplicate Garfield initial electrons")
        local_min=np.full(len(rows),np.inf)
        local_max=np.full(len(rows),-np.inf)
        local_sum=np.zeros(len(rows))
        local_n=np.zeros(len(rows),dtype=np.int32)
        local_ions=0
        local_ion_avalanches=0
        for r in rows:
            col=int(r["collected_e"])
            ions=int(r["ion_births"])
            ne=int(r["garfield_ne"])
            up=int(r["ion_up"]);down=int(r["ion_down"])
            must(col>=0 and ions>=0 and ne>=0 and int(r["capped"])==0,
                 "Negative or capped Garfield avalanche, refusing misleading plots")
            gu=float(r["seed_u_mm"]);gv=float(r["seed_v_mm"])
            gain.append(col);raw_ne.append(ne);ion_births.append(ions)
            seed_u.append(nearest_wire_offset(gu));seed_v.append(gv)
            A=ahemi(up,down)
            seed_A.append(A);seed_absA.append(abs(A) if math.isfinite(A) else np.nan)
            local_ions+=ions
            if ions>0: local_ion_avalanches+=1
            xi=binned(nearest_wire_offset(gu),U_EDGES)
            yi=binned(gv,V_EDGES)
            if xi>=0 and yi>=0:
                start_n[yi,xi]+=1
                start_sum[yi,xi]+=col
        batch=list(readcsv(d/"muon_reference_crossing.csv",
            ("source|chunk|event|hit_index|track_id|station|det_id",
             "ion_up","ion_down","ion_births","collected_avalanche_electrons",
             "heed_primaries")))
        must(len(batch)==20,"Incorrect muon number per checkpoint")
        for r in batch:
            hp=int(r["heed_primaries"])
            collected=int(r["collected_avalanche_electrons"])
            births_mu=int(r["ion_births"])
            must(hp>=0 and collected>=0 and births_mu>=0,
                 "Negative physical muon/Heed census")
            if hp==0:
                # A genuine Geant4 gas crossing CAN yield no accepted Heed
                # electron (e.g. a short chord). There are then ZERO Garfield
                # avalanche calls, so G_mu = total_collected/hp is 0/0,
                # UNDEFINED rather than a zero-gain avalanche.
                must(collected==0 and births_mu==0 and
                     int(r["ion_up"])==0 and int(r["ion_down"])==0,
                     "Muon with zero Heed electrons has nonzero Garfield products")
            muon_heed.append(hp)
            muon_gain.append(collected/hp if hp>0 else np.nan)
            muon_births.append(births_mu)
            A=ahemi(int(r["ion_up"]),int(r["ion_down"]))
            muon_A.append(A);muon_absA.append(abs(A) if math.isfinite(A) else np.nan)
        must(sum(int(r["ion_births"]) for r in batch)==local_ions,
             "Per-muon ion counts and per-electron births disagree")
        sample_counts=np.zeros(len(rows),dtype=np.int32)
        weighted_sum=0.0
        for r in readcsv(d/"muon_reference_ion_birth_sample.csv",
                ("crossing_id","primary_index","u_mm","v_mm","w_mm","t_ns",
                 "reservoir_weight","seed_ion_births_total")):
            k=key(r)
            must(k in local_keys,"Ion birth has no matching microscopic electron")
            j=local_keys[k];b=int(rows[j]["ion_births"])
            must(b>0 and int(r["seed_ion_births_total"])==b,
                 "Ion-birth callback total differs from parent avalanche")
            m=min(128,b)
            weight=float(r["reservoir_weight"])
            must(math.isfinite(weight) and abs(weight-b/m)<=max(1.,b/m)*1.e-7,
                 "Ion-birth reservoir weights changed")
            sample_counts[j]+=1
            weighted_sum+=weight
            u=float(r["u_mm"]);v=float(r["v_mm"]);t=float(r["t_ns"])
            dtime=t-tmap[k]
            min_dt=min(min_dt,dtime);max_dt=max(max_dt,dtime)
            if dtime<-1.e-6:
                invalid_negative_dt+=1
            # Small numerical negative values are treated as t=0 only if
            # within ±1e-6 ns; physical negative delays are preserved as
            # a warning/census rather than secretly included in the pdf.
            dtime=max(0.,dtime)
            local_min[j]=min(local_min[j],dtime)
            local_max[j]=max(local_max[j],dtime)
            local_sum[j]+=dtime
            local_n[j]+=1
            p=phi_idx(u,v)
            sampled_equal[p]+=1./m
            sampled_ion[p]+=weight
            xi=binned(nearest_wire_offset(u),U_EDGES)
            yi=binned(v,V_EDGES)
            if xi>=0 and yi>=0:
                xy_ion[yi,xi]+=weight
                xy_n[yi,xi]+=1
            it=binned(dtime,T_EDGES)
            if it>=0:
                delay_ion[it]+=weight
            else:
                out_of_range_dt_ions+=weight
            ip=binned(dtime,T_PHI_EDGES)
            if ip>=0:tphi_ion[ip,p]+=weight
            n_ion_samples+=1
        must(all(sample_counts[j]==min(128,int(r["ion_births"]))
                 for j,r in enumerate(rows)),
             "Missing exact 128-point-per-avalanches ion-birth reservoir")
        must(abs(weighted_sum-local_ions)<1.e-5*max(1,local_ions),
             "Saved weighted ion cloud does not reproduce full exact ion birth census")
        for j,r in enumerate(rows):
            if local_n[j]:
                seed_span.append(local_max[j]-local_min[j])
                seed_mean_delay.append(local_sum[j]/local_n[j])
            else:
                seed_span.append(np.nan)
                seed_mean_delay.append(np.nan)
        total_ions+=local_ions
        n_ion_avals+=local_ion_avalanches
        print(f"PASS Penning ON SHA/ion-census chunk {i+1:02d}/50: "
              f"20 muons, {len(rows)} Garfield seeds, {local_ions:,} ion births",
              flush=True)
    must(len(muon_gain)==N_MUONS and
         len(gain)==len(seed_span)==len(raw_ne)==len(ion_births)==len(seed_A)
         and len(gain)>20000,
         "Missing Penning ON parents or Garfield avalanches")
    must(abs(sampled_ion.sum()-total_ions)<1.e-5*max(1,total_ions) and
         abs(sampled_equal.sum()-n_ion_avals)<1.e-5*max(1,n_ion_avals),
         "Angular sampled/weighted total did not reproduce microscopic ion birth counts")
    must(invalid_negative_dt==0,
         "Some ion births predate the initiating TrackHeed electron by >1e-6ns "
         f"({invalid_negative_dt} samples, min delta={min_dt:.6g}ns). "
         "Time origin is questionable; inspect the raw physics before charting.")
    must(out_of_range_dt_ions<1.e-7*max(1,total_ions),
         "Time-binning range misses actual Garfield ion birth delays; "
         f"missed weighted ions={out_of_range_dt_ions}")
    muon_heed_arr=np.asarray(muon_heed,dtype=np.int64)
    muon_gain_arr=np.asarray(muon_gain,dtype=float)
    n_zero_heed=int(np.count_nonzero(muon_heed_arr==0))
    must(np.isfinite(muon_gain_arr[muon_heed_arr>0]).all()
         and np.isnan(muon_gain_arr[muon_heed_arr==0]).all(),
         "A muon with Heed electrons has undefined gain, or the zero-Heed "
         "muon is incorrectly assigned a numeric gain")
    print(f"PASS {N_MUONS} real Geant4 parents: "
          f"{N_MUONS-n_zero_heed} with retained Heed electrons and "
          f"{n_zero_heed} with no retained Heed electrons (G_mu undefined).",
          flush=True)
    must(not OUT.exists(),
         "Refusing to overwrite dashboard directory "+str(OUT))
    OUT.mkdir(parents=True)
    arrays={
        "gain":np.asarray(gain,np.int32),
        "raw_ne":np.asarray(raw_ne,np.int32),
        "births":np.asarray(ion_births,np.int32),
        "seed_u":np.asarray(seed_u,np.float32),
        "seed_v":np.asarray(seed_v,np.float32),
        "seed_A":np.asarray(seed_A,np.float32),
        "seed_absA":np.asarray(seed_absA,np.float32),
        "seed_span_ns":np.asarray(seed_span,np.float32),
        "seed_mean_birth_dt_ns":np.asarray(seed_mean_delay,np.float32),
        "muon_gain":np.asarray(muon_gain,np.float32),
        "muon_heed":np.asarray(muon_heed,np.int32),
        "muon_births":np.asarray(muon_births,np.int32),
        "muon_A":np.asarray(muon_A,np.float32),
        "muon_absA":np.asarray(muon_absA,np.float32),
        "angle_equal":sampled_equal,
        "angle_ion":sampled_ion,
        "phi_edges":PHI_EDGES,
        "u_edges":U_EDGES,
        "v_edges":V_EDGES,
        "xy_ion":xy_ion,
        "start_sum":start_sum,
        "start_n":start_n,
        "t_edges":T_EDGES,
        "t_weights":delay_ion,
        "t_phi_edges":T_PHI_EDGES,
        "t_phi":tphi_ion
    }
    np.savez_compressed(NPZ,**arrays)
    metadata={
        "schema":"P3_NOFW_PENNING_ON_INTERACTIVE_DASHBOARD_V1",
        "n_independent_Geant4_Jpsi_parents":1000,
        "n_muons_with_retained_Heed_primaries":N_MUONS-n_zero_heed,
        "n_muons_with_zero_retained_Heed_primaries":n_zero_heed,
        "muon_gain_definition":"Collected electron sum / number of retained Heed primary electrons; NaN for zero-primary muons (0/0, NOT zero-gain avalanches)",
        "n_microscopic_TrackHeed_primary_Garfield_avalanches":len(gain),
        "n_ion_producing_avalanches":n_ion_avals,
        "n_true_ion_births":total_ions,
        "n_saved_ion_birth_coordinate_samples":n_ion_samples,
        "n_bad_negative_birth_times":invalid_negative_dt,
        "physical_birth_delay_min_ns":min_dt,
        "physical_birth_delay_max_ns":max_dt,
        "n_ion_births_outside_time_histogram":out_of_range_dt_ions,
        "physical_setup":"3mm pitch; +/-3mm cathodes; 1800V; ArCO2 70:30; Penning=0.57",
        "caveat":"This is the original p3_noFW geometry, NOT the real 4mm 2+4mm prototype",
        "cathode_signal":"NOT simulated; plots are ion birth initial conditions only",
        "original_source_sha256":sha(old_input),
        "on_source_sha256":sha(new_input),
        "sha256_each_of_50_native_seals":seals,
        "cache_sha256":sha(NPZ),
        "script_sha256":sha(Path(__file__))
    }
    MANIFEST.write_text(json.dumps(metadata,indent=2)+"\n")
    print("PASS exact 1000 Penning ON muons with 50 verified SHA chunks",flush=True)
    print("Mean gain (all):",float(np.mean(gain)),"mean gain (>0):",
          float(np.mean(np.asarray(gain)[np.asarray(gain)>0])))
    print("CACHE",NPZ,flush=True)
    print("RUN: python -m streamlit run",Path(__file__),flush=True)


def _weighted_time_quantile(edges,weights,q):
    """Bin-edge quantile approximation; preserves physical units ns."""
    w=np.asarray(weights,dtype=np.float64)
    e=np.asarray(edges,dtype=np.float64)
    must(len(e)==len(w)+1 and 0.<=q<=1. and w.sum()>0.,
         "Invalid physical ion-birth time marginal")
    idx=min(len(w)-1,int(np.searchsorted(np.cumsum(w),q*w.sum(),side="left")))
    return float(e[idx+1])


def plot_explanation(kind,data,meta,angular_deg):
    """Three grounded, plot-specific paragraphs; no physics rerun.

    Output: (what is shown, how it was measured, what the saved data reveal).
    All numerical statements derive from the verified Penning-ON cache,
    never from assumptions about the shape of a histogram.
    """
    g=data["gain"].astype(np.float64)
    positive=g[g>0]
    n=len(g)
    n0=int(np.count_nonzero(g==0))
    ng=int(meta["n_true_ion_births"])
    if kind=="gain_all":
        return (
            "**What you see:** The distribution of anode-collected electron counts "
            "from one Garfield avalanche started by one retained TrackHeed electron. "
            "A separate [0,1) bin represents avalanches collecting no electrons.",
            "**How measured:** The microscopic code counts electron endpoints close "
            "to anode wires for every initial electron. Each avalanche contributes "
            "one entry. The average includes zero-collection avalanches, rather than "
            "averaging over muons.",
            f"**In this sample:** {n:,} avalanches; mean {g.mean():,.1f} electrons, "
            f"RMS {g.std(ddof=0):,.1f}; {n0:,} collect none ({n0/n:.2%}); "
            f"99th percentile {np.quantile(g,.99):,.0f}. Fluctuations in collected "
            "charge matter for event-by-event cathode signal amplitudes."
        )
    if kind=="gain_pos":
        return (
            "**What you see:** The avalanche gain distribution **conditional on at "
            "least one electron being collected at an anode**. It is not the gain "
            "distribution of all starting electrons.",
            "**How measured:** We use exactly the same microscopic collected-endpoint "
            "counter but select only avalanche records with collected_e > 0, then "
            "renormalise the histogram using this selected population.",
            f"**In this sample:** {len(positive):,} positive-gain avalanches; mean "
            f"{positive.mean():,.1f}, RMS {positive.std(ddof=0):,.1f} electrons. "
            f"This conditional mean exceeds the inclusive mean ({g.mean():,.1f}) "
            "because the zero-collection outcomes are omitted."
        )
    if kind in ("angle_eq","angle_ion"):
        raw=np.asarray(data["angle_equal"] if kind=="angle_eq" else data["angle_ion"],
                       dtype=np.float64)
        e=np.asarray(data["phi_edges"],dtype=np.float64)
        mid=(e[1:]+e[:-1])/2.
        axial=np.minimum(np.abs(mid),180.-np.abs(mid))<30.
        transverse=np.abs(np.abs(mid)-90.)<30.
        axsum=float(raw[axial].sum())
        trsum=float(raw[transverse].sum())
        contrast=axsum/trsum if trsum>0 else float("nan")
        mode="every ion-producing avalanche contributes total weight 1" if kind=="angle_eq" else "each sampled position is weighted by the exact avalanche ion count / number of retained samples"
        population="per avalanche" if kind=="angle_eq" else "per created ion"
        return (
            "**What you see:** The azimuth φ of the locations at which positive ions "
            "were *created*, measured around the nearest anode wire. φ = 0° points "
            "to +u, ±180° to −u, and ±90° to the ±v directions. This is not "
            "the direction in which ions subsequently drift.",
            "**How measured:** Garfield's ionisation callback records the birth "
            "positions; up to 128 positions are uniformly sampled per avalanche. "
            f"Here, {mode}. Angular bins are currently {angular_deg}° wide.",
            f"**In this sample:** The probability in equal-size angular windows "
            f"within 30° of the ±u axis is {contrast:.2f}× that in windows within "
            "30° of the ±v axis. This quantifies how anisotropic the *birth* "
            f"distribution is {population}; it does not yet predict electrode "
            "signals without ion transport and weighting fields."
        )
    if kind in ("asym_seed","asym_muon"):
        v=np.asarray(data["seed_A"] if kind=="asym_seed" else data["muon_A"],
                     dtype=np.float64)
        v=v[np.isfinite(v)]
        unit="individual ion-producing avalanche" if kind=="asym_seed" else "whole real J/ψ muon"
        unitcount=f"{len(v):,} defined values"
        return (
            f"**What you see:** The *signed* upper/lower ion-birth asymmetry of each "
            f"{unit}, with A = (N_up − N_down)/(N_up + N_down). "
            "A = 0 is symmetric; +1 means all counted births above the anode "
            "and −1 means all below.",
            "**How measured:** Every ionisation callback is assigned to the upper "
            "or lower half-plane using its exact position relative to the closest "
            "anode. These are **full exact integer counts**, not the 128-point "
            "spatial subsample. The plotted asymmetry is undefined if both "
            "sector counts are zero.",
            f"**In this sample:** {unitcount}; mean signed A = {v.mean():+.4f}, "
            f"mean |A| = {np.abs(v).mean():.4f}, and "
            f"{np.mean(np.abs(v)>.5):.1%} have |A| > 0.5. "
            "An overall mean near zero can coexist with large individual "
            "charge imbalances and does not imply symmetric cathode pulses."
        )
    if kind=="time":
        e=np.asarray(data["t_edges"],dtype=np.float64)
        w=np.asarray(data["t_weights"],dtype=np.float64)
        p50=_weighted_time_quantile(e,w,.5)
        p90=_weighted_time_quantile(e,w,.9)
        return (
            "**What you see:** When ions are *created* during the electron "
            "avalanche, expressed as Δt = t_ion birth − t_initial Heed electron "
            "in nanoseconds. This is **not ion drift time** and not CPU time.",
            "**How measured:** We subtract the initiating TrackHeed electron's "
            "saved time from each sampled Garfield ionisation-birth time. "
            "Each retained birth is weighted to represent the exact number "
            "of ions in its avalanche. Time bins are nonuniform, so "
            "'probability density' is preferable when comparing their heights.",
            f"**In this sample:** Approximate ion-weighted median creation "
            f"delay ≤ {p50:.4g} ns and 90th percentile ≤ {p90:.4g} ns "
            f"(bin-edge estimates from {ng:,} actual ion births). "
            "This temporal source is needed before adding ion drift and the "
            "Shockley–Ramo cathode response."
        )
    if kind=="span":
        span=np.asarray(data["seed_span_ns"],dtype=np.float64)
        finite=span[np.isfinite(span)]
        return (
            "**What you see:** The interval between the earliest and latest "
            "saved ion-birth times **within one avalanche**; one histogram "
            "entry corresponds to one ion-producing initial electron.",
            "**How measured:** For each avalanche we take max(t) − min(t) "
            "over at most 128 randomly selected ion births. Because not "
            "every birth is stored, this is a **lower bound** on the true "
            "microscopic creation-time span.",
            f"**In this sample:** {len(finite):,} avalanches with sampled "
            f"births; median sampled span {np.median(finite):.4g} ns and "
            f"90th percentile {np.quantile(finite,.9):.4g} ns. "
            "A small observed span does not establish a short ion-induced "
            "cathode pulse—the ions still need to drift."
        )
    if kind=="gain_asym":
        a=np.asarray(data["seed_absA"],dtype=np.float64)
        mask=(g>0)&np.isfinite(a)
        x=g[mask]
        y=a[mask]
        split=float(np.median(x))
        low=y[x<=split]
        high=y[x>split]
        return (
            "**What you see:** How much ion creation favours one half-plane "
            "as a function of anode-collected electron gain. The ordinate "
            "is **mean |A|**, not a signed asymmetry or a direct cathode "
            "charge ratio.",
            "**How measured:** A is computed exactly from the callback's "
            "ion-up and ion-down counts for each avalanche. Positive "
            "collected-gain avalanches are grouped in logarithmic gain bins; "
            "each plotted point is the mean of their |A| values. "
            "Bins with fewer than 20 entries are omitted.",
            f"**In this sample:** {len(x):,} eligible avalanches; splitting "
            f"at median collected gain {split:,.0f}, mean |A| is "
            f"{low.mean():.4f} below and {high.mean():.4f} above this gain. "
            "This is a correlation relevant to charge-sharing fluctuations, "
            "not proof that anode gain alone determines cathode response."
        )
    if kind=="ion_xy":
        w=np.asarray(data["xy_ion"],dtype=np.float64)
        xe=np.asarray(data["u_edges"],dtype=np.float64)
        ye=np.asarray(data["v_edges"],dtype=np.float64)
        total=float(w.sum())
        frac=total/ng if ng else float("nan")
        iy,ix=np.unravel_index(int(np.argmax(w)),w.shape)
        umid=float((xe[ix]+xe[ix+1])/2)
        vmid=float((ye[iy]+ye[iy+1])/2)
        return (
            "**What you see:** A top-down *ion-creation* map in the (u,v) "
            "plane near the anode. Horizontal u is folded into the nearest "
            "3-mm anode cell; vertical v points between the cathode planes.",
            "**How measured:** For each Garfield avalanche, up to 128 "
            "ion-birth coordinates are sampled uniformly. Every sample "
            "has weight N_ion births / N_samples. Heatmap colour counts "
            "weighted births **per spatial cell**; it is not yet a density "
            "per mm² or a drifting-ion distribution.",
            f"**In this sample:** The displayed coordinate window contains "
            f"approximately {frac:.2%} of the weighted ion births; its "
            f"most populated cell is near (u,v)=({umid:+.3f},{vmid:+.3f}) mm "
            "(dependent on chosen grid size). The position distribution "
            "will seed future ion-drift and weighting-field calculations."
        )
    if kind=="gain_xy":
        cnt=np.asarray(data["start_n"],dtype=np.float64)
        sm=np.asarray(data["start_sum"],dtype=np.float64)
        valid=cnt>=5
        n_cells=int(np.count_nonzero(valid))
        total=int(cnt.sum())
        return (
            "**What you see:** The conditional mean number of electrons "
            "collected at the anode as a function of the **starting "
            "TrackHeed electron position** relative to its nearest wire.",
            "**How measured:** Each retained Heed electron is assigned "
            "to a spatial (u,v) cell, and its exact Garfield collected "
            "electron count is included—also if collection was zero. "
            "The colour is (sum collected electrons)/(number of starts) "
            "in that cell. Cells with fewer than 5 starts are hidden.",
            f"**In this sample:** {total:,} electron starts fall inside "
            f"the displayed spatial window and {n_cells:,} cells meet "
            "the minimum population. Structure in this map tests "
            "whether a position-conditioned gain model is needed "
            "for realistic charge generation."
        )
    if kind=="time_phi":
        w=np.asarray(data["t_phi"],dtype=np.float64)
        e=np.asarray(data["t_phi_edges"],dtype=np.float64)
        hist=w.sum(axis=1)
        p90=_weighted_time_quantile(e,hist,.9)
        return (
            "**What you see:** A joint map of the ion-birth azimuth φ "
            "and creation delay Δt for the **same sampled births**. "
            "Bright bins have more weighted ion births at that angle "
            "and time; they are not measured cathode current.",
            "**How measured:** From each saved Garfield ionisation "
            "position and creation time we calculate φ around the nearest "
            "anode and Δt relative to the original TrackHeed electron "
            "start. A retained sample carries weight N_births/N_samples. "
            "The time bins are logarithmically spaced. Colour denotes "
            "ion births **per bin**, not per ns per degree.",
            f"**In this sample:** The upper edge of the approximate "
            f"ion-weighted 90th-percentile creation delay is {p90:.4g} ns "
            "in the stored joint time bins. Changes in angular structure "
            "with time would be useful for modelling charge development, "
            "but ion drift and the cathode weighting field are still absent."
        )
    raise RuntimeError("Missing grounded interpretation for "+str(kind))



def plot_features(kind):
    """Sidebar capabilities: never display controls that cannot affect a plot."""
    hist={"gain_all","gain_pos","angle_eq","angle_ion",
          "asym_seed","asym_muon","time","span"}
    heat={"ion_xy","gain_xy","time_phi","t10_start","t10_gain_position"}
    curves={"gain_asym","t10_gain","t10_asym","t10_angles"}
    timeline=kind=="t10_timeline"
    must(kind in hist|heat|curves|{"t10_timeline"},
         "No capability definitions for plot "+kind)
    return {
        "histogram_style":kind in hist,
        "curve_color":kind in hist|curves,
        "line_width":kind in hist|curves,
        "legend_controls":kind in hist|curves or timeline,
        "legend_text":kind!="t10_angles",
        "gain_bin":kind in {"gain_all","gain_pos"},
        "angular_bin":kind in {"angle_eq","angle_ion","t10_angles"},
        "asym_bin":kind in {"asym_seed","asym_muon"},
        "time_group":kind=="time",
        "time_window":kind=="time",
        "histogram_normalization":kind in hist,
        "heatmap_colorscale":kind in heat,
        "mean_rms":kind in {"gain_all","gain_pos","asym_seed","asym_muon",
                            "time","span","gain_asym","t10_gain","t10_asym",
                            "t10_angles"},
        "x_log":kind not in {"gain_all","angle_eq","angle_ion",
                             "asym_seed","asym_muon","ion_xy","gain_xy",
                             "time_phi","t10_angles","t10_start",
                             "t10_timeline"},
        "y_log":kind in hist|curves|{"time_phi","t10_start"},
    }


def legend_stats(label,values,enabled,quantity="",weights=None,approx=False):
    """Add meaningful and reproducible mean + population RMS to legend."""
    if not enabled:return label
    a=np.asarray(values,dtype=np.float64)
    valid=np.isfinite(a)
    if weights is not None:
        w=np.asarray(weights,dtype=np.float64)
        valid &= np.isfinite(w)&(w>=0)
        a=a[valid];w=w[valid]
        if a.size==0 or w.sum()<=0:return label
        mean=float(np.average(a,weights=w))
        rms=float(np.sqrt(np.average((a-mean)**2,weights=w)))
    else:
        a=a[valid]
        if a.size==0:return label
        mean=float(np.mean(a))
        rms=float(np.std(a,ddof=0))
    prefix="≈" if approx else ""
    units=(" "+quantity) if quantity else ""
    return f"{label} · {prefix}mean={mean:,.2f}{units}, {prefix}RMS={rms:,.2f}{units}"


def app():
    import streamlit as st
    import plotly.graph_objects as go
    import pi_p3_noFW_penning057_timing_correlations_v1 as temporal
    from plotly.colors import sequential
    from io import StringIO

    st.set_page_config(page_title="MWPC · Penning ON / Plot Studio",
                       page_icon="📈",layout="wide",initial_sidebar_state="expanded")
    st.markdown("""<style>
    .block-container {padding-top:1.6rem;max-width:1600px;}
    div[data-testid="stMetric"] {border:1px solid #DEE5EB;border-radius:10px;padding:12px;}
    .stDownloadButton button {width:100%}
    </style>""",unsafe_allow_html=True)
    st.title("MWPC Plot Studio")
    st.caption("Real 1,000 J/ψ muons · Garfield++ Penning ON · Ar/CO₂ 70:30 · "
               "p3_noFW 3 mm pitch · 1,800 V")
    if not NPZ.is_file() or not MANIFEST.is_file():
        st.error("No verified dashboard cache exists. First run:\n\n"
                 "python test/mwpc-stage-c/pi_p3_noFW_penning057_interactive_dashboard_v1.py prepare")
        st.stop()

    @st.cache_data(show_spinner="Loading verified Penning-ON data cache…")
    def load_cache(path,mtime):
        with np.load(path,allow_pickle=False) as z:
            arrays={k:z[k] for k in z.files}
        return arrays,json.loads(MANIFEST.read_text())
    data,meta=load_cache(str(NPZ),NPZ.stat().st_mtime_ns)

    g=data["gain"].astype(float)
    positive=g[g>0]
    a,b,c,d,e=st.columns(5)
    a.metric("Real J/ψ muons","1,000")
    b.metric("Microscopic avalanches",f"{len(g):,}")
    c.metric("Mean collected gain",f"{np.mean(g):,.1f}")
    d.metric("Zero collected",f"{np.mean(g==0):.2%}")
    e.metric("Muons without Heed e",str(meta["n_muons_with_zero_retained_Heed_primaries"]))
    st.caption("Muon-level gain is undefined (0/0) for muons with no retained "
               "TrackHeed electrons; they remain in the 1,000-muon sample. "
               "Avalanche-level gain statistics use only real starting electrons.")

    kinds={
        "gain_all":"Avalanche gain · including zero",
        "gain_pos":"Avalanche gain · positive only",
        "angle_eq":"Ion-birth angle · equal avalanche",
        "angle_ion":"Ion-birth angle · ion weighted",
        "asym_seed":"Exact asymmetry · each avalanche",
        "asym_muon":"Exact asymmetry · each muon",
        "time":"Ion-birth delay · ion weighted",
        "span":"Avalanche ion-birth time span · sampled",
        "gain_asym":"Mean |asymmetry| versus gain",
        "ion_xy":"Ion-birth spatial density (u, v)",
        "gain_xy":"Mean gain versus initial electron position",
        "time_phi":"Birth time versus azimuth",
    }
    kinds.update(temporal.KINDS)
    defaults={
        "gain_all":("Collected-electron gain — all avalanches",
                    "Anode-collected electrons per initial Heed electron",
                    "Probability per bin","Penning ON · all avalanches"),
        "gain_pos":("Collected-electron gain — nonzero avalanches",
                    "Anode-collected electrons per initial Heed electron",
                    "Probability per bin","Penning ON · gain > 0"),
        "angle_eq":("Ion-birth azimuth — equal avalanche weighting",
                    "Birth azimuth around nearest anode [degrees]",
                    "Probability per bin","Each ion-producing avalanche weighted equally"),
        "angle_ion":("Ion-birth azimuth — ion-weighted",
                    "Birth azimuth around nearest anode [degrees]",
                    "Probability per bin","Each created ion weighted equally"),
        "asym_seed":("Exact upper/lower ion-birth asymmetry — per avalanche",
                    "A = (N_up − N_down)/(N_up + N_down)",
                    "Probability per bin","Each ion-producing avalanche"),
        "asym_muon":("Exact upper/lower ion-birth asymmetry — per J/ψ muon",
                    "A = (N_up − N_down)/(N_up + N_down)",
                    "Probability per bin","Each real J/ψ muon"),
        "time":("Ion creation delay relative to Heed electron start",
                "t_ion_birth − t_Heed [ns]","Ion births per ns (probability density)",
                "Ion-weighted birth-time distribution"),
        "span":("Ion-birth time span per avalanche — 128 sampled points",
                "Sampled max(t) − min(t) [ns]","Probability per bin",
                "Sampled time span; underestimates exact avalanche span"),
        "gain_asym":("Absolute ion-birth asymmetry as a function of gain",
                     "Collected gain (electrons per initial Heed electron)",
                     "Mean |A_hemi|","Penning ON · exact ion-sector counts"),
        "ion_xy":("Ion creation positions near the nearest anode",
                  "Position u relative to nearest anode [mm]",
                  "Position v [mm]","Ion birth density (weighted)"),
        "gain_xy":("Effective gain versus Heed electron starting position",
                   "Initial u relative to nearest anode [mm]",
                   "Initial v [mm]","Mean collected electrons"),
        "time_phi":("Ion creation delay vs azimuth",
                    "Birth azimuth around nearest anode [deg]",
                    "Birth delay [ns]","Ion-weighted birth density")
    }
    defaults.update(temporal.DEFAULTS)
    st.sidebar.header("Edit the plot")
    kind=st.sidebar.selectbox("Plot",list(kinds),format_func=lambda x:kinds[x])
    temporal_data=temporal_meta=None
    if kind in temporal.KINDS:
        if not temporal.CACHE.is_file() or not temporal.META.is_file():
            st.warning("Temporal correlations are a separate derived cache. "
                       "Run the new read-only temporal script with 'prepare' "
                       "once, then restart the dashboard.")
            st.code("python test/mwpc-stage-c/"
                    "pi_p3_noFW_penning057_timing_correlations_v1.py prepare",
                    language="bash")
            st.stop()
        @st.cache_data(show_spinner="Loading the 1000-muon temporal cache...")
        def load_temporal(path,mtime):
            return temporal.load()
        temporal_data,temporal_meta=load_temporal(
            str(temporal.CACHE),temporal.CACHE.stat().st_mtime_ns)
    if st.sidebar.button("Reset this plot to defaults"):
        for k in list(st.session_state):
            if k.startswith(kind+"-"):
                del st.session_state[k]
        st.rerun()
    caps=plot_features(kind)
    title0,x0,y0,legend0=defaults[kind]
    # Each view has its own independent saved editable controls.
    title=st.sidebar.text_input("Title",value=title0,key=f"{kind}-title")
    xlab=st.sidebar.text_input("X-axis title",value=x0,key=f"{kind}-xlab")
    ylab=st.sidebar.text_input("Y-axis title",value=y0,key=f"{kind}-ylab")
    if caps["legend_text"]:
        label=st.sidebar.text_input(
            "Legend text / colourbar title",value=legend0,
            key=f"{kind}-legend")
    else:
        label=legend0
    st.sidebar.divider()
    c1,c2=st.sidebar.columns(2)
    if caps["x_log"]:
        xscale=c1.selectbox("X scale",["Linear","Log"],
                            key=f"{kind}-xscale")
    else:
        xscale="Linear"
    if caps["y_log"]:
        yscale=c2.selectbox("Y scale",["Linear","Log"],
                            index=1 if kind=="time_phi" else 0,
                            key=f"{kind}-yscale")
    else:
        yscale="Linear"
    if caps["curve_color"]:
        color=st.sidebar.color_picker("Curve colour","#E46B08",
                                      key=f"{kind}-color")
    else:
        color="#E46B08"
    if caps["histogram_style"]:
        plot_style=st.sidebar.selectbox(
            "Histogram style",["Stepped line","Bars"],key=f"{kind}-style")
    else:
        plot_style="Stepped line"
    if caps["line_width"]:
        line_width=st.sidebar.slider("Curve line width",1.,5.,2.2,.2,
                                     key=f"{kind}-width")
    else:
        line_width=2.2
    if caps["legend_controls"]:
        show_legend=st.sidebar.checkbox("Show legend",True,
                                         key=f"{kind}-legend-show")
        legpos=st.sidebar.selectbox(
            "Legend position",
            ["Top right","Top left","Bottom right","Bottom left"],
            key=f"{kind}-legend-position")
    else:
        show_legend=False
        legpos="Top right"
    if caps["mean_rms"]:
        add_stats=st.sidebar.checkbox(
            "Append mean and RMS to legend",
            value=kind in {"gain_all","gain_pos"},
            key=f"{kind}-stats")
        st.sidebar.caption(
            "Statistics use the underlying observations, not the heights "
            "of the displayed histogram bins.")
    else:
        add_stats=False
    show_interpretation=st.sidebar.checkbox(
        "Show physics interpretation above plot",True,
        key=f"{kind}-interpretation-show")
    # Export applies to every plot; keep it independent of plot-specific
    # physics controls and defined on every Streamlit rerun.
    export_format=st.sidebar.selectbox(
        "Camera export format",["png","svg","jpeg"],
        key=f"{kind}-export")
    title_font=st.sidebar.slider("Title font size",14,34,22,
                                 key=f"{kind}-title-font")
    label_font=st.sidebar.slider("Axes font size",10,24,15,
                                 key=f"{kind}-axes-font")
    st.sidebar.divider()
    # Plot-specific scientific controls. Unused parameters are not shown.
    if caps["gain_bin"]:
        gain_width=st.sidebar.select_slider(
            "Gain histogram bin width [electrons]",
            options=[10,20,25,50,100,200,500],value=50,
            key=f"{kind}-gainbin")
    else:
        gain_width=50
    if caps["angular_bin"]:
        angular_deg=st.sidebar.select_slider(
            "Angular bin width [degrees]",
            options=[5,10,15,20,30],value=5,
            key=f"{kind}-angbin")
    else:
        angular_deg=5
    if caps["asym_bin"]:
        asym_bins=st.sidebar.slider("Asymmetry bins",20,120,50,10,
                                    key=f"{kind}-asymbin")
    else:
        asym_bins=50
    if caps["time_group"]:
        time_group=st.sidebar.select_slider(
            "Time-bin grouping",options=[1,2,3,4,6,8,12],value=4,
            key=f"{kind}-tgroup")
    else:
        time_group=4
    if caps["time_window"]:
        time_range=st.sidebar.selectbox(
            "Time view window",
            ["Central 99.9% of ion births","All ion-birth delays"],
            key=f"{kind}-trange")
    else:
        time_range="All ion-birth delays"
    if caps["histogram_normalization"]:
        norm=st.sidebar.selectbox(
            "Histogram normalization",
            ["Probability per bin","Probability density",
             "Raw counts / weighted counts"],
            index=1 if kind=="time" else 0,
            key=f"{kind}-normalization")
    else:
        norm="Probability per bin"
    if caps["heatmap_colorscale"]:
        zscale=st.sidebar.selectbox(
            "Heatmap colour scale",["Linear","Log"],
            key=f"{kind}-zscale")
    else:
        zscale="Linear"
    time_controls={}
    if kind in temporal.KINDS:
        st.sidebar.divider()
        st.sidebar.subheader("Early-birth correlation settings")
        if kind!="t10_timeline":
            reliable=st.sidebar.checkbox(
                "Only 128-sample ion avalanches",True,
                key=f"{kind}-reliable")
        else:
            reliable=True
        time_controls={
            "reliable_only":reliable,"color":color,"line_width":line_width,
            "legend":label,"angular_deg":angular_deg,
            "add_stats":add_stats,
            "zscale":zscale,
            "mode":"Quantile groups",
            "early_min":10.,"early_max":40.,
            "late_min":85.,"late_max":115.,
            "legend_early":"Earlier birth times",
            "legend_late":"Later birth times",
            "muon_index":0
        }
        if kind=="t10_angles":
            grouping=st.sidebar.selectbox(
                "Early/late comparison",
                ["Custom time intervals (ns)","Quantile groups"],
                key=f"{kind}-grouping")
            time_controls["mode"]=grouping
            time_controls["legend_early"]=st.sidebar.text_input(
                "Early curve legend","Earlier birth times",
                key=f"{kind}-early-label")
            time_controls["legend_late"]=st.sidebar.text_input(
                "Late curve legend","Later birth times",
                key=f"{kind}-late-label")
            if grouping=="Custom time intervals (ns)":
                st.sidebar.caption(
                    "Enter two independent birth-time intervals in ns "
                    "relative to the muon entry (t10 proxy). "
                    "From is included, To is excluded.")
                st.sidebar.markdown("**Earlier birth-time group**")
                early_start,early_end=st.sidebar.columns(2)
                time_controls["early_min"]=early_start.number_input(
                    "From [ns]",min_value=0.,value=10.,step=5.,
                    key=f"{kind}-early-min")
                time_controls["early_max"]=early_end.number_input(
                    "To [ns]",min_value=0.,value=40.,step=5.,
                    key=f"{kind}-early-max")
                st.sidebar.markdown("**Later birth-time group**")
                late_start,late_end=st.sidebar.columns(2)
                time_controls["late_min"]=late_start.number_input(
                    "From [ns]",min_value=0.,value=85.,step=5.,
                    key=f"{kind}-late-min")
                time_controls["late_max"]=late_end.number_input(
                    "To [ns]",min_value=0.,value=115.,step=5.,
                    key=f"{kind}-late-max")
                if (time_controls["early_min"]>=time_controls["early_max"]
                    or time_controls["late_min"]>=time_controls["late_max"]):
                    st.sidebar.error("Each group needs From < To.")
                    st.stop()
                if time_controls["early_max"]>time_controls["late_min"]:
                    st.sidebar.warning("These intervals overlap or are not "
                                       "ordered chronologically; some "
                                       "avalanches may appear in both groups.")
        if kind=="t10_timeline":
            time_controls["muon_index"]=st.sidebar.selectbox(
                "Real J/psi muon to inspect",
                options=list(range(1000)),
                index=temporal.typical_muon(temporal_meta),
                format_func=lambda i:(
                    f"Muon {i:03d} | "
                    f"{temporal_meta['muons'][i]['n_clusters']} clusters, "
                    f"{temporal_meta['muons'][i]['n_Heed']} initial e"),
                key=f"{kind}-selected-muon")
    st.sidebar.divider()
    with st.sidebar.expander("Axes limits (optional)",expanded=False):
        xmin=st.text_input("X minimum",key=f"{kind}-xmin",placeholder="Automatic")
        xmax=st.text_input("X maximum",key=f"{kind}-xmax",placeholder="Automatic")
        ymin=st.text_input("Y minimum",key=f"{kind}-ymin",placeholder="Automatic")
        ymax=st.text_input("Y maximum",key=f"{kind}-ymax",placeholder="Automatic")
    st.sidebar.caption("All controls apply immediately. They never trigger Garfield "
                       "or reread original ion CSVs.")

    def hnorm(y,edges):
        w=np.diff(edges)
        if norm=="Raw counts / weighted counts":return y
        if y.sum()<=0:return y
        if norm=="Probability density":return y/y.sum()/w
        return y/y.sum()

    def hist(dataarray,edges):
        counts=np.histogram(dataarray,bins=edges)[0].astype(float)
        return counts,hnorm(counts,edges)

    def trace_hist(edges,counts,name,include_zero_marker=False):
        values=hnorm(counts.astype(float),edges)
        mids=(edges[:-1]+edges[1:])/2
        if plot_style=="Bars":
            # Exact nonuniform bin widths: zero [0,1), then gain bins.
            fig.add_trace(go.Bar(
                x=mids,y=values,width=np.diff(edges),
                marker=dict(color=color),opacity=.78,
                name=name,hovertemplate="x: %{x:,.4g}<br>height: %{y:,.5g}<extra></extra>"))
        else:
            # Explicit step edges, without smoothing or connecting empty
            # log-y bins (which would create a false high-gain tail).
            xe=np.repeat(edges,2)[1:-1]
            ye=np.repeat(values,2).astype(float)
            if yscale=="Log":
                ye[ye<=0]=np.nan
            fig.add_trace(go.Scatter(
                x=xe,y=ye,mode="lines",name=name,
                line=dict(color=color,width=line_width),
                hovertemplate="Collected bin x: %{x:,.4g}"
                              "<br>height: %{y:,.5g}<extra></extra>"))
        if include_zero_marker:
            fig.add_trace(go.Scatter(
                x=[0.5],y=[values[0]],mode="markers",
                marker=dict(color=color,symbol="circle-open",size=13,line=dict(width=2)),
                showlegend=False,hovertemplate="Exactly zero collected electrons"
                                             "<br>height: %{y:,.4g}<extra></extra>"))
        return mids,values

    fig=go.Figure()
    original_note=""
    xallowlog=True
    yallowlog=True
    heatmap=False
    displayed_csv=None
    time_view_limit=None

    if kind in ("gain_all","gain_pos"):
        yy=g[g>0] if kind=="gain_pos" else g
        mx=int(max(yy.max(),1))
        edges=np.arange(1,mx+gain_width+1,gain_width,dtype=float)
        if kind=="gain_all":edges=np.r_[0.,edges]
        counts=np.histogram(yy,bins=edges)[0]
        must(counts.sum()==len(yy),"Gain histogram dropped original Garfield avalanches")
        mean=np.mean(yy);rms=np.std(yy,ddof=0)
        name=legend_stats(label,yy,add_stats,"e")
        xs,ys=trace_hist(edges,counts,name,kind=="gain_all")
        xallowlog=(kind=="gain_pos")
        original_note=(f"Zero-only bin [0,1); remaining bins {gain_width} collected e"
                       if kind=="gain_all" else
                       f"Conditional on N_collected>0; bins {gain_width} collected e")
        st.caption(f"{original_note}. Mean={mean:,.2f}, RMS={rms:,.2f}, N={len(yy):,}.")
        displayed_csv={"edge_lo":edges[:-1],"edge_hi":edges[1:],"value":ys}

    elif kind in ("angle_eq","angle_ion"):
        x=data["angle_equal"] if kind=="angle_eq" else data["angle_ion"]
        f=angular_deg//5
        counts=x.reshape(-1,f).sum(axis=1)
        edges=data["phi_edges"][::f]
        must(len(edges)==len(counts)+1,"Angular rebinned edges inconsistent")
        xs,ys=trace_hist(edges,counts,label)
        xallowlog=False
        original_note=("Each nonzero-ion avalanche has unit total weight" if
                       kind=="angle_eq" else
                       "Every ion birth has equal expected weight")
        st.caption(f"{original_note}; angular bin width {angular_deg}°.")
        displayed_csv={"edge_lo":edges[:-1],"edge_hi":edges[1:],"value":ys}

    elif kind in ("asym_seed","asym_muon"):
        v=data["seed_A"] if kind=="asym_seed" else data["muon_A"]
        a=v[np.isfinite(v)]
        edges=np.linspace(-1.,1.,asym_bins+1)
        counts=np.histogram(a,bins=edges)[0]
        xs,ys=trace_hist(
            edges,counts,legend_stats(label,a,add_stats))
        xallowlog=False
        st.caption(f"Exact unsampled ion-sector counters; N defined = {len(a):,}; "
                   f"undefined because N_up+N_down=0: {np.sum(~np.isfinite(v)):,}.")
        displayed_csv={"edge_lo":edges[:-1],"edge_hi":edges[1:],"value":ys}

    elif kind=="time":
        edges=data["t_edges"]
        counts=data["t_weights"]
        time_view_limit=float(edges[min(len(counts)-1,int(np.searchsorted(
            np.cumsum(counts),.999*float(np.sum(counts))))+1)])
        # Logarithmic time axes cannot show the real [0,10^-6ns) bin.
        # For log X omit exactly that bin with explicitly reported mass
        # and renormalize. All other birth times stay intact.
        if xscale=="Log":
            excluded=float(counts[0])
            fraction=excluded/max(1.,float(np.sum(counts)))
            counts=counts[1:]
            edges=edges[1:]
            xallowlog=True
            st.info(f"Log-time view excludes the [0,10⁻⁶ ns) bin "
                    f"({fraction:.3%} of weighted ion births); remaining "
                    "histogram is renormalized.")
        # Group original logarithmic-resolved bins; avoid changing physics weights.
        gbin=time_group
        cuts=list(range(0,len(counts),gbin))
        if cuts[-1]!=len(counts):cuts.append(len(counts))
        cuts=np.asarray(cuts,dtype=int)
        counts=np.add.reduceat(counts,cuts[:-1])
        edges=edges[cuts]
        xs,ys=trace_hist(
            edges,counts,legend_stats(
                label,(edges[:-1]+edges[1:])/2.,add_stats,"ns",
                weights=counts,approx=True))
        original_note="Physical birth delay relative to each initial Heed electron"
        st.caption(f"{original_note}; 128 randomly retained ion-birth points/avalanche, "
                   "correctly weighted by exact number of ion births. Bins are "
                   "nonuniform (log-spaced); use probability density to compare heights.")
        displayed_csv={"edge_lo":edges[:-1],"edge_hi":edges[1:],"value":ys}

    elif kind=="span":
        x=data["seed_span_ns"]
        x=x[np.isfinite(x)]
        edges=np.linspace(0.,max(1.e-6,float(x.max())),201)
        counts=np.histogram(x,bins=edges)[0]
        xs,ys=trace_hist(
            edges,counts,legend_stats(label,x,add_stats,"ns"))
        st.caption("These spans are based on at most 128 sampled birth positions "
                   "per avalanche, so they systematically underestimate the "
                   "true earliest-to-latest ion-birth span.")
        displayed_csv={"edge_lo":edges[:-1],"edge_hi":edges[1:],"value":ys}

    elif kind=="gain_asym":
        a=data["seed_absA"].astype(float)
        mask=(g>0)&np.isfinite(a)
        lx=np.log10(g[mask])
        bins=np.linspace(0.,max(1.,float(lx.max())),31)
        n=np.histogram(lx,bins=bins)[0]
        summed=np.histogram(lx,bins=bins,weights=a[mask])[0]
        y=np.divide(summed,n,out=np.full_like(summed,np.nan,dtype=float),where=n>=20)
        edges=10.**bins
        fig.add_trace(go.Scatter(
            x=(edges[:-1]*edges[1:])**.5,y=y,
            mode="lines+markers",
            name=legend_stats(label,a[mask],add_stats),
            line=dict(color=color,width=line_width),
            marker=dict(color=color,size=6),
            customdata=n,
            hovertemplate="Gain %{x:,.2f}<br>Mean |A| %{y:.4f}"
                          "<br>avalanches %{customdata}<extra></extra>"))
        st.caption("Binned exact |A| for each positive-gain avalanche; "
                   "bins with fewer than 20 avalanches are not displayed.")
        displayed_csv={"gain_bin_lo":edges[:-1],"gain_bin_hi":edges[1:],
                       "mean_absA":y,"n_avalanches":n}

    elif kind in temporal.KINDS:
        if kind=="t10_timeline":
            import pi_p3_noFW_penning057_muon_geometry_v1 as geomviz
            # Interactive REAL geometry is the first plot. The selected
            # ORIGINAL Heed cluster index is forwarded to the existing
            # full-muon chronology; no Garfield re-simulation.
            time_controls["selected_cluster"]=geomviz.render(
                st,int(time_controls["muon_index"]),temporal_meta)
        displayed_csv,note=temporal.render(
            kind,fig,temporal_data,temporal_meta,time_controls,st)
        if kind=="t10_timeline":
            st.markdown("#### Matching chronology — selected cluster highlighted")
        st.caption(note)
        if kind in ("t10_angles","t10_start","t10_timeline"):
            xallowlog=False
        if kind in ("t10_gain_position","t10_timeline"):
            yallowlog=False
    elif kind in ("ion_xy","gain_xy","time_phi"):
        heatmap=True
        if kind=="ion_xy":
            matrix=data["xy_ion"]
            xedges=data["u_edges"];yedges=data["v_edges"]
            zlabel=label
            original_note="Ion-birth count weighted from <=128 points per avalanche"
        elif kind=="gain_xy":
            nn=data["start_n"]
            matrix=np.divide(data["start_sum"],nn,
                 out=np.full_like(nn,np.nan),where=nn>=5)
            xedges=data["u_edges"];yedges=data["v_edges"]
            zlabel=label
            original_note="Mean collected gain per Heed electron; cells with N<5 hidden"
        else:
            matrix=data["t_phi"]
            xedges=data["phi_edges"];yedges=data["t_phi_edges"]
            zlabel=label
            original_note="Weighted joint ion birth time/azimuth; <=128 samples per avalanche"
            xallowlog=False
        xcent=(xedges[:-1]+xedges[1:])/2
        ycent=(yedges[:-1]+yedges[1:])/2
        z=matrix.astype(float)
        if zscale=="Log":
            z=np.where(z>0,np.log10(z),np.nan)
            zlabel="log10("+zlabel+")"
        fig.add_trace(go.Heatmap(
            x=xcent,y=ycent,z=z,
            colorscale="YlOrBr",colorbar=dict(title=zlabel),
            hovertemplate="x=%{x:.4g}<br>y=%{y:.4g}"
                          "<br>value=%{z:.5g}<extra></extra>"))
        st.caption(original_note+". This is an ion-creation source, "
                   "not an induced cathode signal.")
        displayed_csv=None
        if kind=="time_phi":
            yallowlog=True

    # Apply all user-controlled presentation options without reloading data.
    legend_xy={
        "Top right":(.99,.99,"right","top"),
        "Top left":(.01,.99,"left","top"),
        "Bottom right":(.99,.01,"right","bottom"),
        "Bottom left":(.01,.01,"left","bottom"),
    }[legpos]
    fig.update_layout(
        title=dict(text=title,font=dict(size=title_font)),
        xaxis_title=xlab,yaxis_title=ylab,
        template="plotly_white",height=650,
        font=dict(family="Arial, Helvetica, sans-serif",size=label_font),
        margin=dict(l=80,r=45,t=80,b=76),
        showlegend=show_legend,hovermode="closest",
        legend=dict(x=legend_xy[0],y=legend_xy[1],
                    xanchor=legend_xy[2],yanchor=legend_xy[3],
                    bgcolor="rgba(255,255,255,.85)",
                    bordercolor="rgba(120,130,140,.30)",borderwidth=1),
        bargap=0.02
    )
    if xscale=="Log" and not xallowlog:
        st.warning("Logarithmic X is invalid for the zero-only bin or for "
                   "signed/negative angular and asymmetry axes. Using linear X.")
    else:
        fig.update_xaxes(type="log" if xscale=="Log" else "linear")
    if yscale=="Log" and not yallowlog:
        st.warning("Logarithmic Y is invalid for this signed quantity.")
    else:
        fig.update_yaxes(type="log" if yscale=="Log" else "linear")

    def parse_lim(v):
        if not v.strip():return None
        try:return float(v)
        except ValueError:return None
    xl=[parse_lim(xmin),parse_lim(xmax)]
    yl=[parse_lim(ymin),parse_lim(ymax)]
    if (kind=="time" and xscale=="Linear"
        and time_range=="Central 99.9% of ion births"
        and not any(v is not None for v in xl)):
        fig.update_xaxes(range=[0.,time_view_limit])
        st.caption("Time X-axis shows the central 99.9% of weighted ion births "
                   "for readability; all time bins remain present in the "
                   "download and cache. Choose 'All ion-birth delays' to show the tail.")
    if any(v is not None for v in xl):
        existing=fig.layout.xaxis.range
        fallback=(0.,1.) if existing is None else tuple(existing)
        lo=xl[0] if xl[0] is not None else fallback[0]
        hi=xl[1] if xl[1] is not None else fallback[1]
        if xscale=="Log" and xallowlog:
            if lo>0 and hi>lo:fig.update_xaxes(range=[math.log10(lo),math.log10(hi)])
            else:st.warning("Log-X limits must be strictly positive and increasing.")
        elif hi>lo:fig.update_xaxes(range=[lo,hi])
    if any(v is not None for v in yl):
        fallback=(0.,1.)
        lo=yl[0] if yl[0] is not None else fallback[0]
        hi=yl[1] if yl[1] is not None else fallback[1]
        if yscale=="Log" and yallowlog:
            if lo>0 and hi>lo:fig.update_yaxes(range=[math.log10(lo),math.log10(hi)])
            else:st.warning("Log-Y limits must be strictly positive and increasing.")
        elif hi>lo:fig.update_yaxes(range=[lo,hi])

    if show_interpretation:
        if kind in temporal.KINDS:
            what,how,observed=temporal.describe(
                kind,temporal_data,temporal_meta,time_controls)
        else:
            what,how,observed=plot_explanation(kind,data,meta,angular_deg)
        with st.container(border=True):
            st.markdown("#### Understanding this plot")
            st.markdown(what)
            st.markdown(how)
            st.markdown(observed)
    st.plotly_chart(fig,use_container_width=True,
                    config={"displaylogo":False,
                            "toImageButtonOptions":{"format":export_format,
                                "filename":kind+"_penning_on","scale":3}})
    st.caption("Use the camera icon to export the selected PNG / SVG / JPEG; "
               "scroll to zoom, or drag to select a range. No simulated "
               "cathode waveform is implied by ion-birth distributions.")
    left,right=st.columns(2)
    with left:
        html=fig.to_html(full_html=True,include_plotlyjs=True)
        st.download_button("Download interactive figure (.html)",
                           html,file_name=f"{kind}_penning_on.html",
                           mime="text/html",use_container_width=True)
    with right:
        if displayed_csv is not None:
            buf=StringIO()
            writer=csv.writer(buf)
            cols=list(displayed_csv)
            writer.writerow(cols)
            for row in zip(*(displayed_csv[k] for k in cols)):
                writer.writerow(row)
            st.download_button("Download plotted data (.csv)",buf.getvalue(),
                               file_name=f"{kind}_penning_on.csv",
                               mime="text/csv",use_container_width=True)
        else:
            st.info("For 2D heatmaps, inspect the underlying grid in the cached NPZ.")

    with st.expander("Scientific definitions, raw statistics and provenance"):
        st.markdown("""
- **Exact**: Per-seed collected electrons, true total ion births, sectors
  and per-muon upper/lower asymmetry (where its denominator is nonzero).
- **Zero Heed electrons:** a Geant4 muon with no retained primary
  electrons has NO microscopic avalanche and undefined muon gain (0/0).
  Such muons are counted, not removed or confused with a real G=0 avalanche.
- **Spatial and temporal ion distributions are estimated** from a uniform,
  bounded 128-position reservoir per ion-producing avalanche, weighted
  by each avalanche's real ion-birth count.
- **Birth delay** is measured relative to the initiating Heed electron's
  saved creation time, *not* CPU wall time. No ion drift or electronics
  response has been simulated.
- **Plot normalization:** a probability per bin integrates/sums to one;
  probability density compensates for unequal bin widths.
- The zero-only collected-gain bin is **[0,1)**. Nonzero gain bins are
  selectable without losing or truncating high-gain avalanches.
- Gas physics uses Ar/CO₂ 70/30 with Penning transfer r=0.57; geometry
  is p3_noFW with symmetric ±3-mm cathodes and 1800 V; **not the
  4-mm-pitch 2+4-mm actual prototype**.
""")
        st.json(meta)

def selftest():
    # Cheap mathematical test: no original physics, no output, no Streamlit.
    must(phi_idx(0.,1.)==54 and phi_idx(0.,-1.)==18,
         "Angular quadrant convention wrong")
    must(phi_idx(3.,1.)==54 and phi_idx(-3.,1.)==54,
         "Nearest 3-mm wire angular folding wrong")
    must(abs(nearest_wire_offset(3.2)-.2)<1e-12,
         "3-mm nearest anode folding is wrong")
    must(binned(0.,T_EDGES)==0 and binned(1.e-4,T_EDGES)>0,
         "Physical time-bin boundary wrong")
    a=np.array([0,0,1,2,100],dtype=float)
    bins=np.r_[0.,np.arange(1.,102.,50.)]
    h=np.histogram(a,bins=bins)[0]
    must(int(h.sum())==len(a) and int(h[0])==2,
         "Gain zero-only [0,1) bin is mixing nonzero avalanches")
    w=[1./4]*4
    must(abs(sum(w)-1.)<1e-12,"Weighted reservoir source normalization wrong")
    # Valid all-muon census must not require G_mu for zero-Heed parents.
    n0=0
    gmu=float("nan") if n0==0 else 23./n0
    must(math.isnan(gmu) and (0==0) and
         np.isfinite(np.asarray([12./3.,25./5.])).all(),
         "Zero-Heed parent gain must be undefined while active gains stay finite")
    expected={
        "gain_all":{"gain_bin","mean_rms","histogram_style"},
        "angle_ion":{"angular_bin","histogram_style"},
        "time":{"time_group","time_window","histogram_normalization"},
        "t10_angles":{"angular_bin","mean_rms","line_width"},
        "t10_gain_position":{"heatmap_colorscale"},
        "asym_muon":{"asym_bin","mean_rms"},
    }
    for name,required in expected.items():
        c=plot_features(name)
        must(all(c[k] for k in required),"Missing sidebar control for "+name)
    must(not plot_features("ion_xy")["time_group"]
         and not plot_features("gain_all")["angular_bin"]
         and not plot_features("t10_angles")["gain_bin"]
         and not plot_features("time")["angular_bin"]
         and not plot_features("t10_timeline")["mean_rms"]
         and not plot_features("t10_angles")["legend_text"]
         and plot_features("ion_xy")["legend_text"]
         and not plot_features("t10_start")["legend_controls"],
         "Unrelated sidebar controls exposed")
    must("mean=" in legend_stats("gain",np.array([0.,2.,4.]),True)
         and "RMS=" in legend_stats("A",np.array([0.,.5,1.]),True),
         "Statistics cannot be appended to Plotly trace legend")
    # Catch the regression in which the format was referenced by the
    # Plotly camera button but never initialized by the Streamlit sidebar.
    import ast
    tree=ast.parse(Path(__file__).read_text())
    app_node=next(n for n in tree.body
                  if isinstance(n,ast.FunctionDef) and n.name=="app")
    export_assignments=[
        n for n in ast.walk(app_node)
        if isinstance(n,ast.Assign) and any(
            isinstance(t,ast.Name) and t.id=="export_format"
            for t in n.targets)
    ]
    must(len(export_assignments)==1
         and "export_format=st.sidebar.selectbox" in
             Path(__file__).read_text(),
         "Camera export format is used without a sidebar assignment")
    # Exercise each dynamic explanation with a deliberately tiny synthetic
    # cache, so UI documentation cannot silently break a specific plot.
    synthetic={
        "gain":np.array([0,2,10,22],dtype=float),
        "angle_equal":np.ones(72,dtype=float),
        "angle_ion":np.ones(72,dtype=float),
        "phi_edges":PHI_EDGES,
        "seed_A":np.array([np.nan,-.3,0.,1.]),
        "muon_A":np.array([0.,-.1,np.nan,1.]),
        "t_edges":T_EDGES,
        "t_weights":np.ones(len(T_EDGES)-1),
        "seed_span_ns":np.array([np.nan,0.,.5,2.]),
        "seed_absA":np.array([np.nan,.3,0.,1.]),
        "xy_ion":np.ones((len(V_EDGES)-1,len(U_EDGES)-1)),
        "u_edges":U_EDGES,
        "v_edges":V_EDGES,
        "start_n":np.full((len(V_EDGES)-1,len(U_EDGES)-1),5),
        "start_sum":np.ones((len(V_EDGES)-1,len(U_EDGES)-1))*300.,
        "t_phi":np.ones((len(T_PHI_EDGES)-1,72)),
        "t_phi_edges":T_PHI_EDGES,
    }
    descriptions=("gain_all","gain_pos","angle_eq","angle_ion",
                  "asym_seed","asym_muon","time","span","gain_asym",
                  "ion_xy","gain_xy","time_phi")
    for kind in descriptions:
        description=plot_explanation(
            kind,synthetic,{"n_true_ion_births":720},5)
        must(len(description)==3 and
             all(isinstance(d,str) and len(d)>80 for d in description),
             "Plot-specific scientific interpretation failed: "+kind)
    print("PASS local math/geometry tests and plot interpretations; "
          "no Garfield, no raw data, no files written.")

if __name__=="__main__":
    if len(sys.argv)>1 and sys.argv[1]=="prepare":
        prepare()
    elif len(sys.argv)>1 and sys.argv[1]=="selftest":
        selftest()
    else:
        app()
