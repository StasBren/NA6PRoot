#!/usr/bin/env python3
"""Read-only Penning-ON wire-response *data export* for independent review.

Requires the already-verified 1000-muon temporal dashboard cache; does
not rerun Geant4, TrackHeed, Garfield, or modify original physics files.

Creates ONE shareable ZIP with:
 - exact per-avalanche gain and ion counts, A, original start u/v/w,
   cluster membership, physical timing, sampled w/angular/radial statistics
 - a 4D joint ion-birth histogram v0 x collected_gain x phi x delta_w
   (both ion and equal-avalanche weight), plus phi/time and phi/radius maps
 - selected ORIGINAL 3D+time Garfield birth sample clouds (<=128 each)
 - ATLAS-inspired white-background figures, provenance and README
 - every 20-muon source chunk SHA verified; exact ion conservation.
No full unsampled ion cloud or independent anode-arrival measurement exists.

Usage:
  python this_script.py selftest
  python -u this_script.py export
  python this_script.py inspect
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import shutil
import tempfile
import zipfile
from collections import defaultdict
from pathlib import Path

import numpy as np

import pi_p3_noFW_penning057_timing_correlations_v1 as temporal

ROOT=temporal.ROOT
ROOT_OUT=ROOT/"wire_model_export_v1"
OUTPUT=ROOT_OUT/"NA60_PENNING_ON_1000_WIRE_MODEL_REVIEW_V1.zip"
RUN=ROOT/"native_chunks"
PITCH=3.0
R_WIRE=.015
V0_EDGES=np.linspace(-3.,3.,25)
GAIN_EDGES=np.asarray([0,1,10,25,50,100,250,500,1000,2000,
                       5000,10000,20000,40000,100000,200001],float)
PHI_EDGES=np.linspace(-180.,180.,37)  # 10-degree for exported 4D matrix
W_EDGES=np.asarray([-100,-25,-10,-5,-2,-1,-.5,-.25,-.1,-.05,
                    -.02,-.01,0,.01,.02,.05,.1,.25,.5,1,2,5,
                    10,25,100],float)
R_EDGES=np.asarray([0.,.015,.02,.03,.05,.1,.2,.4,.8,1.6,
                    3.2,6.4,12.,100.],float)
T_EDGES=np.asarray([0.,1.,2.,5.,10.,15.,20.,25.,30.,40.,50.,
                    65.,80.,100.,125.,160.,200.,300.,500.,1000.,
                    1.e6],float)
BINS=(len(V0_EDGES)-1,len(GAIN_EDGES)-1,
      len(PHI_EDGES)-1,len(W_EDGES)-1)
# Csv values are all mm, ns or counted electron/ion numbers as labelled.
COLUMNS=(
    "parent_muon","crossing_id","cluster_index","primary_index",
    "seed_u_mm","u0_nearest_wire_mm","v0_mm","w0_mm",
    "heed_creation_t_ns","collected_gain_e","exact_ion_births",
    "exact_ion_up","exact_ion_down","A_exact",
    "n_saved_births","sampled_number_wire_cells",
    "sampled_fraction_births_near_seed_wire",
    "phi_circular_mean_deg","phi_resultant_R",
    "phi_circular_dispersion_1_minus_R",
    "phi_n_occupied_30deg_sectors","phi_arc90_deg",
    "phi_largest_empty_gap_deg",
    "delta_w_mean_mm","delta_w_rms_mm",
    "delta_w_q05_mm","delta_w_q50_mm","delta_w_q95_mm",
    "delta_w_q95_minus_q05_mm",
    "radius_mean_mm","radius_rms_mm",
    "radius_q50_mm","radius_q90_mm","radius_q99_mm",
    "birth_delay_q10_ns","birth_delay_q50_ns","birth_delay_q90_ns",
    "birth_delay_rms_ns","birth_delay_q90_minus_q10_ns",
    "ion_births_on_opposite_side_fraction_exact"
)

def need(cond,message):
    if not cond:raise RuntimeError(message)
def checksum(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for buf in iter(lambda:f.read(1<<20),b""):
            h.update(buf)
    return h.hexdigest()
def csvreader(path,required):
    with Path(path).open(newline="") as f:
        rd=csv.DictReader(f)
        need(rd.fieldnames is not None and
             set(required).issubset(rd.fieldnames) and
             len(set(rd.fieldnames))==len(rd.fieldnames),
             "Unexpected original native CSV schema: "+str(path))
        yield from rd
def key(row):
    return row["crossing_id"],int(row["primary_index"])
def cpp_round(x):
    return math.floor(x+.5) if x>=0 else math.ceil(x-.5)
def wire_u(u):
    return PITCH*cpp_round(u/PITCH)
def local_u(u):
    return u-wire_u(u)
def index(values,edges):
    values=np.asarray(values,dtype=np.float64)
    ids=np.searchsorted(edges,values,side="right")-1
    ids=np.where(values==edges[-1],len(edges)-2,ids)
    return ids
def moments(x):
    a=np.asarray(x,dtype=np.float64)
    if not a.size:return [float("nan")]*2
    return [float(a.mean()),float(a.std(ddof=0))]
def minimal_arc90(phi):
    """Minimum contiguous circular interval containing >=90% of draws.

    Distinguishes one angular lobe from widely separated lobes.
    Finite reservoir sampling; NOT a guarantee of a continuous shell.
    """
    theta=np.sort(np.mod(np.asarray(phi,dtype=float),360.))
    n=len(theta)
    if n<2:return float("nan")
    k=max(2,int(math.ceil(.9*n)))
    extended=np.r_[theta,theta+360.]
    return float(np.min(extended[k-1:k-1+n]-extended[:n]))
def angular_summary(phi):
    a=np.asarray(phi,dtype=np.float64)
    if not a.size:return [float("nan")]*4+[0]
    theta=np.deg2rad(a)
    z=np.mean(np.exp(1j*theta))
    R=float(np.abs(z))
    avg=float(np.rad2deg(np.angle(z))) if R>1.e-9 else float("nan")
    theta0=np.sort(np.mod(a,360.))
    largest=float(np.diff(np.r_[theta0,theta0[0]+360.]).max()) if len(a)>1 else float("nan")
    occupied=int(np.unique(np.floor(theta0/30.).astype(int)).size)
    return [avg,R,1.-R,minimal_arc90(a),largest,occupied]

def write_figures(dst,summary,records,joint):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    def style(ax):
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.tick_params(direction="in",top=False,right=False)
        ax.grid(False)
    def save(fig,name):
        fig.tight_layout()
        fig.savefig(dst/name,dpi=170,facecolor="white")
        plt.close(fig)

    phi_c=(PHI_EDGES[1:]+PHI_EDGES[:-1])/2
    v_c=(V0_EDGES[1:]+V0_EDGES[:-1])/2
    ion=np.sum(joint["births_ion_weighted_v0_gain_phi_dw"],axis=(1,3))
    row=ion.sum(axis=1)
    cond=np.divide(ion,row[:,None],
                    out=np.full_like(ion,np.nan),where=row[:,None]>0)
    fig,ax=plt.subplots(figsize=(10.5,5.4))
    heat=ax.pcolormesh(PHI_EDGES,V0_EDGES,cond,
                       cmap="YlOrBr",shading="flat")
    fig.colorbar(heat,ax=ax,label="Ion-weighted P(azimuth bin | starting v0 bin)")
    ax.set_xlabel("Ion birth azimuth around nearest anode [deg]")
    ax.set_ylabel("Initial TrackHeed electron v0 [mm]")
    ax.set_title("Penning ON · ion-birth azimuth conditioned on v0")
    style(ax)
    save(fig,"01_ion_azimuth_given_v0.png")

    v=np.asarray([r["v0_mm"] for r in records],float)
    A=np.asarray([r["A_exact"] for r in records],float)
    n=np.array([r["exact_ion_births"] for r in records])
    fine=np.isfinite(A)&(n>0)
    bc=np.bincount(np.maximum(0,np.minimum(len(v_c)-1,
        index(v[fine],V0_EDGES))),minlength=len(v_c))
    s=np.bincount(np.maximum(0,np.minimum(len(v_c)-1,
        index(v[fine],V0_EDGES))),weights=A[fine],minlength=len(v_c))
    avg=np.divide(s,bc,out=np.full_like(s,np.nan),where=bc>=40)
    fig,ax=plt.subplots(figsize=(9.5,4.9))
    ax.plot(v_c,avg,color="#E46B08",marker="o",lw=1.8,ms=4)
    ax.axhline(0,lw=1,color="#555555",linestyle=":")
    ax.axvline(0,lw=1,color="#555555",linestyle=":")
    ax.set_xlabel("Original Heed electron v0 [mm]")
    ax.set_ylabel("Mean exact ion-birth asymmetry A")
    ax.set_title("Penning ON · signed A versus initial v0")
    style(ax);save(fig,"02_exact_mean_asymmetry_vs_v0.png")

    w_rms=np.asarray([r["delta_w_rms_mm"] for r in records],float)
    count=np.asarray([r["n_saved_births"] for r in records],int)
    fig,ax=plt.subplots(figsize=(9.5,4.9))
    cm=(count==128)&np.isfinite(w_rms)
    v_bins=index(v[cm],V0_EDGES)
    med=np.full(len(v_c),np.nan)
    q20=np.full(len(v_c),np.nan)
    q80=np.full(len(v_c),np.nan)
    for b in range(len(v_c)):
        vals=w_rms[cm][v_bins==b]
        if len(vals)>=40:
            med[b]=np.median(vals)
            q20[b],q80[b]=np.quantile(vals,[.2,.8])
    ax.plot(v_c,med,color="#E46B08",lw=2,label="Median")
    ax.fill_between(v_c,q20,q80,color="#E46B08",alpha=.16,
                    label="20–80% avalanche spread")
    ax.legend(frameon=False)
    ax.set_xlabel("Original Heed electron v0 [mm]")
    ax.set_ylabel("Sampled per-avalanche RMS birth spread along wire [mm]")
    ax.set_title("Penning ON · longitudinal ion-birth extent vs v0")
    style(ax);save(fig,"03_wire_longitudinal_spread_vs_v0.png")

    arc=np.asarray([r["phi_arc90_deg"] for r in records],float)
    gain=np.asarray([r["collected_gain_e"] for r in records],float)
    m=(count==128)&np.isfinite(arc)&(gain>0)
    ab=np.logspace(0,max(1.,np.log10(gain[m].max())),18)
    ac=(ab[1:]*ab[:-1])**.5
    median=np.full(len(ac),np.nan)
    for i in range(len(ac)):
        vals=arc[m][(gain[m]>=ab[i])&(gain[m]<ab[i+1])]
        if len(vals)>=40:median[i]=np.median(vals)
    fig,ax=plt.subplots(figsize=(9.5,4.9))
    ax.semilogx(ac,median,color="#E46B08",marker="o",lw=2,ms=4)
    ax.set_xlabel("Collected gain [electrons per original Heed primary]")
    ax.set_ylabel("Median minimum angular arc containing 90% ion births [deg]")
    ax.set_title("Penning ON · sampled azimuthal extent vs collected gain")
    style(ax);save(fig,"04_angular_extent_vs_gain.png")

    wmid=(W_EDGES[:-1]+W_EDGES[1:])/2.
    for ii,(lo,hi,label) in enumerate([
        (-3.,-.5,"below"),(.5,3.,"above")
    ]):
        vsel=(v_c>=lo)&(v_c<hi)
        selection=np.sum(joint["births_ion_weighted_v0_gain_phi_dw"][vsel],
                         axis=(0,1))  # [phi,delta_w]
        selection=np.where(selection>0,selection,np.nan)
        selection=selection/np.nansum(selection)
        fig,ax=plt.subplots(figsize=(10.5,5.4))
        pic=ax.pcolormesh(W_EDGES,PHI_EDGES,selection,
                          cmap="YlOrBr",shading="flat")
        fig.colorbar(pic,ax=ax,label="Ion-weighted joint probability per bin")
        ax.set_xlim(-1.,1.)
        ax.set_xlabel("Ion-birth displacement along wire Δw [mm]")
        ax.set_ylabel("Ion-birth azimuth around nearest wire [degrees]")
        ax.set_title(f"Penning ON · P(φ,Δw) for initial v0 {label} anode")
        style(ax);save(fig,f"0{5+ii}_phi_delta_w_{label}_anode.png")

def export():
    # The existing temporal cache checks original 1000-muon starting
    # positions, Heed-cluster identities and exact per-electron gain.
    saved,meta=temporal.load()
    need(meta["n_independent_true_Geant4_muons"]==1000,
         "Not an original complete 1000-parent Penning-ON campaign")
    need(not ROOT_OUT.exists(),
         "Refusing to overwrite an earlier physical result. "
         "Keep original ZIP and choose another exporter version: "+str(ROOT_OUT))
    input_path=ROOT/"input/native_muon_crossings_input.csv"
    manifest=ROOT/"input/input_manifest.json"
    need(input_path.is_file() and manifest.is_file(),
         "Missing original 1000 real Geant4 crossings")
    source=json.loads(manifest.read_text())
    need(source["n_parents"]==1000 and
         checksum(input_path)==source["selected_input_sha256"],
         "Geant4 original muon identity / source file SHA changed")
    # No directory is created until the schema/input validation passes.
    arrays={
        "births_ion_weighted_v0_gain_phi_dw":np.zeros(BINS,dtype=np.float64),
        "births_equal_avalanche_v0_gain_phi_dw":np.zeros(BINS,dtype=np.float64),
        "births_ion_weighted_v0_phi_creation_delay":np.zeros(
            (len(V0_EDGES)-1,len(PHI_EDGES)-1,len(T_EDGES)-1),dtype=np.float64),
        "births_ion_weighted_v0_phi_radius":np.zeros(
            (len(V0_EDGES)-1,len(PHI_EDGES)-1,len(R_EDGES)-1),dtype=np.float64),
    }
    rows_all=[]
    exemplars=[]
    source_info=[]
    raw_total=0
    samples_total=0
    ion_creating=0
    full128=0
    selected_clouds=0
    for chunk in range(50):
        folder=RUN/f"chunk_{chunk:04d}"
        rec=json.loads((folder/"complete.json").read_text())
        need(rec["chunk"]==chunk and rec["n_parents"]==20
             and rec["stats"]["n_capped_seeds"]==0 and
             rec["stats"]["same_Heed_electron_positions_and_Garfield_seeds"] is True,
             "Original Penning ON 20-muon checkpoint failed physics gates")
        files=("muon_reference_seed.csv",
               "muon_reference_heed_primaries.csv",
               "muon_reference_ion_birth_sample.csv",
               "muon_reference_crossing.csv",
               "physics.log")
        for filename in files:
            need(rec["file_sha256"].get(filename)==checksum(folder/filename),
                 "Source-file SHA altered; refusing export: "+str(folder/filename))
        need("PHYSICS_PENNING_AR_CO2_70_30_R=0.57"
             in (folder/"physics.log").read_text(errors="replace"),
             "Original Garfield physics log does not confirm Penning ON")
        source_info.append({
            "chunk":chunk,"complete_sha256":checksum(folder/"complete.json"),
            "file_sha256":{name:rec["file_sha256"][name] for name in files}
        })
        primaries={}
        for h in csvreader(folder/"muon_reference_heed_primaries.csv",
                           ("crossing_id","primary_index","cluster_index",
                            "t_ns")):
            k=key(h)
            need(k not in primaries,"Duplicate initial TrackHeed electron")
            primaries[k]=(int(h["cluster_index"]),float(h["t_ns"]))
        seeds=list(csvreader(folder/"muon_reference_seed.csv",
            ("crossing_id","primary_index","seed_u_mm","seed_v_mm",
             "seed_w_mm","collected_e","ion_births","ion_up","ion_down",
             "capped")))
        need(len(seeds)==len(primaries),
             "Each original Garfield avalanche must map to exactly one Heed electron")
        loc={}
        for j,r in enumerate(seeds):
            k=key(r)
            need(k not in loc and k in primaries,
                 "Original native Garfield seed has no exact Heed identity")
            loc[k]=j
        sample_by_key=defaultdict(list)
        weighted_by_key=defaultdict(float)
        for ion in csvreader(folder/"muon_reference_ion_birth_sample.csv",
            ("crossing_id","primary_index","u_mm","v_mm","w_mm",
             "t_ns","seed_ion_births_total","reservoir_weight")):
            k=key(ion)
            need(k in loc,"Orphan ion birth from unknown Garfield primary")
            n=int(seeds[loc[k]]["ion_births"])
            m=min(128,n)
            w=float(ion["reservoir_weight"])
            need(n>0 and int(ion["seed_ion_births_total"])==n and
                 abs(w-n/m)<1.e-7*max(1.,n/m),
                 "Original microscopic ion-birth sample weight invalid")
            sample_by_key[k].append((
                float(ion["u_mm"]),float(ion["v_mm"]),
                float(ion["w_mm"]),float(ion["t_ns"])
            ))
            weighted_by_key[k]+=w
        per_chunk=[]
        sampled_chunk=0
        ion_chunk=0
        for j,r in enumerate(seeds):
            k=key(r)
            mu=chunk*20+next((p for p in range(20)
                if meta["muons"][chunk*20+p]["id"]==k[0]),-1000)
            need(chunk*20<=mu<chunk*20+20,"Garfield avalanche parent identity incorrect")
            cluster,th=primaries[k]
            gain=int(r["collected_e"])
            n=int(r["ion_births"])
            up,down=int(r["ion_up"]),int(r["ion_down"])
            need(min(gain,n,up,down)>=0 and up+down<=n and
                 int(r["capped"])==0,
                 "Invalid/capped native avalanche ion/electron census")
            srcu=float(r["seed_u_mm"])
            u0=local_u(srcu)
            v0=float(r["seed_v_mm"])
            w0=float(r["seed_w_mm"])
            A=(up-down)/(up+down) if up+down else float("nan")
            vv=sample_by_key.get(k,[])
            m=min(128,n)
            need(len(vv)==m and
                 abs(weighted_by_key.get(k,0.)-n)<1.e-6*max(1.,n),
                 "Incomplete original bounded Garfield ion birth reservoir")
            # Compare untouched exact cached collected counts/primary origins.
            gi=len(rows_all)
            need(gi<len(saved["gain"]) and
                 mu==int(saved["parent"][gi]) and
                 int(r["primary_index"])==int(saved["primary"][gi]) and
                 gain==int(saved["gain"][gi]) and
                 n==int(saved["ions"][gi]) and
                 abs(u0-float(saved["u0"][gi]))<1.e-5 and
                 abs(v0-float(saved["v0"][gi]))<1.e-5,
                 "Exported source no longer matches existing verified timing cache")
            if m:
                p=np.asarray(vv,dtype=np.float64)
                need(np.isfinite(p).all(),"Invalid Garfield physical 3D/time ion coordinates")
                near=np.array([wire_u(z) for z in p[:,0]])
                du=p[:,0]-near
                phi=np.rad2deg(np.arctan2(p[:,1],du))
                # Force atan2's +180 upper-bound into leftmost periodic bin.
                phi=np.where(phi>=180.,phi-360.,phi)
                dw=p[:,2]-w0
                rad=np.hypot(du,p[:,1])
                delay=p[:,3]-th
                need(np.min(delay)>=-1.e-6,
                     "Original ion birth time precedes initiating Heed creation")
                delay=np.maximum(delay,0.)
                ang=angular_summary(phi)
                wmean,wrms=moments(dw)
                rmean,rrms=moments(rad)
                tmean,trms=moments(delay)
                qdw=np.quantile(dw,[.05,.5,.95])
                qr=np.quantile(rad,[.5,.9,.99])
                qt=np.quantile(delay,[.1,.5,.9])
                ns=int(np.unique(np.round(near/PITCH).astype(int)).size)
                fraction_start=float(np.mean(np.isclose(near,wire_u(srcu),atol=1.e-7)))
                n_v=index([v0],V0_EDGES)[0]
                n_g=index([gain],GAIN_EDGES)[0]
                n_phi=index(phi,PHI_EDGES)
                n_w=index(dw,W_EDGES)
                n_t=index(delay,T_EDGES)
                n_r=index(rad,R_EDGES)
                need(0<=n_v<len(V0_EDGES)-1 and
                     0<=n_g<len(GAIN_EDGES)-1 and
                     np.all((n_phi>=0)&(n_phi<len(PHI_EDGES)-1)) and
                     np.all((n_w>=0)&(n_w<len(W_EDGES)-1)) and
                     np.all((n_t>=0)&(n_t<len(T_EDGES)-1)) and
                     np.all((n_r>=0)&(n_r<len(R_EDGES)-1)),
                     "A physical birth lies outside an exported histogram axis; "
                     "refusing to silently drop ion births")
                weight=n/m
                np.add.at(arrays["births_ion_weighted_v0_gain_phi_dw"],
                          (n_v,n_g,n_phi,n_w),weight)
                np.add.at(arrays["births_equal_avalanche_v0_gain_phi_dw"],
                          (n_v,n_g,n_phi,n_w),1./m)
                np.add.at(arrays["births_ion_weighted_v0_phi_creation_delay"],
                          (n_v,n_phi,n_t),weight)
                np.add.at(arrays["births_ion_weighted_v0_phi_radius"],
                          (n_v,n_phi,n_r),weight)
                ion_creating+=1
                full128+=(m==128)
                ion_chunk+=n
                sampled_chunk+=m
                metrics=dict(
                    phi_circular_mean_deg=ang[0],
                    phi_resultant_R=ang[1],
                    phi_circular_dispersion_1_minus_R=ang[2],
                    phi_arc90_deg=ang[3],
                    phi_largest_empty_gap_deg=ang[4],
                    phi_n_occupied_30deg_sectors=ang[5],
                    delta_w_mean_mm=wmean,
                    delta_w_rms_mm=wrms,
                    delta_w_q05_mm=qdw[0],
                    delta_w_q50_mm=qdw[1],
                    delta_w_q95_mm=qdw[2],
                    delta_w_q95_minus_q05_mm=qdw[2]-qdw[0],
                    radius_mean_mm=rmean,
                    radius_rms_mm=rrms,
                    radius_q50_mm=qr[0],
                    radius_q90_mm=qr[1],
                    radius_q99_mm=qr[2],
                    birth_delay_q10_ns=qt[0],
                    birth_delay_q50_ns=qt[1],
                    birth_delay_q90_ns=qt[2],
                    birth_delay_rms_ns=trms,
                    birth_delay_q90_minus_q10_ns=qt[2]-qt[0],
                    sampled_number_wire_cells=ns,
                    sampled_fraction_births_near_seed_wire=fraction_start,
                )
            else:
                metrics={fld:float("nan") for fld in (
                    "phi_circular_mean_deg","phi_resultant_R",
                    "phi_circular_dispersion_1_minus_R",
                    "phi_arc90_deg","phi_largest_empty_gap_deg",
                    "delta_w_mean_mm","delta_w_rms_mm",
                    "delta_w_q05_mm","delta_w_q50_mm","delta_w_q95_mm",
                    "delta_w_q95_minus_q05_mm","radius_mean_mm",
                    "radius_rms_mm","radius_q50_mm","radius_q90_mm",
                    "radius_q99_mm","birth_delay_q10_ns",
                    "birth_delay_q50_ns","birth_delay_q90_ns",
                    "birth_delay_rms_ns","birth_delay_q90_minus_q10_ns",
                    "sampled_fraction_births_near_seed_wire")}
                metrics["phi_n_occupied_30deg_sectors"]=0
                metrics["sampled_number_wire_cells"]=0
            row=dict(parent_muon=mu,crossing_id=k[0],
                     cluster_index=cluster,primary_index=k[1],
                     seed_u_mm=srcu,u0_nearest_wire_mm=u0,
                     v0_mm=v0,w0_mm=w0,heed_creation_t_ns=th,
                     collected_gain_e=gain,exact_ion_births=n,
                     exact_ion_up=up,exact_ion_down=down,
                     A_exact=A,n_saved_births=m,
                     ion_births_on_opposite_side_fraction_exact=(
                         (1.-np.sign(v0)*A)/2. if abs(v0)>.025 and
                         np.isfinite(A) else float("nan")),
                     **metrics)
            need(set(row)==set(COLUMNS),
                 "Internal exporter CSV schema changed unexpectedly")
            per_chunk.append(row)
            rows_all.append(row)
        need(sum(int(r["ion_births"]) for r in seeds)==ion_chunk,
             "Weighted original ion birth total not conserved by sampled inputs")
        raw_total+=ion_chunk
        samples_total+=sampled_chunk
        # Select illustrative original avalanche clouds STRATIFIED BY
        # different physical features; never pretend 128 points is full.
        candidates=[r for r in per_chunk if r["n_saved_births"]==128]
        chosen=[]
        if candidates:
            def choose(feature,reverse=True):
                ordered=sorted(candidates,key=lambda r:r[feature],
                               reverse=reverse)
                for z in ordered:
                    k=(z["crossing_id"],int(z["primary_index"]))
                    if k not in chosen:
                        chosen.append(k)
                        break
            choose("collected_gain_e",True)
            choose("delta_w_rms_mm",True)
            choose("phi_arc90_deg",True)
            choose("phi_largest_empty_gap_deg",True)
            choose("A_exact",True)
            choose("A_exact",False)
            choose("v0_mm",True)
            choose("v0_mm",False)
            for k in chosen:
                z=next(row for row in per_chunk
                       if row["crossing_id"]==k[0] and row["primary_index"]==k[1])
                for u,v,w,t in sample_by_key[k]:
                    exemplars.append((
                        z["parent_muon"],k[0],z["primary_index"],
                        z["cluster_index"],z["u0_nearest_wire_mm"],
                        z["v0_mm"],z["collected_gain_e"],
                        z["exact_ion_births"],u,v,w,t,
                        local_u(u),w-z["w0_mm"],
                        t-z["heed_creation_t_ns"]
                    ))
                selected_clouds+=1
        print(f"PASS SHA/census #{chunk+1:02d}/50: "
              f"{len(seeds)} Garfield avalanches, "
              f"{ion_chunk:,} ions, {sampled_chunk:,} sampled births",
              flush=True)
    need(len(rows_all)==len(saved["gain"]) and
         raw_total==meta["n_true_ion_births"] and
         samples_total==meta["n_sampled_ion_birth_positions"],
         "Full source 1000-muon microscopic census disagrees with sealed cache")
    for k in arrays:
        expected=(float(raw_total) if k.startswith("births_ion_weighted")
                  else float(ion_creating))
        need(abs(arrays[k].sum()-expected)<1.e-6*max(1.,expected),
             "Four-dimensional ion-birth probability mass lost: "+k)
    out_summary={
        "schema":"NA60_PENNING_ON_1000_WIRE_MODEL_REVIEW_V1",
        "n_independent_Geant4_Jpsi_muons":1000,
        "n_retained_TrackHeed_electrons_and_Garfield_avalanches":len(rows_all),
        "n_ion_creating_avalanches":ion_creating,
        "n_avalanches_with_128_sampled_births":full128,
        "n_true_created_positive_ions_exact":raw_total,
        "n_retained_uniform_ion_birth_samples":samples_total,
        "n_exported_example_avalanches":selected_clouds,
        "n_exported_example_ion_birth_points":len(exemplars),
        "gas":"Ar/CO2 (70:30), Penning transfer r=0.57",
        "geometry":"p3_noFW: pitch 3 mm, 30um wire diameter, +/-3mm cathodes, 1800V",
        "quantities_note":"A and G are EXACT per starting electron; spatial/time moments are from <=128 uniform ion-birth samples per avalanche",
        "physical_limitations":[
            "No independently measured drift/arrival angle, time or local electron energy at anode",
            "Delta w measured from initial Heed electron w0 includes electron drift diffusion and multiplication",
            "Ion birth is not positive-ion drift, Shockley-Ramo current or cathode collected charge",
            "Angle around EACH locally nearest birth wire; multiwire seeds are flagged",
            "A near +/-1 can arise from small avalanches; use full-128 subsample for morphology",
            "Station-dependent real wire numbering is not inferred from local p3 periodic birth coordinates",
            "Original 1000 Jpsi muon selection is one SHA-random crossing per Geant4 parent, not a flux-unbiased sample",
        ],
        "source_temporal_cache_sha256":checksum(temporal.CACHE),
        "original_geant_source_sha256":checksum(input_path),
        "native_provenance_50_chunks":source_info,
    }
    positive=[r for r in rows_all if r["n_saved_births"]==128]
    for name,field in (("sampled_sigma_w_mm","delta_w_rms_mm"),
                       ("sampled_angular_arc90_deg","phi_arc90_deg"),
                       ("sampled_angular_largest_gap_deg","phi_largest_empty_gap_deg"),
                       ("sampled_radius_median_mm","radius_q50_mm")):
        a=np.asarray([z[field] for z in positive],float)
        out_summary[name+"_full128"]={
            "median":float(np.median(a)),
            "q10":float(np.quantile(a,.1)),
            "q90":float(np.quantile(a,.9))}
    A=np.asarray([r["A_exact"] for r in rows_all],float)
    v0=np.asarray([r["v0_mm"] for r in rows_all],float)
    good=np.isfinite(A)
    out_summary["corr_exact_A_vs_Heed_v0_Pearson"]={
        "value":float(np.corrcoef(A[good],v0[good])[0,1]),
        "n":int(good.sum()),
        "note":"Descriptive association; individual avalanches from a common muon are correlated"
    }
    out_summary["per_seed_mean_collected_gain_e"]=float(np.mean(
        [r["collected_gain_e"] for r in rows_all]))
    # Construct outputs in a private sibling directory, then publish once.
    # Preexisting experiment files and previous dashboard cache are untouched.
    # Build everything in an isolated temporary directory UNDER the
    # established physics campaign root. No final output directory
    # appears until every CSV/NPZ/PNG and its SHA has been produced.
    # A failed export can be retried without deleting physics or artifacts.
    temp=Path(tempfile.mkdtemp(prefix=".assembling_wire_model_v1_",dir=ROOT))
    try:
        with (temp/"per_avalanche_morphology.csv").open("x",newline="") as f:
            writer=csv.DictWriter(f,fieldnames=COLUMNS)
            writer.writeheader()
            for row in rows_all:
                writer.writerow({k:("" if isinstance(row[k],float) and
                    not np.isfinite(row[k]) else row[k]) for k in COLUMNS})
        example_headers=("parent_muon","crossing_id","primary_index",
                        "cluster_index","u0_mm","v0_mm","gain_e",
                        "ions_exact","birth_u_mm","birth_v_mm",
                        "birth_w_mm","birth_t_ns","birth_u_wire_offset_mm",
                        "birth_delta_w_mm","birth_delta_t_ns")
        with (temp/"representative_original_birth_clouds.csv").open(
                "x",newline="") as f:
            writer=csv.writer(f)
            writer.writerow(example_headers)
            writer.writerows(exemplars)
        np.savez_compressed(temp/"joint_ion_source_histograms.npz",
            **arrays,
            v0_edges_mm=V0_EDGES,gain_edges_e=GAIN_EDGES,
            birth_phi_edges_deg=PHI_EDGES,birth_delta_w_edges_mm=W_EDGES,
            birth_radius_edges_mm=R_EDGES,birth_delay_edges_ns=T_EDGES)
        (temp/"summary.json").write_text(json.dumps(
            out_summary,indent=2,allow_nan=False)+"\n")
        (temp/"README.txt").write_text(
            "PENNING ON p3_noFW REAL 1000 MUONS — CATHODE SOURCE MODEL DATA\n"
            "========================================================\n\n"
            "per_avalanche_morphology.csv: exact G and ion counts + per-seed\n"
            "signed A, start u/v/w, sample moments, 90% angular arc, largest\n"
            "empty arc, radial and longitudinal (w) spread, and timing.\n"
            "Undefined morphology on 0-ion avalanches has empty fields.\n\n"
            "joint_ion_source_histograms.npz: four-dimensional arrays\n"
            "[start v0, collected gain bin, phi, birth delta-w].\n"
            "births_ion_weighted_v0_gain_phi_dw sums to EXACT total ion\n"
            "births; each sample carries true ions/sample count.\n"
            "births_equal_avalanche_v0_gain_phi_dw sums to NUMBER OF\n"
            "AVALANCHES CREATING IONS (each contributing total one).\n"
            "There are also ion-weighted [v0,phi,time] and [v0,phi,r]\n"
            "joint arrays, all with matching *_edges_* keys.\n\n"
            "representative_original_birth_clouds.csv: a stratified set\n"
            "of ACTUAL sampled original ion-birth positions from the same\n"
            "immutable Penning-ON Garfield runs (<=128 points per seed).\n\n"
            "Important: A measures where ions were BORN, not a Shockley\n"
            "Ramo cathode charge or ion drift. w spread includes preceding\n"
            "electron drift diffusion; it is not isolated multiplication.\n"
            "Do not equate sampled angular spread with a continuous sheath\n"
            "or current reaching anode/cathode electrodes.\n\n"
            "Source provenance and exact run limitations: summary.json.\n"
            "Share only the generated ZIP, not the huge raw native folders.\n")
        write_figures(temp,out_summary,rows_all,arrays)
        hashes={p.name:checksum(p) for p in sorted(temp.iterdir())
                if p.is_file()}
        (temp/"file_sha256.json").write_text(
            json.dumps(hashes,indent=2,sort_keys=True)+"\n")
        need(not ROOT_OUT.exists(),"An earlier export directory already exists")
        ROOT_OUT.mkdir(parents=False,exist_ok=False)
        with zipfile.ZipFile(OUTPUT,"x",compression=zipfile.ZIP_DEFLATED,
                             compresslevel=6,allowZip64=True) as z:
            for p in sorted(temp.iterdir()):
                if p.is_file():
                    z.write(p,arcname=p.name)
        with zipfile.ZipFile(OUTPUT) as z:
            corrupt=z.testzip()
            need(corrupt is None,"Physical morphology ZIP integrity failed: "+str(corrupt))
        print("PASS exact 1000 Geant4 parents and 50 SHA-sealed Penning ON chunks")
        print("Garfield avalanches",len(rows_all),
              "| ion-producing",ion_creating,
              "| 128-point",full128)
        print("True ion births",raw_total,
              "| bounded ion coordinates",samples_total)
        print("Export ZIP",OUTPUT,"bytes",OUTPUT.stat().st_size)
        print("Upload this ZIP to ChatGPT for quantitative joint morphology analysis.")
    except Exception:
        if OUTPUT.exists():
            OUTPUT.unlink()
        if ROOT_OUT.exists() and not any(ROOT_OUT.iterdir()):
            ROOT_OUT.rmdir()
        raise
    finally:
        shutil.rmtree(temp,ignore_errors=True)
def inspect():
    need(OUTPUT.exists(),"No generated export archive: run 'export' first")
    with zipfile.ZipFile(OUTPUT) as z:
        need(z.testzip() is None,"ZIP failed CRC verification")
        print("ZIP files:",*z.namelist(),sep="\n - ")
        summary=json.loads(z.read("summary.json"))
        print(json.dumps({k:v for k,v in summary.items()
                          if k!="native_provenance_50_chunks"},
                         indent=2))
def selftest():
    p=np.array([-179.,-160.,-159.,-90.,0.,90.,160.,179.])
    need(minimal_arc90(p)<360 and minimal_arc90(p)>100,
         "Circular 90% arc algorithm broken")
    a=angular_summary([0.,0.,0.,0.])
    need(abs(a[1]-1.)<1e-12 and abs(a[2])<1e-12 and a[-1]==1,
         "Circular resultant / angular occupancy broken")
    need(cpp_round(.5)==1 and cpp_round(-.5)==-1 and
         abs(local_u(3.2)-.2)<1e-12,
         "Nearest-wire phase geometry incorrect")
    need(index([0.,1.,200000.],GAIN_EDGES).tolist()==[
         0,1,len(GAIN_EDGES)-2],
         "Export gain includes the zero-only bin incorrectly")
    need(index([180.,-180.,0.],PHI_EDGES).tolist()==[
         len(PHI_EDGES)-2,0,len(PHI_EDGES)//2],
         "Periodic phi histogram boundary invalid")
    need(np.all(index([0.,.015,.03],R_EDGES)>=0),
         "30um physical wire radius bins invalid")
    fake=np.array([[.0,.015,.0,10.],[.02,0.,.01,11.],
                   [-.02,0.,-.01,12.]])
    wmean,wrms=moments(fake[:,2])
    need(abs(wmean)<1e-12 and wrms>0,
         "Longitudinal sampled ion-cloud spread mathematics broken")
    print("PASS exporter: periodic 90% angular extent, exact gain-zero bins, "
          "30um wire radius, longitudinal moments and wire phase. "
          "No output and no Garfield runs.")
def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode",choices=("selftest","export","inspect"))
    opts=parser.parse_args()
    if opts.mode=="selftest":selftest()
    elif opts.mode=="export":export()
    else:inspect()
if __name__=="__main__":
    main()
