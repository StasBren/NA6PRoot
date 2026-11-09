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
            muon_heed.append(hp)
            muon_gain.append(int(r["collected_avalanche_electrons"])/hp if hp else np.nan)
            muon_births.append(int(r["ion_births"]))
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
    must(np.isfinite(muon_gain).all(),
         "A true Geant4 muon yielded zero initial electrons")
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

def app():
    import streamlit as st
    import plotly.graph_objects as go
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
    a,b,c,d=st.columns(4)
    a.metric("Real J/ψ muons","1,000")
    b.metric("Microscopic avalanches",f"{len(g):,}")
    c.metric("Mean collected gain",f"{np.mean(g):,.1f}")
    d.metric("Zero collected",f"{np.mean(g==0):.2%}")

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
    st.sidebar.header("Edit the plot")
    kind=st.sidebar.selectbox("Plot",list(kinds),format_func=lambda x:kinds[x])
    if st.sidebar.button("Reset this plot to defaults"):
        for k in list(st.session_state):
            if k.startswith(kind+"-"):
                del st.session_state[k]
        st.rerun()
    title0,x0,y0,legend0=defaults[kind]
    # Each view has its own independent saved editable controls.
    title=st.sidebar.text_input("Title",value=title0,key=f"{kind}-title")
    label=st.sidebar.text_input("Legend text / colourbar title",value=legend0,key=f"{kind}-legend")
    xlab=st.sidebar.text_input("X-axis title",value=x0,key=f"{kind}-xlab")
    ylab=st.sidebar.text_input("Y-axis title",value=y0,key=f"{kind}-ylab")
    st.sidebar.divider()
    c1,c2=st.sidebar.columns(2)
    xscale=c1.selectbox("X scale",["Linear","Log"],key=f"{kind}-xscale")
    yscale=c2.selectbox("Y scale",["Linear","Log"],key=f"{kind}-yscale")
    color=st.sidebar.color_picker("Curve colour","#E46B08",key=f"{kind}-color")
    show_legend=st.sidebar.checkbox("Show legend",True,key=f"{kind}-legend-show")
    add_stats=st.sidebar.checkbox("Append mean / RMS to legend",
                                  kind.startswith("gain"),key=f"{kind}-stats")
    st.sidebar.divider()
    gain_width=st.sidebar.select_slider("Gain histogram bin width [electrons]",
        options=[10,20,25,50,100,200,500],value=50,key=f"{kind}-gainbin")
    angular_deg=st.sidebar.select_slider("Angular bin width [degrees]",
        options=[5,10,15,20,30],value=5,key=f"{kind}-angbin")
    asym_bins=st.sidebar.slider("Asymmetry bins",20,120,50,10,key=f"{kind}-asymbin")
    time_group=st.sidebar.select_slider("Time-bin grouping",
        options=[1,2,3,4,6,8,12],value=4,key=f"{kind}-tgroup")
    norm=st.sidebar.selectbox("Histogram normalization",
        ["Probability per bin","Probability density","Raw counts / weighted counts"],
        index=1 if kind=="time" else 0,
        key=f"{kind}-normalization")
    zscale=st.sidebar.selectbox("Heatmap colour scale",["Linear","Log"],key=f"{kind}-zscale")
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
        # Use bars for correct nonuniform bin widths and honest density interpretation.
        fig.add_trace(go.Bar(
            x=mids,y=values,width=np.diff(edges),
            marker=dict(color=color),opacity=.75,
            name=name,hovertemplate="x: %{x:,.4g}<br>height: %{y:,.5g}<extra></extra>"))
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

    if kind in ("gain_all","gain_pos"):
        yy=g[g>0] if kind=="gain_pos" else g
        mx=int(max(yy.max(),1))
        edges=np.arange(1,mx+gain_width+1,gain_width,dtype=float)
        if kind=="gain_all":edges=np.r_[0.,edges]
        counts=np.histogram(yy,bins=edges)[0]
        must(counts.sum()==len(yy),"Gain histogram dropped original Garfield avalanches")
        mean=np.mean(yy);rms=np.std(yy,ddof=0)
        name=(f"{label} — mean={mean:,.2f}, RMS={rms:,.2f}"
              if add_stats else label)
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
        xs,ys=trace_hist(edges,counts,label)
        xallowlog=False
        st.caption(f"Exact unsampled ion-sector counters; N defined = {len(a):,}; "
                   f"undefined because N_up+N_down=0: {np.sum(~np.isfinite(v)):,}.")
        displayed_csv={"edge_lo":edges[:-1],"edge_hi":edges[1:],"value":ys}

    elif kind=="time":
        edges=data["t_edges"]
        counts=data["t_weights"]
        # Group original logarithmic-resolved bins; avoid changing physics weights.
        gbin=time_group
        cuts=list(range(0,len(counts),gbin))
        if cuts[-1]!=len(counts):cuts.append(len(counts))
        cuts=np.asarray(cuts,dtype=int)
        counts=np.add.reduceat(counts,cuts[:-1])
        edges=edges[cuts]
        xs,ys=trace_hist(edges,counts,label)
        xallowlog=False # zero-delay bin [0,1e-6) is physically real
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
        xs,ys=trace_hist(edges,counts,label)
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
            mode="lines+markers",name=label,
            line=dict(color=color,width=2.5),
            marker=dict(color=color,size=6),
            customdata=n,
            hovertemplate="Gain %{x:,.2f}<br>Mean |A| %{y:.4f}"
                          "<br>avalanches %{customdata}<extra></extra>"))
        st.caption("Binned exact |A| for each positive-gain avalanche; "
                   "bins with fewer than 20 avalanches are not displayed.")
        displayed_csv={"gain_bin_lo":edges[:-1],"gain_bin_hi":edges[1:],
                       "mean_absA":y,"n_avalanches":n}

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
    fig.update_layout(
        title=dict(text=title,font=dict(size=22)),
        xaxis_title=xlab,yaxis_title=ylab,
        template="plotly_white",height=650,
        font=dict(family="Arial, Helvetica, sans-serif",size=14),
        margin=dict(l=80,r=45,t=80,b=76),
        showlegend=show_legend,hovermode="closest",
        legend=dict(x=.99,y=.99,xanchor="right",yanchor="top",
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

    st.plotly_chart(fig,use_container_width=True,
                    config={"displaylogo":False,
                            "toImageButtonOptions":{"format":"png",
                                "filename":kind+"_penning_on","scale":3}})
    st.caption("Use the camera icon above the plot to download PNG; "
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
  and per-muon upper/lower asymmetry.
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

if __name__=="__main__":
    if len(sys.argv)>1 and sys.argv[1]=="prepare":
        prepare()
    else:
        app()
