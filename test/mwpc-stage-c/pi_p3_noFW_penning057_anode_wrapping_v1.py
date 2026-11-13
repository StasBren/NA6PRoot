#!/usr/bin/env python3
"""Penning-ON MWPC: exact ion-birth A vs initial electron v0 and wire wrapping.

Read-only extension to the existing 1000-muon Plot Studio.
NO new Garfield runs and NO new full-campaign cache. Uses the already
SHA-checked temporal correlation NPZ for population statistics.
The single-avalanche 2D/3D view reads ONLY a selected original
SHA-checked 20-muon native chunk (<=128 saved ion-birth positions/seed).

'Wrapping' here means *ion-birth spatial extent around a wire*, not
electron arrival angle, positive-ion drift, or induced cathode charge.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np

import pi_p3_noFW_penning057_timing_correlations_v1 as temporal

ROOT=temporal.ROOT
N_BOOTSTRAP_PARENTS=1000
PITCH_MM=3.
WIRE_RADIUS_MM=.015
V_EDGES=np.linspace(-3.,3.,25)
U_START_EDGES=np.linspace(-1.5,1.5,25)
ABS_V_EDGES=np.linspace(0.,3.,13)
THETA=(np.asarray(temporal.PHI_EDGES[:-1])+np.asarray(temporal.PHI_EDGES[1:]))/2.

KINDS={
    "wrap_v0_A":"Signed exact ion asymmetry vs Heed starting v₀",
    "wrap_v0_opposite":"Ion births on opposite side vs |v₀|",
    "wrap_v0_angle":"Ion-birth azimuth vs starting v₀ (heatmap)",
    "wrap_u0_angle":"Ion-birth azimuth vs starting u₀ (heatmap)",
    "wrap_dispersion":"Avalanche azimuth spread vs collected gain",
    "wrap_single":"Individual avalanche: ion cloud around anode (2D + 3D)",
}
DEFAULTS={
    "wrap_v0_A":("Signed ion-birth asymmetry vs initial electron v₀",
                 "Original Heed electron v₀ [mm]","Mean signed A_ion,hemi",
                 "Each ion-producing avalanche · exact A"),
    "wrap_v0_opposite":("Fraction of ion births opposite initial electron side",
                        "Initial electron distance |v₀| [mm]",
                        "Mean opposite-hemisphere fraction",
                        "Exact ion sectors · each avalanche equal weight"),
    "wrap_v0_angle":("Ion-birth angle conditioned on initial electron v₀",
                     "Ion-birth azimuth φ around closest wire [degrees]",
                     "Initial Heed electron v₀ [mm]",
                     "Ion-weighted probability per angular bin"),
    "wrap_u0_angle":("Ion-birth angle conditioned on initial electron u₀",
                     "Ion-birth azimuth φ around closest wire [degrees]",
                     "Initial Heed electron u₀ [mm]",
                     "Ion-weighted probability per angular bin"),
    "wrap_dispersion":("Sampled angular spread vs avalanche collected gain",
                       "Collected gain per initial Heed electron [electrons]",
                       "Mean circular dispersion (1−R)",
                       "128 sampled ion births per avalanche"),
    "wrap_single":("Actual ion-birth positions around one anode",
                   "u − selected wire centre [mm]",
                   "v relative to anode wire plane [mm]",
                   "Sampled ion-birth times relative to Heed start [ns]"),
}

def must(ok,reason):
    if not ok:raise RuntimeError(reason)
def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for b in iter(lambda:f.read(1<<20),b""):
            h.update(b)
    return h.hexdigest()
def readcsv(path,fields):
    with Path(path).open(newline="") as f:
        r=csv.DictReader(f)
        must(r.fieldnames is not None and
             set(fields).issubset(r.fieldnames) and
             len(set(r.fieldnames))==len(r.fieldnames),
             "Native Garfield CSV schema mismatch: "+str(path))
        yield from r
def cpp_round(x):
    return math.floor(x+.5) if x>=0 else math.ceil(x-.5)

def exact_eligible(d):
    """Exact signed A is valid even with <128 saved positions."""
    return np.isfinite(d["A"])&np.isfinite(d["v0"])&(d["ions"]>0)

def grouped_parent_bootstrap(x,y,parent,edges,minimum=40,boot=350):
    """Confidence bands resample independent real Geant4 muon parents.

    One electron is NOT one independent parent: each real muon generates
    many correlated ionisation electrons/avalanches.
    """
    x=np.asarray(x,float);y=np.asarray(y,float);parent=np.asarray(parent,int)
    edges=np.asarray(edges,float)
    nparents=N_BOOTSTRAP_PARENTS
    valid=np.isfinite(x)&np.isfinite(y)&(parent>=0)&(parent<nparents)
    xi=np.searchsorted(edges,x[valid],side="right")-1
    xi=np.where(x[valid]==edges[-1],len(edges)-2,xi)
    p=parent[valid];v=y[valid]
    good=(xi>=0)&(xi<len(edges)-1)
    xi=xi[good];p=p[good];v=v[good]
    nbin=len(edges)-1
    cnt=np.zeros((nparents,nbin),dtype=np.float64)
    summ=np.zeros((nparents,nbin),dtype=np.float64)
    np.add.at(cnt,(p,xi),1.)
    np.add.at(summ,(p,xi),v)
    counts=cnt.sum(axis=0)
    sums=summ.sum(axis=0)
    mean=np.divide(sums,counts,out=np.full(nbin,np.nan),
                   where=counts>=minimum)
    rng=np.random.default_rng(20261110)
    # Each bootstrap resamples complete real-muon parents, not electrons.
    draws=rng.multinomial(nparents,
                          np.full(nparents,1./nparents),size=boot)
    bcnt=draws@cnt
    bsum=draws@summ
    with np.errstate(divide="ignore",invalid="ignore"):
        bm=np.where(bcnt>=minimum,bsum/bcnt,np.nan)
    lower=np.full(nbin,np.nan)
    upper=np.full(nbin,np.nan)
    for j in range(nbin):
        vals=bm[:,j][np.isfinite(bm[:,j])]
        if len(vals)>=int(boot*.8) and counts[j]>=minimum:
            lower[j],upper[j]=np.quantile(vals,[.025,.975])
    return (edges[:-1]+edges[1:])/2,mean,lower,upper,counts

def circular_dispersion(fractions):
    """1-R; zero for a narrow single angular direction.

    1-R does NOT prove continuous encirclement. Two opposite lobes can
    also produce R≈0. Therefore inspect the 2D ion cloud as well.
    """
    f=np.asarray(fractions,float)
    must(f.ndim==2 and f.shape[1]==72,"Need per-avalanche 5-degree angles")
    theta=np.deg2rad(THETA)
    z=f@np.exp(1j*theta)
    return np.clip(1.-np.abs(z),0.,1.)

def plot_binned(fig,x,y,p,edges,label,color,line_width):
    import plotly.graph_objects as go
    centres,means,low,hi,count=grouped_parent_bootstrap(x,y,p,edges)
    usable=np.isfinite(means)
    fig.add_trace(go.Scatter(
        x=centres[usable],y=means[usable],mode="lines+markers",
        name=label,line=dict(color=color,width=line_width),
        marker=dict(size=7),customdata=count[usable],
        hovertemplate="x=%{x:.3f}<br>Mean=%{y:.4f}"
                      "<br>N avalanches=%{customdata:.0f}<extra></extra>"))
    valid=np.isfinite(low)&np.isfinite(hi)
    if np.any(valid):
        xx=centres[valid]
        fig.add_trace(go.Scatter(
            x=np.r_[xx,xx[::-1]],y=np.r_[hi[valid],low[valid][::-1]],
            fill="toself",fillcolor="rgba(228,107,8,.15)",
            line=dict(width=0),hoverinfo="skip",
            showlegend=False))
    return {"bin_centre":centres,"mean":means,"ci95_low":low,
            "ci95_high":hi,"n_avalanches":count}

def options(st,kind,d,meta,allow_lut=True):
    """Only physically relevant controls for each specific view."""
    out={"min_entries":40,"ion_weighted":True,"reliable":True}
    if kind in {"wrap_v0_A","wrap_v0_opposite"}:
        st.sidebar.caption(f"Shaded band: 95% bootstrap CI from "
                           f"{N_BOOTSTRAP_PARENTS:,} independent Geant4 muon parents.")
        st.sidebar.caption("Every avalanche with defined exact A is included; "
                           "no 128-point sampling requirement.")
    if kind in {"wrap_v0_angle","wrap_u0_angle"}:
        out["ion_weighted"]=st.sidebar.checkbox(
            "Weight by number of ion births",True,
            key=f"{kind}-ion-weighting")
        out["reliable"]=st.sidebar.checkbox(
            "Require all 128 sampled ion births",True,
            key=f"{kind}-reliable")
        if kind=="wrap_u0_angle":
            other="v₀"
            lo,hi,step=-3.,3.,.25
        else:
            other="u₀"
            lo,hi,step=-1.5,1.5,.125
        out["other_range"]=st.sidebar.slider(
            f"Condition on initial {other} [mm]",
            min_value=lo,max_value=hi,value=(lo,hi),step=step,
            key=f"{kind}-other-range")
        out["compare_drift"]=(allow_lut and st.sidebar.checkbox(
            "Compare with electron drift LUT500",False,
            key=f"{kind}-drift-compare"))
        if out["compare_drift"]:
            out["drift_file"]=st.sidebar.text_input(
                "Original microscopic p3_noFW LUT500 CSV",
                value=str(Path(__file__).resolve().parents[2]/
                "test_runs/mwpc_stage_c/pi_pitch_fieldwire_500_p3_v1/"
                "pilot/p3p0_anode_only/garfield_electrons.csv"),
                key=f"{kind}-drift-path")
            st.sidebar.caption(
                "Both maps use the same initial-position bins and 5° angles. "
                "Arrival angles refer to the collecting wire, but birth "
                "angles to the nearest wire to the original electron. "
                "These are independent microscopic populations, not paired trajectories.")
    if kind=="wrap_dispersion":
        st.sidebar.caption("1−R measures circular angular dispersion. "
                           "It is not proof of continuous wrapping.")
    if kind=="wrap_single":
        out["muon_index"]=st.sidebar.selectbox(
            "Choose real Geant4 muon",list(range(1000)),
            index=temporal.typical_muon(meta),
            format_func=lambda i:(
                f"Muon {i:03d} · {meta['muons'][i]['n_Heed']} Heed e"),
            key="wrap_single-muon")
        same=np.flatnonzero(d["parent"]==out["muon_index"])
        if same.size:
            seq=np.argsort(-d["ions"][same],kind="stable")
            indices=d["primary"][same[seq]].astype(int).tolist()
            out["seed_index"]=st.sidebar.selectbox(
                "Individual avalanche / Heed electron",indices,
                format_func=lambda j:(
                    f"Seed {j} · cluster {int(d['cluster'][same[np.where(d['primary'][same]==j)[0][0]]])}"
                    f" · gain {int(d['gain'][same[np.where(d['primary'][same]==j)[0][0]]]):,}"
                    f" · ions {int(d['ions'][same[np.where(d['primary'][same]==j)[0][0]]]):,}"
                ),
                key=f"wrap_single-seed-{out['muon_index']}")
        else:
            out["seed_index"]=None
        out["zoom"]=st.sidebar.select_slider(
            "Cross-section half-width [mm]",
            options=[.05,.1,.25,.5,1.5,3.],
            value=.1,key="wrap_single-zoom")
        out["view_3d"]=st.sidebar.checkbox(
            "Also show 3D ion-birth positions",True,
            key="wrap_single-3d")
    return out

def one_seed(index,seed,meta):
    """Select ONE muon/avalanche, read-only, with SHA raw-file verification."""
    must(0<=index<1000 and meta["muons"][index]["n_Heed"]>0,
         "A real muon without Heed electrons has no individual avalanche")
    parent=meta["muons"][index]["id"]
    d=ROOT/"native_chunks"/f"chunk_{index//20:04d}"
    receipt=json.loads((d/"complete.json").read_text())
    must(receipt["chunk"]==index//20 and receipt["n_parents"]==20
         and receipt["stats"]["n_capped_seeds"]==0,
         "Selected chunk is not a completed validated Penning-on run")
    for name in ("muon_reference_seed.csv",
                 "muon_reference_ion_birth_sample.csv",
                 "muon_reference_heed_primaries.csv"):
        must(sha(d/name)==receipt["file_sha256"][name],
             "Raw ion-birth source failed original SHA "+str(d/name))
    target=None
    for r in readcsv(d/"muon_reference_seed.csv",
                ("crossing_id","primary_index","seed_u_mm","seed_v_mm",
                 "ion_births","collected_e","ion_up","ion_down")):
        if r["crossing_id"]==parent and int(r["primary_index"])==seed:
            must(target is None,"Duplicate microscopic avalanche seed")
            target=r
    must(target is not None,"Selected seed not among real muon Garfield entries")
    th=None;cluster=None
    for r in readcsv(d/"muon_reference_heed_primaries.csv",
                     ("crossing_id","primary_index","cluster_index","t_ns")):
        if r["crossing_id"]==parent and int(r["primary_index"])==seed:
            must(th is None,"Duplicate initiating Heed electron")
            th=float(r["t_ns"]);cluster=int(r["cluster_index"])
    must(th is not None,"Original microscopic electron lost Heed identity")
    samples=[]
    n=int(target["ion_births"])
    for r in readcsv(d/"muon_reference_ion_birth_sample.csv",
                     ("crossing_id","primary_index","u_mm","v_mm",
                      "w_mm","t_ns","seed_ion_births_total","reservoir_weight")):
        if r["crossing_id"]==parent and int(r["primary_index"])==seed:
            must(int(r["seed_ion_births_total"])==n,
                 "Ion reservoir does not match exact avalanches")
            samples.append([float(r["u_mm"]),float(r["v_mm"]),
                            float(r["w_mm"]),float(r["t_ns"])])
    must(len(samples)==min(128,n),
         "Individual avalanche has incomplete sampled ion birth coordinates")
    rawseed=float(target["seed_u_mm"])
    wire=cpp_round(rawseed/PITCH_MM)*PITCH_MM
    pts=np.asarray(samples,dtype=float).reshape((-1,4))
    if len(pts):pts[:,0]-=wire
    return {
        "seed":seed,"parent":index,"cluster":cluster,
        "start_u_wire_mm":rawseed-wire,"start_v_mm":float(target["seed_v_mm"]),
        "t_heed":th,"collected_e":int(target["collected_e"]),
        "ion_births":n,
        "A":(int(target["ion_up"])-int(target["ion_down"]))/
             (int(target["ion_up"])+int(target["ion_down"]))
             if int(target["ion_up"])+int(target["ion_down"]) else float("nan"),
        "wire_u_local_mm":wire,
        "points":pts,
    }

def three_d_view(observation,halfwidth,colour_mode="heed"):
    import plotly.graph_objects as go
    pts=observation["points"]
    zmin,zmax=float(pts[:,2].min()),float(pts[:,2].max())
    extension=max(.05,.2*(zmax-zmin))
    z0=zmin-extension;z1=zmax+extension
    phi=np.linspace(0.,2*np.pi,45)
    cylinder_u=np.tile(WIRE_RADIUS_MM*np.cos(phi),(2,1))
    cylinder_v=np.tile(WIRE_RADIUS_MM*np.sin(phi),(2,1))
    cylinder_w=np.tile(np.linspace(z0,z1,2)[:,None],(1,len(phi)))
    if colour_mode not in {"heed","first"}:
        raise ValueError("Unknown 3D time-colour mode")
    t0=(float(np.min(pts[:,3])) if colour_mode=="first"
        else observation["t_heed"])
    dt=pts[:,3]-t0
    label=("Time from first ion [ns]" if colour_mode=="first"
           else "Birth delay [ns]")
    fig=go.Figure()
    fig.add_trace(go.Surface(
        x=cylinder_u,y=cylinder_v,z=cylinder_w,
        surfacecolor=np.zeros_like(cylinder_u),
        colorscale=[[0,"#666666"],[1,"#666666"]],
        opacity=.65,showscale=False,name="Real 30-μm wire",
        hoverinfo="skip"))
    fig.add_trace(go.Scatter3d(
        x=pts[:,0],y=pts[:,1],z=pts[:,2],
        mode="markers",marker=dict(size=4,color=dt,
                                colorscale="YlOrRd",showscale=True,
                                colorbar=dict(title=label)),
        name="Saved Garfield ion births",
        hovertemplate="Δu=%{x:.4g} mm<br>v=%{y:.4g} mm"
                      "<br>w=%{z:.4g} mm<extra></extra>"))
    fig.update_layout(
        title="3D sampled ion-birth cloud around anode wire",
        template="plotly_white",height=690,
        scene=dict(xaxis=dict(title="Δu [mm]",range=[-halfwidth,halfwidth]),
                   yaxis=dict(title="v [mm]",range=[-halfwidth,halfwidth]),
                   zaxis=dict(title="w [mm]",range=[z0,z1]),
                   aspectmode="data"),
        margin=dict(l=0,r=0,t=70,b=0))
    return fig


def conditional_birth_angular_map(d,kind,controls):
    """Row-conditional birth-angle probabilities, with optional orthogonal slice.

    Efficient bounded chunks avoid a temporary N_avalanches x 72 float64
    array when all 30,000 muons are selected.
    """
    axis="u0" if kind=="wrap_u0_angle" else "v0"
    other="v0" if axis=="u0" else "u0"
    edges=U_START_EDGES if axis=="u0" else V_EDGES
    x=np.asarray(d[axis],float)
    y=np.asarray(d[other],float)
    ions=np.asarray(d["ions"])
    lo,hi=controls.get("other_range",(-3.,3.) if axis=="u0" else (-1.5,1.5))
    good=(ions>0)&np.isfinite(x)&np.isfinite(y)&(y>=lo)&(y<=hi)
    if controls["reliable"]:
        good &= np.asarray(d["sample_count"])==128
    inds=np.flatnonzero(good)
    nb=len(edges)-1
    bins=np.zeros((nb,72),dtype=np.float64)
    n=np.zeros(nb,dtype=np.int64)
    angle_bins=np.arange(72,dtype=np.int64)[None,:]
    for batch in np.array_split(inds,max(1,(len(inds)+8191)//8192)):
        if not len(batch):continue
        ix=np.searchsorted(edges,x[batch],side="right")-1
        ix=np.where(x[batch]==edges[-1],nb-1,ix)
        keep=(ix>=0)&(ix<nb)
        if not np.any(keep):continue
        idx=ix[keep]
        rows=batch[keep]
        np.add.at(n,idx,1)
        weights=(np.asarray(d["angle_fraction"][rows],dtype=np.float64)*
                 (ions[rows,None] if controls["ion_weighted"] else 1.))
        flat=(idx[:,None]*72+angle_bins).ravel()
        bins+=np.bincount(flat,weights=weights.ravel(),
                          minlength=nb*72).reshape(nb,72)
    denom=bins.sum(axis=1)
    probs=np.divide(bins,denom[:,None],
                    out=np.full_like(bins,np.nan),
                    where=denom[:,None]>0)
    return edges,probs,n,int(n.sum())


DRIFT_LUT500_SHA="8b00c43680c8fb63abb0f83e1da6600d48c83cce211bd205bca8c9688740cdbf"


def load_verified_drift_angles(path,mtime,size):
    """Read-only raw frozen drift LUT. Different trajectory from full avalanches."""
    import pandas as pd
    file=Path(path)
    if not file.is_file() or file.stat().st_mtime_ns!=mtime or file.stat().st_size!=size:
        raise ValueError("Frozen drift LUT was removed or changed")
    if sha(file)!=DRIFT_LUT500_SHA:
        raise ValueError("Frozen p3_noFW drift LUT500 SHA mismatch; refusing comparison")
    cfg=json.loads((file.parent/"frozen_config.json").read_text())
    expected={"pitch_mm":3.,"field_wire_diam_um":0.,
              "n_u":30,"n_v":30,"samples_per_box":500,
              "gap_minus_mm":3.,"gap_plus_mm":3.,
              "anode_voltage_v":1800.,"magnetic_field_t":0.}
    for key,val in expected.items():
        if cfg.get(key)!=val:
            raise ValueError(f"Drift LUT geometry mismatch: {key}={cfg.get(key)!r}")
    columns=["iu","iv","sample","u_start_mm","v_start_mm",
             "collected","phi_hit_deg","wire_index"]
    table=pd.read_csv(file,usecols=columns)
    if len(table)!=450000:
        raise ValueError("Frozen drift LUT500 must contain 450000 trials")
    iu=table["iu"].to_numpy(dtype=np.int64)
    iv=table["iv"].to_numpy(dtype=np.int64)
    trial=table["sample"].to_numpy(dtype=np.int64)
    if np.any((iu<0)|(iu>=30)|(iv<0)|(iv>=30)|
              (trial<0)|(trial>=500)):
        raise ValueError("Drift LUT cell/trial indices are invalid")
    idx=(iu*30+iv)*500+trial
    if len(np.unique(idx))!=450000:
        raise ValueError("Drift LUT cell/trial identities are incomplete or duplicated")
    u=table["u_start_mm"].to_numpy(dtype=float)
    v=table["v_start_mm"].to_numpy(dtype=float)
    if not np.all(np.isfinite(u)&np.isfinite(v)&
                  (u>=-1.5+.1*iu-1e-8)&(u<=-1.5+.1*(iu+1)+1e-8)&
                  (v>=-3.+.2*iv-1e-8)&(v<=-3.+.2*(iv+1)+1e-8)):
        raise ValueError("Drift LUT starts do not match frozen 30x30 grid")
    collected=table["collected"].to_numpy(dtype=np.int64)==1
    phi=table["phi_hit_deg"].to_numpy(dtype=float)
    if np.any(collected & (~np.isfinite(phi)|(phi < -180.)|(phi > 180.))):
        raise ValueError("Collected drift electron without physical arrival azimuth")
    return u[collected],v[collected],phi[collected],int(collected.sum())


def drift_angle_comparison(st,kind,controls,birth_probs,birth_counts):
    """Show original drift-arrival phi below Garfield ion-birth conditional phi."""
    import plotly.graph_objects as go
    file=Path(controls["drift_file"]).expanduser()
    if not file.is_file():
        st.warning("Original microscopic drift LUT not found: "+str(file))
        return None
    try:
        u,v,phi,nall=st.cache_data(show_spinner="Verifying and loading the original LUT500…")(
            load_verified_drift_angles)(str(file.resolve()),file.stat().st_mtime_ns,
                                        file.stat().st_size)
    except (OSError,ValueError,KeyError) as ex:
        st.error("Drift comparison blocked: "+str(ex))
        return None
    axis="u0" if kind=="wrap_u0_angle" else "v0"
    x=u if axis=="u0" else v
    other=v if axis=="u0" else u
    edges=U_START_EDGES if axis=="u0" else V_EDGES
    lo,hi=controls.get("other_range",(-3.,3.) if axis=="u0" else (-1.5,1.5))
    ok=(other>=lo)&(other<=hi)
    hist,_,_=np.histogram2d(x[ok],phi[ok],bins=(edges,temporal.PHI_EDGES))
    den=hist.sum(axis=1)
    probs=np.divide(hist,den[:,None],out=np.full_like(hist,np.nan),
                    where=den[:,None]>0)
    comparable=(birth_counts>=40)&(den>=40)&np.isfinite(birth_probs).all(axis=1)
    if np.any(comparable):
        tv=.5*np.abs(birth_probs[comparable]-probs[comparable]).sum(axis=1)
        st.caption(f"Start-coordinate bins with at least 40 full Garfield avalanches "
                   f"and 40 drift endpoints: {int(comparable.sum())}. "
                   f"Median conditional total-variation distance = {np.median(tv):.3f} "
                   "(0 identical; 1 disjoint). These are NOT paired histories "
                   "or a proof of exact drift→avalanche matching.")
    fig=go.Figure(go.Heatmap(
        x=THETA,y=(edges[:-1]+edges[1:])/2.,z=probs,colorscale="Blues",
        zmin=0.,zmax=None,
        colorbar=dict(title="P(arrival φ | start bin)"),
        hovertemplate="Electron arrival φ=%{x:.1f}°<br>"
                      f"initial {axis}=%{{y:.3f}} mm"
                      "<br>conditional probability=%{z:.4g}<extra></extra>"))
    fig.update_layout(
        title="Electron arrival azimuth · frozen drift-only LUT500",
        xaxis_title="Azimuth around collecting anode [degrees]",
        yaxis_title=f"Initial electron {axis} [mm]",
        template="plotly_white",height=640,
        margin=dict(l=70,r=110,t=80,b=70),
        font=dict(family="Arial,Helvetica,sans-serif",size=15))
    fig.update_xaxes(range=[-180.,180.])
    return fig

def render(kind,fig,d,meta,controls,st):
    """Return (histogram columns or None, annotation, optional second Figure)."""
    import plotly.graph_objects as go
    x=d["v0"].astype(float)
    A=d["A"].astype(float)
    p=d["parent"].astype(int)
    mask=exact_eligible(d)
    c=controls["color"];lw=controls["line_width"]
    if kind=="wrap_v0_A":
        out=plot_binned(fig,x[mask],A[mask],p[mask],V_EDGES,
                        controls["legend"],c,lw)
        fig.add_hline(y=0,line_dash="dot",line_color="#8B939E",line_width=1)
        fig.add_vline(x=0,line_dash="dot",line_color="#8B939E",line_width=1)
        return out,(
            f"Exact ion-sector asymmetry; {int(mask.sum()):,} valid individual "
            "avalanches (all ion multiplicities), each given equal weight. "
            f"Shaded interval resamples {N_BOOTSTRAP_PARENTS:,} complete "
            "muons, not electrons."),None
    if kind=="wrap_v0_opposite":
        good=mask&(np.abs(x)>=.025)
        # Exact fractions relative to UP/DOWN BIRTHS (boundary excluded).
        # v0 sign is start hemisphere, NOT actual anode-arrival side.
        other=(1.-np.sign(x[good])*A[good])/2.
        out=plot_binned(fig,np.abs(x[good]),other,p[good],
                        ABS_V_EDGES,controls["legend"],c,lw)
        fig.update_yaxes(range=[0.,1.])
        return out,(
            f"{int(good.sum()):,} avalanches with |v₀|≥0.025 mm. "
            "Opposite means the hemisphere opposite the INITIAL Heed "
            "electron position, not measured electron arrival side. "
            "Boundary births are excluded from exact A's denominator."),None
    if kind in {"wrap_v0_angle","wrap_u0_angle"}:
        edges,matrix,n,total=conditional_birth_angular_map(d,kind,controls)
        axis="u₀" if kind=="wrap_u0_angle" else "v₀"
        z=matrix
        title=controls["legend"]
        if controls.get("zscale")=="Log":
            with np.errstate(divide="ignore",invalid="ignore"):
                z=np.where(matrix>0,np.log10(matrix),np.nan)
            title="log10("+title+")"
        fig.add_trace(go.Heatmap(
            x=THETA,y=(edges[:-1]+edges[1:])/2.,
            z=z,colorscale="YlOrBr",
            colorbar=dict(title=title),
            hovertemplate=f"Ion birth φ=%{{x:.1f}}°<br>initial {axis}=%{{y:.3f}} mm"
                          "<br>conditional probability=%{z:.4g}<extra></extra>"))
        fig.update_xaxes(range=[-180.,180.])
        other="v₀" if kind=="wrap_u0_angle" else "u₀"
        lo,hi=controls.get("other_range",(np.nan,np.nan))
        aux=None
        if controls.get("compare_drift"):
            aux=drift_angle_comparison(st,kind,controls,matrix,n)
        return None,(
            f"{total:,} avalanche ion clouds. Each row normalised independently; "
            f"ion CREATION azimuth conditioned on initial Heed {axis}, "
            f"with {other} in [{lo:.3g},{hi:.3g}] mm. "
            +("Ion weighted." if controls["ion_weighted"] else
              "Each avalanche equally weighted.")
            +" Birth azimuth refers to the initial nearest wire; "
            "drift LUT arrival azimuth refers to its collecting wire. "
            "A different wire is a real but relatively rare exception."),aux
    if kind=="wrap_dispersion":
        good=(d["sample_count"]==128)&(d["ions"]>0)&(d["gain"]>0)
        spread=circular_dispersion(d["angle_fraction"][good])
        gg=d["gain"][good].astype(float)
        edges=np.geomspace(1.,max(2.,float(gg.max())),23)
        # Bootstrapped mean dispersion by independent parent.
        out=plot_binned(fig,gg,spread,p[good],edges,
                        controls["legend"],c,lw)
        return out,(
            f"{int(good.sum()):,} ion-producing avalanches with all 128 "
            "samples and gain>0. R=|sum p(φ)e^{iφ}|; 1−R=0 for "
            "a narrow angular direction and approaches 1 for balanced "
            "azimuthal populations. An opposite two-lobe distribution "
            "can also give 1−R≈1, so **this does not prove continuous "
            "encirclement**. Inspect individual 2D/3D birth clouds."),None
    if kind=="wrap_single":
        mi=controls["muon_index"];seed=controls["seed_index"]
        if seed is None:
            st.info("The selected real muon has no retained Heed electrons.")
            return None,"No microscopic avalanche available for this muon.",None
        @st.cache_data(show_spinner="Reading one sealed 20-muon Garfield chunk...")
        def cached_single(muon_idx,primary_idx,source_sha):
            return one_seed(muon_idx,primary_idx,meta)
        src=ROOT/"native_chunks"/f"chunk_{mi//20:04d}"/"complete.json"
        obs=cached_single(mi,seed,src.stat().st_mtime_ns)
        matches=(d["parent"]==mi)&(d["primary"]==seed)
        must(int(matches.sum())==1 and
             int(d["gain"][matches][0])==obs["collected_e"] and
             int(d["ions"][matches][0])==obs["ion_births"],
             "The selected microscopic Garfield avalanche no longer "
             "matches the verified temporal cache")
        pts=obs["points"];half=controls["zoom"]
        wirecol="#555555"
        fig.add_shape(type="circle",xref="x",yref="y",
                      x0=-WIRE_RADIUS_MM,y0=-WIRE_RADIUS_MM,
                      x1=WIRE_RADIUS_MM,y1=WIRE_RADIUS_MM,
                      line=dict(color=wirecol,width=1),
                      fillcolor="#A0A0A0")
        for j in (-1,1):
            if abs(j*PITCH_MM)<=half:
                fig.add_shape(type="circle",xref="x",yref="y",
                              x0=j*PITCH_MM-WIRE_RADIUS_MM,
                              x1=j*PITCH_MM+WIRE_RADIUS_MM,
                              y0=-WIRE_RADIUS_MM,y1=WIRE_RADIUS_MM,
                              line=dict(color=wirecol,width=1),
                              fillcolor="#A0A0A0")
        if len(pts):
            dt=pts[:,3]-obs["t_heed"]
            fig.add_trace(go.Scatter(
                x=pts[:,0],y=pts[:,1],mode="markers",
                marker=dict(size=8,color=dt,colorscale="YlOrRd",
                            showscale=True,
                            colorbar=dict(title="Birth delay [ns]"),
                            line=dict(width=.5,color="#444444")),
                name=controls["legend"],
                customdata=np.c_[pts[:,2],pts[:,3]],
                hovertemplate="Δu=%{x:.5f} mm<br>v=%{y:.5f} mm"
                              "<br>w=%{customdata[0]:.4f} mm"
                              "<br>absolute birth t=%{customdata[1]:.4f} ns"
                              "<extra></extra>"))
        fig.add_trace(go.Scatter(
            x=[obs["start_u_wire_mm"]],y=[obs["start_v_mm"]],
            mode="markers",marker=dict(symbol="star",size=12,color="#1676B3"),
            name="Original Heed electron start",
            hovertemplate="Starting electron Δu=%{x:.4g} mm"
                          "<br>v=%{y:.4g} mm<extra></extra>"))
        fig.update_xaxes(range=[-half,half],showgrid=False)
        fig.update_yaxes(range=[-half,half],showgrid=False,
                         scaleanchor="x",scaleratio=1)
        # The drawn anode circle is physically 30um in diameter.
        aux=three_d_view(obs,half) if len(pts)>0 and controls["view_3d"] else None
        if len(pts)>0:
            r=np.hypot(pts[:,0],pts[:,1])
            frac50=float(np.mean(r<.05))
            radial=(f"{frac50:.1%} of the sampled ion births lie "
                    "within 50 μm of this anode axis")
        else:
            radial="No ion births or sampled coordinates in this avalanche"
        return None,(
            f"Actual Geant4 muon #{mi}, TrackHeed cluster "
            f"#{obs['cluster']}, initial electron #{seed}; "
            f"**exact collected electrons={obs['collected_e']:,}**, "
            f"exact ion births={obs['ion_births']:,}, "
            f"exact A={obs['A']:+.3f}; "
            f"{len(pts)} uniformly retained ion-birth positions. "
            f"Wire centre u={obs['wire_u_local_mm']:.3f} mm in "
            "Garfield's local 3-mm wire lattice; wire radius **15 μm** "
            "(drawn to scale). "+radial+
            ". The blue star is the original Heed starting position, "
            "which may lie outside the chosen zoom. Neither panel "
            "simulates post-creation ion drift."),aux
    raise RuntimeError("No wire-wrapping plot renderer for "+kind)

def describe(kind,d,meta,settings):
    exact=exact_eligible(d)
    if kind=="wrap_v0_A":
        good=exact
        x=d["v0"][good];A=d["A"][good]
        neg=A[x<-.25];pos=A[x>.25]
        return (
            "**What you see:** Mean **signed** ion-birth hemisphere asymmetry "
            "A as a function of the original electron starting coordinate "
            "v₀ relative to the anode-wire plane. Negative and positive "
            "A mean predominantly below and above the anode.",
            "**How measured:** A is an exact full ion-birth callback ratio "
            "per electron avalanche. Every ion-producing avalanche has "
            "equal weight. The independent-muon bootstrap 95% interval "
            "accounts for many correlated electrons coming from one muon.",
            f"**This Penning-ON sample:** {len(A):,} defined avalanches. "
            f"Mean A for v₀<−0.25 mm = {neg.mean():+.3f}; "
            f"for v₀>+0.25 mm = {pos.mean():+.3f}. "
            "This tests initial-position predictability, not a measured "
            "electron arrival azimuth or actual cathode charge ratio."
        )
    if kind=="wrap_v0_opposite":
        good=exact&(np.abs(d["v0"])>=.025)
        f=(1.-np.sign(d["v0"][good])*d["A"][good])/2.
        return (
            "**What you see:** The fraction of counted ion births in the "
            "**opposite hemisphere from where the initial electron was born**, "
            "averaged across avalanches at a given |v₀|.",
            "**How measured:** Exact sector counters give "
            "f_opp=(1−sign(v₀)A)/2; the undefined v₀≈0 region is omitted. "
            "The denominator excludes sector-boundary counts. "
            "Parent-level bootstrap handles correlated electron starts.",
            f"**This Penning-ON sample:** {int(good.sum()):,} eligible "
            f"avalanches, mean opposite-side ion fraction {f.mean():.2%}. "
            "Even a small opposite-side fraction does not prove that "
            "ions create a continuous shell wrapping the wire."
        )
    if kind in {"wrap_v0_angle","wrap_u0_angle"}:
        axis="u₀" if kind=="wrap_u0_angle" else "v₀"
        return (
            "**What you see:** A conditional angular map. Each horizontal "
            f"row fixes the starting electron {axis}; horizontal position is "
            "where the avalanche's ions were born in azimuth around a wire.",
            "**How measured:** Up to 128 sampled ion-birth positions per "
            "avalanche become a 5° histogram. In ion-weighted mode, each "
            "sample represents exact N_ions/N_samples births; each v₀ "
            "row is separately normalised. The angle is not the incoming "
            "electron's drift/arrival angle.",
            "**Interpretation:** A bright angular band over one hemisphere "
            "means preferential birth there; bands spanning both sides "
            "suggest multiplication around a wider portion of the wire. "
            "The joint distribution is still not a cathode-current model."
        )
    if kind=="wrap_dispersion":
        return (
            "**What you see:** Mean sampled circular angular dispersion "
            "1−R of individual ion-birth clouds versus collected gain.",
            "**How measured:** The 128-sample 5° azimuth histogram gives "
            "R=|Σp(φ)exp(iφ)|. R≈1 describes tightly concentrated angles; "
            "a low R gives higher dispersion. Only full-128-sample "
            "avalanches enter the comparison.",
            "**Caution:** Opposite narrow lobes can also have R≈0, "
            "so this is an angular-spread diagnostic—not a measurement "
            "that a continuous ring surrounds the wire. Confirm with "
            "individual real 2D and 3D ion-birth positions."
        )
    if kind=="wrap_single":
        return (
            "**What you see:** The real **sampled ion-birth points for "
            "one specific Garfield electron avalanche** around its "
            "nearest 30-μm-diameter anode wire, with optional 3D "
            "rotation along the wire axis.",
            "**How measured:** Original, SHA-verified Garfield callback "
            "coordinates (u,v,w,t) are loaded for just this starting "
            "electron, no more than 128 uniformly selected positions. "
            "The centre is the nearest wire in the same p3_noFW "
            "analytic lattice used by the simulation. Point colour is "
            "physical birth time relative to the original Heed electron.",
            "**Interpretation:** Births on both sides of the wire "
            "indicate geometric spreading. A narrow arc means one-sided "
            "multiplication. The wire is drawn to scale, but only "
            "a random sample of ion births is shown. No ion drift, "
            "weighting-field response, or post-avalanche cathode charge "
            "sharing is represented."
        )
    raise RuntimeError("No physical explanation for "+kind)

def selftest():
    import plotly.graph_objects as go
    a=np.array([[1./72]*72,[1./72]*72])
    assert np.allclose(circular_dispersion(a),[1.,1.],atol=1e-12)
    hot=np.zeros((1,72));hot[0,36]=1
    assert circular_dispersion(hot)[0]<1e-12
    m=np.array([.5,-.5]);v=np.array([1.,-1.])
    assert np.allclose((1-np.sign(v)*m)/2,[.25,.25])
    assert cpp_round(.5)==1 and cpp_round(-.5)==-1
    parents=np.tile(np.arange(1000),2)
    xx=np.tile(np.array([-1.,1.]),1000)
    yy=xx.copy()
    centres,mean,lo,hi,n=grouped_parent_bootstrap(
        xx,yy,parents,np.array([-2.,0.,2.]),minimum=20,boot=30)
    assert np.allclose(mean,[-1.,1.]) and (n==1000).all()
    assert np.nanmax(np.abs(hi-mean))<1e-12
    synthetic={
       "gain":np.tile(np.array([120,950]),1000),
       "ions":np.full(2000,150),
       "A":np.tile(np.array([-.8,.8]),1000),
       "v0":xx,
       "parent":parents,
       "sample_count":np.full(2000,128),
       "angle_fraction":np.tile(np.eye(1,72,36,dtype=np.float32),(2000,1)),
    }
    class Mock:
        def warning(self,_):pass
    options={"color":"#E46B08","line_width":2.,"legend":"test",
             "ion_weighted":True,"reliable":True,"zscale":"Linear"}
    synthetic["u0"]=xx/2.
    for mode in ("wrap_v0_A","wrap_v0_opposite",
                 "wrap_v0_angle","wrap_u0_angle","wrap_dispersion"):
        explanation=describe(mode,synthetic,{},options)
        assert len(explanation)==3
        fig=go.Figure()
        chart,note,aux=render(mode,fig,synthetic,{},options,Mock())
        assert len(fig.data)>0 and aux is None and note
    # Test true wire radius, 3D cylindrical surface and physical axes
    # without opening any real muon or needing a Streamlit runtime.
    obs={"points":np.array([[.02,.03,.001,21.],[.04,-.02,.003,22.]]),
         "t_heed":19.}
    cylinder=three_d_view(obs,.1)
    assert len(cylinder.data)==2 and cylinder.data[0].type=="surface"
    assert np.allclose(
        np.sqrt(cylinder.data[0].x**2+cylinder.data[0].y**2),
        WIRE_RADIUS_MM,atol=1.e-12)
    print("PASS exact vs sampled geometry semantics, 1000-parent "
          "confidence bands, true 30-um wire cylinder, circular spread, "
          "population chart renderers; NO Garfield and no files modified.")

if __name__=="__main__":
    selftest()
