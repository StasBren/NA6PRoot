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
PITCH_MM=3.
WIRE_RADIUS_MM=.015
V_EDGES=np.linspace(-3.,3.,25)
ABS_V_EDGES=np.linspace(0.,3.,13)
THETA=(np.asarray(temporal.PHI_EDGES[:-1])+np.asarray(temporal.PHI_EDGES[1:]))/2.

KINDS={
    "wrap_v0_A":"Signed exact ion asymmetry vs Heed starting v₀",
    "wrap_v0_opposite":"Ion births on opposite side vs |v₀|",
    "wrap_v0_angle":"Ion-birth azimuth vs starting v₀ (heatmap)",
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
    """Confidence bands resample 1000 independent real Geant4 muons.

    One electron is NOT one independent parent: each real muon generates
    many correlated ionisation electrons/avalanches.
    """
    x=np.asarray(x,float);y=np.asarray(y,float);parent=np.asarray(parent,int)
    edges=np.asarray(edges,float)
    valid=np.isfinite(x)&np.isfinite(y)&(parent>=0)&(parent<1000)
    xi=np.searchsorted(edges,x[valid],side="right")-1
    xi=np.where(x[valid]==edges[-1],len(edges)-2,xi)
    p=parent[valid];v=y[valid]
    good=(xi>=0)&(xi<len(edges)-1)
    xi=xi[good];p=p[good];v=v[good]
    nbin=len(edges)-1
    cnt=np.zeros((1000,nbin),dtype=np.float64)
    summ=np.zeros((1000,nbin),dtype=np.float64)
    np.add.at(cnt,(p,xi),1.)
    np.add.at(summ,(p,xi),v)
    counts=cnt.sum(axis=0)
    sums=summ.sum(axis=0)
    mean=np.divide(sums,counts,out=np.full(nbin,np.nan),
                   where=counts>=minimum)
    rng=np.random.default_rng(20261110)
    # Each bootstrap resamples complete real-muon parents, not electrons.
    draws=rng.multinomial(1000,np.full(1000,.001),size=boot)
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

def options(st,kind,d,meta):
    """Only physically relevant controls for each specific view."""
    out={"min_entries":40,"ion_weighted":True,"reliable":True}
    if kind in {"wrap_v0_A","wrap_v0_opposite"}:
        st.sidebar.caption("Shaded band: 95% bootstrap CI from "
                           "1000 independent Geant4 muon parents.")
        st.sidebar.caption("Every avalanche with defined exact A is included; "
                           "no 128-point sampling requirement.")
    if kind=="wrap_v0_angle":
        out["ion_weighted"]=st.sidebar.checkbox(
            "Weight by number of ion births",True,
            key="wrap_v0_angle-ion-weighting")
        out["reliable"]=st.sidebar.checkbox(
            "Require all 128 sampled ion births",True,
            key="wrap_v0_angle-reliable")
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

def three_d_view(observation,halfwidth):
    import plotly.graph_objects as go
    pts=observation["points"]
    zmin,zmax=float(pts[:,2].min()),float(pts[:,2].max())
    extension=max(.05,.2*(zmax-zmin))
    z0=zmin-extension;z1=zmax+extension
    phi=np.linspace(0.,2*np.pi,45)
    cylinder_u=np.tile(WIRE_RADIUS_MM*np.cos(phi),(2,1))
    cylinder_v=np.tile(WIRE_RADIUS_MM*np.sin(phi),(2,1))
    cylinder_w=np.tile(np.linspace(z0,z1,2)[:,None],(1,len(phi)))
    dt=pts[:,3]-observation["t_heed"]
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
                                colorbar=dict(title="Birth delay [ns]")),
        name="Saved Garfield ion births",
        hovertemplate="Δu=%{x:.4g} mm<br>v=%{y:.4g} mm"
                      "<br>w=%{z:.4g} mm<extra></extra>"))
    fig.update_layout(
        title="3D sampled ion-birth cloud around anode wire",
        template="plotly_white",height=540,
        scene=dict(xaxis=dict(title="Δu [mm]",range=[-halfwidth,halfwidth]),
                   yaxis=dict(title="v [mm]",range=[-halfwidth,halfwidth]),
                   zaxis=dict(title="w [mm]",range=[z0,z1]),
                   aspectmode="data"),
        margin=dict(l=0,r=0,t=70,b=0))
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
            "Shaded interval resamples 1,000 complete muons, not electrons."),None
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
    if kind=="wrap_v0_angle":
        good=(d["ions"]>0)&np.isfinite(x)&np.isfinite(A)
        if controls["reliable"]:good&=(d["sample_count"]==128)
        v=x[good];fraction=d["angle_fraction"][good]
        weight=d["ions"][good].astype(float) if controls["ion_weighted"] else np.ones(len(v))
        ix=np.searchsorted(V_EDGES,v,side="right")-1
        ix=np.where(v==V_EDGES[-1],len(V_EDGES)-2,ix)
        valid=(ix>=0)&(ix<len(V_EDGES)-1)
        bins=np.zeros((len(V_EDGES)-1,72))
        np.add.at(bins,ix[valid],fraction[valid]*weight[valid,None])
        # Condition on given INITIAL v0; rows sum to one.
        denom=bins.sum(axis=1)
        matrix=np.divide(bins,denom[:,None],
                         out=np.full(bins.shape,np.nan),
                         where=denom[:,None]>0)
        z=matrix
        title=controls["legend"]
        if controls.get("zscale")=="Log":
            with np.errstate(divide="ignore",invalid="ignore"):
                z=np.where(matrix>0,np.log10(matrix),np.nan)
            title="log10("+title+")"
        fig.add_trace(go.Heatmap(
            x=THETA,y=(V_EDGES[:-1]+V_EDGES[1:])/2,
            z=z,colorscale="YlOrBr",
            colorbar=dict(title=title),
            hovertemplate="φ=%{x:.1f}°<br>initial v₀=%{y:.2f} mm"
                          "<br>conditional bin prob=%{z:.4g}"
                          "<extra></extra>"))
        return None,(
            f"{int(np.sum(good)):,} avalanches. Heatmap rows normalised "
            "independently; estimates of **where ions were CREATED** "
            "around a wire, conditioned on initial Heed v₀. "
            +("Ion weighted." if controls["ion_weighted"] else
              "Each ion-producing avalanche equally weighted.")+
            " 128-point sampling is finite."),None
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
    if kind=="wrap_v0_angle":
        return (
            "**What you see:** A conditional angular map. Each horizontal "
            "row fixes the starting electron v₀; horizontal position is "
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
    for mode in ("wrap_v0_A","wrap_v0_opposite",
                 "wrap_v0_angle","wrap_dispersion"):
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
