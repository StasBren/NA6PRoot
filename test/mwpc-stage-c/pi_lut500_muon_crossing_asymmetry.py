#!/usr/bin/env python3
"""LUT20/LUT500 per-MUON-GAS-CROSSING anode hemisphere and polar-sector asymmetry.

Read-only fold of original native Geant4 J/psi gas-crossing CSV, original
stable-seeded Heed birth sampler and immutable p=3 / p=4
20- or 500-microscopic-trials/cell Garfield tables (explicit selection). Every crossing has its own
(source, chunk, event, hit_index, track_id, station, DET) provenance.
No Geant4/Garfield/Heed library re-simulation and NO cathode signal model.

Hemispheres: up 0<phi_hit<180; down -180<phi_hit<0. Exact phi=0,+/-180
unclassified. Distinct illustrative sectors: top [45,135), bottom
[-135,-45), side = remaining surface. A values with zero denominators
are UNDEFINED, never set to zero. Side hits remain in hemisphere A,
but NOT in the separate top-vs-bottom sector A.
LUT trials are repeatedly reused and Heed siblings/parent crossings are
correlated; event population RMS is NOT a standard error.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

from pi_2m_native_six_geometry_fold import crossing_provenance
from pi_2m_zero_path_compatible_resume import coincident_zero_chord
from pi_native_centered_timing_pilot import ANGLES
from pi_six_geometry_native_timing import (
    HEED_TRACKS, HEED_CLUSTERS, birth_sample, load_map,
    response_for_births, stable_seed, MAPS as ROOT20,
)
from pi_lut_20_vs_500_convergence import ROOT500 as ROOT3
from pi_pitch_fieldwire_pilot import tag
from pi_event_asymmetry_500_ms0 import up_down
from run_pi_per_chamber_response import read_heed
from run_pi_geant4_campaign import REPO
from mwpc_checkerboard_projection import build_chambers

HERE = Path(__file__).resolve().parent
BASE = REPO / "test_runs/mwpc_stage_c"
ROOT4 = BASE / "pi_pitch_fieldwire_500_p4_v1/pilot"
OUT = BASE / "pi_lut500_muon_crossing_asymmetry_v1"
SOURCE = "jpsi40-bore"
BINS = np.linspace(-1., 1., 42)
DELTA_BINS = np.linspace(-2., 2., 42)
SCHEMA = "PI_JPSI_NATIVE_MUON_GAS_CROSSING_ASYM_P3_P4_LUT500_V1"
RAW_FIELDS = (
    "source","chunk","event","hit_index","track_id","station","det_id",
    "theta_u_deg","theta_w_deg","u_in_mm","u_out_mm","pitch_mm",
    "lut_trials_per_cell",
    "field_wire_um","n_primary","n_raw_anode","n_anode_finite",
    "n_edge_lost","up","down","hemisphere_boundary","top","bottom","side",
    "A_hemisphere","A_top_bottom","f_top_bottom",
)

def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()

def ratio(a,b):
    return a / b if b else None

def classify(samples, finite):
    """Count collected primary-electron endpoints, NOT cathode charge."""
    up, down, boundary = up_down(samples, finite)
    phi = samples[finite, 4]
    if not np.all(np.isfinite(phi) & (phi >= -180.) & (phi <= 180.)):
        raise ValueError("Bad collected anode endpoint azimuth")
    top = int(np.count_nonzero((phi >= 45.) & (phi < 135.)))
    bottom = int(np.count_nonzero((phi >= -135.) & (phi < -45.)))
    side = len(phi) - top - bottom
    if top + bottom + side != up + down + boundary or side < 0:
        raise AssertionError("Azimuth partition not conserved")
    return dict(up=up,down=down,hemisphere_boundary=boundary,
                top=top,bottom=bottom,side=side,
                A_hemisphere=ratio(up-down,up+down),
                A_top_bottom=ratio(top-bottom,top+bottom),
                f_top_bottom=ratio(top+bottom,len(phi)))

class Moments:
    """Equal-weight physical crossings; no electron-as-event pseudo-statistics."""
    def __init__(self):
        self.total=0
        self.primaries=0
        self.collected=0
        self.upper=self.lower=self.boundary=0
        self.top=self.bottom=self.side=0
        self.hist=np.zeros(41,dtype=np.int64)
        self.hist_tb=np.zeros(41,dtype=np.int64)
        self.valid=self.valid_tb=0
        self.sum=self.sumsq=self.sum_tb=self.sumsq_tb=0.
        self.f_n=self.f_sum=self.f_sumsq=0.

    def add(self,c,primaries,ncollected):
        self.total+=1; self.primaries+=int(primaries)
        self.collected+=int(ncollected)
        self.upper+=c["up"];self.lower+=c["down"]
        self.boundary+=c["hemisphere_boundary"]
        self.top+=c["top"];self.bottom+=c["bottom"];self.side+=c["side"]
        for key,attr,ss,sq,h in (
                ("A_hemisphere","valid","sum","sumsq","hist"),
                ("A_top_bottom","valid_tb","sum_tb","sumsq_tb","hist_tb")):
            a=c[key]
            if a is not None:
                setattr(self,attr,getattr(self,attr)+1)
                setattr(self,ss,getattr(self,ss)+a)
                setattr(self,sq,getattr(self,sq)+a*a)
                j=min(40,max(0,int(np.searchsorted(BINS,a,side="right")-1)))
                getattr(self,h)[j]+=1
        if c["f_top_bottom"] is not None:
            f=c["f_top_bottom"]
            self.f_n+=1;self.f_sum+=f;self.f_sumsq+=f*f

    @staticmethod
    def stats(n,s,ss):
        if not n:return {"mean":None,"rms":None}
        mean=s/n
        return {"mean":mean,"rms":math.sqrt(max(0.,ss/n-mean*mean))}

    def result(self):
        if (self.upper+self.lower+self.boundary != self.collected or
            self.top+self.bottom+self.side != self.collected or
            self.valid != int(self.hist.sum()) or
            self.valid_tb != int(self.hist_tb.sum()) or
            self.valid > self.total or self.valid_tb > self.total):
            raise AssertionError("Crossing electron and histogram conservation failure")
        return {
            "crossings":self.total,"n_primaries":self.primaries,
            "n_anode_collected":self.collected,
            "n_upper":self.upper,"n_lower":self.lower,
            "n_boundary":self.boundary,
            "n_top":self.top,"n_bottom":self.bottom,"n_side":self.side,
            "top_bottom_electron_fraction":ratio(self.top+self.bottom,self.collected),
            "hemisphere_A":{"valid":self.valid,"undefined":self.total-self.valid,
               "hist":self.hist.tolist(),**self.stats(self.valid,self.sum,self.sumsq)},
            "top_bottom_A":{"valid":self.valid_tb,"undefined":self.total-self.valid_tb,
               "hist":self.hist_tb.tolist(),
               **self.stats(self.valid_tb,self.sum_tb,self.sumsq_tb)},
            "f_top_bottom_by_crossing":{
                "valid":self.f_n,"undefined":self.total-self.f_n,
                **self.stats(self.f_n,self.f_sum,self.f_sumsq)}
        }

class Paired:
    def __init__(self):
        self.total=self.valid=0
        self.delta=np.zeros(41,dtype=np.int64)
        self.sum=self.sumsq=0.

    def add(self,a0,a1):
        self.total+=1
        if a0 is None or a1 is None:return
        d=a1-a0
        self.valid+=1;self.sum+=d;self.sumsq+=d*d
        j=min(40,max(0,int(np.searchsorted(DELTA_BINS,d,side="right")-1)))
        self.delta[j]+=1

    def result(self):
        if int(self.delta.sum())!=self.valid:
            raise AssertionError("FW paired difference not conserved")
        return {"both_defined":self.valid,"not_both_defined":self.total-self.valid,
                "FW_minus_noFW_hemisphere_A":
                    {**Moments.stats(self.valid,self.sum,self.sumsq),
                     "hist":self.delta.tolist()}}



def facing_fraction_plot(summary,out,pitch):
    """Mean fraction of endpoints in ±45° cathode-facing sectors, per crossing."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(9.2,5.7))
    for fw,color,ls,label in (
            (0.,"#2158a5","-","No field wire"),
            (75.,"#c52b2a","--","Field wire 75 μm")):
        cfg=tag(pitch,fw)
        values=[summary["stations"][str(st)]["configs"][cfg][
               "f_top_bottom_by_crossing"]["mean"] for st in range(6)]
        ax.plot(range(6),values,color=color,ls=ls,marker="o",
                markersize=5,lw=1.65,label=label)
    ax.set_xticks(range(6),[f"MS{i}" for i in range(6)])
    ax.set_ylim(0,1.)
    ax.set_xlabel("Muon spectrometer station")
    ax.set_ylabel("Mean fraction of anode endpoints in top/bottom ±45° sectors")
    ax.tick_params(direction="in",top=True,right=True)
    ax.legend(frameon=False)
    ax.set_title("NA60+ Simulation · Work in progress\n"
                 f"J/ψ, p={pitch:g} mm, LUT{summary['lut_trials_per_cell']} · per-crossing facing-sector fraction",
                 fontsize=12)
    fig.text(.5,.02,"One physical crossing has equal weight (when ≥1 anode endpoint). "
             "Geometric proxy only; NOT charge on cathodes.",
             ha="center",fontsize=8)
    fig.tight_layout(rect=(0,.055,1,1))
    name=f"f_top_bottom_by_crossing_p{int(pitch)}.png"
    fig.savefig(out/name,dpi=165,facecolor="white")
    plt.close(fig)
    return name

def up_down_scatter(crossings_by_pair, out, pitch, fw, summary):
    """Two-dimensional joint event-count distribution; NO parent pooling."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    name=tag(pitch,fw)
    fig,axes=plt.subplots(2,3,figsize=(13,7.8),sharex=True,sharey=True)
    all_counts=[]
    all_pop=[]
    for st in range(6):
        for (u,d),n in crossings_by_pair[st,name].items():
            all_counts.extend((u,d))
            all_pop.append(n)
    max_count=max(all_counts,default=1)
    max_pop=max(all_pop,default=1)
    lim=math.log10(max_count+1.)*1.035
    upper=max(1.,lim)
    scatter_last=None
    for st,ax in enumerate(axes.flat):
        counter=crossings_by_pair[st,name]
        coords=list(counter)
        if coords:
            xx=np.log10(1.+np.asarray([p[0] for p in coords],dtype=float))
            yy=np.log10(1.+np.asarray([p[1] for p in coords],dtype=float))
            population=np.asarray([counter[p] for p in coords])
            order=np.argsort(population)[::-1]
            scatter_last=ax.scatter(xx[order],yy[order],c=population[order],
                       cmap="viridis",norm=LogNorm(vmin=1,vmax=max(2,max_pop)),
                       s=8,marker="s",edgecolors="none",alpha=.85,rasterized=True)
        ax.plot([0,upper],[0,upper],ls="--",color="black",lw=1.1)
        ax.set_xlim(0,upper);ax.set_ylim(0,upper)
        ticks=[n for n in (0,1,3,10,30,100,300,1000,3000) if
               math.log10(n+1)<=upper]
        ax.set_xticks([math.log10(n+1) for n in ticks],labels=[str(n) for n in ticks])
        ax.set_yticks([math.log10(n+1) for n in ticks],labels=[str(n) for n in ticks])
        ax.tick_params(direction="in",top=True,right=True)
        ax.set_title(f"MS{st} · {sum(counter.values()):,} crossings",fontsize=10)
    for ax in axes[1]:ax.set_xlabel("Upper-hemisphere anode electrons")
    for ax in axes[:,0]:ax.set_ylabel("Lower-hemisphere anode electrons")
    fig.suptitle("NA60+ Simulation · Work in progress\n"
                 f"J/ψ, p={pitch:g} mm · {'FW 75 μm' if fw else 'No FW'}"
                 f" · per-crossing N(up) vs N(down) · LUT{summary['lut_trials_per_cell']}",
                 fontsize=13,y=.99)
    if scatter_last is not None:
        fig.colorbar(scatter_last,ax=axes.ravel().tolist(),fraction=.015,pad=.02,
                     label="Crossings per electron-count pair")
    fig.text(.5,.015,"Axes plotted in log10(N+1), labels are actual electron counts. "
             "One physical gas crossing per data point; colour encodes repeated count pairs.",
             ha="center",fontsize=8)
    fig.subplots_adjust(left=.09,right=.88,top=.88,bottom=.13,wspace=.24,hspace=.25)
    filename=f"up_down_joint_p{int(pitch)}_{'FW75' if fw else 'noFW'}.png"
    fig.savefig(out/filename,dpi=165,facecolor="white")
    plt.close(fig)
    return filename


def conditional_angle_curves(summary,out,pitch,axis):
    """Binned conditional RMS and mean A; x positions denote bins, not angles."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig,axs=plt.subplots(2,3,figsize=(14,8.1),sharex=True,sharey=True)
    edges=np.asarray(summary["angle_edges_deg"],dtype=float)
    xs=np.arange(len(edges)-1)
    names=[f"{edges[i]:g}..{edges[i+1]:g}" for i in xs]
    ymax=0.
    for st,ax in enumerate(axs.flat):
        station=summary["stations"][str(st)]
        for fw,color,ls,label in (
                (0.,"#2158a5","-","No FW"),
                (75.,"#c52b2a","--","FW 75 μm")):
            bins=station["angle"].get(tag(pitch,fw),{}).get(axis,{})
            y=np.asarray([
                bins[str(i)]["hemisphere_A"]["rms"]
                if str(i) in bins and bins[str(i)]["hemisphere_A"]["valid"]>=30
                else np.nan for i in xs],dtype=float)
            if np.isfinite(y).any():
                ymax=max(ymax,float(np.nanmax(y)))
            ax.plot(xs,y,color=color,linestyle=ls,marker="o",
                    markersize=3.5,label=label,lw=1.4)
        ax.set_title(f"MS{st} · {station['gas_crossings']:,} crossings",fontsize=11)
        ax.tick_params(direction="in",top=True,right=True)
        ax.set_xlim(-.3,len(xs)-.7)
        ax.set_xticks(xs,labels=names,rotation=48,ha="right",fontsize=7)
    for ax in axs[1]:ax.set_xlabel(f"Signed muon θ_{axis} interval [deg]")
    for ax in axs[:,0]:ax.set_ylabel("Population RMS(A), dimensionless")
    for ax in axs.flat:ax.set_ylim(0,max(.2,ymax*1.12))
    axs[0,0].legend(frameon=False,fontsize=8)
    fig.suptitle("NA60+ Simulation · Work in progress\n"
                 f"J/ψ · p={pitch:g} mm · LUT{summary['lut_trials_per_cell']} · "
                 f"RMS of per-crossing hemisphere A vs θ_{axis}",
                 fontsize=13,y=.995)
    fig.text(.5,.005,"Each marker represents an ENTIRE signed entrance-angle bin "
             "(not continuous theta). Bins with fewer than 30 valid crossings hidden; "
             "RMS is not standard error.",fontsize=8,ha="center")
    fig.tight_layout(rect=(0,.04,1,.90),h_pad=2.5)
    filename=f"RMS_A_vs_theta_{axis}_p{int(pitch)}.png"
    fig.savefig(out/filename,dpi=165,facecolor="white")
    plt.close(fig)
    return filename

def curves(summary, out, pitch, key, *, paired=False):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        "font.family":"DejaVu Sans","font.size":10,
        "axes.edgecolor":"black","axes.labelcolor":"black",
        "text.color":"black","xtick.color":"black","ytick.color":"black",
        "figure.facecolor":"white","axes.facecolor":"white",
        "axes.linewidth":1.05,
    })
    fig, axs=plt.subplots(2,3,figsize=(13,7.7),sharex=True,sharey=True)
    maxval=0.
    for st,ax in enumerate(axs.flat):
        station=summary["stations"][str(st)]
        if paired:
            row=station["paired"][str(pitch)]["FW_minus_noFW_hemisphere_A"]
            n=station["paired"][str(pitch)]["both_defined"]
            counts=np.asarray(row["hist"],dtype=float)
            yy=counts/n if n else counts
            ax.stairs(yy,DELTA_BINS,color="#214b83",lw=1.65)
            maxval=max(maxval,float(yy.max()))
            ax.set_title(f"MS{st} · paired crossings {n:,}",fontsize=11)
        else:
            for fw,color,style,label in (
                (0., "#2158a5","-","No field wire"),
                (75.,"#c52b2a","--","Field wire 75 μm")):
                cfg=tag(pitch,fw)
                row=station["configs"][cfg][key]
                n=row["valid"]
                yy=np.asarray(row["hist"],dtype=float)/max(n,1)
                ax.stairs(yy,BINS,label=label,lw=1.6,linestyle=style,color=color)
                maxval=max(maxval,float(yy.max()))
            ax.set_title(f"MS{st} · {station['gas_crossings']:,} crossings",fontsize=11)
        ax.tick_params(direction="in",top=True,right=True,which="both")
        ax.set_xlim((-2.,2.) if paired else (-1.,1.))
        ax.grid(False)
    for ax in axs[1]:
        ax.set_xlabel("Paired ΔA (dimensionless)" if paired
                      else ("A (upper/lower hemisphere)" if key=="hemisphere_A"
                            else "A (top/bottom ±45° sectors)"))
    for ax in axs[:,0]:
        ax.set_ylabel("Fraction of valid crossings / "+
                      ("0.098" if paired else "0.049"))
    for ax in axs.flat:
        ax.set_ylim(0,maxval*1.09 if maxval else 1.)
    if not paired:axs[0,0].legend(loc="upper right",frameon=False,fontsize=8)
    fig.suptitle(f"NA60+ Simulation · Work in progress\n"
                 f"J/ψ · 100k parents maximum · p = {pitch:g} mm · "
                 f"LUT{summary['lut_trials_per_cell']} · B = 0 T · {'paired ΔA' if paired else key}",
                 fontsize=13,y=.995)
    fig.text(.5,.005,
             "One physical muon gas crossing is one observation. "
             "Population RMS is not uncertainty on the mean. "
             "No gain or cathode signal.",
             ha="center",fontsize=8)
    fig.tight_layout(rect=(0,.035,1,.90))
    name=(f"paired_delta_A_p{int(pitch)}.png" if paired
          else f"{key}_p{int(pitch)}.png")
    fig.savefig(out/name,dpi=165,facecolor="white")
    plt.close(fig)
    return name

def run(chunks,stations,max_crossings,pitches,out,plots=True,*,lut_trials=500):
    if not (1<=chunks<=40) or max_crossings<0 or not stations or not pitches or lut_trials not in (20,500):
        raise ValueError("Invalid exposure selection")
    if any(st not in range(6) for st in stations):
        raise ValueError("Invalid station")
    out=Path(out)
    if any((out/name).exists() for name in
           ("crossings.csv.gz","summary.json","manifest.json")):
        raise FileExistsError(f"Non-destructive output already exists: {out}")
    if (out/"crossings.csv.gz.incomplete").exists():
        raise FileExistsError("Previous incomplete output must be inspected")
    configs={}
    provenance={}
    for pitch in pitches:
        root=ROOT20 if lut_trials==20 else (ROOT3 if pitch==3. else ROOT4)
        for fw in (0.,75.):
            cfg=tag(pitch,fw)
            pool,details=load_map(pitch,fw,root)
            if pool.shape!=(round(pitch/.1)*30,lut_trials,5) or details["samples_per_box"]!=lut_trials:
                raise ValueError(f"Expected p={pitch:g}mm immutable LUT{lut_trials}")
            configs[cfg]=(pitch,fw,pool)
            provenance[cfg]=details
    layout={c.det_id:c.station for group in build_chambers() for c in group}
    lam,marks=read_heed(HEED_TRACKS,HEED_CLUSTERS)
    stats=defaultdict(Moments)
    by_det=defaultdict(Moments)
    angle=defaultdict(Moments)
    paired=defaultdict(Paired)
    joint=defaultdict(Counter)  # (station, cfg) -> {(up,down): crossings}
    exposures=Counter()
    provenance_native=[]
    out.mkdir(parents=True,exist_ok=True)
    tmp=out/"crossings.csv.gz.incomplete"
    with gzip.open(tmp,"wt",encoding="utf-8",newline="",compresslevel=4) as f:
        writer=csv.DictWriter(f,fieldnames=RAW_FIELDS)
        writer.writeheader()
        for chunk in range(1,chunks+1):
            path,meta=crossing_provenance("Jpsi",chunk)
            provenance_native.append({"chunk":chunk,"csv_sha256":meta["csv_sha256"]})
            done=Counter()
            with path.open(newline="",encoding="utf-8") as src:
                for row in csv.DictReader(src):
                    st=int(row["station"])
                    if st not in stations or abs(int(row["pdg"]))!=13:
                        continue
                    if max_crossings and done[st]>=max_crossings:
                        continue
                    det=int(row["det_id"])
                    if layout.get(det)!=st:raise ValueError("Invalid real DET ID")
                    theta=(float(row["theta_u_in_deg"]),
                           float(row["theta_w_in_deg"]))
                    if not all(math.isfinite(t) and -180<=t<=180 for t in theta):
                        raise ValueError("Invalid signed muon angles")
                    length=float(row["path_length_mm"])
                    if length==0.:
                        if not coincident_zero_chord(row):
                            raise ValueError("Invalid zero-chord crossing")
                        pos=np.empty((0,3));q=np.empty(0)
                    else:
                        seed=stable_seed("native-six-pitch-shared-Heed-v1",
                             SOURCE,f"chunk{chunk:04d}",row["event"],
                             row["hit_index"],row["track_id"])
                        pos,q=birth_sample(row,lam,marks,seed)
                    per_crossing={}
                    for cfg,(pitch,fw,pool) in configs.items():
                        finite,raw,wires,samples=response_for_births(
                            pos[:,0],pos[:,1],q,pitch,st,pool)
                        c=classify(samples,finite)
                        collected=int(np.count_nonzero(finite))
                        rawcollected=int(np.count_nonzero(raw))
                        if collected>rawcollected:
                            raise AssertionError("Physical collected > anode collected")
                        stats[st,cfg].add(c,len(q),collected)
                        joint[st,cfg][(c['up'],c['down'])]+=1
                        by_det[st,det,cfg].add(c,len(q),collected)
                        for axis,t in (("u",theta[0]),("w",theta[1])):
                            bi=min(len(ANGLES)-2,max(0,int(
                               np.searchsorted(ANGLES,t,side="right")-1)))
                            angle[st,cfg,axis,bi].add(c,len(q),collected)
                        writer.writerow(dict(
                            source=SOURCE,chunk=chunk,event=row["event"],
                            hit_index=row["hit_index"],track_id=row["track_id"],
                            station=st,det_id=det,theta_u_deg=theta[0],
                            theta_w_deg=theta[1],
                            u_in_mm=row["u_in_mm"],u_out_mm=row["u_out_mm"],
                            pitch_mm=pitch,field_wire_um=fw,
                            lut_trials_per_cell=lut_trials,
                            n_primary=len(q),n_raw_anode=rawcollected,
                            n_anode_finite=collected,
                            n_edge_lost=rawcollected-collected,
                            **{k:("" if v is None else v) for k,v in c.items()}))
                        per_crossing[(pitch,fw)]=c["A_hemisphere"]
                    for pitch in pitches:
                        paired[st,pitch].add(per_crossing[pitch,0.],
                                             per_crossing[pitch,75.])
                    done[st]+=1
                    exposures[st]+=1
            print(f"Finished chunk {chunk}/{chunks}: "+
                  ", ".join(f"MS{st}={done[st]:,}" for st in stations),
                  flush=True)
    os.replace(tmp,out/"crossings.csv.gz")
    summary={
      "schema":SCHEMA.replace("LUT500",f"LUT{lut_trials}"),"source":SOURCE,
      "lut_trials_per_cell":lut_trials,
      "generated_parent_slots":chunks*50_000,
      "source_chunks":chunks,"partial_exposure":bool(max_crossings),
      "limit_per_station_per_chunk":max_crossings,
      "stations_selected":stations,"pitches_mm":pitches,
      f"LUT{lut_trials}":provenance,"native_source_csv":provenance_native,
      "heed":{"tracks":str(HEED_TRACKS),"clusters":str(HEED_CLUSTERS),
              "mean_clusters_per_mm":float(lam)},
      "definitions":{
          "unit":"ONE native physical muon gas crossing; source/chunk/parent/track/hit/DET preserved",
          "hemisphere":"A=(N_up-N_down)/(N_up+N_down), phi_hit >0 vs <0; exclude exact 0,+/-180",
          "sectors":"top: [45,135) deg; bottom: [-135,-45) deg; side: rest",
          "sector_A":"(N_top-N_bottom)/(N_top+N_bottom), undefined when denominator=0",
          "f_top_bottom":"(N_top+N_bottom)/N_anode_collected",
          "undefined":"missing denominator => null/blank, NEVER A=0",
          "weighting":"equal-weight crossing observations (not electron-weighted)",
          "limitations":f"Primary electron endpoint counts only; no avalanche, induced charge, cathode/electronics; reused correlated LUT{lut_trials}."
      },
      "angle_edges_deg":ANGLES.tolist(),
      "stations":{}
    }
    for st in stations:
        row={
          "gas_crossings":int(exposures[st]),
          "configs":{cfg:stats[st,cfg].result() for cfg in configs},
          "paired":{str(p):paired[st,p].result() for p in pitches},
          "det":{},"angle":{}
        }
        for (station,det,cfg),v in by_det.items():
            if station==st:row["det"].setdefault(str(det),{})[cfg]=v.result()
        for (station,cfg,axis,bi),v in angle.items():
            if station==st:
                row["angle"].setdefault(cfg,{}).setdefault(axis,{})[
                    str(bi)]=v.result()
        for cfg in configs:
            if sum(joint[st,cfg].values()) != exposures[st]:
                raise AssertionError('Crossing joint-count conservation failure')
            if sum(d[cfg]["n_anode_collected"] for d in row["det"].values()
                   if cfg in d) != row["configs"][cfg]["n_anode_collected"]:
                raise AssertionError("DET/station finite-anode conservation error")
            if sum(d[cfg]["crossings"] for d in row["det"].values()
                   if cfg in d) != exposures[st]:
                raise AssertionError("DET/station crossing conservation error")
        summary["stations"][str(st)]=row
    if not bool(max_crossings) and stations==list(range(6)):
        if any(exposures[st]==0 for st in stations):
            raise ValueError("Missing full station data")
    (out/"summary.json").write_text(
        json.dumps(summary,indent=2,allow_nan=False,ensure_ascii=False)+"\n",
        encoding="utf-8")
    figures=[]
    if plots and stations==list(range(6)):
        for pitch in pitches:
            for key in ("hemisphere_A","top_bottom_A"):
                figures.append(curves(summary,out,pitch,key))
            figures.append(curves(summary,out,pitch,"",paired=True))
            figures.append(facing_fraction_plot(summary,out,pitch))
            for fw in (0.,75.):
                figures.append(up_down_scatter(joint,out,pitch,fw,summary))
            for axis in ("u","w"):
                figures.append(conditional_angle_curves(summary,out,pitch,axis))
    manifest={
        "schema":summary["schema"],
        "complete_six_station":max_crossings==0 and stations==list(range(6)),
        "generated_parents":chunks*50_000,
        "sha256_outputs":{p.name:sha(p) for p in
            [out/"crossings.csv.gz",out/"summary.json"]+
            [out/n for n in figures]},
        "source_csv_sha256":provenance_native,
        "garfield_lut_sha256":provenance
    }
    (out/"manifest.json").write_text(
        json.dumps(manifest,indent=2)+"\n",encoding="utf-8")
    print(f"PASS: {out}, {sum(exposures.values()):,} physical crossings, "
          f"{len(configs)} Garfield LUT{lut_trials} variants. "
          "Plots are descriptive, no false independent-electron error bars.")
    return summary

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--chunks",type=int,default=1)
    p.add_argument("--stations",type=int,nargs="+",default=list(range(6)))
    p.add_argument("--pitches",type=float,nargs="+",choices=(3.,4.),default=[3.])
    p.add_argument("--max-crossings",type=int,default=200,
                   help="Per station per 50k-parent chunk; 0 = full")
    p.add_argument("--out",type=Path,default=OUT)
    p.add_argument("--lut-trials",type=int,choices=(20,500),default=500,
                   help="Frozen independent microscopic drifts per cell; LUT20 uses the original six-geometry maps.")
    p.add_argument("--no-plots",action="store_true")
    a=p.parse_args()
    run(a.chunks,sorted(set(a.stations)),a.max_crossings,
        sorted(set(a.pitches)),a.out,not a.no_plots,
        lut_trials=a.lut_trials)
if __name__=="__main__":
    main()
