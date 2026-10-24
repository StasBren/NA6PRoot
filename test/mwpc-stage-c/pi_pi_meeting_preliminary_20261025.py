#!/usr/bin/env python3
"""PI meeting 2026-10-25: two PRELIMINARY questions from existing 100k data.

1. Electron drift-time PDF: all six independent 3/3.5/4 mm anode pitch
   geometries with/without 75 um field wire, chamber and station.
2. P(electron approach angle | signed muon entrance angle): all six
   geometries, with theta_u and theta_w shown SEPARATELY. Incidence angle
   is last Garfield microscopic step relative to LOCAL ANODE SURFACE NORMAL,
   not muon entrance angle, not macroscopic E-field direction, and not
   an avalanche/ASIC arrival-time measurement.

This is a READ-ONLY presentation postprocessor. Uses only saved 100k
per-source (2*50k) native merged Geant4/Heed/Garfield data, never reruns
transport, Heed or Garfield, never touches running 6M production.
Includes exact 10 ns >=300 ns overflow; normalizes collected known-time
electrons only. No independent-electron statistical error bars are inferred
(parent-event bootstrap not available yet).

Examples:
 python3 test/mwpc-stage-c/pi_pi_meeting_preliminary_20261025.py --source Jpsi --station 4 --det 143
 python3 test/mwpc-stage-c/pi_pi_meeting_preliminary_20261025.py --source Phi --station 0 --det 1
"""
from __future__ import annotations
import argparse
import json
import math
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from pi_2m_native_six_geometry_fold import OUTPUT
from pi_six_geometry_native_dashboard_pages import read_groups_stream
from pi_six_geometry_native_timing import n_anodes
from pi_native_centered_timing_pilot import ANGLES
from pi_pitch_fieldwire_pilot import tag

CONFIGS=[tag(p,d) for p in (3.,3.5,4.) for d in (0.,75.)]
SOURCE_TAGS={"Jpsi":"jpsi40-bore","Phi":"phi40-bore","Omega":"omega40-bore"}
STYLE={
    "p3p0_anode_only":"#007f86",
    "p3p0_field75p0um":"#22a7a0",
    "p3p5_anode_only":"#af4c91",
    "p3p5_field75p0um":"#e48bba",
    "p4p0_anode_only":"#a56219",
    "p4p0_field75p0um":"#e3b14d",
}
REPORT_ROOT=OUTPUT.parent/"pi_pi_meeting_20261025_v1"
TIME_EDGES=np.arange(0,301,10,dtype=float)
INC_EDGES=np.arange(0,181,10,dtype=float)


def config_label(cfg):
    pitch=cfg.split("_")[0][1:].replace("p",".")
    return f'{pitch} mm · '+("75 µm field wire" if "field" in cfg else "no field wire")


def validate_group(row):
    if len(row["hist"])!=31 or sum(row["hist"])!=row["known"]:
        raise ValueError("Time histogram lost collected known-time electrons")
    if row["collected"]!=row["known"]+row["unknown"]:
        raise ValueError("Collected != known + unknown")
    if len(row["angle_hist"])!=18 or sum(row["angle_hist"])!=row["angle_n"]:
        raise ValueError("Garfield polar approach-angle count mismatch")
    if len(row["phi_hist"])!=36 or sum(row["phi_hist"])!=row["phi_n"]:
        raise ValueError("Garfield approach-azimuth count mismatch")


def selected_data(filename,source,station,det):
    # Stream huge input; keep only station, chosen chamber, its anodes,
    # and conditional muon-angle station/chamber bins in bounded RAM.
    iterrows=iter(read_groups_stream(filename))
    header=next(iterrows)
    meta=header["meta"]
    if (meta.get("n_generated_parent_events_per_source")!=100000 or
        meta.get("source_station")!=f"{SOURCE_TAGS[source]} MS{station}" or
        meta.get("geometries")!=6 or
        set(meta.get("garfield_luts",{}))!=set(CONFIGS)):
        raise ValueError("Wrong/mixed six-geometry 100k source-station provenance")
    edges=list(meta.get("angle_bins_deg",[]))
    if edges!=ANGLES.tolist():
        raise ValueError("Muon signed entrance-angle bins differ from frozen model")
    output={
      "schema":"NA60_PI_MEETING_20261025_SIX_GEOMETRY_READONLY_V1",
      "source":source,"station":station,"chosen_detector":det,
      "exposure":"100,000 generated parents for this one source; "
                 "2 completed 50k native Geant4 chunks",
      "angle_semantics":"theta_u,w=atan2(local momentum transverse component,"
                        " local gas normal component); "
                        "electron incidence_normal_deg is the last "
                        "Garfield microscopic step relative to local anode normal; "
                        "phi_hit_deg is electron incidence azimuth",
      "muon_angle_bin_edges_deg":edges,
      "electron_incidence_bin_edges_deg":INC_EDGES.tolist(),
      "electron_phi_bin_edges_deg":list(range(-180,181,10)),
      "time_bin_width_ns":10,
      "time_overflow":"hist[30] = t>=300ns (included in denominator)",
      "normalization":"Conditional on collected electrons with known drift "
            "time for timing; conditional on collected electrons with finite "
            "Garfield approach angle in each signed muon angle bin for angles.",
      "limitations":"No parent bootstrap/statistical error bars; B=0; "
            "frozen sampled Garfield LUT; no avalanche, ion or ASIC time",
      "data":{cfg:{"station":None,"chamber":None,"anodes":[],
                   "theta_u_station":{},"theta_w_station":{},
                   "theta_u_chamber":{},"theta_w_chamber":{}} for cfg in CONFIGS},
    }
    seen=set()
    for row in iterrows:
        if row["station"]!=station or row["source"]!=SOURCE_TAGS[source]:
            raise ValueError("Mixed source/station detected in merged file")
        cfg=row["config"]
        if cfg not in output["data"]:raise ValueError(f"Unexpected LUT {cfg}")
        did=int(row["det_id"]);wire=int(row["anode"])
        axis=row["axis"]; abin=int(row["anglebin"])
        target=None
        if did==-1 and wire==-1:
            if axis=="all" and abin==-1:target="station"
            elif axis in ("u","w") and 0<=abin<len(edges)-1:
                target=f"theta_{axis}_station"
        elif did==det:
            if wire==-1:
                if axis=="all" and abin==-1:target="chamber"
                elif axis in ("u","w") and 0<=abin<len(edges)-1:
                    target=f"theta_{axis}_chamber"
            elif wire>=0 and axis=="all" and abin==-1:
                n=n_anodes(station,float(cfg.split("_")[0][1:].replace("p",".")))
                if wire>=n:raise ValueError(f"Anode {wire} beyond {n}")
                target="anodes"
        if target is None:continue
        validate_group(row)
        ident=(cfg,did,wire,axis,abin)
        if ident in seen:raise ValueError(f"Repeated six-LUT native group {ident}")
        seen.add(ident)
        if target=="anodes":
            output["data"][cfg][target].append({
                "anode":wire,"known":row["known"],"hist":row["hist"],
                "mean_ns":row.get("mean_ns"),"rms_ns":row.get("rms_ns"),
            })
        elif target.startswith("theta_"):
            output["data"][cfg][target][str(abin)]={
                "segments":row["segments"],"angle_n":row["angle_n"],
                "angle_hist":row["angle_hist"],"phi_n":row["phi_n"],
                "phi_hist":row["phi_hist"],
            }
        else:
            output["data"][cfg][target]={
                k:row[k] for k in ("segments","primaries","collected","known",
                    "unknown","hist","mean_ns","rms_ns","angle_n","angle_hist",
                    "phi_n","phi_hist")
            }
    for cfg in CONFIGS:
        data=output["data"][cfg]
        if not data["station"] or not data["chamber"]:
            raise ValueError(f"Missing MS{station} or DET{det} summary for {cfg}")
        if len(data["theta_u_station"])!=len(edges)-1 or \
           len(data["theta_w_station"])!=len(edges)-1:
            # Zero-hit theta bins are absent from sparse native group data.
            for axis in ("u","w"):
                for i in range(len(edges)-1):
                    data[f"theta_{axis}_station"].setdefault(str(i),{
                        "segments":0,"angle_n":0,"angle_hist":[0]*18,
                        "phi_n":0,"phi_hist":[0]*36})
        chamber_hist=np.zeros(31,dtype=np.int64)
        chamber_known=0
        for anode in data["anodes"]:
            chamber_hist+=np.asarray(anode["hist"],dtype=np.int64)
            chamber_known+=anode["known"]
        if chamber_known!=data["chamber"]["known"] or \
           chamber_hist.tolist()!=data["chamber"]["hist"]:
            raise ValueError(f"Per-anode chamber histogram conservation FAILED: {cfg}")
        for axis in ("u","w"):
            rowset=data[f"theta_{axis}_station"]
            if sum(x["angle_n"] for x in rowset.values())!= \
               data["station"]["angle_n"]:
                raise ValueError(f"Muon {axis} conditional Garfield angles "
                                 f"do not sum to full station: {cfg}")
            if sum(x["segments"] for x in rowset.values())!= \
               data["station"]["segments"]:
                raise ValueError(f"Muon {axis} conditioned crossing exposure "
                                 f"does not sum to station: {cfg}")
        data["anodes"].sort(key=lambda x:x["anode"])
    output["input_file"]=str(filename)
    output["input_scope_validated"]=True
    return output


def axes_time(ax,row,cfg):
    known=row["known"]
    h=np.array(row["hist"][:30],dtype=float)
    if known<=0:raise ValueError("No collected known-time electrons")
    ax.bar(TIME_EDGES[:-1],100*h/known,width=9.4,align="edge",
           color=STYLE[cfg],alpha=.82)
    ax.axvline(100,ls="--",lw=1,color="#d4772b")
    ax.set_xlim(0,145)
    ax.set_title(f"{config_label(cfg)}\n"
                 f"μ={row['mean_ns']:.2f} ns · RMS={row['rms_ns']:.2f} ns",
                 fontsize=10)
    ax.set_xlabel("Electron drift time since birth [ns]")
    ax.set_ylabel("Collected electrons / 10 ns [%]")
    ax.grid(axis="y",alpha=.15)
    ax.set_axisbelow(True)


def plot_six_times(report,level,out):
    fig,axes=plt.subplots(2,3,figsize=(14,7.7),sharex=True,sharey=True)
    rows=[report["data"][cfg][level] for cfg in CONFIGS]
    high=max(max(r["hist"][:30])/r["known"]*100 for r in rows)
    for cfg,ax,r in zip(CONFIGS,axes.ravel(),rows):
        axes_time(ax,r,cfg)
        ax.set_ylim(0,high*1.17)
    marker=(f"DET{report['chosen_detector']}" if level=="chamber" else
            f"all MS{report['station']} chambers")
    fig.suptitle(f"Electron drift-time PDFs · {report['source']} · MS{report['station']} · {marker}\n"
        "Same native muon crossings + shared Heed births · preliminary 100k generated parents",
        fontsize=14,weight="bold")
    fig.text(.5,.008,"Each panel normalised by its own collected known-time electrons; "
        "t≥300 ns overflow included in denominator, not in drawn 0–150 ns range.",
        ha="center",fontsize=9)
    fig.tight_layout(rect=(0,.035,1,.94))
    fig.savefig(out/f"01_time_six_{level}.png",dpi=200)
    plt.close(fig)


def plot_overlay_times(report,out):
    fig,ax=plt.subplots(figsize=(12,6.3))
    for cfg in CONFIGS:
        row=report["data"][cfg]["station"]
        h=np.asarray(row["hist"][:30],dtype=float)
        ax.step(TIME_EDGES[:-1]+5,100*h/row["known"],where="mid",
                lw=2.0,color=STYLE[cfg],label=config_label(cfg))
    ax.axvline(100,ls=":",color="black",alpha=.6,label="100 ns")
    ax.set(xlim=(0,130),xlabel="Electron drift time [ns]",
           ylabel="Collected known-time electrons / 10 ns [%]",
           title=f"Six geometries · station MS{report['station']} · {report['source']}")
    ax.grid(alpha=.18);ax.legend(fontsize=9,ncol=2,frameon=False)
    fig.tight_layout()
    fig.savefig(out/"02_time_overlay_station.png",dpi=200)
    plt.close(fig)


def plot_numbered_anode_means(report,out):
    fig,axs=plt.subplots(2,1,figsize=(12,8),sharex=True)
    for cfg in CONFIGS:
        pitch=float(cfg.split("_")[0][1:].replace("p","."))
        number=n_anodes(report["station"],pitch)
        anodes=report["data"][cfg]["anodes"]
        if not anodes:continue
        x=np.array([2*r["anode"]/(number-1)-1 for r in anodes])
        count=np.array([r["known"] for r in anodes])
        mean=np.array([r["mean_ns"] if r["mean_ns"] is not None else
                       np.nan for r in anodes])
        rms=np.array([r["rms_ns"] if r["rms_ns"] is not None else
                      np.nan for r in anodes])
        # Some DET anodes have very few collected electrons. Omit <=30.
        ok=count>=30
        axs[0].plot(x[ok],mean[ok],linewidth=1.15,alpha=.85,
                    color=STYLE[cfg],label=config_label(cfg))
        axs[1].plot(x[ok],rms[ok],linewidth=1.15,alpha=.85,color=STYLE[cfg])
    axs[0].set_ylabel("Mean drift time [ns]")
    axs[1].set_ylabel("Centered RMS [ns]")
    axs[1].set_xlabel("Physical anode position in DET gas width (-1 edge; +1 edge)")
    for ax in axs:ax.grid(alpha=.15)
    axs[0].legend(ncol=2,fontsize=8)
    fig.suptitle(f"True numbered anodes · DET{report['chosen_detector']} / "
                 f"MS{report['station']} · {report['source']}\n"
                 "Horizontal axis normalized separately for each pitch; "
                 "anodes with <30 collected known-time e⁻ omitted",fontsize=12)
    fig.tight_layout(rect=(0,0,1,.92))
    fig.savefig(out/"03_numbered_anodes_mean_rms_chamber.png",dpi=200)
    plt.close(fig)


def angle_array(report,cfg,axis,level,angle_field="angle_hist"):
    rows=report["data"][cfg][f"theta_{axis}_{level}"]
    bins=len(report["muon_angle_bin_edges_deg"])-1
    values=np.array([rows.get(str(i),{}).get(angle_field,[
          0]*(18 if angle_field=="angle_hist" else 36)) for i in range(bins)],
          dtype=float)
    counts=values.sum(axis=1)
    prob=np.divide(values,counts[:,None],out=np.full_like(values,np.nan),
                   where=counts[:,None]!=0)
    return prob,counts


def plot_angle_heatmaps(report,axis,out):
    muedges=report["muon_angle_bin_edges_deg"]
    labels=[f"{a:g}…{b:g}" for a,b in zip(muedges,muedges[1:])]
    matrices=[]
    for cfg in CONFIGS:
        pp,nn=angle_array(report,cfg,axis,"station")
        # Sparse source bins: blank bins with fewer than 200 electrons.
        pp[nn<200,:]=np.nan
        matrices.append(pp)
    top=max([np.nanmax(x) for x in matrices if np.isfinite(x).any()]+[.01])
    fig,axes=plt.subplots(2,3,figsize=(17,9.2),sharex=True,sharey=True)
    im=None
    for cfg,ax,mat in zip(CONFIGS,axes.ravel(),matrices):
        im=ax.imshow(100*mat.T,origin="lower",aspect="auto",
                     vmin=0,vmax=100*top,cmap="viridis")
        ax.set_title(config_label(cfg),fontsize=11)
        ax.set_xticks(np.arange(len(labels)));ax.set_xticklabels(
            labels,rotation=60,ha="right",fontsize=7)
        ax.set_yticks(np.arange(0,18,2))
        ax.set_yticklabels([f"{10*x}-{10*(x+1)}" for x in range(0,18,2)],
                           fontsize=8)
        ax.set_ylabel("Electron incidence relative to local normal [degrees]")
        ax.set_xlabel(f"Signed muon entrance θ_{axis} bin [degrees]")
    fig.colorbar(im,ax=axes.ravel().tolist(),label="Collected electron fraction per muon-angle bin [%]",
                 fraction=.02,pad=.015)
    fig.suptitle(f"PI question: P(e⁻ approach angle | muon θ_{axis}) · "
                 f"{report['source']} · MS{report['station']}\n"
                 "Six Garfield geometries, native 100k-parent subset · "
                 "white = <200 collected angle-known e⁻",fontsize=14)
    fig.subplots_adjust(top=.88,bottom=.12,left=.07,right=.89,
                        hspace=.41,wspace=.26)
    fig.savefig(out/f"04_electron_approach_vs_muon_theta_{axis}_heatmaps.png",dpi=190)
    plt.close(fig)


def plot_angle_overlays(report,axis,out):
    cfg=tag(3.5,0)
    prob,counts=angle_array(report,cfg,axis,"station")
    edges=report["muon_angle_bin_edges_deg"]
    fig,ax=plt.subplots(figsize=(11,6.5))
    palette=plt.get_cmap("tab10")
    shown=0
    for i in (2,3,4,5,6,7,8):
        if counts[i]<200:continue
        ax.step(INC_EDGES[:-1]+5,100*prob[i],where="mid",
           lw=2,color=palette(i%10),label=
           f"θ_{axis} [{edges[i]:g},{edges[i+1]:g})° · {int(counts[i]):,} e⁻")
        shown+=1
    if not shown:raise ValueError("No sufficiently populated muon incidence bin")
    ax.set(xlim=(0,180),xlabel="Electron incidence to anode normal [degrees]",
           ylabel="Collected electrons per 10° [%]",
           title=f"Electron approach-angle distributions conditional on muon θ_{axis}\n"
                 f"{report['source']} · MS{report['station']} · {config_label(cfg)}")
    ax.legend(fontsize=8,ncol=2);ax.grid(alpha=.18)
    fig.tight_layout()
    fig.savefig(out/f"05_electron_approach_angle_slices_theta_{axis}.png",dpi=200)
    plt.close(fig)


def plot_station_chamber_comparison(report,out):
    fig,ax=plt.subplots(figsize=(10.5,6.0))
    for level,ls in (("station","-"),("chamber","--")):
        row=report["data"][tag(3.5,0)][level]
        ax.step(TIME_EDGES[:-1]+5,
            100*np.asarray(row["hist"][:30])/row["known"],
            where="mid",lw=2.3,ls=ls,
            label=f"{level}: N_e={row['known']:,}")
    ax.set(xlim=(0,130),xlabel="Drift time [ns]",
           ylabel="Known collected electrons / 10 ns [%]",
           title=f"Chamber-vs-station check · MS{report['station']} · "
                 f"DET{report['chosen_detector']}\n"
                 "3.5 mm pitch, no field wire")
    ax.grid(alpha=.2);ax.legend()
    fig.tight_layout();fig.savefig(out/"06_chamber_vs_station.png",dpi=200)
    plt.close(fig)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source",choices=tuple(SOURCE_TAGS),default="Jpsi")
    p.add_argument("--station",type=int,choices=range(6),default=4)
    p.add_argument("--det",type=int,default=143)
    p.add_argument("--input",type=Path,default=None,
        help="100k native pre-merged source/station JSON")
    p.add_argument("--out",type=Path,default=REPORT_ROOT)
    args=p.parse_args()
    path=args.input or (OUTPUT/"full"/"merged"/args.source.lower()/
                        f"MS{args.station}.json")
    if not path.is_file():raise SystemExit(
        f"Missing validated 100k source station file: {path}")
    report=selected_data(path,args.source,args.station,args.det)
    args.out.mkdir(parents=True,exist_ok=True)
    for level in ("station","chamber"):
        plot_six_times(report,level,args.out)
    plot_overlay_times(report,args.out)
    plot_numbered_anode_means(report,args.out)
    for axis in ("u","w"):
        plot_angle_heatmaps(report,axis,args.out)
        plot_angle_overlays(report,axis,args.out)
    plot_station_chamber_comparison(report,args.out)
    status={
        "source":args.source,"station":args.station,"det_id":args.det,
        "plots":[x.name for x in sorted(args.out.glob("*.png"))],
        "physics_definition":report["angle_semantics"],
        "input_scope":report["exposure"],
    }
    tmp=args.out/"selected_data.json.incomplete"
    tmp.write_text(json.dumps(report,separators=(",",":"),
          ensure_ascii=False,allow_nan=False),encoding="utf-8")
    os.replace(tmp,args.out/"selected_data.json")
    (args.out/"README.txt").write_text(
      "PI preliminary: 100k generated parents / Jpsi or selected source.\n"
      "Six pitches+field wire compared with SAME Geant4/Heed input.\n"
      "Muon signed theta_u and theta_w treated separately.\n"
      "Electron incidence: last microscopic step relative anode normal.\n"
      "All PDFs conditional on collected Garfield electrons with finite "
      "observable. No parent bootstrap uncertainties yet.\n"
      "DO NOT label last step alpha as absolute angle to anode wire axis.\n",
      encoding="utf-8")
    print(json.dumps(status,indent=2,ensure_ascii=False))
    print("PASS: PI meeting plots and compact data ready at",args.out)


if __name__=="__main__":main()
