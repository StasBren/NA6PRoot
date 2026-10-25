#!/usr/bin/env python3
"""PI plots: azimuth of near-anode Garfield endpoints around anode circumference.

Six genuine 20-trial Garfield LUT response geometries, matched Heed births,
station and DET chamber levels, conditional on signed Geant4 muon entrance
theta_u and theta_w. Read only an existing validated selected_data.json.

CRITICAL: phi_hit_deg is atan2(v_end, u_end-u_wire), the near-anode ENDPOINT
POSITION azimuth. It is NOT the electron-velocity direction, and it is NOT
the angle incidence_normal_deg relative to local wire-surface normal.

Input from precomputed 2M source:
  test_runs/mwpc_stage_c/pi_full_channel_2m_PI_v1/Jpsi_MS4_DET143/selected_data.json
Or preliminary 100k:
  test_runs/mwpc_stage_c/pi_pi_meeting_20261025_v1/selected_data.json

The 2M physics file is created by pi_full_channel_2m_PI_plots.py after
the frozen 40 Jpsi chunk-station fold files have been verified.
"""
from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CONFIGS=[(f'p{str(p).replace(".","p")}_'+
          ("anode_only" if fw==0 else "field75p0um"),p,fw)
         for p in (3.,3.5,4.) for fw in (0,75)]
BINS=(2,4,5,7)  # populated theta_mu intervals at Jpsi MS4


def provenance(src):
    n=src.get("generated_parent_events")
    if n is None:
        if "100,000" in src.get("exposure",""): n=100_000
        else: raise ValueError("Source generated-parent exposure not validated")
    if (n not in (100_000,2_000_000) or src.get("source")!="Jpsi" or
        src.get("station") not in range(6) or not isinstance(src.get("chosen_detector"),int)
        or set(src.get("data",{}))!={item[0] for item in CONFIGS} or
        len(src.get("muon_angle_bin_edges_deg",[]))!=11 or
        len(src.get("electron_phi_bin_edges_deg",[]))!=37):
        raise ValueError("Unexpected non-Jpsi source/statistics/geometry")
    return n


def hist(r):
    h=np.asarray(r["phi_hist"],dtype=np.float64)
    if len(h)!=36 or h.sum()!=r["phi_n"]:
        raise ValueError("Garfield phi endpoint-count conservation failed")
    return h


def setup(ax,title,maxpct):
    r0=maxpct*.42
    angles=np.linspace(-np.pi,np.pi,361)
    ax.fill_between(angles,0,r0,facecolor="lightgray",edgecolor="dimgray",
                    linewidth=1.2,zorder=4)
    ax.set_theta_zero_location("E");ax.set_theta_direction(1)
    ax.set_ylim(0,r0+maxpct*1.18)
    ax.set_thetagrids(range(0,360,45),
        labels=["0° (+u)","45°","90° (+v)","135°",
                "180° (−u)","225°","270° (−v)","315°"],fontsize=8)
    ax.set_yticklabels([]);ax.grid(alpha=.18)
    ax.text(0,0,"ANODE",ha="center",va="center",fontsize=8,zorder=5)
    ax.set_title(title,fontsize=11,pad=19)
    return r0


def l(p,fw):
    return f"{p:.1f} mm · "+("field 75 µm" if fw else "no field wire")


def summary(src,level,out,nparents):
    arr=[]
    for cfg,p,fw in CONFIGS:
        row=src["data"][cfg][level];h=hist(row)
        if not h.sum():raise ValueError("No electron azimuths in "+cfg)
        arr.append((p,fw,h/h.sum()*100,int(row["phi_n"])))
    peak=max(v.max() for _,_,v,_ in arr)
    fig,axes=plt.subplots(2,3,subplot_kw={"projection":"polar"},
                          figsize=(15.5,10))
    theta=np.deg2rad(np.arange(-175,180,10))
    for ax,(p,fw,v,n) in zip(axes.flat,arr):
        r0=setup(ax,l(p,fw)+f"\nN={n:,} anode electrons",peak)
        ax.bar(theta,v,bottom=r0,width=np.deg2rad(9.5),
               align="center",zorder=3,alpha=.84)
    title=f"J/ψ · MS{src['station']} · "+(
         "whole station" if level=="station" else f"DET{src['chosen_detector']}")
    fig.suptitle(f"{title} · {nparents:,} GENERATED parents\n"
                 "Near-anode electron endpoint POSITION azimuth φhit, "
                 "six pitch / field-wire geometries",fontsize=14,y=.98)
    fig.text(.5,.017,"10° azimuth bins. Radial bar height = count / all "
      "known-φ anode electrons (per geometry). Common radial scale. "
      "This is not drift-direction azimuth.",ha="center",fontsize=9)
    fig.tight_layout(rect=(0,.037,1,.935),h_pad=3,w_pad=2)
    name=f"01_wraparound_phi_six_{level}.png"
    fig.savefig(out/name,dpi=200,bbox_inches="tight");plt.close(fig)
    return name


def conditional(src,level,axis,out,nparents):
    limits=src["muon_angle_bin_edges_deg"]
    series=[]
    for cfg,p,fw in CONFIGS:
        grouped=src["data"][cfg][f"theta_{axis}_{level}"]
        subset=[]
        for j in BINS:
            v=grouped.get(str(j))
            if v and v["phi_n"]>=200:
                h=hist(v)
                subset.append((j,100*h/h.sum(),v["phi_n"]))
        series.append((p,fw,subset))
    peak=max([float(v.max()) for _,_,subset in series for _,v,_ in subset]+[0.1])
    fig,axes=plt.subplots(2,3,subplot_kw={"projection":"polar"},
                          figsize=(15.5,10))
    theta=np.deg2rad(np.arange(-175,180,10))
    for ax,(p,fw,subset) in zip(axes.flat,series):
        r0=setup(ax,l(p,fw),peak)
        for j,v,n in subset:
            ax.plot(np.r_[theta,theta[0]+2*np.pi],
                    np.r_[r0+v,r0+v[0]],lw=1.7,
                    label=f"θ{axis} [{limits[j]:g},{limits[j+1]:g})°")
    title=f"J/ψ · MS{src['station']} · "+(
          "whole station" if level=="station" else f"DET{src['chosen_detector']}")
    fig.suptitle(f"{title} · {nparents:,} generated parents\n"
         f"Near-anode endpoint P(φhit | θμ,{axis}) · six geometries · "
         "each entrance-angle bin normalised separately",fontsize=14,y=.98)
    h,lbl=axes.flat[0].get_legend_handles_labels()
    if h:fig.legend(h,lbl,loc="lower center",bbox_to_anchor=(.5,.01),
                    ncol=len(h),fontsize=9,frameon=False)
    fig.tight_layout(rect=(0,.075,1,.935),h_pad=3,w_pad=2)
    name=f"02_wraparound_phi_muon_theta_{axis}_{level}.png"
    fig.savefig(out/name,dpi=190,bbox_inches="tight");plt.close(fig)
    return name


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input",type=Path,required=True)
    parser.add_argument("--out",type=Path,required=True)
    a=parser.parse_args()
    src=json.loads(a.input.read_text(encoding="utf-8"))
    n=provenance(src)
    a.out.mkdir(parents=True,exist_ok=True)
    images=[]
    comparisons=[]
    for level in ("station","chamber"):
        images.append(summary(src,level,a.out,n))
        for axis in ("u","w"):
            images.append(conditional(src,level,axis,a.out,n))
    for pitch in (3.,3.5,4.):
        p0=next(cfg for cfg,p,fw in CONFIGS if p==pitch and fw==0)
        pf=next(cfg for cfg,p,fw in CONFIGS if p==pitch and fw==75)
        h0=hist(src["data"][p0]["station"])
        hf=hist(src["data"][pf]["station"])
        d=np.abs(h0/h0.sum()-hf/hf.sum()).sum()/2
        comparisons.append({
          "pitch_mm":pitch,
          "total_variation_FW_vs_noFW":d,
          "n_electrons_noFW":int(h0.sum()),"n_electrons_withFW":int(hf.sum())})
    with (a.out/"phi_hit_FW_comparison.csv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=comparisons[0].keys())
        w.writeheader();w.writerows(comparisons)
    # Explicit denominator and muon-crossing exposures for every plotted curve.
    bounds=src["muon_angle_bin_edges_deg"]
    with (a.out/"conditional_angle_statistics.csv").open(
            "w",newline="",encoding="utf-8") as output:
        w=csv.writer(output)
        w.writerow(("scope","config","pitch_mm","field_wire_um",
                    "muon_angle_axis","theta_bin_index",
                    "theta_low_deg","theta_high_deg",
                    "muon_gas_crossings","anode_electrons_known_phi"))
        for level in ("station","chamber"):
            for cfg,p,fw in CONFIGS:
                for axis in ("u","w"):
                    rows=src["data"][cfg][f"theta_{axis}_{level}"]
                    for j in range(len(bounds)-1):
                        row=rows.get(str(j),{})
                        if row:
                            hist(row)
                        w.writerow((level,cfg,p,fw,axis,j,bounds[j],bounds[j+1],
                                    row.get("segments",0),row.get("phi_n",0)))
    (a.out/"README.txt").write_text(
      f"Jpsi, MS{src['station']}, DET{src['chosen_detector']}, "
      f"{n:,} generated parents.\\n"
      "φ_hit is endpoint position around wire circumference from Garfield, "
      "NOT arrival velocity direction.\\n"
      "Conditional muon θ_u and θ_w use precomputed 20-trial Garfield LUT. "
      "Per-wire conditional distributions are NOT saved in compact source; "
      "chamber and station PDFs are.\\n"
      "Each curve has its own collected-electron denominator; no naive "
      "independent-electron uncertainty bands.\\n",
      encoding="utf-8")
    print(f"PASS J/ψ {n:,} parents, six geometries, "
          f"{len(images)} real circular plots in {a.out}")
    for r in comparisons:print("pitch",r["pitch_mm"],
          "FW/noFW azimuth TV distance",f"{r['total_variation_FW_vs_noFW']:.3f}")


if __name__=="__main__":
    main()
