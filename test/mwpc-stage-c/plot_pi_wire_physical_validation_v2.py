#!/usr/bin/env python3
"""Physically interpretable PI chamber/anode validation atlas V2.

Every channel/station panel compares:
  1) actual XY gas crossing midpoints, retaining both spatial dimensions
  2) geometric anode-cell intersections, using the SAME fixed physical wires
  3) Garfield model >=1 collected primary e- on anode
  4) Garfield model >=5 collected primary e-
  5) >=1 collected primary e- by 100 ns
  6) collected minus geometric anode activation, same physical wire index.

Wire maps are deliberately constant along X: anode-wire occupancy alone
cannot reconstruct the muon's longitudinal position along the anode!
Color in station wire maps = hits per 10k ALL GAS CHAMBER CROSSINGS in
that station/source; this avoids spurious "bright" rare chambers from
their unreliable per-module conditional fraction.
Per-chamber plots retain conditional rates and compare to geometry.
"""
from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize, TwoSlopeNorm, PowerNorm
from matplotlib.patches import Rectangle

from mwpc_checkerboard_projection import (
    SOURCES, NX, NY, GAS_X_CM, GAS_Y_CM, MS0_GAS_Y_CM,
    build_chambers,
)
from pi_fixed_wire_occupancy import n_wires, positions
from plot_pi_fixed_wire_atlas import load_wire_csv

CHAMBERS=build_chambers()
BY_DET={ch.det_id: ch for station in CHAMBERS for ch in station}
CONDITION_LABELS=(
    "Geometry: wire cells crossed",
    "Garfield: >=1 collected electron",
    "Garfield: >=5 collected electrons",
    "Garfield: >=1 electron by 100 ns",
)
CHANNEL_LABELS={"jpsi40-bore":r"J/$\psi$",
                "phi40-bore":r"$\phi$",
                "omega40-bore":r"$\omega$"}


def load_geometric(path,exposure):
    result=np.zeros((3,298,128),dtype=np.int64)
    seen=np.zeros((3,298,128),dtype=bool)
    with path.open(newline="",encoding="utf-8") as f:
        for row in csv.DictReader(f):
            isrc=SOURCES.index(row["source"])
            det=int(row["det_id"])
            j=int(row["wire_index"])
            if det not in BY_DET or int(row["station"])!=BY_DET[det].station:
                raise ValueError("Incorrect geometry control chamber label")
            if not (0<=j<n_wires(BY_DET[det].station)):
                raise ValueError("Incorrect nominal wire index")
            n=int(row["N_crossings"])
            if n!=exposure[isrc,det]:
                raise ValueError(
                    f"Geometry/response crossing exposure differs for source "
                    f"{SOURCES[isrc]}, DET{det}: {n} vs {exposure[isrc,det]}"
                )
            if seen[isrc,det,j]:raise ValueError("Duplicate geometry anode")
            seen[isrc,det,j]=True
            count=int(row["n_geom"])
            if not 0<=count<=n:raise ValueError("Bad geometric wire count")
            result[isrc,det,j]=count
    for i in range(3):
        for det,ch in BY_DET.items():
            if not np.all(seen[i,det,:n_wires(ch.station)]):
                raise ValueError(f"Missing geometry anode source={i} DET{det}")
    return result


def chamber_rect(ax,ch,alpha=1):
    gy=MS0_GAS_Y_CM if ch.station==0 else GAS_Y_CM
    ax.add_patch(Rectangle(
        (ch.x-GAS_X_CM/2.,ch.y-gy/2.),GAS_X_CM,gy,
        fill=False,lw=.65,edgecolor="#697586",alpha=alpha,zorder=5,
    ))


def physical_bounds(st):
    chamber=CHAMBERS[st]
    gy=MS0_GAS_Y_CM if st==0 else GAS_Y_CM
    xlo=min(c.x-GAS_X_CM/2 for c in chamber)
    xhi=max(c.x+GAS_X_CM/2 for c in chamber)
    ylo=min(c.y-gy/2 for c in chamber)
    yhi=max(c.y+gy/2 for c in chamber)
    # Same station-specific frame for all three channels.
    mx=max(4.,.04*(xhi-xlo))
    my=max(4.,.04*(yhi-ylo))
    return (xlo-mx,xhi+mx,ylo-my,yhi+my)


def stats_over_station(exposure,i,st):
    return int(sum(exposure[i,c.det_id] for c in CHAMBERS[st]))


def make_wire_lines(i,st,counts_per_wire,per_station_crossings,group=4):
    segs, vals=[],[]
    denom=max(per_station_crossings,1)
    for ch in CHAMBERS[st]:
        n=n_wires(st)
        wire_positions=positions(st)
        v=counts_per_wire[i,ch.det_id,:n]
        for start in range(0,n,group):
            last=min(start+group,n)
            # Rate averaged over neighboring wires, not sum.
            rate=10000.*float(np.mean(v[start:last]))/denom
            yy=ch.y+float(np.mean(wire_positions[start:last]))/10.
            # Full wire geometry has no X position information along it.
            segs.append(((ch.x-GAS_X_CM/2.+1.,yy),
                         (ch.x+GAS_X_CM/2.-1.,yy)))
            vals.append(rate)
    return segs,np.asarray(vals)


def station_figure(folder,src_index,st,exposure,geom,wire_counts,xy,
                   normal_max,difference_max):
    src=SOURCES[src_index]
    nstation=stats_over_station(exposure,src_index,st)
    fig,axes=plt.subplots(2,3,figsize=(20,11.8),layout="constrained")
    xy_key=f"{src}_MS{st}"
    image_grid=xy[xy_key]
    xe,ye=xy[xy_key+"_xedges"],xy[xy_key+"_yedges"]
    ax=axes[0,0]
    image=ax.imshow(
        image_grid,origin="lower",
        extent=(xe[0],xe[-1],ye[0],ye[-1]),
        interpolation="nearest",aspect="equal",cmap="Greys",
        norm=PowerNorm(gamma=.5,vmin=0,vmax=max(1,int(np.quantile(image_grid[image_grid>0],.995)))
                      if np.any(image_grid>0) else 1),
    )
    fig.colorbar(image,ax=ax,shrink=.7,label="Crossings / 25 cm$^2$")
    ax.set_title("Muon XY crossings (actual gas-box midpoints)",fontsize=11)
    wire_data=[
        geom,
        wire_counts[:,:,:,0],
        wire_counts[:,:,:,1],
        wire_counts[:,:,:,2],
        wire_counts[:,:,:,0]-geom,
    ]
    norm=Normalize(0.,normal_max,clip=False)
    signed=TwoSlopeNorm(vcenter=0,vmin=-difference_max,vmax=difference_max)
    headings=[
        "Geometric anode cells touched",
        "Garfield anode collection (>=1e)",
        "Garfield anode collection (>=5e)",
        "Garfield anode collection (by 100 ns)",
        "Garfield >=1e minus geometric",
    ]
    for j,metric in enumerate(wire_data):
        ax=axes.flat[j+1]
        segments, values=make_wire_lines(src_index,st,metric,nstation)
        diverging=j==4
        # Every mappable has an independent normalization to avoid old
        # distro matplotlib Colorbar weakref callback issues.
        localnorm=(TwoSlopeNorm(vcenter=0,vmin=-difference_max,vmax=difference_max)
                   if diverging else Normalize(0.,normal_max,clip=False))
        cmap="coolwarm" if diverging else "inferno"
        collection=LineCollection(segments,cmap=cmap,norm=localnorm,
                                   linewidths=3.1,zorder=3)
        collection.set_array(values)
        ax.add_collection(collection)
        fig.colorbar(collection,ax=ax,shrink=.7,
                     extend="both" if diverging else "max",
                     label="Activations / 10k station chamber crossings")
        ax.set_title(headings[j],fontsize=11)
    bounds=physical_bounds(st)
    for ax in axes.flat:
        for ch in CHAMBERS[st]:chamber_rect(ax,ch)
        ax.set_xlim(bounds[:2]);ax.set_ylim(bounds[2:])
        ax.set_aspect("equal",adjustable="box")
        ax.set_xlabel("Global X [cm] (along wires)")
        ax.set_ylabel("Global Y [cm] (across wires)")
        ax.grid(False)
    fig.suptitle(
        f"{CHANNEL_LABELS[src]} | MS{st} | {nstation:,} gas-module crossings | "
        f"{len(CHAMBERS[st])} chambers\n"
        "Muon XY footprint vs fixed physical-anode activity. "
        "Anode X is not resolved. Four wires per drawn stripe.",
        fontsize=14,
    )
    fig.text(.5,.008,
             "Crossing density uses gas-box midpoints. "
             "Anode activity uses collected PRIMARY electrons (not front-end hits). "
             "Contours are nominal gas extents; no final edge-field simulation.",
             ha="center",fontsize=9,color="#555555")
    path=folder/f"physical_validation_{src}_MS{st}.png"
    fig.savefig(path,dpi=135)
    plt.close(fig)


def detailed_figure(folder,source_index,ch,exp,geom,wire_counts):
    src=SOURCES[source_index]
    det=ch.det_id
    n=n_wires(ch.station)
    exposure=int(exp[source_index,det])
    g=geom[source_index,det,:n]
    a=wire_counts[source_index,det,:n,0]
    b=wire_counts[source_index,det,:n,1]
    t=wire_counts[source_index,det,:n,2]
    if exposure:
        p=lambda q:q.astype(float)/exposure
        pg,pa,pb,pt=[p(x) for x in (g,a,b,t)]
    else:
        pg=pa=pb=pt=np.full(n,np.nan)
    fig,axes=plt.subplots(3,1,figsize=(9.2,8.4),sharex=True,layout="constrained")
    idx=np.arange(n)
    labels=[
        "1. Validation: geometrically crossed wire cell vs collected >=1e",
        "2. Threshold sensitivity: collected >=1e vs >=5e",
        "3. Timing sensitivity: all collected >=1e vs by 100 ns",
    ]
    series=[[(pg,"Geometry","#727e89"),(pa,"Garfield >=1e","#c45813")],
            [(pa,"Garfield >=1e","#c45813"),(pb,"Garfield >=5e","#2674a5")],
            [(pa,"All collected >=1e","#c45813"),(pt,"By 100 ns","#257b63")]]
    for k,ax in enumerate(axes):
        if exposure:
            for yy,label,col in series[k]:
                ax.step(idx,yy,where="mid",label=label,color=col,lw=1.25)
            ax.set_ylim(0,max(.015,float(np.nanmax([x for x,_,_ in series[k]]))*1.14))
            ax.legend(loc="upper right",fontsize=8)
        else:
            ax.text(.5,.5,"NO CHAMBER CROSSINGS",
                    transform=ax.transAxes,ha="center")
            ax.set_ylim(0,1)
        ax.set_ylabel("P per chamber crossing")
        ax.set_title(labels[k],fontsize=11)
        ax.grid(alpha=.18)
    axes[-1].set_xlabel(
        f"Fixed physical wire index along global Y (0..{n-1}; 4 mm pitch)"
    )
    label=(f"{CHANNEL_LABELS[src]} | MS{ch.station} | DET{det:03d} "
           f"(row {ch.row}, col {ch.col}, layer {ch.layer}) | N={exposure:,}")
    if 0<exposure<100:
        label+=" | LOW STATISTICS — interpret cautiously"
    fig.suptitle(label+"\nNominal mechanical wire grid; primary-electron response only",
                 fontsize=12)
    path=folder/src/f"MS{ch.station}"/f"wire_DET{det:03d}.png"
    path.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(path,dpi=125)
    plt.close(fig)


def contact_figure(folder,i,st,exposure,geom,wire_counts):
    src=SOURCES[i]
    nx,ny=NX[st],NY[st]
    fig,axes=plt.subplots(ny,nx,figsize=(nx*2.9,ny*1.9),
                          squeeze=False,layout="constrained")
    for ch in CHAMBERS[st]:
        ax=axes[ny-1-ch.row,ch.col]
        det=ch.det_id
        n=n_wires(st)
        g=geom[i,det,:n]
        a=wire_counts[i,det,:n,0]
        b=wire_counts[i,det,:n,1]
        hits=int(exposure[i,det])
        if hits:
            ax.fill_between(np.arange(n),g/hits,step="mid",
                            color="#8892a2",alpha=.35,lw=0)
            ax.plot(g/hits,color="#8892a2",lw=.65)
            ax.plot(a/hits,color="#c45813",lw=.85)
            ax.plot(b/hits,color="#2674a5",lw=.65)
            ax.set_ylim(0,max(.02,max(np.max(g),np.max(a))*1.15/hits))
        else:
            ax.set_facecolor("#eaeaea")
            ax.text(.5,.5,"NO HITS",ha="center",va="center",
                    transform=ax.transAxes,fontsize=7,color="#555")
        if hits<100:
            for spine in ax.spines.values():
                spine.set_color("#be501f")
                spine.set_linewidth(1.6)
        ax.set_title(f"{det}:{ch.layer}  N={hits}",fontsize=7)
        ax.set_xlim(0,n-1)
        ax.set_xticks([])
        ax.set_yticks([])
    fig.suptitle(
        f"{CHANNEL_LABELS[src]} MS{st}: all physical chamber wire profiles\n"
        "gray = geometric cell crossing; orange = Garfield >=1e; "
        "blue = Garfield >=5e; red rim = fewer than 100 crossings",
        fontsize=14,
    )
    path=folder/f"per_chamber_profiles_{src}_MS{st}.png"
    fig.savefig(path,dpi=120)
    plt.close(fig)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--input",type=Path,
        default=Path("test_runs/mwpc_stage_c/pi_fixed_wire_maps"),
        help="Completed wire folding output (no new Garfield or Heed jobs)",
    )
    ap.add_argument(
        "--controls",type=Path,default=None,
        help="Geometric per-wire control and XY grid; defaults to input/physical_validation",
    )
    ap.add_argument("--no-details",action="store_true")
    args=ap.parse_args()
    controls=args.controls or args.input/"physical_validation"
    arr,exposure,total=load_wire_csv(args.input/"pi_fixed_wire_counts.csv")
    geom=load_geometric(controls/"pi_geom_wire_counts.csv",exposure)
    # integer counts from CSV aggregator? load_wire_csv returns per-det SUM;
    # reconstruct exact per-wire integer counts from probabilities times N.
    collected=np.zeros((3,298,128,3),dtype=np.int64)
    for i in range(3):
        for det,ch in BY_DET.items():
            n=int(exposure[i,det])
            if n:
                vals=arr[i,det,:n_wires(ch.station),:]
                collected[i,det,:n_wires(ch.station),:]=np.rint(vals*n).astype(np.int64)
    xy=np.load(controls/"pi_crossing_density_xy.npz")
    output=controls/"atlas_v2"
    for item in ("station_validation","chamber_profiles","contact_sheets"):
        (output/item).mkdir(parents=True,exist_ok=True)
    # Shared per-channel color ceiling in activations per 10k station chamber
    # crossings; station counts are NOT conditional module probabilities.
    for i,src in enumerate(SOURCES):
        all_rates=[]
        all_differences=[]
        for st in range(6):
            den=max(1,stats_over_station(exposure,i,st))
            for ch in CHAMBERS[st]:
                n=n_wires(st);d=ch.det_id
                all_rates.append((10000./den)*geom[i,d,:n])
                for k in range(3):
                    all_rates.append((10000./den)*collected[i,d,:n,k])
                all_differences.append(
                    (10000./den)*(collected[i,d,:n,0]-geom[i,d,:n])
                )
        positives=np.concatenate(all_rates)
        positives=positives[positives>0]
        normal_max=float(np.quantile(positives,.998)) if len(positives) else 1.
        normal_max=max(normal_max,.01)
        difference_max=max(
            float(np.quantile(np.abs(np.concatenate(all_differences)),.998)),
            .01,
        )
        print(
            f"{src}: one common six-station rate scale "
            f"max={normal_max:.3g} /10k station crossings; "
            f"|difference| max={difference_max:.3g}",
            flush=True,
        )
        for st in range(6):
            station_figure(output/"station_validation",i,st,exposure,
                           geom,collected,xy,normal_max,difference_max)
            contact_figure(output/"contact_sheets",i,st,exposure,geom,collected)
            if not args.no_details:
                for ch in CHAMBERS[st]:
                    detailed_figure(output/"chamber_profiles",i,ch,
                                    exposure,geom,collected)
            print(f"  {src} MS{st} completed",flush=True)
    print(
        f"PI physical validation V2 finished. "
        f"18 six-panel station figures; 18 geometric/Garfield contacts; "
        f"{0 if args.no_details else 894} geometric/Garfield chamber profiles.",
        flush=True,
    )
    print(output,flush=True)


if __name__=="__main__":
    main()
