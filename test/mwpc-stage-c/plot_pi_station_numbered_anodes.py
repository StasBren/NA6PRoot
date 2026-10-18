#!/usr/bin/env python3
"""PI-explainable station diagram and EVERY numbered anode of one MWPC module.

The station layout and numbered wire positions are generated from the SAME
Stage-C chamber geometry and fixed-anode definition as the simulations.
The example wire activity is read from the complete per-wire CSV, never faked.

Creates:
  station_module_and_wire_legend_MS<station>_<source>_DET<id>.png
  all_numbered_wires_MS<station>_<source>_DET<id>.png
All physical locations are nominal until exact mechanical wire origins and
final number of wires have been confirmed with the PI.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.patches import Rectangle
import numpy as np

from mwpc_checkerboard_projection import (
    SOURCES, build_chambers, GAS_X_CM, GAS_Y_CM, MS0_GAS_Y_CM, NX, NY,
)
from pi_fixed_wire_occupancy import n_wires, positions


CHAMBERS = build_chambers()
INDEX = {ch.det_id: ch for station in CHAMBERS for ch in station}


def load_one(path, source, det_id, station):
    result = {}
    with path.open(newline="", encoding="utf-8") as file:
        for row in csv.DictReader(file):
            if row["source"] != source or int(row["det_id"]) != det_id:
                continue
            idx = int(row["wire_index"])
            if idx in result:
                raise ValueError(f"Repeated wire #{idx} for DET{det_id}")
            if int(row["station"]) != station:
                raise ValueError("Wrong station for requested chamber")
            result[idx] = row
    nw = n_wires(station)
    if set(result) != set(range(nw)):
        raise ValueError(f"Expected numbered physical wires 0..{nw-1}")
    n = int(result[0]["chamber_crossings"])
    if any(int(row["chamber_crossings"]) != n for row in result.values()):
        raise ValueError("Per-wire denominators differ in one chamber")
    return result, n


def make_station(out, src, chamber, rows, n):
    st = chamber.station
    gy = MS0_GAS_Y_CM if st == 0 else GAS_Y_CM
    group = CHAMBERS[st]
    fig, axs = plt.subplots(
        1, 2, figsize=(20, 11.7), gridspec_kw={"width_ratios": [1.35, 1]},
        layout="constrained",
    )
    ax, zoom = axs
    for ch in group:
        chosen = ch.det_id == chamber.det_id
        ax.add_patch(Rectangle(
            (ch.x-GAS_X_CM/2, ch.y-gy/2), GAS_X_CM, gy,
            facecolor="#f3b985" if chosen else "#e8eef2",
            edgecolor="#ba4b20" if chosen else "#738697",
            linewidth=3.2 if chosen else 1, alpha=.86,
            zorder=4 if chosen else 2,
        ))
        ax.text(
            ch.x, ch.y, f"DET{ch.det_id:03d}\n{ch.layer}",
            ha="center", va="center", fontsize=8,
            fontweight="bold" if chosen else "normal", zorder=5,
        )

    x_max=max(abs(ch.x)+GAS_X_CM/2 for ch in group)
    y_max=max(abs(ch.y)+gy/2 for ch in group)
    ax.set(xlim=(-x_max-7,x_max+7),ylim=(-y_max-7,y_max+7),
           xlabel="Global X [cm]: along anode wires",
           ylabel="Global Y [cm]: across anode wires",
           title=f"MS{st}: {NX[st]} columns x {NY[st]} rows = "
                 f"{len(group)} modules\nOrange = selected DET{chamber.det_id:03d}")
    ax.set_aspect("equal")
    ax.grid(alpha=.1,zorder=0)

    # Enlarged geometry: all wires physically horizontal along X.
    ypos = positions(st)/10.  # mm -> cm, module-local
    wire_rates=np.array(
        [float(rows[j]["P_ge1e"]) for j in range(len(ypos))],dtype=float
    )
    wire_colors=plt.get_cmap("YlOrRd")(
        np.clip(np.nan_to_num(wire_rates) / max(np.nanmax(wire_rates),.001),0,1)
    )
    for j,y in enumerate(ypos):
        zoom.plot([-GAS_X_CM/2,GAS_X_CM/2],[y,y],
                  color=wire_colors[j],lw=1.4,alpha=.9,zorder=2)
    zoom.add_patch(Rectangle(
        (-GAS_X_CM/2,-gy/2),GAS_X_CM,gy,
        fill=False,lw=1.8,edgecolor="#3d5365",
    ))
    samples=sorted(set((0,(len(ypos)-1)//4,(len(ypos)-1)//2,
                        3*(len(ypos)-1)//4,len(ypos)-1)))
    for j in samples:
        y=ypos[j]
        zoom.scatter([-GAS_X_CM/2],[y],s=35,color="#ba4b20",zorder=8)
        zoom.annotate(
            f"wire {j}",xy=(-GAS_X_CM/2,y),
            xytext=(-GAS_X_CM/2-3.5,y),
            textcoords="data",ha="right",va="center",fontsize=10,
            fontweight="bold",color="#943d25",
        )
    zoom.set(
        xlim=(-GAS_X_CM/2-21,GAS_X_CM/2+4),
        ylim=(-gy/2-1,gy/2+1),
        xlabel="Local X [cm] (physical wire length)",
        ylabel="Local Y [cm] (number increases upward)",
        title=f"DET{chamber.det_id:03d} layer {chamber.layer}, "
              f"N={n:,} muon chamber crossings\n"
              f"{len(ypos)} real-indexed nominal wires (4 mm pitch)",
    )
    zoom.set_aspect("equal")
    zoom.grid(alpha=.1)
    fig.suptitle(
        f"NA60+ MWPC Stage C — station MS{st} and module DET{chamber.det_id:03d}\n"
        f"{src}; projected geometry and simulated primary-electron anode activity",
        fontsize=15,
    )
    fig.text(.5,.015,
             "Nominal wire origin and count; full mechanical specification "
             "and electronics threshold NOT confirmed. Colors on wires "
             "illustrate relative primary-electron activity, not X position.",
             ha="center",fontsize=9)
    name=out / (
        f"station_module_and_wire_legend_MS{st}_{src}_DET{chamber.det_id:03d}.png"
    )
    fig.savefig(name,dpi=170)
    plt.close(fig)
    print("Wrote:",name,flush=True)


def make_all_wires(out,src,chamber,rows,n):
    st=chamber.station
    nw=n_wires(st)
    width=GAS_X_CM
    gy=MS0_GAS_Y_CM if st==0 else GAS_Y_CM
    y=positions(st)/10
    p1=np.array([float(rows[j]["P_ge1e"]) for j in range(nw)])
    p5=np.array([float(rows[j]["P_ge5e"]) for j in range(nw)])
    p100=np.array([float(rows[j]["P_le100ns"]) for j in range(nw)])
    fig,ax=plt.subplots(figsize=(15,22),layout="constrained")
    # Enough physical height for each index to be printed legibly.
    for j in range(nw):
        yy=y[j]
        color="#c85b20" if np.isfinite(p1[j]) and p1[j]>.0 else "#8f9ea7"
        ax.plot([-width/2,width/2],[yy,yy],color=color,lw=.85)
        ax.text(-width/2-2, yy, f"{j:03d}",
                fontsize=7.0,va="center",ha="right",color="#1f2d3b",
                fontfamily="monospace")
        if np.isfinite(p1[j]):
            ax.text(
                width/2+2,yy,
                f"{100*p1[j]:5.2f}%   {100*p5[j]:5.2f}%   {100*p100[j]:5.2f}%",
                fontsize=6.5,va="center",ha="left",
                color="#263847",fontfamily="monospace",
            )
        if j in (0,nw//4,(nw-1)//2,3*(nw-1)//4,nw-1):
            ax.scatter([0],[yy],s=13,color="#205b89",zorder=6)
    ax.add_patch(Rectangle(
        (-width/2,-gy/2),width,gy,fill=False,color="#405b6f",lw=1.8,
    ))
    ax.set_xlim(-width/2-12,width/2+62)
    ax.set_ylim(-gy/2-1.0,gy/2+2.8)
    ax.set_aspect("auto")
    ax.grid(False)
    ax.set_xticks([-width/2,0,width/2])
    ax.set_xlabel("Local X [cm] — anode wire extends horizontally")
    ax.set_ylabel("Local Y [cm] — increasing wire index")
    ax.set_title(
        f"Every physical anode of MS{st}, DET{chamber.det_id:03d}, layer {chamber.layer}\n"
        f"{src} | N={n:,} projected chamber crossings | 4-mm pitch\n"
        f"LEFT: wire index (000..{nw-1:03d}); RIGHT: P(>=1 e-)  P(>=5 e-)  "
        f"P(>=1 e- in 100 ns)\n"
        "Right columns are probabilities PER CROSSING OF THIS MODULE",
        fontsize=12,pad=12,
    )
    ax.text(
        .5,.005,
        "Nominal centered fixed wire grid. One muon can activate >1 anode. "
        "No avalanche/electronics gain and no X localization on each wire.",
        ha="center",va="bottom",transform=fig.transFigure,fontsize=10,
    )
    name=out / f"all_numbered_wires_MS{st}_{src}_DET{chamber.det_id:03d}.png"
    fig.savefig(name,dpi=145)
    plt.close(fig)
    print("Wrote:",name,flush=True)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--counts",type=Path,
        default=Path("test_runs/mwpc_stage_c/pi_fixed_wire_maps/pi_fixed_wire_counts.csv"),
    )
    ap.add_argument(
        "--out",type=Path,
        default=Path("test_runs/mwpc_stage_c/pi_fixed_wire_maps/pi_explanation"),
    )
    ap.add_argument("--det-id",type=int,default=17)
    ap.add_argument("--source",choices=SOURCES,default="jpsi40-bore")
    args=ap.parse_args()
    if args.det_id not in INDEX:
        ap.error("--det-id must be an existing 0..297 detector ID")
    ch=INDEX[args.det_id]
    data,n=load_one(args.counts,args.source,ch.det_id,ch.station)
    args.out.mkdir(parents=True,exist_ok=True)
    make_station(args.out,args.source,ch,data,n)
    make_all_wires(args.out,args.source,ch,data,n)


if __name__=="__main__":
    main()
