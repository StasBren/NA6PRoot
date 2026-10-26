#!/usr/bin/env python3
"""PI question 2: eventwise up/down collection asymmetry using p=3mm LUT500.

READ ONLY on original Geant4 crossings, original Heed library, immutable
Garfield microscopic LUT500. Stream native J/psi MS0 gas crossings; sample
the EXACT original stable-seeded Poisson-Heed births for each crossing.
Apply same birth positions / lookup quantiles to each 500-trial geometry.
Never run new Geant4 or Garfield physics or change frozen 720-shard folds.

Two distinct observables:
  * parent-event A: pool collected electrons from ALL MS0 muon crossings
    belonging to the same (source, chunk, local parent event).
  * crossing A: one native chamber gas crossing, with its own signed Geant4
    muon entrance angles theta_u / theta_w for the angle-dependence study.
The parent-event distribution and crossing-angle RMS are NOT the same
statistical population. A is undefined when no classified electrons exist.

"Up" means v>0 on the anode circumference, inferred from the sign of
Garfield phi_hit=atan2(v_end,u_end-u_anode) in (-180,180). Exactly 0/+/-180
degree hits are classified as boundary and excluded from A denominators;
edge-lost and non-anode-collected electrons are excluded too. No gain or
electronics are simulated. Finite spatial LUT, not independent fresh Garfield
physics per Geant4 electron, still limits predictive precision.

Outputs: new isolated JSON and two PNG summaries. Default small pilot.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import time
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from pi_2m_native_six_geometry_fold import crossing_provenance
from pi_2m_zero_path_compatible_resume import coincident_zero_chord
from pi_six_geometry_native_timing import (
    HEED_TRACKS, HEED_CLUSTERS, birth_sample, load_map, n_anodes,
    response_for_births, stable_seed, ANGLES,
)
from pi_lut_20_vs_500_convergence import ROOT500
from pi_pitch_fieldwire_pilot import tag
from run_pi_per_chamber_response import read_heed
from run_pi_geant4_campaign import REPO

OUT=REPO/"test_runs/mwpc_stage_c/pi_event_asymmetry_500_ms0_v1"
SOURCE="jpsi40-bore"
PITCH=3.
CONFIGS=(0.,75.)
N_PARENTS_CHUNK=50_000
A_EDGES=np.linspace(-1.,1.,42)  # 41 bins: include +/-1 in endpoint bins
VERSION="PI_NATIVE_JPSI_MS0_PARENT_AND_CROSSING_ASYM_P3_LUT500_V1"


def up_down(samples, finite):
    """Three-way exclusive split of accepted Garfield anode hits."""
    if len(samples)!=len(finite):
        raise ValueError("Response and collection flags differ")
    if not np.any(finite):
        return 0,0,0
    phi=samples[finite,4]
    if not np.all(np.isfinite(phi)&(phi>=-180.)&(phi<=180.)):
        raise ValueError("Collected electron has undefined Garfield phi_hit")
    # In the current B=0 anode geometry v_end sign is sign(sin(phi)).
    # phi=0 or +/-180 is a geometric boundary, not 'up' or 'down'.
    boundary=np.isclose(phi,0.,atol=1e-10,rtol=0.) | \
             np.isclose(np.abs(phi),180.,atol=1e-10,rtol=0.)
    up=int(np.count_nonzero((phi>0.)&~boundary))
    down=int(np.count_nonzero((phi<0.)&~boundary))
    n_boundary=int(boundary.sum())
    if up+down+n_boundary!=int(np.count_nonzero(finite)):
        raise AssertionError("Electron azimuthal hemisphere conservation")
    return up,down,n_boundary


def asymmetry(up,down):
    """Finite-defined A, otherwise None; do not replace zero counts by A=0."""
    if up<0 or down<0:
        raise ValueError("Negative electron count")
    return (up-down)/(up+down) if up+down else None


@dataclass
class Moments:
    hits: int=0
    valid: int=0
    empty: int=0
    up: int=0
    down: int=0
    sum_a: float=0.
    sum_a2: float=0.
    sum_n: int=0

    def add(self,u,d):
        self.hits+=1
        self.up+=u
        self.down+=d
        a=asymmetry(u,d)
        if a is None:
            self.empty+=1
            return
        self.valid+=1
        self.sum_a+=a
        self.sum_a2+=a*a
        self.sum_n+=u+d

    def result(self):
        if self.hits!=self.valid+self.empty:
            raise AssertionError("Crossing exposure conservation failure")
        if not self.valid:
            return dict(hits=self.hits,valid=0,empty=self.empty,
                        mean=None,rms=None,mean_collected=None,
                        electron_weighted_asymmetry=None)
        mean=self.sum_a/self.valid
        return dict(hits=self.hits,valid=self.valid,empty=self.empty,
            mean=mean,
            rms=math.sqrt(max(0.,self.sum_a2/self.valid-mean**2)),
            mean_collected=self.sum_n/self.valid,
            electron_weighted_asymmetry=((self.up-self.down)/
                (self.up+self.down)) if self.up+self.down else None)


def summarize_parent_counts(all_chunks):
    """Counts from each source+chunk+local event, with no fake zero-A events."""
    vals=[]
    n_with_crossing=0
    n_with_classified=0
    total_parent_slots=0
    up=down=0
    for item in all_chunks:
        pair=item["pair"]
        valid=(pair[:,0]+pair[:,1])>0
        total_parent_slots+=len(pair)
        n_with_crossing+=int(np.count_nonzero(item["has_crossing"]))
        n_with_classified+=int(np.count_nonzero(valid))
        up+=int(pair[:,0].sum())
        down+=int(pair[:,1].sum())
        v=((pair[valid,0]-pair[valid,1])/
           (pair[valid,0]+pair[valid,1]))
        vals.append(v)
    a=np.concatenate(vals) if vals else np.empty(0,dtype=float)
    if not len(a):
        return dict(parents_generated=total_parent_slots,
                    parents_with_MS0_crossing=n_with_crossing,
                    parents_with_classified_electrons=0,
                    n_with_undefined_asymmetry=total_parent_slots,
                    mean=None,rms=None,hist=[0]*41,
                    histogram_edges=A_EDGES.tolist(),
                    total_up=up,total_down=down)
    hist,_=np.histogram(a,bins=A_EDGES)
    if hist.sum()!=len(a) or np.any(a < -1) or np.any(a > 1):
        raise AssertionError("Parent-event A histogram conservation")
    return dict(parents_generated=total_parent_slots,
        parents_with_MS0_crossing=n_with_crossing,
        parents_with_classified_electrons=n_with_classified,
        n_with_undefined_asymmetry=total_parent_slots-n_with_classified,
        mean=float(np.mean(a)),rms=float(np.std(a,ddof=0)),
        median=float(np.median(a)),
        histogram_counts=hist.tolist(),
        histogram_edges=A_EDGES.tolist(),total_up=up,total_down=down,
        electron_weighted_asymmetry=(up-down)/(up+down))


def plot_results(result, output):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise RuntimeError("Install matplotlib in the slide venv") from exc

    fig,ax=plt.subplots(figsize=(10.2,5.5),dpi=155)
    edges=np.asarray(A_EDGES)
    for key,color,label in (
        (tag(PITCH,0.),"#404040","No field wire"),
        (tag(PITCH,75.),"#176f9b","Field wire 75 μm")):
        row=result["parent_level"][key]
        counts=np.asarray(row.get("histogram_counts",row.get("hist",[])),
                          dtype=float)
        if len(counts)!=41:
            raise ValueError("No parent-level A histogram")
        density=counts/(max(1,counts.sum())*np.diff(edges))
        ax.stairs(density,edges,color=color,label=label,linewidth=1.8)
    ax.set_xlabel(r"Parent-event asymmetry $A=(N_{\rm up}-N_{\rm down})/(N_{\rm up}+N_{\rm down})$")
    ax.set_ylabel("Probability density")
    ax.set_title("J/ψ, MS0: parent-event electron collection asymmetry (LUT500)")
    ax.legend(frameon=False)
    ax.grid(alpha=.18)
    fig.tight_layout()
    fig.savefig(output/"parent_event_A_distribution.png")
    plt.close(fig)

    fig,axes=plt.subplots(1,2,figsize=(14.5,5.1),dpi=155,sharey=True)
    for ax,axis in zip(axes,("u","w")):
        edges=np.asarray(ANGLES,dtype=float)
        for key,color,label in (
            (tag(PITCH,0.),"#404040","No field wire"),
            (tag(PITCH,75.),"#176f9b","Field wire 75 μm")):
            rows=result["crossing_by_angle"][key][axis]
            xx=[j for j,r in enumerate(rows) if r["rms"] is not None
                and r["valid"]>=30]
            yy=[rows[j]["rms"] for j in xx]
            ax.plot(xx,yy,marker="o",markersize=4,linewidth=1.5,
                    color=color,label=label)
        ax.set_xticks(range(len(edges)-1))
        ax.set_xticklabels(
            [f"{edges[i]:g}..{edges[i+1]:g}" for i in range(len(edges)-1)],
            rotation=55,ha="right",fontsize=8)
        ax.set_xlabel(r"Muon signed entrance $\theta_"+axis+r"$ bin [deg]")
        ax.set_title(r"Crossing-level RMS(A) vs $\theta_"+axis+"$")
        ax.grid(alpha=.2)
    axes[0].set_ylabel("RMS(A) (crossing-level)")
    axes[0].legend(frameon=False)
    fig.suptitle("Muon-angle dependence (individual gas crossings; RMS shown only for N >= 30)")
    fig.tight_layout()
    fig.savefig(output/"crossing_RMS_vs_muon_angles.png")
    plt.close(fig)


def run(chunks,max_hits,out,make_plots=True):
    if not 1<=chunks<=40 or max_hits<0:
        raise ValueError("--chunks 1..40 and --max-hits >= 0 required")
    t0=time.perf_counter()
    lut={}
    lut_meta={}
    for diameter in CONFIGS:
        name=tag(PITCH,diameter)
        lut[name],lut_meta[name]=load_map(PITCH,diameter,ROOT500)
        if lut[name].shape!=(900,500,5):
            raise ValueError("Expected complete 900-cell, 500-drift Garfield LUT")
    lam,marks=read_heed(HEED_TRACKS,HEED_CLUSTERS)
    crossing={key:{axis:[Moments() for _ in range(len(ANGLES)-1)]
                   for axis in ("u","w")} for key in lut}
    parent={key:[] for key in lut}
    counts={key:dict(primary=0,anode_collected=0,edge_collected=0,
                     up=0,down=0,boundary=0)
            for key in lut}
    processed=zero_paths=0
    crossing_prov=[]
    for chunk in range(1,chunks+1):
        path,meta=crossing_provenance("Jpsi",chunk)
        crossing_prov.append(meta["csv_sha256"])
        has_crossing=np.zeros(N_PARENTS_CHUNK,dtype=bool)
        event_counts={key:np.zeros((N_PARENTS_CHUNK,2),dtype=np.int64)
                      for key in lut}
        n_chunk=0
        with path.open(newline="",encoding="utf-8") as handle:
            for row in csv.DictReader(handle):
                if int(row["station"])!=0 or abs(int(row["pdg"]))!=13:
                    continue
                if max_hits and processed>=max_hits:
                    break
                event=int(row["event"])
                if not 0<=event<N_PARENTS_CHUNK:
                    raise ValueError("Bad local parent index in native Geant4")
                has_crossing[event]=True
                length=float(row["path_length_mm"])
                if length==0.:
                    if not coincident_zero_chord(row):
                        raise ValueError("Noncompliant exact-zero gas path")
                    zero_paths+=1
                    positions=np.empty((0,3),dtype=float)
                    q=np.empty(0,dtype=float)
                else:
                    seed=stable_seed("native-six-pitch-shared-Heed-v1",
                        SOURCE,f"chunk{chunk:04d}",row["event"],
                        row["hit_index"],row["track_id"])
                    positions,q=birth_sample(row,lam,marks,seed)
                for name,pool in lut.items():
                    finite,rawcol,wires,samples=response_for_births(
                        positions[:,0],positions[:,1],q,PITCH,0,pool)
                    u,d,b=up_down(samples,finite)
                    nfinite=int(np.count_nonzero(finite))
                    nprim=len(q)
                    counts[name]["primary"]+=nprim
                    counts[name]["anode_collected"]+=nfinite
                    counts[name]["edge_collected"]+=int(
                        np.count_nonzero(rawcol&~finite))
                    counts[name]["up"]+=u
                    counts[name]["down"]+=d
                    counts[name]["boundary"]+=b
                    event_counts[name][event,0]+=u
                    event_counts[name][event,1]+=d
                    for axis,field in (("u","theta_u_in_deg"),
                                       ("w","theta_w_in_deg")):
                        theta=float(row[field])
                        if not math.isfinite(theta) or not -180<=theta<=180:
                            raise ValueError("Bad signed native muon incidence")
                        bid=int(np.clip(np.searchsorted(ANGLES,theta,
                                             side="right")-1,0,len(ANGLES)-2))
                        crossing[name][axis][bid].add(u,d)
                n_chunk+=1
                processed+=1
                if processed%10000==0:
                    print(f"Progress: {processed:,} MS0 muon crossings "
                          f"across {chunk} Geant4 chunk(s)",flush=True)
        for key in lut:
            parent[key].append({"pair":event_counts[key],
                                "has_crossing":has_crossing})
        print(f"PASS J/ψ chunk {chunk}/{chunks}: "
              f"{n_chunk:,} MS0 gas crossings, "
              f"cumulative={processed:,}",flush=True)
        if max_hits and processed>=max_hits:
            break
    if not processed:
        raise ValueError("No MS0 crossings selected")
    parent_level={key:summarize_parent_counts(value)
                  for key,value in parent.items()}
    for key in lut:
        c=counts[key]
        if c["up"]+c["down"]+c["boundary"]!=c["anode_collected"]:
            raise AssertionError("Total Garfield anode-hit conservation failure")
        if c["up"]!=parent_level[key]["total_up"] or \
           c["down"]!=parent_level[key]["total_down"]:
            raise AssertionError("Crossing vs parent count conservation failure")
        if sum(m.hits for m in crossing[key]["u"])!=processed:
            raise AssertionError("Muon crossing angle exposure mismatch")
        if sum(m.hits for m in crossing[key]["w"])!=processed:
            raise AssertionError("Muon entrance-angle exposure mismatch")
    result={
        "schema":VERSION,
        "selection":dict(source=SOURCE,station=0,pitch_mm=PITCH,
                         chunks_requested=chunks,max_hits=max_hits,
                         sampled_muon_gas_crossings=processed,
                         zero_path_crossings=zero_paths,
                         partial_exposure=bool(max_hits),
                         parent_identifier="source+chunk+local parent event",
                         crossing_definition="one native chamber gas crossing",
                         entrance_angles="signed theta_u_in_deg and theta_w_in_deg",
                         up_down_definition="Garfield phi_hit_deg sign, "
                             "strictly excluding phi=0,+/-180 and missing "
                             "or finite-edge-lost anode responses",
                         zero_collected_policy="A undefined, exclude from "
                             "P(A)/mean/RMS; retain exposure counts",
                         no_charge_gain=True,
                         original_Heed_seed="native-six-pitch-shared-Heed-v1",
                         measured_wall_seconds=time.perf_counter()-t0),
        "garfield_LUT500":lut_meta,
        "native_csv_sha256":crossing_prov,
        "heed_library":dict(tracks=str(HEED_TRACKS),
                            clusters=str(HEED_CLUSTERS),
                            clusters_per_mm=lam,
                            multiplicity_mean=float(np.mean(marks))),
        "collection_counts":counts,
        "parent_level":parent_level,
        "crossing_by_angle":{
            key:{axis:[m.result() for m in vals]
                 for axis,vals in bucket.items()}
            for key,bucket in crossing.items()},
    }
    out.mkdir(parents=True,exist_ok=True)
    target=out/"asymmetry_summary.json"
    tmp=target.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(result,indent=2,ensure_ascii=False,
                              allow_nan=False)+"\n",encoding="utf-8")
    tmp.replace(target)
    if make_plots:
        plot_results(result,out)
    for key in lut:
        v=parent_level[key]
        mean_text=f"{v['mean']:+.5f}" if v["mean"] is not None else "undefined"
        rms_text=f"{v['rms']:.5f}" if v["rms"] is not None else "undefined"
        print(f"{key}: {v['parents_with_classified_electrons']:,} "
              f"valid parent events, <A>={mean_text}, "
              f"RMS(A)={rms_text}",flush=True)
    print("PASS: independent per-parent and per-crossing LUT500 asymmetry "
          f"diagnostics: {out}",flush=True)
    return result


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--chunks",type=int,default=1)
    ap.add_argument("--max-hits",type=int,default=5000,
                    help="At most this many MS0 crossings; 0 uses all")
    ap.add_argument("--out",type=Path,default=OUT)
    ap.add_argument("--no-plots",action="store_true")
    args=ap.parse_args()
    run(args.chunks,args.max_hits,args.out,not args.no_plots)


if __name__=="__main__":
    main()
