#!/usr/bin/env python3
"""Read-only joint angle/time and up-down checks for attachment-free fast MC.

Matched birth starts: micro was simulated WITH attachment and fast MC
WITHOUT it. Only the *common-birth anode-survivor* cohort is compared.
This is NOT a pure microscopic attachment-off reference. Disjoint fit
and validation samples are kept, and 5 time-quantile edges are defined
using microscopic FIT only, frozen before validation.
"""
from __future__ import annotations

import argparse
import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np

from pi_fast_lut_existing_starts_report import GEOMS, read_micro, read_mc, sha256

PHI_EDGES=np.linspace(-180.,180.,13)

def phi_time(row,model):
    if model=="micro":
        phi=float(row["reference_phi_hit_deg"])
        time=float(row["reference_drift_time_ns"])
    else:
        phi=float(row["phi_hit_deg"])
        time=float(row["drift_time_ns"])
    if not math.isfinite(phi) or not -180.<=phi<=180.:
        raise ValueError("Invalid common-survivor phi")
    if not math.isfinite(time) or time<0.:
        raise ValueError("Invalid common-survivor drift time")
    return phi,time

def common(pairs):
    return [(r,m) for r,m in pairs
            if int(r["reference_collected"])==1
               and int(m["anode_collected"])==1]

def joint_hist(entries,which,t_edges):
    angles,times=zip(*(phi_time(r if which=="micro" else m,which)
                     for r,m in entries))
    counts=np.histogram2d(angles,times,
                 bins=(PHI_EDGES,t_edges))[0].astype(int)
    if int(counts.sum())!=len(entries):
        raise ValueError("Joint PDF lost collected electron records")
    return counts

def normhist(c):
    x=np.asarray(c,dtype=float)
    if x.sum()==0: raise ValueError("Empty joint PDF")
    return x/x.sum()

def total_variation(x,y):
    return float(.5*np.abs(normhist(x)-normhist(y)).sum())

def angle_hemi(pairs,which):
    vals=[phi_time(r if which=="micro" else m,which)
          for r,m in pairs]
    n=len(vals)
    if not n: raise ValueError("Empty hemisphere sample")
    upper=[t for phi,t in vals if phi>=0.]
    lower=[t for phi,t in vals if phi<0.]
    mean_upper=float(np.mean(upper)) if upper else None
    mean_lower=float(np.mean(lower)) if lower else None
    return {
        "n":n,
        "P_phi_ge_0":len(upper)/n,
        "upper_minus_lower_count_over_n":(len(upper)-len(lower))/n,
        "mean_time_upper_ns":mean_upper,
        "mean_time_lower_ns":mean_lower,
        "mean_time_upper_minus_lower_ns":
          mean_upper-mean_lower if upper and lower else None,
        "p95_time_ns":float(np.quantile([v[1] for v in vals],.95)),
        "p99_time_ns":float(np.quantile([v[1] for v in vals],.99))
    }

def hemi_cluster_bootstrap(pairs,replicates,seed):
    cells=defaultdict(list)
    for r,m in pairs:
        cells[(int(r["iu"]),int(r["iv"]))].append((r,m))
    # Some of 900 cells may have 0 micro survivors; they must be
    # represented as zero counts rather than omitted from bootstrap.
    data=[]
    for _,rows in sorted(cells.items()):
        n=len(rows)
        n_micro=sum(1 for r,_ in rows if phi_time(r,"micro")[0]>=0.)
        n_mc=sum(1 for _,m in rows if phi_time(m,"mc")[0]>=0.)
        data.append((n_micro,n_mc,n))
    a=np.asarray(data,dtype=float)
    rng=np.random.default_rng(seed)
    samples=rng.integers(len(a),size=(replicates,len(a)))
    sums=a[samples].sum(axis=1)
    diff=(sums[:,1]-sums[:,0])/sums[:,2]
    return [float(v) for v in np.quantile(diff,[.025,.975])]

def grouped_stats(pairs):
    out={}
    for key,groupby in (
        ("u_block",lambda r:int(r["iu"])//5),
        ("v_block",lambda r:int(r["iv"])//5),
        ("u_edge_vs_inner",lambda r:
              "edge" if int(r["iu"])//5 in (0,5) else "inner"),
    ):
        out[key]={}
        labels={groupby(r) for r,_ in pairs}
        for label in sorted(labels,key=str):
            selected=[(r,m) for r,m in pairs if groupby(r)==label]
            if len(selected)<10:
                raise ValueError("Too few surviving samples in region")
            micro=angle_hemi(selected,"micro")
            fast=angle_hemi(selected,"mc")
            out[key][str(label)]={
                "n":len(selected),
                "micro_P_upper":micro["P_phi_ge_0"],
                "mc_P_upper":fast["P_phi_ge_0"],
                "delta_P_upper":fast["P_phi_ge_0"]-micro["P_phi_ge_0"],
                "micro_mean_time_upper_minus_lower_ns":
                  micro["mean_time_upper_minus_lower_ns"],
                "mc_mean_time_upper_minus_lower_ns":
                  fast["mean_time_upper_minus_lower_ns"],
            }
    return out

def report(ref_path,mc_path,out,replicates,seed):
    if out.exists():raise FileExistsError("Refusing to overwrite "+str(out))
    reference,provenance=read_micro(ref_path)
    mc,step,attachment=read_mc(mc_path,reference)
    if attachment!="off" or step!=5.:
        raise ValueError("Expected attachment-off 5um fast MC output")
    result={
        "schema":"PI_MWPC_TRANSPORT_JOINT_PHI_TIME_V1",
        "reference_sha256":sha256(ref_path),
        "mc_sha256":sha256(mc_path),
        "original_microscopic_sha256s":{
             k:v["sha256"] for k,v in
             provenance["source_provenance"].items()},
        "selection":"Common birth sites where MICRO and MC both collected. "
               "Microscopic still includes attachment selection.",
        "time_edges_ns_by_geometry_defined_on_micro_fit_only":{},
        "geometries":{},
        "limitations":[
            "The microscopic reference is attachment-on; this is a "
            "survivor-selected transport diagnostic, not validation of "
            "truly attachment-free microscopic trajectories.",
            "12 phi bins x 5 time bins have appreciable finite-sample "
            "noise with ~3000 common electron starts per split.",
            "Small time-tail fluctuations, bin edges, or spatial cells "
            "must not be declared physics failures without uncertainty.",
            "The binning is frozen using reference microscopic FIT; "
            "VALIDATION is not used in defining bin edges.",
            "Sign phi >= 0 is a convenient upper/lower impact-angle "
            "proxy, not an induced strip charge or full event asymmetry."
        ]
    }
    for j,geom in enumerate(GEOMS):
        cohorts={}
        for split in ("fit","validation"):
            pairs=[(r,mc[key]) for key,r in reference.items()
                   if key[0]==geom and r["split"]==split]
            if len(pairs)!=3600:
                raise ValueError("Expected 3600 matched birth coordinates")
            cohorts[split]=common(pairs)
        fit_times=[phi_time(r,"micro")[1] for r,_ in cohorts["fit"]]
        # 0..infinity lets us detect early and late errors, with
        # interior 20/40/60/80% thresholds defined by FIT-only.
        edges=np.asarray([0.,*np.quantile(fit_times,[.2,.4,.6,.8]),
                          np.inf])
        if not np.all(np.diff(edges)>0.):
            raise ValueError("Time edges not strictly increasing")
        result["time_edges_ns_by_geometry_defined_on_micro_fit_only"][geom]=[
            float(x) for x in edges]
        gs={}
        ref_f=joint_hist(cohorts["fit"],"micro",edges)
        ref_v=joint_hist(cohorts["validation"],"micro",edges)
        ref_noise=total_variation(ref_f,ref_v)
        gs["microscopic_fit_vs_validation_joint_TV"]=ref_noise
        for split in ("fit","validation"):
            same=cohorts[split]
            micro_joint=joint_hist(same,"micro",edges)
            mc_joint=joint_hist(same,"mc",edges)
            micro_h=angle_hemi(same,"micro")
            mc_h=angle_hemi(same,"mc")
            diff=mc_h["P_phi_ge_0"]-micro_h["P_phi_ge_0"]
            interval=hemi_cluster_bootstrap(same,replicates,
                       seed+j*10+(split=="validation"))
            x={
                "n_common_birth_survivors":len(same),
                "joint_phi_time_TV_12_by_5":total_variation(
                                              micro_joint,mc_joint),
                "joint_hist_micro_12_by_5":micro_joint.tolist(),
                "joint_hist_mc_12_by_5":mc_joint.tolist(),
                "micro_hemisphere_time":micro_h,
                "mc_hemisphere_time":mc_h,
                "MC_minus_micro_P_upper":diff,
                "P_upper_delta_cell_cluster_bootstrap_95ci":interval,
                "conditional_spatial_regions":grouped_stats(same)
              }
            gs[split]=x
            print(f"{geom} {split}: common N={len(same)}; "
                  f"joint phi-time TV={x['joint_phi_time_TV_12_by_5']:.4f}; "
                  f"micro-fit/val joint noise={ref_noise:.4f}; "
                  f"delta P_upper={diff:+.4f} (CI {interval[0]:+.4f}.."
                  f"{interval[1]:+.4f}); "
                  f"Q99 micro={micro_h['p99_time_ns']:.2f} ns,"
                  f" MC={mc_h['p99_time_ns']:.2f} ns",flush=True)
        result["geometries"][geom]=gs
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print("PASS joint transport-only distribution audit:",out)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--reference",type=Path,required=True)
    p.add_argument("--mc",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    p.add_argument("--bootstrap",type=int,default=500)
    p.add_argument("--seed",type=int,default=261027)
    a=p.parse_args()
    if not 100<=a.bootstrap<=2000:
        p.error("bootstrap must be in [100,2000]")
    report(a.reference,a.mc,a.out,a.bootstrap,a.seed)

if __name__=="__main__":
    main()
