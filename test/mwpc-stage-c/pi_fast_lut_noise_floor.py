#!/usr/bin/env python3
"""Read-only intrinsic microscopic angular-noise and MC step-convergence check.

The microscopic reference has fit/validation trials with independent
original track IDs in EACH of 900 cells. Cell-stratified permutation
estimates the noise floor of a 36-bin phi total-variation (TV)
comparison due to finite reference statistics alone.

Important: MC attachment-off has different survivors than micro.
This reference-only null is a NOISE BENCHMARK, NOT a formal p-value for
model equivalence. Matching only this TV threshold cannot certify a
fast transport model.
"""
from __future__ import annotations

import argparse
import json
import math
import random
from collections import defaultdict
from pathlib import Path

import numpy as np

from pi_fast_lut_existing_starts_report import (
    GEOMS, read_micro, read_mc, sha256, phi_pdf
)


def tv_from_counts(first,second):
    a=np.asarray(first,dtype=np.float64)
    b=np.asarray(second,dtype=np.float64)
    if not a.sum() or not b.sum():
        return None
    return float(.5*np.abs(a/a.sum()-b/b.sum()).sum())


def angle_bin(row,reference):
    is_hit=int(row["reference_collected"] if reference else row["anode_collected"])
    if not is_hit:
        return -1
    phi=float(row["reference_phi_hit_deg"] if reference else row["phi_hit_deg"])
    if not math.isfinite(phi) or phi < -180 or phi > 180:
        raise ValueError("Invalid collected electron phi angle")
    return max(0,min(35,int(math.floor((phi+180.)/10.))))


def hist_for(rows,reference=True):
    counts=np.zeros(36,dtype=np.int64)
    for r in rows:
        j=angle_bin(r,reference)
        if j>=0:counts[j]+=1
    return counts


def reference_null(rows,permutations,seed):
    cells=defaultdict(list)
    fit,validation=[],[]
    for r in rows:
        key=(int(r["iu"]),int(r["iv"]))
        cells[key].append(r)
        if r["split"]=="fit":
            fit.append(r)
        elif r["split"]=="validation":
            validation.append(r)
        else:
            raise ValueError("Invalid fit/validation split label")
    if len(cells)!=900:
        raise ValueError(f"Expected 900 populated microscopic cells, got {len(cells)}")
    if not fit or not validation or len(fit)!=len(validation):
        raise ValueError("Microscopic reference splits must be balanced")
    sizes={len(v) for v in cells.values()}
    if len(sizes)!=1 or next(iter(sizes))%2:
        raise ValueError("Need equal, even number of trial starts per LUT cell")
    nfit=next(iter(sizes))//2
    for values in cells.values():
        if sum(r["split"]=="fit" for r in values)!=nfit:
            raise ValueError("Unbalanced microscopic splits within spatial cell")
    f=hist_for(fit)
    v=hist_for(validation)
    observed=tv_from_counts(f,v)
    rng=random.Random(seed)
    histogram_cells=[
        [angle_bin(r,True) for r in rows] for _,rows in sorted(cells.items())
    ]
    null=[]
    for _ in range(permutations):
        h1=np.zeros(36,dtype=np.int64)
        h2=np.zeros(36,dtype=np.int64)
        for angles in histogram_cells:
            shuffled=rng.sample(angles,len(angles))
            for bin_index in shuffled[:nfit]:
                if bin_index>=0:h1[bin_index]+=1
            for bin_index in shuffled[nfit:]:
                if bin_index>=0:h2[bin_index]+=1
        null.append(tv_from_counts(h1,h2))
    return {
      "microscopic_fit_n_start":len(fit),
      "microscopic_validation_n_start":len(validation),
      "microscopic_fit_survivors":int(f.sum()),
      "microscopic_validation_survivors":int(v.sum()),
      "observed_fit_vs_validation_phi_TV_36bins":observed,
      "cell_stratified_permutation_count":permutations,
      "permutation_TV_p05_p50_p95":
          [float(np.quantile(null,q)) for q in (.05,.5,.95)],
      "notes":"Within-cell fixed-exposure permutation null conditional on "
              "the existing 8 microscopic trials per cell. "
              "Not a formal acceptance threshold for MC with different survival."
    }


def matched_step_comparison(rows_a,rows_b):
    """Compare MC stochastic outputs at same starts across two steps."""
    if set(rows_a)!=set(rows_b):
        raise ValueError("Two MC step outputs do not share all trial IDs")
    a=list(rows_a.values())
    b=list(rows_b.values())
    a_hist=hist_for(a,reference=False)
    b_hist=hist_for(b,reference=False)
    a_survive=int(a_hist.sum())
    b_survive=int(b_hist.sum())
    mismatches=sum(
      int(rows_a[key]["anode_collected"])!=int(rows_b[key]["anode_collected"])
      for key in rows_a)
    drift_a=[
      float(r["drift_time_ns"]) for r in a if int(r["anode_collected"])]
    drift_b=[
      float(r["drift_time_ns"]) for r in b if int(r["anode_collected"])]
    if not all(math.isfinite(t) for t in drift_a+drift_b):
        raise ValueError("Nonfinite collected MC drift time")
    return {
      "n_same_birth_sites":len(a),
      "anode_a_fraction":a_survive/len(a),
      "anode_b_fraction":b_survive/len(b),
      "individual_anode_outcome_mismatches":mismatches,
      "phi_total_variation_36bins":tv_from_counts(a_hist,b_hist),
      "mean_drift_a_ns":sum(drift_a)/len(drift_a) if drift_a else None,
      "mean_drift_b_ns":sum(drift_b)/len(drift_b) if drift_b else None,
      "drift_us_per_start_a":sum(float(r["drift_wall_us"]) for r in a)/len(a),
      "drift_us_per_start_b":sum(float(r["drift_wall_us"]) for r in b)/len(b)
    }


def analyze(reference_path,mc_paths,out,permutations,seed):
    if out.exists():
        raise FileExistsError("Refusing to overwrite existing noise floor report "+str(out))
    ref,info=read_micro(reference_path)
    report={
        "study":"Fast LUT intrinsic micro statistical noise and MC step convergence",
        "reference_sha256":sha256(reference_path),
        "frozen_micro_input_provenance":{
          k:v["sha256"] for k,v in info["source_provenance"].items()},
        "reference_noise":{},
        "mc_steps":{},
        "mc_step_pairwise":{},
        "warnings":[
          "The micro split TV and its stratified null quantify sampling variation; "
          "they are not an equivalence test for MC.",
          "MC no-attachment and micro include different anode survivor populations.",
          "Same MC RNG seeds across different step sizes do not imply "
          "identical stochastic drift/diffusion trajectories.",
          "Phi TV and collection alone do not validate drift-time tails or "
          "eventwise NA60+ detector response."
        ]
    }
    mc_by_step={}
    for geom in GEOMS:
        rows=[v for key,v in ref.items() if key[0]==geom]
        result=reference_null(rows,permutations,seed+(1 if geom==GEOMS[0] else 2))
        report["reference_noise"][geom]=result
        print(geom,": micro fit/validation intrinsic phi TV=",
              f"{result['observed_fit_vs_validation_phi_TV_36bins']:.4f}",
              ", permuted 5/50/95% TV=",result["permutation_TV_p05_p50_p95"],
              flush=True)
    for path in mc_paths:
        data,step,mode=read_mc(path,ref)
        if step in mc_by_step:
            raise ValueError("Duplicate MC step")
        if mode!="off":
            raise ValueError("Noise diagnostic intentionally requires MC attachment OFF")
        mc_by_step[step]=data
        report["mc_steps"][str(step)]={
            "source":str(path),"sha256":sha256(path),"attachment":mode,
            "mean_transport_us_per_electron":sum(
               float(r["drift_wall_us"]) for r in data.values())/len(data)}
        for geom in GEOMS:
            for split in ("fit","validation"):
                group_ref=[v for k,v in ref.items() if k[0]==geom and v["split"]==split]
                group_mc=[data[k] for k,v in ref.items()
                          if k[0]==geom and v["split"]==split]
                ref_angle=hist_for(group_ref)
                mc_angle=hist_for(group_mc,False)
                tv=tv_from_counts(ref_angle,mc_angle)
                report["mc_steps"][str(step)][geom+"/"+split]={
                    "mc_vs_micro_phi_TV_36bins":tv,
                    "micro_anode_survivors":int(ref_angle.sum()),
                    "mc_anode_survivors":int(mc_angle.sum())}
                print(f"{step:g}um {geom} {split}: "
                      f"MC vs micro phi TV={tv:.4f}, micro anode="
                      f"{int(ref_angle.sum())}/{len(group_ref)}, MC anode="
                      f"{int(mc_angle.sum())}/{len(group_mc)}",flush=True)
    steps=sorted(mc_by_step)
    for i,step_a in enumerate(steps):
        for step_b in steps[i+1:]:
            for geom in GEOMS:
                for split in ("fit","validation"):
                    data_a={k:v for k,v in mc_by_step[step_a].items()
                        if k[0]==geom and v["split"]==split}
                    data_b={k:v for k,v in mc_by_step[step_b].items()
                        if k[0]==geom and v["split"]==split}
                    key=f"{step_a:g}um_vs_{step_b:g}um/{geom}/{split}"
                    report["mc_step_pairwise"][key]=matched_step_comparison(
                        data_a,data_b)
                    z=report["mc_step_pairwise"][key]
                    print(f"MC step {key}: phi TV={z['phi_total_variation_36bins']:.4f}"
                          f", mean drift {z['mean_drift_a_ns']:.2f}ns"
                          f" vs {z['mean_drift_b_ns']:.2f}ns",flush=True)
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print("PASS noise-floor and step diagnostic:",out)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--reference",type=Path,required=True)
    p.add_argument("--mc",type=Path,nargs="+",required=True)
    p.add_argument("--out",type=Path,required=True)
    p.add_argument("--permutations",type=int,default=300)
    p.add_argument("--seed",type=int,default=271026)
    a=p.parse_args()
    if not 50<=a.permutations<=5000:
        p.error("permutations must be 50..5000")
    analyze(a.reference,a.mc,a.out,a.permutations,a.seed)


if __name__=="__main__":
    main()
