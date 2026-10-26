#!/usr/bin/env python3
"""Calibrate k*eta(E) on FIT only, audit validation and geometry transfer.

No Garfield drifts or microscopic calculations. This reads the immutable
selected reference and the matching un-killed MC *hazard* CSV.
The calibrated exp(-k*H) is a probabilistic response surrogate, not a
drop-in replacement for Garfield's attachment/multiplication algorithm.

H_vdt = sum eta(E_mid)*|v_d(E_mid)|*delta_t is the PRIMARY convention.
H_chord = sum eta(E_mid)*|delta_x_MC| is a diffusion-length sensitivity
control. Their difference diagnoses whether the meaning of ds matters.
"""
from __future__ import annotations

import argparse,csv,hashlib,json,math
from collections import defaultdict
from pathlib import Path

import numpy as np

from pi_fast_lut_existing_starts_report import (
    GEOMS,read_micro,read_mc,sha256
)

NEEDED={"hazard_eta_vdt","hazard_eta_chord","drift_length_cm",
        "chord_length_cm","path_steps","eta_E_min_Vcm","eta_E_max_Vcm"}


def load_hazards(path,reference,baseline=None):
    rows,step,mode=read_mc(path,reference)
    if mode!="off":
        raise ValueError("Un-killed hazard extraction requires MC attachment off")
    if baseline is not None:
        old,oldstep,oldmode=read_mc(baseline,reference)
        if oldstep!=step or oldmode!="off":
            raise ValueError("Baseline MC setting mismatch")
        # Detect any RNG or stepping changes when extraction was added.
        # Exact starts and seeds are already checked by read_mc.
        for key,row in rows.items():
            other=old[key]
            for prop in ("status","anode_collected","field_wire_collected",
                         "wire_index"):
                if int(row[prop])!=int(other[prop]):
                    raise ValueError(f"Changed un-killed MC endpoint {key}: {prop}")
            for prop in ("u_end_mm","v_end_mm","drift_time_ns","phi_hit_deg"):
                a=float(row[prop]);b=float(other[prop])
                if math.isnan(a) and math.isnan(b):
                    continue
                if not math.isclose(a,b,rel_tol=0,abs_tol=2e-8):
                    raise ValueError(f"Changed MC trajectory by hazard extraction {key}: {prop}")
    count=0
    for key,row in rows.items():
        if not NEEDED.issubset(row):
            raise ValueError("Hazard fields absent")
        for name in ("hazard_eta_vdt","hazard_eta_chord","drift_length_cm",
                     "chord_length_cm"):
            x=float(row[name])
            if not math.isfinite(x) or x<0:
                raise ValueError(f"Invalid {name} at {key}: {x}")
        if int(row["path_steps"])<1:
            raise ValueError(f"No recorded drift steps {key}")
        if float(row["chord_length_cm"])+1e-10<float(row["drift_length_cm"])*.1:
            raise ValueError(f"Implausibly short stochastic path at {key}")
        count+=1
    if count!=len(reference):
        raise ValueError("Incomplete hazard table")
    return rows,step


def vectorize(reference,rows,hazard):
    grouped=defaultdict(list)
    for key,ref in reference.items():
        mc=rows[key]
        e={
          "geom":ref["geometry"],"split":ref["split"],
          "cell":(int(ref["iu"]),int(ref["iv"])),
          "block":(int(ref["iu"])//5,int(ref["iv"])//5),
          "reference_hit":int(ref["reference_collected"]),
          "mc_hit":int(mc["anode_collected"]),
          "hazard":float(mc[hazard]),
          "ref_phi":float(ref["reference_phi_hit_deg"]),
          "mc_phi":float(mc["phi_hit_deg"]),
          "ref_time":float(ref["reference_drift_time_ns"]),
          "mc_time":float(mc["drift_time_ns"]),
          "start_v":float(mc["v_start_mm"]),
          "start_u":float(mc["u_start_mm"]),
        }
        if (e["mc_hit"] or e["reference_hit"]) and e["hazard"]<0:
            raise ValueError("Negative hazard")
        if not e["mc_hit"] and e["hazard"]!=0:
            # Since hazard is integrated for non-anode endpoints, this is fine.
            pass
        grouped[(e["geom"],e["split"])].append(e)
    return grouped


def predicted_weight(e,k):
    return float(e["mc_hit"])*math.exp(-k*e["hazard"])


def mean_pred(items,k):
    return sum(predicted_weight(e,k) for e in items)/len(items)


def fit_global_k(items):
    if not items:raise ValueError("No calibration records")
    target=sum(e["reference_hit"] for e in items)/len(items)
    maxsurvive=sum(e["mc_hit"] for e in items)/len(items)
    if maxsurvive<target-1e-9:
        raise ValueError("Even no-attachment MC has fewer collected electrons than reference")
    if abs(maxsurvive-target)<1.e-12:
        return 0.
    floor=sum(e["mc_hit"] and e["hazard"]==0 for e in items)/len(items)
    if target<=floor+1e-10:
        raise ValueError(f"Cannot fit k: micro collected fraction {target} "
                         f"is at/below hazard-zero floor {floor}")
    lo,hi=0.,1.
    while mean_pred(items,hi)>target:
        hi*=2.
        if hi>1e7:raise ValueError("Unable to bracket global k")
    for _ in range(65):
        mid=(lo+hi)*.5
        if mean_pred(items,mid)>target:
            lo=mid
        else:
            hi=mid
    return .5*(lo+hi)


def one_geom(items,k):
    if not items:raise ValueError("No validation group")
    obs=np.array([e["reference_hit"] for e in items],dtype=np.float64)
    prob=np.array([predicted_weight(e,k) for e in items])
    n=len(items)
    nref=int(obs.sum())
    # This uncertainty is approximate: model is a second stochastic
    # process and ignores correlation induced by shared birth cells.
    naive_se=float(math.sqrt((obs.mean()*(1.-obs.mean())
       +np.mean(prob*(1.-prob)))/n))
    delta=float(prob.mean()-obs.mean())
    counts_ob=np.zeros(36,dtype=float)
    counts_mc=np.zeros(36,dtype=float)
    ref_time=[]
    mc_time=[]
    weighted_time=[]
    for e,w in zip(items,prob):
        if e["reference_hit"]:
            phi=e["ref_phi"]
            if not math.isfinite(phi) or not math.isfinite(e["ref_time"]):
                raise ValueError("Reference survivor has invalid phi/time")
            counts_ob[max(0,min(35,int((phi+180.)/10.)))]+=1.
            ref_time.append(e["ref_time"])
        if e["mc_hit"]:
            phi=e["mc_phi"]
            if not math.isfinite(phi) or not math.isfinite(e["mc_time"]):
                raise ValueError("MC survivor has invalid phi/time")
            counts_mc[max(0,min(35,int((phi+180.)/10.)))]+=w
            mc_time.append(e["mc_time"])
            weighted_time.append(w)
    pred_sum=sum(weighted_time)
    pobs=counts_ob/counts_ob.sum() if counts_ob.sum()>0 else counts_ob
    pmc=counts_mc/counts_mc.sum() if counts_mc.sum()>0 else counts_mc
    # Diagnostics by grouped 5x5 spatial supercells (100 starts each in fit
    # or validation). Do not assess tight spatial gates on 4e / tiny cells.
    blocks=defaultdict(list)
    for e in items:
        blocks[e["block"]].append(e)
    regional=[]
    for block,bb in sorted(blocks.items()):
        actual=sum(e["reference_hit"] for e in bb)/len(bb)
        fitted=mean_pred(bb,k)
        regional.append({"block_uv":list(block),"n":len(bb),
            "microscopic_fraction":actual,"predicted_fraction":fitted,
            "residual":fitted-actual})
    rms_reg=math.sqrt(sum(z["residual"]**2 for z in regional)/len(regional))
    return {
      "n":n,"micro_anode_fraction":float(obs.mean()),
      "mc_no_attachment_anode_fraction":float(np.mean([e["mc_hit"] for e in items])),
      "predicted_survival_fraction":float(prob.mean()),
      "delta_predicted_minus_micro":delta,
      "approx_95ci_delta_iid":[delta-1.96*naive_se,delta+1.96*naive_se],
      "phi_tv_36bins_weighted_survivors":float(.5*np.abs(pobs-pmc).sum()),
      "phi_pdf_micro_36bins":pobs.tolist(),
      "phi_pdf_model_36bins":pmc.tolist(),
      "mean_drift_micro_ns":float(np.mean(ref_time)) if ref_time else None,
      "mean_drift_model_weighted_ns":
          (float(np.average(mc_time,weights=weighted_time))
           if pred_sum>0 else None),
      "coarse_block_rms_collection_residual":rms_reg,
      "worst_6_abs_coarse_block_residuals":sorted(regional,
          key=lambda z:abs(z["residual"]),reverse=True)[:6],
      "region_n":len(regional),
      "warnings":"Global fit is constrained by efficiency; regional "
          "and phi/time residuals are independent diagnostic checks."
    }


def evaluate(reference_path,hazards_path,out,baseline):
    if out.exists():raise FileExistsError(f"Refusing calibrated output overwrite: {out}")
    reference,manifest=read_micro(reference_path)
    rows,step=load_hazards(hazards_path,reference,baseline)
    result={
      "schema":"PI_FAST_LUT_ATTACHMENT_K_FIRST_PASS_V1",
      "source_micro_selection_sha256":sha256(reference_path),
      "hazard_csv_sha256":sha256(hazards_path),
      "source_frozen_micro_shas":{
        k:v["sha256"] for k,v in manifest["source_provenance"].items()},
      "mc_step_um":step,"scope":"p3, B=0, Ar/CO2 70:30, 293.15K, 760Torr",
      "fitting":"Deterministic bisection to match FIT microscopic global "
        "anode fraction using mean(MC_anode * exp(-k*hazard)).",
      "models":{},
      "cautions":[
         "Exploratory first calibration, not approved production physics.",
         "Global fit does NOT guarantee local collection or angular/time agreement.",
         "Reference and MC stochastic trajectories are independent even "
         "for matched birth sites.",
         "v_d*dt and chord ds integration differ due to diffusion; "
         "neither is yet demonstrated identical to Garfield internal attachment.",
         "Error intervals shown here are naive iid approximations, not "
         "cell-clustered confidence intervals; final approval requires stronger uncertainty tests.",
         "A predicted survival weight is not an event-level Bernoulli draw. "
         "Do not reuse weighted mean as a substitute for eventwise fluctuations."
       ]
    }
    for hazard in ("hazard_eta_vdt","hazard_eta_chord"):
        grouped=vectorize(reference,rows,hazard)
        fit_a=grouped[(GEOMS[0],"fit")]
        fit_b=grouped[(GEOMS[1],"fit")]
        ks={"joint_both_geometries":fit_global_k(fit_a+fit_b),
            "fit_noFW_only":fit_global_k(fit_a),
            "fit_FW75_only":fit_global_k(fit_b)}
        hazard_result={"estimated_k":ks,"assessments":{}}
        print("\nHazard convention:",hazard)
        for label,k in ks.items():
            print(f"  fitted {label} k={k:.6g}")
            hazard_result["assessments"][label]={}
            for geom in GEOMS:
                for split in ("fit","validation"):
                    result_here=one_geom(grouped[(geom,split)],k)
                    hazard_result["assessments"][label][geom+"/"+split]=result_here
                    if split=="validation":
                        print(f"    {geom} validation: micro="
                              f"{result_here['micro_anode_fraction']:.4f}; "
                              f"pred={result_here['predicted_survival_fraction']:.4f}; "
                              f"phiTV={result_here['phi_tv_36bins_weighted_survivors']:.4f}; "
                              f"region RMS={result_here['coarse_block_rms_collection_residual']:.4f}")
        result["models"][hazard]=hazard_result
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print("PASS exploratory attachment calibration:",out)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--reference",type=Path,required=True)
    p.add_argument("--hazards",type=Path,required=True)
    p.add_argument("--baseline-mc",type=Path,default=None)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    evaluate(a.reference,a.hazards,a.out,a.baseline_mc)


if __name__=="__main__":
    main()
