#!/usr/bin/env python3
"""Transport-only (attachment OFF) audit against immutable microscopic LUT500.

Run entirely on existing selected reference + MC CSVs. In the microscopic
reference electron attachment was ON; hence this comparison is a
SURVIVOR-CONDITIONED transport diagnostic, not validation of an
attachment-free microscopic benchmark and not an electronics-response
prediction.

Two cohorts are reported separately:
  'each_own_anode_survivors': endpoints conditional on each model's own hit;
  'same_birth_both_anode_survivors': same birth coordinates, microscopic
       and MC both hit anode; their collision histories remain independent.

Outputs: JSON with provenance, tables of quantile/CDF differences,
36-bin angles, wire indices, 6x6 spatial mean-time residuals, cell-bootstrap
CIs for mean-time bias on common survivors; two optional figures.
No Garfield, table generation, or file overwrite.
"""
from __future__ import annotations

import argparse
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

from pi_fast_lut_existing_starts_report import GEOMS, read_micro, read_mc, sha256

ANGLES=np.arange(-180,181,10,dtype=float)
QUANTILES=(.1,.5,.9,.99)
MODEL_STATUS="transport only; attachment intentionally off"


def arr_quantiles(xs):
    return {str(int(q*100)):float(np.quantile(xs,q)) for q in QUANTILES}


def angles(rows,key):
    output=[]
    for row in rows:
        h=int(row["reference_collected" if key=="micro" else "anode_collected"])
        if not h:continue
        phi=float(row["reference_phi_hit_deg" if key=="micro" else "phi_hit_deg"])
        if not math.isfinite(phi) or not -180<=phi<=180:
            raise ValueError("Invalid phi of anode endpoint")
        output.append(phi)
    return np.asarray(output,dtype=float)


def times(rows,key):
    output=[]
    for row in rows:
        h=int(row["reference_collected" if key=="micro" else "anode_collected"])
        if not h:continue
        t=float(row["reference_drift_time_ns" if key=="micro" else "drift_time_ns"])
        if not math.isfinite(t) or t<0:
            raise ValueError("Invalid collected-electron drift time")
        output.append(t)
    return np.asarray(output,dtype=float)


def hist36(a):
    counts=np.histogram(a,bins=ANGLES)[0]
    if counts.sum()==0: raise ValueError("No angle endpoints")
    return counts/counts.sum()


def tv(a,b):
    return float(.5*np.abs(hist36(a)-hist36(b)).sum())


def ecdf_ks(a,b):
    """Two-sample KS distance on empirical distributions, O(N log N)."""
    if not len(a) or not len(b):raise ValueError("Empty distribution")
    aa=np.sort(a)
    bb=np.sort(b)
    grid=np.unique(np.concatenate((aa,bb)))
    return float(np.max(np.abs(
        np.searchsorted(aa,grid,side="right")/len(aa)
      - np.searchsorted(bb,grid,side="right")/len(bb))))


def summary_distribution(ref,mc):
    if not len(ref) or not len(mc):
        raise ValueError("Missing collected electrons in distribution")
    return {
       "microscopic": {"n":int(len(ref)),"mean":float(ref.mean()),
                       "rms":float(ref.std()),"quantiles_ns":arr_quantiles(ref)},
       "fast_mc": {"n":int(len(mc)),"mean":float(mc.mean()),
                   "rms":float(mc.std()),"quantiles_ns":arr_quantiles(mc)},
       "mean_difference_mc_minus_micro_ns":float(mc.mean()-ref.mean()),
       "rms_difference_mc_minus_micro_ns":float(mc.std()-ref.std()),
       "quantile_difference_mc_minus_micro_ns":{
           label: float(np.quantile(mc,int(label)/100)
                       - np.quantile(ref,int(label)/100))
           for label in ("10","50","90","99")},
       "KS_distribution_distance":ecdf_ks(ref,mc)
    }


def wire_probabilities(ref_rows,mc_rows):
    def prob(rows,flag,wire):
        counts=Counter()
        collected=0
        for row in rows:
            if not int(row[flag]):continue
            j=int(row[wire])
            if not -4<=j<=4: raise ValueError("Invalid anode wire number")
            counts[j]+=1
            collected+=1
        if not collected:raise ValueError("No collected wire index")
        return {str(j):counts[j]/collected for j in range(-4,5)}
    r=prob(ref_rows,"reference_collected","reference_wire_index")
    m=prob(mc_rows,"anode_collected","wire_index")
    same=sum(int(a["reference_wire_index"])==int(b["wire_index"])
             for a,b in zip(ref_rows,mc_rows)
             if int(a["reference_collected"]) and int(b["anode_collected"]))
    both=sum(int(a["reference_collected"]) and int(b["anode_collected"])
             for a,b in zip(ref_rows,mc_rows))
    return {
       "micro_wire_pdf":r,"mc_wire_pdf":m,
       "wire_pdf_TV":float(.5*sum(abs(r[k]-m[k]) for k in r)),
       "same_wire_on_same_birth_and_both_survived_fraction":same/both if both else None,
       "same_wire_sample_n":both,
       "caveat":"Same-wire agreement is a stochastic diagnostic, not an "
                "electron-by-electron identity requirement."
    }


def bootstrap_cell_bias(pairs,seed,bootstrap=500):
    """Distribution of mean time differences on *common-survivor starts*.

    Cluster resampling 900 source cells, preserving all per-cell sampled
    microscopic + MC pairs. The two random processes are independent;
    paired *birth sites* are not paired collision paths. This CI only
    quantifies sampling error of the conditioned mean-time contrast.
    """
    bycell=defaultdict(list)
    for a,b in pairs:
        if not int(a["reference_collected"]) or not int(b["anode_collected"]):
            continue
        cell=(int(a["iu"]),int(a["iv"]))
        diff=float(b["drift_time_ns"])-float(a["reference_drift_time_ns"])
        if not math.isfinite(diff): raise ValueError("Nonfinite paired time")
        bycell[cell].append(diff)
    if not bycell: raise ValueError("No jointly collected start sites")
    agg=np.asarray([(sum(v),len(v)) for _,v in sorted(bycell.items())],
                   dtype=float)
    total=agg.sum(axis=0)
    actual=float(total[0]/total[1])
    rng=np.random.default_rng(seed)
    indices=rng.integers(0,len(agg),size=(bootstrap,len(agg)))
    bsum=agg[indices].sum(axis=1)
    results=bsum[:,0]/bsum[:,1]
    return {
      "mean_time_bias_ns":actual,
      "n_joint_survivor_starts":int(total[1]),
      "n_cells_with_survivors":len(agg),
      "cell_bootstrap_95pct_interval_ns":[
         float(x) for x in np.quantile(results,[.025,.975])],
      "bootstrap_replicates":bootstrap
    }


def blocks(pairs):
    result=[]
    for bu in range(6):
        for bv in range(6):
            loc=[(r,m) for r,m in pairs
               if int(r["iu"])//5==bu and int(r["iv"])//5==bv]
            common=[(r,m) for r,m in loc
                 if int(r["reference_collected"]) and int(m["anode_collected"])]
            if len(loc)==0:raise ValueError("Missing 6x6 region")
            if common:
                r=times([row for row,_ in common],"micro")
                m=times([row for _,row in common],"mc")
                bias=float(m.mean()-r.mean())
                micromean=float(r.mean())
                mcmean=float(m.mean())
            else: bias=micromean=mcmean=None
            result.append({
              "block_iu":bu,"block_iv":bv,
              "n_starts":len(loc),"n_joint_survivors":len(common),
              "fraction_micro_hit":sum(int(r["reference_collected"])
                           for r,_ in loc)/len(loc),
              "fraction_mc_hit":sum(int(m["anode_collected"])
                           for _,m in loc)/len(loc),
              "common_survivors_mean_time_mc_minus_micro_ns":bias,
              "common_survivors_mean_time_micro_ns":micromean,
              "common_survivors_mean_time_mc_ns":mcmean
            })
    return result


def group_report(pairs,seed,bootstrap):
    ref=[r for r,_ in pairs]
    mc=[m for _,m in pairs]
    if any(m["mc_attachment"]!="off" for m in mc):
        raise ValueError("Transport report requires MC attachment OFF")
    joint=[(a,b) for a,b in pairs
           if int(a["reference_collected"]) and int(b["anode_collected"])]
    if not joint:raise ValueError("No common survivor cohort")
    joint_micro=[a for a,_ in joint]
    joint_mc=[b for _,b in joint]
    base=summary_distribution(times(ref,"micro"),times(mc,"mc"))
    common=summary_distribution(times(joint_micro,"micro"),times(joint_mc,"mc"))
    return {
       "n_birth_sites":len(pairs),
       "n_micro_anode_survivors":sum(int(x["reference_collected"]) for x in ref),
       "n_mc_anode_survivors":sum(int(x["anode_collected"]) for x in mc),
       "n_common_survivor_birth_sites":len(joint),
       "own_survivor_time_ns":base,
       "common_birth_survivor_time_ns":common,
       "common_birth_time_cell_bootstrap":bootstrap_cell_bias(pairs,seed,bootstrap),
       "own_survivor_angular_pdf_36bins":{
          "micro":hist36(angles(ref,"micro")).tolist(),
          "mc":hist36(angles(mc,"mc")).tolist(),
          "TV":tv(angles(ref,"micro"),angles(mc,"mc"))},
       "common_birth_survivor_angular_pdf_36bins":{
          "micro":hist36(angles(joint_micro,"micro")).tolist(),
          "mc":hist36(angles(joint_mc,"mc")).tolist(),
          "TV":tv(angles(joint_micro,"micro"),angles(joint_mc,"mc"))},
       "wire_index":wire_probabilities(ref,mc),
       "coarse_spatial_blocks":blocks(pairs)
    }


def report(reference_path,mc_path,out,bootstrap,seed):
    if out.exists() or out.with_suffix(".png").exists():
        raise FileExistsError("Refusing to overwrite transport-only results")
    ref,provenance=read_micro(reference_path)
    mc,step,mode=read_mc(mc_path,ref)
    if mode!="off":raise ValueError("Only attachment-off MC is permitted")
    if step!=5.:raise ValueError("First detailed transport audit expects 5um MC")
    result={
      "schema":"PI_MWPC_TRANSPORT_ONLY_5UM_AUDIT_V1",
      "model":MODEL_STATUS,
      "scope":"p3mm, B=0, Ar/CO2 70:30, 1800V, 30um anode, 75um FW optional",
      "micro_reference_selection_sha256":sha256(reference_path),
      "fast_mc_sha256":sha256(mc_path),
      "original_microscopic_sha256s":{
           name:data["sha256"] for name,data in provenance["source_provenance"].items()},
      "methodological_limits":[
         "Microscopic data were collected with REAL electron attachment ON.",
         "Fast MC uses attachment OFF and nearly all birth electrons survive.",
         "Neither own-survivor nor common-birth comparison is a clean "
         "attachment-free microscopic transport benchmark.",
         "Microscopic and fast MC are independent stochastic collisions/"
         "diffusions even from identical birth coordinates.",
         "Same-wire agreement is not an electron-by-electron validation criterion.",
         "Spatial blocks have only 100 birth sites/split; never use their "
         "10-degree phi PDFs as high-precision truths.",
         "This report covers electron drift only; no avalanche, signal, "
         "thresholds, event asymmetry or gain validation."
      ],
      "mc_step_um":step,
      "geometries":{}
    }
    plotdata=[]
    for ig,geom in enumerate(GEOMS):
        entries={}
        for split in ("fit","validation"):
            pairs=[(row,mc[key]) for key,row in ref.items()
                   if key[0]==geom and row["split"]==split]
            if len(pairs)!=3600:raise ValueError("Expected 3600 independent starts")
            stats=group_report(pairs,seed+ig*10+(0 if split=="fit" else 1),
                               bootstrap)
            entries[split]=stats
            v=stats["common_birth_survivor_time_ns"]
            print(f"{geom} {split}: N_common="
                  f"{stats['n_common_survivor_birth_sites']}, "
                  f"dt_common={v['mean_difference_mc_minus_micro_ns']:+.3f} ns, "
                  f"KS_time={v['KS_distribution_distance']:.4f}, "
                  f"phiTV_common={stats['common_birth_survivor_angular_pdf_36bins']['TV']:.4f}, "
                  f"wirePDF_TV={stats['wire_index']['wire_pdf_TV']:.4f}",
                  flush=True)
        result["geometries"][geom]=entries
        plotdata.append((geom,entries["validation"]))
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("NOTE: matplotlib unavailable; JSON complete")
    else:
        fig,axes=plt.subplots(2,3,figsize=(16,9))
        for i,(geom,data) in enumerate(plotdata):
            a=data["common_birth_survivor_angular_pdf_36bins"]
            axes[i,0].stairs(a["micro"],ANGLES,label="Microscopic | both collected",
                             linewidth=1.7,color="#3d3d3d")
            axes[i,0].stairs(a["mc"],ANGLES,label="MC no attachment | both",
                             linewidth=1.7,color="#087eac")
            axes[i,0].set_xlabel("phi_hit [degrees]")
            axes[i,0].set_ylabel("Probability / 10 degrees")
            axes[i,0].set_title(f"{geom}: common-birth angle PDF")
            axes[i,0].legend(frameon=False,fontsize=8)
            entries=[(r,mc[key]) for key,r in ref.items()
                     if key[0]==geom and r["split"]=="validation"
                     and int(r["reference_collected"]) and
                         int(mc[key]["anode_collected"])]
            microt=times([r for r,_ in entries],"micro")
            mct=times([m for _,m in entries],"mc")
            for vals,label,color in ((microt,"Microscopic","#3d3d3d"),
                     (mct,"Fast MC","#087eac")):
                t=np.sort(vals)
                axes[i,1].step(t,np.arange(1,len(t)+1)/len(t),
                               label=label,where="post",color=color)
            axes[i,1].set_xlabel("Drift time [ns]")
            axes[i,1].set_ylabel("Cumulative probability")
            axes[i,1].set_title(f"{geom}: common-birth time CDF")
            axes[i,1].legend(frameon=False,fontsize=8)
            hm=np.full((6,6),np.nan)
            for b in data["coarse_spatial_blocks"]:
                hm[b["block_iv"],b["block_iu"]]=(
                    b["common_survivors_mean_time_mc_minus_micro_ns"])
            im=axes[i,2].imshow(hm,origin="lower",aspect="auto",cmap="coolwarm")
            fig.colorbar(im,ax=axes[i,2],label="MC minus microscopic [ns]")
            axes[i,2].set_xlabel("u block (0-5)")
            axes[i,2].set_ylabel("v block (0-5)")
            axes[i,2].set_title(f"{geom}: common-birth mean-time bias")
        fig.suptitle("Transport-only diagnostic, attachment OFF in MC; "
                     "Microscopic survivors include attachment selection")
        fig.tight_layout()
        fig.savefig(out.with_suffix(".png"),dpi=150)
        plt.close(fig)
        print("PASS plot:",out.with_suffix(".png"))
    print("PASS transport-only audit:",out)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--reference",type=Path,required=True)
    p.add_argument("--mc",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    p.add_argument("--bootstrap",type=int,default=500)
    p.add_argument("--seed",type=int,default=261026)
    args=p.parse_args()
    if not 100<=args.bootstrap<=5000:
        p.error("bootstrap must be 100..5000")
    report(args.reference,args.mc,args.out,args.bootstrap,args.seed)


if __name__=="__main__":
    main()
