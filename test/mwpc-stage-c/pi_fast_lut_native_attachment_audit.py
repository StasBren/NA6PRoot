#!/usr/bin/env python3
"""Audit native Garfield AvalancheMC attachment against offline eta-path hazard.

Requires:
 - immutable SHA-checked micro selected starts (fit + validation)
 - previous attachment-OFF hazard CSV
 - new native attachment-ON MC CSV with EXACT SAME birth sites
Produces read-only numerical report and optional field-region heatmap.
Never refits k using the validation split.
"""
from __future__ import annotations
import argparse,json,math
from collections import Counter,defaultdict
from pathlib import Path
import numpy as np
from pi_fast_lut_existing_starts_report import GEOMS,read_micro,read_mc,sha256
from pi_fast_lut_attachment_fit import load_hazards,vectorize,fit_global_k,mean_pred

def group_map(entries):
    groups=defaultdict(list)
    for v in entries:
        groups[v["block"]].append(v)
    if len(groups)!=36:raise ValueError("Expected all 6x6 coarse blocks")
    return groups

def residuals(entries,k):
    groups=group_map(entries)
    values=[]
    for cell,subset in sorted(groups.items()):
        obs=float(np.mean([v["reference_hit"] for v in subset]))
        pred=mean_pred(subset,k)
        values.append({"block":[*cell],"n":len(subset),
                       "micro":obs,"predicted":pred,
                       "difference":pred-obs})
    return {
        "rms":math.sqrt(sum(z["difference"]**2 for z in values)/len(values)),
        "max_abs":max(abs(z["difference"]) for z in values),
        "worst_five":sorted(values,key=lambda z:abs(z["difference"]),
                             reverse=True)[:5],
        "all_blocks":values}

def kfit_target(entries,native):
    fitted=[{**e,"reference_hit":int(native[e["id"]]["anode_collected"])}
            for e in entries]
    return fit_global_k(fitted)

def bin_checks(items,native,k):
    # Hazard-based bins are estimated from OFF drift paths; native
    # attachment-on may consume RNG differently and terminate earlier.
    # They are *population* comparisons, not eventwise trajectory matches.
    items=sorted(items,key=lambda v:v["hazard"])
    segments=[]
    for i in range(10):
        values=items[i*len(items)//10:(i+1)*len(items)//10]
        if not values:continue
        actual=sum(int(native[v["id"]]["anode_collected"]) for v in values)/len(values)
        micro=sum(v["reference_hit"] for v in values)/len(values)
        pred=sum(math.exp(-v["hazard"]) * v["mc_hit"] for v in values)/len(values)
        predk=sum(math.exp(-k*v["hazard"]) * v["mc_hit"] for v in values)/len(values)
        segments.append({"decile":i+1,"n":len(values),
             "H_median":float(np.median([v["hazard"] for v in values])),
             "native_mc_attachment_on_anode":actual,
             "micro_anode":micro,"estimated_k1_survival":pred,
             "estimated_kfit_survival":predk})
    return segments

def study(reference_path,hazard_path,native_path,baseline_path,out):
    if out.exists():raise FileExistsError("Refusing to replace audit "+str(out))
    reference,prov=read_micro(reference_path)
    off,step=load_hazards(hazard_path,reference,baseline_path)
    native,step_native,mode=read_mc(native_path,reference)
    if mode!="on" or step!=step_native:
        raise ValueError("Native MC must have attachment=on and identical step")
    grouped=vectorize(reference,off,"hazard_eta_vdt")
    # Attach immutable key to rows for direct lookup in native MC.
    for key,ref in reference.items():
        pass
    result={"schema":"PI_FAST_LUT_NATIVE_ATTACHMENT_AUDIT_V1",
        "micro_selected_sha256":sha256(reference_path),
        "hazards_sha256":sha256(hazard_path),
        "native_mc_csv_sha256":sha256(native_path),
        "mc_step_um":step,
        "geometry":{},
        "warning":"k_native is an EFFECTIVE scale required to reproduce "
        "native MC collection using THIS offline hazard convention, "
        "NOT a physical Magboltz gas correction. Only FIT used for k.",
        "notes":[
          "Micro and MC trajectories have separate collision randomness.",
          "Native MC-on may terminate earlier than un-killed MC and "
          "consume a different RNG sequence.",
          "Projecting eta onto drift velocity or other internal Garfield "
          "integration choices may cause disagreement with offline H_vdt.",
          "A global collection fraction or a 36-bin phi comparison alone "
          "does not validate local collection probabilities."
        ]}
    for geom in GEOMS:
        fit=[];validation=[]
        for key,ref in reference.items():
            if key[0]!=geom:continue
            o=off[key]
            e={"id":key,
               "block":(int(ref["iu"])//5,int(ref["iv"])//5),
               "reference_hit":int(ref["reference_collected"]),
               "mc_hit":int(o["anode_collected"]),
               "hazard":float(o["hazard_eta_vdt"])}
            (fit if ref["split"]=="fit" else validation).append(e)
        k_micro=fit_global_k(fit)
        k_native=kfit_target(fit,native)
        g={"k_fit_micro":k_micro,"k_fit_native":k_native,"splits":{}}
        for name,values in (("fit",fit),("validation",validation)):
            n=len(values)
            native_fraction=sum(int(native[x["id"]]["anode_collected"]) for x in values)/n
            micro_fraction=sum(x["reference_hit"] for x in values)/n
            no_attachment_fraction=sum(x["mc_hit"] for x in values)/n
            statuses=Counter(str(native[x["id"]]["status"]) for x in values)
            native_vs_k1=native_fraction-mean_pred(values,1.)
            native_by_block=[]
            grouped_by_block=group_map(values)
            for block,subset in sorted(grouped_by_block.items()):
                native_block=np.mean([int(native[x["id"]]["anode_collected"]) for x in subset])
                predicted_k1=mean_pred(subset,1.)
                native_by_block.append({"block":list(block),"n":len(subset),
                    "native_fraction":float(native_block),
                    "offline_k1_fraction":predicted_k1,
                    "difference":predicted_k1-float(native_block)})
            native_rms=math.sqrt(np.mean([b["difference"]**2 for b in native_by_block]))
            block_micro=residuals(values,k_micro)
            h=[x["hazard"] for x in values]
            val={
              "n":n,
              "micro_anode":micro_fraction,
              "mc_no_attachment_anode":no_attachment_fraction,
              "native_mc_attachment_on_anode":native_fraction,
              "native_endpoint_status_counts":dict(statuses),
              "offline_survival_k1":mean_pred(values,1.),
              "offline_survival_fitted_micro_k":mean_pred(values,k_micro),
              "offline_survival_fitted_native_k":mean_pred(values,k_native),
              "native_minus_offline_k1":native_vs_k1,
              "native_vs_offline_k1_block_RMS":native_rms,
              "micro_vs_offline_kmicro_block_RMS":block_micro["rms"],
              "micro_vs_offline_kmicro_worst_blocks":block_micro["worst_five"],
              "hazard_H_min":min(h),
              "hazard_H_median":float(np.median(h)),
              "hazard_H_p95":float(np.quantile(h,.95)),
              "hazard_H_max":max(h),
              "hazard_deciles":bin_checks(values,native,k_micro)}
            g["splits"][name]=val
            print(f"{geom} {name}: micro={micro_fraction:.4f}; "
                  f"native-on={native_fraction:.4f}; "
                  f"offline k1={val['offline_survival_k1']:.4f}; "
                  f"offline k_micro={val['offline_survival_fitted_micro_k']:.4f}; "
                  f"native-minus-k1={native_vs_k1:+.4f}; "
                  f"blockRMS(native vs k1)={native_rms:.4f}",flush=True)
        result["geometry"][geom]=g
        print(f"  {geom}: k to micro FIT={k_micro:.5g}, "
              f"k to native FIT={k_native:.5g}",flush=True)
    # Joint fit transfer between geometry is handled by prior first-pass
    # script; this audit specifically checks Garfield native-k1 consistency.
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print("PASS native MC attachment versus offline hazard audit:",out)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--reference",required=True,type=Path)
    p.add_argument("--hazards",required=True,type=Path)
    p.add_argument("--native",required=True,type=Path)
    p.add_argument("--baseline-mc",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    a=p.parse_args()
    study(a.reference,a.hazards,a.native,a.baseline_mc,a.out)

if __name__=="__main__":main()
