#!/usr/bin/env python3
"""Read-only, SHA-audited ON/OFF review for 1000 MATCHED real J/psi muons.

No TrackHeed, Garfield or modification of reference archives.
Exact avalanche charges and EXACT per-muon ion-birth hemisphere counts.
Angular positions from <=128 unbiased Garfield ion-birth samples/avalanche.
Parent-level bootstrap uncertainties, never naive per-avalanche bootstrap.
"""
from __future__ import annotations
import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import pi_p3_noFW_jpsi1000_penning057_campaign_v1 as camp

BASE=camp.base
ROOT=camp.ROOT
OLD=camp.ORIG
OUTPUT=ROOT/"paired_1000_penning_review_v1"
NBINS=72
PHI=np.linspace(-180.,180.,NBINS+1)
CENTRES=(PHI[:-1]+PHI[1:])/2
REQUIRED_SEED=("crossing_id","primary_index","collected_e","garfield_ne",
               "ion_births","successful","capped")
REQUIRED_EVENT=("source|chunk|event|hit_index|track_id|station|det_id",
                "heed_primaries","collected_avalanche_electrons","ion_births",
                "ion_up","ion_down","ion_boundary")
REQUIRED_BIRTH=("crossing_id","primary_index","u_mm","v_mm","ion_births",
                "seed_ion_births_total","reservoir_weight")

def need(value,why):
    if not value: raise RuntimeError(why)
def to_i(x):
    z=int(x)
    need(z>=0,"Negative Garfield microscopic census")
    return z
def reader(path,req):
    with Path(path).open(newline="") as f:
        rd=csv.DictReader(f)
        need(rd.fieldnames and len(set(rd.fieldnames))==len(rd.fieldnames)
             and set(req).issubset(rd.fieldnames),
             "Missing or duplicate native CSV column "+str(path))
        for row in rd: yield row
def ident(r):
    return r["crossing_id"],int(r["primary_index"])
def hemisphere(u,d):
    n=u+d
    return (u-d)/n if n else float("nan")
def stat(x):
    a=np.asarray(x,dtype=np.float64)
    need(a.size>0,"Empty microscopic Garfield charges")
    v=a[a>0]
    need(len(v)>0,"No nonzero Garfield avalanche")
    return {
        "n":int(a.size),
        "mean_all":float(a.mean()),
        "rms_all":float(a.std(ddof=0)),
        "median_all":float(np.median(a)),
        "mean_positive":float(v.mean()),
        "rms_positive":float(v.std(ddof=0)),
        "n_zero":int(np.sum(a==0)),
        "p_zero":float(np.mean(a==0)),
        "p90":float(np.quantile(a,.90)),
        "p99":float(np.quantile(a,.99)),
        "p999":float(np.quantile(a,.999)),
        "max":int(a.max())
    }
def parent_stat(a):
    a=np.asarray(a,dtype=np.float64)
    finite=a[np.isfinite(a)]
    need(finite.size>=990,"Too many zero-ion parent muons")
    return {
        "n_defined":len(finite),
        "mean":float(finite.mean()),
        "std":float(finite.std(ddof=1)),
        "mean_abs":float(np.abs(finite).mean()),
        "p_abs_gt_05":float(np.mean(np.abs(finite)>.5)),
        "p_abs_gt_08":float(np.mean(np.abs(finite)>.8)),
        "n_undefined":int(np.sum(~np.isfinite(a)))
    }
def verify_sources():
    need(BASE.sha(OLD/"input/native_muon_crossings_input.csv")==
         BASE.sha(ROOT/"input/native_muon_crossings_input.csv"),
         "Penning ON and OFF original 1000 true Geant4 input differ")
    origins=list(BASE.csvread(OLD/"input/native_muon_crossings_input.csv",BASE.FIELDS))
    new,listmeta=BASE.check_input()
    need(len(origins)==len(new)==1000 and
         [BASE.id_of(r) for r in origins]==[BASE.id_of(r) for r in new],
         "1000 Geant4 parent identities changed")
    phy=BASE.checked_physics()
    chunks=BASE.chunks(new)
    return origins,new,phy,chunks

def load_chunks(with_angles):
    _,new,phy,chunks=verify_sources()
    out={side:{
        "muon_charge":np.zeros(1000,dtype=np.int64),
        "muon_heed":np.zeros(1000,dtype=np.int32),
        "muon_births":np.zeros(1000,dtype=np.int64),
        "A":np.full(1000,np.nan),
        "charges":[],"garfield_ne":[],"births":[],
        "seed_unsuccessful":0,"seed_capped":0,
        "angel_ions":np.zeros(NBINS),
        "angel_equal":np.zeros(NBINS),
        "ion_birth_total":0,"ion_sample_points":0,"n_ion_avalanches":0
    } for side in ("off","on")}
    for chunk_i,original_batch in enumerate(chunks):
        dnew=BASE.folder(chunk_i)
        receipt=BASE.audit_chunk(chunk_i,original_batch,phy,verify_files=True)
        need(receipt["stats"]["same_Heed_electron_positions_and_Garfield_seeds"] is True,
             "Paired Penning ON/OFF Heed and seed audit not confirmed")
        dold,oldreceipt=camp.original_sealed_chunk(chunk_i)
        if with_angles:
            name="muon_reference_ion_birth_sample.csv"
            need(BASE.sha(dold/name)==oldreceipt["file_sha256"][name],
                 "Penning OFF angular sample SHA differs from frozen source")
        # Compare the exact native parent AND per-electron order.
        paired={}
        for side,d in (("off",dold),("on",dnew)):
            parents=list(reader(d/"muon_reference_crossing.csv",REQUIRED_EVENT))
            need(len(parents)==20,"Not 20 complete real muons in native chunk")
            seedlist=list(reader(d/"muon_reference_seed.csv",REQUIRED_SEED))
            seen=set()
            for row in parents:
                key=row["source|chunk|event|hit_index|track_id|station|det_id"]
                loc=20*chunk_i+len(seen)
                need(key==BASE.id_of(original_batch[len(seen)]) and
                     key not in seen,"Native real muon ordering/source mismatch")
                seen.add(key)
                U,D,B=to_i(row["ion_up"]),to_i(row["ion_down"]),to_i(row["ion_boundary"])
                total=to_i(row["ion_births"])
                need(U+D+B==total,"Exact hemisphere full ion census incorrect")
                out[side]["muon_charge"][loc]=to_i(row["collected_avalanche_electrons"])
                out[side]["muon_heed"][loc]=to_i(row["heed_primaries"])
                out[side]["muon_births"][loc]=total
                out[side]["A"][loc]=hemisphere(U,D)
            out[side]["charges"].extend([to_i(r["collected_e"]) for r in seedlist])
            out[side]["garfield_ne"].extend([to_i(r["garfield_ne"]) for r in seedlist])
            out[side]["births"].extend([to_i(r["ion_births"]) for r in seedlist])
            out[side]["seed_unsuccessful"]+=sum(int(r["successful"])==0 for r in seedlist)
            out[side]["seed_capped"]+=sum(to_i(r["capped"]) for r in seedlist)
            need(sum(to_i(r["collected_e"]) for r in seedlist)==
                 sum(int(x["collected_avalanche_electrons"]) for x in parents),
                 "Parent collected charge differs from seed census")
            need(sum(to_i(r["ion_births"]) for r in seedlist)==
                 sum(int(x["ion_births"]) for x in parents),
                 "Parent ion birth counts differ from seed census")
            paired[side]={"seed":seedlist,"parents":parents}
            if with_angles:
                angles_for_one(d,seedlist,out[side])
        x,y=paired["off"]["seed"],paired["on"]["seed"]
        need(len(x)==len(y),"The paired sample lost real starting electrons")
        for a,b in zip(x,y):
            need(ident(a)==ident(b),
                 "Penning ON/OFF original Garfield primary identity differs")
        need(np.array_equal(out["off"]["muon_heed"][20*chunk_i:20*(chunk_i+1)],
                            out["on"]["muon_heed"][20*chunk_i:20*(chunk_i+1)]),
             "Penning changed original number of TrackHeed electrons per muon")
        if (chunk_i+1)%10==0:
            print(f"PASS SHA + exact ON/OFF source and parent ion census: "
                  f"{chunk_i+1}/50 chunks",flush=True)
    need(sum(x["muon_heed"].sum() for x in out.values())==
         len(out["on"]["charges"])+len(out["off"]["charges"]),
         "Complete exact TrackHeed electron count inconsistent")
    for s in ("off","on"):
        need(int(out[s]["muon_charge"].sum())==sum(out[s]["charges"]) and
             int(out[s]["muon_births"].sum())==sum(out[s]["births"]),
             "All 1000 parent-to-seed microscopic counts inconsistent")
        need(out[s]["seed_capped"]==0,"An avalanche was capped; cannot trust gain")
    return out

def angles_for_one(d,seedlist,info):
    expected={ident(x):to_i(x["ion_births"]) for x in seedlist}
    count={}
    birth_total=0
    for x in seedlist:
        b=to_i(x["ion_births"])
        if b>0:
            info["n_ion_avalanches"]+=1
            birth_total+=b
    info["ion_birth_total"]+=birth_total
    for row in reader(d/"muon_reference_ion_birth_sample.csv",
                      ("crossing_id","primary_index","u_mm","v_mm",
                       "seed_ion_births_total","reservoir_weight")):
        k=ident(row)
        need(k in expected and expected[k]>0,
             "Source has an ion not associated with a real Garfield avalanche")
        b=expected[k]; m=min(128,b)
        need(int(row["seed_ion_births_total"])==b,
             "Weighted ion reservoir has incorrect true ion count")
        w=float(row["reservoir_weight"])
        need(abs(w-b/m)<1e-7*max(1.,b/m),
             "Ion reservoir sampling weight is not uniform")
        phi=math.degrees(math.atan2(
            float(row["v_mm"]),
            float(row["u_mm"])-3*camp.cpp_round(float(row["u_mm"])/3)))
        idx=max(0,min(NBINS-1,int((phi+180)*NBINS/360.)))
        info["angel_ions"][idx]+=w
        info["angel_equal"][idx]+=1/m
        count[k]=count.get(k,0)+1
        info["ion_sample_points"]+=1
    need(all(count.get(k,0)==min(128,b)
             for k,b in expected.items() if b>0),
         "A nonzero ion avalanche lost reservoir positions")
def bootstrap(on,off,iterations=2000,seed=20261109):
    rng=np.random.default_rng(seed)
    N=len(on["muon_charge"])
    need(N==1000,"Parent bootstrap expects 1000 paired independent muons")
    gainratio=[]
    abs_diff=[]
    for _ in range(iterations):
        ids=rng.integers(0,N,N)
        num=on["muon_charge"][ids].sum()
        den=off["muon_charge"][ids].sum()
        gainratio.append(num/den if den>0 else float("nan"))
        aon=on["A"][ids];aoff=off["A"][ids]
        mask=np.isfinite(aon)&np.isfinite(aoff)
        abs_diff.append(np.mean(np.abs(aon[mask])-np.abs(aoff[mask])))
    gainratio=np.asarray(gainratio)
    abs_diff=np.asarray(abs_diff)
    need(np.isfinite(gainratio).all() and np.isfinite(abs_diff).all(),
         "Paired parent bootstrap failed")
    return {"iterations":iterations,"seed":seed,
            "ratio_95ci":np.quantile(gainratio,[.025,.975]).tolist(),
            "delta_mean_abs_A_95ci":np.quantile(abs_diff,[.025,.975]).tolist()}

def create():
    need(not OUTPUT.exists(),
         "Refusing to overwrite previous SHA-audited data review "+str(OUTPUT))
    results=load_chunks(with_angles=True)
    on,off=results["on"],results["off"]
    aon,aoff=on["A"],off["A"]
    defined=np.isfinite(aon)&np.isfinite(aoff)
    ratio=sum(on["charges"])/sum(off["charges"])
    summary={
        "schema":"P3_NOFW_PENNING057_1000_PAIRED_GARFIELD_REVIEW_V1",
        "independent_true_Jpsi_parents":1000,
        "n_garfield_avalanches":len(on["charges"]),
        "physics":"3mm-pitch, 30um-anodes, +/-3mm-cathodes, 1800V, Ar/CO2 70:30",
        "Penning_ON":"Ar Penning r=0.57",
        "Penning_OFF":"Original no-penning reference",
        "same_Geant4_Heed_primaries_Garfield_seeds":True,
        "comparison":{
            "collected_gain_ratio_on_off":ratio,
            "total_created_ions_on_off":sum(on["births"])/sum(off["births"]),
            "mean_absolute_hemisphere_A_difference_on_minus_off":
                float(np.mean(np.abs(aon[defined])-np.abs(aoff[defined]))),
            "pearson_parent_A_on_vs_off":float(np.corrcoef(aon[defined],aoff[defined])[0,1]),
            "parents_with_both_defined":int(defined.sum())
        },
        "parent_bootstrap":bootstrap(on,off),
    }
    for k,x in results.items():
        summary[k]={
            "gain_collected_e":stat(x["charges"]),
            "garfield_ne":stat(x["garfield_ne"]),
            "ion_births_per_primary":stat(x["births"]),
            "A_ion_hemi_per_parent":parent_stat(x["A"]),
            "total_ion_births":int(sum(x["births"])),
            "n_failed_Garfield_calls":x["seed_unsuccessful"],
            "n_capped_Garfield_calls":x["seed_capped"],
            "angular_n_points":x["ion_sample_points"],
            "angular_n_ion_producing_avalanches":x["n_ion_avalanches"]
        }
        need(abs(x["angel_ions"].sum()-x["ion_birth_total"])<
             1e-5*x["ion_birth_total"],"Ion weighted angle conservation failed")
        need(abs(x["angel_equal"].sum()-x["n_ion_avalanches"])<
             1e-5*x["n_ion_avalanches"],"Equal avalanche angular conservation failed")
    OUTPUT.mkdir(parents=True)
    phi= CENTRES
    with (OUTPUT/"paired_birth_angles_5deg.csv").open("x",newline="") as f:
        wr=csv.writer(f)
        wr.writerow(("phi_deg","Penning_ON_equal_avalanche",
                     "Penning_OFF_equal_avalanche",
                     "Penning_ON_ion_weighted",
                     "Penning_OFF_ion_weighted"))
        for j,v in enumerate(phi):
            wr.writerow((v,
                on["angel_equal"][j]/on["n_ion_avalanches"],
                off["angel_equal"][j]/off["n_ion_avalanches"],
                on["angel_ions"][j]/on["ion_birth_total"],
                off["angel_ions"][j]/off["ion_birth_total"]))
    with (OUTPUT/"paired_exact_parent_hemisphere.csv").open("x",newline="") as f:
        wr=csv.writer(f)
        wr.writerow(("parent_index","n_Heed","gain_on_e","gain_off_e",
                     "ion_births_on","ion_births_off","A_hemi_ON","A_hemi_OFF"))
        for i in range(1000):
            wr.writerow([i,int(on["muon_heed"][i]),
                int(on["muon_charge"][i]),int(off["muon_charge"][i]),
                int(on["muon_births"][i]),int(off["muon_births"][i]),
                on["A"][i],off["A"][i]])
    fig,axes=plt.subplots(1,2,figsize=(13,4.7))
    m0=np.isfinite(aoff);m1=np.isfinite(aon)
    edges=np.linspace(-1,1,51)
    axes[0].hist(aoff[m0],bins=edges,density=True,histtype="step",
                 label="Penning OFF",color="#737373",lw=1.7)
    axes[0].hist(aon[m1],bins=edges,density=True,histtype="step",
                 label="Penning ON",color="#e46b08",lw=1.7)
    axes[0].set_xlabel("Exact A_hemi per real J/psi muon")
    axes[0].set_ylabel("Density / 0.04 bin")
    axes[0].legend();axes[0].grid(alpha=.15)
    axes[1].scatter(aoff[defined],aon[defined],s=5,alpha=.3,c="#e46b08")
    axes[1].plot([-1,1],[-1,1],ls="--",color="#777",lw=1)
    axes[1].set_xlim(-1,1);axes[1].set_ylim(-1,1)
    axes[1].set_xlabel("A_hemi Penning OFF")
    axes[1].set_ylabel("A_hemi Penning ON")
    axes[1].grid(alpha=.15)
    fig.suptitle("Same 1,000 real J/psi muons: exact (unsampled) ion-birth hemisphere")
    fig.tight_layout()
    fig.savefig(OUTPUT/"paired_exact_hemisphere.png",dpi=180);plt.close(fig)
    for weight,title in (("angel_equal","Per ion-producing avalanche"),
                         ("angel_ions","Per created ion")):
        fig,ax=plt.subplots(figsize=(10.8,4.8))
        ax.step(phi,off[weight]/(off["n_ion_avalanches"] if weight=="angel_equal" else off["ion_birth_total"]),
                where="mid",color="#777",linestyle="--",lw=1.5,label="Penning OFF (matched)")
        ax.step(phi,on[weight]/(on["n_ion_avalanches"] if weight=="angel_equal" else on["ion_birth_total"]),
                where="mid",color="#e46b08",lw=1.8,label="Penning ON r=0.57")
        ax.set_xlabel("Ion-birth azimuth around closest 3mm anode [degrees]")
        ax.set_ylabel("Probability per 5-degree bin")
        ax.set_xlim(-180,180);ax.legend();ax.grid(alpha=.15)
        ax.set_title(f"Matched Garfield ion-birth angle — {title}\n"
                     "Sampled <=128 births/avalanche with exact count weights")
        fig.tight_layout()
        outname=("paired_angular_equal_avalanche.png" if weight=="angel_equal"
                 else "paired_angular_ion_weighted.png")
        fig.savefig(OUTPUT/outname,dpi=180);plt.close(fig)
    summary["source_sha256"]=BASE.sha(Path(__file__))
    summary["output_sha256"]={p.name:BASE.sha(p) for p in OUTPUT.iterdir()}
    BASE.jwrite(OUTPUT/"review_summary.json",summary)
    print("PASS 1000 matched true Geant4 parents; full per-chunk SHA audit")
    print("Microscopic Garfield avalanche calls:",len(on["charges"]))
    for key in ("off","on"):
        s=summary[key]["gain_collected_e"]
        a=summary[key]["A_ion_hemi_per_parent"]
        print(f"{key.upper()} G(all)={s['mean_all']:.2f}, "
              f"G(>0)={s['mean_positive']:.2f}, RMS={s['rms_all']:.2f}, "
              f"zero={s['p_zero']:.3%}; |A|={a['mean_abs']:.4f}, "
              f"sd(A)={a['std']:.4f}")
    print("Paired gain ON/OFF ratio",round(ratio,4),
          "parent bootstrap 95% interval",summary["parent_bootstrap"]["ratio_95ci"])
    print("Paired mean |A| ON-OFF",
          summary["comparison"]["mean_absolute_hemisphere_A_difference_on_minus_off"],
          "parent 95% interval",summary["parent_bootstrap"]["delta_mean_abs_A_95ci"])
    print("OUTPUT",OUTPUT)
def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("mode",choices=["review","status"])
    a=p.parse_args()
    if a.mode=="status":
        camp.status(full=True)
    else:
        create()
if __name__=="__main__":
    main()
