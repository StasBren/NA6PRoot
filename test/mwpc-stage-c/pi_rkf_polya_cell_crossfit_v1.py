#!/usr/bin/env python3
"""Read-only selected-cell crossfit of RKF + zero-inflated Polya versus A20.

Reuses sealed 16k matched RKF/microscopic starts. NO Garfield physics.
This is an exploratory count-only test, not a spatial cloud or cathode test.
Cells were A100-selected on older A20/120-muon information, so no claim of
independent event-level or flux-level validation.
"""
from __future__ import annotations
import argparse
import csv
import json
import math
from collections import defaultdict
import numpy as np
from scipy.stats import gamma

import pi_audit_rkf_v3_16k_status as audited
import pi_run_rkf_ion_census_full16k_v3 as v3
from pi_run_lut500_overnight import ROOT as A20

OUT = v3.BASE / "pi_RKF_Polya_cell_crossfit_v1"
Q = (.05, .5, .95, .99)

def must(ok, why):
    if not ok: raise ValueError(why)

def csvrows(path):
    with path.open(newline="") as f: return list(csv.DictReader(f))

def get_data():
    audited.check_archive_seal()  # source, binary, paired reference, manifests SHA
    m=json.loads((A20/"manifest.json").read_text())
    must(m.get("status")=="complete" and
         m.get("n_completed_avalanches")==84000, "Incomplete original A20")
    local=defaultdict(list)
    for g in v3.GEOMS:
        p=A20/("trials_"+g+".csv")
        must(v3.sha(p)==m["files_sha256"][p.name], "A20 SHA mismatch "+g)
        for r in csvrows(p):
            local[(g,int(r["box_id"]))].append((int(r["sample"]),int(r["ion_births"])))
    must(len(local)==4200 and all(len(v)==20 for v in local.values()),"A20 cell census")
    local={k:(np.array([y for x,y in sorted(v)],float),{x for x,y in v}) for k,v in local.items()}
    a100=json.loads((v3.A100/"manifest.json").read_text())
    truth={}
    for g in v3.GEOMS:
        p=v3.A100/("trials_extra_"+g+".csv")
        must(v3.sha(p)==a100["files_sha256"][p.name],"A100 SHA mismatch "+g)
        for r in csvrows(p):
            key=(g,int(r["box_id"]),int(r["sample"]))
            must(key not in truth,"Duplicated A100 micro start")
            truth[key]=int(r["ion_births"])
    must(len(truth)==16000,"Not all 16k real microscopic outcomes")
    data=defaultdict(list);found=set()
    for p in v3.read_csv(v3.PAIRS,v3.PAIRCOLS):
        g,b,s=p["geometry"],int(p["box_id"]),int(p["sample"])
        key=(g,b,s)
        must(key in truth and key not in found,"RKF/micro pair mismatch")
        found.add(key)
        must(int(p["micro_ion_births"])==truth[key],"Microscopic census changed")
        must(p["rkf_transport_ok"]=="1" and p["rkf_metal"]=="0",
             "Nontransport RKF case needs separate physical treatment "+str(key))
        rkf=float(p["rkf_ni"])
        must(math.isfinite(rkf) and rkf>=0,"Invalid RKF estimate")
        must((g,b) in local and s not in local[g,b][1],
             "Original A20 reused in independent new80")
        data[g,b].append((s,rkf,truth[key]))
    must(found==set(truth) and len(data)==200 and
         all(len(v)==80 for v in data.values()),"Incomplete independent cell data")
    data={k:(np.array([x for s,x,y in sorted(rows)],float),
             np.array([y for s,x,y in sorted(rows)],float))
          for k,rows in data.items()}
    return data,local

def fit(train,constant=False):
    x=np.concatenate([a for a,b in train])
    y=np.concatenate([b for a,b in train])
    if constant: x=np.ones_like(x)
    p0=float((np.sum(y==0)+.5)/(len(y)+1))
    c=float(np.sum(y)/np.sum(x)) if np.sum(x)>0 else float("nan")
    must(math.isfinite(c) and c>0,"Invalid mean-gain normalization")
    positive=y>0
    must(np.all(x[positive]>0),"RKF predicts zero at a real positive-ion start")
    z=y[positive]/(x[positive]*c/(1-p0))
    var=float(np.var(z,ddof=1))
    must(math.isfinite(var) and var>0,"Gamma fluctuation not identifiable")
    raw=float(np.mean(z)**2/var)
    shape=float(np.clip(raw,.1,100))
    return {"zero_prob":p0,"mean_ratio":c,"gamma_shape":shape,
            "raw_shape":raw,"n":int(len(y)),
            "note":"Gamma shape clipped to [.1,100] if needed"}

def quantiles(x,par):
    x=np.asarray(x,float)
    p0=par["zero_prob"]
    mean=x*par["mean_ratio"]
    must(np.all(np.isfinite(x)) and np.all(x>=0),"Invalid heldout RKF values")
    q={}
    for a in Q:
        if a<=p0:q[a]=np.zeros_like(x)
        else:
            z=(a-p0)/(1-p0)
            scale=np.maximum(mean,1e-30)/((1-p0)*par["gamma_shape"])
            q[a]=np.where(x>0,gamma.ppf(z,par["gamma_shape"],scale=scale),0.)
    return mean,q

def metric(y,mean,q,pzero):
    pin=lambda lev:float(np.mean(np.maximum(lev*(y-q[lev]),(lev-1)*(y-q[lev]))))
    return {
      "cell_mean_absolute_error":abs(float(y.mean()-mean.mean())),
      "one_draw_MAE":float(np.mean(np.abs(y-mean))),
      "zero_Brier":float(np.mean(((y==0).astype(float)-pzero)**2)),
      "PI90_coverage":float(np.mean((y>=q[.05])&(y<=q[.95]))),
      "upper95_exceedance":float(np.mean(y>q[.95])),
      "upper99_exceedance":float(np.mean(y>q[.99])),
      "q95_pinball":pin(.95),"q99_pinball":pin(.99),
    }

def bootstrap(vals,reps,rng):
    vals=np.array(vals,float)
    idx=rng.integers(len(vals),size=(reps,len(vals)))
    means=vals[idx].mean(axis=1)
    return {"value":float(vals.mean()),"CI95_cell_bootstrap":
            list(map(float,np.quantile(means,[.025,.975])))}

def analyze(reps,seed):
    must(not OUT.exists(),"Result already exists; refusing overwrite "+str(OUT))
    data,old=get_data()
    rng=np.random.default_rng(seed)
    results=[];folds=[]
    for g in v3.GEOMS:
        cells=np.array(sorted(b for h,b in data if h==g))
        must(len(cells)==50,"Expected 50 selected cells per geometry")
        rng.shuffle(cells)
        halves=(list(map(int,cells[:25])),list(map(int,cells[25:])))
        for fi in (0,1):
            train=halves[1-fi];test=halves[fi]
            params=fit([data[g,b] for b in train])
            control=fit([data[g,b] for b in train],constant=True)
            folds.append({"geometry":g,"fold":fi,"train":sorted(train),
                          "test":sorted(test),"RKF_Polya":params,
                          "geometry_only_Polya":control})
            for b in test:
                x,y=data[g,b];n=len(y)
                mean,q=quantiles(x,params)
                cm,cq=quantiles(np.ones(n),control)
                old20=old[g,b][0]
                baseline={a:np.full(n,np.quantile(old20,a)) for a in Q}
                rec={"geometry":g,"box_id":b,"n":n,
                     "micro_mean":float(y.mean()),
                     "pzero_micro":float(np.mean(y==0)),
                     "n_RKF_zero_microscopic_positive":int(np.sum((x==0)&(y>0))),
                     "RKF_Polya":metric(y,mean,q,params["zero_prob"]),
                     "geom_Polya":metric(y,cm,cq,control["zero_prob"]),
                     "A20":metric(y,np.full(n,old20.mean()),baseline,
                        float((np.sum(old20==0)+.5)/(len(old20)+1)))}
                results.append(rec)
    must(len(results)==200,"A cell was omitted in the crossfit")
    summary={}
    for g in v3.GEOMS:
        arr=[r for r in results if r["geometry"]==g]
        must(len(arr)==50,"Missing cell crossfit prediction "+g)
        s={"n_validation_cells":50,"n_validation_avalanches":4000}
        for tag in ("RKF_Polya","geom_Polya","A20"):
            s[tag]={name:bootstrap([r[tag][name] for r in arr],reps,rng)
                    for name in arr[0][tag]}
        for tag in ("geom_Polya","A20"):
            s["RKF_minus_"+tag]={
              name:bootstrap([r["RKF_Polya"][name]-r[tag][name]
                               for r in arr],reps,rng)
              for name in ("cell_mean_absolute_error","PI90_coverage",
                           "q95_pinball","q99_pinball","zero_Brier")}
        summary[g]=s
    payload={"schema":"NA60_RKF_POLYA_CELLS_V1","status":"complete",
             "total_selected_cells":200,"total_validation_avalanches":16000,
             "type":"Two-fold crossfit of 50 whole spatial cells per geometry",
             "warning":"SELECTED CELL INTERNAL validation only; no new Geant4 muon blind test",
             "warning2":"90% central intervals may overcover because pzero ~17%; check upper-tail pinball and missing rate",
             "limitation":"Only ion_births. Does not validate drift-time, avalanche angle, (u,v,w,t) ion births or signal.",
             "summary":summary}
    OUT.mkdir(parents=True)
    files={"summary.json":payload,"fold_parameters.json":folds,
           "per_cell_results.json":results}
    for n,obj in files.items():
        (OUT/n).write_text(json.dumps(obj,indent=2,sort_keys=True,allow_nan=False)+"\n")
    seal={"schema":"NA60_RKF_POLYA_CELLS_SHA_V1","source_SHA256":{
       str(p):v3.sha(p) for p in (v3.RECEIPT,v3.OUTPUT,v3.PAIRS,v3.SEAL,
                                 A20/"manifest.json",v3.A100/"manifest.json")},
       "output_SHA256":{n:v3.sha(OUT/n) for n in files},"seed":seed,
       "bootstrap":reps}
    (OUT/"manifest.json").write_text(json.dumps(seal,indent=2,sort_keys=True)+"\n")
    for g in v3.GEOMS:
        s=summary[g]
        print(g,"RKF+Polya/A20/geomOnly",
              "cell MAE",*[round(s[t]["cell_mean_absolute_error"]["value"],2)
                            for t in ("RKF_Polya","A20","geom_Polya")],
              "PI90%",*[round(100*s[t]["PI90_coverage"]["value"],2)
                          for t in ("RKF_Polya","A20","geom_Polya")],
              "q95 pinball",*[round(s[t]["q95_pinball"]["value"],2)
                              for t in ("RKF_Polya","A20","geom_Polya")],flush=True)
    print("PASS immutable physical inputs, 200-cell crossfit, SHA output",OUT,flush=True)

def selftest():
    rng=np.random.default_rng(17)
    train=[]
    for _ in range(25):
        x=rng.uniform(100,500,80)
        y=np.where(rng.random(80)<.17,0,
                   rng.gamma(1.5,x*(.8/.83)/1.5))
        train.append((x,y))
    p=fit(train);m,q=quantiles(np.array([100.,500.]),p)
    must(.12<p["zero_prob"]<.22 and .7<p["mean_ratio"]<.9,
         "Synthetic fit wrong")
    must(np.all(q[.95]>q[.5]) and np.all(q[.5]>=q[.05]),
         "Invalid predictive quantiles")
    scores=metric(train[0][1],*quantiles(train[0][0],p),p["zero_prob"])
    must(0<=scores["PI90_coverage"]<=1 and scores["q95_pinball"]>=0,
         "Invalid scoring")
    print("PASS synthetic zero-atom+Polya fit and 95/99 quantiles; NO physics",flush=True)

def main():
    a=argparse.ArgumentParser(description=__doc__)
    a.add_argument("mode",choices=("selftest","audit"))
    a.add_argument("--bootstrap",type=int,default=3000)
    a.add_argument("--seed",type=int,default=20261104)
    z=a.parse_args()
    must(200<=z.bootstrap<=10000,"Use 200..10000 bootstrap replicas")
    if z.mode=="selftest":selftest()
    else:analyze(z.bootstrap,z.seed)
if __name__=="__main__":main()
