#!/usr/bin/env python3
"""Physics-only strict-parent holdout of conditional p3_noFW avalanche models.
TRAIN: 10,000 additional REAL J/psi Geant4/TrackHeed/Garfield parents.
TEST: original 1,000 disjoint REAL J/psi parents.
Compare local expectation, *joint* empirical (ion_up,ion_down) resampling,
and independent zero-inflated Gamma gain + empirical individual A.
No new Garfield, no ion cloud coordinates, no cathode signals. Read-only.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,math
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO=Path(__file__).resolve().parents[2]
BASE=REPO/"test_runs/mwpc_stage_c"
TRAIN=BASE/"pi_p3_noFW_jpsi_extra10000_native_v1"
TEST=BASE/"pi_p3_noFW_jpsi1000_native_v1"
OUT=BASE/"pi_p3_noFW_conditional_gain_asym_v1"
FIRST="source|chunk|event|hit_index|track_id|station|det_id"
SOURCE=("source","chunk","event","hit_index","track_id","station","det_id","u_in_mm")
SEED=("crossing_id","configuration","primary_index","seed_u_mm","seed_v_mm",
      "ion_births","ion_up","ion_down","ion_boundary","capped")
EVENT=(FIRST,"configuration","heed_primaries","ion_births","ion_up","ion_down","ion_boundary")
NU,NV=8,12

def must(ok,msg):
    if not ok:raise RuntimeError(msg)
def sha(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda:f.read(1<<20),b""):h.update(b)
    return h.hexdigest()
def load(p,need):
    with p.open(newline="") as f:
        reader=csv.DictReader(f)
        must(set(need).issubset(reader.fieldnames or []),"Bad columns: "+str(p))
        yield from reader
def idof(r):
    return "|".join(r[k] for k in SOURCE[:7])
def parent(r):
    return "|".join(r[k] for k in SOURCE[:3])
def round_cpp(x):
    return math.floor(x+.5) if x>=0 else math.ceil(x-.5)
def entryphase(u,st):
    n=math.floor((307.2 if st==0 else 512.)/3.+.5)
    p=1.5 if n%2==0 else 0.
    return u-p-3*round_cpp((u-p)/3)
def cells(u,v):
    w=u/3.
    k=np.where(w>=0,np.floor(w+.5),np.ceil(w-.5))
    local=u-3*k
    bu=np.clip(np.floor((local+1.5)*NU/3).astype(int),0,NU-1)
    bv=np.clip(np.floor((v+3)*NV/6).astype(int),0,NV-1)
    return bu*NV+bv

def campaign(root,n):
    meta=json.loads((root/"input/input_manifest.json").read_text())
    src=root/"input/native_muon_crossings_input.csv"
    must(meta["n_parents"]==n and meta["geometry"]=="p3_noFW"
         and meta["group_size"]==20 and sha(src)==meta["selected_input_sha256"],
         "Not a frozen 20-parent p3-noFW campaign")
    original=list(load(src,SOURCE))
    ids=[idof(r) for r in original]
    parents={parent(r) for r in original}
    must(len(ids)==len(set(ids))==len(parents)==n,"Duplicate Geant4 parent")
    idx={x:i for i,x in enumerate(ids)}
    phase=np.asarray([entryphase(float(r["u_in_mm"]),int(r["station"]))
                      for r in original])
    up=np.zeros(n,dtype=np.int64);down=np.zeros(n,dtype=np.int64)
    heed=np.zeros(n,dtype=np.int64)
    seen=np.zeros(n,dtype=bool);si=[];su=[];sv=[];nu=[];nd=[]
    receipts=[]
    for ch in range(n//20):
        d=root/"native_chunks"/f"chunk_{ch:04d}"
        receipt=json.loads((d/"complete.json").read_text())
        must(receipt["chunk"]==ch and receipt["n_parents"]==20,
             "Bad SHA-sealed 20-muon checkpoint "+str(d))
        layout=receipt.get("recovery_layout")
        must(layout in (None,"high_heed_replay_v1"),"Unknown replay layout")
        phys=d/layout if layout else d
        for name in ("native_muon_crossings_input.csv",
                     "muon_reference_crossing.csv","muon_reference_seed.csv"):
            must(sha(phys/name)==receipt["file_sha256"][name],
                 "Corrupted frozen Garfield CSV: "+str(phys/name))
        receipts.append([ch,sha(d/"complete.json")])
        got=[idof(r) for r in load(phys/"native_muon_crossings_input.csv",SOURCE)]
        must(got==ids[20*ch:20*(ch+1)],"Geant4 source identity changed")
        for r in load(phys/"muon_reference_crossing.csv",EVENT):
            must(r["configuration"]=="p3_noFW","Wrong physical geometry")
            i=idx[r[FIRST]]
            U,D,B=(int(r[k]) for k in ("ion_up","ion_down","ion_boundary"))
            must(not seen[i] and U+D+B==int(r["ion_births"]),
                 "Muon complete ion births invalid or duplicate")
            up[i],down[i],heed[i]=U,D,int(r["heed_primaries"])
            seen[i]=True
        for r in load(phys/"muon_reference_seed.csv",SEED):
            must(r["configuration"]=="p3_noFW" and int(r["capped"])==0,
                 "Unexpected geometry or capped physical avalanche")
            i=idx[r["crossing_id"]]
            U,D,B=(int(r[k]) for k in ("ion_up","ion_down","ion_boundary"))
            must(U+D+B==int(r["ion_births"]),"Invalid full per-seed ion census")
            si.append(i);su.append(float(r["seed_u_mm"]))
            sv.append(float(r["seed_v_mm"]));nu.append(U);nd.append(D)
    must(seen.all(),"Incomplete real parent physics census")
    i=np.asarray(si,dtype=np.int32)
    U=np.asarray(nu,dtype=np.float64);D=np.asarray(nd,dtype=np.float64)
    must(np.array_equal(np.bincount(i,minlength=n),heed),"Heed seed census mismatch")
    must(np.array_equal(np.bincount(i,weights=U,minlength=n).astype(np.int64),up)
         and np.array_equal(np.bincount(i,weights=D,minlength=n).astype(np.int64),down),
         "True whole-muon sector counts differ from per-electron Garfield")
    su=np.asarray(su);sv=np.asarray(sv)
    must(np.isfinite(su).all() and np.isfinite(sv).all() and (np.abs(sv)<3.02).all(),
         "Invalid starting electron positions")
    print("PASS",root.name,n,"Geant4 parents,",len(i),"full microscopic avalanches",flush=True)
    return dict(parents=parents,idx=i,U=U,D=D,cell=cells(su,sv),
                phase=phase,up=up,down=down,ids=ids,receipts=receipts,
                source_sha=sha(src))

def train_cells(tr, minN=80):
    label=tr["cell"];gain=tr["U"]+tr["D"]
    bank=[];diagnostics=[]
    for u in range(NU):
      for v in range(NV):
        c=u*NV+v
        for rad in range(max(NU,NV)+1):
            grid=[a*NV+b for a in range(max(0,u-rad),min(NU,u+rad+1))
                         for b in range(max(0,v-rad),min(NV,v+rad+1))]
            pool=np.flatnonzero(np.isin(label,grid))
            if len(pool)>=minN:break
        must(len(pool)>=minN,"Insufficient real Garfield donors")
        positive=pool[gain[pool]>0]
        must(len(positive)>0,"All avalanches have zero positive ions")
        g=gain[positive];mn=float(g.mean());var=float(g.var())
        shape=np.clip(mn**2/var,1e-3,1e5) if var>0 else 1e5
        a=(tr["U"][positive]-tr["D"][positive])/g
        rho=float(np.corrcoef(g,a)[0,1]) if len(g)>2 and var>0 and np.std(a)>0 else None
        bank.append((pool,positive,mn,shape,1-len(positive)/len(pool)))
        diagnostics.append(dict(u_bin=u,v_bin=v,
           n_direct=int(np.sum(label==c)),n_pool=len(pool),radius=rad,
           p_zero=1-len(positive)/len(pool),gamma_shape=float(shape),
           gamma_mean_positive=mn,empirical_gain_asym_corr=rho))
    return bank,diagnostics

def simulate(hold,tr,bank,R,seed):
    rng=np.random.default_rng(seed)
    ix=hold["idx"];labels=hold["cell"];N=len(hold["ids"])
    U=tr["U"];D=tr["D"];G=U+D
    valid=(hold["up"]+hold["down"])>0
    true=np.divide(hold["up"]-hold["down"],hold["up"]+hold["down"],
                   out=np.zeros(N),where=valid)
    # Nonstochastic conditional mean is a deliberately weak baseline.
    mU=np.zeros(len(ix));mD=np.zeros(len(ix))
    for c,(pool,positive,mn,k,p0) in enumerate(bank):
        j=np.flatnonzero(labels==c)
        mU[j]=np.mean(U[pool]);mD[j]=np.mean(D[pool])
    up0=np.bincount(ix,weights=mU,minlength=N)
    dn0=np.bincount(ix,weights=mD,minlength=N)
    mean=np.divide(up0-dn0,up0+dn0,out=np.zeros(N),where=up0+dn0>0)
    models={"paired_empirical":[],"factorized_gamma":[]}
    for rep in range(R):
        u=np.zeros(len(ix));d=np.zeros(len(ix))
        fg=np.zeros(len(ix));fa=np.zeros(len(ix))
        for c,(pool,positive,mn,k,p0) in enumerate(bank):
            j=np.flatnonzero(labels==c)
            if not len(j):continue
            picks=pool[rng.integers(len(pool),size=len(j))]
            u[j]=U[picks];d[j]=D[picks]
            active=rng.random(len(j))>=p0
            gain=rng.gamma(k,mn/k,size=len(j))*active
            angle=positive[rng.integers(len(positive),size=len(j))]
            fg[j]=gain
            fa[j]=(U[angle]-D[angle])/G[angle]
        ou=np.bincount(ix,weights=u,minlength=N)
        od=np.bincount(ix,weights=d,minlength=N)
        pa=np.divide(ou-od,ou+od,out=np.zeros(N),where=ou+od>0)
        gg=np.bincount(ix,weights=fg,minlength=N)
        da=np.bincount(ix,weights=fg*fa,minlength=N)
        ga=np.divide(da,gg,out=np.zeros(N),where=gg>0)
        models["paired_empirical"].append(pa)
        models["factorized_gamma"].append(ga)
    return true,valid,mean,{x:np.asarray(y) for x,y in models.items()}

def summary(a):
    return dict(mean=float(a.mean()),sd=float(a.std(ddof=1)),abs_mean=float(np.abs(a).mean()),
                p_abs_gt_05=float(np.mean(np.abs(a)>.5)),
                p_abs_gt_08=float(np.mean(np.abs(a)>.8)))
def phase_bin(phi,a,valid,lo,hi):
    j=valid&(phi>=lo)&(phi<hi)
    return (float(np.mean(np.abs(a[j]))) if np.sum(j)>0 else None,int(j.sum()))
def render(truth,valid,mean,models,phi):
    edges=np.linspace(-1,1,51)
    fig,axes=plt.subplots(1,2,figsize=(12,4.6))
    for x,label in [(truth,"Full Garfield truth"),(mean,"Only E[U,D]"),
                    (models["paired_empirical"][0],"Paired empirical (replica 1)"),
                    (models["factorized_gamma"][0],"Independent Gamma+A (replica 1)")]:
        axes[0].hist(x[valid],bins=edges,density=True,histtype="step",lw=1.5,label=label)
    axes[0].set_xlabel("Muon whole-event A_hemi");axes[0].set_ylabel("Density")
    axes[0].legend(fontsize=8);axes[0].grid(alpha=.15)
    ph_edges=np.linspace(-1.5,1.5,13)
    centres=(ph_edges[:-1]+ph_edges[1:])/2
    for x,label in [(truth,"True Garfield"),(mean,"Conditional mean"),
                    (np.mean(np.abs(models["paired_empirical"]),axis=0),"Paired empirical |A| mean"),
                    (np.mean(np.abs(models["factorized_gamma"]),axis=0),"Factorized |A| mean")]:
        if label in ("Factorized |A| mean","Paired empirical |A| mean"):
            vals=[float(np.mean(x[valid&(phi>=lo)&(phi<hi)]))
                  for lo,hi in zip(ph_edges[:-1],ph_edges[1:])]
        else:
            vals=[phase_bin(phi,x,valid,lo,hi)[0] for lo,hi in zip(ph_edges[:-1],ph_edges[1:])]
        axes[1].plot(centres,vals,marker="o",ms=3,label=label)
    axes[1].set_xlabel("True muon entrance phase relative anode [mm]")
    axes[1].set_ylabel("Mean |A_hemi|");axes[1].legend(fontsize=8);axes[1].grid(alpha=.15)
    fig.suptitle("Real J/psi parents: 10k training / strictly disjoint 1k holdout")
    fig.tight_layout();fig.savefig(OUT/"strict_heldout_conditional_gain_asym.png",dpi=180);plt.close(fig)

def selftest():
    assert round_cpp(.5)==1 and round_cpp(-.5)==-1
    assert entryphase(0,0)==1.5 and entryphase(0,1)==0.
    assert cells(np.array([0.,3.,-3.,1.5]),np.zeros(4)).min()>=0
    assert cells(np.array([0.,3.,-3.,1.5]),np.zeros(4)).max()<96
    assert summary(np.array([-.5,.5]))["mean"]==0.
    print("PASS C++ wire-grid parity, zero-safe exact A, 96 conditioning bins")
    print("No Garfield calls; no source or output files touched.")

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("mode",choices=("selftest","analyze"))
    p.add_argument("--replicas",type=int,default=24)
    p.add_argument("--seed",type=int,default=20261108)
    a=p.parse_args()
    if a.mode=="selftest":selftest();return
    must(not OUT.exists(),"Refusing to overwrite immutable analysis: "+str(OUT))
    must(a.replicas>=3,"Use at least 3 independent whole-muon surrogate replicas")
    tr=campaign(TRAIN,10000);hold=campaign(TEST,1000)
    must(tr["parents"].isdisjoint(hold["parents"]),"TRAIN/TEST original Geant4 parent leakage")
    bank,diagnostics=train_cells(tr)
    truth,valid,mean,models=simulate(hold,tr,bank,a.replicas,a.seed)
    report=dict(schema="P3_NOFW_CONDITIONAL_FULL_GARFIELD_HOLDOUT_V1",
      train_parents=10000,holdout_parents=1000,
      train_avalanches=len(tr["idx"]),holdout_avalanches=len(hold["idx"]),
      true=summary(truth[valid]),
      mean_only_r2=1-np.mean((truth[valid]-mean[valid])**2)/np.var(truth[valid]),
      source_sha256=sha(Path(__file__)),replicas=a.replicas,seed=a.seed,
      train_checkpoint_sha=tr["receipts"],holdout_checkpoint_sha=hold["receipts"],
      limitations="Ion BIRTH hemispheres only, not Gaussian full clouds or cathode signal",
      models={},phase_bins=[])
    for name,x in models.items():
        stats=[summary(a[valid]) for a in x]
        report["models"][name]={
            k:dict(mean=float(np.mean([v[k] for v in stats])),
                   between_replica_sd=float(np.std([v[k] for v in stats],ddof=1)),
                   bias_vs_true=float(np.mean([v[k] for v in stats])-report["true"][k]))
            for k in stats[0]}
    phase_edges=np.linspace(-1.5,1.5,13)
    for lo,hi in zip(phase_edges[:-1],phase_edges[1:]):
        v,n=phase_bin(hold["phase"],truth,valid,lo,hi)
        row=dict(lo=float(lo),hi=float(hi),n=n,true_mean_abs_A=v)
        for name,data in models.items():
            row[name+"_mean_abs_A"]=float(np.mean([
                phase_bin(hold["phase"],x,valid,lo,hi)[0] for x in data]))
        report["phase_bins"].append(row)
    OUT.mkdir(parents=True,exist_ok=False)
    render(truth,valid,mean,models,hold["phase"])
    (OUT/"spatial_conditioning_cells.json").write_text(json.dumps(diagnostics,indent=2)+"\n")
    report["output_sha256"]={p.name:sha(p) for p in OUT.iterdir()}
    (OUT/"summary.json").write_text(json.dumps(report,indent=2)+"\n")
    print("PASS disjoint true-muon surrogate comparison",flush=True)
    print("TRUE",report["true"],flush=True)
    for name,st in report["models"].items():print(name,st,flush=True)
    print("OUTPUT",OUT,flush=True)
if __name__=="__main__":main()
