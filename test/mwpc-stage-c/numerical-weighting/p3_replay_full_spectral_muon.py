#!/usr/bin/env python3
"""Complete weighted muon, direct RKF trajectories + spectral potential LUT.
No trajectory-bank interpolation yet. Time-grid compression is tested against
native 10 ns event-clock endpoints on a deterministic independent row subset.
The output remains provisional ion-only, ideal strip geometry/mobility.
"""
import argparse
import csv
import json
import time
from pathlib import Path
import numpy as np
from p3_spectral_weighting_lut import WeightingLUT, P, A, R
from p3_fourier_mfs_mode_pilot import fit

EFC=1.602176634e-4
DT=10.
NBINS=10000
EDGES=np.arange(NBINS+1)*DT

def read_export(directory):
    seal=dict(line.strip().split("=",1) for line in
              (directory/"complete.txt").read_text().splitlines() if line.strip())
    with (directory/"ions.csv").open(newline="") as f:
        rows=list(csv.DictReader(f))
    n=int(seal["rows"])
    if len(rows)!=n or n==0 or float(seal["step_mm"])!=.005:
        raise ValueError("invalid archive seal")
    vv=int(seal["vertices"])
    v=np.memmap(directory/"vertices.f64",dtype="<f8",mode="r",shape=(vv,5))
    if (directory/"vertices.f64").stat().st_size!=40*vv:
        raise ValueError("vertex size mismatch")
    offset=0
    for j,r in enumerate(rows):
        if (int(r["row"])!=j or int(r["first_vertex"])!=offset or
                int(r["vertices"])<2):
            raise ValueError("invalid vertex index")
        offset+=int(r["vertices"])
    weights=np.array([float(r["weight"]) for r in rows])
    if offset!=vv or np.any(weights<=0) or not np.all(np.isfinite(weights)):
        raise ValueError("archive census mismatch")
    exact=int(seal["exact_ions"])
    if not np.isclose(weights.sum(),exact,rtol=1e-9,atol=1e-6):
        raise ValueError("archive weight census mismatch")
    return rows,v,weights,exact

def path_for(r,v):
    start=int(r["first_vertex"]); stop=start+int(r["vertices"])
    p=np.asarray(v[start:stop])
    if (not np.all(np.isfinite(p)) or np.any(np.diff(p[:,0])<=0) or
            abs(p[0,0]-float(r["birth_ns"]))>1e-7 or
            abs(p[-1,0]-float(r["end_ns"]))>1e-7 or p[0,0]<0 or
            abs(abs(p[-1,2])-3.)>1e-12):
        raise ValueError("invalid path at row "+r["row"])
    if np.max(np.abs(p[:,3:]-p[0,3:]))>1e-7:
        raise ValueError("unexpected longitudinal motion")
    rho=np.hypot(p[:,1]-3*np.rint(p[:,1]/3),p[:,2])
    if np.any(rho<=R) or np.any(abs(p[:,2])>3+1e-10):
        raise ValueError("non-gas point in archived RKF path")
    return p

def query(lut,p,t,offsets):
    x=np.column_stack([np.interp(t,p[:,0],p[:,i]) for i in (1,2)])
    out=[]
    for side in range(2):
        xx=x if side==0 else -x
        w=p[0,3+side]-P*offsets
        phi=lut.upper(np.repeat(xx,len(offsets),axis=0),
                      np.tile(w,len(t))).reshape(len(t),len(offsets))
        out.append(phi)
    a=np.stack(out,axis=1)
    if not np.all(np.isfinite(a)) or a.min()<-1e-5 or a.max()>1+1e-5:
        raise ValueError("LUT maximum-principle excursion beyond tolerance")
    return a

def knots(t0,t1,stride):
    # Native event-clock edges around birth and collection. Coarser
    # absolute-clock edges in between. Keep exact birth and end.
    left=max(0,int(np.floor(t0/DT))+1)
    right=min(NBINS,int(np.ceil(t1/DT))-1)
    e=EDGES[left:right+1] if right>=left else np.empty(0)
    select=((np.arange(left,right+1)%stride==0)|
            (e-t0<=1000.)|(t1-e<=1000.)) if len(e) else np.empty(0,dtype=bool)
    # If the ion extends past the recorded window, include window end;
    # the endpoint is kept separately to calculate the unrecorded tail.
    return np.unique(np.r_[t0,e[select],t1])

def waveform(times,phi):
    # Integrate Q by differences at the original event-clock bin endpoints.
    # Piecewise linear cumulative potential, NOT shifted current templates.
    # Only active time bins are materialized.
    first=min(NBINS,max(0,int(np.floor(times[0]/DT))))
    stop=min(NBINS,max(first,int(np.ceil(times[-1]/DT))))
    edges=EDGES[first:stop+1]
    shape=(stop-first,)+phi.shape[1:]
    if stop==first:
        return first,stop,np.zeros(shape)
    flat=phi.reshape(len(times),-1)
    q=flat-flat[0]
    ee=np.column_stack([np.interp(edges,times,q[:,j],left=0.,right=q[-1,j])
                        for j in range(q.shape[1])])
    return first,stop,np.diff(ee,axis=0).reshape(shape)

def report(a,b):
    d=a-b
    return float(np.linalg.norm(d)/max(np.linalg.norm(b),1e-300)),float(np.max(abs(d)))

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("table",type=Path)
    ap.add_argument("export",type=Path)
    ap.add_argument("out",type=Path)
    ap.add_argument("--half-strips",type=int,default=5)
    ap.add_argument("--late-step-ns",type=int,default=100)
    ap.add_argument("--reference-rows",type=int,default=32)
    args=ap.parse_args()
    if not 2<=args.half_strips<=6 or not 10<=args.late_step_ns<=200 or args.late_step_ns%10:
        raise ValueError("invalid strip or time controls")
    rows,vertices,weights,exact=read_export(args.export)
    n=len(rows)
    if not 1<=args.reference_rows<=n:
        raise ValueError("invalid reference sample")
    offsets=np.arange(-args.half_strips,args.half_strips+1)
    nc=len(offsets)
    bases=np.array([int(rows[0]["upper_index"]),int(rows[0]["lower_index"])])
    subset=set(((2*np.arange(args.reference_rows)+1)*n//(2*args.reference_rows)).tolist())
    lut=WeightingLUT(args.table)
    stride=args.late_step_ns//10
    args.out.mkdir(parents=True,exist_ok=False)
    total=np.zeros((NBINS,2,nc))
    coarse_sample=np.zeros_like(total); fine_sample=np.zeros_like(total)
    all_endpoint=np.zeros((2,nc)); independent_endpoint=np.zeros_like(all_endpoint)
    window_check=np.zeros_like(all_endpoint)
    landing=np.zeros_like(all_endpoint)
    shape=(2,nc)
    sample_stats=[]
    check_points=[]; check_sums=[]
    total_queries=0; late_rows=0; late_weight=0.; started=time.perf_counter()
    for j,(r,weight) in enumerate(zip(rows,weights)):
        if [int(r["upper_index"]),int(r["lower_index"])]!=bases.tolist():
            raise ValueError("inconsistent physical strip frame")
        p=path_for(r,vertices)
        t0,t1=p[[0,-1],0]; q=weight*EFC
        t=knots(t0,t1,stride)
        phi=query(lut,p,t,offsets)
        total_queries+=len(t)*2*nc
        # Independently assign ideal contiguous landing electrode.
        side=0 if p[-1,2]>0 else 1
        eta=p[-1,3+side]-(1 if side==0 else -1)*A*p[-1,1]
        k=int(np.floor(eta/P+.5))
        if abs(abs(eta-k*P)-P/2)<1e-7:
            raise ValueError("ambiguous ideal landing at strip edge")
        if k not in offsets:
            raise ValueError("selected channels omit the landing strip; widen or recenter")
        endpoint=np.zeros(shape); endpoint[side,k+args.half_strips]=1.
        if np.max(abs(phi[-1]-endpoint))>1e-8:
            raise ValueError("LUT endpoint does not match physical strip")
        landing+=q*endpoint
        all_endpoint+=q*(phi[-1]-phi[0])
        independent_endpoint+=q*(endpoint-phi[0])
        win=np.clip(EDGES[-1],t0,t1)
        flat=phi.reshape(len(t),-1)
        at_win=np.array([np.interp(win,t,flat[:,ch]) for ch in range(2*nc)]).reshape(shape)
        window_check+=q*(at_win-phi[0])
        a,b,wave=waveform(t,phi)
        total[a:b]+=q*wave
        if t1>EDGES[-1]:
            late_rows+=1;late_weight+=weight
        if j in subset:
            tf=knots(t0,t1,1)
            pf=query(lut,p,tf,offsets)
            total_queries+=len(tf)*2*nc
            af,bf,wf=waveform(tf,pf)
            if (af,bf)!=(a,b):
                raise RuntimeError("reference window alignment")
            coarse_sample[a:b]+=q*wave
            fine_sample[a:b]+=q*wf
            for si in range(2):
                if len(wave):
                    l2,mx=report(wave[:,si],wf[:,si])
                    sample_stats.append((j,si,l2,mx))
            # Independent k=0 strip-sum check at a few path points.
            tx=np.unique(np.clip(t0+np.array([0.,100.,1000.,10000.,t1-t0]),t0,t1))
            val=query(lut,p,tx,offsets)
            xx=np.column_stack([np.interp(tx,p[:,0],p[:,i]) for i in (1,2)])
            check_points.append(xx);check_sums.append(val.sum(axis=2))
        if (j+1)%128==0 or j+1==n:
            print(f"replay {j+1}/{n} rows, queries={total_queries}, "
                  f"elapsed={time.perf_counter()-started:.1f}s",flush=True)
    window=total.sum(axis=0)
    closure=float(np.max(abs(window-window_check)))
    endpoint_error=float(np.max(abs(all_endpoint-independent_endpoint)))
    if closure>1e-7 or endpoint_error>1e-7:
        raise RuntimeError(f"charge conservation failure: {closure}, {endpoint_error}")
    if abs(landing.sum()-exact*EFC)>1e-7:
        raise RuntimeError("collected ion census mismatch")
    # Same numerical reference, independently solved normalized k=0 mode.
    fn,_=fit(0.,24,0)
    xx=np.concatenate(check_points)
    sums=np.concatenate(check_sums)
    plane=np.column_stack((fn(xx).real,fn(-xx).real))
    missing=plane-sums
    strip_tail_max=float(np.max(np.abs(missing)))
    print(f"subset whole-plane minus selected-strip sum: max_abs={strip_tail_max:.8g}, "
          f"min_signed={missing.min():.8g} (finite channel truncation + lookup error)")
    print("Time compression test: same RKF + LUT, weighted subset")
    errors={}
    for side,label in enumerate(("U","V")):
        for end in (100,1000,10000,100000):
            k=int(end/DT)
            rel,mx=report(coarse_sample[:k,side],fine_sample[:k,side])
            errors[f"{label}_{end}ns"]=rel
            print(f"  {label} 0..{end:6d}ns: L2={100*rel:.6g}% max_abs_bin_fC={mx:.8g}")
    sample_array=np.array(sample_stats)
    if len(sample_array):
        print(f"Individual-trajectory full-window L2: median={100*np.median(sample_array[:,2]):.6g}% "
              f"p95={100*np.percentile(sample_array[:,2],95):.6g}% "
              f"max={100*sample_array[:,2].max():.6g}%")
    # Explicit separate outputs; no hidden rescaling of charge or tails.
    labels=[f"{side}_{base+int(k)}" for side,base in zip(("U","V"),bases) for k in offsets]
    np.savetxt(args.out/"waveforms.csv",
               np.column_stack((EDGES[:-1]+DT/2,total.reshape(NBINS,-1)/DT)),
               delimiter=",",header="time_center_ns,"+",".join(labels),comments="")
    with (args.out/"channels.csv").open("x",newline="") as f:
        wr=csv.writer(f);wr.writerow(["side","strip_index","offset","window_charge_fC",
            "endpoint_charge_fC","unrecorded_tail_fC","collected_ion_charge_fC",
            "first100ns_charge_fC","waveform_L2_fC_per_ns"])
        for si,side in enumerate(("U","V")):
            for ci,k in enumerate(offsets):
                wr.writerow([side,int(bases[si]+k),int(k),window[si,ci],
                    all_endpoint[si,ci],all_endpoint[si,ci]-window[si,ci],
                    landing[si,ci],total[:10,si,ci].sum(),
                    np.linalg.norm(total[:,si,ci]/DT)])
    with (args.out/"subset_comparison.npz").open("xb") as f:
        np.savez(f,coarse=coarse_sample,fine=fine_sample,rows=sorted(subset),
                 sample_stats=sample_array)
    worst=max(errors.values())
    summary={"rows":n,"exact_ions":exact,"weight_sum":float(weights.sum()),
             "late_rows":late_rows,"late_weight":late_weight,"n_channels":len(labels),
             "queries":total_queries,"wall_s":time.perf_counter()-started,
             "dt_ns":DT,"late_step_ns":args.late_step_ns,"subset_rows":len(subset),
             "subset_relative_L2":errors,"window_closure_fC":closure,
             "endpoint_closure_fC":endpoint_error,"max_subset_strip_sum_gap":strip_tail_max,
             "time_grid_pilot_pass":bool(worst<.005),"field_accuracy_detector_validated":False,
             "electron_current_included":False,"strip_origin_assumed_zero":True,
             "mobility_provisional":True}
    (args.out/"summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    print(json.dumps(summary,indent=2))
    print("Provisional timing gate: <0.5% weighted subset L2 for every reported window.")
    if worst>=.005:
        print("REVIEW: refine late-step-ns in a new output directory; preserve paths.")
        return 2
    print("PASS complete weighted muon replay timing pilot; not final fast model or physical certification.")
    return 0
if __name__=="__main__":
    raise SystemExit(main())
