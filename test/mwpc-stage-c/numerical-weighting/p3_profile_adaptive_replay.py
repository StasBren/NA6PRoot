#!/usr/bin/env python3
"""S0/S1 development benchmark: profile fixed RKF paths and compare age-adaptive
time grids against native 10 ns cumulative-charge replay on disjoint subsets
of developmental muon 0301. Does not modify paths, LUT, or holdout archives.
No change to spatial interpolation, transport, or strip geometry.
"""
import argparse
import json
import time
from pathlib import Path
import numpy as np
from p3_replay_full_spectral_muon import (read_export, path_for, query, knots,
    waveform, report, WeightingLUT, EDGES, DT, NBINS, EFC, P, A)

def age_knots(t0,t1,ratio,max_step):
    # Build in event-clock edges. Preserve native bin endpoints within
    # the first 100ns of ion life and last 200ns before collection.
    last=min(t1,float(EDGES[-1]))
    if last<=t0:
        return np.array([t0,t1])
    a=max(0,int(np.floor(t0/DT))+1)
    b=min(NBINS,int(np.ceil(last/DT))-1)
    edges=EDGES[a:b+1]
    # Forward and backward logarithmic sampling, rounded to native edges.
    # Anchor exact birth, collection, and recorded window end.
    chosen=[t0,last,t1]
    if len(edges):
        chosen.extend(edges[(edges-t0<=100.)|(t1-edges<=200.)])
    for origin,sign,initial in ((t0,1.,100.),(t1,-1.,200.)):
        distance=initial
        duration=t1-t0
        # Cap far-time step; geometric grid near each endpoint.
        while distance<duration:
            t=origin+sign*distance
            if t0<t<last:
                e=DT*np.round(t/DT)
                if t0<e<last:
                    chosen.append(e)
            distance+=min(max_step,max(DT,ratio*distance))
    return np.unique(chosen)

def ramps(times,phi):
    """Slope-change accumulation of piecewise-linear cumulative charge.
    Same physics and same knot interpolation as existing waveform().
    Local time origin reduces cancellation in S0.
    """
    first=min(NBINS,max(0,int(np.floor(times[0]/DT))))
    stop=min(NBINS,max(first,int(np.ceil(times[-1]/DT))))
    shape=(stop-first,)+phi.shape[1:]
    if stop==first:
        return first,stop,np.zeros(shape)
    origin=first*DT
    tau=times-origin
    flat=phi.reshape(len(times),-1)
    slopes=np.diff(flat,axis=0)/np.diff(times)[:,None]
    delta=np.diff(np.vstack((np.zeros_like(slopes[:1]),slopes,
                            np.zeros_like(slopes[:1]))),axis=0)
    idx=np.ceil(tau/DT).astype(np.int64)
    valid=(idx>=0)&(idx<=stop-first)
    idx=idx[valid]; d=delta[valid]
    s1=np.zeros((stop-first+1,flat.shape[1]))
    s0=np.zeros_like(s1)
    np.add.at(s1,idx,d)
    np.add.at(s0,idx,-d*tau[valid,None])
    edges=np.arange(stop-first+1)*DT
    c=np.cumsum(s1,axis=0)*edges[:,None]+np.cumsum(s0,axis=0)
    return first,stop,np.diff(c,axis=0).reshape(shape)

def stats(waves,ref):
    out={}
    for side,label in enumerate(("U","V")):
        for end in (100,1000,10000,100000):
            stop=int(end/DT)
            a=waves[:stop,side]; b=ref[:stop,side]
            l2,mx=report(a,b)
            out[f"{label}_{end}ns"]=float(l2)
        # Channel-wise L2 error scaled to strongest reference channel,
        # avoiding ratios to nearly zero bipolar channels.
        d=np.linalg.norm(waves[:,side]-ref[:,side],axis=0)
        norm=max(np.max(np.linalg.norm(ref[:,side],axis=0)),1e-300)
        out[f"{label}_max_channel_error_over_strongest"]=float(np.max(d)/norm)
        for end in (10,30,100,1000,100000):
            k=int(end/DT)
            r=ref[:k,side].sum(axis=0)
            delta=waves[:k,side].sum(axis=0)-r
            out[f"{label}_Q{end}ns_abs_fC"]=float(np.max(abs(delta)))
    return out

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("table",type=Path)
    ap.add_argument("export",type=Path)
    ap.add_argument("output",type=Path)
    ap.add_argument("--rows",type=int,default=64)
    args=ap.parse_args()
    if args.rows<8 or args.rows%2:
        raise ValueError("rows must be even and >=8")
    started=time.perf_counter()
    rows,vertices,weights,exact=read_export(args.export)
    lut=WeightingLUT(args.table)
    print(f"load/census/LUT time={time.perf_counter()-started:.4f}s")
    if args.rows>len(rows):
        raise ValueError("sample too large")
    ids=(2*np.arange(args.rows)+1)*len(rows)//(2*args.rows)
    # Both subsets are developmental, not independent muon holdouts.
    # Interleave sample rows for similar coverage of archived primaries.
    offsets=np.arange(-5,6)
    configs=[("current100",None,None),("age025",.025,200.),
             ("age05",.05,400.),("age10",.10,800.)]
    names=["native10"]+[c[0] for c in configs]
    shape=(NBINS,2,len(offsets))
    sums={name:np.zeros((2,)+shape) for name in names}
    metrics={name:{"knots":[],"query_s":0.,"binning_s":0.,"grid_s":0.,
                   "individual_l2":[[],[]],"max_ramp_delta_Q_over_e":0.,
                   "max_ramp_closure_Q_over_e":0.} for name in names}
    path_s=0.; total_started=time.perf_counter()
    for i,j in enumerate(ids):
        j=int(j); split=i%2
        t=time.perf_counter()
        p=path_for(rows[j],vertices)
        path_s+=time.perf_counter()-t
        t0,t1=p[[0,-1],0]
        w=weights[j]*EFC
        reference=None
        for name,ratio,step in [("native10",None,None)]+configs:
            m=metrics[name]
            t=time.perf_counter()
            ts=(knots(t0,t1,1) if name=="native10" else
                knots(t0,t1,10) if name=="current100" else
                age_knots(t0,t1,ratio,step))
            m["grid_s"]+=time.perf_counter()-t
            t=time.perf_counter()
            phi=query(lut,p,ts,offsets)
            m["query_s"]+=time.perf_counter()-t
            m["knots"].append(len(ts))
            t=time.perf_counter()
            if name in ("native10","current100"):
                first,stop,q=waveform(ts,phi)
            else:
                first,stop,q=ramps(ts,phi)
            m["binning_s"]+=time.perf_counter()-t
            # Independently verify accumulation on identical knots.
            if name!="native10":
                a,b,check=(ramps(ts,phi) if name=="current100" else waveform(ts,phi))
                if (first,stop)!=(a,b):
                    raise RuntimeError("ramp time alignment failure")
                if q.size:
                    err=float(np.max(abs(q-check)))
                    close=float(np.max(abs(q.sum(axis=0)-check.sum(axis=0))))
                    m["max_ramp_delta_Q_over_e"]=max(m["max_ramp_delta_Q_over_e"],err)
                    m["max_ramp_closure_Q_over_e"]=max(m["max_ramp_closure_Q_over_e"],close)
                    if err>1e-9 or close>1e-8:
                        raise RuntimeError("ramp identity failed")
            sums[name][split,first:stop]+=w*q
            if name=="native10":
                reference=(first,stop,q)
            elif q.size:
                a,b,ref=reference
                if (a,b)!=(first,stop):
                    raise RuntimeError("native reference alignment failure")
                for side in range(2):
                    m["individual_l2"][split].append(float(report(q[:,side],ref[:,side])[0]))
        if (i+1)%8==0:
            print(f"development rows {i+1}/{args.rows}",flush=True)
    result={"sample_rows":ids.tolist(),"path_read_check_s":path_s,
            "benchmark_wall_s":time.perf_counter()-total_started,
            "vertices":int(len(vertices)),"vertex_file_MiB":(args.export/"vertices.f64").stat().st_size/2**20,
            "models":{}}
    for name in names:
        m=metrics[name]
        entry={k:v for k,v in m.items() if k not in ("knots","individual_l2")}
        entry["mean_knots"]=float(np.mean(m["knots"]))
        entry["max_knots"]=int(np.max(m["knots"]))
        entry["potential_queries"]=int(sum(m["knots"])*22)
        for s,label in enumerate(("tuning","dev_check")):
            e=stats(sums[name][s],sums["native10"][s])
            individual=m["individual_l2"][s]
            if individual:
                e["individual_p95"]=float(np.percentile(individual,95))
                e["individual_max"]=float(max(individual))
            entry[label]=e
            main_errors=[v for k,v in e.items() if k.endswith("ns") and "_Q" not in k]
            print(f"{name:12s} {label:10s} max_window_L2={100*max(main_errors):.6g}% "
                  f"U_full={100*e['U_100000ns']:.6g}% V_full={100*e['V_100000ns']:.6g}%")
        print(f"  mean_knots={entry['mean_knots']:.2f} query_s={m['query_s']:.4f} "
              f"binning_s={m['binning_s']:.4f} grid_s={m['grid_s']:.4f} "
              f"ramp_delta={m['max_ramp_delta_Q_over_e']:.3g}",flush=True)
        result["models"][name]=entry
    args.output.mkdir(parents=True,exist_ok=False)
    (args.output/"profile.json").write_text(json.dumps(result,indent=2)+"\n")
    with (args.output/"subset_waveforms.npz").open("xb") as f:
        np.savez(f,**sums)
    print("No production setting selected automatically. Compare both developmental subsets,"
          " per-channel and per-ion errors, timing, and closure before freezing.")
    print("Sealed muons untouched. RKF transport timing not included in this profile.")

if __name__=="__main__":
    main()
