#!/usr/bin/env python3
"""Numerical sensitivity of U/V induced charge for 16 frozen archive rows.
Ion ages (not event-clock bins). Five strips per side. Exact ideal endpoint
when collected. Independent field verification and fast operator still pending.
"""
import argparse
import csv
from pathlib import Path
import numpy as np
from p3_fourier_strip_reconstruction_pilot import Correction, P, A

AGES=np.array([*range(0,101,10),1000,10000,100000],dtype=float)
OFFSETS=np.arange(-2,3)

def compare(a,b):
    d=b-a
    return (np.linalg.norm(d)/max(np.linalg.norm(b),1e-300),
            np.max(np.abs(d)))

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("paths_csv",type=Path)
    ap.add_argument("output_dir",type=Path)
    args=ap.parse_args()
    with args.paths_csv.open(newline="") as f:
        rows=list(csv.DictReader(f))
    ids=sorted({int(r["row"]) for r in rows})
    if len(ids)!=16 or len(rows)!=16*len(AGES):
        raise ValueError("expected exactly 16 sampled archived ions, 14 ages each")
    rows.sort(key=lambda r:(int(r["row"]),float(r["age_ns"])))
    if not np.array_equal(np.array([float(r["age_ns"]) for r in rows]),
                          np.tile(AGES,16)):
        raise ValueError("age grid mismatch")
    x=np.array([[float(r["u_mm"]),float(r["v_mm"])] for r in rows])
    w=np.array([[float(r["upper_w0_mm"]),float(r["lower_w0_mm"])] for r in rows])
    xx=np.vstack((x,-x))
    if not np.all(np.isfinite(xx)) or not np.all(np.isfinite(w)):
        raise ValueError("nonfinite path")
    if np.any(abs(xx[:,1])>3.+1e-12):
        raise ValueError("path outside cathodes")
    phase=np.r_[w[:,0]-A*x[:,0],w[:,1]+A*x[:,0]]
    if np.max(abs(phase))+2*P>=19.:
        raise ValueError("shifted strip phase too far for this period test")
    args.output_dir.mkdir(exist_ok=False)
    configs=[("L48_N24",48.,6.,24),("L96_N24",96.,6.,24),
             ("L96_N32",96.,6.,32)]
    predicts={}
    for name,L,K,N in configs:
        print(f"RUN {name}: L={L:g} cutoff={K:g} N={N}",flush=True)
        model=Correction(xx,L,K,N,4)
        phi=np.stack([model.phi(np.r_[w[:,0],w[:,1]]-j*P)
                      for j in OFFSETS],axis=-1)
        phi=np.stack([phi[:len(rows)],phi[len(rows):]],axis=1)
        phi=phi.reshape(16,len(AGES),2,5)
        if not np.all(np.isfinite(phi)):
            raise RuntimeError("nonfinite weighting potential")
        print(f"  sampled range={phi.min():.9g}..{phi.max():.9g}",flush=True)
        if phi.min() < -1e-5 or phi.max()>1+1e-5:
            raise RuntimeError("sampled maximum-principle violation")
        # Exact ideal electrode endpoint, independently computed from strip phase.
        terminal=np.array([int(r["at_cathode"]) for r in rows],dtype=bool)
        for i in np.flatnonzero(terminal):
            for side in range(2):
                sign=1 if side==0 else -1
                eta=w[i,side]-sign*A*x[i,0]-OFFSETS*P
                if np.any(abs(abs(eta)-P/2)<1e-7):
                    raise RuntimeError("ambiguous strip-edge landing")
                expect=(abs(eta)<P/2).astype(float) if x[i,1]*sign>0 else np.zeros(5)
                if np.max(abs(phi.reshape(len(rows),2,5)[i,side]-expect))>1e-7:
                    raise RuntimeError("terminal strip BC mismatch")
        predicts[name]=phi
        del model
    print("Numerical convergence, per unit ion charge, unweighted samples:")
    for a,b in zip(configs[:-1],configs[1:]):
        p,q=predicts[a[0]],predicts[b[0]]
        pc=p-p[:,0:1]; qc=q-q[:,0:1]
        for side,label in enumerate(("U","V")):
            vals=[("10ns bins first100ns",np.diff(p[:,:11,side],axis=1),
                   np.diff(q[:,:11,side],axis=1))]
            for age in (100.,1000.,10000.,100000.):
                idx=int(np.flatnonzero(AGES==age)[0])
                vals.append((f"Q({age:g}ns)",pc[:,idx,side],qc[:,idx,side]))
            for what,pa,qb in vals:
                rel,mx=compare(pa,qb)
                print(f"  {a[0]}->{b[0]} {label} {what:22s} "
                      f"relL2={100*rel:.6g}% max_abs_Q/e={mx:.8g}")
    with (args.output_dir/"charge_samples.csv").open("x",newline="") as f:
        writer=csv.writer(f)
        names=[c[0] for c in configs]
        writer.writerow(["row","primary","weight","birth_ns","age_ns",
                         "side","strip_offset","strip_index","end_age_ns"]+
                        [name+"_phi" for name in names]+
                        [name+"_Q_over_e" for name in names])
        for i,r in enumerate(rows):
            trial,age=divmod(i,len(AGES))
            for side,label in enumerate(("U","V")):
                for j,offset in enumerate(OFFSETS):
                    base=int(r["upper_index" if side==0 else "lower_index"])
                    ph=[predicts[n][trial,age,side,j] for n in names]
                    ch=[predicts[n][trial,age,side,j]-predicts[n][trial,0,side,j]
                        for n in names]
                    writer.writerow([r[k] for k in ("row","primary","weight",
                                     "birth_ns","age_ns")]+
                                    [label,offset,base+offset,r["end_age_ns"]]+ph+ch)
    print("COMPLETED numerical sensitivity pilot only; not a full archive,"
          " waveform validation or independent electrostatic certification.")

if __name__=="__main__": main()
