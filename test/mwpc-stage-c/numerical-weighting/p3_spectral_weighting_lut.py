#!/usr/bin/env python3
"""Offline potential LUT for the ideal P3 nine-wire Fourier-MFS reference.

Reference settings: L=48 mm, kmax=6/mm, N=24 sources per wire.
Upper strip only; U/V symmetry and negative-u reflection are exact for this
ideal model. Near wire tabulate phi/log(r/R), bulk tabulate wire correction.
Use exact bare-strip solution (including cathode landing indicator) in bulk.
No approximate metal geometry and no new microscopic data.
This is a potential interpolation pilot; NOT the final flow-map replay bank.
"""
import argparse
import csv
import itertools
import json
import math
import time
from pathlib import Path
import numpy as np
from p3_fourier_mfs_mode_pilot import A, R, GAP, H, fit, slab
from p3_fourier_strip_reconstruction_pilot import bare, P

L=48.; KMAX=6.; N=24; IMAGE_EXTRA=4
PHASE_MAX=12.; NR_MAX=.75; BLEND0=.50; BLEND1=.65
PITCH=3.; UMAX=13.5

def mode_samples(x):
    """All *correction* Fourier coefficients at x, with cached solve per mode.
    Remove the known phase exp(-ik A u) so lookup uses eta=w-A u.
    Grounded metal extension: total potential is zero inside/at wire.
    """
    j=np.clip(np.rint(x[:,0]/PITCH),-4,4)
    inside=np.hypot(x[:,0]-j*PITCH,x[:,1])<=R*(1+1e-12)
    ks=2*np.pi/L*np.arange(int(np.floor(KMAX*L/(2*np.pi)))+1)
    result=np.empty((len(ks),len(x)),dtype=np.complex128)
    for m,k in enumerate(ks):
        images=0 if m==0 else max(8,int(np.ceil(24/(2*H*k))))+IMAGE_EXTRA
        fn,_=fit(k,N,images)
        scale=P/L if m==0 else 2*np.sin(k*P/2)/(k*L)
        for start in range(0,len(x),384):
            sl=slice(start,start+384)
            z=x[sl]
            b=slab(z,k)
            vals=-b.copy()
            gas=~inside[sl]
            if np.any(gas):
                vals[gas]=fn(z[gas])-b[gas]
            result[m,sl]=scale*vals*np.exp(1j*k*A*z[:,0])
        if m==0 or m==len(ks)-1 or m%5==0:
            print(f"  retained mode {m}/{len(ks)-1}, k={k:.6g}/mm",flush=True)
    return ks,result

def build(output,nr,nt,step,nphase):
    if output.exists():
        raise FileExistsError(output)
    if nr<12 or nt<16 or nphase<121 or step<=0 or step>.2:
        raise ValueError("insufficient table grid")
    nu=int(round(UMAX/step))+1; nv=int(round(2*GAP/step))+1
    if abs((nu-1)*step-UMAX)>1e-9 or abs((nv-1)*step-2*GAP)>1e-9:
        raise ValueError("bulk step must divide 13.5 and 6 mm exactly")
    eta=np.linspace(-PHASE_MAX,PHASE_MAX,nphase)
    logs=np.linspace(0.,np.log(NR_MAX/R),nr)
    theta=np.arange(nt)*2*np.pi/nt
    j,r,t=np.meshgrid(np.arange(5),R*np.exp(logs),theta,indexing="ij")
    near_x=np.column_stack(((PITCH*j+r*np.cos(t)).ravel(),(r*np.sin(t)).ravel()))
    u,v=np.meshgrid(np.linspace(0,UMAX,nu),np.linspace(-GAP,GAP,nv),indexing="ij")
    bulk_x=np.column_stack((u.ravel(),v.ravel()))
    x=np.vstack((near_x,bulk_x))
    print(f"BUILD {len(near_x)} near + {len(bulk_x)} bulk cross-section points;", 
          f"{nphase} strip phases",flush=True)
    started=time.perf_counter()
    ks,co=mode_samples(x)
    weights=np.exp(1j*ks[:,None]*eta[None,:])
    weights[1:]*=2
    corr=np.real(co.T@weights)
    n=len(near_x)
    # Use exact slab term at each point and phase. The near table stores
    # full phi, anchored at R. Correct zero at the grounded wire.
    full=corr[:n].copy()
    for start in range(0,n,256):
        stop=min(n,start+256)
        z=near_x[start:stop]
        zz=np.repeat(z,nphase,axis=0)
        ww=(A*z[:,0,None]+eta[None,:]).ravel()
        full[start:stop]+=bare(zz,ww).reshape(stop-start,nphase)
    near=full.reshape(5,nr,nt,nphase)
    # Grounded surface row is exactly zero by construction (not division
    # by zero). The normalized limit is extrapolated from the first two
    # nonzero log-r rows and tested directly at held-out near-wire points.
    near[:,0]=0.
    near[:,1:]/=logs[None,1:,None,None]
    near[:,0]=2*near[:,1]-near[:,2]
    bulk=corr[n:].reshape(nu,nv,nphase)
    bulk[:,0]=0.; bulk[:,-1]=0.
    if not np.all(np.isfinite(near)) or not np.all(np.isfinite(bulk)):
        raise RuntimeError("nonfinite table")
    meta={"model":"P3_9wire_infinite_plates_ideal_contiguous_strips",
          "L_mm":L,"cutoff_per_mm":KMAX,"sources":N,"image_extra":IMAGE_EXTRA,
          "phase_max_mm":PHASE_MAX,"radius_mm":R,"pitch_mm":PITCH,"strip_mm":P,
          "gap_mm":GAP,"stereo_tan":A,"near_max_mm":NR_MAX,
          "blend_mm":[BLEND0,BLEND1],"u_max_mm":UMAX,
          "near_shape":list(near.shape),"bulk_shape":list(bulk.shape),
          "field_accuracy_detector_validated":False}
    # Exclusive file creation; no historical data modified.
    with output.open("xb") as f:
        np.savez(f,near=near,bulk=bulk,meta=json.dumps(meta))
    print(f"BUILT {output} ({output.stat().st_size/2**20:.1f} MiB), "
          f"wall={time.perf_counter()-started:.2f}s",flush=True)

def stencil(t,n,periodic=False):
    start=np.floor(t).astype(np.int64)-1
    if not periodic:
        start=np.clip(start,0,n-4)
    x=t-start
    indices=start[:,None]+np.arange(4)
    if periodic:
        indices%=n
    weights=np.column_stack((-(x-1)*(x-2)*(x-3)/6,
                              x*(x-2)*(x-3)/2,
                             -x*(x-1)*(x-3)/2,
                              x*(x-1)*(x-2)/6))
    return indices,weights

def cubic(table,coords,periodic_axis=-1):
    """Local tensor-product cubic Lagrange interpolation, no spline fit."""
    if not len(coords[0]):
        return np.empty(0,dtype=float)
    st=[stencil(t,table.shape[d],d==periodic_axis)
        for d,t in enumerate(coords)]
    y=np.zeros(len(coords[0]))
    for i,j,k in itertools.product(range(4),repeat=3):
        y+=table[st[0][0][:,i],st[1][0][:,j],st[2][0][:,k]]*(
            st[0][1][:,i]*st[1][1][:,j]*st[2][1][:,k])
    return y

class WeightingLUT:
    def __init__(self,path):
        with np.load(path,allow_pickle=False) as f:
            self.near=f["near"]
            self.bulk=f["bulk"]
            self.meta=json.loads(f["meta"].item())
        for key,val in (("L_mm",L),("cutoff_per_mm",KMAX),("sources",N),
                        ("radius_mm",R),("strip_mm",P),("phase_max_mm",PHASE_MAX)):
            if self.meta[key]!=val:
                raise ValueError(f"table {key} mismatch")
        if list(self.near.shape)!=self.meta["near_shape"] or list(self.bulk.shape)!=self.meta["bulk_shape"]:
            raise ValueError("table shape mismatch")
    def upper(self,x,w):
        x=np.array(x,dtype=float,copy=True)
        w=np.broadcast_to(np.asarray(w,dtype=float),(len(x),)).copy()
        if x.shape!=(len(w),2) or not np.all(np.isfinite(x)) or not np.all(np.isfinite(w)):
            raise ValueError("bad query")
        flipped=x[:,0]<0
        x[flipped,0]*=-1
        w[flipped]*=-1
        eta=w-A*x[:,0]
        if np.any(x[:,0]>UMAX+1e-10) or np.any(abs(x[:,1])>GAP+1e-10):
            raise ValueError("query outside validated cross-section")
        if np.any(abs(eta)>PHASE_MAX+1e-10):
            raise ValueError("query outside declared strip phase range")
        j=np.clip(np.rint(x[:,0]/PITCH).astype(int),0,4)
        du=x[:,0]-j*PITCH
        r=np.hypot(du,x[:,1])
        phi=np.zeros(len(x))
        onwall=np.abs(x[:,1])>=GAP-1e-12
        if np.any(onwall):
            phi[onwall]=bare(x[onwall],w[onwall])
        gas=(r>R)&(~onwall)
        nsel=gas&(r<BLEND1)
        bsel=gas&(r>BLEND0)
        pn=np.zeros(len(x)); pb=np.zeros(len(x))
        nr,nt,ns=self.near.shape[1:]
        phase=(eta+PHASE_MAX)*(ns-1)/(2*PHASE_MAX)
        for wire in range(5):
            take=nsel&(j==wire)
            if not np.any(take):
                continue
            rr=np.log(r[take]/R)
            t=np.mod(np.arctan2(x[take,1],du[take]),2*np.pi)*nt/(2*np.pi)
            pn[take]=rr*cubic(self.near[wire],
                             [rr*(nr-1)/np.log(NR_MAX/R),t,phase[take]],
                             periodic_axis=1)
        if np.any(bsel):
            nu,nv,nb=self.bulk.shape
            pb[bsel]=bare(x[bsel],w[bsel])+cubic(self.bulk,
                [x[bsel,0]*(nu-1)/UMAX,
                 (x[bsel,1]+GAP)*(nv-1)/(2*GAP),
                 (eta[bsel]+PHASE_MAX)*(nb-1)/(2*PHASE_MAX)])
        z=np.clip((r-BLEND0)/(BLEND1-BLEND0),0.,1.)
        z=z*z*(3-2*z)
        phi[gas]=(1-z[gas])*pn[gas]+z[gas]*pb[gas]
        # Do not silently clip overshooting interpolants in this pilot.
        return phi

def read_paths(path):
    with path.open(newline="") as f:
        rows=sorted(csv.DictReader(f),key=lambda r:(int(r["row"]),float(r["age_ns"])))
    if len(rows)!=224:
        raise ValueError("expected 224 existing 16-ion trajectory points")
    return rows

def validate(path,paths,reference):
    lut=WeightingLUT(path)
    rows=read_paths(paths)
    with reference.open(newline="") as f:
        ref={(int(r["row"]),float(r["age_ns"]),r["side"],int(r["strip_offset"])):
             float(r["L96_N32_phi"]) for r in csv.DictReader(f)}
    if len(ref)!=2240:
        raise ValueError("unexpected reference count")
    x=np.array([[float(r["u_mm"]),float(r["v_mm"])] for r in rows])
    rp=np.empty((224,2,5)); pp=rp.copy()
    started=time.perf_counter()
    for side,label in enumerate(("U","V")):
        xx=x if side==0 else -x
        w=np.array([float(r["upper_w0_mm" if side==0 else "lower_w0_mm"]) for r in rows])
        for col,k in enumerate(range(-2,3)):
            pp[:,side,col]=lut.upper(xx,w-k*P)
            rp[:,side,col]=[ref[(int(row["row"]),float(row["age_ns"]),label,k)] for row in rows]
    elapsed=time.perf_counter()-started
    p=pp.reshape(16,14,2,5); q=rp.reshape(16,14,2,5)
    print(f"2240 tested potential queries (Python, incl. CSV comparisons): {elapsed:.6f}s")
    print(f"max_abs_phi_error={np.max(abs(p-q)):.9g}; sampled_range={p.min():.9g}..{p.max():.9g}")
    def report(label,pa,qa):
        diff=pa-qa
        rel=np.linalg.norm(diff)/max(np.linalg.norm(qa),1e-300)
        print(f"  {label:25s} relative_L2={100*rel:.6g}% max_abs_Q/e={np.max(abs(diff)):.9g}")
        return rel
    worst=0.
    for side,label in enumerate(("U","V")):
        worst=max(worst,report(label+" early 10ns bins",np.diff(p[:,:11,side],axis=1),
                              np.diff(q[:,:11,side],axis=1)))
        for idx,age in ((10,100),(11,1000),(12,10000),(13,100000)):
            worst=max(worst,report(f"{label} cumulative {age}ns",
                    p[:,idx,side]-p[:,0,side],q[:,idx,side]-q[:,0,side]))
    terminal=np.array([int(r["at_cathode"]) for r in rows],dtype=bool)
    if np.max(abs(pp[terminal]-rp[terminal]))>1e-8:
        raise RuntimeError("cathode endpoint mismatch")
    print("Provisional table pilot target: <0.5% L2 on each tested observable.")
    if worst>=.005 or p.min()<-1e-5 or p.max()>1+1e-5:
        print("REVIEW: interpolation pilot outside target")
        return 2
    print("PASS archived-path table pilot; independent spatial and full-waveform tests still required.")
    return 0

def main():
    ap=argparse.ArgumentParser()
    sub=ap.add_subparsers(dest="command",required=True)
    b=sub.add_parser("build")
    b.add_argument("output",type=Path)
    b.add_argument("--radial",type=int,default=32)
    b.add_argument("--angles",type=int,default=48)
    b.add_argument("--bulk-step",type=float,default=.15)
    b.add_argument("--phases",type=int,default=241)
    v=sub.add_parser("validate")
    v.add_argument("table",type=Path)
    v.add_argument("paths",type=Path)
    v.add_argument("reference",type=Path)
    args=ap.parse_args()
    if args.command=="build":
        build(args.output,args.radial,args.angles,args.bulk_step,args.phases)
        return 0
    return validate(args.table,args.paths,args.reference)
if __name__=="__main__":
    raise SystemExit(main())
