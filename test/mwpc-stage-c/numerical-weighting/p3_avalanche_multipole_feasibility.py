#!/usr/bin/env python3
"""T1: test local 2D grounded-cylinder moments against the *direct* 3D
Fourier-MFS weighting reference, not a fitted spatial LUT.

All lengths mm. Fixed longitudinal coordinate on each ring, with phase
xi = w - A*u_wire - strip_center; NOT w-A*u at the moving evaluation point.
Positive wire classes 0..4 suffice for this symmetric ideal geometry.
Fit angular coefficients on an anchor circle; independently test other
radii and offset angles. Also report same-radius angular truncation to
distinguish longitudinal/radial-model error from angular truncation.
No archive read, no transport, no front-end model, no holdout access.
"""
import argparse
import csv
import json
import time
from pathlib import Path
import numpy as np
from p3_spectral_weighting_lut import mode_samples, L, KMAX, N, A, R, PITCH
from p3_fourier_strip_reconstruction_pilot import bare

RADII=np.array([1.01*R,1.1*R,1.5*R,2*R,.05,.075,.1,.15,.2,.3,.45,.6,.75,1.,1.2,1.4])
ANCHORS=np.array([.1,.2,.3,.45,.6,.75,1.,1.2,1.4])
PHASES=np.linspace(-11.5,11.5,49)
NTH=48

def direct_reference():
    # Train and test angular grids do not share nodes. Phase is held
    # fixed relative to the wire center for each angular projection.
    angles=(np.arange(NTH)[None,:]+np.array([[0.],[.371]]))*2*np.pi/NTH
    j,rad,ang=np.broadcast_arrays(np.arange(5)[:,None,None,None],
        RADII[None,:,None,None],angles[None,None,:,:])
    du=rad*np.cos(ang)
    x=np.column_stack(((j*PITCH+du).ravel(),(rad*np.sin(ang)).ravel()))
    eta_shift=-A*du.ravel()
    if np.max(np.abs(eta_shift))+np.max(abs(PHASES))>12.:
        raise RuntimeError("outside spectral budget phase range")
    print(f"Direct reference: {len(x)} transverse points x {len(PHASES)} phases",
          flush=True)
    ks,co=mode_samples(x)
    # mode_samples removes exp(-ik A*u); restore local ring dependence.
    co*=np.exp(1j*ks[:,None]*eta_shift[None,:])
    factors=np.exp(1j*ks[:,None]*PHASES[None,:])
    factors[1:]*=2
    result=np.real(co.T@factors)
    # Restore isolated bare strip (not periodic). w fixed on a ring.
    for i in range(0,len(x),256):
        z=x[i:i+256]
        center=j.ravel()[i:i+len(z)]*PITCH
        w=(A*center[:,None]+PHASES[None,:]).ravel()
        result[i:i+len(z)]+=bare(np.repeat(z,len(PHASES),axis=0),w).reshape(len(z),-1)
    if not np.all(np.isfinite(result)):
        raise RuntimeError("nonfinite reference")
    if result.min()<-1e-6 or result.max()>1+1e-6:
        raise RuntimeError("reference maximum principle excursion")
    ref=result.reshape(5,len(RADII),2,NTH,len(PHASES))
    return angles,ref

def basis(theta,order):
    return np.stack([np.ones(len(theta))]+[
        f(n*theta) for n in range(1,order+1) for f in (np.cos,np.sin)],axis=-1)

def radial_ratio(r,b,order):
    # Normalization at anchor radius avoids ill-scaled powers of mm.
    vals=[np.log(r/R)/np.log(b/R)]
    for n in range(1,order+1):
        v=( (r/b)**n * (-np.expm1(-2*n*np.log(r/R))) /
                      (-np.expm1(-2*n*np.log(b/R))) )
        vals.extend([v,v])
    return np.asarray(vals)

def metrics(pred,ref):
    d=pred-ref
    return (float(np.max(abs(d))),float(np.linalg.norm(d)/
            max(np.linalg.norm(ref),1e-300)))

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("out",type=Path)
    args=ap.parse_args()
    if args.out.exists():
        raise FileExistsError(args.out)
    started=time.perf_counter()
    theta,ref=direct_reference()
    print(f"Reference range [{ref.min():.9g}, {ref.max():.9g}]",flush=True)
    args.out.mkdir(parents=True,exist_ok=False)
    all_rows=[]; summary=[]
    train=ref[:,:,0]; test=ref[:,:,1]
    print("Errors against independent offset angles; relL2 uses all wires/angles/phases per shell.")
    for order in range(5):
        bt=basis(theta[0],order)
        bv=basis(theta[1],order)
        coeff=np.einsum("ma,jrap->jrmp",np.linalg.pinv(bt),train)
        # Angular-only fit at each radius, no imposed radial dependence.
        angular=np.einsum("am,jrmp->jrap",bv,coeff)
        for idx,r in enumerate(RADII):
            ma,l2=metrics(angular[:,idx],test[:,idx])
            all_rows.append(["angular_only",order,r,r,ma,l2])
        for anchor in ANCHORS:
            anchor_idx=int(np.flatnonzero(RADII==anchor)[0])
            c=coeff[:,anchor_idx]
            use=np.flatnonzero(RADII<=anchor)
            shell_errors=[]
            for idx in use:
                r=RADII[idx]
                factor=radial_ratio(r,anchor,order)
                pred=np.einsum("am,jmp,m->jap",bv,c,factor)
                ma,l2=metrics(pred,test[:,idx])
                all_rows.append(["grounded_2d",order,anchor,r,ma,l2])
                shell_errors.append((ma,l2))
            worst_abs=max(e[0] for e in shell_errors)
            worst_l2=max(e[1] for e in shell_errors)
            # Screening only, NOT a detector waveform or physics gate.
            screened=worst_abs<1e-4 and worst_l2<.005
            summary.append(dict(order=order,rank=2*order+1,anchor_mm=float(anchor),
                 max_abs=worst_abs,max_shell_relL2=worst_l2,screen_pass=screened))
            print(f"order={order} rank={2*order+1} rmax={anchor:.3f}mm "
                  f"max_abs={worst_abs:.6g} max_shell_L2={100*worst_l2:.6g}% "
                  f"screen={'PASS' if screened else 'REVIEW'}")
    with (args.out/"radial_errors.csv").open("x",newline="") as f:
        wr=csv.writer(f)
        wr.writerow(["model","order","anchor_mm","test_radius_mm","max_abs_phi","relative_L2"])
        wr.writerows(all_rows)
    with (args.out/"reference_samples.npz").open("xb") as f:
        np.savez(f,phi=ref,radii_mm=RADII,theta_rad=theta,phase_mm=PHASES,
                 wire_classes=np.arange(5))
    meta={"reference":{"period_mm":L,"cutoff_per_mm":KMAX,"sources":N},
          "phase_convention":"w-A*u_wire-strip_center",
          "screen_abs":1e-4,"screen_relL2":.005,
          "wall_s":time.perf_counter()-started,"models":summary}
    (args.out/"summary.json").write_text(json.dumps(meta,indent=2)+"\n")
    print("Maximum radial/angle error screening only. Phase interpolation,"
          " waveform/charge error, longitudinal ion spread, and runtime untested.")
    print("Angular-only errors in CSV separate angular truncation from the 2D radial-law assumption.")
    print("No holdout access; no additional FEM or microscopic simulation.")

if __name__=="__main__":
    main()
