#!/usr/bin/env python3
"""Pilot: isolated strip in an infinite parallel-plate, finite nine-wire model.
Exact wire-free slab + periodic Fourier approximation to the wire correction.
The artificial longitudinal period, mode cutoff, source count and image count
are numerical controls, not hardware parameters. No archived particle data.
"""
import argparse
import csv
from pathlib import Path
import numpy as np
from p3_fourier_mfs_mode_pilot import A, R, GAP, H, ring, fit, slab

P = 1.7  # strip pitch AND width in w at fixed u, mm

def bare(x, w):
    """Infinite isolated upper strip, phase w - A*u, exact slab Poisson integral."""
    x = np.asarray(x)
    y = x[:, 1] + GAP
    phase = np.asarray(w) - A*x[:, 0]
    scale = np.sqrt(1+A*A)
    d = P/2
    z = np.pi*y/(2*H)
    t = np.tan(z)
    hi = np.pi*(phase+d)/(2*H*scale)
    lo = np.pi*(phase-d)/(2*H*scale)
    val = (np.arctan(t*np.tanh(hi))-np.arctan(t*np.tanh(lo)))/np.pi
    val = np.where(y==0.,0.,val)
    val = np.where(y==H,np.where(abs(phase)<d,1.,
                               np.where(abs(phase)==d,.5,0.)),val)
    return val

class Correction:
    def __init__(self, x, period, cutoff, n, extra_images):
        self.x=np.asarray(x)
        self.period=period
        dk=2*np.pi/period
        self.ks=dk*np.arange(1,int(np.floor(cutoff/dk))+1)
        f,_=fit(0.,n,0)
        self.zero=(f(self.x)-slab(self.x,0.)).real*P/period
        self.coeff=[]
        for k in self.ks:
            images=max(8,int(np.ceil(24/(2*H*k))))+extra_images
            f,_=fit(k,n,images)
            # Fourier coefficient of strip indicator, periodic over L.
            coefficient=2*np.sin(k*P/2)/(k*period)
            self.coeff.append(coefficient*(f(self.x)-slab(self.x,k)))
        self.coeff=np.array(self.coeff)
    def phi(self,w):
        w=np.broadcast_to(np.asarray(w,dtype=float),(len(self.x),))
        return bare(self.x,w)+self.zero+2*np.real(np.sum(
            self.coeff*np.exp(1j*self.ks[:,None]*w[None,:]),axis=0))

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("output", type=Path)
    args=ap.parse_args()
    args.output.mkdir(parents=True,exist_ok=False)
    obs=[]
    labels=[]
    for j in (-4,0,4):
        for ratio in (1.01,1.1,2.,5.,20.):
            for theta in (.37,1.94,3.51,5.08):
                obs.append([3*j+R*ratio*np.cos(theta),R*ratio*np.sin(theta)])
                labels.append(f"wire{j}_r{ratio}")
    for u in (-12.,0.,12.):
        for v in (-2.95,-1.,1.,2.95):
            obs.append([u,v])
            labels.append("bulk")
    # Independent circular boundary points, offset from all source/collocation grids.
    wires=ring(19,R,.281)
    walls=np.array([[u,v] for u in (-12.,0.,12.) for v in (-GAP,GAP)])
    x=np.vstack([np.array(obs),wires,walls])
    no=len(obs); nw=len(wires)
    # Values around both strip edges and far from the selected strip.
    phase=np.array([-5.,-1.7,-.86,-.4,0.,.4,.86,1.7,5.])
    ws=[A*x[:,0]+s for s in phase]
    # Vary one parameter at a time, starting with an inexpensive base.
    configs=[
        ("base",24.,4.,16,0),
        ("angular",24.,4.,24,0),
        ("images",24.,4.,24,4),
        ("cutoff",24.,6.,24,4),
        ("period",48.,6.,24,4),
    ]
    pred={}
    final=None
    for name,L,kmax,n,images in configs:
        print(f"RUN {name}: L={L:g}mm cutoff={kmax:g}/mm sources={n} "
              f"extra_images={images}",flush=True)
        model=Correction(x,L,kmax,n,images)
        val=np.stack([model.phi(w) for w in ws],axis=1)
        if not np.all(np.isfinite(val)):
            raise RuntimeError("nonfinite reconstructed strip potential")
        pred[name]=val
        wall_expected=np.stack([bare(x,w) for w in ws],axis=1)[no+nw:]
        print(f"  modes={len(model.ks)} wire_max_abs="
              f"{np.max(abs(val[no:no+nw])):.8g} "
              f"wall_max_abs={np.max(abs(val[no+nw:]-wall_expected)):.8g} "
              f"gas_range=[{val[:no].min():.8g},{val[:no].max():.8g}]",
              flush=True)
        final=model
    print("Independent parameter sensitivity (gas points only):")
    for a,b in zip(configs[:-1],configs[1:]):
        old=pred[a[0]][:no]; new=pred[b[0]][:no]; diff=new-old
        print(f"  {a[0]}->{b[0]} max_abs={np.max(abs(diff)):.8g} "
              f"relL2={np.linalg.norm(diff)/np.linalg.norm(new):.8g}")
    # Sum of 17 contiguous individual, translated strips; finite sum misses tails.
    # Last model has L=48mm; sampled phases plus +/-13.6mm remain well within L/2.
    f,_=fit(0.,24,0)
    phi0=f(x[:no]).real
    total=np.zeros(len(x))
    for j in range(-8,9):
        total+=final.phi(A*x[:,0]-j*P)
    missing=phi0-total[:no]
    print(f"17-strip sum vs k=0: max_abs_difference={np.max(abs(missing)):.8g} "
          f"min_signed_difference={missing.min():.8g} "
          f"relative_L2={np.linalg.norm(missing)/np.linalg.norm(phi0):.8g}")
    with (args.output/"strip_samples.csv").open("x",newline="") as f:
        wr=csv.writer(f)
        wr.writerow(["region","u_mm","v_mm","w_mm","phase_mm"]+
                    [c[0] for c in configs])
        for i in range(no):
            for j,s in enumerate(phase):
                wr.writerow([labels[i],*x[i],ws[j][i],s]+
                            [pred[c[0]][i,j] for c in configs])
    print("COMPLETED PILOT ONLY: sampled upper strip, not production certification.")
    print("Need independent reference for nonzero modes, tighter parameter test,"
          " and full trajectory/multi-strip validation.")
if __name__=="__main__":
    main()
