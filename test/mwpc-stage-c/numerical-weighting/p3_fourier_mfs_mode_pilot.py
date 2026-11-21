#!/usr/bin/env python3
"""Finite nine-wire, infinite parallel-plate Fourier-mode MFS pilot.
Coordinates mm. Normalized upper cathode mode exp(-i*k*0.1*u).
No finite sidewalls, full inverse Fourier quadrature, ion transport, or
strip-response validation. Sources are inside metal, not physical gas.
Requires numpy and scipy. Angular and image truncations still need
independent refinement studies and external benchmarks.
"""
import numpy as np
from scipy.special import k0
from scipy.linalg import solve

R=.015; GAP=3.; H=2*GAP; PITCH=3.; A=.1
CENTERS=np.arange(-4,5)*PITCH

def ring(n,r,offset=0.):
    t=(np.arange(n)+offset)*2*np.pi/n
    return np.column_stack(((CENTERS[:,None]+r*np.cos(t)).ravel(),
        np.broadcast_to(r*np.sin(t),(9,n)).ravel()))

def green(x,s,k,images):
    dx=x[:,0,None]-s[None,:,0]
    y=x[:,1,None]+GAP; ys=s[None,:,1]+GAP
    dy=y-ys; sy=y+ys
    # Dirichlet strip Green function of -Laplace+k^2.
    if k==0:
        d=np.sinh(np.pi*dx/(2*H))**2
        return np.log((d+np.sin(np.pi*sy/(2*H))**2)/
                      (d+np.sin(np.pi*dy/(2*H))**2))/(4*np.pi)
    z=np.zeros_like(dx)
    for n in range(-images,images+1):
        z+=k0(k*np.hypot(dx,dy+2*n*H))-k0(k*np.hypot(dx,sy+2*n*H))
    return z/(2*np.pi)

def slab(x,k):
    q=abs(k)*np.hypot(1,A)
    v=x[:,1]; u=x[:,0]
    z=(v+GAP)/H if q==0 else np.exp(-q*(GAP-v))*(
        -np.expm1(-2*q*(v+GAP)))/(-np.expm1(-2*q*H))
    return z*np.exp(-1j*k*A*u)

def fit(k,n,images):
    src=ring(n,.65*R); col=ring(n,R)
    mat=green(col,src,k,images)
    c=solve(mat,-slab(col,k))
    return lambda x: slab(x,k)+green(x,src,k,images)@c, np.linalg.cond(mat)

def main():
    check=ring(73,R,.371)
    obs=np.vstack([ring(13,r,.27) for r in
                   [1.0066667*R,1.1*R,2*R,5*R,20*R]])
    walls=np.array([[u,v] for u in np.linspace(-13.4,13.4,41)
                    for v in [-GAP,GAP]])
    for k in [0.,.25,1.,3.]:
        images=0 if k==0 else max(8,int(np.ceil(24/(2*H*k))))
        f1,c1=fit(k,16,images)
        f2,c2=fit(k,24,images+4 if k else 0)
        a,b=f1(obs),f2(obs)
        wall_err=np.max(abs(f2(walls)-slab(walls,k)))
        print(f'k={k:g} /mm, N=16->24 sources/wire, '
              f'independent wire max residual={np.max(abs(f1(check))):.3e}'
              f'->{np.max(abs(f2(check))):.3e}; '
              f'observed mode delta={np.max(abs(a-b)):.3e}; '
              f'relative L2={np.linalg.norm(a-b)/np.linalg.norm(b):.3e}; '
              f'cathode correction={wall_err:.3e}; cond={c2:.3e}',flush=True)
    print("Mode-level numerical pilot only; no full strip or physical accuracy certification.")

if __name__=='__main__': main()
