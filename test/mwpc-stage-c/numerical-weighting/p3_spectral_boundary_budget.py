#!/usr/bin/env python3
"""Final pre-table diagnostic: sampled boundary error budget of Fourier-MFS.
No FEM, no ion transport, no archive writes.

The periodic full solution has periodic bare-strip boundary data. Our
isolated-strip evaluator uses the isolated bare strip and periodic wire
correction. On an exact periodic solution the maximum principle bounds
the potential error by the wire boundary residual, and the periodic copy
correction can be bounded by the corresponding wire-free strip tails.
IMPORTANT: wire angle maxima here are sampled, not certified continuous
suprema. The resulting total is a diagnostic estimate, not a rigorous
global error certificate. The Fourier and copy tail formulas are analytic
bounds for the stated ideal geometry.
"""
import argparse
import math
import numpy as np
from p3_fourier_mfs_mode_pilot import A, R, GAP, H, ring, fit
from p3_fourier_strip_reconstruction_pilot import P

def cutoff_tail(L, first_k, dk):
    # |f_n| <= P/L; upper slab factor at actual wire surface
    # <= exp(-q*(GAP-R))/(1-exp(-2*q*H)).
    # The omitted periodic full solution has a wire correction equal
    # to minus the omitted slab modes on the grounded wires.
    scale=math.hypot(1,A)
    d=scale*(GAP-R)
    return (2*P/L*math.exp(-first_k*d)/
            ((-math.expm1(-dk*d))*(-math.expm1(-2*first_k*scale*H))))

def copy_tail(L, phase_max):
    # For a nonselected periodic copy the wire-induced correction has
    # magnitude <= the corresponding positive wire-free strip potential.
    # Upper slab Poisson kernel <= 1/[2H(cosh(pi*dist/H)-1)].
    # Bound both positive and negative copies by geometric series.
    scale=math.hypot(1,A)
    width=P/scale
    d1=math.pi*(L-phase_max-P/2)/(H*scale)
    step=math.pi*L/(H*scale)
    if d1<=0:
        raise ValueError("phase range overlaps adjacent periodic copies")
    return (2*width/H*math.exp(-d1)/
            ((-math.expm1(-d1))**2*(-math.expm1(-step))))

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--period",type=float,default=48.)
    p.add_argument("--cutoff",type=float,default=6.)
    p.add_argument("--sources",type=int,default=24)
    p.add_argument("--angles",type=int,default=128)
    p.add_argument("--phase-max",type=float,default=12.)
    args=p.parse_args()
    L,K,N,na,S=args.period,args.cutoff,args.sources,args.angles,args.phase_max
    if not (L>0 and K>0 and N>=8 and na>=64 and 0<=S<L/2):
        raise ValueError("invalid numerical controls")
    dk=2*math.pi/L
    nm=int(math.floor(K/dk))
    # Two independently offset angular grids, all nine wires.
    # Interleaving grids detects sensitivity of sampled maxima.
    points=np.vstack([ring(na,R,.173),ring(na,R,.673)])
    m=9*na
    sums=np.zeros(2)
    maxmode=0.
    for i in range(nm+1):
        k=i*dk
        images=0 if i==0 else max(8,int(math.ceil(24/(2*H*k))))+4
        f,cond=fit(k,N,images)
        residual=np.abs(f(points))
        if not np.all(np.isfinite(residual)):
            raise RuntimeError("nonfinite mode boundary residual")
        coeff=P/L if i==0 else 4*abs(math.sin(k*P/2))/(k*L)
        v=np.array([np.max(residual[:m]),np.max(residual[m:])])
        sums+=coeff*v
        maxmode=max(maxmode,float(v.max()))
        if i==0 or i==nm or i%10==0:
            print(f"mode {i:3d}/{nm} k={k:.5g} "
                  f"max_sampled_residual={v.max():.8g} cond={cond:.6g}",flush=True)
    tail=cutoff_tail(L,(nm+1)*dk,dk)
    copies=copy_tail(L,S)
    sampled=float(max(sums))
    budget=sampled+tail+copies
    print(f"period={L:g}mm sources={N} cutoff={K:g}/mm "
          f"sampled_angles_per_wire={2*na} phase_limit={S:g}mm")
    print(f"max_normalized_mode_boundary_residual={maxmode:.10g}")
    print(f"weighted_mode_residual_sum_grid_A={sums[0]:.10g}")
    print(f"weighted_mode_residual_sum_grid_B={sums[1]:.10g}")
    print(f"analytic_omitted_fourier_tail_bound={tail:.10g}")
    print(f"analytic_periodic_copy_bound={copies:.10g}")
    print(f"estimated_absolute_potential_error_budget={budget:.10g}")
    print(f"estimated_absolute_bin_charge_error_per_e={2*budget:.10g}")
    print("Sampled angular residual estimate, NOT a certified continuous supremum.")
    print("Applies to ideal infinite parallel plates, nine wires, and declared phase range.")
    print("Provisional pre-table diagnostic target: 1e-6 absolute potential error.")
    if budget>1e-6:
        print("REVIEW: estimated budget exceeds provisional target")
        return 2
    print("PASS sampled numerical budget (not detector-physics certification)")
    return 0

if __name__=="__main__":
    raise SystemExit(main())
