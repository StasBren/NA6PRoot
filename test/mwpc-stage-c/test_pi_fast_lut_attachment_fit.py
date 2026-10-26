#!/usr/bin/env python3
"""Synthetic attachment-k fitting regression, no Garfield or full LUT data."""
from __future__ import annotations
import math
from pi_fast_lut_attachment_fit import fit_global_k,mean_pred,one_geom

def item(hit,hazard,angle=45.,time=20.,ref=1,block=(0,0)):
    return {
      "geom":"p3p0_anode_only","split":"fit","cell":(0,0),
      "block":block,"reference_hit":ref,"mc_hit":hit,
      "hazard":hazard,"ref_phi":angle,"mc_phi":angle,
      "ref_time":time,"mc_time":time,"start_v":0.,
      "start_u":0.}

def main():
    # Analytically k=1 for eta integrated hazard log(2) and
    # half of the microscopic training electrons surviving.
    h=math.log(2)
    starts=[item(1,h,ref=int(i%2==0),
                 block=(i//10,0)) for i in range(100)]
    k=fit_global_k(starts)
    assert abs(k-1.)<1e-10,k
    assert abs(mean_pred(starts,k)-.5)<1e-10
    result=one_geom(starts,k)
    assert abs(result["predicted_survival_fraction"]-.5)<1e-10
    assert result["phi_tv_36bins_weighted_survivors"]==0.
    assert abs(result["mean_drift_model_weighted_ns"]-20.)<1e-10

    # A H=0 path can never undergo attachment under k*eta.
    # Calibrator must reject a reference efficiency below H=0 floor.
    impossible=[item(1,0.,ref=0) for _ in range(30)]
    impossible += [item(1,1.,ref=0) for _ in range(30)]
    try:
        fit_global_k(impossible)
    except ValueError as e:
        assert "hazard-zero floor" in str(e)
    else:
        raise AssertionError("Calibration silently fitted unreachable target")

    # Un-killed MC cannot be less efficient than reference if we want
    # to use attachment alone to reproduce the missing electron losses.
    impossible=[item(0,1.,ref=1) for _ in range(10)]
    try:
        fit_global_k(impossible)
    except ValueError as e:
        assert "fewer collected electrons" in str(e)
    else:
        raise AssertionError("Calibration increased survival above MC max")
    print("PASS attachment-k synthetic: known k=1 recovered, "
          "weighted phi/time, impossible-survival fit rejected")

if __name__=="__main__":
    main()
