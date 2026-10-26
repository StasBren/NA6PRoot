#!/usr/bin/env python3
"""No-Garfield smoke checks for native/on versus offline attachment audit."""
from __future__ import annotations
import math
from pi_fast_lut_native_attachment_audit import (
    bin_checks,kfit_target,residuals
)

def entry(i,p):
    return {
        "id":("geometry",0,0,i),
        "block":(i//100,0),
        "reference_hit":1 if i%2==0 else 0,
        "mc_hit":1,"hazard":math.log(2),
    }

def main():
    entries=[entry(i,1) for i in range(3600)]
    native={r["id"]:{"anode_collected":(1 if j%2==0 else 0)}
        for j,r in enumerate(entries)}
    k=kfit_target(entries,native)
    assert abs(k-1.)<1e-9,k
    bins=bin_checks(entries,native,k)
    assert len(bins)==10 and all(x["n"]==360 for x in bins)
    assert all(abs(x["estimated_k1_survival"]-.5)<1e-9 for x in bins)
    assert all(abs(x["native_mc_attachment_on_anode"]-.5)<1e-9 for x in bins)
    blocks=residuals(entries,k)
    assert blocks["rms"]<1e-9
    assert blocks["max_abs"]<1e-9
    print("PASS native-attachment audit synthetic: effective k, "
          "eta deciles, 6x6-block residual formulas")

if __name__=="__main__":main()
