#!/usr/bin/env python3
"""Synthetic fast Heed expectation tests: no ROOT, Geant4, or Garfield."""
from __future__ import annotations

import math
import numpy as np

from pi_lut_20_vs_500_expected import segment_fractions
from pi_six_geometry_native_timing import first_indices


def run():
    n=102  # Actual p=3mm MS0 finite anode count
    # Normal incidence: every 0.2mm v layer gets 1/30 of the 6mm chord.
    boxes,wires,f=segment_fractions(-.75,-3.,-.75,3.,n)
    assert len(f)==30 and np.allclose(f,np.full(30,1./30.),atol=1e-12)
    assert len(np.unique(wires))==1
    assert np.array_equal(boxes%30,np.arange(30))

    # Across the central Voronoi boundary: exactly half on each anode.
    b,w,f=segment_fractions(-.05,0.,.05,0.,n)
    assert len(f)==2 and np.allclose(f,[.5,.5])
    assert np.array_equal(w,[n//2,n//2+1])
    assert np.array_equal(b%30,[15,15])
    assert np.array_equal(b//30,[29,0])

    # Reversal symmetry is a consequence of integration, not RNG.
    for u0,v0,u1,v1 in (
        (-1.4,-2.9,1.4,2.9),(-.05,-1.1,.31,1.1),
        (1.37,2.81,-1.31,-2.5),(.03,-2.99,.13,2.96)):
        x,y,f=segment_fractions(u0,v0,u1,v1,n)
        xr,yr,fr=segment_fractions(u1,v1,u0,v0,n)
        z=np.zeros((900,n+2),dtype=float)
        zr=np.zeros_like(z)
        np.add.at(z,(x,y),f)
        np.add.at(zr,(xr,yr),fr)
        assert np.allclose(z,zr,atol=1e-12,rtol=0.),(u0,v0,u1,v1)
        assert math.isclose(float(z.sum()),1.,abs_tol=1e-12)

    # Agreement with uniform sampling along the exact Geant4 chord.
    # This is test-only numerical quadrature, not the runtime method.
    u0,v0,u1,v1=-1.37,-2.43,.84,2.39
    b,w,f=segment_fractions(u0,v0,u1,v1,n)
    exact=np.zeros((900,n+2))
    np.add.at(exact,(b,w),f)
    t=(np.arange(200000,dtype=float)+.5)/200000
    u=u0+(u1-u0)*t
    v=v0+(v1-v0)*t
    nearest,cell=first_indices(u,3.,n)
    iu=np.clip(np.floor((cell+1.5)/.1).astype(int),0,29)
    iv=np.clip(np.floor((v+3.)/.2).astype(int),0,29)
    numeric=np.zeros_like(exact)
    np.add.at(numeric,(iu*30+iv,nearest+1),1/len(t))
    assert np.max(np.abs(exact-numeric))<1e-5, np.max(np.abs(exact-numeric))
    print("PASS: 30 v layers, finite-wire split, reversal, "
          "spatial occupancy conservation, direct segment quadrature")


if __name__=="__main__":
    run()
