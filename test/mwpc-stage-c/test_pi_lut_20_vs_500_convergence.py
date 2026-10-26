#!/usr/bin/env python3
"""Synthetic no-Geant4 test of LUT finite-wire weighting and correlated SE."""
from __future__ import annotations
import math
import numpy as np
from pi_lut_20_vs_500_convergence import microscopic_pdf, NCELL

def run():
    occupancy=np.zeros((NCELL,4),dtype=np.int64) # nwire=2; nearest -1,0,1,2
    occupancy[0,1]=10
    occupancy[0,2]=20
    pool=np.full((NCELL,3,5),np.nan,dtype=float)
    pool[:,:,0]=0.
    pool[:,:,1]=0.
    pool[0,:,0]=1.
    pool[0,:,4]=[45.,45.,-45.]
    d=microscopic_pdf(occupancy,pool)
    p=np.array(d["phi_probabilities"])
    se=np.array(d["garfield_mc_standard_errors"])
    assert math.isclose(d["collection_fraction"],1)
    assert math.isclose(d["expected_collected_heed_electrons"],30)
    # Only 3 statistically independent Garfield drifts, not 90 correlated
    # "Heed samples" after the two physical nearest wires reuse the same LUT.
    plus=int(np.searchsorted(np.arange(-180,181,10),45.,side="right")-1)
    minus=int(np.searchsorted(np.arange(-180,181,10),-45.,side="right")-1)
    assert math.isclose(p[plus],2/3)
    assert math.isclose(p[minus],1/3)
    assert math.isclose(se[plus],1/3)
    assert math.isclose(se[minus],1/3)
    assert math.isclose(p.sum(),1)

    # A primary beyond the first finite anode must be rejected if Garfield
    # reports zero wire migration and accepted if it migrates +1 wire.
    edge_occ=np.zeros((NCELL,4),dtype=np.int64)
    edge_occ[0,0]=12 # physical nearest wire=-1
    edge_pool=np.full((NCELL,3,5),np.nan,dtype=float)
    edge_pool[:,:,0]=0.
    edge_pool[:,:,1]=0.
    edge_pool[0,:,0]=1.
    edge_pool[0,:,1]=[0.,1.,1.]
    edge_pool[0,:,4]=[45.,45.,-45.]
    edge=microscopic_pdf(edge_occ,edge_pool)
    assert math.isclose(edge["collection_fraction"],2/3)
    assert math.isclose(edge["expected_collected_heed_electrons"],8)
    assert math.isclose(sum(edge["phi_probabilities"]),1)
    print("PASS: finite anode edge, shared microscopic trial covariance, "
          "probability normalization and Garfield-only uncertainty")

if __name__=="__main__":
    run()
