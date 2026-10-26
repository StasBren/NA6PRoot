#!/usr/bin/env python3
"""Small synthetic regression for eventwise Garfield LUT500 asymmetry logic."""
from __future__ import annotations

import math
import numpy as np
from pi_event_asymmetry_500_ms0 import (
    up_down,asymmetry,Moments,summarize_parent_counts,A_EDGES
)


def run():
    # Accepted anode hits with positive/negative phi and exact boundaries.
    # Rejected hits must never affect up/down or event A.
    sample=np.full((9,5),np.nan)
    sample[:,4]=[90.,-90.,30.,-170.,0.,180.,-180.,-40.,80.]
    mask=np.array([True,True,True,True,True,True,True,False,True])
    u,d,b=up_down(sample,mask)
    assert (u,d,b)==(3,2,3), (u,d,b)
    assert math.isclose(asymmetry(u,d),.2)
    assert asymmetry(0,0) is None
    assert asymmetry(0,5)==-1.
    assert asymmetry(5,0)==1.

    stat=Moments()
    stat.add(10,0)
    stat.add(0,10)
    stat.add(0,0)
    got=stat.result()
    assert got["hits"]==3 and got["valid"]==2 and got["empty"]==1
    assert got["mean"]==0. and got["rms"]==1.
    assert got["mean_collected"]==10.
    assert got["electron_weighted_asymmetry"]==0.

    # Event index zero appears in TWO different Geant4 chunks.
    # It must remain two independent parent events, never merged.
    ch1=np.zeros((4,2),dtype=np.int64)
    ch2=np.zeros((4,2),dtype=np.int64)
    ch1[0]=[8,0]
    ch2[0]=[0,8]
    r=summarize_parent_counts([
        {"pair":ch1,"has_crossing":np.array([1,0,0,0],dtype=bool)},
        {"pair":ch2,"has_crossing":np.array([1,0,0,0],dtype=bool)},
    ])
    assert r["parents_generated"]==8
    assert r["parents_with_MS0_crossing"]==2
    assert r["parents_with_classified_electrons"]==2
    assert r["n_with_undefined_asymmetry"]==6
    assert r["mean"]==0. and r["rms"]==1.
    assert sum(r["histogram_counts"])==2 and len(r["histogram_counts"])==41
    assert r["histogram_edges"]==A_EDGES.tolist()
    assert r["histogram_counts"][0]==1 and r["histogram_counts"][-1]==1
    print("PASS: exclusive up/down/boundary classification, empty event "
          "policy, correct RMS and separate source-chunk-local parents")


if __name__=="__main__":
    run()
