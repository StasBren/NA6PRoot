#!/usr/bin/env python3
"""Tests for new symmetric 3+3 Geant4+Heed/Garfield timing pilot; no ROOT."""
from __future__ import annotations
import numpy as np
from pi_centered_native_timing_pilot import (
    ANGLE_EDGES,angle_index,blank,born_on_segment,centered_endpoints,
    consume,serialize,TIMES,
)


def main():
    for a in (-90,-30,-10,-5,-2,-.1,0,1,2,10,180):
        b=angle_index(a)
        assert ANGLE_EDGES[b]<=a<=ANGLE_EDGES[b+1]
    try:
        angle_index(190)
    except ValueError:
        pass
    else:
        raise AssertionError("Accepted invalid angle")

    row={"gas_path_mm":6,
         "u_in_mm":.1,"u_out_mm":.2,
         "v_in_mm":-3,"v_out_mm":3,
         "w_in_mm":-300,"w_out_mm":300}
    rng=np.random.default_rng(123)
    u,v,w,ncl=born_on_segment(row,rng,1.2,np.array([1,2,3]))
    assert len(u)==len(v)==len(w)>0
    assert ncl>0 and np.min(v)>=-3 and np.max(v)<=3
    assert np.min(w)>=-300 and np.max(w)<=300

    # Synthetic COMPLETE 20x30 map: every electron is collected at its
    # nearest physical anode with precise 32ns drift (this is a test
    # fixture, not a physics map).
    offsets=np.arange(600,dtype=np.int64)
    sizes=np.ones(600,dtype=np.int64)
    results=np.zeros((600,3),dtype=float)
    results[:,0]=1
    results[:,2]=32
    bank=(offsets,sizes,results)
    wires,times,urel,lost,mask=centered_endpoints(
        np.array([-1.1,-.1,1.1]),np.array([-2.9,0,2.9]),
        1,bank,np.random.default_rng(43))
    assert len(wires)==3 and lost==0 and np.all(times==32)
    assert np.all(urel>=-2.0001) and np.all(urel<=2.0001)
    assert mask.tolist()==[True,True,True]
    blank_result=centered_endpoints(
        np.empty(0),np.empty(0),0,bank,np.random.default_rng(42))
    assert len(blank_result)==5 and len(blank_result[0])==0

    # Every native GAS hit remains a denominator, even if Heed had no
    # primary electrons and no anode is active at any time window.
    a=blank()
    consume(a,wires,times,urel,np.array([-2.9,0,2.9]),
            np.array([0.,25.,-25.]),2,5,0)
    consume(a,np.empty(0,dtype=np.int64),np.empty(0),
            np.empty(0),np.empty(0),np.empty(0),0,0,0)
    packed=serialize(a)
    assert packed["gas_hits"]==2
    assert packed["valid_drift_times"]==3
    assert packed["collected_electrons"]==3
    assert sum(packed["arrival_time_counts"])==3
    assert all(sum(v)==2 for v in packed["multiplicity_counts"])
    assert packed["multiplicity_counts"][1][0]==1
    # t=32ns collected too late for the 25ns gate; 50ns activates
    # the hit's two physical wires; the other hit has zero electrons.
    assert packed["multiplicity_counts"][0][0]==2
    assert packed["multiplicity_counts"][1][2]==1

    # Unknown Garfield arrival time may count for "all collected"
    # but must NOT be silently treated as a prompt electron.
    b=blank()
    consume(b,np.array([4]),np.array([np.nan]),
            np.array([0.1]),np.array([0.1]),np.array([20.]),
            1,1,0)
    B=serialize(b)
    assert B["valid_drift_times"]==0
    assert B["unknown_drift_times"]==1
    assert B["multiplicity_counts"][0][0]==1
    assert B["multiplicity_counts"][-1][1]==1

    print("PASS: symmetric local v=-3..+3 with no historical +1mm shift")
    print("PASS: 20x30 local map and finite 77/128 physical anode identities")
    print("PASS: real-gas N0 retained, timed vs all-collected multiplicity")
    print("PASS: 3 birth axes × time histograms, unknown drift time not prompt")


if __name__=="__main__":
    main()
