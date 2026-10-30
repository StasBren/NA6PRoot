#!/usr/bin/env python3
"""SYNTHETIC/OFFLINE ONLY: targeted 20->100 full-avalanche trial design.
Does NOT run any original Geant4, TrackHeed or Garfield physics.
"""
from collections import Counter
from pi_A100_targeted_50cell_campaign import (
    TOP,ADD,TOTAL,CHUNK,N_CHUNKS,GEOMS,EXPORT_FIELDS,spatial_cell,
)
from pi_run_lut500_overnight import selected_j
from pi_run_muon_lut_validation import field_wire_blocked_zero_start
from pi_compare_A100_targeted_vs_A20 import N_REP

def main():
    assert len(GEOMS)==4 and TOP==50 and ADD==80
    assert TOTAL==16000 and N_CHUNKS==40 and CHUNK==400
    assert N_REP==200
    assert EXPORT_FIELDS[:3]==(
        "lut_trials","geometry","original_lut_sha256")
    for g in GEOMS:
        for box in (0,15,215,880 if g.startswith("p3_") else 1199):
            original=[selected_j(g,box,idx) for idx in range(20)]
            new=[selected_j(g,box,idx) for idx in range(20,100)]
            assert len(original)==20 and len(set(original))==20
            assert len(new)==80 and len(set(new))==80
            assert set(original).isdisjoint(new)
            assert all(0<=x<500 for x in original+new)
    # Repeated physical electrons in ONE cell have one shared
    # uncertain 20-sample cell mean: uncertainty scales as n**2.
    count=Counter({11:20,12:1})
    sigma={11:850.,12:400.}
    training=sum(n*n*sigma[b]**2/20 for b,n in count.items())
    wrong_independent=sum(n*sigma[b]**2/20 for b,n in count.items())
    assert training>15*wrong_independent
    # Confirm one real metal start is counted zero only in wire geometry.
    r={"configuration":"p4_FW75","successful":"0","capped":"0",
       "garfield_ne":"1","garfield_ni":"0","ion_births":"0",
       "collected_e":"0","attached":"0","other":"1",
       "ion_up":"0","ion_down":"0","ion_top":"0","ion_bottom":"0",
       "seed_u_mm":"-2.00535685975886",
       "seed_v_mm":"-0.00391288042673045"}
    assert field_wire_blocked_zero_start(r)
    r["configuration"]="p4_noFW"
    assert not field_wire_blocked_zero_start(r)
    print("PASS SYNTHETIC/OFFLINE ONLY: 16,000 extra FULL Garfield starts "
          "from disjoint original 500 D500 indices, 40x400 immutable "
          "chunk plan, n_cell^2 finite A20 mean uncertainty, strict metal "
          "zero-ion classification; NO physics execution.")

if __name__=="__main__":
    main()
