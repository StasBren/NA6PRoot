#!/usr/bin/env python3
"""Synthetic test ONLY: matched Poisson-Heed clusters, finite binomial
cell fractions and identical six-Garfield LUT input electrons.

No native Geant4/Heed/Garfield executables or experiment data required.
"""
from __future__ import annotations

import math
import numpy as np
from pi_geometric_vs_lut_benchmark import (
    shared_heed_births,geometric_wire_indices,
    two_cell_probabilities,compute,THRESHOLDS,
)
from pi_six_geometry_native_timing import (
    birth_sample,n_anodes,first_indices,stable_seed,
)
from test_pi_six_geometry_native_timing import dummy_map
from pi_pitch_fieldwire_pilot import PITCHES,tag


def main():
    row={
        "event":"45","hit_index":"3","track_id":"7",
        "station":"0","det_id":"0","pdg":"13",
        "path_length_mm":"7.1",
        "u_in_mm":"-1.5","u_out_mm":"1.5",
        "v_in_mm":"-3","v_out_mm":"3",
        "w_in_mm":"0","w_out_mm":"0",
        "theta_u_in_deg":"21","theta_w_in_deg":"0",
    }
    marks=np.array([1,1,2,3,5,8],dtype=np.int64)
    key=("native-six-pitch-shared-Heed-v1","jpsi40-bore",
         "chunk0001",row["event"],row["hit_index"],row["track_id"])
    seed=stable_seed(*key)
    pos,q,tc,n_e,length,endpoints=shared_heed_births(row,5.4,marks,seed)
    expected_pos,expected_q=birth_sample(row,5.4,marks,seed)
    assert np.array_equal(pos,expected_pos)
    assert np.array_equal(q,expected_q)
    assert len(pos)==int(n_e.sum())
    assert len(tc)==len(n_e)
    assert np.all(n_e>=1)
    assert np.array_equal(
        np.repeat(endpoints[0,0]+(endpoints[1,0]-endpoints[0,0])*tc,n_e),
        pos[:,0])
    # MS0 p=3.0 has an even number of anodes, central boundary u=0:
    # our symmetric crossing MUST have exactly equal geometric path shares.
    idx,segments=geometric_wire_indices(-1.5,1.5,tc,3.,0)
    assert n_anodes(0,3.)==102
    assert len(segments)==2 and np.allclose(segments,[.5,.5])
    assert set(idx.tolist()).issubset({50,51})
    conditional,unconditional,pminor=two_cell_probabilities(
        segments,len(n_e),5.4*length)
    assert math.isclose(conditional,
        1.-2.*(.5**len(n_e)),abs_tol=1.e-12)
    assert math.isclose(unconditional,
        (1-math.exp(-5.4*length*.5))**2,abs_tol=1.e-12)
    assert pminor==.5
    # Poisson/binomial model is about whole CLUSTERS, not independent
    # primary electrons. All m_i electrons have the same assigned wire.
    counts=np.bincount(idx,minlength=102)
    assert sum(counts)==len(n_e)
    assert sum(np.repeat(idx,n_e)==50)+sum(np.repeat(idx,n_e)==51)==len(pos)
    # Same born electron population and same q for each of SIX LUTs.
    pools={tag(p,d):(p,d,dummy_map(p)) for p in PITCHES for d in (0.,75.)}
    rows=[row,{**row,"event":"46","u_in_mm":"-.4",
               "u_out_mm":".4"}]
    comparison,clock,common,two=compute(
        rows,"Jpsi",1,0,5.4,marks,pools)
    assert len(comparison)==6
    for cfg,stat in comparison.items():
        assert stat["gas_hits"]==2
        assert stat["primary_electrons"]>0
        assert stat["lut_collected_electrons"]==stat["primary_electrons"]
        assert stat["lut_lost_electrons"]==0
        assert stat["lut_migrated_collected_electrons"]==0
        assert clock[cfg]["geometry_ns"]>0
        assert clock[cfg]["lut_ns"]>0
        assert two[cfg]["n"]<=2
        for threshold in THRESHOLDS:
            assert 0<=stat[f"geo_multiwire_atleast_{threshold}e"]<=2
            assert 0<=stat[f"lut_multiwire_atleast_{threshold}e"]<=2
    assert common>0
    print("PASS: exact unchanged original Heed electron seed/positions/q; "
          "cluster multiplicity survives binomial two-wire splitting; "
          "analytic two-cell Poisson+Binomial sharing; matched six-LUT "
          "geometry vs response clock scopes and finite-anode conservation.")


if __name__=="__main__":
    main()
