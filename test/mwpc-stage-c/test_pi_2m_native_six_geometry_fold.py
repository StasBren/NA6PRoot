#!/usr/bin/env python3
"""No Geant4/Garfield jobs: validate source+chunk seed and exact moment merging."""
from __future__ import annotations
import math
import numpy as np
from pi_six_geometry_native_timing import (
    Count,aggregate,validate_conservation,stable_seed,PITCHES,tag,
)
from pi_2m_native_six_geometry_fold import (
    as_row,key,restore,snapshot,
)
from test_pi_six_geometry_native_timing import dummy_map


def main():
    # Controlled independent source+chunk+event seed namespace.
    start=("native-six-pitch-shared-Heed-v1","jpsi40-bore")
    seed_1=stable_seed(*start,"chunk0001",5,3,7)
    seed_2=stable_seed(*start,"chunk0002",5,3,7)
    seed_phi=stable_seed("native-six-pitch-shared-Heed-v1",
                         "phi40-bore","chunk0001",5,3,7)
    assert seed_1!=seed_2 and seed_1!=seed_phi
    first=Count();second=Count()
    t1=np.array([1.,21.,41.])
    t2=np.array([11.,31.])
    first.counts(segments=1,primaries=4,collected=3,known=3,edge=0)
    first.add_electrons(t1,np.full(3,30.),np.full(3,-25.))
    second.counts(segments=1,primaries=3,collected=2,known=2,edge=0)
    second.add_electrons(t2,np.full(2,40.),np.full(2,15.))
    single=Count()
    restore(single,snapshot(first))
    restore(single,snapshot(second))
    out=single.done()
    truth=np.concatenate((t1,t2))
    assert out["segments"]==2 and out["primaries"]==7
    assert out["collected"]==out["known"]==5
    assert sum(out["hist"])==5
    assert math.isclose(out["mean_ns"],float(np.mean(truth)),abs_tol=1.e-12)
    assert math.isclose(out["rms_ns"],float(np.std(truth)),abs_tol=1.e-12)
    assert out["angle_n"]==5 and out["phi_n"]==5
    assert out["angle_hist"][3]==3 and out["angle_hist"][4]==2
    # Arithmetic mean of individual chunk RMS is NOT true merged RMS.
    assert abs(out["rms_ns"]-(np.std(t1)+np.std(t2))/2)>1.e-2
    # Original first-chunk pipeline still supports source overview by
    # default; production opts out to avoid doubling memory.
    row={
        "event":"1","hit_index":"3","track_id":"7",
        "station":"0","det_id":"0","pdg":"13",
        "path_length_mm":"6","u_in_mm":"0","u_out_mm":"0",
        "v_in_mm":"-3","v_out_mm":"3",
        "w_in_mm":"0","w_out_mm":"6",
        "theta_u_in_deg":"-12","theta_w_in_deg":"+7",
        "edep_eV":"2000",
    }
    maps={tag(p,d):(p,d,dummy_map(p))
          for p in PITCHES for d in (0.,75.)}
    rows=[("jpsi40-bore",row),
          ("jpsi40-bore",{**row,"station":"1","det_id":"1"})]
    groups,processed=aggregate(rows,maps,4.,np.array([1,2],dtype=int),
           source_overview=False,station_filter=0,chunk_key="chunk0002")
    assert processed=={"jpsi40-bore":1}
    assert all(k[0]=="jpsi40-bore" and k[2]==0 for k in groups)
    assert len(groups)>12
    validate_conservation(groups)
    key0=next(k for k in groups if k[3]==-1 and k[5]=="all")
    assert key(as_row(key0,groups[key0]))==key0
    print("PASS: unique native source/chunk Heed seeds, no all-source "
          "RAM duplication, selected-station fold conserves actual "
          "anode electrons, exact two-chunk time mean and centered RMS")


if __name__=="__main__":
    main()
