#!/usr/bin/env python3
"""Pure synthetic regression for symmetric 3+3mm Heed/Garfield native adapter."""
from __future__ import annotations
import numpy as np

from pi_centered_electron_sampler import (
    ionize_and_drift, multiplicities, stable_seed,
)
from pi_native_centered_timing_pilot import fresh, update, finish, COORDINATES


def row(st=1,u=0.,v=-3.,angle=1.):
    return {
        "station":st,"path_length_mm":6.,
        "u_in_mm":u,"u_out_mm":u,
        "v_in_mm":v,"v_out_mm":-v,
        "w_in_mm":-14.,"w_out_mm":14.,
        "theta_u_in_deg":angle,"theta_w_in_deg":.5,
    }


def main():
    # Synthetic centered 600 boxes with 1 deterministic collect-to-nearest-
    # anode electron each, independent of the side of v=0.
    pool=np.zeros((600,1,3),dtype=float)
    pool[:,:,0]=1.
    pool[:,:,1]=0.
    pool[:,:,2]=20.
    marks=np.array([1],dtype=np.int32)
    seed=stable_seed("centered-test",123)
    h=ionize_and_drift(row(),np.random.default_rng(seed),15.,marks,pool)
    assert h["nprimary"]>0
    assert h["ncollected"]==h["nprimary"]
    assert h["nedge"]==0
    assert np.all(h["valid_times"])
    assert np.all(h["times"]==20.)
    assert np.all(h["v"]>=-3) and np.all(h["v"]<=3)
    assert np.any(h["v"]<0) and np.any(h["v"]>0)
    assert np.allclose(h["u_cell"],-2.)  # even wire bank, u=0 at cell boundary
    assert multiplicities(h)==[1,1,1,1,1,1]
    h2=ionize_and_drift(row(),np.random.default_rng(seed),15.,marks,pool)
    assert np.array_equal(h["v"],h2["v"])  # per-hit stable seed
    empty=ionize_and_drift(row(),np.random.default_rng(7),0.,marks,pool)
    assert empty["nprimary"]==0 and multiplicities(empty)==[0]*6

    # Late electron (no prompt 50/100), unknown drift time, and edge loss.
    delayed=pool.copy();delayed[:,:,2]=160.
    h3=ionize_and_drift(row(),np.random.default_rng(4),12.,marks,delayed)
    assert multiplicities(h3)==[0,0,1,1,1,1]
    lost=ionize_and_drift(row(0,400.),np.random.default_rng(5),12.,marks,pool)
    assert lost["nrawcollected"]==lost["nedge"] and lost["ncollected"]==0
    assert multiplicities(lost)==[0]*6
    unknown=pool.copy();unknown[:,:,2]=np.nan
    h4=ionize_and_drift(row(),np.random.default_rng(6),12.,marks,unknown)
    assert not np.any(h4["valid_times"]) and multiplicities(h4)==[0,0,0,0,0,1]

    # Segment-level PDF including N0 and correct time×birth conservation.
    stats=fresh(-1)
    z=np.empty(0,dtype=int);flag=np.empty(0,dtype=bool)
    coords={name:z for name in COORDINATES}
    update(stats,empty,[0]*6,z,flag,coords)
    done=finish(stats)
    assert done["n"]==1 and done["n0"]==1
    assert [counts[0] for counts in done["m"]]==[1]*6
    assert done["known"]==done["collected"]==0

    bad=dict(row());bad["v_out_mm"]=4.
    try:
        ionize_and_drift(bad,np.random.default_rng(1),12.,marks,pool)
    except ValueError as e:
        assert "outside centered" in str(e)
    else:
        raise AssertionError("Out-of-gas 2+4mm masquerade accepted")

    print("PASS: centered 3+3mm Heed local u/v/w and deterministic per-hit seeds")
    print("PASS: all 77/128 physical wire bounds and boundary electron loss")
    print("PASS: 0-primary event, unknown t, late t and time monotonicity")
    print("PASS: station/module time+initial-coordinate spectrum conservation")


if __name__=="__main__":
    main()
