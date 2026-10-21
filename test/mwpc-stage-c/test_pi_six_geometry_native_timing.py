#!/usr/bin/env python3
"""Synthetic tests of the six-geometry true Geant4/Heed/Garfield fold.

No ROOT, Garfield, Heed C++, or native data needed.
"""
from __future__ import annotations
import json
import math
import tempfile
from pathlib import Path

import numpy as np
from pi_six_geometry_native_timing import (
    WIDTHS, DT, NBINS, Count, PITCHES, first_indices,
    n_anodes, signed_bin, response_for_births, aggregate, write_report,
)
from pi_six_geometry_native_timing_html import build


def dummy_map(pitch):
    # n_u=30/35/40, n_v=30; two independent microscopic trials/cell.
    arr=np.full((round(pitch/.1)*30,2,5),np.nan)
    arr[:,:,0]=1  # every test electron collected on nearest anode
    arr[:,:,1]=0
    arr[:,:,2]=25.
    arr[:,:,3]=35.
    arr[:,:,4]=45.
    return arr


def main():
    assert WIDTHS==(307.2,512.,512.,512.,512.,512.)
    assert [n_anodes(0,p) for p in PITCHES]==[102,88,77]
    assert [n_anodes(1,p) for p in PITCHES]==[171,146,128]
    assert NBINS==31 and DT==10
    assert signed_bin(-1.) != signed_bin(+1.)
    assert signed_bin(-10.0) != signed_bin(+10.0)

    # Exact centered local cell and installed p=4 station-specific phase.
    j,cell=first_indices(np.array([-255.,0.,+255.]),4.,128)
    assert np.all(np.abs(cell)<=2.+1.e-12)
    assert list(j)==[0,64,127]
    # A point outside physical gas can land on a nonexistent anode,
    # but finite collection must reject it.
    n=n_anodes(0,4.)
    fine,raw,wire,samples=response_for_births(
        np.array([0.,400.]),np.array([1.,0.]),np.array([.1,.5]),
        4.,0,dummy_map(4.))
    assert raw.tolist()==[True,True]
    assert fine.tolist()==[True,False]
    assert np.all(samples[:,2]==25.)

    # [0,10), [10,20), ... [290,300), and [300,infinity)
    bucket=Count()
    t=np.array([0.,9.9,10.,299.99,300.,500.])
    a=np.full(6,35.)
    bucket.counts(segments=1,primaries=7,collected=6,
                  known=6,unknown=0,edge=1)
    # Keep electron incidence (35 deg) and signed wire-hit azimuth
    # (45 deg) DISTINCT: their histogram definitions are different.
    phi=np.full(6,45.)
    bucket.add_electrons(t,a,phi)
    summary=bucket.done()
    assert summary["hist"][0]==2
    assert summary["hist"][1]==1
    assert summary["hist"][29]==1
    assert summary["hist"][30]==2
    assert sum(summary["hist"])==6
    assert summary["collected"]==6
    assert summary["angle_hist"][3]==6
    assert sum(summary["angle_hist"])==summary["angle_n"]==6
    assert summary["phi_hist"][22]==6  # 45 degrees = bin [40,50)
    assert sum(summary["phi_hist"])==summary["phi_n"]==6
    assert abs(summary["mean_ns"]-np.mean(t))<1.e-9
    assert abs(summary["rms_ns"]-np.std(t))<1.e-9

    # Real first-chunk dictionary field names; no assumptions about Edep
    # to Heed electrons. Exactly the same seed per source/event/segment.
    row={
        "event":"1","hit_index":"3","track_id":"7",
        "station":"0","det_id":"0","pdg":"13",
        "path_length_mm":"6",
        "u_in_mm":"0","u_out_mm":"0",
        "v_in_mm":"-3","v_out_mm":"3",
        "w_in_mm":"0","w_out_mm":"6",
        "theta_u_in_deg":"-12","theta_w_in_deg":"+7",
        "edep_eV":"2000",
    }
    maps={}
    for p in PITCHES:
        for fw in (0,75):
            name=f"p{str(p).replace('.','p')}_"+(
                "anode_only" if fw==0 else "field75p0um")
            maps[name]=(p,fw,dummy_map(p))
    assert len(maps)==6
    # Zero-Heed case must still increment true GAS-HIT denominators
    # for stations and chamber, for all six configurations.
    groups,processed=aggregate([("jpsi40-bore",row)],maps,0.,
                                np.array([1],dtype=np.int32))
    assert processed=={"jpsi40-bore":1}
    key=("jpsi40-bore","p4p0_anode_only",0,-1,-1,"all",-1)
    assert groups[key].segments==1
    assert groups[key].known==0 and groups[key].primaries==0
    # Positive Heed: electron births shared across all six maps.
    groups,processed=aggregate([("jpsi40-bore",row)],maps,8.,
                                np.array([1,2],dtype=np.int32))
    first=("jpsi40-bore","p4p0_anode_only",0,-1,-1,"all",-1)
    nprim=groups[first].primaries
    assert nprim>0
    for name in maps:
        k=("jpsi40-bore",name,0,-1,-1,"all",-1)
        assert groups[k].primaries==nprim
        assert groups[k].collected==nprim
        assert groups[k].known==nprim
        assert sum(groups[k].hist)==nprim
        assert groups[k].hist[2]==nprim
    mu_angle=signed_bin(-12.)
    assert groups[("jpsi40-bore","p4p0_anode_only",0,-1,-1,
                   "u",mu_angle)].known==nprim
    station_anode=groups[("jpsi40-bore","p4p0_anode_only",0,0,38,"all",-1)]
    assert station_anode.known==nprim

    with tempfile.TemporaryDirectory(prefix="six_native_") as tmp:
        path=Path(tmp)
        meta={
            "garfield_luts":{name:{} for name in maps},
            "gas_hit_counts":processed,
            "angle_bins_deg":[-180,-10,0,10,180],
            "input_model":"synthetic",
            "data_scope":"SYNTHETIC ONLY",
            "statistical_errors":"unavailable",
            "finite_anodes_by_pitch":{},
        }
        out=write_report(path,groups,processed,meta)
        target=build(path/"pi_six_geometry_native_timing.json",
                     path/"dashboard.html")
        page=target.read_text(encoding="utf-8")
        assert "Дрейф электронов" in page and "application/json" in page
        assert "SYNTHETIC ONLY" in page
        assert len(out["groups"])>10
    print("PASS: finite 3/3.5/4mm physical lattice, true 10ns overflow,")
    print("PASS: gas-hit N0, common Heed births, signed muon angles,")
    print("PASS: per-anode spectrum conservation and offline dashboard.")


if __name__=="__main__":
    main()
