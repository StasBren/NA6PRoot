#!/usr/bin/env python3
"""Pure synthetic parent-bootstrap errors and exact station-matching tests."""
import math

import numpy as np

from pi_six_geometry_event_bootstrap import (
    N_PARENT,NCFG,NST,NF,T0,T1,AN,AS,AQ,AH0,AH1,PN,PH0,
    ensemble_errors,scalar_metrics,finite_moments,
    source_bootstrap,reference_check,
)


def test_exact_moments_and_errors():
    assert finite_moments(4.,60.,100.)==(15.,5.)
    assert finite_moments(0.,0.,0.)==(None,None)
    p=np.zeros(NF)
    p[0]=10;p[1]=8;p[2]=7;p[3]=70;p[4]=875
    p[T0+0]=3;p[T0+1]=4
    p[AN]=5;p[AS]=100;p[AQ]=2200
    p[AH0+2]=5
    p[PN]=7;p[PH0+21]=7
    s=scalar_metrics(p)
    assert s["n_primaries"]==10 and s["n_anode_collected"]==8
    assert math.isclose(s["anode_collection_fraction"],.8)
    assert s["mean_drift_ns"]==10
    assert math.isclose(s["centered_rms_drift_ns"],5)
    assert s["time_hist"][0]==3 and s["time_hist"][1]==4
    assert math.isclose(sum(s["time_pdf_fraction"]),1)
    assert s["angle_hist"][2]==5 and s["phi_hist"][21]==7
    boot=np.stack([p,p*2,p*3,p*4])
    e=ensemble_errors(boot)
    assert e["mean_drift_ns"]==0
    assert e["centered_rms_drift_ns"]==0
    assert e["anode_collection_fraction"]==0
    assert all(x==0 for x in e["time_pdf_fraction"])
    # Change a single generated parent's physical electron response,
    # retaining cluster-preserving weights, to produce NONZERO SD.
    b2=np.stack([p,p.copy(),p.copy(),p.copy()])
    b2[1,T0]=0
    b2[1,T0+1]=7
    b2[1,3]=120
    b2[1,4]=2450
    err=ensemble_errors(b2)
    assert err["mean_drift_ns"]>0 and err["time_pdf_fraction"][0]>0


def test_parent_resampling():
    cube=np.zeros((N_PARENT,NST,NCFG,NF),dtype=np.float64)
    # Parent 0 emits 20 correlated electrons in one gas chamber.
    cube[0,:, :,0]=20
    cube[0,:, :,1]=20
    cube[0,:, :,2]=20
    cube[0,:, :,3]=200
    cube[0,:, :,4]=2500
    cube[0,:,:,T0+1]=20
    # Parent 1 emits a different cloud, in the same station.
    cube[1,:,:,0]=4
    cube[1,:,:,1]=3
    cube[1,:,:,2]=3
    cube[1,:,:,3]=60
    cube[1,:,:,4]=1500
    cube[1,:,:,T0+2]=3
    names=[f"config_{k}" for k in range(NCFG)]
    out=source_bootstrap("synthetic",cube,names,100,23)
    assert len(out)==NST*NCFG
    g=out[0]
    assert g["n_generated_parents"]==3000
    assert g["statistic"]["n_known_t"]==23
    assert g["statistic"]["time_hist"][1]==20
    assert g["statistic"]["time_hist"][2]==3
    assert g["bootstrap_sd"]["mean_drift_ns"]>0
    assert g["bootstrap_sd"]["time_pdf_fraction"][1]>0


def test_exact_reference_comparison():
    # Frozen source/station/config test reference: all counts and spectral
    # distributions MUST match the independent event-indexed fold.
    totals=np.zeros((NST,NCFG,NF))
    rows=[]
    configs=[f"config_{i}" for i in range(NCFG)]
    for st in range(NST):
        for ci,config in enumerate(configs):
            totals[st,ci,0]=10
            totals[st,ci,1]=8
            totals[st,ci,2]=7
            totals[st,ci,3]=70
            totals[st,ci,4]=875
            totals[st,ci,T0]=7
            totals[st,ci,AN]=5
            totals[st,ci,AS]=100
            totals[st,ci,AQ]=2200
            totals[st,ci,AH0+2]=5
            totals[st,ci,PN]=4
            totals[st,ci,PH0+21]=4
            mean,rms=finite_moments(7,70,875)
            rows.append(dict(source="synthetic",station=st,
                config=config,primaries=10,collected=8,known=7,
                hist=[int(x) for x in totals[st,ci,T0:T1]],
                angle_n=5,
                angle_hist=[int(x) for x in totals[st,ci,AH0:AH1]],
                phi_n=4,phi_hist=[int(x) for x in totals[st,ci,PH0:]],
                mean_ns=mean,rms_ns=rms))
    report={"groups":rows,"meta":{"garfield_luts":{k:{} for k in configs}}}
    assert reference_check("synthetic",totals,report)
    rows[0]["known"]=8
    try:
        reference_check("synthetic",totals,report)
    except AssertionError as exc:
        assert "Full native fold differs" in str(exc)
    else:
        raise AssertionError("Mismatched frozen full-fold check was accepted")


def main():
    test_exact_moments_and_errors()
    test_parent_resampling()
    test_exact_reference_comparison()
    print("PASS: paired-parent-event bootstrap, 3000 generated parents "
          "including zeros, exact full-fold histogram audit")


if __name__=="__main__":
    main()
