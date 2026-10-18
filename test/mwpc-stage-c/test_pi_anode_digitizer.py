#!/usr/bin/env python3
"""Synthetic deterministic tests: anode-only stochastic gain and digitization.

No native Geant4, new Garfield, cathode strips or ROOT imports required.
"""
from dataclasses import replace
import csv
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

from pi_anode_digitizer import (
    AnodeConfig, digitize_anode, digitize_hit, draw_gains,
    empirical_gains, threshold_features, ELEMENTARY_CHARGE_FC
)


def main():
    cfg=AnodeConfig(noise_rms_mV=0.,threshold_mV=0.1,
                    mean_gain=10000.,zero_gain_probability=0.,
                    observation_ns=400.,sampling_ns=10.)
    a=digitize_anode([0.,30.,50.],np.random.default_rng(42),cfg)
    b=digitize_anode([0.,30.,50.],np.random.default_rng(42),cfg)
    assert a==b, "fixed seed should give identical waveforms"
    assert a["nPrimaryCollected"]==3
    assert a["nTimed"]==3 and a["nUnknownOrOutOfWindow"]==0
    assert a["nZeroAvalanches"]==0
    assert abs(a["avalancheCharge_fC"] -
               a["effectiveGainSum_e"]*ELEMENTARY_CHARGE_FC)<1.e-4
    assert len(a["samples_mV"])==41
    assert a["signalOnlyPeak_mV"]>=0.
    assert all(np.isfinite(a["samples_mV"]))
    assert (threshold_features(a["samples_mV"],cfg,threshold=1000.)["fired"]
            is False)
    assert (threshold_features(a["samples_mV"],cfg,threshold=0.00001)["fired"]
            is True)
    dead=replace(cfg,zero_gain_probability=0.999999)
    e=digitize_anode([10.,30.],np.random.default_rng(1),dead)
    assert e["nZeroAvalanches"]==2
    assert not e["fired"]
    assert e["avalancheCharge_fC"]==0
    missing=digitize_anode([float("nan"),-2.,5.,600.],
                           np.random.default_rng(1),cfg)
    assert missing["nTimed"]==1
    assert missing["nUnknownOrOutOfWindow"]==3

    hits=digitize_hit([(7,15.),(7,35.),(8,60.)],
                      np.random.default_rng(987),cfg)
    assert set(hits)=={"7","8"}
    assert hits["7"]["nPrimaryCollected"]==2
    assert hits["8"]["nPrimaryCollected"]==1

    with TemporaryDirectory(prefix="pi_anode_test_") as temp:
        p=Path(temp)/"phaseA_gain.csv"
        with p.open("w",newline="",encoding="utf-8") as stream:
            writer=csv.writer(stream)
            writer.writerow(["event","collected_on_wire"])
            for k in range(40):
                writer.writerow([k,0 if k%4==0 else 12000+k])
        gains=empirical_gains(p)
        assert len(gains)==40 and np.count_nonzero(gains==0)==10
        sampled=draw_gains(np.random.default_rng(5),100,cfg,gains)
        assert all(v in set(gains) for v in sampled)
        empirical=digitize_anode([5.]*30,np.random.default_rng(5),
                                cfg,gains)
        assert 0<=empirical["nZeroAvalanches"]<=30
    try:
        digitize_anode([0.],np.random.default_rng(1),
                       replace(cfg,fast_charge_fraction=2.))
        raise AssertionError("Invalid physical parameters accepted")
    except ValueError:
        pass
    print("PASS anode-only: deterministic waveforms, gain charge, timing,")
    print("PASS 100% zero gain => no electronics fired, threshold scan,")
    print("PASS 40-row true zero-inclusive empirical gain library, no cathodes")


if __name__=="__main__":
    main()
