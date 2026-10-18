#!/usr/bin/env python3
"""Anode-only experimental fast digitizer. NO CATHODE STRIPS.

A real per-electron Garfield transport endpoint (wire, drift time) feeds:
- independent stochastic gas gains: empirical UNCONDITIONED Phase-A gain
  histogram if available, otherwise explicitly UNCALIBRATED Polya+zero model;
- approximate per-avalanche fast electron and slow-ion INDUCED charge streams;
- nominal charge-sensitive amplifier CR-RC impulse response (user parameters);
- fixed once-per-event Gaussian time-sample noise, threshold and time-over-
  threshold digitization for each ANODE that collected primary electrons.

The waveform/current split is a PARAMETRIC SURROGATE, NOT a simulation of
the actual anode weighting field and amplifier. True Shockley-Ramo templates
and experimental calibration are REQUIRED for predictive efficiency.
Saved waveform samples can be rethresholded instantly in an offline viewer.

Does not simulate noise-only hits on unvisited anodes or spatial crosstalk.
No cathode/stereo/strip signals or reconstruction are included.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import csv
import math
from pathlib import Path

import numpy as np

ELEMENTARY_CHARGE_FC = 1.602176634e-4


@dataclass(frozen=True)
class AnodeConfig:
    mean_gain: float = 10000.0
    polya_shape: float = 1.5
    zero_gain_probability: float = 0.05
    fast_charge_fraction: float = 0.2
    fast_decay_ns: float = 5.0
    ion_decay_ns: float = 1200.0
    peaking_ns: float = 80.0
    preamp_mV_per_fC: float = 0.35
    noise_rms_mV: float = 0.35
    threshold_mV: float = 5.0
    sampling_ns: float = 10.0
    observation_ns: float = 3000.0
    deadtime_ns: float = 0.0

    def check(self):
        if not (all(math.isfinite(v) for v in asdict(self).values()) and
                self.mean_gain > 0 and self.polya_shape > 0 and
                0 <= self.zero_gain_probability < 1 and
                0 <= self.fast_charge_fraction <= 1 and
                self.fast_decay_ns > 0 and self.ion_decay_ns > 0 and
                self.peaking_ns > 0 and self.preamp_mV_per_fC > 0 and
                self.noise_rms_mV >= 0 and self.threshold_mV > 0 and
                self.sampling_ns > 0 and self.observation_ns >= 2*self.sampling_ns and
                self.deadtime_ns >= 0):
            raise ValueError("Invalid anode gain/shaping/noise/time parameters")
        n=int(math.floor(self.observation_ns/self.sampling_ns))+1
        if n > 10000:
            raise ValueError("Anode samples/window too large")
        return n


def empirical_gains(path: Path) -> np.ndarray:
    """Import true per-seed zero-inclusive collected-on-wire gains from A.

    Do not use min-avalanche-ion selected Stage-B example libraries here:
    their gain distribution is conditioned on nonzero amplification.
    """
    gains=[]
    with path.open(newline="",encoding="utf-8") as stream:
        rd=csv.DictReader(stream)
        if "collected_on_wire" not in (rd.fieldnames or []):
            raise ValueError("Need UNCONDITIONED Stage-A gain CSV with collected_on_wire")
        for row in rd:
            g=float(row["collected_on_wire"])
            if not math.isfinite(g) or g < 0:
                raise ValueError("Invalid empirical avalanche gain")
            gains.append(g)
    if len(gains)<20:
        raise ValueError("Empirical gain bank needs at least 20 seeds")
    return np.asarray(gains,dtype=np.float64)


def draw_gains(rng, n, cfg, empirical=None):
    if empirical is not None:
        return rng.choice(empirical,size=n,replace=True).astype(np.float64)
    # Do not resample zero-gain avalanches; that would bias the efficiency.
    gains=rng.gamma(cfg.polya_shape,cfg.mean_gain/cfg.polya_shape,size=n)
    gains[rng.random(n)<cfg.zero_gain_probability]=0.
    return gains


def threshold_features(wave_mV, cfg, threshold=None):
    threshold=cfg.threshold_mV if threshold is None else float(threshold)
    if not (math.isfinite(threshold) and threshold > 0):
        raise ValueError("Invalid discriminator threshold")
    above=np.flatnonzero(np.asarray(wave_mV) >= threshold)
    return {
        "fired":bool(above.size),
        "time_ns":round(float(above[0])*cfg.sampling_ns,3) if above.size else None,
        "tot_ns":round(float(above.size)*cfg.sampling_ns,3),
        "peak_mV":round(float(np.max(wave_mV)),4),
        "threshold_mV":threshold,
    }


def digitize_anode(endpoints_ns, rng, cfg=AnodeConfig(), empirical=None):
    """Digitize ONE wire, with all its collected primary-electron arrival times.

    Returns full waveform and physically labelled intermediate charge.
    Invalid Garfield drift times remain unmodelled and are EXPLICITLY
    counted, never silently recast as prompt arrivals.
    """
    nbins=cfg.check()
    times=np.asarray(endpoints_ns,dtype=np.float64)
    if times.ndim!=1:
        raise ValueError("Expected flat electron arrival list")
    finite=np.isfinite(times)&(times>=0.)&(times<=cfg.observation_ns)
    discarded=int(times.size-np.count_nonzero(finite))
    arrival=times[finite]
    gain=draw_gains(rng,len(arrival),cfg,empirical)
    charge_fC=float(np.sum(gain)*ELEMENTARY_CHARGE_FC)
    zeros=int(np.count_nonzero(gain==0))
    # Binary charge contributions on discrete timing bins, in fC.
    rapid=np.zeros(nbins,dtype=np.float64)
    ions=np.zeros(nbins,dtype=np.float64)
    bins=np.minimum(np.floor(arrival/cfg.sampling_ns).astype(int),nbins-1)
    qq=gain*ELEMENTARY_CHARGE_FC
    np.add.at(rapid,bins,qq*cfg.fast_charge_fraction)
    np.add.at(ions,bins,qq*(1.-cfg.fast_charge_fraction))

    # Fast electron charge fraction and slow positive-ion charge stream.
    # Each packet has exponential shape, integrated over finite bins.
    # This is NOT yet the true Ramo induced current of the MWPC anode.
    induced=np.zeros(nbins,dtype=np.float64)
    reservoir_fast=reservoir_ions=0.
    af=math.exp(-cfg.sampling_ns/cfg.fast_decay_ns)
    ai=math.exp(-cfg.sampling_ns/cfg.ion_decay_ns)
    for k in range(nbins):
        reservoir_fast=af*reservoir_fast+rapid[k]
        reservoir_ions=ai*reservoir_ions+ions[k]
        induced[k]=(1-af)*reservoir_fast+(1-ai)*reservoir_ions

    # Nominal CR-RC response with a unit-peak impulse for an incoming
    # charge packet, multiplied by the nominal voltage conversion.
    t=np.arange(nbins,dtype=float)*cfg.sampling_ns
    h=np.e*(t/cfg.peaking_ns)*np.exp(-t/cfg.peaking_ns)
    h[0]=0.
    signal=np.convolve(induced,h,mode="full")[:nbins]
    signal*=cfg.preamp_mV_per_fC
    noisy=signal+rng.normal(0.,cfg.noise_rms_mV,size=nbins)
    discr=threshold_features(noisy,cfg)
    return {
        "nPrimaryCollected":int(times.size),
        "nTimed":int(arrival.size),
        "nUnknownOrOutOfWindow":discarded,
        "nZeroAvalanches":zeros,
        "avalancheCharge_fC":round(charge_fC,5),
        "effectiveGainSum_e":round(float(np.sum(gain)),2),
        "firstElectron_ns":round(float(np.min(arrival)),3) if arrival.size else None,
        "baselineNoise_rms_mV":cfg.noise_rms_mV,
        "samples_mV":[round(float(v),3) for v in noisy],
        "signalOnlyPeak_mV":round(float(np.max(signal)),3),
        **discr,
    }


def digitize_hit(endpoints, rng, cfg=AnodeConfig(), empirical=None):
    """endpoints: per collected primary electron (anodeIndex, driftTime_ns).

    Deliberately digitizes only collected-electron anodes; noise-only
    candidates among the other finite anodes are NOT evaluated.
    """
    by_wire={}
    for wire,t in endpoints:
        by_wire.setdefault(int(wire),[]).append(float(t))
    return {
        str(wire):digitize_anode(tlist,rng,cfg,empirical)
        for wire,tlist in sorted(by_wire.items())
    }
