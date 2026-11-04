#!/usr/bin/env python3
"""Read-only higher-cap audit of SHA-sealed P3 full Garfield ion-BIRTH clouds.

No new Garfield simulation, ion drift, weighting field or cathode signal.
This is a paired within-existing-avalanche STORAGE-compression experiment.
It complements, never overwrites, sealed 32/128/256/512/1024 analysis.
"""
from __future__ import annotations

import argparse
import math
import random
from pathlib import Path

import numpy as np

import pi_p3_extended_full_cloud_campaign_v1 as campaign
import pi_p3_extended_cloud_readout_v1 as readout

core = campaign.core
TAG = 'NA60-P3-EXTENDED-HIGHCAP-20261104-V1'
DEFAULT_CAPS = (2048, 4096)
ALLOWED_CAPS = (2048, 4096, 8192)
METRICS = tuple(key for key, _ in readout.METRICS)


def must(condition, message):
    if not condition:
        raise ValueError(message)


def score_cloud(points, caps, reps, geom, box, sample):
    """Repeated reservoirs of ONE exact birth cloud; no extra physics trials."""
    must(len(points) > 0, 'Do not score a zero-ion birth cloud')
    full = core.cloud_stats(points)
    baseline = core.joint_hist(full, float(np.min(full['t'])),
                               float(np.max(full['t'])))
    result = {}
    for cap in caps:
        k = min(cap, len(points))
        rng = random.Random(core.stable(TAG, geom, box, sample, cap))
        errors = []
        if k == len(points):
            errors = [core.compression_errors(points, points, full, baseline)]
        else:
            for _ in range(reps):
                indices = rng.sample(range(len(points)), k)
                errors.append(core.compression_errors(points, points[indices],
                                                       full, baseline))
        values = {}
        for metric in METRICS:
            finite = [float(e[metric]) for e in errors if e[metric] is not None]
            must(all(math.isfinite(v) and v >= 0 for v in finite),
                 'Nonfinite/negative reconstructed compression error')
            values[metric] = ({'mean': float(np.mean(finite)),
                               'q90': float(np.quantile(finite, .9)),
                               'maximum': float(np.max(finite))}
                              if finite else None)
        result[str(cap)] = values
    return result


def augment_exact_births(caps, reps):
    """SHA-audit 120 physical chunks, stream all raw births, preserve originals."""
    prior, source, _summary = readout.checked(campaign.ROOT)
    expected = {(r['geometry'], int(r['box_id']), int(r['sample'])): r
                for r in prior}
    must(len(expected) == len(prior) == source['n_selected'] == 960,
         'Unexpected full-cloud source census')
    annotated = []
    seen = set()
    for chunk_id in range(120):
        folder = core.chunk_path(chunk_id)
        prefix = core.prefix(folder)
        trials = list(core.csv_records(Path(str(prefix)+core.SUFFIX[0])))
        moments = list(core.csv_records(Path(str(prefix)+core.SUFFIX[3])))
        must(len(trials) == len(moments) == core.GROUP,
             'Native full-cloud metadata malformed in chunk '+str(chunk_id))
        bytrial = {core.trial_key(r): r for r in trials}
        bymoment = {core.trial_key(r): r for r in moments}
        must(len(bytrial) == len(bymoment) == core.GROUP and
             set(bytrial) == set(bymoment),
             'Missing/duplicated trial or moments '+str(chunk_id))
        nonzero = set()
        for label, points in core.full_cloud_iterator(
                Path(str(prefix)+core.SUFFIX[1])):
            must(label in bytrial and label not in nonzero,
                 'Unknown or repeated ion birth source')
            trial = bytrial[label]
            must(int(trial['ion_births']) > 0,
                 'Unexpected birth CSV for zero-ion native trial')
            core.verify_moments(points, bymoment[label], trial)
            key = (trial['configuration'], int(trial['box_id']),
                   int(trial['sample']))
            must(key in expected and key not in seen,
                 'Source identity was not preserved: '+str(key))
            old = expected[key]
            must(len(points) == old['n_ion_births'] and
                 key[0] in campaign.GEOMS and '1024' in old['cap_study'],
                 'Baseline/full-cloud birth mismatch: '+str(key))
            new_caps = score_cloud(points, caps, reps, *key)
            annotated.append({**old, 'cap_study': {**old['cap_study'], **new_caps}})
            seen.add(key)
            nonzero.add(label)
        for label, trial in bytrial.items():
            if label in nonzero:
                continue
            must(int(trial['ion_births']) == 0,
                 'Missing full births for nonzero avalanche')
            key = (trial['configuration'], int(trial['box_id']),
                   int(trial['sample']))
            must(key in expected and key not in seen and
                 expected[key]['n_ion_births'] == 0,
                 'Zero-ion trial not represented correctly')
            annotated.append(expected[key])
            seen.add(key)
    must(seen == set(expected) and len(annotated) == 960,
         'Not all 960 exact-start avalanches examined')
    return annotated


def emit(rows, caps, bootstrap, seed, reps):
    rng = np.random.default_rng(seed)
    print('PASS all 120 native full-cloud chunks SHA-audited; old analyses read-only')
    print('P3 HIGH-CAP ION-BIRTH STORAGE COMPRESSION, NOT A TRAINED LUT OR SIGNAL')
    print('New caps:', ', '.join(map(str, caps)), '; reservoir draws per avalanche:', reps)
    print('Uncertainty unit: independent D500 spatial cells, NOT ions/draws')
    for geom in campaign.GEOMS:
        active = [r for r in rows if r['geometry'] == geom and r['n_ion_births'] > 0]
        n_total = sum(r['geometry'] == geom for r in rows)
        print('\n'+'='*68)
        print(geom, 'total starts', n_total, 'nonzero', len(active),
              'zero', n_total-len(active), 'full ion births',
              sum(r['n_ion_births'] for r in active))
        for cap in (1024,)+tuple(caps):
            print('CAP', cap, 'genuinely compressed',
                  sum(r['n_ion_births'] > cap for r in active), '/', len(active))
        for regime in (('near_anode','transition','bulk') if geom == 'p3_noFW'
                       else ('near_anode','near_field','transition','bulk')):
            subset = [r for r in active if r['regime'] == regime]
            print('REGIME', regime, 'n', len(subset), 'spatial cells',
                  len(set(r['box_id'] for r in subset)))
            for cap in (1024,)+tuple(caps):
                all_ = readout.clustered(subset,cap,'joint_phi_time_TV',bootstrap,rng)
                compressed = readout.clustered(
                    subset,cap,'joint_phi_time_TV',bootstrap,rng,compressed_only=True)
                print(' TV',cap,readout.fmt(all_),
                      '; compressed-only',readout.fmt(compressed))
        print('ALL REGIMES: cell-bootstrap errors')
        for metric in METRICS:
            print(' METRIC',metric)
            for cap in (1024,)+tuple(caps):
                x=readout.clustered(active,cap,metric,bootstrap,rng)
                print('  ',cap,readout.fmt(x))
        print('ION-BIRTH WEIGHTED (sensitivity only)')
        for metric in ('joint_phi_time_TV','mean_radius_rel_abs','abs_q90_t_ns'):
            print(' METRIC',metric)
            for cap in (1024,)+tuple(caps):
                x=readout.clustered(active,cap,metric,bootstrap,rng,
                                    ion_birth_weight=True)
                print('  ',cap,readout.fmt(x))
        for low,high in zip((1024,)+tuple(caps[:-1]),caps):
            print('PAIRED IMPROVEMENT',low,'->',high,'(same full avalanche)')
            for metric in ('joint_phi_time_TV','mean_radius_rel_abs','abs_q90_t_ns'):
                overall=readout.paired(active,low,high,metric,bootstrap,rng)
                compressed=readout.paired(active,low,high,metric,bootstrap,rng,
                                          compressed_only=True)
                print(' ',metric,readout.fmt(overall),
                      '; only N_birth >',high,readout.fmt(compressed))
    print('\nCAUTION: precision threshold unspecified; no cap accepted by default.')
    print('No new full Garfield, joint LUT training, ion drift or cathode weighting.')
    print('PASS high-cap read-only study; frozen files untouched')


def selftest():
    toy=np.asarray([[.1+.002*j,(-1)**j*.05,0.,j*.001]
                    for j in range(300)],dtype=float)
    for caps in ((2048,),(2048,4096)):
        out=score_cloud(toy,caps,8,'p3_noFW',123,124)
        must(set(out)==set(map(str,caps)),'Missing toy high-cap output')
        for cap in caps:
            must(all(v is None or abs(v['mean'])<1e-10
                     for v in out[str(cap)].values()),
                 'Exact no-compression identity failed')
    larger=np.asarray([[.1+(.002*j%1.2),(-1)**j*(.05+j*.0002),
                        .001*j,j*.0001] for j in range(2200)],dtype=float)
    a=score_cloud(larger,(2048,),8,'p3_FW75',200,100)
    b=score_cloud(larger,(2048,),8,'p3_FW75',200,100)
    must(a==b and a['2048']['joint_phi_time_TV']['mean']>0,
         'Deterministic actual-compression selftest failed')
    print('PASS high-cap exact identity and deterministic real compression; '
          'no Garfield/files written')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode',choices=('selftest','report'))
    parser.add_argument('--caps',nargs='+',type=int,default=list(DEFAULT_CAPS))
    parser.add_argument('--reps',type=int,default=32)
    parser.add_argument('--bootstrap',type=int,default=3000)
    parser.add_argument('--seed',type=int,default=20261104)
    args=parser.parse_args()
    caps=tuple(args.caps)
    must(caps and caps==tuple(sorted(set(caps))) and set(caps).issubset(ALLOWED_CAPS),
         'Choose increasing unique caps from 2048,4096,8192')
    must(8<=args.reps<=128,'Require 8..128 reservoirs per avalanche')
    must(200<=args.bootstrap<=10000,'Require 200..10000 cell bootstrap draws')
    if args.mode=='selftest':
        selftest()
    else:
        data=augment_exact_births(caps,args.reps)
        emit(data,caps,args.bootstrap,args.seed,args.reps)


if __name__=='__main__':
    main()
