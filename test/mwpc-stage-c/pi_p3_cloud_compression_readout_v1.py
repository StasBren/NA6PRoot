#!/usr/bin/env python3
"""Read-only forensic report for TWO p=3mm complete-cloud pilot geometries.

Checks complete analysis-v2 SHA, six physical checkpoint seals, source identity,
then reports 32/128/512 birth-reservoir fidelity and parent-level bootstrap CIs.
Never calls Garfield; never overwrites or creates physics/analysis data.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from collections import Counter
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parents[2]
ROOT = HERE / 'test_runs/mwpc_stage_c/pi_p3_full_cloud_convergence_v1'
GEOMS = ('p3_noFW', 'p3_FW75')
CAPS = (32, 128, 512)
METRICS = (
    ('joint_phi_time_TV', '12x6 phi-time total-variation'),
    ('abs_A_TB', 'absolute top-bottom ion-birth asymmetry error'),
    ('abs_top_frac', 'absolute top-sector fraction error'),
    ('mean_radius_rel_abs', 'relative mean ion-birth radius error'),
    ('abs_q90_t_ns', 'absolute 90% birth-time quantile error [ns]'),
)


def must(cond, msg):
    if not cond:
        raise ValueError(msg)


def load(path):
    with path.open(encoding='utf-8') as fh:
        return json.load(fh)


def sha(path):
    hasher = hashlib.sha256()
    with path.open('rb') as fh:
        for buf in iter(lambda: fh.read(1 << 20), b''):
            hasher.update(buf)
    return hasher.hexdigest()


def checked(root):
    analysis = root / 'analysis_complete_v2'
    receipt = load(analysis / 'manifest.json')
    must(receipt.get('schema') == 'NA60_P3_FULL_CLOUD_ANALYSIS_SHA_V2',
         'Not the fixed, sealed complete-cloud analysis v2')
    must(set(receipt.get('output_sha256', {})) == {'summary.json', 'per_avalanche.json'},
         'Unexpected analysis output set')
    for name, expected in receipt['output_sha256'].items():
        must(sha(analysis / name) == expected, 'Analysis output SHA mismatch: ' + name)
    source = load(root / 'input/selection_manifest.json')
    must(source['schema'] == 'NA60_P3_CLOUD_INPUT_V1' and
         source['per_geometry'] == 24 and source['n_selected'] == 48,
         'Unexpected full-cloud pilot source geometry/cardinality')
    must(sha(root / 'input/exact_D500_p3_unseen.csv') == source['input_sha256'],
         'Frozen D500 electron start CSV changed')
    must(receipt['selection_sha256'] == source['input_sha256'],
         'Analysis and frozen D500 selection disagree')
    must(set(receipt['input_chunk_sha256']) ==
         {f'chunk_{i:04d}/complete.json' for i in range(6)},
         'Expected exactly six immutable full-cloud checkpoints')
    total = 0
    for i in range(6):
        name = f'chunk_{i:04d}'
        cd = root / 'full_cloud_chunks' / name
        seal = cd / 'complete.json'
        must(sha(seal) == receipt['input_chunk_sha256'][name + '/complete.json'],
             'Physical checkpoint receipt was replaced: ' + name)
        record = load(seal)
        must(record['schema'] == 'NA60_P3_FULL_CLOUD_SHA_CHUNK_V1' and
             record['chunk'] == i and record['n_avalanches'] == 8,
             'Malformed sealed native Garfield checkpoint: ' + name)
        for file_name, hexdigest in record['sha256_files'].items():
            must(sha(cd / file_name) == hexdigest,
                 'ORIGINAL full Garfield birth-cloud SHA mismatch: ' + name + '/' + file_name)
        total += record['n_ion_births']
    summary = load(analysis / 'summary.json')
    must(summary.get('schema') == 'NA60_P3_FULL_BIRTH_COMPRESS_ANALYSIS_V2' and
         summary['selected_starts'] == 48 and summary['completed_chunks'] == 6,
         'Wrong/partial cloud compression analysis')
    must(set(summary['summary_by_geometry']) == set(GEOMS),
         'A non-p3 geometry entered pilot')
    cloud = load(analysis / 'per_avalanche.json')
    must(len(cloud) == 48, 'Not all 48 independent electron-start avalanches analyzed')
    must(len({(v['geometry'], v['box_id'], v['sample']) for v in cloud}) == 48,
         'Reused D500 microscopic electron-start ID')
    must(sum(v['n_ion_births'] for v in cloud) == total,
         'Missing exact full microscopic ion births in compression analysis')
    for g in GEOMS:
        sub = [v for v in cloud if v['geometry'] == g]
        must(len(sub) == 24 and
             sum(v['n_ion_births'] for v in sub) ==
             summary['summary_by_geometry'][g]['total_exact_ion_births'],
             'Geometry count/full-cloud mismatch: ' + g)
        for r in sub:
            must(r['n_ion_births'] >= 0, 'Negative ion birth count')
            if not r['n_ion_births']:
                must(r['exact_A_TB'] is None and not r['cap_study'],
                     'Zero-ion avalanche has fabricated shape')
            else:
                must(set(r['cap_study']) == {'32','128','512'},
                     'Nonzero avalanche missing reservoir-size comparison')
    return cloud, summary, total


def stats(values, bootstrap, rng):
    """Bootstrap across INDEPENDENT microscopic avalanches, not 64 subsamples."""
    vals = np.asarray(values, dtype=float)
    if len(vals) == 0:
        return None
    must(bool(np.all(np.isfinite(vals))), 'Nonfinite cap-compression measure')
    out = {'n_avalanches':len(vals), 'mean':float(vals.mean()),
           'median':float(np.median(vals)), 'q90':float(np.quantile(vals, .9))}
    if len(vals) >= 2:
        idx = rng.integers(len(vals), size=(bootstrap,len(vals)))
        out['bootstrap95_across_avalanches'] = list(map(float,
                      np.quantile(vals[idx].mean(axis=1), [.025,.975])))
    else:
        out['bootstrap95_across_avalanches'] = None
    return out


def values_for(rows, cap, metric, *, informative=False):
    result=[]
    for r in rows:
        n=r['n_ion_births']
        if not n or (informative and n <= cap):
            continue
        measure = r['cap_study'][str(cap)][metric]
        if measure is not None:
            result.append(float(measure['mean']))
    return result


def report(root, bootstrap, seed):
    must(200 <= bootstrap <= 10000, 'Require 200-10000 bootstrap draws')
    cloud, meta, n_births = checked(root)
    rng=np.random.default_rng(seed)
    print('PASS SHA: all six immutable full-cloud physics checkpoints and analysis_complete_v2')
    print('FULL MICROSCOPIC: 48 independent D500 electron starts,', n_births,
          'EXACT ionization births; NO signal response')
    print('TV=0 indicates identical coarse angular/time probabilities; 1 means disjoint distributions.')
    print('All scores below are PER AVALANCHE, not weighted by its ion count.')
    print('64 subsamples are repeated measurements of the SAME avalanche, not 64 independent avalanches.')
    print('\nBASIC AND BY-REGIME:')
    for g in GEOMS:
        arr=[v for v in cloud if v['geometry']==g]
        active=[v for v in arr if v['n_ion_births']>0]
        missing=sum(v['exact_A_TB'] is None for v in active)
        print('\n'+('='*77))
        print(g+f': 24 starts; {len(active)} nonzero, {24-len(active)} zero; '
              f'undefined A_TB among nonzero = {missing}')
        print('births=',sum(x['n_ion_births'] for x in arr),
              '  births/avalanche mean=',round(float(np.mean([r['n_ion_births'] for r in arr])),2))
        for cap in CAPS:
            n_full=sum(v['n_ion_births']<=cap for v in active)
            n_undef=sum(v['exact_A_TB'] is not None and
                        v['cap_study'][str(cap)]['abs_A_TB'] is None
                        for v in active)
            print(f'  CAP {cap}: {n_full}/{len(active)} nonzero avalanches do not need any compression; '
                  f'{len(active)-n_full} genuinely compressed; '
                  f'{n_undef} originally-defined A_TB lost in compressed draw(s)')
        for regime in ('near_anode','near_field','transition','bulk'):
            rr=[r for r in arr if r['regime']==regime]
            if not rr:continue
            print('\n  '+regime+f': {len(rr)} independent avalanches, '
                  f'{sum(v["n_ion_births"]>0 for v in rr)} nonzero')
            for cap in CAPS:
                vals=values_for(rr,cap,'joint_phi_time_TV')
                v=stats(vals,bootstrap,rng)
                vals_st=values_for(rr,cap,'joint_phi_time_TV',informative=True)
                inf=stats(vals_st,bootstrap,rng)
                print(f'    cap={cap:3d} TV mean={v["mean"] if v else None!s:>19} '
                      f'CI95={v["bootstrap95_across_avalanches"] if v else None} '
                      f'[compressed-only n={len(vals_st)}, mean={round(inf["mean"],4) if inf else None}]')
        print('\n  ALL REGIMES PER CAP (unweighted across independent nonzero avalanches):')
        for label,metric in METRICS:
            print('    '+label)
            for cap in CAPS:
                value=stats(values_for(arr,cap,metric),bootstrap,rng)
                compressed=stats(values_for(arr,cap,metric,informative=True),bootstrap,rng)
                if value:
                    ci=value['bootstrap95_across_avalanches']
                    print(f'      {cap:3d}: {value["mean"]:.6g} [95% CI {ci[0]:.6g},{ci[1]:.6g}] '
                          f'n={value["n_avalanches"]}; truly compressed {compressed["n_avalanches"] if compressed else 0}')
                else: print(f'      {cap:3d}: no defined observations')
        # Paired benefit from more ion samples; each input avalanche is one
        # independent unit, not each random subsampling repetition.
        for metric in ('joint_phi_time_TV','abs_A_TB','abs_q90_t_ns'):
            for only_compressed in (False, True):
                pairs=[]
                for r in active:
                    if only_compressed and r['n_ion_births']<=512:
                        continue
                    a,b=r['cap_study']['128'][metric],r['cap_study']['512'][metric]
                    if a is not None and b is not None:
                        pairs.append(a['mean']-b['mean'])
                v=stats(pairs,bootstrap,rng)
                if v:
                    print(f'  Improvement 512 relative to 128, {metric} '
                          f'({"Nbirth>512 only" if only_compressed else "all defined"}): '
                          f'+{v["mean"]:.6g}, 95% CI {v["bootstrap95_across_avalanches"]}, '
                          f'paired n={len(pairs)}')
    print('\nCAUTION: 18 or 19 nonzero starts/geometry are far too few for p99 tails,')
    print('fine-grained spatial generalization, or a definitive 128-versus-512 decision.')
    print('All angle/time binning was fixed using the FULL cloud of each avalanche;')
    print('this test checks compression only, NOT a new conditioned generator.')
    print('PASS read-only SHA-sealed p3 cloud compression report. No files created or overwritten.')


def selftest():
    import tempfile
    z=np.random.default_rng(0)
    raw=[{'n_ion_births':100,'cap_study':{'32':{'joint_phi_time_TV':{'mean':.15}},
                                              '128':{'joint_phi_time_TV':{'mean':0.}},
                                              '512':{'joint_phi_time_TV':{'mean':0.}}}},
         {'n_ion_births':0,'cap_study':{}},
         {'n_ion_births':600,'cap_study':{'32':{'joint_phi_time_TV':{'mean':.25}},
                                             '128':{'joint_phi_time_TV':{'mean':.11}},
                                             '512':{'joint_phi_time_TV':{'mean':.02}}}}]
    must(values_for(raw,32,'joint_phi_time_TV')==[.15,.25], 'Nonzero compression rows')
    must(values_for(raw,128,'joint_phi_time_TV',informative=True)==[.11],
         'Uncompressed avalanche must not count as informative 128 comparison')
    a=stats([0.,1.,2.],300,z)
    must(a['mean']==1. and a['n_avalanches']==3, 'Bootstrap independence')
    with tempfile.TemporaryDirectory() as temp:
        p=Path(temp)/'file.json'
        p.write_text('{"x":1}\n')
        must(sha(p)==hashlib.sha256(p.read_bytes()).hexdigest(),
             'SHA checksum')
    print('PASS synthetic reservoir cap strata, per-avalanche (NOT per-ion) bootstrap, and SHA. No Garfield physics.')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('mode',choices=('selftest','report'))
    p.add_argument('--root',type=Path,default=ROOT)
    p.add_argument('--bootstrap',type=int,default=3000)
    p.add_argument('--seed',type=int,default=20261104)
    a=p.parse_args()
    if a.mode=='selftest':selftest()
    else:report(a.root,a.bootstrap,a.seed)

if __name__=='__main__':main()
