#!/usr/bin/env python3
"""P3 complete-ion-BIRTH bank: independent-cell and within-cell diagnostics.

Read-only and SHA-verified. Does not change native physics or existing bank.
Run selftest then report. Both geometries independently; ion count and 7D
joint whole-cloud feature energy score are NOT a full 4D point-process test.
"""
from __future__ import annotations
import argparse
import hashlib
from collections import defaultdict

import numpy as np

import pi_p3_joint_fullcloud_bank_v1 as bankmod

MODELS = ('nearest', 'guarded')
METRICS = ('count_crps', 'pzero_brier', 'joint_energy')
CV_TAG = 'NA60-P3-FULLCLOUD-5FOLD-AND-WITHINCELL-20261104-V1'


def require(ok, why):
    if not ok:
        raise ValueError(why)


def fold(g, cell):
    h = hashlib.sha256(f'{CV_TAG}|spatial|{g}|{cell}'.encode('utf-8')).digest()
    return int.from_bytes(h[:8], 'big') % 5


def start_order(g, row):
    h = hashlib.sha256(
        f'{CV_TAG}|within|{g}|{row["box_id"]}|{row["sample"]}'.encode('utf-8')
    ).digest()
    return int.from_bytes(h[:8], 'big')


def feature_cache(b, g):
    out = {}
    for r in b.items[g]:
        key = (r['box_id'], r['sample'])
        require(key not in out, 'Duplicate full-cloud metadata')
        out[key] = bankmod.cloud_features(r, b.points[g])
    return out


def training_scale(training, features):
    require(training, 'Empty training')
    a = np.stack([features[(r['box_id'], r['sample'])] for r in training])
    width = np.quantile(a, .75, axis=0) - np.quantile(a, .25, axis=0)
    return np.where(width > 1.e-5, width, 1.)


def predict_one(b, g, truth, train, feats, scale, mode):
    choices = b.candidates(truth['u0_mm'], truth['v0_mm'], g,
                           mode=mode, train=train)
    if not choices:
        return None
    support = [(r['n_ion_births'], p) for r, p in choices]
    zprob = sum(p for r, p in choices if r['n_ion_births'] == 0)
    y = (truth['n_ion_births'] == 0)
    joint = [
        (feats[(r['box_id'], r['sample'])]/scale, p)
        for r, p in choices
    ]
    target = feats[(truth['box_id'], truth['sample'])]/scale
    scores = [
        bankmod.crps_discrete(truth['n_ion_births'], support),
        (zprob - float(y))**2,
        bankmod.energy_score(target, joint),
    ]
    require(np.all(np.isfinite(scores)), 'Undefined conditional scores')
    return tuple(float(x) for x in scores)


def summarize(label, percell, rng, boots):
    """Every test CELL one unit. A cell with multiple independent starts
    first contributes its average loss; no ion-level pseudo-replication."""
    models = sorted({m for cell in percell.values() for m in cell})
    allcells = sorted(percell)
    print(f'  {label} test cells={len(allcells)}')
    for mode in models:
        cells = [c for c in allcells if mode in percell[c] and percell[c][mode]]
        if not cells:
            print('   ', mode, 'no supported cells')
            continue
        z = np.array([
            np.mean(percell[c][mode], axis=0) for c in cells], dtype=float)
        print('   ', mode, 'supported_cells', len(cells),
              'supported_avalanches', sum(len(percell[c][mode]) for c in cells),
              'abstained_cells', len(allcells)-len(cells),
              'mean', {k:float(v) for k,v in zip(METRICS,z.mean(axis=0))})
    for left, right in [('nearest', 'guarded'), ('local', 'guarded')]:
        common = [c for c in allcells if
                  left in percell[c] and right in percell[c] and
                  percell[c][left] and percell[c][right]]
        if not common:
            continue
        differences = np.stack([
            np.mean(percell[c][left], axis=0) -
            np.mean(percell[c][right], axis=0) for c in common])
        ix = rng.integers(0, len(common), size=(boots,len(common)))
        boot = differences[ix].mean(axis=1)
        lo = np.quantile(boot, .025, axis=0)
        hi = np.quantile(boot, .975, axis=0)
        print('    paired',left,'minus',right,'on',len(common),'cells')
        for i, metric in enumerate(METRICS):
            print('      ',metric,'delta',float(differences[:,i].mean()),
                  'CI95',[float(lo[i]),float(hi[i])])


def independent_spatial_folds(b, g, feats, boots, rng):
    print(' SPATIAL 5-FOLD: each box held out exactly once, no same-cell training')
    allcells = {r['box_id'] for r in b.items[g]}
    require(allcells, 'No cells')
    percell = defaultdict(lambda: defaultdict(list))
    n_by_fold = []
    for f in range(5):
        train = [r for r in b.items[g] if fold(g,r['box_id']) != f]
        test = [r for r in b.items[g] if fold(g,r['box_id']) == f]
        require(train and test, 'Empty held-out fold')
        require(not ({x['box_id'] for x in train} &
                     {x['box_id'] for x in test}), 'Spatial leakage')
        scale = training_scale(train, feats)
        n_by_fold.append((f, len(train), len(test),
                          len({r['box_id'] for r in test})))
        for truth in test:
            for mode in MODELS:
                val = predict_one(b, g, truth, train, feats, scale, mode)
                if val is not None:
                    percell[truth['box_id']][mode].append(val)
                else:
                    # Explicitly retain unsupported cells in denominator.
                    _ = percell[truth['box_id']]
    require(set(percell) == allcells, 'Holdout did not cover every spatial cell')
    print('  fold_counts',n_by_fold)
    print('  Note: CV training folds overlap; CIs condition on these partitions,')
    print('  and do not fully include the uncertainty from constructing the training bank.')
    summarize('SPATIAL', percell, rng, boots)


def within_cell(b,g,feats,boots,rng):
    print(' WITHIN-CELL: one independent trial held out in each cell with >=2 starts')
    groups = defaultdict(list)
    for r in b.items[g]:
        groups[r['box_id']].append(r)
    chosen = {
        cell: min(rows, key=lambda r: start_order(g,r))
        for cell,rows in groups.items() if len(rows) >= 2
    }
    print('  source cells total',len(groups),
          'cells with >=2 independent starts',len(chosen),
          'heldout avalanches',len(chosen))
    if not chosen:
        return
    excluded = {(r['box_id'],r['sample']) for r in chosen.values()}
    train = [r for r in b.items[g] if (r['box_id'],r['sample']) not in excluded]
    scale = training_scale(train, feats)
    bycell = defaultdict(lambda: defaultdict(list))
    for cell, truth in chosen.items():
        # Historical local-only baseline vs candidate guarded whose local
        # 25% component is actually present; unlike the spatial CV test.
        for mode in ('local','guarded'):
            outcome = predict_one(b,g,truth,train,feats,scale,mode)
            if outcome is not None:
                bycell[cell][mode].append(outcome)
            else:
                _ = bycell[cell]
    summarize('WITHIN-CELL',bycell,rng,boots)


def report(boots):
    require(100 <= boots <= 20000, 'Choose 100..20000 cell bootstrap resamples')
    bankmod.audit()
    b = bankmod.Bank()
    print('P3 COMPLETE BIRTH CLOUD LUT: 5 SPATIAL FOLDS + WITHIN CELL')
    print('Scores: lower better, independently for each geometry, no U/V signal.')
    print('No new microscopic calculation; bank and source manifests verified.')
    rng = np.random.default_rng(20261104)
    for g in bankmod.GEOMS:
        print('\n===',g,'===')
        feats = feature_cache(b,g)
        independent_spatial_folds(b,g,feats,boots,rng)
        within_cell(b,g,feats,boots,rng)
    print('\nPASS read-only 5-fold spatial & independent-within-cell diagnostics')


def selftest():
    assert fold('p3_noFW',250) == fold('p3_noFW',250)
    assert 0 <= fold('p3_FW75',100) < 5
    assert start_order('p3_noFW', {'box_id':3,'sample':7}) == \
           start_order('p3_noFW', {'box_id':3,'sample':7})
    rng = np.random.default_rng(20261104)
    data = {1:{'nearest':[(2.,.1,3.)], 'guarded':[(1.,.05,2.)]},
            2:{'nearest':[(1.,.1,1.)], 'guarded':[(2.,.2,2.)]},
            3:{'nearest':[], 'guarded':[]}}
    summarize('TOY',data,rng,200)
    print('PASS P3 fullcloud crossval synthetic selftest')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('action',choices=('selftest','report'))
    p.add_argument('--bootstrap',type=int,default=3000)
    a=p.parse_args()
    (selftest if a.action == 'selftest' else lambda:report(a.bootstrap))()


if __name__ == '__main__':
    main()
