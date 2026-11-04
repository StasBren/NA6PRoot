#!/usr/bin/env python3
"""P3 read-only full-ion-BIRTH library and spatially-disjoint sampling pilot.

Separate from SHA-sealed native Garfield archives and the validated A20 LUT.
Stores whole (u,v,w,t) clouds without subsampling and never changes old data.
Modes: selftest, build, audit, evaluate. No new Garfield, Heed, ion drift,
weighting field, or cathode signal. First evaluation is exploratory only.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
import random
import shutil
import tempfile
from collections import defaultdict
from pathlib import Path

import numpy as np

import pi_p3_extended_full_cloud_campaign_v1 as campaign
import pi_p3_extended_cloud_readout_v1 as readout

CORE = campaign.core
GEOMS = campaign.GEOMS
DEST = campaign.ROOT.parent / 'pi_p3_joint_fullcloud_bank_v1'
TAG = 'NA60-P3-JOINT-FULLCLOUD-SPATIAL-HOLDOUT-20261104-V1'


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    return campaign.sha(path)


def dump(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n')


def read(path):
    return json.loads(path.read_text())


def sealed_sources():
    """Audit physical and analysis SHA seals before interpreting any points."""
    cloud, source, analysis = readout.checked(campaign.ROOT)
    # readout.checked also calls configure_core and audits all 120 chunks.
    rows, tags, selection, chunks = CORE.completed()
    require(len(rows) == len(cloud) == len(tags) == 960 and len(chunks) == 120,
            'Unexpected incomplete P3 native source')
    return rows, tags, cloud, selection


def build():
    rows, tags, analysis_rows, selection = sealed_sources()
    require(not DEST.exists(), 'Output bank already exists: ' + str(DEST))
    tmp = DEST.with_name(DEST.name + '.INCOMPLETE')
    require(not tmp.exists(), 'Preserve unsealed partial bank for inspection: ' + str(tmp))
    require(shutil.disk_usage(DEST.parent).free > 1024**3,
            'Less than 1 GiB free for immutable binary bank')
    tmp.mkdir(parents=True)
    lookup = {(r['geometry'], int(r['box_id']), int(r['sample'])): r
              for r in analysis_rows}
    require(len(lookup) == len(analysis_rows), 'Duplicate analyzed avalanche identity')
    indices = {g: [] for g in GEOMS}
    mapping = {}
    totals = {g: 0 for g in GEOMS}
    for src, label in zip(rows, tags):
        g = src['geometry']
        key = (g, int(src['box_id']), int(src['sample']))
        require(g in GEOMS and key in lookup, 'Invalid D500 start identity')
        outcome = lookup[key]
        n = int(outcome['n_ion_births'])
        require(n >= 0, 'Negative number of ionization births')
        item = {
            'geometry': g, 'box_id': key[1], 'sample': key[2],
            'regime': label['regime'], 'v_sign': label['v_sign'],
            'd500_seed': int(src['random_seed']),
            'u0_mm': float(src['u_start_mm']),
            'v0_mm': float(src['v_start_mm']),
            'w0_mm': float(src['w_start_mm']),
            'n_ion_births': n, 'offset': totals[g],
            'crossing_id': f'LUT500|{g}|{key[1]}|{key[2]}',
        }
        require(CORE.stratum(item['u0_mm'], item['v0_mm'], g) ==
                (item['regime'], item['v_sign']), 'D500 geometry tag changed')
        indices[g].append(item)
        mapping[key] = item
        totals[g] += n
    require(len(mapping) == 960, 'Unexpected source identities')
    require(sum(totals.values()) == 2981533, 'Exact physical ion census changed')
    shapes = {g: (totals[g], 4) for g in GEOMS}
    arrays = {g: np.lib.format.open_memmap(
        tmp / ('clouds_' + g + '.npy'), mode='w+', dtype='<f8', shape=shapes[g])
        for g in GEOMS}
    seen = set()
    for j in range(120):
        folder = campaign.CHUNKS / f'chunk_{j:04d}'
        prefix = folder / 'cloud'
        trials = list(CORE.csv_records(Path(str(prefix) + '_trials.csv')))
        moments = list(CORE.csv_records(Path(str(prefix) + '_moments.csv')))
        bytrial = {CORE.trial_key(x): x for x in trials}
        bymoment = {CORE.trial_key(x): x for x in moments}
        require(len(bytrial) == len(bymoment) == len(trials) == len(moments) == 8
                and set(bytrial) == set(bymoment), 'Duplicate/missing trial or moment')
        for native_key, points in CORE.full_cloud_iterator(
                Path(str(prefix) + '_birth_reservoir.csv')):
            require(native_key in bytrial, 'Unknown native ion source')
            trial = bytrial[native_key]
            key = (trial['configuration'], int(trial['box_id']), int(trial['sample']))
            require(key in mapping and key not in seen, 'Duplicate or unselected cloud')
            item = mapping[key]
            require(native_key[0] == item['crossing_id'] and
                    len(points) == item['n_ion_births'] > 0 and
                    int(trial['primary_index']) == 0,
                    'Native avalanche identity/count does not agree with selection')
            CORE.verify_moments(points, bymoment[native_key], trial)
            start = item['offset']
            arrays[key[0]][start:start + len(points)] = points
            seen.add(key)
        for native_key, trial in bytrial.items():
            key = (trial['configuration'], int(trial['box_id']), int(trial['sample']))
            require(key in mapping and native_key[0] == mapping[key]['crossing_id'],
                    'Mismatched trial/crossing ID')
            n = mapping[key]['n_ion_births']
            require(n == int(trial['ion_births']) and
                    ((key in seen) == (n > 0)),
                    'Unaccounted zero/nonzero Garfield avalanche')
        if (j+1) % 20 == 0:
            print('VERIFIED native chunks', j+1, '/120', flush=True)
    require(len(seen) == sum(r['n_ion_births'] > 0 for g in GEOMS
                             for r in indices[g]), 'Nonzero cloud missing')
    for a in arrays.values():
        a.flush()
    arrays.clear()
    for g in GEOMS:
        stored = np.load(tmp / ('clouds_' + g + '.npy'), mmap_mode='r')
        require(stored.shape == shapes[g] and np.all(np.isfinite(stored)),
                'Nonfinite or truncated binary ion cloud')
    dump(tmp / 'index.json', {'schema': 'NA60_P3_JOINT_FULLCLOUD_INDEX_V1',
                             'geometries': indices,
                             'source_note': 'Independent D500-start Garfield full ion-BIRTH clouds, NOT ion drift, not muon flux'})
    data_files = ['index.json'] + ['clouds_' + g + '.npy' for g in GEOMS]
    dump(tmp / 'manifest.json', {
        'schema': 'NA60_P3_JOINT_FULLCLOUD_SHA_BANK_V1',
        'n_avalanches': 960, 'n_ion_births': totals,
        'source_selection_sha256': selection['input_sha256'],
        'source_analysis_manifest_sha256': sha(campaign.ROOT / 'analysis_complete_v2/manifest.json'),
        'files_sha256': {name: sha(tmp / name) for name in data_files},
        'note': 'Full clouds, zero clouds kept. Timing relative to Garfield start t=0. No new physics.'})
    tmp.rename(DEST)
    print('PASS bank built without compression', DEST, 'births=', totals, flush=True)


def audit():
    rows, tags, physical, selection = sealed_sources()
    m = read(DEST / 'manifest.json')
    idx = read(DEST / 'index.json')
    require(m['schema'] == 'NA60_P3_JOINT_FULLCLOUD_SHA_BANK_V1' and
            idx['schema'] == 'NA60_P3_JOINT_FULLCLOUD_INDEX_V1' and
            m['source_selection_sha256'] == selection['input_sha256'] and
            m['source_analysis_manifest_sha256'] ==
            sha(campaign.ROOT / 'analysis_complete_v2/manifest.json'),
            'Bank source manifest changed')
    for name, digest in m['files_sha256'].items():
        require(sha(DEST / name) == digest, 'Changed bank data: ' + name)
    expect = {(r['geometry'], int(r['box_id']), int(r['sample'])): r
              for r in physical}
    require(set(idx['geometries']) == set(GEOMS), 'Geometry set changed')
    found = set()
    for g in GEOMS:
        arr = np.load(DEST / ('clouds_' + g + '.npy'), mmap_mode='r')
        offset = 0
        for r in idx['geometries'][g]:
            key = (g, r['box_id'], r['sample'])
            require(key in expect and key not in found and
                    r['n_ion_births'] == expect[key]['n_ion_births']
                    and r['offset'] == offset, 'Changed cloud index/source count')
            found.add(key)
            offset += r['n_ion_births']
        require(arr.shape == (offset, 4) and offset == m['n_ion_births'][g]
                and np.all(np.isfinite(arr)), 'Corrupt array or cloud offsets')
        print('SHA PASS', g, len(idx['geometries'][g]),
              'complete avalanches /', offset, 'births', flush=True)
    require(found == set(expect) and
            sum(map(len, idx['geometries'].values())) == len(expect) == 960,
            'Incomplete bank')
    print('PASS physical SHA + full-cloud bank audit', flush=True)
    return idx


class Bank:
    """Empirical whole-cloud sampler. Never rescale ion count or birth times."""
    def __init__(self, path=DEST):
        p = Path(path)
        self.items = read(p / 'index.json')['geometries']
        self.points = {g: np.load(p / ('clouds_' + g + '.npy'), mmap_mode='r')
                       for g in self.items}

    @staticmethod
    def holdout(g, cell):
        k = hashlib.sha256(f'{TAG}|{g}|{cell}'.encode()).digest()
        return int.from_bytes(k[:8], 'big') % 5 == 0

    @staticmethod
    def guards(u, v, g):
        require(g in GEOMS, 'Unknown P3 geometry')
        require(-1.5 <= u < 1.5 and -3 <= v < 3, 'Outside 900-cell P3 local domain')
        region = CORE.classify(u, v, g)
        if region is None:
            return None  # physical metal: must not create a gas avalanche
        side = 'positive' if v >= 0 else 'negative'
        # Near wires: protect opposite azimuthal sides, not just v hemisphere.
        anode_side = 1 if u >= 0 else -1
        fw_x = 1.5 if u >= 0 else -1.5
        fw_side = 1 if u >= fw_x else -1
        return region, side, anode_side, fw_x, fw_side

    def candidates(self, u, v, g, mode='guarded', train=None):
        q = self.guards(u, v, g)
        if q is None:
            return []
        require(mode in ('local', 'nearest', 'guarded'), 'Unknown sampler')
        source = self.items[g] if train is None else train
        iu = int(math.floor((u + 1.5) / .1))
        iv = int(math.floor((v + 3.0) / .2))
        current = iu * 30 + iv
        groups = defaultdict(list)
        for r in source:
            if r['regime'] != q[0] or r['v_sign'] != q[1]:
                continue
            du = int(r['box_id']) // 30 - iu
            dv = int(r['box_id']) % 30 - iv
            if abs(du) > 2 or abs(dv) > 2:
                continue
            if q[0] == 'near_anode' and (1 if r['u0_mm'] >= 0 else -1) != q[2]:
                continue
            if q[0] == 'near_field':
                sx = 1.5 if r['u0_mm'] >= 0 else -1.5
                ss = 1 if r['u0_mm'] >= sx else -1
                if sx != q[3] or ss != q[4]:
                    continue
            groups[int(r['box_id'])].append(r)
        if mode == 'local':
            groups = {current: groups[current]} if current in groups else {}
        if mode == 'nearest' and groups:
            best = min(groups, key=lambda box: (
                (box // 30 - iu)**2 + (box % 30 - iv)**2, box))
            groups = {best: groups[best]}
        if not groups:
            return []
        if mode in ('nearest', 'local'):
            weights = {box: 1.0 / len(groups) for box in groups}
        else:
            bybox = {box: math.exp(-.5 * ((box//30-iu)**2 + (box%30-iv)**2))
                     for box in groups}
            if current in groups and len(groups) > 1:
                local_weight = .25
                other = sum(v for b, v in bybox.items() if b != current)
                weights = {b: (local_weight if b == current else
                               (1-local_weight)*bybox[b]/other) for b in groups}
            else:
                total = sum(bybox.values())
                weights = {b: x/total for b, x in bybox.items()}
        result = []
        for box, rs in groups.items():
            for r in rs:
                result.append((r, weights[box]/len(rs)))
        z = sum(p for _, p in result)
        require(z > 0 and math.isfinite(z), 'Invalid conditional bank weights')
        return [(r, p/z) for r, p in result]

    def sample(self, u0, v0, w0, t0, g, rng, mode='guarded', train=None):
        choices = self.candidates(u0, v0, g, mode=mode, train=train)
        if not choices:
            return None  # ABSTAIN: missing coverage or physical metal
        picked = rng.choices([r for r, _ in choices],
                             weights=[p for _, p in choices], k=1)[0]
        n = picked['n_ion_births']
        start = picked['offset']
        cloud = np.array(self.points[g][start:start+n], dtype=np.float64, copy=True)
        # Native full clouds start at t=0 and w0=0. Do NOT add D500 drift time.
        cloud[:, 2] += w0 - picked['w0_mm']
        cloud[:, 3] += t0
        return {'n_ion_births': n, 'births_uvwt_mm_ns': cloud,
                'source_cell': picked['box_id'],
                'source_sample': picked['sample'],
                'mode': mode, 'physical_basis': 'replayed one complete Garfield birth cloud'}


def crps_discrete(outcome, support):
    """Exact proper CRPS for discrete predictive support with weights."""
    if not support:
        return None
    x = np.asarray([float(v) for v, p in support])
    p = np.asarray([float(p) for v, p in support])
    p = p / np.sum(p)
    return float(np.dot(p, np.abs(x - float(outcome))) -
                 .5 * np.dot(p, np.abs(x[:, None] - x[None, :]) @ p))


def cloud_features(record, points):
    """Event-level joint features of a WHOLE cloud; zeros remain physical zeros."""
    n = record['n_ion_births']
    if n == 0:
        return np.zeros(7, dtype=float)
    a = points[record['offset']:record['offset'] + n]
    stats = CORE.cloud_stats(a)
    return np.array([math.log1p(n), stats['mean_radius_mm'],
                     stats['sd_radius_mm'], stats['mean_t_ns'],
                     stats['q90_t_ns'], stats['mean_w_mm'],
                     stats['top_frac']], dtype=float)


def energy_score(y, support):
    """Proper multivariate energy score on pre-scaled whole-cloud features."""
    if not support:
        return None
    x = np.stack([v for v, p in support])
    p = np.asarray([p for v, p in support], dtype=float)
    p /= p.sum()
    first = float(np.dot(p, np.linalg.norm(x - y, axis=1)))
    distances = np.linalg.norm(x[:, None] - x[None, :], axis=2)
    return first - .5 * float(p @ distances @ p)


def evaluate():
    audit()
    b = Bank()
    print('READ-ONLY spatial-cell holdout; sparse, stratified D500 source, NOT muon flux')
    print('This is NOT a check of full 4D joint distribution or U/V signal.')
    for g in GEOMS:
        train = [r for r in b.items[g] if not b.holdout(g, r['box_id'])]
        test = [r for r in b.items[g] if b.holdout(g, r['box_id'])]
        require(train and test and not ({r['box_id'] for r in train} &
                                        {r['box_id'] for r in test}),
                'Leaking a spatial cell into both train and test')
        # Scales fixed using TRAINING CELLS ONLY; no held-out cloud leakage.
        features = {(r['box_id'], r['sample']): cloud_features(r, b.points[g])
                    for r in b.items[g]}
        train_features = np.stack([features[r['box_id'], r['sample']]
                                   for r in train])
        scale = np.quantile(train_features, .75, axis=0) - \
                np.quantile(train_features, .25, axis=0)
        scale = np.where(scale > 1.e-5, scale, 1.)
        print('\nGEOMETRY', g, 'train', len(train), 'test', len(test),
              'held-out spatial cells', len({r['box_id'] for r in test}))
        preds = {mode: defaultdict(list) for mode in ('nearest', 'guarded')}
        missing = {mode: 0 for mode in preds}
        for truth in test:
            for mode in preds:
                cand = b.candidates(truth['u0_mm'], truth['v0_mm'], g,
                                    mode=mode, train=train)
                if not cand:
                    missing[mode] += 1
                    continue
                n = truth['n_ion_births']
                support = [(r['n_ion_births'], p) for r, p in cand]
                pzero = sum(p for r, p in cand if r['n_ion_births'] == 0)
                truth_z = features[truth['box_id'], truth['sample']] / scale
                support_z = [(features[r['box_id'], r['sample']] / scale, p)
                             for r, p in cand]
                preds[mode][truth['box_id']].append((
                    crps_discrete(n, support),
                    (pzero - float(n == 0))**2,
                    energy_score(truth_z, support_z)))
        common = set(preds['nearest']) & set(preds['guarded'])
        for mode in preds:
            cells = preds[mode]
            n_ok = sum(map(len, cells.values()))
            if not n_ok:
                print(mode, 'NO SUPPORTED QUERIES; no score')
                continue
            crps = float(np.mean([np.mean([x[0] for x in arr])
                                  for arr in cells.values()]))
            brier = float(np.mean([np.mean([x[1] for x in arr])
                                   for arr in cells.values()]))
            energy = float(np.mean([np.mean([x[2] for x in arr])
                                    for arr in cells.values()]))
            print(mode, 'coverage', n_ok, '/', len(test),
                  'tested cells', len(cells), 'abstained', missing[mode],
                  'cell-balanced count CRPS', round(crps, 5),
                  'zero-Brier', round(brier, 6),
                  'joint-feature energy score', round(energy, 6))
        if common:
            delta = np.array([
                np.mean([x[0] for x in preds['nearest'][cell]]) -
                np.mean([x[0] for x in preds['guarded'][cell]])
                for cell in sorted(common)])
            energy_delta = np.array([
                np.mean([x[2] for x in preds['nearest'][cell]]) -
                np.mean([x[2] for x in preds['guarded'][cell]])
                for cell in sorted(common)])
            rng = np.random.default_rng(20261104)
            sample = rng.integers(0, len(delta), size=(2000, len(delta)))
            ci = list(map(float, np.quantile(delta[sample].mean(axis=1), [.025, .975])))
            energy_ci = list(map(float, np.quantile(
                energy_delta[sample].mean(axis=1), [.025, .975])))
            print('PAIRED nearest-minus-guarded count CRPS',
                  float(np.mean(delta)), 'common held-out cells', len(common),
                  '(positive favors guarded; exploratory)',
                  'cell-bootstrap CI95', ci)
            print('PAIRED nearest-minus-guarded joint-feature energy score',
                  float(np.mean(energy_delta)), 'CI95', energy_ci)
        else:
            print('No cells supported by both samplers')
    print('\nPASS read-only evaluation; NO Garfield calls, NO LUT acceptance')


def selftest():
    with tempfile.TemporaryDirectory(prefix='p3_bank_selftest_') as temp:
        p = Path(temp)
        np.save(p/'clouds_p3_noFW.npy', np.array([
            [0.01, 0.02, 0., 4.], [0.02, 0.02, 0., 4.5]], dtype='<f8'))
        np.save(p/'clouds_p3_FW75.npy', np.empty((0, 4), dtype='<f8'))
        item = {'geometry': 'p3_noFW', 'box_id': 465, 'sample': 101,
                'regime': 'near_anode', 'v_sign': 'positive',
                'u0_mm': .02, 'v0_mm': .01, 'w0_mm': 0.,
                'offset': 0, 'n_ion_births': 2}
        dump(p/'index.json', {'geometries': {
            'p3_noFW': [item], 'p3_FW75': []}})
        bank = Bank(p)
        sample = bank.sample(.02, .01, 3., 20., 'p3_noFW', random.Random(42),
                             mode='local')
        require(sample is not None and sample['n_ion_births'] == 2,
                'Selftest whole-cloud selection failed')
        require(np.allclose(sample['births_uvwt_mm_ns'][:, 2], 3.) and
                np.allclose(sample['births_uvwt_mm_ns'][:, 3], [24., 24.5]),
                'Birth clocks or w translation incorrect')
        require(bank.sample(.02, .01, 0., 0., 'p3_noFW', random.Random(42),
                            mode='local', train=[]) is None,
                'No-data query must abstain')
        require(bank.sample(0., 0., 0., 0., 'p3_noFW', random.Random(42)) is None,
                'Metal starting point must not generate gas ionization')
        require(abs(energy_score(np.array([1.]),
                     [(np.array([0.]), .5), (np.array([2.]), .5)])-.5)<1e-12,
                'Energy-score selftest wrong')
        require(abs(crps_discrete(1, [(0, .5), (2, .5)])-.5)<1e-12,
                'Discrete proper scoring formula wrong')
        require(Bank.holdout('p3_noFW', 50) ==
                Bank.holdout('p3_noFW', 50), 'Spatial split not deterministic')
    print('PASS P3 bank selftest: full cloud, clock, w, metal, abstain, CRPS')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('selftest', 'build', 'audit', 'evaluate'))
    args = parser.parse_args()
    {'selftest': selftest, 'build': build,
     'audit': audit, 'evaluate': evaluate}[args.action]()


if __name__ == '__main__':
    main()
