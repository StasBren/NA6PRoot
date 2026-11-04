#!/usr/bin/env python3
"""P3-only complete ion-BIRTH cloud learning curve: FIRST 5 starts per cell.

Isolated NEW 9000-start Garfield++ campaign (900 cells * 5 * 2 p3 fields).
Original A20, A100, D500, 48- and 960-cloud archives are all SHA-frozen
read-only. Exact initial D500 sites, full ion (u,v,w,t) output, never drift
endpoints or subsampled/weighted ion cloud. Status tolerates one active chunk.
Modes: selftest, prepare, preflight, run, status, audit.
No gamma fitting, no cathode signal, no implicit automatic large run.
"""
from __future__ import annotations

import argparse
import collections
import hashlib
import json
import math
from pathlib import Path

import pi_p3_extended_full_cloud_campaign_v1 as old
from pi_audit_avalanche_lut_reuse import BASE, CONFIGS, EXPORT_FIELDS, locations, sha
from pi_run_lut500_overnight import source_maps, offset_for, STRIDE

core = old.core
ROOT = BASE / 'pi_p3_fullcloud_learning5_v1'
INPUT = ROOT / 'input'
CHUNKS = ROOT / 'full_cloud_chunks'
GEOMS = old.GEOMS
GROUP = old.GROUP
FIELDS = tuple(EXPORT_FIELDS)
SCHEMA = 'NA60_P3_FULL_CLOUD_LEARNING5_SOURCE_V1'
TAG = 'NA60-P3-FULLCLOUD-LEARNING5-PASS140-219-20261104-V1'
N_PER_CELL = 5
N_CELLS = 900
PASSES = range(140, 220)  # never A20/A100 [0,100) or extended pilot [100,140)
FUTURE_PHYSICS_HOLDOUT = range(460, 500)  # RESERVE; do not train on these


def must(ok, why):
    if not ok:
        raise ValueError(why)


def choose_seed(*parts):
    msg = '|'.join(map(str, (TAG,) + parts))
    return int.from_bytes(hashlib.sha256(msg.encode()).digest()[:8], 'big')


def sample_for_pass(g, box, pass_index):
    """Original 500-per-cell D500 permutation, valid for ANY pass 0..499.

    The earlier extended-cloud module's old.sample_for() intentionally
    restricts its API to passes 100..139. DO NOT call that helper for this
    separate campaign's passes 140..219 or its reserved blind holdouts.
    """
    must(g in GEOMS and 0 <= box < N_CELLS and 0 <= pass_index < 500,
         'Invalid original D500 geometry, cell or permutation pass')
    return (offset_for(g, box) + STRIDE * pass_index) % 500


def reserved(g, box):
    """80 unused original D500 permutation indices, no electron rerun."""
    must(g in GEOMS and 0 <= box < N_CELLS, 'Invalid original D500 cell')
    mapping = {sample_for_pass(g, box, p): p for p in PASSES}
    must(len(mapping) == len(PASSES), 'D500 reserved source passes collide')
    return mapping


def immutable_sources():
    """Cross-check source identities of the 48 and 960 old complete clouds."""
    old_keys, old48sha = old.old_pilot_keys()
    older_rows, _, older_meta = old.check_plan()
    old_ext = {(r['geometry'], int(r['box_id']), int(r['sample']))
               for r in older_rows}
    must(len(old_ext) == len(older_rows) == 960 and not
         old_ext.intersection(old_keys), 'Overlapping old full-cloud source')
    return old_keys, old_ext, old48sha, older_meta['input_sha256']


def source_selection():
    old_keys, ext_keys, old48sha, extsha = immutable_sources()
    index, index_sha = source_maps()
    selected = {g: collections.defaultdict(list) for g in GEOMS}
    source_shas = {}
    stats = {}
    for g, pitch, fw in CONFIGS:
        if g not in GEOMS:
            continue
        _, filename = locations(500, pitch, fw)
        must(sha(filename) == index[g]['lut_csv_sha256'],
             'Original SHA-frozen D500 source changed: ' + g)
        source_shas[g] = sha(filename)
        candidate_samples = [reserved(g, box) for box in range(N_CELLS)]
        seen = set()
        rejected_metal = 0
        rejected_oldpilot = 0
        with filename.open(newline='') as stream:
            reader = core.csv.DictReader(stream)
            expected_fields = {
                'box_id', 'sample', 'iu', 'iv', 'u_start_mm', 'v_start_mm',
                'w_start_mm', 'random_seed', 'drift_ok', 'collected',
                'field_wire_collected', 'wire_index', 'phi_hit_deg',
                'drift_time_ns', 'last_step_energy_eV'}
            must(expected_fields.issubset(reader.fieldnames or []),
                 'Missing original D500 source columns')
            for r in reader:
                box = int(r['box_id'])
                sample = int(r['sample'])
                must(0 <= box < N_CELLS, 'Invalid original D500 source box')
                passage = candidate_samples[box].get(sample)
                if passage is None:
                    continue
                k = (box, passage)
                must(k not in seen, 'Duplicate original D500 source pass')
                seen.add(k)
                u = float(r['u_start_mm'])
                v = float(r['v_start_mm'])
                w = float(r['w_start_mm'])
                iu = int(r['iu'])
                iv = int(r['iv'])
                must(box == iu*30+iv and math.isfinite(u) and
                     math.isfinite(v) and math.isfinite(w) and
                     abs(w) < 1.e-8 and
                     -1.5 + iu*.1 - 1.e-6 <= u <= -1.5 + (iu+1)*.1 + 1.e-6 and
                     -3.0 + iv*.2 - 1.e-6 <= v <= -3.0 + (iv+1)*.2 + 1.e-6,
                     'Original D500 site outside its claimed spatial cell')
                key = (g, box, sample)
                if key in old_keys or key in ext_keys:
                    rejected_oldpilot += 1
                    continue
                label = core.stratum(u, v, g)
                if label is None:
                    rejected_metal += 1
                    continue
                data = {
                    'lut_trials': 500, 'geometry': g,
                    'original_lut_sha256': index[g]['lut_csv_sha256'],
                    'box_id': box, 'iu': iu, 'iv': iv, 'sample': sample,
                    'random_seed': int(r['random_seed']),
                    'u_start_mm': r['u_start_mm'],
                    'v_start_mm': r['v_start_mm'],
                    'w_start_mm': r['w_start_mm'],
                    'old_drift_ok': r['drift_ok'],
                    'old_collected': r['collected'],
                    'old_field_wire_collected': r['field_wire_collected'],
                    'old_wire_index': r['wire_index'],
                    'old_phi_hit_deg': r['phi_hit_deg'],
                    'old_drift_time_ns': r['drift_time_ns'],
                    'old_last_step_energy_eV': r['last_step_energy_eV']
                }
                must(set(data) == set(FIELDS), 'Unexpected D500 source schema')
                selected[g][box].append((choose_seed(g, box, passage),
                                         passage, data, label))
        must(len(seen) == N_CELLS * len(PASSES),
             'Missing one or more reserved D500 start sites: ' + g)
        for box in range(N_CELLS):
            must(len(selected[g][box]) >= N_PER_CELL,
                 'Insufficient original gas starts in ' + g + '/' + str(box))
            selected[g][box].sort(key=lambda x: (x[0], x[1]))
            selected[g][box] = selected[g][box][:N_PER_CELL]
        stats[g] = {
            'candidate_starts': len(seen),
            'excluded_oldpilot': rejected_oldpilot,
            'excluded_inside_metal': rejected_metal,
            'chosen': N_CELLS*N_PER_CELL
        }
    must(set(selected) == set(GEOMS), 'Some geometry missing')
    # For every round, include all 900 cells and both geometries;
    # SHA ordering distributes source locations, including near-wire.
    box_order = sorted(range(N_CELLS), key=lambda b: (choose_seed('box-order',b),b))
    rows, tags = [], []
    for rep in range(N_PER_CELL):
        for box in box_order:
            for g in GEOMS:
                _, passage, data, label = selected[g][box][rep]
                rows.append(data)
                tags.append({
                    'geometry': g,
                    'box_id': box,
                    'sample': int(data['sample']),
                    'pass_index': passage,
                    'replica': rep,
                    'regime': label[0],
                    'v_sign': label[1]
                })
    must(len(rows) == len(tags) == 2*N_CELLS*N_PER_CELL and
         len(rows)%GROUP == 0,
         'Learning5 source count is not 9000 or multiple of chunk size')
    return rows, tags, stats, source_shas, index_sha, old48sha, extsha


def check_plan():
    rows = list(core.csv_records(INPUT / 'exact_D500_p3_unseen.csv', FIELDS))
    tags = core.read_json(INPUT / 'strata.json')
    m = core.read_json(INPUT / 'selection_manifest.json')
    must(m.get('schema') == SCHEMA and
         len(rows) == len(tags) == m.get('n_selected') == 9000 and
         m.get('per_geometry') == 4500 and
         m.get('group_size') == GROUP and
         m.get('source_passes') == [140, 219] and
         tuple(m.get('geometries', [])) == GEOMS and
         sha(INPUT / 'exact_D500_p3_unseen.csv') == m['input_sha256'] and
         sha(INPUT / 'strata.json') == m['strata_sha256'],
         'Learning5 source selection or SHA receipt changed')
    old_keys, ext_keys, old48sha, extsha = immutable_sources()
    must(old48sha == m['original_48pilot_input_sha256'] and
         extsha == m['original_960pilot_input_sha256'],
         'An old complete-cloud source identity changed')
    unique = set()
    counter = collections.Counter()
    for i, (r, label) in enumerate(zip(rows,tags)):
        g = r['geometry']
        box = int(r['box_id'])
        s = int(r['sample'])
        passage = int(label['pass_index'])
        k = (g,box,s)
        must(g == GEOMS[i%2] and label['geometry'] == g and
             label['box_id'] == box and label['sample'] == s and
             label['replica'] in range(N_PER_CELL) and passage in PASSES and
             sample_for_pass(g,box,passage) == s and
             core.stratum(float(r['u_start_mm']),float(r['v_start_mm']),g) ==
                 (label['regime'], label['v_sign']) and
             k not in unique and k not in old_keys and k not in ext_keys,
             'Learning5 reused, embedded-metal, or mislabelled D500 start ' + str(k))
        unique.add(k)
        counter[g,box] += 1
    must(len(unique) == 9000 and
         all(counter[g,box] == N_PER_CELL for g in GEOMS
             for box in range(N_CELLS)),
         'Not exactly five independent source starts per D500 cell')
    return rows, tags, m


def configure_core():
    # Reuse original full-Garfield immutable chunk/seal/preflight code,
    # but make all mutable paths point ONLY to a new independent folder.
    core.ROOT = ROOT
    core.INPUT = INPUT
    core.CHUNKS = CHUNKS
    core.FIELDS = FIELDS
    core.RNG_TAG = TAG
    core.check_plan = check_plan


def prepare():
    must(not ROOT.exists(), 'Campaign directory already exists; refuse overwrite: '+str(ROOT))
    rows, tags, stats, shas, indexsha, old48sha, extsha = source_selection()
    INPUT.mkdir(parents=True)  # only now create new campaign files
    core.write_csv(INPUT / 'exact_D500_p3_unseen.csv', rows)
    core.write_json(INPUT / 'strata.json', tags)
    core.write_json(INPUT / 'selection_manifest.json', {
        'schema': SCHEMA, 'n_selected': len(rows),
        'per_geometry': len(rows)//len(GEOMS),
        'n_per_cell': N_PER_CELL, 'cells_per_geometry': N_CELLS,
        'group_size': GROUP, 'geometries': GEOMS,
        'source_passes': [min(PASSES), max(PASSES)],
        'reserved_future_validation_passes': [460,499],
        'original_48pilot_input_sha256': old48sha,
        'original_960pilot_input_sha256': extsha,
        'original_D500_sha256': shas,
        'original_D500_index_manifest_sha256': indexsha,
        'selection_note': 'Full original D500 gas-phase electron start; 5 per cell. '
            'A20/A100 0-99, ext 100-139 excluded; 140-219 candidates; '
            'future blind validation 460-499 reserved, not simulated. '
            'Exact ion BIRTH cloud (u,v,w,t) and ion census; no cathode signal.',
        'selection_counts': stats,
        'B_T': 0, 'anode_voltage_V': 1800, 'cathode_gap_mm': 3,
        'gas': 'Ar70CO2_30_293.15K_760T',
        'input_sha256': sha(INPUT / 'exact_D500_p3_unseen.csv'),
        'strata_sha256': sha(INPUT / 'strata.json')
    })
    check_plan()
    print('PASS FROZEN new full-cloud P3 learning5 inputs: 9000 D500 original starts')
    for g in GEOMS:
        print('GEOMETRY',g,'exactly 5 fresh source electrons in EVERY 900 spatial cells')
        print('SOURCE', stats[g])
    print('NO Garfield physics run by prepare')


def selftest():
    must(GEOMS == ('p3_noFW','p3_FW75') and GROUP == 8 and
         len(PASSES) == 80 and N_PER_CELL == 5,
         'Unexpected phase-B full-cloud driver interface')
    for g in GEOMS:
        for box in (0,1,450,899):
            current = set(reserved(g,box))
            older = {sample_for_pass(g,box,p) for p in range(140)}
            future = {sample_for_pass(g,box,p) for p in FUTURE_PHYSICS_HOLDOUT}
            all500 = {sample_for_pass(g,box,p) for p in range(500)}
            must(len(all500) == 500 and len(current) == 80 and
                 len(older) == 140 and len(future) == 40 and
                 not current.intersection(older) and
                 not current.intersection(future) and
                 not older.intersection(future),
                 'D500 permutation/source-pass collision or reuse')
            for pass_index in (100,111,125,139):
                must(sample_for_pass(g,box,pass_index) ==
                     old.sample_for(g,box,pass_index),
                     'D500 absolute permutation differs from old SHA-frozen helper')
    for g in GEOMS:
        must(core.classify(0.,0.,g) is None,'Anode metal incorrectly accepted')
        must(core.classify(1.5,0.,'p3_FW75') is None,
             'Field-wire metal incorrectly accepted')
    must(len(set(range(N_CELLS))) == N_CELLS,
         'Invalid uniform cell enumeration')
    print('PASS P3 learning5 disjoint original D500 pass ranges and metal guards')
    print('NO source directory created, NO Garfield physics')


def preflight():
    core.preflight()


def audit():
    rows, _, _, done = core.completed()
    require_n = len(rows)//GROUP
    must(len(done) == require_n,
         'Cannot audit incomplete physical campaign: '+
         str(len(done))+'/'+str(require_n))
    must(sum(x['n_avalanches'] for x in done) == len(rows) == 9000,
         'Missing true full Garfield avalanches')
    print('PASS all 1125 physics chunks SHA-sealed, 9000 full Garfield avalanches')
    print('TOTAL EXACT ION BIRTHS', sum(x['n_ion_births'] for x in done))


def status():
    rows, _, _ = check_plan()
    binary = core.verified_binary()
    groups = core.batch(rows)
    ids, pending = old.scan_chunk_state([core.chunk_path(i) for i in range(len(groups))])
    # Audits the SHA of every SEALED chunk, but not an unsealed writer chunk.
    done = [core.audit_chunk(i,groups[i],binary) for i in ids]
    completed = sum(x['n_avalanches'] for x in done)
    must(completed == GROUP*len(ids), 'Non-contiguous sealed source census')
    print('PASS new learning5 P3 physics status: SHA-sealed avalanches',
          completed,'/',len(rows))
    print('VERIFIED CHUNKS',len(ids),'/',len(groups))
    print('EXACT ION BIRTH POINTS',sum(x['n_ion_births'] for x in done))
    print('BALANCED per geometry',completed//2)
    if pending is not None:
        print('UNSEALED CHUNK',pending,'possibly ACTIVE or INTERRUPTED; NOT counted')
        print('Check current writer process and physics.log before re-running.')
        print('Never start two writers, delete a partial folder or change the compiled binary.')
    elif len(ids) < len(groups):
        print('NEXT UNUSED CHUNK',len(ids))
    else:
        print('PASS ALL 1125 SHA-SEALED FULL CLOUD CHUNKS')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('action', choices=('selftest','prepare','preflight','run','status','audit'))
    p.add_argument('--max-chunks',type=int,default=1)
    p.add_argument('--hours',type=float,default=2.)
    a = p.parse_args()
    configure_core()
    if a.action == 'selftest':selftest()
    elif a.action == 'prepare':prepare()
    elif a.action == 'preflight':preflight()
    elif a.action == 'run':core.run(a.max_chunks,a.hours)
    elif a.action == 'status':status()
    else:audit()


if __name__ == '__main__':
    main()
