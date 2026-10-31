#!/usr/bin/env python3
"""Strict recovery and continuation of the NEW 120-parent Garfield holdout.

Never edit existing microscopic CSVs, logs, old 120 reference, LUT or C++.
Accept only geometrically proven zero-ion failures starting within the
30um anode (or old permitted 75um field wire).
Modes: selftest, audit, recover, run, evaluate, status.
"""
from __future__ import annotations
import argparse
import json
import math
import types
from collections import Counter
from pathlib import Path

import pi_external_muon120_guarded_lut_v1 as ext
from pi_run_muon_lut_validation import field_wire_blocked_zero_start

legacy = ext.legacy
CHUNK = 7
ANODE_R_MM = .015
INNER_MARGIN_MM = 1.e-6


def anode_blocked_zero_start(r):
    """Inside one of nine 30um-diameter anodes, with a true no-ion endpoint."""
    if r['configuration'] not in ext.GEOMS:
        return False
    if (int(r['successful']) != 0 or int(r['capped']) != 0 or
        int(r['garfield_ne']) != 1 or int(r['garfield_ni']) != 0 or
        int(r['ion_births']) != 0 or int(r['collected_e']) != 1 or
        int(r['attached']) != 0 or int(r['other']) != 0 or
        any(int(r[k]) != 0 for k in
            ('ion_up', 'ion_down', 'ion_top', 'ion_bottom'))):
        return False
    x, y = float(r['seed_u_mm']), float(r['seed_v_mm'])
    if not (math.isfinite(x) and math.isfinite(y)):
        return False
    pitch = 3. if r['configuration'].startswith('p3_') else 4.
    nearest = legacy.round_away(x/pitch)
    return (-4 <= nearest <= 4 and
            math.hypot(x-nearest*pitch, y) < ANODE_R_MM-INNER_MARGIN_MM)


def strict_metal_start(r):
    return field_wire_blocked_zero_start(r) or anode_blocked_zero_start(r)


def isolated_validator():
    """Reuse original bytecode and all conservation checks, with LOCAL metal test.

    Crucially: no mutation of canonical validator globals, any original source
    file, or the frozen old120/A20/A100 protocol.
    """
    original = legacy.validate_reference
    names = dict(original.__globals__)
    if 'field_wire_blocked_zero_start' not in names:
        raise RuntimeError('Original validator does not have expected metal hook')
    names['field_wire_blocked_zero_start'] = strict_metal_start
    strict = types.FunctionType(original.__code__, names, original.__name__,
                                original.__defaults__, original.__closure__)

    def validate(path, selected):
        stat = strict(path, selected)
        seeds = legacy.read_csv(Path(path)/'muon_reference_seed.csv',
                                legacy.SEED_FIELDS)
        fw, anode = Counter(), Counter()
        for r in seeds:
            if int(r['successful']) == 1:
                continue
            if field_wire_blocked_zero_start(r):
                fw[r['configuration']] += 1
            elif anode_blocked_zero_start(r):
                anode[r['configuration']] += 1
            else:
                raise ValueError('Unexpected nonmetal Garfield unsuccessful seed')
        if sum(fw.values())+sum(anode.values()) != stat[
                'n_geometrically_verified_field_wire_zero_starts']:
            raise ValueError('Metal census disagrees with original conservation validator')
        stat.update({
            'n_geometrically_verified_field_wire_zero_starts':sum(fw.values()),
            'field_wire_zero_starts_by_geometry':dict(fw),
            'n_geometrically_verified_anode_zero_starts':sum(anode.values()),
            'anode_zero_starts_by_geometry':dict(anode),
            'note_on_unsuccessful':
                'Only geometrically demonstrated zero-ion starts strictly inside '
                '30um anode or 75um field-wire metal are permitted; original '
                'successful=0 and all ion/endpoint counters are unchanged.'
        })
        return stat
    return validate


def inspect_chunk():
    rows, meta = ext.check()
    binary = legacy.check_executable()
    if (ext.REF/'manifest.json').exists():
        raise FileExistsError('Final external full reference already sealed')
    cfgfile = ext.REF/'run_config.json'
    if not cfgfile.is_file():
        raise FileNotFoundError('Missing original full Garfield run config')
    cfg = json.loads(cfgfile.read_text())
    if (cfg.get('input_sha256') != meta['sha256_input'] or
        cfg.get('cpp_source_sha256') != legacy.sha(legacy.CXX) or
        cfg.get('cpp_executable_sha256') != binary):
        raise ValueError('Frozen full native input/source/binary SHA mismatch')
    batches = legacy.groups(rows)
    for i in range(CHUNK):
        legacy.read_seal(legacy.chunk_dir(i), i, batches[i], binary)
    c = legacy.chunk_dir(CHUNK)
    if not c.is_dir() or (c/'complete.json').exists():
        raise RuntimeError('Expected unsealed completed chunk_07 missing/already sealed')
    inp = c/'native_muon_crossings_input.csv'
    actual = legacy.read_csv(inp, ext.FIELDS)
    if ([legacy.id_of(x) for x in actual] !=
            [legacy.id_of(x) for x in batches[CHUNK]] or
        legacy.sha(inp) != legacy.sha_bytes_for_subset(ext.FIELDS,
                                                       batches[CHUNK])):
        raise ValueError('Original eighth group physical muon input changed')
    raw = [inp, c/'physics.log'] + [
        c/name for name in legacy.FIELDS_TO_MERGE]
    if not all(p.is_file() and p.stat().st_size for p in raw):
        raise RuntimeError('Original native physical outputs missing/empty')
    physics_log = (c/'physics.log').read_text(errors='replace')
    for marker in ('PASS: muon-driven B1 all crossing/geometry rows written',
                   'PASS: immutable same-Heed primary/cluster cache'):
        if marker not in physics_log:
            raise RuntimeError('Full native C++ physics did NOT complete: '+marker)
    stat = isolated_validator()(c,batches[CHUNK])
    if (stat['n_native_muons'] != 10 or stat['n_heed_primaries'] < 1 or
        stat['n_full_avalanche_trials'] != 4*stat['n_heed_primaries'] or
        any(int(x['n_capped_seeds']) for x in stat['by_geometry'].values())):
        raise ValueError('Native full Garfield census inconsistent or capped')
    if stat['n_geometrically_verified_anode_zero_starts'] < 1:
        raise ValueError('Documented anode-start example not found in physical CSV')
    return c, raw, stat, binary, batches[CHUNK]


def audit():
    c,raw,stat,_,_ = inspect_chunk()
    print('PASS forensic read-only audit: all native physics finished, '
          'original total/sector/cluster and metal census valid',flush=True)
    print(json.dumps({
        'chunk':'chunk_07','n_real_muons':stat['n_native_muons'],
        'n_Heed_primaries':stat['n_heed_primaries'],
        'n_full_micro_avalanches':stat['n_full_avalanche_trials'],
        'anode_metal_zero_starts':stat['anode_zero_starts_by_geometry'],
        'field_wire_metal_zero_starts':stat['field_wire_zero_starts_by_geometry'],
        'unmeasured_wall_time':'UNKNOWN after original Python validation abort',
        'source_file_sha256':{p.name:legacy.sha(p) for p in raw},
    },indent=2),flush=True)


def recover():
    c,raw,stat,binary,_ = inspect_chunk()
    record = {
        'schema':'PI_10_REAL_GEANT4_MUON_GARFIELD_SHA_CHECKPOINT_V1',
        'status':'complete','chunk':CHUNK,'n_muons':10,
        'cpp_executable_sha256':binary,'input_sha256':legacy.sha(raw[0]),
        'stats':stat,
        'recovery':{
            'schema':'NA60_EXTERNAL_120_STRICT_ANODE_FIELD_METAL_V1',
            'original_microscopic_physics_rerun':False,
            'original_raw_files_modified':False,
            'physics_wall_seconds':'UNKNOWN: validator aborted before timing seal',
            'accepted':'Only physically inside metal, zero ion, uncapped, '
                       'complete original microscopic conservation',
        },
        'outputs_sha256':{p.name:legacy.sha(p) for p in raw},
    }
    with (c/'complete.json').open('x') as f:
        f.write(json.dumps(record,indent=2,sort_keys=True)+'\n')
    rows,_=ext.check()
    legacy.read_seal(c,CHUNK,legacy.groups(rows)[CHUNK],binary)
    print('PASS RECOVERED eighth block: new SHA seal only, '
          'original Garfield files unchanged; no physics rerun',flush=True)
    print('anode starts',stat['anode_zero_starts_by_geometry'],
          'field-wire starts',stat['field_wire_zero_starts_by_geometry'],flush=True)


def selftest():
    r={'configuration':'p3_noFW','seed_u_mm':'.00237593484574722',
       'seed_v_mm':'.0107657602328344','successful':'0','capped':'0',
       'garfield_ne':'1','garfield_ni':'0','ion_births':'0',
       'collected_e':'1','attached':'0','other':'0',
       'ion_up':'0','ion_down':'0','ion_top':'0','ion_bottom':'0'}
    assert anode_blocked_zero_start(r)
    assert not field_wire_blocked_zero_start(r)
    assert not anode_blocked_zero_start(dict(r,ion_births='1'))
    assert not anode_blocked_zero_start(dict(r,seed_u_mm='0.02'))
    assert not anode_blocked_zero_start(dict(r,seed_v_mm='0.02'))
    assert not anode_blocked_zero_start(dict(r,capped='1'))
    assert not anode_blocked_zero_start(dict(r,collected_e='0'))
    assert not anode_blocked_zero_start(dict(r,successful='1'))
    original=legacy.validate_reference
    old=original.__globals__['field_wire_blocked_zero_start']
    assert callable(isolated_validator())
    assert original.__globals__['field_wire_blocked_zero_start'] is old
    print('PASS strict 30um anode geometry/failed ion-census '
          'and isolated-original-validator selftests; no physics',flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('mode',choices=['selftest','audit','recover',
                                  'run','evaluate','status'])
    p.add_argument('--hours',type=float,default=8.)
    p.add_argument('--max-chunks',type=int,default=12)
    p.add_argument('--reps',type=int,default=3000)
    p.add_argument('--seed',type=int,default=20261101)
    a=p.parse_args()
    if a.mode=='selftest':selftest()
    elif a.mode=='audit':audit()
    elif a.mode=='recover':recover()
    elif a.mode=='status':ext.status()
    elif a.mode=='run':
        ext.check()
        legacy.validate_reference=isolated_validator()
        legacy.run(a.hours,a.max_chunks)
    elif a.mode=='evaluate':
        ext.check()
        legacy.validate_reference=isolated_validator()
        ext.evaluate(a.reps,a.seed)


if __name__=='__main__':main()
