#!/usr/bin/env python3
"""Isolated Phase-C V2: compare native RKF ESTIMATED ion population with 100
existing, matched full microscopic Garfield++ ion-birth trials.

V1 source, gas, exact-start manifest and all original physics are READ ONLY.
No Geant4, Heed, or AvalancheMicroscopic computation is performed here.
This is NOT independent muon validation or a microscopic ion-census replica.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
import subprocess
from collections import defaultdict
from pathlib import Path
from typing import Optional

from pi_prepare_muon_avalanche_b1 import REPO, sha
from pi_A100_targeted_50cell_campaign import ROOT as A100, GEOMS
from pi_run_rkf_polya_pilot import (
    ROOT as V1_ROOT, SOURCE as V1_SOURCE,
    GAS as V1_GAS, PLAN as V1_PLAN, OUTPUT as V1_OUTPUT,
    INPUT_SUBSET as V1_INPUT, RKF_COLUMNS as V1_COLUMNS,
    read_csv, verify_plan as verify_v1_plan,
)

ROOT = A100.parent / 'pi_RKF_ion_census_B0_pilot_v2'
CXX = REPO / 'test/mwpc-avalanche/phaseC_rkf_ion_census_pilot_v2.cxx'
EXE = REPO / 'test/mwpc-avalanche/build/mwpc_phase_c_rkf_ion_census_pilot_v2'
OUTPUT = ROOT / 'rkf_estimated_ions_100_exact_starts.csv'
RECEIPT = ROOT / 'run_receipt.json'
PER_SEED = ROOT / 'matched_100_micro_vs_RKF_ions.csv'
SUMMARY = ROOT / 'summary.json'
MANIFEST = ROOT / 'comparison_manifest.json'

COLUMNS = V1_COLUMNS[:10] + ('rkf_electrons_estimate', 'rkf_ions_estimate') + V1_COLUMNS[10:]
METADATA_FIELDS = ("geometry", "box_id", "sample")
COMPARE_FIELDS = (
    'geometry', 'box_id', 'sample', 'start_u_mm', 'start_v_mm',
    'micro_ion_births', 'micro_garfield_ni', 'micro_collected_e',
    'rkf_transport_ok', 'physical_metal', 'rkf_electrons_estimate',
    'rkf_ions_estimate', 'rkf_mean_gain', 'rkf_loss',
    'rkf_gain_v1', 'rkf_ions_minus_micro_births', 'rkf_drift_seconds'
)


def write_json(p: Path, data):
    p.write_text(json.dumps(data, sort_keys=True, indent=2, allow_nan=False) + '\n')


def identity(r):
    return r['geometry'], int(r['box_id']), int(r['sample'])


def verified_v1():
    rows, plan = verify_v1_plan()  # SHA of V1 source, D500 selection, A100 manifest
    if not V1_GAS.is_file() or not V1_OUTPUT.is_file():
        raise FileNotFoundError('V1 gas/physics output not present. Never regenerate or replace it.')
    receipt_file = V1_ROOT / 'rkf_run_receipt.json'
    summary_file = V1_ROOT / 'summary.json'
    seal_file = V1_ROOT / 'report_manifest.json'
    if not all(p.is_file() for p in (receipt_file, summary_file, seal_file)):
        raise FileNotFoundError('V1 physics receipt/summary/manifest missing')
    rec = json.loads(receipt_file.read_text())
    summary = json.loads(summary_file.read_text())
    seal = json.loads(seal_file.read_text())
    expected = {
        'pilot_input_sha256': sha(V1_INPUT),
        'source_cpp_sha256': sha(V1_SOURCE),
        'Magboltz_table_sha256': sha(V1_GAS),
        'rkf_file_sha256': sha(V1_OUTPUT),
    }
    if any(rec.get(k) != v for k, v in expected.items()):
        raise ValueError('Original V1 RKF receipt SHA/provenance mismatch')
    if (rec.get('status') != 'complete' or summary.get('status') != 'complete'
        or summary.get('n_exact_pairs') != 100 or
        summary.get('rkf_output_sha256') != sha(V1_OUTPUT) or
        summary.get('gas_table_sha256') != sha(V1_GAS)):
        raise ValueError('V1 original run/summary not complete or SHA altered')
    for name, digest in seal.get('files_sha256', {}).items():
        if sha(V1_ROOT / name) != digest:
            raise ValueError('V1 comparison sealed file changed: ' + name)
    if len(seal.get('files_sha256', {})) != 5:
        raise ValueError('V1 archived file seal incomplete')
    v1 = read_csv(V1_OUTPUT, V1_COLUMNS)
    if len(v1) != 100 or {identity(x) for x in v1} != {identity(x) for x in rows}:
        raise ValueError('V1 source/result exact-start identities mismatch')
    return rows, plan, {identity(x): x for x in v1}, rec


def check_native():
    rows, plan, v1, rec = verified_v1()
    if not CXX.is_file() or not EXE.is_file():
        raise FileNotFoundError('Build new target mwpc_phase_c_rkf_ion_census_pilot_v2 first')
    result = subprocess.run([str(EXE), '--check', str(V1_INPUT)],
                            capture_output=True, text=True, check=False)
    if result.returncode or 'PASS C++ ORIGINAL D500 exact-start parser (100 rows)' not in result.stdout:
        raise RuntimeError('V2 compiled binary parser FAIL:\n' + result.stdout + '\n' + result.stderr)
    print(result.stdout.strip(), flush=True)
    return rows, plan, v1, rec


def run():
    rows, plan, v1, rec = check_native()
    if ROOT.exists():
        raise FileExistsError('V2 output exists: no overwrites allowed: ' + str(ROOT))
    ROOT.mkdir(parents=True)
    # SAME random proxy parameters as frozen V1; only new deterministic RKF
    # estimated electron and ion counts are added in V2 C++.
    theta = rec['theta_proxy']
    p0 = rec['zero_gate_proxy']
    subprocess.run([str(EXE), '--run', str(V1_INPUT), str(V1_GAS),
                    str(OUTPUT), str(theta), str(p0)], check=True)
    result = read_csv(OUTPUT, COLUMNS)
    if len(result) != 100 or len({identity(x) for x in result}) != 100 or \
            {identity(x) for x in result} != {identity(x) for x in rows}:
        raise ValueError('V2 output row count/physical starting coordinate IDs wrong')
    for x in result:
        old = v1[identity(x)]
        for field in ('source_u_mm', 'source_v_mm', 'source_w_mm',
                      'physical_metal', 'rkf_transport_ok', 'rkf_end_status'):
            if x[field] != old[field]:
                raise ValueError(f'V2 vs V1 physical coordinate/transport changed at {identity(x)} field {field}')
        for field in ('rkf_mean_gain', 'rkf_loss'):
            if not math.isclose(float(x[field]), float(old[field]), rel_tol=1.e-8, abs_tol=1.e-9):
                raise ValueError(f'V2 vs V1 RKF transport/gain not reproduced at {identity(x)} field {field}')
        if x['rkf_transport_ok'] == '1':
            ne = float(x['rkf_electrons_estimate'])
            ni = float(x['rkf_ions_estimate'])
            if not (math.isfinite(ne) and math.isfinite(ni) and ne >= 0 and ni >= 0):
                raise ValueError('V2 did not return valid estimated particle counts')
        elif x['physical_metal'] != '1':
            raise ValueError('Unexpected V2 transport failure not observed in V1')
    write_json(RECEIPT, {
        'schema': 'MWPC_PHASEC_RKF_ION_CENSUS_V2_RUN_V1',
        'status': 'complete', 'n_exact_starts': len(result),
        'v1_plan_sha256': sha(V1_PLAN),
        'v1_gas_sha256': sha(V1_GAS),
        'v1_output_sha256': sha(V1_OUTPUT),
        'v2_cpp_sha256': sha(CXX), 'v2_executable_sha256': sha(EXE),
        'v2_output_sha256': sha(OUTPUT),
        'original_A100_manifest_sha256': plan['original_A100_manifest_sha256'],
        'interpretation': ('GetAvalancheSize reports deterministic/mean RKF estimates. '
                           'No original microscopic Garfield avalanche recomputed; '
                           'physical random micro outcome is not expected to equal RKF expectation.')
    })
    print('PASS 100/100 exact original starts; V1 coordinate/gain reproduction;'
          ' V2 RKF estimated ne/ni saved; no original physics changed', flush=True)


def correlated_pearson(xs, ys) -> Optional[float]:
    if len(xs) < 3: return None
    mx, my = statistics.mean(xs), statistics.mean(ys)
    xx = sum((v - mx)**2 for v in xs)
    yy = sum((v - my)**2 for v in ys)
    return sum((x - mx)*(y - my) for x, y in zip(xs, ys))/math.sqrt(xx*yy) if xx > 0 and yy > 0 else None


def stats(items):
    good = [r for r in items if r['rkf_transport_ok'] == '1']
    metals = [r for r in items if r['physical_metal'] == '1']
    failures = [r for r in items if r['rkf_transport_ok'] != '1' and r['physical_metal'] != '1']
    micro = [float(r['micro_ion_births']) for r in items]
    ion_est = [float(r['rkf_ions_estimate']) for r in good]
    micro_good = [float(r['micro_ion_births']) for r in good]
    residuals = [float(r['rkf_ions_minus_micro_births']) for r in good]
    return {
        'n': len(items), 'rkf_transport_ok': len(good), 'in_metal': len(metals),
        'nonmetal_failures': len(failures),
        'mean_microscopic_ion_births': statistics.mean(micro),
        'mean_RKF_estimated_ion_population_good': statistics.mean(ion_est) if ion_est else None,
        'mean_RKF_estimated_electron_population_good': statistics.mean(
            float(r['rkf_electrons_estimate']) for r in good) if good else None,
        'mean_RKF_estimate_to_microscopic_ratio_pilot_only':
            statistics.mean(ion_est)/statistics.mean(micro_good)
            if ion_est and statistics.mean(micro_good)>0 else None,
        'mean_signed_RKF_estimate_minus_microscopic': statistics.mean(residuals) if residuals else None,
        'MAE_single_microscopic_realization_not_model_accuracy':
            statistics.mean(abs(v) for v in residuals) if residuals else None,
        'pearson_single_micro_realizations_not_model_calibration':
            correlated_pearson(ion_est, micro_good) if good else None,
        'n_RKF_estimated_ion_population_zero_among_good':sum(v == 0 for v in ion_est),
        'n_microscopic_zero_ion_births':sum(v == 0 for v in micro),
        'median_RKF_seconds_drift':statistics.median(float(r['rkf_drift_seconds']) for r in items),
    }


def compare():
    rows, plan, v1, rec = verified_v1()
    if not RECEIPT.is_file():
        raise FileNotFoundError('First run V2 with new native executable; no run receipt')
    receipt = json.loads(RECEIPT.read_text())
    expected = {
        'v1_plan_sha256': sha(V1_PLAN),
        'v1_gas_sha256': sha(V1_GAS),
        'v1_output_sha256': sha(V1_OUTPUT),
        'v2_cpp_sha256': sha(CXX),
        'v2_executable_sha256': sha(EXE),
        'v2_output_sha256': sha(OUTPUT),
        'original_A100_manifest_sha256': plan['original_A100_manifest_sha256'],
    }
    if receipt.get('status') != 'complete' or any(receipt.get(k) != v for k, v in expected.items()):
        raise ValueError('V2 receipt/source/provenance SHA mismatch')
    if PER_SEED.exists() or SUMMARY.exists() or MANIFEST.exists():
        raise FileExistsError('Immutable V2 comparison already exists; no rewriting')
    observed = read_csv(OUTPUT, COLUMNS)
    if {identity(x) for x in observed} != {identity(x) for x in rows}:
        raise ValueError('V2 matched original starts changed')
    a100_path = A100 / 'manifest.json'
    if sha(a100_path) != receipt['original_A100_manifest_sha256']:
        raise ValueError('Frozen A100 source manifest changed')
    a100 = json.loads(a100_path.read_text())
    if a100.get('status') != 'complete' or a100.get('new_full_Garfield_avalanches') != 16000:
        raise ValueError('Original A100 physics incomplete')
    wanted = {identity(x) for x in rows}
    truth = {}
    for geom in GEOMS:
        source = A100 / f'trials_extra_{geom}.csv'
        if sha(source) != a100['files_sha256'][source.name]:
            raise ValueError('Original A100 trial source SHA mismatch: '+source.name)
        with source.open(newline='') as f:
            for r in csv.DictReader(f):
                key = identity(r)
                if key in wanted:
                    if key in truth: raise ValueError('Duplicate microscopic physics at '+str(key))
                    if r['configuration'] != geom:raise ValueError('Microscopic geometry mismatch')
                    truth[key] = r
    if set(truth) != wanted:
        raise ValueError('Not all 100 original starts found in frozen A100 physics')
    matched = []
    for row in observed:
        k = identity(row)
        phys = truth[k]
        prev = v1[k]
        if any(not math.isclose(float(row[f]),float(prev[f]),rel_tol=1.e-8,abs_tol=1.e-9)
               for f in ('rkf_mean_gain','rkf_loss')):
            raise ValueError('V1/v2 drift field changed: '+str(k))
        if row['physical_metal'] != '0' or row['rkf_transport_ok'] != '1':
            raise ValueError('Unexpected V2 metal/failure in previously 100/100 successful pilot')
        ne = float(row['rkf_electrons_estimate'])
        ni = float(row['rkf_ions_estimate'])
        n_micro = int(phys['ion_births'])
        matched.append({
            'geometry': k[0], 'box_id': k[1], 'sample': k[2],
            'start_u_mm': row['source_u_mm'], 'start_v_mm': row['source_v_mm'],
            'micro_ion_births': n_micro, 'micro_garfield_ni':int(phys['garfield_ni']),
            'micro_collected_e': int(phys['collected_e']),
            'rkf_transport_ok': row['rkf_transport_ok'],
            'physical_metal': row['physical_metal'],
            'rkf_electrons_estimate':ne, 'rkf_ions_estimate':ni,
            'rkf_mean_gain':float(row['rkf_mean_gain']),
            'rkf_loss':float(row['rkf_loss']),
            'rkf_gain_v1':float(prev['rkf_mean_gain']),
            'rkf_ions_minus_micro_births':ni-n_micro,
            'rkf_drift_seconds':float(row['seconds_drift']),
        })
    if len(matched) != 100:
        raise ValueError('Wrong number of matched outcomes')
    summary = {
        'schema': 'MWPC_PHASEC_RKF_ION_CENSUS_MATCHED100_V2', 'status': 'complete',
        'n_exact_d500_pairs': len(matched),
        'by_geometry': {g:stats([x for x in matched if x['geometry']==g]) for g in GEOMS},
        'overall':stats(matched),
        'v1_summary_sha256':sha(V1_ROOT/'summary.json'),
        'v2_output_sha256':sha(OUTPUT),
        'source_A100_manifest_sha256':sha(a100_path),
        'interpretation':('RKF ne/ni are macroscopic estimates of electron/ion population, '
             'not actual microscopic collision histories. One real Garfield realization '
             'per exact position is not an expectation and has large heavy-tailed variance. '
             'This 100-start source-selected pilot cannot calibrate event-by-event physics '
             'or establish performance on independent full muons. No charge induction.'),
    }
    with PER_SEED.open('w',newline='') as f:
        writer = csv.DictWriter(f, fieldnames=COMPARE_FIELDS)
        writer.writeheader();writer.writerows(matched)
    write_json(SUMMARY, summary)
    write_json(MANIFEST, {
        'status':'complete', 'schema':'MWPC_PHASEC_RKF_ION_CENSUS_MATCHED100_V2_SEAL',
        'sha256':{x.name:sha(x) for x in (RECEIPT, PER_SEED, SUMMARY)},
        'preserved_v1_sha256':{
            'gas':sha(V1_GAS),'pilot_input':sha(V1_INPUT),
            'source_cpp':sha(V1_SOURCE),'pilot_output':sha(V1_OUTPUT)},
    })
    print(json.dumps(summary, indent=2, ensure_ascii=False), flush=True)
    print('PASS V2 matched 100 actual microscopic Garfield starts vs RKF ion estimates; V1 unchanged', flush=True)


def selftest():
    a = [
        {'geometry':'p3_noFW','box_id':1,'sample':20,'rkf_transport_ok':'1','physical_metal':'0',
         'micro_ion_births':10,'rkf_ions_estimate':12,'rkf_electrons_estimate':13,'rkf_ions_minus_micro_births':2,'rkf_drift_seconds':.001},
        {'geometry':'p3_noFW','box_id':2,'sample':21,'rkf_transport_ok':'1','physical_metal':'0',
         'micro_ion_births':0,'rkf_ions_estimate':5,'rkf_electrons_estimate':6,'rkf_ions_minus_micro_births':5,'rkf_drift_seconds':.002},
    ]
    s = stats(a)
    assert s['n']==2 and s['rkf_transport_ok']==2
    assert s['mean_microscopic_ion_births']==5 and s['mean_RKF_estimated_ion_population_good']==8.5
    assert s['n_microscopic_zero_ion_births']==1 and s['nonmetal_failures']==0
    assert identity(a[0])==('p3_noFW',1,20)
    assert COLUMNS[10:12] == ('rkf_electrons_estimate','rkf_ions_estimate')
    print('PASS V2 read-only statistical selftest; no Garfield physics',flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('mode', choices=['selftest','buildcheck','run','compare','status'])
    a = p.parse_args()
    if a.mode == 'selftest': selftest()
    elif a.mode == 'buildcheck': check_native()
    elif a.mode == 'run': run()
    elif a.mode == 'compare':compare()
    else:
        print('V2_ROOT',ROOT, 'EXISTS',ROOT.exists(),
              'OUTPUT',OUTPUT.exists(),'RECEIPT',RECEIPT.exists(),
              'SUMMARY',SUMMARY.exists(),'MANIFEST',MANIFEST.exists())
        print('V1_ROOT',V1_ROOT,'V1_PLAN',V1_PLAN.exists(),
              'V1_GAS',V1_GAS.exists(),'V1_OUTPUT',V1_OUTPUT.exists())

if __name__ == '__main__': main()