#!/usr/bin/env python3
"""Read-only independent audit of completed and SHA-sealed RKF V3 16k results.

Repairs the V3 *reporting* type bug without overwriting or re-running any
Garfield++/Magboltz/Geant4/TrackHeed calculations. Recompute validated status,
ion means, cell correlations and a conservative paired-cell bootstrap in NEW
sibling output. Old incorrect source report and its SHA seal remain intact.
"""
from __future__ import annotations
import argparse
import csv
import json
import math
import statistics
from collections import Counter, defaultdict
from pathlib import Path

import pi_run_rkf_ion_census_full16k_v3 as v3

OUT = v3.BASE / 'pi_RKF_ion_census_B0_A100all16k_v3_status_audit_v1'

def require(condition, reason):
    if not condition:
        raise ValueError('V3 audit: '+reason)

def close(a, b, field, tol=1e-8):
    require(math.isclose(float(a), float(b), rel_tol=tol, abs_tol=tol),
            'mismatch in '+field+f': {a} != {b}')

def check_archive_seal():
    v3.verify_source()  # all old 16k A100 full Garfield, gas, V1 & V2 SHA
    for p in (v3.RECEIPT, v3.OUTPUT, v3.PAIRS, v3.CELLS, v3.SUMMARY, v3.SEAL):
        require(p.is_file(), 'missing immutable V3 original: '+str(p))
    r=json.loads(v3.RECEIPT.read_text())
    require(r.get('status')=='complete' and r.get('n')==16000, 'V3 run receipt incomplete')
    for name,path in [('input_sha256',v3.INPUT),('gas_sha256',v3.GAS),
                      ('v3_cpp_sha256',v3.CXX),('v3_executable_sha256',v3.EXE),
                      ('v3_output_sha256',v3.OUTPUT)]:
        require(r.get(name)==v3.sha(path), 'original V3 receipt SHA mismatch: '+name)
    m=json.loads(v3.SEAL.read_text())
    require(m.get('status')=='complete' and
            set(m.get('sha256',{}))=={p.name for p in (v3.RECEIPT,v3.PAIRS,v3.CELLS,v3.SUMMARY)},
            'incomplete original comparison seal')
    for name,h in m['sha256'].items():
        require(v3.sha(v3.ROOT/name)==h,'original comparison SHA mismatch: '+name)
    return r

def row_id(row):
    return row['geometry'], int(row['box_id']), int(row['sample'])

def analyze():
    if OUT.exists():
        raise FileExistsError('Immutable audit output already exists (will not overwrite): '+str(OUT))
    receipt=check_archive_seal()
    raw=v3.read_csv(v3.OUTPUT,v3.RKFCOLS)
    pairs=v3.read_csv(v3.PAIRS,v3.PAIRCOLS)
    cells=v3.read_csv(v3.CELLS,v3.CELLCOLS)
    original=json.loads(v3.SUMMARY.read_text())
    require(len(raw)==len(pairs)==16000 and len(cells)==200, 'original row count')
    lookup={row_id(x):x for x in raw}
    require(len(lookup)==16000, 'duplicate RKF raw starts')
    micro_lookup={}
    for geom in v3.GEOMS:
        for x in v3.read_csv(v3.A100/f'trials_extra_{geom}.csv'):
            require(x['configuration']==geom and x['original_geometry']==geom,
                    'original microscopic geometry mismatch')
            k=(geom,int(x['box_id']),int(x['sample']))
            require(k not in micro_lookup,'duplicate microscopic source')
            micro_lookup[k]=int(x['ion_births'])
    require(len(micro_lookup)==16000, 'original A100 trials incomplete')
    groups=defaultdict(list)
    rows=[]
    for p in pairs:
        k=row_id(p)
        require(k in lookup and k in micro_lookup, 'missing raw or microscopic start: '+repr(k))
        r=lookup.pop(k)
        require(p['rkf_transport_ok'] in ('0','1') and p['rkf_metal'] in ('0','1'),
                'invalid normalized RKF status flags')
        require(r['rkf_transport_ok']==p['rkf_transport_ok'] and
                r['physical_metal']==p['rkf_metal'], 'status differs raw vs matched: '+str(k))
        require(int(p['micro_ion_births'])==micro_lookup[k], 'micro ion census differs source')
        require(int(p['micro_zero'])==(micro_lookup[k]==0),'micro zero indicator changed')
        for a,b in (('rkf_ni','rkf_ions_estimate'),('rkf_ne','rkf_electrons_estimate'),
                    ('rkf_gain','rkf_mean_gain'),('rkf_loss','rkf_loss'),
                    ('rkf_wall_seconds','seconds_drift')):
            if p[a] or r[b]: close(p[a],r[b],a)
        if p['rkf_transport_ok']=='1':
            ni=float(p['rkf_ni'])
            require(math.isfinite(ni) and ni>=0, 'invalid successful RKF ion estimate')
            require(p['rkf_metal']=='0','metal transport unexpectedly succeeded')
        rows.append(p);groups[(k[0],k[1])].append(p)
    require(not lookup, 'unpaired RKF raw start(s)')
    require(len(groups)==200 and all(len(x)==80 for x in groups.values()),
            'incorrect selected cell group sizes')
    recorded_cells={(x['geometry'],int(x['box_id'])):x for x in cells}
    require(set(groups)==set(recorded_cells),'cell ids mismatch')
    for k,p in groups.items():
        source=recorded_cells[k]
        m=statistics.mean(int(x['micro_ion_births']) for x in p)
        good=[x for x in p if x['rkf_transport_ok']=='1']
        require(len(good)==80 and not any(x['rkf_metal']=='1' for x in p),
                'a failed/metal RKF start cannot be silently used in mean: '+str(k))
        ni=statistics.mean(float(x['rkf_ni']) for x in good)
        close(source['micro_mean_ions'],m,'cell micro mean '+str(k))
        close(source['rkf_mean_ions'],ni,'cell RKF mean '+str(k))
    by_geom={}
    for g in v3.GEOMS:
        rr=[x for x in rows if x['geometry']==g]
        cc=[dict(box_id=int(x['box_id']),micro_mean_ions=float(x['micro_mean_ions']),
                 rkf_mean_ions=float(x['rkf_mean_ions'])) for x in cells if x['geometry']==g]
        require(len(rr)==4000 and len(cc)==50,'wrong geometry size '+g)
        st=v3.summary_for(rr) # fixed integer-vs-string normalization
        st['n_cells']=50
        st['cell_mean_pearson']=v3.pearson([x['rkf_mean_ions'] for x in cc],
                                           [x['micro_mean_ions'] for x in cc])
        st['ratio_rkf_over_micro_cluster_bootstrap95']=v3.bootstrap_ratio(cc)
        st['crossfit_25cells_x2']=v3.crossfit_scale(cc)
        # These remain scientifically interpretable despite former flag bug.
        old=original['by_geometry'][g]
        close(st['cell_mean_pearson'],old['cell_mean_pearson'],'original cell Pearson '+g)
        for name in ('mean_heldout_abs_cell_error_unscaled',
                     'mean_heldout_abs_cell_error_scaled'):
            close(st['crossfit_25cells_x2'][name],old['crossfit_25cells_x2'][name],
                  'original cell crossfit '+g+' '+name)
        by_geom[g]=st
    overall=v3.summary_for(rows)
    require(overall['n']==16000 and overall['n_transport_ok']==16000 and
            overall['n_nonmetal_fail']==0 and overall['n_in_metal']==0,
            'physical paired status count contradicts 16k successful cells')
    out={'schema':'MWPC_RKF_V3_16K_ORIGINAL_SHA_AUDIT_CORRECTED_STATUS_V1',
         'status':'complete','n_original_pairs':16000,
         'original_raw_v3_sha256':v3.sha(v3.OUTPUT),
         'original_paired_v3_sha256':v3.sha(v3.PAIRS),
         'original_incorrect_summary_sha256':v3.sha(v3.SUMMARY),
         'original_v3_seal_sha256':v3.sha(v3.SEAL),
         'root_cause':('Original V3 summary_for compared integer paired rkf_transport_ok/rkf_metal '
                       'to string flags, reporting all 16000 events failed even though the '
                       'cell-level code correctly selected 80 successful starts per cell.'),
         'method':('All 16000 raw RKF and matched Garfield rows independently cross-joined '
                    'and SHA verified; all 200 x 80 cell averages verified against archived cells. '
                    'Original source summary left unchanged and sealed.'),
         'corrected_overall':overall, 'corrected_by_geometry':by_geom,
         'original_bad_status':{'overall_n_transport_ok':original['overall']['n_transport_ok'],
                                'overall_n_nonmetal_fail':original['overall']['n_nonmetal_fail']},
         'caveats':('100% successful RKF transport is numerical, not avalanche survival. '
                   'RKF estimated ni is a macroscopic mean, not a stochastic microscopic event. '
                   'Cells were A20/120-source selected. Twofold cell calibration is internal, '
                   'not unseen Geant4 muon validation. N=80 separate starts/cell, '
                   'not N=80 repeated trials at identical position. No induced strip signals.')}
    OUT.mkdir(parents=True)
    v3.write_json(OUT/'corrected_summary.json',out)
    report=['# RKF V3 corrected status audit (read-only)', '',
        'Original 16,000-row RKF and Garfield datasets and original comparison SHA seal were verified. '
        'The original report remains intact but had a status integer/string comparison bug.', '',
        '| Geometry | RKF transport | Micro mean ions | RKF mean estimated ions | Ratio | Cell Pearson | Cell MAE before → after calibration |',
        '|---|---:|---:|---:|---:|---:|---:|']
    for g,s in by_geom.items():
        x=s['crossfit_25cells_x2']
        report.append(f"| {g} | {s['n_transport_ok']}/{s['n']} | "
                      f"{s['mean_micro_ions']:.2f} | {s['mean_rkf_ions_on_transport_ok']:.2f} | "
                      f"{s['rkf_micro_ratio_on_good']:.3f} | {s['cell_mean_pearson']:.3f} | "
                      f"{x['mean_heldout_abs_cell_error_unscaled']:.1f} → "
                      f"{x['mean_heldout_abs_cell_error_scaled']:.1f} |")
    report+=['',f"All 16,000: RKF successes {overall['n_transport_ok']}, "
             f"nonmetal failures {overall['n_nonmetal_fail']}, "
             f"original micro zeros {overall['n_micro_zero']}.",
             '', '**Do not use the original V3 summary\'s 0/16,000 transport statistic.**',
             'Use the tabulated cell correlations and calibration errors to assess spatial predictivity; '
             'not a validated stochastic replacement of microscopic Garfield.',
             'No new microscopic or macroscopic physics was run.']
    (OUT/'audit_note.md').write_text('\n'.join(report)+'\n')
    v3.write_json(OUT/'audit_manifest.json',{
       'schema':'MWPC_RKF_V3_STATUS_AUDIT_SEAL_V1','status':'complete',
       'sha256':{p.name:v3.sha(p) for p in [OUT/'corrected_summary.json',OUT/'audit_note.md']},
       'preserved_original_SHA':{'original_v3_seal':v3.sha(v3.SEAL),
                                 'original_v3_output':v3.sha(v3.OUTPUT)}})
    print('\n'.join(report),flush=True)
    print('PASS SHA-VERIFIED V3 READ-ONLY STATUS AUDIT; 16k raw/pairs and all 200 cells; original output untouched',flush=True)

def selftest():
    good={'geometry':'p3_noFW','rkf_transport_ok':1,'rkf_metal':0,
          'micro_ion_births':200,'micro_zero':0,'rkf_ni':'280',
          'rkf_wall_seconds':'.001'}
    same=dict(good,rkf_transport_ok='1',rkf_metal='0')
    for xs in ([good],[same],[good,same]):
        t=v3.summary_for(xs)
        require(t['n_transport_ok']==len(xs) and t['n_nonmetal_fail']==0,
                'regression: integer and string flags')
    bad=dict(good,rkf_transport_ok=0,rkf_metal=0)
    s=v3.summary_for([good,bad])
    require(s['n_transport_ok']==1 and s['n_nonmetal_fail']==1,
            'regression: nonmetal failure classification')
    print('PASS V3 audit synthetic flag/string/integer selftest; no physics',flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('mode',choices=['selftest','audit'])
    a=p.parse_args()
    selftest() if a.mode=='selftest' else analyze()