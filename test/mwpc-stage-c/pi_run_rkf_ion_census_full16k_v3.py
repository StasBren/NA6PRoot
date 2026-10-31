#!/usr/bin/env python3
"""Separate, SHA-verified all-16k A100 original D500-start RKF ion census.

Runs the user's *already compiled* Garfield++ macroscopic DriftLineRKF on
16,000 ORIGINAL starting positions from the 200 A100-selected cells, with
one existing and immutable microscopic Garfield avalanche per start.
This does not run any new microscopic avalanches, Geant4, or TrackHeed.

The output is a diagnostic of group-mean ion response and its spatial
variation, NOT a proof of microscopic avalanche reproducibility, and the
A100 target cells are position/exposure-selected on the original 120 muons.

V3 inherits a Pólya proxy column from V2 C++; it is deliberately excluded
from every physics statistic, since the proxy is not physically calibrated.

Modes: selftest | buildcheck | run | compare | status
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
import random
import statistics
import subprocess
from collections import defaultdict
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
BASE = REPO / 'test_runs' / 'mwpc_stage_c'
A100 = BASE / 'pi_A100_targeted_50cells_per_geometry_v1'
V1 = BASE / 'pi_RKF_Polya_separate_B0_pilot_v1'
V2 = BASE / 'pi_RKF_ion_census_B0_pilot_v2'
ROOT = BASE / 'pi_RKF_ion_census_B0_A100all16k_v3'
CXX = REPO / 'test/mwpc-avalanche/phaseC_rkf_ion_census_full16k_v3.cxx'
EXE = REPO / 'test/mwpc-avalanche/build/mwpc_phase_c_rkf_ion_census_full16k_v3'
INPUT = A100 / 'selected_exact_old_D500_starts.csv'
GAS = V1 / 'Ar70_CO2_30_293K_760T_B0_Magboltz40E.gas'
OUTPUT = ROOT / 'rkf_estimated_ions_16k_exact_starts.csv'
RECEIPT = ROOT / 'run_receipt.json'
PAIRS = ROOT / 'matched_16k_micro_vs_rkf_ions.csv'
CELLS = ROOT / 'cell_80_micro_vs_80_rkf.csv'
SUMMARY = ROOT / 'summary.json'
SEAL = ROOT / 'comparison_manifest.json'
GEOMS = ('p3_noFW', 'p3_FW75', 'p4_noFW', 'p4_FW75')
SHA_A100 = '2a8520c61dc38dc2b32aa31319c3665e36d9ed86f9119e8f3737f32849d7da82'
SHA_V1_GAS = 'f6f6a8a561d52c788c7fe8168256d3add78bdca131172fcf9526159707e1e693'
SHA_V1_OUTPUT = 'a6464df3a5ba52b5077b1994f4963477fa186f583335d5836b56a04a866f45a4'
SHA_V2_OUTPUT = 'e5dd149a2d5f02d0fd2e07f9a4ced9ed386eee0442bd374b6706995a81d9cee7'
TRIALCOLS = ('lut_trials','geometry','original_lut_sha256','box_id','iu','iv','sample',
 'random_seed','u_start_mm','v_start_mm','w_start_mm','old_drift_ok','old_collected',
 'old_field_wire_collected','old_wire_index','old_phi_hit_deg','old_drift_time_ns',
 'old_last_step_energy_eV')
RKFCOLS = ('geometry','box_id','sample','source_u_mm','source_v_mm','source_w_mm',
 'physical_metal','rkf_transport_ok','rkf_mean_gain','rkf_loss',
 'rkf_electrons_estimate','rkf_ions_estimate','rkf_final_u_mm',
 'rkf_final_v_mm','rkf_final_w_mm','rkf_end_status',
 'polya_proxy_ion_births','polya_theta','p0_proxy','seconds_drift')
PAIRCOLS = ('geometry','box_id','sample','u_start_mm','v_start_mm',
 'micro_ion_births','micro_ni','micro_electrons','micro_zero','rkf_ne','rkf_ni',
 'rkf_gain','rkf_loss','rkf_metal','rkf_transport_ok',
 'rkf_minus_micro','rkf_wall_seconds')
CELLCOLS = ('geometry','box_id','n','micro_mean_ions','rkf_mean_ions',
 'micro_zero_fraction','rkf_to_micro_mean_ratio', 'mean_signed_error',
 'rkf_MAE_vs_single_micro_diagnostics_only','rkf_ions_sd_across_starts',
 'micro_ions_sd_across_starts','n_metal','n_nonmetal_failed')


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):
            h.update(b)
    return h.hexdigest()


def require_sha(path,expected):
    if not path.is_file() or sha(path)!=expected:
        raise ValueError('Missing or wrong SHA for immutable source: '+str(path))


def read_csv(path,fields=None):
    with Path(path).open(newline='') as f:
        reader=csv.DictReader(f)
        if fields is not None and tuple(reader.fieldnames or ()) != fields:
            raise ValueError(f'CSV schema mismatch: {path}: {reader.fieldnames}')
        return list(reader)


def write_csv(path,fields,rows):
    with path.open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=fields)
        writer.writeheader();writer.writerows(rows)


def write_json(path,obj):
    path.write_text(json.dumps(obj,sort_keys=True,indent=2,allow_nan=False)+'\n')


def key(x,field='geometry'):
    return x[field],int(x['box_id']),int(x['sample'])


def pearson(a,b):
    if len(a)<3:return None
    ma=statistics.mean(a);mb=statistics.mean(b)
    xx=sum((v-ma)**2 for v in a);yy=sum((v-mb)**2 for v in b)
    return sum((x-ma)*(y-mb) for x,y in zip(a,b))/math.sqrt(xx*yy) if xx>0 and yy>0 else None


def summary_for(items):
    good=[a for a in items if a['rkf_transport_ok']=='1']
    metals=[a for a in items if a['rkf_metal']=='1']
    failed=[a for a in items if a['rkf_transport_ok']!='1' and a['rkf_metal']!='1']
    micro=[float(a['micro_ion_births']) for a in items]
    rkf=[float(a['rkf_ni']) for a in good]
    micro_good=[float(a['micro_ion_births']) for a in good]
    return {'n':len(items),'n_transport_ok':len(good),'n_in_metal':len(metals),
        'n_nonmetal_fail':len(failed), 'n_micro_zero':sum(int(a['micro_zero']) for a in items),
        'mean_micro_ions':statistics.mean(micro),
        'mean_rkf_ions_on_transport_ok':statistics.mean(rkf) if rkf else None,
        'rkf_micro_ratio_on_good':(statistics.mean(rkf)/statistics.mean(micro_good))
             if rkf and statistics.mean(micro_good)>0 else None,
        'pearson_single_draw_uninformative_for_conditional_mean':pearson(rkf,micro_good),
        'rkf_ion_zero_on_good':sum(v==0 for v in rkf),
        'median_drift_seconds':statistics.median(float(a['rkf_wall_seconds']) for a in items)}


def verify_source():
    require_sha(A100/'manifest.json',SHA_A100)
    man=json.loads((A100/'manifest.json').read_text())
    if man.get('status')!='complete' or man.get('new_full_Garfield_avalanches')!=16000:
        raise ValueError('A100 source not complete')
    require_sha(INPUT,man['files_sha256'][INPUT.name])
    require_sha(GAS,SHA_V1_GAS)
    require_sha(V1/'rkf_polya_proxy_100_starts.csv',SHA_V1_OUTPUT)
    require_sha(V2/'rkf_estimated_ions_100_exact_starts.csv',SHA_V2_OUTPUT)
    v1plan=json.loads((V1/'pilot_manifest.json').read_text())
    if v1plan['original_A100_manifest_sha256']!=SHA_A100:
        raise ValueError('V1 pilot source identity changed')
    v2rec=json.loads((V2/'run_receipt.json').read_text())
    if (v2rec.get('status')!='complete' or v2rec.get('v1_gas_sha256')!=SHA_V1_GAS
        or v2rec.get('v2_output_sha256')!=SHA_V2_OUTPUT):
        raise ValueError('V2 run receipt changed')
    starts=read_csv(INPUT,TRIALCOLS)
    expected=16000
    if len(starts)!=expected or len({key(x) for x in starts})!=expected:
        raise ValueError('A100 exact start 16k incomplete or duplicated')
    for g in GEOMS:
        group=[x for x in starts if x['geometry']==g]
        if len(group)!=4000 or len({x['box_id'] for x in group})!=50:
            raise ValueError('A100 selected cells not 50x80 for '+g)
        counts=defaultdict(int)
        for x in group:counts[int(x['box_id'])]+=1
        if any(v!=80 for v in counts.values()):raise ValueError('Not 80 starts per cell')
        p=A100/f'trials_extra_{g}.csv'
        require_sha(p,man['files_sha256'][p.name])
    return starts,man


def check_native():
    starts,man=verify_source()
    if not CXX.is_file() or not EXE.is_file():
        raise FileNotFoundError('Compile separate RKF full16k V3 target first')
    proc=subprocess.run([str(EXE),'--check',str(INPUT)],text=True,capture_output=True)
    if proc.returncode or 'PASS C++ ORIGINAL D500 exact-start parser (16000 rows)' not in proc.stdout:
        raise RuntimeError('V3 native CSV preflight failed:\n'+proc.stdout+'\n'+proc.stderr)
    print(proc.stdout.strip(),flush=True)
    return starts,man


def run():
    starts,man=check_native()
    if ROOT.exists():raise FileExistsError('Immutable V3 output directory already exists: '+str(ROOT))
    ROOT.mkdir(parents=True)
    # Proxy is inherited from V2 and excluded from all V3 scientific conclusions.
    cmd=[str(EXE),'--run',str(INPUT),str(GAS),str(OUTPUT),'0.4','0.17']
    subprocess.run(cmd,check=True)
    out=read_csv(OUTPUT,RKFCOLS)
    if len(out)!=16000 or {key(x) for x in out}!={key(x) for x in starts}:
        raise ValueError('V3 output incomplete or mismatched exact start identifiers')
    start_idx={key(x):x for x in starts}
    for r in out:
        x=start_idx[key(r)]
        for col,outcol in [('u_start_mm','source_u_mm'),('v_start_mm','source_v_mm'),('w_start_mm','source_w_mm')]:
            if not math.isclose(float(x[col]),float(r[outcol]),abs_tol=1e-8):
                raise ValueError('V3 source coordinates changed: '+str(key(r)))
        if r['rkf_transport_ok']=='1' and (not math.isfinite(float(r['rkf_ions_estimate'])) or float(r['rkf_ions_estimate'])<0):
            raise ValueError('Invalid V3 ion population: '+str(key(r)))
    write_json(RECEIPT,{'schema':'MWPC_RKF_V3_16000_RUN','status':'complete',
        'n':16000,'input_sha256':sha(INPUT),'gas_sha256':sha(GAS),
        'a100_manifest_sha256':SHA_A100,'v1_result_sha256':SHA_V1_OUTPUT,
        'v2_result_sha256':SHA_V2_OUTPUT,'v3_cpp_sha256':sha(CXX),
        'v3_executable_sha256':sha(EXE),'v3_output_sha256':sha(OUTPUT),
        'warning':'Proxy Pólya field inherited from prior pilot; NOT used for any physics comparisons.'})
    print('PASS V3 exactly 16,000 RKF macroscopic estimates, receipt SHA sealed; no microscopic physics recomputed',flush=True)


def crossfit_scale(cells):
    """Spatial cell-holdout calibration diagnostic, not genuine external generalization."""
    ordered=sorted(cells,key=lambda z:int(z['box_id']))
    shuffled=list(ordered)
    random.Random(20261031).shuffle(shuffled)
    folds=(shuffled[:25],shuffled[25:])
    if len(ordered)!=50 or any(len(x)!=25 for x in folds):raise ValueError('Crossfit expects 50 selected cells')
    res=[]
    for train,test in ((folds[0],folds[1]),(folds[1],folds[0])):
        scale=sum(x['micro_mean_ions'] for x in train)/sum(x['rkf_mean_ions'] for x in train)
        errbase=[abs(x['rkf_mean_ions']-x['micro_mean_ions']) for x in test]
        errcal=[abs(scale*x['rkf_mean_ions']-x['micro_mean_ions']) for x in test]
        res.append({'train_cell_count':25,'test_cell_count':25,
             'training_scale_micro_over_rkf':scale,
             'heldout_mean_abs_cell_error_unscaled':statistics.mean(errbase),
             'heldout_mean_abs_cell_error_scaled':statistics.mean(errcal),
             'heldout_micro_mean_ions':statistics.mean(x['micro_mean_ions'] for x in test),
             'heldout_rkf_mean_ions':statistics.mean(x['rkf_mean_ions'] for x in test)})
    return {'two_folds':res,
      'mean_heldout_abs_cell_error_unscaled':statistics.mean(x['heldout_mean_abs_cell_error_unscaled'] for x in res),
      'mean_heldout_abs_cell_error_scaled':statistics.mean(x['heldout_mean_abs_cell_error_scaled'] for x in res)}


def bootstrap_ratio(cells,reps=3000,seed=20261031):
    rnd=random.Random(seed)
    values=[]
    for _ in range(reps):
        chosen=[cells[rnd.randrange(len(cells))] for _ in cells]
        a=sum(x['micro_mean_ions'] for x in chosen)
        b=sum(x['rkf_mean_ions'] for x in chosen)
        values.append(b/a if a>0 else float('nan'))
    values=[x for x in values if math.isfinite(x)]
    if len(values)!=reps:raise ValueError('Invalid finite bootstrap sampled ratios')
    values.sort()
    return [values[round(.025*(reps-1))],values[round(.975*(reps-1))]]


def compare():
    starts,man=verify_source()
    if not RECEIPT.is_file():raise FileNotFoundError('V3 run not complete; never rerun existing output')
    receipt=json.loads(RECEIPT.read_text())
    expected={'input_sha256':sha(INPUT),'gas_sha256':sha(GAS),
        'v3_cpp_sha256':sha(CXX),'v3_executable_sha256':sha(EXE),
        'v3_output_sha256':sha(OUTPUT)}
    if receipt.get('status')!='complete' or receipt.get('n')!=16000 or any(receipt.get(k)!=v for k,v in expected.items()):
        raise ValueError('V3 receipt SHA mismatch or incomplete')
    if any(x.exists() for x in (PAIRS,CELLS,SUMMARY,SEAL)):
        raise FileExistsError('V3 comparison output already exists, preserve immutable results')
    rkf=read_csv(OUTPUT,RKFCOLS)
    if len(rkf)!=16000:raise ValueError('RKF count wrong')
    old={key(x):x for x in rkf}
    if len(old)!=16000:raise ValueError('Duplicate V3 RKF start')
    trial=[]
    for g in GEOMS:
        path=A100/f'trials_extra_{g}.csv'
        vals=read_csv(path)
        if len(vals)!=4000:raise ValueError('Corrupt physical A100 '+g)
        for r in vals:
            if r['original_geometry']!=g or r['configuration']!=g:
                raise ValueError('A100 original & micro geometry mismatch')
            k=key(r,'configuration')
            if k not in old:raise ValueError('A100 start lacks RKF: '+str(k))
            v=old[k]
            for src,dst in (('original_u_start_mm','source_u_mm'),
                            ('original_v_start_mm','source_v_mm'),
                            ('original_w_start_mm','source_w_mm')):
                if not math.isclose(float(r[src]),float(v[dst]),abs_tol=1e-8):
                    raise ValueError('Original microscopic versus RKF coordinates mismatch '+str(k))
            if int(r['successful'])!=1 or int(r['capped'])!=0:
                raise ValueError('Unexpected original microscopic failed/capped event')
            micro=int(r['ion_births'])
            if int(r['garfield_ni'])!=micro:
                raise ValueError('Original Garfield ni differs from microscopic ion birth census '+str(k))
            ok=int(v['rkf_transport_ok']);metal=int(v['physical_metal'])
            value=None if not ok else float(v['rkf_ions_estimate'])
            trial.append({'geometry':g,'box_id':int(r['box_id']),'sample':int(r['sample']),
              'u_start_mm':r['original_u_start_mm'],'v_start_mm':r['original_v_start_mm'],
              'micro_ion_births':micro,'micro_ni':int(r['garfield_ni']),
              'micro_electrons':int(r['collected_e']),'micro_zero':int(micro==0),
              'rkf_ne':v['rkf_electrons_estimate'],'rkf_ni':v['rkf_ions_estimate'],
              'rkf_gain':v['rkf_mean_gain'],'rkf_loss':v['rkf_loss'],
              'rkf_metal':metal,'rkf_transport_ok':ok,
              'rkf_minus_micro':'' if value is None else value-micro,
              'rkf_wall_seconds':v['seconds_drift']})
    if len(trial)!=16000 or len({key(x) for x in trial})!=16000:
        raise ValueError('Unexpected trial count or duplicates in V3 comparison')
    # Existing V2 frozen RKF mean/ion values MUST be exactly matched at 100
    # common original physical start positions. Pólya random proxy is ignored.
    v2=read_csv(V2/'rkf_estimated_ions_100_exact_starts.csv',RKFCOLS)
    if len(v2)!=100:raise ValueError('V2 source not 100 rows')
    for z in v2:
        found=old[key(z)]
        for f in ('source_u_mm','source_v_mm','source_w_mm',
                  'rkf_mean_gain','rkf_loss','rkf_electrons_estimate','rkf_ions_estimate'):
            if not math.isclose(float(z[f]),float(found[f]),rel_tol=1e-8,abs_tol=1e-8):
                raise ValueError('V3 RKF drift disagrees with original V2 '+str(key(z))+' '+f)
        if any(z[f]!=found[f] for f in ('physical_metal','rkf_transport_ok','rkf_end_status')):
            raise ValueError('V3 vs V2 drift outcome mismatch')
    groups=defaultdict(list)
    for x in trial:groups[(x['geometry'],int(x['box_id']))].append(x)
    if len(groups)!=200 or any(len(x)!=80 for x in groups.values()):
        raise ValueError('V3 missing original 200 spatial cell groups of 80 each')
    cells=[]
    for (geom,box),rr in sorted(groups.items()):
        m=[float(x['micro_ion_births']) for x in rr]
        g=[float(x['rkf_ni']) for x in rr if x['rkf_transport_ok']==1]
        failed=sum(x['rkf_transport_ok']==0 and x['rkf_metal']==0 for x in rr)
        if failed or len(g)!=80:
            # Retain failed records but do NOT manufacture false zero or calibration.
            raise ValueError('RKF nonmetal/metal subset prevents fair unconditional cell mean comparison: '+str((geom,box)))
        mm=statistics.mean(m);gg=statistics.mean(g)
        cells.append({'geometry':geom,'box_id':box,'n':80,
            'micro_mean_ions':mm,'rkf_mean_ions':gg,
            'micro_zero_fraction':sum(x==0 for x in m)/80,
            'rkf_to_micro_mean_ratio':gg/mm if mm>0 else '',
            'mean_signed_error':gg-mm,
            'rkf_MAE_vs_single_micro_diagnostics_only':statistics.mean(abs(x-y) for x,y in zip(g,m)),
            'rkf_ions_sd_across_starts':statistics.pstdev(g),
            'micro_ions_sd_across_starts':statistics.pstdev(m),
            'n_metal':sum(x['rkf_metal']==1 for x in rr),
            'n_nonmetal_failed':failed})
    by_geom={}
    for g in GEOMS:
        rr=[x for x in trial if x['geometry']==g]
        cc=[x for x in cells if x['geometry']==g]
        stat=summary_for(rr)
        stat.update({'n_cells':50,'cell_mean_pearson':pearson([x['rkf_mean_ions'] for x in cc],
                              [x['micro_mean_ions'] for x in cc]),
                     'ratio_rkf_over_micro_cluster_bootstrap95':bootstrap_ratio(cc),
                     'crossfit_25cells_x2':crossfit_scale(cc)})
        by_geom[g]=stat
    summary={'schema':'MWPC_RKF_V3_MICRO_A100_SELECTED_200CELLS_x80_V1',
        'status':'complete','n_original_independent_single_micro_trials':16000,
        'n_selected_cells':200,'n_same_starts_as_v2_verified':100,
        'overall':summary_for(trial),'by_geometry':by_geom,
        'source_manifest_sha256':SHA_A100,'v1_gas_sha256':SHA_V1_GAS,
        'v2_rkf_verified_sha256':SHA_V2_OUTPUT,
        'interpretation':('This is a deterministic macroscopic RKF ion estimate compared with ONE '
        'independent microscopic Garfield draw per original D500 start. 80 starts per selected '
        'spatial cell permit cell-mean comparison, and 50 cell bootstrap respects within-cell '
        'dependence. All 200 cells were selected conditional on old A20 variance/120-source '
        'exposure, so this is NOT blind generalization. Crossfit holds out spatial cells '
        'only from fitting illustrative multiplicative calibration; do not call it external '
        'model validation. No microscopic Garfield, Geant4, ion-cloud or induced signal reruns. '
        'Uncalibrated random Pólya proxy column in C++ is not used in any scientific statistic.'),
        'comparison_caveat':'If any physically in-metal or nonmetal failed RKF start is encountered, '
               'this version stops instead of constructing a biased mean; inspect such starts separately.'}
    write_csv(PAIRS,PAIRCOLS,trial)
    write_csv(CELLS,CELLCOLS,cells)
    write_json(SUMMARY,summary)
    write_json(SEAL,{'schema':'MWPC_RKF_V3_16K_SHA_SEAL_V1','status':'complete',
        'sha256':{x.name:sha(x) for x in (RECEIPT,PAIRS,CELLS,SUMMARY)},
        'source_sha256':{'a100':SHA_A100,'gas':SHA_V1_GAS,'v2_rkf':SHA_V2_OUTPUT}})
    print(json.dumps(summary,indent=2,ensure_ascii=False),flush=True)
    print('PASS V3 matched 16k Garfield microscopic starts with deterministic RKF ion estimates',flush=True)


def selftest():
    row={'geometry':'p3_noFW','box_id':'12','sample':'24'}
    assert key(row)==('p3_noFW',12,24)
    assert key({'configuration':'p3_noFW','box_id':'12','sample':'24'},'configuration')==key(row)
    assert pearson([1,2,3],[2,4,6])>0.9999
    assert pearson([1,1,1],[1,2,3]) is None
    dat=[dict(micro_mean_ions=10.,rkf_mean_ions=20.,box_id=i) for i in range(50)]
    cf=crossfit_scale(dat)
    assert all(math.isclose(x['training_scale_micro_over_rkf'],.5) for x in cf['two_folds'])
    assert cf['mean_heldout_abs_cell_error_scaled']==0
    assert bootstrap_ratio(dat,200)==[2.,2.]
    print('PASS V3 offline statistics/crossfit/identity selftest; no Garfield physics',flush=True)


def status():
    for k,v in (('A100',A100),('V1',V1),('V2',V2),('V3',ROOT),
                ('Input16k',INPUT),('Gas',GAS),('V3 RKF output',OUTPUT),
                ('V3 receipt',RECEIPT),('Matched16k',PAIRS),('Summary',SUMMARY),('Manifest',SEAL)):
        print(k,v.exists())


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('mode',choices=('selftest','buildcheck','run','compare','status'))
    args=ap.parse_args()
    if args.mode=='selftest':selftest()
    elif args.mode=='buildcheck':check_native()
    elif args.mode=='run':run()
    elif args.mode=='compare':compare()
    else:status()


if __name__=='__main__':main()