#!/usr/bin/env python3
"""Read-only posthoc diagnostic of SHA-sealed independent new120 LUT evaluation.

No Geant4, Heed, microscopic Garfield, training, hyperparameter search,
new bootstrap training-sample draws, or mutation of original results.
The source 120-parent clustered benchmark is already completed and frozen.
"""
from __future__ import annotations
import argparse
import csv
import json
import math
from collections import defaultdict
import numpy as np
import pi_external_muon120_guarded_lut_v1 as ext

MODES=('A20','A100','guarded_A20','guarded_A100')
GEOMS=ext.GEOMS
PER=ext.EVAL/'per_external_muon_4_models.csv'
SUM=ext.EVAL/'summary.json'
SEAL=ext.EVAL/'manifest.json'
REFSEAL=ext.REF/'manifest.json'
SEEDS=ext.REF/'muon_reference_seed.csv'

def data():
    ext.check()
    for p in (PER,SUM,SEAL,REFSEAL,SEEDS):
        if not p.is_file():
            raise FileNotFoundError(p)
    man=json.loads(SEAL.read_text())
    for p in (PER,SUM):
        if man['files_sha256'].get(p.name)!=ext.sha(p):
            raise ValueError('Frozen evaluation SHA mismatch: '+str(p))
    phys=json.loads(REFSEAL.read_text())
    if phys.get('status')!='complete' or phys.get('n_physically_distinct_muons')!=120:
        raise ValueError('Original full Garfield 120-parent source is not sealed')
    if phys['outputs_sha256'].get(SEEDS.name)!=ext.sha(SEEDS):
        raise ValueError('Original full Garfield electron seeds SHA differs')
    summ=json.loads(SUM.read_text())
    if summ.get('status')!='complete' or not summ.get('external_holdout'):
        raise ValueError('Not a successful independent external evaluation')
    if summ['source']['new_full_Garfield_manifest_sha256']!=ext.sha(REFSEAL):
        raise ValueError('Evaluation was not based on this exact full Garfield truth')
    with PER.open(newline='') as f:
        rows=list(csv.DictReader(f))
    if len(rows)!=480:raise ValueError('Expected exactly 480 linked parent×geometry results')
    for i in range(120):
        block=rows[4*i:4*i+4]
        if tuple(x['geometry'] for x in block)!=GEOMS or len({x['parent'] for x in block})!=1:
            raise ValueError('New120 parent×4 geometry grouping corrupted')
    if len({r['parent'] for r in rows})!=120:
        raise ValueError('Less than 120 distinct parent events')
    for m in MODES:
        if (sum(int(r[m+'_covered']) for r in rows)!=summ['overall'][m]['covered'] or
            not math.isclose(np.mean([float(r[m+'_MAE']) for r in rows]),
                             summ['overall'][m]['mean_abs_error'],rel_tol=1e-12,abs_tol=1e-5) or
            not math.isclose(np.mean([float(r[m+'_width']) for r in rows]),
                             summ['overall'][m]['mean_90PI_width'],rel_tol=1e-12,abs_tol=1e-5)):
            raise ValueError('Stored CSV does not reproduce original external aggregate '+m)
    return rows,summ

def measure(rr):
    out={}
    for m in MODES:
        coverage=[int(x[m+'_covered']) for x in rr]
        out[m]=dict(n=len(rr), covered=sum(coverage),
                    coverage_pct=round(100*np.mean(coverage),3),
                    mean_absolute_error=round(float(np.mean([float(x[m+'_MAE']) for x in rr])),2),
                    mean_90PI_width=round(float(np.mean([float(x[m+'_width']) for x in rr])),2))
    return out

def grouped_bootstrap(rows,reps=12000,seed=20261104):
    """Paired parent-level 95% percentile intervals, not 480 IID comparisons."""
    rng=np.random.default_rng(seed)
    ix=rng.integers(120,size=(reps,120))
    result={}
    for base,improved in [('A20','guarded_A20'),('A100','guarded_A100'),
                          ('guarded_A20','guarded_A100')]:
        vals={}
        for name,suffix,flip in [('coverage_pp','covered',1),
                                  ('MAE_reduction_ions','MAE',-1),
                                  ('PI_width_reduction_ions','width',-1)]:
            # 120 parent means, each over four correlated geometries.
            d=np.asarray([
                np.mean([float(rows[4*i+j][improved+'_'+suffix])-
                         float(rows[4*i+j][base+'_'+suffix])
                         for j in range(4)])
                for i in range(120)])
            d=flip*d*(100 if suffix=='covered' else 1)
            boot=d[ix].mean(axis=1)
            lo,hi=np.quantile(boot,[.025,.975])
            vals[name]={'estimate':round(float(d.mean()),4),
                        'parent_cluster_bootstrap95':[
                            round(float(lo),4),round(float(hi),4)]}
        result[improved+'_vs_'+base]=vals
    return result

def wire_stratification(rows):
    """Exploratory: exact within-muon Heed proximity to anode/FW centres.

    The existing guarded model uses bin centres and fixed rules; this
    diagnostic is purely observational and does not retune those rules.
    """
    near=defaultdict(lambda:[0,0,0])
    with SEEDS.open(newline='') as f:
        rd=csv.DictReader(f)
        for r in rd:
            g=r['configuration'];parent=r['crossing_id']
            if g not in GEOMS:raise ValueError('Unexpected original geometry')
            u=float(r['seed_u_mm']);v=float(r['seed_v_mm'])
            pitch=3. if g.startswith('p3_') else 4.
            nearest=round(u/pitch)*pitch
            da=math.hypot(u-nearest,v)
            dw=math.inf
            if 'FW75' in g:
                k=round(u/pitch-.5)
                dw=math.hypot(u-(k+.5)*pitch,v)
            dist=min(da,dw)
            z=near[(parent,g)]
            z[0]+=1
            z[1]+=int(dist <.35)   # study radius, not an electron-metal veto
            z[2]+=int(da<.015 or dw<.0375)
    if len(near)!=480:raise ValueError('Different parent × geometry count in electron source')
    result={}
    for g in GEOMS:
        original=[r for r in rows if r['geometry']==g]
        with_fractions=[]
        for r in original:
            k=(r['parent'],g)
            counts=near[k]
            if counts[0]!=int(r['n_primaries']):
                raise ValueError('Actual Garfield Heed seeds != evaluation primary count '+str(k))
            with_fractions.append((counts[1]/counts[0],r,counts))
        with_fractions.sort(key=lambda x:(x[0],x[1]['parent']))
        lo=[r for _,r,_ in with_fractions[:60]]
        hi=[r for _,r,_ in with_fractions[60:]]
        result[g]={
            'exact_Heed_primaries_total':sum(z[2][0] for z in with_fractions),
            'fraction_Heed_primaries_within_0_35mm_of_any_wire':round(
                sum(z[2][1] for z in with_fractions)/sum(z[2][0] for z in with_fractions),6),
            'n_primaries_starting_inside_wire_metal':sum(z[2][2] for z in with_fractions),
            'median_parent_fraction_near_wire':round(float(np.median([z[0] for z in with_fractions])),6),
            'lower_60_parents_proximity':measure(lo),
            'higher_60_parents_proximity':measure(hi),
        }
    return result

def selftest():
    # A parent contributes FOUR linked geometries; no IID 480 bootstrap.
    demo=[]
    for i in range(120):
        for g in GEOMS:
            demo.append(dict(parent=str(i),geometry=g,
                A20_covered='0',guarded_A20_covered='1',
                A100_covered='0',guarded_A100_covered='1',
                A20_MAE='10',guarded_A20_MAE='5',
                A100_MAE='10',guarded_A100_MAE='5',
                A20_width='30',guarded_A20_width='20',
                A100_width='30',guarded_A100_width='20'))
    v=grouped_bootstrap(demo,100,1)['guarded_A20_vs_A20']
    assert v['coverage_pp']['estimate']==100
    assert v['MAE_reduction_ions']['estimate']==5
    assert v['PI_width_reduction_ions']['estimate']==10
    assert v['coverage_pp']['parent_cluster_bootstrap95']==[100,100]
    print('PASS posthoc cluster pairing, coverage, MAE/width sign selftest; no physics',flush=True)

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('mode',choices=['selftest','audit'])
    a=ap.parse_args()
    if a.mode=='selftest':
        selftest();return
    rows,_=data()
    report={'schema':'NA60_MWPC_EXTERNAL_NEW120_READ_ONLY_POSTHOC_V1',
            'status':'read_only', 'n_parent_events':120,
            'n_correlated_muon_geometry_rows':480,
            'overall_from_sealed_CSV':measure(rows),
            'by_geometry':{g:measure([r for r in rows if r['geometry']==g]) for g in GEOMS},
            'paired_parent_95percent_uncertainty':grouped_bootstrap(rows),
            'targeted_A100_exposure':{
                'n_rows_in_target_A100_cells':sum(int(r['n_target'])>0 for r in rows),
                'n_rows_without_target_A100_cells':sum(int(r['n_target'])==0 for r in rows),
                'no_A100_target':measure([r for r in rows if int(r['n_target'])==0]),
                'with_A100_target':measure([r for r in rows if int(r['n_target'])>0])},
            'near_wire_observational_groups':wire_stratification(rows),
            'limitations':'Posthoc exploratory stratification, 120 independent parent clusters, fixed original physics and frozen model, not a further model-fitting or new-muon validation. 90% PI describes ion births, NOT induced U/V strip response.'}
    print(json.dumps(report,indent=2,ensure_ascii=False,allow_nan=False),flush=True)
    print('PASS READ-ONLY SHA-VERIFIED new120 diagnostics; original files untouched; no physics rerun',flush=True)

if __name__=='__main__':main()
