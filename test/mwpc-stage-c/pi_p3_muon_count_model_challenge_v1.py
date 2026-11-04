#!/usr/bin/env python3
"""p3 noFW and FW75 separately: archived new120 actual Geant4 muon count challenge.

Compare frozen guarded A20 count emulator against independently fitted
geometry-only zero+Gamma (Pólya) distribution, using EXACT same cached real
Geant4->TrackHeed primary electron populations and one full Garfield truth
per muon/geometry. Original physical, training and holdout data untouched.

This is a retrospective test on PREVIOUSLY OBSERVED external-new120; it is
not a newly blinded model-selection holdout. Neither candidate is a full
spatial-time ion cloud model (A20 preserves exact angular count vectors only).
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
import time
from pathlib import Path
import numpy as np
from scipy.stats import gamma as gammadist

import pi_external_muon120_guarded_lut_v1 as external
from pi_spatial_lut_muon120_crn_v1 import load_native,donors_for
from pi_audit_avalanche_lut_reuse import BASE,sha

GEOMS=('p3_noFW','p3_FW75')
OUT=BASE/'pi_p3_muon_count_model_challenge_v1'


def must(ok,msg):
    if not ok:raise ValueError(msg)
def write(p,obj):p.write_text(json.dumps(obj,sort_keys=True,indent=2,allow_nan=False)+'\n')
def stable(*terms):
    b='|'.join(map(str,terms)).encode()
    return int.from_bytes(hashlib.sha256(b).digest()[:8],'big')

def fit_geom_gamma(old20):
    x=np.asarray(old20,dtype=np.float64)[:,:,0].reshape(-1)
    must(len(x)==18000 and np.isfinite(x).all() and (x>=0).all(),
         'P3 A20 full 900x20 ion count source expected')
    positive=x[x>0]
    a,loc,scale=gammadist.fit(positive,floc=0)
    must(loc==0 and math.isfinite(a*scale) and a>0 and scale>0,
         'Invalid fitted positive-count gamma')
    return {'pzero':float((np.sum(x==0)+.5)/(len(x)+1)),
            'shape':float(a),'scale':float(scale),
            'mean_ion_count_per_primary':float(x.mean()),
            'train_n':len(x)}

def gamma_count(n,par,reps,rng):
    positive=rng.binomial(n,1-par['pzero'],size=reps)
    y=np.zeros(reps,float)
    mask=positive>0
    if mask.any():
        y[mask]=rng.gamma(par['shape']*positive[mask],par['scale'])
    return y

def guarded_count(counts,old20,g,reps,rng):
    # Original fixed 25% local / 75% physically guarded neighbors, no added
    # A100 outcomes and no weights depending on the full Garfield truth.
    vals=np.zeros(reps,float)
    ncell=len(old20)
    must(ncell==900,'Must be p3 geometry only')
    for box,n in counts:
        box=int(box);n=int(n)
        must(0<=box<ncell and n>0,'Invalid physical Heed cell count')
        inds=rng.integers(20,size=(reps,n))
        local=old20[box,inds,0]
        donors,cdf=donors_for(g,box,guard=True)
        if len(donors):
            choice=rng.random((reps,n))
            from_nei=(choice>=.25)
            chosen=rng.random((reps,n))
            target=donors[np.searchsorted(cdf,chosen,side='right').clip(0,len(donors)-1)]
            donoridx=rng.integers(20,size=(reps,n))
            val=np.where(from_nei,old20[target,donoridx,0],local)
        else:val=local
        vals+=np.sum(val,axis=1)
    return vals

def empirical_crps(x,y):
    # Proper scoring for a single held-out true muon draw; sorted O(RlogR).
    x=np.sort(np.asarray(x,float));n=len(x)
    return float(np.mean(abs(x-y))-
                 np.sum((2*np.arange(1,n+1)-n-1)*x)/(n*n))

def score(x,y):
    lo,med,hi,q99=np.quantile(x,[.05,.5,.95,.99])
    return {'truth':int(y),'forecast_mean':float(np.mean(x)),
            'forecast_median':float(med),'interval_low':float(lo),
            'interval_high':float(hi),'forecast_q99':float(q99),
            'MAE_predicted_mean':float(abs(np.mean(x)-y)),
            'CRPS':empirical_crps(x,y),
            'covered_PI90':int(lo<=y<=hi),'above95':int(y>hi),
            'above99':int(y>q99),'width_PI90':float(hi-lo),
            'pinball95':float(max(.95*(y-hi),-.05*(y-hi))),
            'pinball99':float(max(.99*(y-q99),-.01*(y-q99)))}

def cell_boot(vals,reps,rng):
    v=np.asarray(vals,float)
    draw=v[rng.integers(len(v),size=(reps,len(v)))].mean(axis=1)
    return {'value':float(v.mean()),'bootstrap95':list(map(float,np.quantile(draw,[.025,.975])))}

def load_original():
    source,metadata=external.check()
    _,nstats,receipt=external.legacy.verify_full()
    tables,targets,*_=load_native() # archived old20/new80 SHA and exact geometry IDs
    rec,nprim=external.build_external_records(tables,targets,source)
    must(nprim==6193 and len(rec)==480,'Unexpected new120 primary/geometry count')
    reference_manifest=external.REF/'manifest.json'
    return tables,rec,{'external_input_sha256':sha(external.INPUT/'native_muon_crossings_input.csv'),
                       'external_full_garfield_reference_sha256':sha(reference_manifest),
                       'A20_source_sha256':sha(BASE/'pi_lut500_avalanche_overnight_v1/run_20percell_firstpass_v1/manifest.json')}

def run(reps,bootstrap,seed):
    must(1000<=reps<=10000 and 200<=bootstrap<=10000,'Use 1000..10000 draw and 200..10000 bootstrap reps')
    must(not OUT.exists(),'Immutable model challenge exists; refusing overwrite: '+str(OUT))
    tables,rec,provenance=load_original()
    result={};allrows=[]
    for g in GEOMS:
        old,_=tables[g]
        params=fit_geom_gamma(old)
        records=[r for r in rec if r['geom']==g]
        must(len(records)==120 and len({r['parent'] for r in records})==120,
             'Not 120 independent real Geant4 parent identities for '+g)
        rows=[];times={}
        start=time.monotonic()
        for i,r in enumerate(records):
            key=r['parent']; n=r['n_primary']
            must(n==sum(nc for b,nc in r['counts']) and n>0,'Heed population incomplete')
            x=guarded_count(r['counts'],old,g,reps,
                            np.random.default_rng(stable(seed,g,key,'guarded')))
            y=gamma_count(n,params,reps,
                          np.random.default_rng(stable(seed,g,key,'gamma')))
            ref=r['ref_gain']
            rows.append({'geometry':g,'parent':key,'n_Heed_electrons':n,
                         'true_ion_births':ref,
                         'guarded_A20':score(x,ref),
                         'geometry_hurdle_gamma':score(y,ref)})
            if (i+1)%30==0:print('MUON CHALLENGE',g,'done',i+1,'of 120',flush=True)
        times['seconds_models_and_scoring']=time.monotonic()-start
        rng=np.random.default_rng(stable(seed,g,'bootstrap'))
        metrics=['MAE_predicted_mean','CRPS','covered_PI90','above95','above99',
                 'width_PI90','pinball95','pinball99']
        metrics_result={}
        for metric in metrics:
            metrics_result[metric]={
               mode:cell_boot([r[mode][metric] for r in rows],bootstrap,rng)
               for mode in ('guarded_A20','geometry_hurdle_gamma')}
            metrics_result[metric]['gamma_minus_guarded']=cell_boot([
               r['geometry_hurdle_gamma'][metric]-r['guarded_A20'][metric]
               for r in rows],bootstrap,rng)
        result[g]={'n_muons':120,'n_Heed_electrons':sum(r['n_Heed_electrons'] for r in rows),
                   'gamma_fit_from_original_A20':params,
                   'elapsed_seconds':times['seconds_models_and_scoring'],
                   'metric':metrics_result}
        allrows+=rows
    must(len(allrows)==240,'Two p3 independent geometry analyses incomplete')
    payload={'schema':'NA60_P3_MUON_COUNT_CHALLENGE_V1','status':'complete',
             'source':'Original 120 external new parent Geant4 events with SAME TrackHeed primary positions for two p3 geometries',
             'training':'Original pre-existing A20 only; no Gamma parameters fit using any external new120 Garfield outcomes',
             'warning':'The new120 test set was previously analyzed for guarded A20; it is NOT a newly untouched final blind model selection test.',
             'warning2':'This is only count response, neither cathode signals nor joint ion cloud; Gamma has no conditional angular ion distribution.',
             'sampling':'Fixed 120 parents evenly balanced over station and signed angle, NOT unbiased muon flux',
             'reps_per_muon':reps,'bootstrap_independent_parents':bootstrap,
             'provenance':provenance,'by_geometry':result}
    OUT.mkdir(parents=True)
    write(OUT/'summary.json',payload)
    write(OUT/'per_parent.json',allrows)
    write(OUT/'manifest.json',{
        'schema':'NA60_P3_MUON_COUNT_CHALLENGE_SEAL_V1',
        'sha256':{p.name:sha(p) for p in (OUT/'summary.json',OUT/'per_parent.json')},
        'source':provenance,'rng_seed':seed})
    for g in GEOMS:
        print('\n'+g+': independent muon count challenge (NOT 4-geometry pooled)')
        for metric in ('CRPS','MAE_predicted_mean','covered_PI90','above99','pinball99'):
            q=result[g]['metric'][metric]
            print('  '+metric,'guarded',round(q['guarded_A20']['value'],4),
                  'gamma',round(q['geometry_hurdle_gamma']['value'],4),
                  'gamma-minus-guarded',round(q['gamma_minus_guarded']['value'],4),
                  '95% parent CI',*[round(a,4) for a in q['gamma_minus_guarded']['bootstrap95']])
    print('PASS full-reference SHA, p3-only 120-parent independent-group scoring',OUT,flush=True)

def selftest():
    arr=np.zeros((900,20,5),dtype=np.int64)
    r=np.random.default_rng(42)
    arr[:,:,0]=r.choice([0,100,400,800],size=(900,20))
    p=fit_geom_gamma(arr)
    must(0<p['pzero']<1 and p['scale']>0 and p['shape']>0,'Bad Gamma count fit')
    for g in GEOMS:
        gc=guarded_count([(400,3),(450,2)],arr,g,1000,r)
        gg=gamma_count(5,p,1000,r)
        must(gc.shape==gg.shape==(1000,) and
             np.isfinite(gc).all() and np.isfinite(gg).all(),
             'Wrong joint primary count sampler')
        must(empirical_crps(gc,200)>=0,'CRPS must be nonnegative')
    print('PASS synthetic two-geometry count generation, parent-level CRPS, gamma fit; NO Garfield',flush=True)

def main():
    a=argparse.ArgumentParser(description=__doc__)
    a.add_argument('mode',choices=('selftest','audit'))
    a.add_argument('--reps',type=int,default=3000)
    a.add_argument('--bootstrap',type=int,default=4000)
    a.add_argument('--seed',type=int,default=20261104)
    args=a.parse_args()
    if args.mode=='selftest':selftest()
    else:run(args.reps,args.bootstrap,args.seed)
if __name__=='__main__':main()
