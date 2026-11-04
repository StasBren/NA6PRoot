#!/usr/bin/env python3
"""Read-only per-geometry, per-spatial-cell 32/128/256/512/1024 readout.

Verifies ALL physics SHA checkpoints and completed v2 exact-cloud analysis.
This score is compression-of-a-known-cloud only; NOT a trained joint LUT
or cathode signal. Bootstrap resampling unit = independent spatial BOX,
not independent ion and not 32 repeated reservoir draws.
"""
from __future__ import annotations
import argparse
import contextlib
import io
import json
from collections import defaultdict
from pathlib import Path
from unittest.mock import patch

import numpy as np

import pi_p3_extended_full_cloud_campaign_v1 as campaign

GEOMS=campaign.GEOMS
CAPS=campaign.CAPS
ROOT=campaign.ROOT
METRICS=(
 ('joint_phi_time_TV','angular × birth-time distribution: TV'),
 ('abs_top_frac','top-sector birth probability: absolute error'),
 ('mean_radius_rel_abs','mean ion-birth radius: relative error'),
 ('abs_q90_t_ns','90th-percentile ion-birth time: absolute ns'),
 ('abs_A_TB','top-bottom ion-birth asymmetry: absolute error'),
)

def must(x,msg):
    if not x:raise ValueError(msg)
def load(p):return json.loads(p.read_text())

def checked(root):
    must(root==ROOT,'This SHA-verified readout is bound to the extended p3 campaign directory')
    campaign.configure_core()
    rows,tags,source,done=campaign.core.completed()
    must(len(done)==len(rows)//campaign.GROUP and
         len(rows)==2*source['per_geometry'],
         'Extended Garfield physical runs not all SHA sealed')
    target=root/'analysis_complete_v2'
    manifest=load(target/'manifest.json')
    must(manifest['schema']=='NA60_P3_FULL_CLOUD_ANALYSIS_SHA_V2' and
         manifest['selection_sha256']==source['input_sha256'] and
         len(manifest['input_chunk_sha256'])==len(done),
         'Extended full-cloud analysis/checkpoints are incomplete')
    for name,dig in manifest['input_chunk_sha256'].items():
        must(campaign.sha(root/'full_cloud_chunks'/name)==dig,
             'Native full Garfield chunk changed '+name)
    for name,dig in manifest['output_sha256'].items():
        must(campaign.sha(target/name)==dig,'Analysis changed '+name)
    summary=load(target/'summary.json')
    cloud=load(target/'per_avalanche.json')
    must(summary['schema']=='NA60_P3_FULL_BIRTH_COMPRESS_ANALYSIS_V2' and
         summary['selected_starts']==len(rows)==len(cloud) and
         summary['completed_chunks']==len(done) and
         summary['reps_per_cap_per_avalanche']>=8 and
         set(summary['summary_by_geometry'])==set(GEOMS),
         'Unexpected spatial/time compression analysis schema')
    expected={(r['geometry'],int(r['box_id']),int(r['sample'])):t
              for r,t in zip(rows,tags)}
    seen=set();g_counts=defaultdict(int);g_births=defaultdict(int)
    for row in cloud:
        key=(row['geometry'],int(row['box_id']),int(row['sample']))
        must(key in expected and key not in seen, 'Duplicate/wrong D500 avalanche: '+str(key))
        seen.add(key)
        must(row['regime']==expected[key]['regime'] and
             row['v_sign']==expected[key]['v_sign'],'Wrong source stratum')
        n=int(row['n_ion_births'])
        must(n>=0,'Negative full-ion census')
        if n:must(set(row['cap_study'])==set(map(str,CAPS)),
                  'Expected complete five-cap 32/128/256/512/1024 analysis')
        else:must(not row['cap_study'] and row['exact_A_TB'] is None,
                  'Invalid non-physical zero-ion distribution')
        g_counts[row['geometry']]+=1;g_births[row['geometry']]+=n
    must(set(expected)==seen,'Some extended physical starts missing')
    for geom in GEOMS:
        must(g_counts[geom]==source['per_geometry'] and
             g_births[geom]==summary['summary_by_geometry'][geom]['total_exact_ion_births'],
             'Geometry-specific source/birth count mismatch')
    must(sum(g_births.values())==sum(c['n_ion_births'] for c in done),
         'All original Garfield ion births do not sum to the analyzed clouds')
    return cloud,source,summary

def error(row,cap,key):
    if row['n_ion_births']==0:return None
    item=row['cap_study'][str(cap)][key]
    if item is None:return None
    v=float(item['mean'])
    must(np.isfinite(v) and v>=0,'Invalid physical compression error')
    return v

def clustered(rows,cap,key,bootstrap,rng,compressed_only=False,
              ion_birth_weight=False):
    """Bootstrap WHOLE independent spatial cells, retaining their trials."""
    bycell=defaultdict(list)
    for row in rows:
        if row['n_ion_births']<=0:continue
        if compressed_only and row['n_ion_births']<=cap:continue
        z=error(row,cap,key)
        if z is not None:bycell[int(row['box_id'])].append(
            (z,int(row['n_ion_births'])))
    if not bycell:return None
    boxes=sorted(bycell)
    dat=[bycell[b] for b in boxes]
    n=sum(len(v) for v in dat)
    if ion_birth_weight:
        num=np.array([sum(e*weight for e,weight in v) for v in dat],float)
        den=np.array([sum(weight for e,weight in v) for v in dat],float)
    else:
        # Each spatial box gets equal weight even with many independent
        # microscopic starts, to avoid accidentally changing source PDF.
        num=np.array([np.mean([e for e,weight in v]) for v in dat],float)
        den=np.ones(len(dat),float)
    avg=float(np.sum(num)/np.sum(den))
    ci=None
    if len(boxes)>=2:
        ix=rng.integers(0,len(boxes),size=(bootstrap,len(boxes)))
        means=np.sum(num[ix],axis=1)/np.sum(den[ix],axis=1)
        ci=[float(x) for x in np.quantile(means,[.025,.975])]
    return {'mean':avg,'n_spatial_cells':len(boxes),
            'n_independent_avalanches':n,'cluster_CI95':ci}

def paired(rows,small,big,metric,bootstrap,rng,compressed_only=False):
    bycell=defaultdict(list)
    for row in rows:
        if row['n_ion_births']<=0:continue
        if compressed_only and row['n_ion_births']<=big:continue
        a=error(row,small,metric);b=error(row,big,metric)
        if a is not None and b is not None:
            bycell[int(row['box_id'])].append(a-b)
    if not bycell:return None
    z=np.array([np.mean(v) for _,v in sorted(bycell.items())],float)
    ci=None
    if len(z)>=2:
        idx=rng.integers(len(z),size=(bootstrap,len(z)))
        ci=list(map(float,np.quantile(z[idx].mean(axis=1),[.025,.975])))
    return {'mean_improvement':float(np.mean(z)),
            'n_spatial_cells':len(z),'n_avalanches':sum(map(len,bycell.values())),
            'cluster_CI95':ci}

def fmt(x):
    if x is None:return 'no defined avalanche'
    ci=x['cluster_CI95']
    return (f"{x.get('mean',x.get('mean_improvement')):.6g} "
            f"CI95={ci if ci is not None else 'unavailable (one cell)'} "
            f"[boxes={x['n_spatial_cells']}, avalanches={x.get('n_independent_avalanches',x.get('n_avalanches'))}]")

def emit(cloud,bootstrap,seed):
    rng=np.random.default_rng(seed)
    print('P3 EXTENDED FULL-BIRTH CLOUD COMPRESSION (no induced signal)')
    print('Each model geometry and each spatial regime is scored SEPARATELY.')
    print('Independent units for uncertainty: D500 spatial cells; not ions or 32 repeated reservoir samples.')
    for geom in GEOMS:
        allrows=[r for r in cloud if r['geometry']==geom]
        active=[r for r in allrows if r['n_ion_births']>0]
        print('\n'+'='*75)
        print(f'{geom}: {len(allrows)} exact D500 starting electrons, '
              f'{len(active)} nonzero and {len(allrows)-len(active)} zero; '
              f'{sum(r["n_ion_births"] for r in allrows)} total ion births')
        for cap in CAPS:
            print(f' cap={cap}: '+str(sum(r['n_ion_births']>cap for r in active))+
                  f' genuinely compressed nonzero avalanches / {len(active)}')
        regimes=['near_anode','transition','bulk'] if geom=='p3_noFW' else [
                 'near_anode','near_field','transition','bulk']
        for regime in regimes:
            rows=[r for r in active if r['regime']==regime]
            print(f'\n   {regime}: independent nonzero n={len(rows)}, '
                  f'distinct spatial cells={len(set(r["box_id"] for r in rows))}')
            for cap in CAPS:
                s=clustered(rows,cap,'joint_phi_time_TV',bootstrap,rng)
                t=clustered(rows,cap,'joint_phi_time_TV',bootstrap,rng,
                            compressed_only=True)
                print(f'     TV cap={cap}: {fmt(s)}; compressed-only: {fmt(t)}')
        print('\n  ALL REGIMES, SPATIAL-CELL BOOTSTRAP:')
        for key,title in METRICS:
            print('   '+title)
            for cap in CAPS:
                print(f'      {cap}: {fmt(clustered(active,cap,key,bootstrap,rng))}')
        print('\n  ION-BIRTH-WEIGHTED OVER FULL CLOUDS (bootstrap by cells):')
        for key in ('joint_phi_time_TV','mean_radius_rel_abs','abs_q90_t_ns'):
            print('   '+key)
            for cap in CAPS:
                print(f'      {cap}: {fmt(clustered(active,cap,key,bootstrap,rng,ion_birth_weight=True))}')
        print('\n  PAIRED BETTER COMPRESSION (same microscopic avalanche):')
        for small,big in ((128,512),(512,1024)):
            for key in ('joint_phi_time_TV','mean_radius_rel_abs','abs_q90_t_ns'):
                print(f'     {small} -> {big} {key}: '+
                      fmt(paired(active,small,big,key,bootstrap,rng)))
                print('       only Nbirth> '+str(big)+': '+
                      fmt(paired(active,small,big,key,bootstrap,rng,True)))
    print('\nCAUTION: cap-reconstruction error is NOT error of a new spatially conditioned avalanche LUT.')
    print('No information from either p3 geometry is pooled into the other.')

def report(bootstrap,seed):
    cloud,source,summary=checked(ROOT)
    print('PASS SHA: every native physical chunk and separate extended compression analysis')
    print('Compressor used',summary['reps_per_cap_per_avalanche'],
          'random draws PER avalanche; these are NOT independent avalanches.')
    emit(cloud,bootstrap,seed)
    print('PASS p3 extended full-ion-cloud compression report; read-only')

def selftest():
    vals=[]
    for geom in GEOMS:
        regs=['near_anode','transition','bulk'] if geom=='p3_noFW' else [
              'near_anode','near_field','transition','bulk']
        for i in range(24):
            n=0 if i%11==0 else (50,160,350,600,3200)[i%5]
            shape={}
            if n:
                for cap in CAPS:
                    shape[str(cap)]={key:{'mean':(j+1)*.01*max(0,n-cap)/n}
                                     for j,(key,_) in enumerate(METRICS)}
                if i%7==0:
                    for cap in CAPS:shape[str(cap)]['abs_A_TB']=None
            vals.append({'geometry':geom,'box_id':i//3,'sample':100+i,
                         'regime':regs[i%len(regs)],'n_ion_births':n,
                         'cap_study':shape})
    for key,_ in METRICS:
        must(error(vals[1],128,key) is not None,'Metric mapping')
        must(clustered([v for v in vals if v['geometry']=='p3_noFW'],128,
                       key,200,np.random.default_rng(3)) is not None,
             'Spatial cell bootstrap failed')
    out=io.StringIO()
    with contextlib.redirect_stdout(out):emit(vals,200,20261104)
    s=out.getvalue()
    must(all(x in s for x in ('p3_noFW','p3_FW75','512 -> 1024',
          'ION-BIRTH-WEIGHTED','P3 EXTENDED FULL-BIRTH CLOUD')),
         'End-to-end five-cap two-geometry renderer failed')
    print('PASS synthetic 2 independent p3 geometries, full 5-cap report, '
          'physical zeros, undefined sectors, spatial-cell bootstrap; NO Garfield')

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('mode',choices=('selftest','report'))
    p.add_argument('--bootstrap',type=int,default=3000)
    p.add_argument('--seed',type=int,default=20261104)
    a=p.parse_args()
    must(200<=a.bootstrap<=10000,'Require 200..10000 bootstrap draws')
    if a.mode=='selftest':selftest()
    else:report(a.bootstrap,a.seed)
if __name__=='__main__':main()
