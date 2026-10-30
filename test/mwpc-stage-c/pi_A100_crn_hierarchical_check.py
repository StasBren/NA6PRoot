#!/usr/bin/env python3
"""Read-only, source-conditioned A20 vs targeted A100 physics-LUT diagnosis.

Stage C1: coupled common-random-number (CRN) sampling from original 20 and
20+80 empirical Garfield responses.  Unchanged cells are EXACTLY the same.
Stage C2: nonparametric cell-TRAINING bootstrap, one common reweighted cell
response pool per replica and geometry shared by ALL physical muons. This
adds model-training uncertainty to independent electron-response sampling;
it is NOT an independent Geant4 generalization study.

Inputs are immutable 84k A20, 16k A100 and 120 real Geant4/TrackHeed/
Garfield reference. No Geant4, TrackHeed, Garfield++, field modeling,
ion drift or cathode signal is run. Explicitly NO metal-start veto applied.

Normal WSL: python3 pi_A100_crn_hierarchical_check.py run
Published archival cross-check (no original WSL chunk dirs):
 python3 pi_A100_crn_hierarchical_check.py run --published-dir PATH --output PATH
The published mode validates public aggregate SHA digests but cannot check
original per-chunk manifests: use native mode for full provenance.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
import random
import statistics
import sys
from collections import Counter,defaultdict
from pathlib import Path

GEOMS=("p3_noFW","p3_FW75","p4_noFW","p4_FW75")
N_MU=120
DIM=5  # ion births, hemi up/down, top/bottom
SEED=20261030
FIELDS=("ion_births","ion_up","ion_down","ion_top","ion_bottom")


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda:f.read(1024*1024),b''):
            h.update(chunk)
    return h.hexdigest()


def csv_rows(path):
    with Path(path).open(newline='') as f:
        return list(csv.DictReader(f))


def quantile(values,p):
    x=sorted(values)
    if not x:return None
    at=(len(x)-1)*p
    low=int(math.floor(at)); high=int(math.ceil(at))
    return x[low]+(at-low)*(x[high]-x[low])


def asym(a,b):
    return (a-b)/(a+b) if a+b else None


def assert_pool(table,target):
    assert set(table)==set(GEOMS)
    assert len(target)==200
    for geom in GEOMS:
        ncell=900 if geom.startswith('p3_') else 1200
        assert len(table[geom])==ncell, f"wrong {geom} LUT cell cardinality"
        for box in range(ncell):
            vals=table[geom][box]
            want=100 if (geom,box) in target else 20
            if len(vals)!=want:
                raise ValueError(f"Incorrect independent full Garfield pool: {geom} {box} {len(vals)}!={want}")
            for v in vals:
                if len(v)!=DIM or min(v)<0 or v[1]+v[2]>v[0] or v[3]+v[4]>v[0]:
                    raise ValueError(f"Invalid original microscopic ion census {geom}/{box}")


def load_from_published(root):
    root=Path(root)
    old=json.loads((root/'A20_published_manifest.json').read_text())
    new=json.loads((root/'A100_published_manifest.json').read_text())
    cmp=json.loads((root/'compare_published_manifest.json').read_text())
    if (old['status']!='complete' or old['n_completed_avalanches']!=84000
        or new['status']!='complete' or new['new_full_Garfield_avalanches']!=16000
        or cmp['status']!='complete' or cmp['n_true_heldout_Geant4_muons']!=120):
        raise ValueError('Incomplete official physics publications')
    if (cmp['original_A20_manifest_sha256']!=new['original_A20_manifest_sha256']
        or cmp['targeted_A100_full_Garfield_manifest_sha256']!=sha(root/'A100_published_manifest.json')
        or cmp['original_A20_manifest_sha256']!=sha(root/'A20_published_manifest.json')):
        raise ValueError('Official source manifest cross-consistency failure')
    for geom in GEOMS:
        for prefix,manifest in [('trials_',old),('trials_extra_',new)]:
            name=f'{prefix}{geom}.csv'
            if sha(root/name)!=manifest['files_sha256'][name]:
                raise ValueError('Published source SHA mismatch '+name)
    for name in ['target_cells.csv','cell_ranking.csv']:
        if sha(root/name)!=new['files_sha256'][name]:
            raise ValueError('Published A100 target SHA mismatch '+name)
    if sha(root/'muon_A20_vs_targeted_A100.csv')!=cmp['outputs_sha256']['muon_A20_vs_targeted_A100.csv']:
        raise ValueError('Published comparison source SHA mismatch')
    full_output_hash={
        'native_muon_crossings_input.csv':'64338d9af0cff9b96e50b6d44b8f8b71fc5f29929aa89d71568c4d37384be9ba',
        'muon_reference_heed_primaries.csv':'66e0db1366f4995034fe8f32996d0afc0a9ea44ad8b0af1287bdc8280b9130a9',
        'muon_reference_crossing.csv':'45dc019079830dff083c8c148e273c36153fd1e4b1291ff20e287bcf8d5adc35'}
    for name,h in full_output_hash.items():
        if sha(root/name)!=h:
            raise ValueError('Published original full Geant4/Garfield source SHA mismatch '+name)
    target={(r['geometry'],int(r['box_id'])) for r in csv_rows(root/'target_cells.csv')}
    if len(target)!=200 or any(sum(g==geom for g,b in target)!=50 for geom in GEOMS):
        raise ValueError('A100 original top50 selection corrupt')
    table={}
    for geom in GEOMS:
        pools=defaultdict(list)
        for tag,expect in [('trials_',20),('trials_extra_',80)]:
            seen=defaultdict(set)
            for r in csv_rows(root/f'{tag}{geom}.csv'):
                b=int(r['box_id']); k=int(r['sample'])
                if k in seen[b] or int(r['capped'])!=0:
                    raise ValueError('Duplicate trial sample or capped physics')
                seen[b].add(k)
                if tag=='trials_extra_' and (geom,b) not in target:
                    raise ValueError('A100 physics outside frozen selected target')
                pools[b].append(tuple(int(r[k]) for k in FIELDS))
            if len(seen)!=(50 if tag=='trials_extra_' else (900 if geom.startswith('p3_') else 1200)):
                raise ValueError('Wrong number original physics cells')
            if any(len(used)!=expect for used in seen.values()):
                raise ValueError('Bad per-cell original physics count')
        table[geom]=dict(pools)
    assert_pool(table,target)
    input_mu=csv_rows(root/'native_muon_crossings_input.csv')
    primaries=csv_rows(root/'muon_reference_heed_primaries.csv')
    crosses=csv_rows(root/'muon_reference_crossing.csv')
    old_comparison=csv_rows(root/'muon_A20_vs_targeted_A100.csv')
    return table,target,input_mu,primaries,crosses,old_comparison,{
        'input_type':'Google Drive SHA-matched published aggregates (NOT native chunks)',
        'original_A20_manifest_sha256':cmp['original_A20_manifest_sha256'],
        'A100_manifest_sha256':cmp['targeted_A100_full_Garfield_manifest_sha256'],
        'full120_manifest_sha256':cmp['original_Geant4_120_full_Garfield_manifest_sha256'],
        'published_original_comparison_csv_sha256':sha(root/'muon_A20_vs_targeted_A100.csv'),
    }


def load_native():
    # Eager verification includes every exact A100 sample and all 40 original
    # chunk seals, old A20 84k SHA and original 120 full native Garfield chunks.
    from pi_compare_A100_targeted_vs_A20 import (validate_extra,load_extended,
        OUTPUT as CMP_ROOT)
    from pi_holdout120_A20_from_D500 import (verify_full,REF,INPUT)
    from pi_run_muon_lut_validation import (PRIMARY_FIELDS,read_csv)
    from pi_run_native_muon_b1 import CROSSING_COLUMNS
    from pi_A100_targeted_50cell_campaign import ROOT as NEW
    full_a100=validate_extra()
    lut,ha,hnew,target=load_extended(full_a100)
    originals,_,href=verify_full()
    cm=json.loads((CMP_ROOT/'manifest.json').read_text())
    if cm['status']!='complete' or cm['targeted_A100_full_Garfield_manifest_sha256']!=hnew or cm['original_A20_manifest_sha256']!=ha or cm['original_Geant4_120_full_Garfield_manifest_sha256']!=href:
        raise ValueError('Original comparison manifest wrong source SHA')
    for name,digest in cm['outputs_sha256'].items():
        if sha(CMP_ROOT/name)!=digest:
            raise ValueError('Original comparison output edited '+name)
    targetset={(r['geometry'],int(r['box_id'])) for r in full_a100[2]}
    if target != targetset:raise ValueError('Different native target selection')
    assert_pool(lut,target)
    rows=list(originals)
    pri=read_csv(REF/'muon_reference_heed_primaries.csv',PRIMARY_FIELDS)
    crossing=read_csv(REF/'muon_reference_crossing.csv',CROSSING_COLUMNS)
    # Native COMP file exact schema already checked by sha and unique identity below.
    comparison=csv_rows(CMP_ROOT/'muon_A20_vs_targeted_A100.csv')
    return lut,target,rows,pri,crossing,comparison,{
      'input_type':'native WSL per-chunk SHA and Geant4 120 physical validation',
      'original_A20_manifest_sha256':ha,
      'A100_manifest_sha256':hnew,
      'full120_manifest_sha256':href,
      'published_original_comparison_csv_sha256':sha(CMP_ROOT/'muon_A20_vs_targeted_A100.csv'),
    },NEW


def stable_id(mu):
    return '|'.join(str(mu[k]) for k in ('source','chunk','event','hit_index','track_id','station','det_id'))


def round_away(v):
    return math.floor(v+.5) if v>=0 else math.ceil(v-.5)


def bin_for(geom,mu,primary):
    pitch=3. if geom.startswith('p3_') else 4.
    station=int(mu['station'])
    width=307.2 if station==0 else 512.
    nw=math.floor(width/pitch+.5)
    phase=pitch*.5 if nw%2==0 else 0.
    u_in=float(mu['u_in_mm'])
    nearest=round_away((u_in-phase)/pitch)*pitch+phase
    phase=u_in-nearest
    u=10*float(primary['u_rel_cm'])+phase
    local_u=u-round_away(u/pitch)*pitch
    v=10*float(primary['v_abs_cm'])
    if not (math.isfinite(local_u) and math.isfinite(v) and -3<=v<=3):
        raise ValueError('Original Heed out of gas LUT range')
    nu=30 if geom.startswith('p3_') else 40
    iu=min(nu-1,max(0,math.floor((local_u+pitch/2)/.1)))
    iv=min(29,max(0,math.floor((v+3)/.2)))
    return int(iu*30+iv)


def records_from_source(table,target,source_mu,pri,cross,old_compare):
    pr=defaultdict(list)
    for row in pri: pr[row['crossing_id']].append(row)
    crossdict={(r['crossing_id'] if 'crossing_id' in r else r[next(iter(r))],r['configuration']):r for r in cross}
    olddict={(r['crossing_id'],r['configuration']):r for r in old_compare}
    ids=[stable_id(m) for m in source_mu]
    if len(ids)!=N_MU or len(set(ids))!=N_MU or len({(m['source'],m['chunk'],m['event']) for m in source_mu})!=N_MU:
        raise ValueError('Not 120 independent Geant4 parent sources')
    if len(pri)!=5881 or len(olddict)!=480 or len(crossdict)!=480:
        raise ValueError('Wrong original Geant4/Heed/physical comparator row counts')
    records=[]
    for mu,mid in zip(source_mu,ids):
        p=pr[mid]
        if not p:raise ValueError('No physical Heed source for '+mid)
        for geom in GEOMS:
            cellcounts=Counter(bin_for(geom,mu,e) for e in p)
            if any(b not in table[geom] for b in cellcounts):
                raise ValueError('Heed electron outside validated original LUT')
            selected=sum(n for b,n in cellcounts.items() if (geom,b) in target)
            truth=crossdict[mid,geom]
            original=olddict[mid,geom]
            birth=int(truth['ion_births'])
            if (int(original['full_reference_ion_births'])!=birth or
                int(original['n_real_Heed_electrons'])!=len(p) or
                int(original['n_electrons_in_selected_A100_cells'])!=selected):
                raise ValueError(f'Published Heed+Garfield true physics inconsistent for {mid}/{geom}')
            records.append({
                'crossing_id':mid,'geometry':geom,'counts':tuple(sorted(cellcounts.items())),
                'n_prim':len(p),'n_target':selected,
                'full_ion_births':birth,
                'full_A_hemi':asym(int(truth['ion_up']),int(truth['ion_down'])),
                'full_A_TB':asym(int(truth['ion_top']),int(truth['ion_bottom'])),
                'published_A20_200rep_covered':int(original['A20_original_90pct_covered']),
                'published_A100_200rep_covered':int(original['A100_targeted_90pct_covered']),
            })
    if len(records)!=480 or sum(r['n_target']==0 for r in records)!=158:
        raise ValueError('Changed original real muon spatial exposure')
    return records


def sample_crn_cell(arr,n,rng,oldvec,newvec):
    """Valid CRN coupling. At a target, A100 is 20% old + 80% extra;
    for an A100 old choice it is EXACT same physical old draw as A20.
    At non-target cell, both outcomes ALWAYS identical.
    """
    is_target=len(arr)==100
    for _ in range(n):
        orig=arr[rng.randrange(20)]
        if is_target:
            plus=orig if rng.randrange(5)==0 else arr[20+rng.randrange(80)]
        else:
            plus=orig
        for k in range(DIM):
            oldvec[k]+=orig[k]
            newvec[k]+=plus[k]


def pooled_result(vec):
    ion,up,down,top,bot=vec
    return ion,asym(up,down),asym(top,bot)


def crn_replay(records,table,reps,seed):
    out=[]
    for i,r in enumerate(records):
        geom=r['geometry'];mid=r['crossing_id']
        rng=random.Random(int.from_bytes(hashlib.sha256(f'CRN_SOURCE_CONDITIONED_V1\x1f{seed}\x1f{mid}\x1f{geom}'.encode()).digest()[:8],'big'))
        old=[[],[],[]];new=[[],[],[]]
        for j in range(reps):
            v0=[0]*DIM;v1=[0]*DIM
            for box,n in r['counts']:
                sample_crn_cell(table[geom][box],n,rng,v0,v1)
            z0=pooled_result(v0);z1=pooled_result(v1)
            if r['n_target']==0 and z0!=z1:
                raise AssertionError('CRN changed physical sample outside any A100 target')
            for k in range(3):
                if z0[k] is not None:old[k].append(z0[k])
                if z1[k] is not None:new[k].append(z1[k])
        out.append(quantile_summary(r,old,new,'CRN'))
        if (i+1)%40==0:
            print(f'PASS CRN {i+1}/480 physical muon×geometry, {reps} conditional replicas each',flush=True)
    return out


def hierarchical_replay(records,table,reps,seed):
    """For each global replica resample whole cell training pools ONCE; reuse
    for ALL 120 physical muons, respecting correlated finite training error.
    Within resampled cell, independent primary avalanche draws are made.
    Target old20 and A100new100 bootstraps independent *conditional on frozen*
    original empirical data; unchanged old/new pools and event draws identical.
    """
    samples=[([[],[],[]],[[],[],[]]) for r in records]
    rng=random.Random(int.from_bytes(hashlib.sha256(f'HIER_SOURCE_CONDITIONED_V1\x1f{seed}'.encode()).digest()[:8],'big'))
    for t in range(reps):
        boot={}  # (geom,box) -> (20 resampled indices, 100 indices or None)
        for k,r in enumerate(records):
            a=[0]*DIM;b=[0]*DIM
            geom=r['geometry']
            for cell,n in r['counts']:
                key=(geom,cell)
                original=table[geom][cell]
                if key not in boot:
                    sampled20=[rng.randrange(20) for _ in range(20)]
                    sampled100=([rng.randrange(100) for _ in range(100)]
                                if len(original)==100 else None)
                    boot[key]=(sampled20,sampled100)
                pool20,pool100=boot[key]
                for _ in range(n):
                    oldval=original[pool20[rng.randrange(20)]]
                    newval=oldval if pool100 is None else original[pool100[rng.randrange(100)]]
                    for c in range(DIM):
                        a[c]+=oldval[c];b[c]+=newval[c]
            z0=pooled_result(a);z1=pooled_result(b)
            if r['n_target']==0 and z0!=z1:
                raise AssertionError('Hier bootstrap changed unchanged A20/A100 pooled sample')
            for c in range(3):
                if z0[c] is not None:samples[k][0][c].append(z0[c])
                if z1[c] is not None:samples[k][1][c].append(z1[c])
        if (t+1)%250==0:
            print(f'PASS GLOBAL TRAINING BOOTSTRAP {t+1}/{reps}: identical sampled cell pool shared across 120 muons',flush=True)
    return [quantile_summary(r,old,new,'HIER') for r,(old,new) in zip(records,samples)]


def quantile_summary(record,old,new,stage):
    data={'crossing_id':record['crossing_id'],'geometry':record['geometry'],
          'n_Heed_primaries':record['n_prim'],'n_selected_A100_primaries':record['n_target']}
    for name,sample in [('A20',old),('A100',new)]:
        for variable,arr,ref in zip(('gain','A_hemi','A_TB'),sample,
                                     (record['full_ion_births'],record['full_A_hemi'],record['full_A_TB'])):
            if not arr or ref is None:
                # Never mislabel non-defined asymmetry as a coverage miss.
                data.update({f'{stage}_{name}_{variable}_covered':None,
                             f'{stage}_{name}_{variable}_p05':None,
                             f'{stage}_{name}_{variable}_p95':None})
                continue
            lo=quantile(arr,.05);hi=quantile(arr,.95)
            data[f'{stage}_{name}_{variable}_covered']=int(lo<=ref<=hi)
            data[f'{stage}_{name}_{variable}_p05']=lo
            data[f'{stage}_{name}_{variable}_p95']=hi
            data[f'{stage}_{name}_{variable}_mean']=statistics.fmean(arr)
    return data


def bootstrap_coverage(records,crn,hier,nres=15000):
    def values(rows,tag):
        return [r[f'{tag}_A20_gain_covered'] for r in rows],[r[f'{tag}_A100_gain_covered'] for r in rows]
    out={}
    for tag,rs in [('CRN',crn),('HIER',hier)]:
        geom_map={}
        for geom in GEOMS:
            a,b=values([v for v in rs if v['geometry']==geom],tag)
            geom_map[geom]={'old_covered':sum(a),'new_covered':sum(b),
                 'n':len(a),'old_coverage':statistics.mean(a),
                 'new_coverage':statistics.mean(b),
                 'new_minus_old':statistics.mean(y-x for x,y in zip(a,b))}
        # Same 120 parent muons used for four detector configurations, cluster them.
        d=[]
        for j in range(N_MU):
            block=rs[4*j:4*j+4]
            if tuple(v['geometry'] for v in block)!=GEOMS:
                raise ValueError('Broken four-fields-per-parent bootstrap grouping')
            d.append(sum(v[f'{tag}_A100_gain_covered']-v[f'{tag}_A20_gain_covered'] for v in block)/4)
        rng=random.Random(SEED+int(tag=='HIER'))
        boot=sorted(sum(d[rng.randrange(N_MU)] for _ in range(N_MU))/N_MU for _ in range(nres))
        out[tag]={'by_geometry':geom_map,
                  'overall_gain_A20_covered':sum(v[f'{tag}_A20_gain_covered'] for v in rs),
                  'overall_gain_A100_covered':sum(v[f'{tag}_A100_gain_covered'] for v in rs),
                  'overall_gain_new_minus_old':statistics.mean(d),
                  'overall_parent_clustered_bootstrap95':[quantile(boot,.025),quantile(boot,.975)],
                  'n_parent_clusters':N_MU}
    return out


def selftest():
    rng=random.Random(44)
    a0=[(0,0,0,0,0)]*20
    old=[0]*DIM;new=[0]*DIM
    sample_crn_cell(a0,500,rng,old,new)
    assert old==new==[0]*DIM
    target=a0+[(100,50,50,50,50)]*80
    old=[0]*DIM;new=[0]*DIM
    sample_crn_cell(target,20000,rng,old,new)
    assert old[0]==0 and abs(new[0]/20000-80)<2, (old,new)
    # Bootstrapping a shared finite pool produces positive inter-electron
    # covariance; for a muon with many electrons, must enlarge variance.
    table={g:{b:[(0,0,0,0,0)]*10+[(100,50,50,50,50)]*10
              for b in range(900 if g.startswith('p3_') else 1200)} for g in GEOMS}
    for b in range(900):table['p3_noFW'][b]=table['p3_noFW'][b]
    demo={'crossing_id':'synthetic','geometry':'p3_noFW',
          'counts':((0,12),),'n_prim':12,'n_target':0,
          'full_ion_births':600,'full_A_hemi':0.,'full_A_TB':0.}
    x=crn_replay([demo],table,1200,SEED)[0]
    y=hierarchical_replay([demo],table,1200,SEED)[0]
    assert x['CRN_A20_gain_p05']==x['CRN_A100_gain_p05']
    assert y['HIER_A20_gain_p05']==y['HIER_A100_gain_p05']
    assert (y['HIER_A20_gain_p95']-y['HIER_A20_gain_p05'] >
            (x['CRN_A20_gain_p95']-x['CRN_A20_gain_p05'])*1.05)  # synthetic covariance check
    print('PASS SYNTHETIC CRN identity, correct 20/80 mixture and cell-training bootstrap, NO physics',flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('mode',choices=('run','selftest'))
    p.add_argument('--reps',type=int,default=3000)
    p.add_argument('--training-reps',type=int,default=1200)
    p.add_argument('--published-dir',type=Path,default=None)
    p.add_argument('--output',type=Path,default=None)
    args=p.parse_args()
    if args.mode=='selftest':
        selftest();return
    if not 500<=args.reps<=10000 or not 500<=args.training_reps<=10000:
        p.error('Require 500..10000 independent conditional/global bootstrap replicas')
    if args.published_dir:
        table,target,mu,pri,cross,old,provenance=load_from_published(args.published_dir)
        out=args.output
        if out is None:p.error('Explicit --output is REQUIRED for read-only published archive mode')
    else:
        table,target,mu,pri,cross,old,provenance,root=load_native()
        out=args.output or root/'diagnose_CRN_A20_A100_training_bootstrap_v1'
    if out.exists():
        raise FileExistsError('Refusing to replace existing diagnostics '+str(out))
    records=records_from_source(table,target,mu,pri,cross,old)
    # All checks above precede any output write. CRN and bootstrap use only
    # frozen ion births and sectors, never run any microscopic physics.
    crn=crn_replay(records,table,args.reps,SEED)
    hier=hierarchical_replay(records,table,args.training_reps,SEED)
    summary=bootstrap_coverage(records,crn,hier)
    summary.update({'status':'complete','schema':'PI_A100_A20_COMMON_RNG_HIERARCHICAL_V1',
        'conditional_reps_per_muon_geom':args.reps,'global_cell_training_bootstrap_reps':args.training_reps,
        'n_independent_muon_parent_tracks':120,'n_correlated_geometry_comparisons':480,
        'common_random_numbers':'A20 old20 index shared with A100 when A100 chooses one of old20 (20%); unchanged cells bitwise identical. Target new80 chosen with 80% probability.',
        'cell_training_bootstrap':'Nonparametric with-replacement bootstrap of old20 and full100 per cell ONCE per global replica, shared over all muons per geom. Each Heed electron draws independently inside bootstrapped pool; unchanged cells share draws. A20 and A100 target pools bootstrapped separately (not CRN-optimal but valid marginal).',
        'caveats':'Same 120 Geant4/TrackHeed source coordinates informed A100 cell ranking; only one true full Garfield realization per muon and geometry. Training bootstrap is a statistical approximation and not calibrated external validation; bootstrap values derived from frozen original 20/100 trials, no physics and no exact metal-start veto. Most of 4200 cells have only n=20; skew and binning remain.',
        'provenance':provenance})
    rows=[]
    for r,x,y in zip(records,crn,hier):
        row={k:v for k,v in r.items() if k!='counts'}
        if (row['crossing_id'],row['geometry'])!=(x['crossing_id'],x['geometry']) or (row['crossing_id'],row['geometry'])!=(y['crossing_id'],y['geometry']):
            raise ValueError('Lost real source identity during paired analysis')
        row.update({k:v for k,v in x.items() if k not in row})
        row.update({k:v for k,v in y.items() if k not in row})
        rows.append(row)
    out.mkdir(parents=True)
    per=out/'per_muon_C_RNG_and_training_bootstrap.csv'
    with per.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    (out/'summary.json').write_text(json.dumps(summary,indent=2,sort_keys=True)+'\n')
    (out/'manifest.json').write_text(json.dumps({'schema':summary['schema'],'status':'complete',
        'files_sha256':{f.name:sha(f) for f in [per,out/'summary.json']}},sort_keys=True,indent=2)+'\n')
    for mode in ('CRN','HIER'):
        v=summary[mode]
        print(mode,'A20',v['overall_gain_A20_covered'],'/480; A100',v['overall_gain_A100_covered'],
              '/480; delta',round(v['overall_gain_new_minus_old'],6),'bootstrap95',v['overall_parent_clustered_bootstrap95'],flush=True)
        for g,val in v['by_geometry'].items():print('RESULT',mode,g,val['old_covered'],val['new_covered'],flush=True)
    print('PASS SOURCE-CONDITIONED CRN AND HIERARCHICAL TRAINING BOOTSTRAP; immutable original Garfield sources unchanged; SHA-SEALED',out,flush=True)


if __name__=='__main__':main()