#!/usr/bin/env python3
"""NA60+ separate staged full-microscopic Geant4/TrackHeed/Garfield campaign.

Never overwrites old120/new120, old A20/A100, or original native C++ target.
Prepares an angle/station-balanced source from existing SHA-frozen J/psi
Geant4 chunk2, one gas crossing per independent generated parent event.

Stage 80% of selected parents for descriptive precision monitoring; the last
20% are an explicitly RESERVED validation set and are not analyzed by the
ordinary 'report' mode. All 4 geometries use identical TrackHeed electrons.
The separate spacetime_v1 C++ target samples up to 128 birth space/time
points per primary with independent RNG; complete ion count/sector census
uses ALL microscopic births. No induced cathode signal is calculated here.
"""
from __future__ import annotations
import argparse
import csv
import json
import math
import subprocess
import time
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from pi_2m_native_six_geometry_fold import crossing_provenance, SOURCE
from pi_prepare_muon_avalanche_b1 import FIELDS, validate, sha, seed
from pi_garfield_runtime_env import garfield_runtime_env
from pi_run_native_muon_b1 import check_crossings
from pi_external_muon120_wire_geometry_recovery_v1 import isolated_validator
from pi_holdout120_A20_from_D500 import read_csv, write_csv, write_json, id_of

REPO=Path(__file__).resolve().parents[2]
BASE=REPO/'test_runs/mwpc_stage_c'
ROOT=BASE/'pi_adaptive_muon_spacetime_v1'
INPUT=ROOT/'input'
REF=ROOT/'reference_chunks'
REPORTS=ROOT/'reports'
THIS=Path(__file__).resolve()
CXX=REPO/'test/mwpc-avalanche/phaseB_native_muon_3p3_spacetime_v1.cxx'
EXE=REPO/'test/mwpc-avalanche/build/mwpc_phase_b_native_muon_3p3_spacetime_v1'
SOURCE_CHUNK=2
GROUP=10
GEOMS=('p3_noFW','p3_FW75','p4_noFW','p4_FW75')
ST=( (station,side) for station in range(6) for side in (-1,1) )
DATA_FILES=('muon_reference_crossing.csv','muon_reference_seed.csv',
            'muon_reference_ion_birth_sample.csv',
            'muon_reference_heed_primaries.csv','muon_reference_heed_clusters.csv')


def load_json(p):return json.loads(p.read_text())
def check_binary():
    if not CXX.is_file() or not EXE.is_file():
        raise FileNotFoundError('Build NEW spacetime target first (see campaign README): '+str(EXE))
    if 'cacheMuonPhysics?128:100' not in CXX.read_text():
        raise RuntimeError('Expected frozen 128-point spatial/time C++ target absent')
    return {'cpp_source_sha256':sha(CXX),'cpp_executable_sha256':sha(EXE)}


def prior_parents():
    locations=(
       BASE/'pi_muon_holdout120_A20_from_D500_v1/input_chunk2_120',
       BASE/'pi_external_new120_guarded_lut_blind_v1/input_unseen_chunk2_120')
    result=set();verified=[]
    for d in locations:
        p=d/'native_muon_crossings_input.csv';m=d/'input_manifest.json'
        if not p.exists() or not m.exists():
            raise FileNotFoundError('Protected original holdout source missing: '+str(d))
        meta=load_json(m)
        if (meta.get('status')!='complete' or meta.get('n_muons')!=120
            or sha(p)!=meta.get('sha256_input')):
            raise ValueError('Prior source manifest integrity mismatch '+str(d))
        rows=read_csv(p,FIELDS)
        if len(rows)!=120:raise ValueError('Wrong frozen holdout cardinality')
        for r in rows:
            k=(str(r['source']),str(r['chunk']),str(r['event']))
            if k in result:raise ValueError('Old/new holdout parent IDs overlap')
            result.add(k)
        verified.append({'path':str(p),'sha256':sha(p)})
    if len(result)!=240:raise ValueError('Expected exactly 240 excluded parents')
    return result,verified


def parent_key(r):return (str(r['source']),str(r['chunk']),str(r['event']))

def select(pool,per_stratum):
    """Unbiased SHA-ranked within each eligible stratum, one parent total.

    All 12 strata participate in each selection cycle; enough extra
    candidates are retained for cross-station parent conflicts.
    """
    picked=[];used=set();n=Counter()
    for _ in range(per_stratum):
        for station in range(6):
            for side in (-1,1):
                k=(station,side)
                row=next((rec for _,_,rec in pool[k] if parent_key(rec) not in used),None)
                if row is None:
                    raise RuntimeError('Not enough independent eligible Geant4 parents in stratum '+str(k))
                used.add(parent_key(row));picked.append(row);n[k]+=1
    if len(picked)!=12*per_stratum or len(used)!=len(picked):
        raise AssertionError('Independent/stratified selection failed')
    return picked,n


def prepare(n):
    if n%120!=0 or n<600 or n>4800:
        raise ValueError('Use 600..4800 events in multiples of 120; 10-parent chunks, 12 strata')
    if ROOT.exists():raise FileExistsError('New campaign already exists: '+str(ROOT))
    exclude,proof=prior_parents()
    path,srcproof=crossing_provenance('Jpsi',SOURCE_CHUNK)
    if sha(path)!=srcproof['csv_sha256']:
        raise ValueError('Frozen Geant4 source CSV SHA mismatch')
    max_per=n//12
    pools=defaultdict(list);counts=Counter()
    with path.open(newline='') as f:
        reader=csv.DictReader(f)
        required=set(FIELDS)-{'source','chunk','momentum_GeV','heed_seed'}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError('Missing original native Geant4 columns')
        for r in reader:
            if (SOURCE['Jpsi'],str(SOURCE_CHUNK),str(r['event'])) in exclude:continue
            checked=validate(r)
            if checked is None:continue
            p,v=checked
            if abs(v['u_out_mm']-v['u_in_mm'])>8 or abs(v['w_out_mm']-v['w_in_mm'])>40:
                continue
            station=int(r['station'])
            if station not in range(6):continue
            sign=-1 if float(r['theta_u_in_deg'])<0 else 1
            key=(station,sign);counts[key]+=1
            orig=(SOURCE['Jpsi'],SOURCE_CHUNK,r['event'],r['hit_index'],
                  r['track_id'],station,r['det_id'])
            q=seed('NA60-SPACETIME-PARENT-20261104',*orig)
            rec={c:r[c] for c in FIELDS if c not in ('source','chunk','momentum_GeV','heed_seed')}
            rec.update({'source':SOURCE['Jpsi'],'chunk':SOURCE_CHUNK,
                        'momentum_GeV':f'{p:.12g}',
                        'heed_seed':seed('NA60-SPACETIME-HEED-20261104',*orig)})
            cur=pools[key];cur.append((q,tuple(map(str,orig)),rec))
            cur.sort(key=lambda x:(x[0],x[1]))
            if len(cur)>max_per*3+80:cur.pop()
    selected,per=select(pools,max_per)
    holdout=12*max(1,max_per//5)
    # Preselect train before a completely uninspected tail of 20% test parents.
    # Since selection cycles each include 12 strata, split stays balanced.
    assert holdout%12==0 and (n-holdout)%GROUP==0
    ROOT.mkdir(parents=True)
    INPUT.mkdir()
    out=INPUT/'native_muon_crossings_input.csv'
    write_csv(out,FIELDS,selected)
    meta={'schema':'NA60_STAGED_MICROSCOPIC_SPACETIME_V1','status':'prepared',
          'n_parents':n,'dev_parents':n-holdout,'reserved_holdout_parents':holdout,
          'n_per_station_signed_angle':max_per,
          'n_dev_per_stratum':(n-holdout)//12,
          'n_holdout_per_stratum':holdout//12,
          'source_chunk':SOURCE_CHUNK,'source_id':SOURCE['Jpsi'],
          'source_geant4_sha256':srcproof['csv_sha256'],
          'source_root_sha256':srcproof['immutable_input']['input_root_sha256'],
          'input_sha256':sha(out),'excluded_old_new_240':proof,
          'eligible_source_counts':{f'{a}_{b}':counts[a,b] for a,b in ST},
          'selection':'SHA ranking within station x signed theta-u; 1 parent event total; '
                      'valid finite 9-wire analytical patch only; separate 20pct untouched holdout',
          'note':'Balanced strata are not a natural Jpsi flux sample; reweight '
                 'by eligible source stratum counts to estimate source-weighted moments.',
          'script_source_sha256':sha(THIS)}
    write_json(INPUT/'input_manifest.json',meta)
    check_input()
    print('PASS PREPARED',n,'unseen Geant4 parents;',n-holdout,
          'development /',holdout,'reserved untouched validation;',n//GROUP,'chunks',flush=True)


def check_input():
    m=load_json(INPUT/'input_manifest.json')
    p=INPUT/'native_muon_crossings_input.csv'
    if (m['schema']!='NA60_STAGED_MICROSCOPIC_SPACETIME_V1'
        or m['status']!='prepared' or sha(p)!=m['input_sha256']
        or sha(THIS)!=m['script_source_sha256']):
        raise ValueError('Campaign input/protocol changed since frozen prepare')
    forbidden,proof=prior_parents()
    if proof!=m['excluded_old_new_240']:raise ValueError('Prior protected source changed')
    path,s=crossing_provenance('Jpsi',SOURCE_CHUNK)
    if s['csv_sha256']!=m['source_geant4_sha256'] or sha(path)!=m['source_geant4_sha256']:
        raise ValueError('Frozen Geant4 source changed')
    rr=read_csv(p,FIELDS)
    if len(rr)!=m['n_parents'] or len({parent_key(x) for x in rr})!=len(rr):
        raise ValueError('Duplicate parent or missing input')
    if any(parent_key(x) in forbidden or validate(x) is None for x in rr):
        raise ValueError('Invalid or previously used source parent')
    return rr,m


def chunk_dir(i):return REF/f'chunk_{i:04d}'
def source_split(rows):return [rows[i:i+GROUP] for i in range(0,len(rows),GROUP)]

def preflight():
    rows,m=check_input();binary=check_binary()
    first=rows[:GROUP]
    prefix=INPUT/'preflight10_native_input.csv'
    if not prefix.exists():write_csv(prefix,FIELDS,first)
    if [id_of(x) for x in read_csv(prefix,FIELDS)]!=[id_of(x) for x in first]:
        raise ValueError('Existing preflight subset differs from frozen selection')
    for mode in ('--validate-input','--validate-heed'):
        proc=subprocess.run([str(EXE),mode,str(prefix)],
                            capture_output=True,text=True,env=garfield_runtime_env())
        if proc.returncode or 'PASS' not in proc.stdout:
            raise RuntimeError('Preflight failed '+mode+'\n'+proc.stdout[-4000:]+'\n'+proc.stderr[-4000:])
        print('PASS',mode,proc.stdout.strip().splitlines()[-1],flush=True)
    record={'schema':'NA60_SPACETIME_PREFLIGHT_V1','input_sha256':m['input_sha256'],
            'preflight_csv_sha256':sha(prefix),**binary}
    fn=INPUT/'preflight.json'
    if fn.exists() and load_json(fn)!=record:raise ValueError('Existing preflight changed')
    if not fn.exists():write_json(fn,record)
    print('PASS PRE-FLIGHT 10 actual Heed gas tracks; no full avalanches',flush=True)


def expected_binary():
    m=load_json(INPUT/'preflight.json');b=check_binary()
    if any(m[k]!=v for k,v in b.items()):
        raise ValueError('Native physics source or binary changed after preflight')
    return b


def check_birth_reservoir(folder):
    """Check count/weight/coordinate conservation of the NEW 128-point sample."""
    seeds={};counts=Counter()
    with (folder/'muon_reference_seed.csv').open(newline='') as f:
        for r in csv.DictReader(f):
            k=r['crossing_id'],r['configuration'],int(r['primary_index'])
            if k in seeds:raise ValueError('Duplicate microscopic seed')
            seeds[k]=int(r['ion_births'])
    with (folder/'muon_reference_ion_birth_sample.csv').open(newline='') as f:
        for r in csv.DictReader(f):
            k=r['crossing_id'],r['configuration'],int(r['primary_index'])
            if k not in seeds:raise ValueError('Sampled ion birth has no microscopic seed')
            n=seeds[k]
            if n!=int(r['seed_ion_births_total']) or n<=0:
                raise ValueError('Sampled ion weight/total inconsistent')
            for c in ('u_mm','v_mm','w_mm','t_ns'):
                if not math.isfinite(float(r[c])):raise ValueError('Nonfinite birth coordinate')
            counts[k]+=1
            if counts[k]>min(n,128):raise ValueError('More ion births sampled than reservoir cap')
            w=float(r['reservoir_weight'])
            want=n/min(n,128)
            if abs(w-want)>1e-8*max(1,want):
                raise ValueError('Weighted birth sample normalization wrong')
    for k,n in seeds.items():
        if counts[k]!=min(n,128):
            raise ValueError('Incomplete per-avalanche weighted space/time birth reservoir')
    return sum(counts.values())


def audit_chunk(i,subset,binary):
    d=chunk_dir(i);seal=d/'complete.json'
    if not seal.is_file():
        raise RuntimeError('Exists but unsealed checkpoint: '+str(d))
    m=load_json(seal)
    if (m['status']!='complete' or m['chunk']!=i or m['n_parents']!=len(subset)
        or m['physics']!=binary):raise ValueError('Checkpoint identity/binary changed')
    for filename,digest in m['sha256_files'].items():
        if sha(d/filename)!=digest:raise ValueError('Checkpoint file SHA changed: '+str(d/filename))
    if [id_of(x) for x in read_csv(d/'native_muon_crossings_input.csv',FIELDS)] != [id_of(x) for x in subset]:
        raise ValueError('Checkpoint source muons changed')
    return m


def run(max_chunks,hours):
    rr,meta=check_input();binary=expected_binary()
    batches=source_split(rr)
    REF.mkdir(parents=True,exist_ok=True)
    configfile=REF/'run_config.json'
    config={'schema':'NA60_SPACETIME_FULL_MICRO_RUN_V1',
            'input_sha256':meta['input_sha256'],'physics':binary,
            'max_birth_sample_per_primary':128,'chunks':len(batches)}
    if configfile.exists():
        if load_json(configfile)!=config:raise ValueError('Frozen campaign run configuration changed')
    else:write_json(configfile,config)
    original_completed=sum((chunk_dir(i)/'complete.json').is_file() for i in range(len(batches)))
    attempted=0;started=time.monotonic()
    print('Original completed chunks',original_completed,'/',len(batches),flush=True)
    for i,batch in enumerate(batches):
        d=chunk_dir(i)
        if d.exists():
            audit_chunk(i,batch,binary)
            continue
        if attempted>=max_chunks or time.monotonic()-started>=hours*3600:
            print('SAFE STOP:',original_completed+attempted,'/',len(batches),
                  'physics checkpoints complete',flush=True)
            return
        d.mkdir() # atomic; never overwrite an unfinished native run
        write_csv(d/'native_muon_crossings_input.csv',FIELDS,batch)
        log=d/'physics.log'
        print('REAL GARFIELD:',i,'10 independent parents x 4 geometries',flush=True)
        t0=time.monotonic()
        with log.open('x') as f:
            cp=subprocess.run([str(EXE),'--muon-cache-physics',
                               str(d/'native_muon_crossings_input.csv'),
                               str(d/'muon_reference')],stdout=f,stderr=subprocess.STDOUT,
                              env=garfield_runtime_env(),cwd=d)
        elapsed=time.monotonic()-t0
        if cp.returncode or 'PASS: immutable same-Heed primary/cluster cache' not in log.read_text(errors='replace'):
            raise RuntimeError('Native full Garfield incomplete: '+str(d)+'; raw files kept')
        # Verified by original full-real-muon checks, amended for *strict*
        # 30um anode / 75um FW in-metal zero-ion starts only.
        stats=isolated_validator()(d,batch)
        if stats['n_full_avalanche_trials']!=4*stats['n_heed_primaries']:
            raise ValueError('Lost one or more exact Heed primary avalanches')
        sampled=check_birth_reservoir(d)
        names=['native_muon_crossings_input.csv','physics.log',*DATA_FILES]
        record={'schema':'NA60_SPACETIME_10PARENT_SHA_CHUNK_V1','status':'complete',
                'chunk':i,'n_parents':len(batch),'physics':binary,
                'elapsed_wall_seconds':elapsed,
                'n_full_Garfield_avalanches':stats['n_full_avalanche_trials'],
                'n_Heed_primaries':stats['n_heed_primaries'],
                'n_weighted_birth_samples':sampled,
                'sha256_files':{name:sha(d/name) for name in names}}
        write_json(d/'complete.json',record)
        audit_chunk(i,batch,binary)
        attempted+=1
        print('PASS SHA:',i,'full avalanches',stats['n_full_avalanche_trials'],
              'saved birth points',sampled,'wall s',round(elapsed,2),flush=True)
    print('ALL SELECTED PARENTS SHA-VERIFIED',len(batches),'chunks',flush=True)


def read_completed():
    rr,m=check_input();b=expected_binary();all_batches=source_split(rr)
    finished=[];openchunk=None
    for i,rows in enumerate(all_batches):
        d=chunk_dir(i)
        if d.exists():
            if openchunk is not None:raise ValueError('Gap in completed sequential checkpoint list')
            audit_chunk(i,rows,b)
            finished.append(i)
        elif openchunk is None:openchunk=i
    return rr,m,finished


def status():
    rr,m,done=read_completed()
    print('STATUS total parents',len(rr),'verified blocks',len(done),
          'out of',len(rr)//GROUP,
          'dev-parent limit',m['dev_parents'],
          'reserved independent parents',m['reserved_holdout_parents'])
    if done:
        n=sum(load_json(chunk_dir(i)/'complete.json')['n_full_Garfield_avalanches'] for i in done)
        print('REAL completed microscopic AvalancheMicroscopic calls:',n)


def quantile_interval(x):
    return [round(float(np.quantile(x,.025)),6),round(float(np.quantile(x,.975)),6)]


def moments(x):
    q=np.array([t for t in x if math.isfinite(t)],float)
    return (float(q.mean()),float(q.std(ddof=0))) if len(q)>0 else (float('nan'),float('nan'))


def report(reps,rms_target,gain_target,paired_target,final_holdout=False):
    rr,m,done=read_completed()
    if len(done)<12:raise ValueError('Require at least 120 full physical parent events for first report')
    if final_holdout and len(done)!=len(rr)//GROUP:
        raise ValueError('Do not unblind final holdout before ALL preselected physics chunks are complete')
    n=min(len(done)*GROUP,m['dev_parents'])
    begin=m['dev_parents'] if final_holdout else 0
    stop=m['n_parents'] if final_holdout else n
    if stop<=begin:raise ValueError('Not enough development data')
    rows=[]
    for i in range(begin//GROUP,stop//GROUP):
        with (chunk_dir(i)/'muon_reference_crossing.csv').open(newline='') as f:
            rows.extend(csv.DictReader(f))
    source=rr[begin:stop]
    ix={id_of(x):j for j,x in enumerate(source)}
    raw={g:{'gain':np.full(len(source),np.nan),
            'hemi':np.full(len(source),np.nan),
            'TB':np.full(len(source),np.nan)} for g in GEOMS}
    collected={g:0 for g in GEOMS};primaries=0
    for h in rows:
        g=h['configuration'];j=ix[h['source|chunk|event|hit_index|track_id|station|det_id']]
        if g not in raw:raise ValueError('Unknown geometry')
        for metric,col in (('gain','G_eff_per_Heed_primary'),
                            ('hemi','A_ion_hemi'),('TB','A_ion_TB')):
            raw[g][metric][j]=float(h[col]) if h[col] else float('nan')
        collected[g]+=int(h['collected_avalanche_electrons'])
        if g==GEOMS[0]:primaries+=int(h['heed_primaries'])
    strata=defaultdict(list)
    for j,s in enumerate(source):
        k=(int(s['station']),-1 if float(s['theta_u_in_deg'])<0 else 1)
        strata[k].append(j)
    if len(strata)!=12:raise ValueError('Missing station/signed-angle validation stratum')
    # Reweighting is a SOURCE-STRATUM PROXY, not an exact unstratified flux:
    # a single source parent is selected once, while the source can generate
    # multiple correlated crossings at different stations.
    src_counts=m['eligible_source_counts']
    denom=sum(int(v) for v in src_counts.values())
    source_weights=np.zeros(len(source),float)
    for (station,sign),ii in strata.items():
        source_weights[ii]=int(src_counts[f'{station}_{sign}'])/denom/len(ii)
    if abs(float(source_weights.sum())-1)>1e-10:
        raise ValueError('Geant4 source stratum weight normalization failed')
    # Stable predeclared stratified parent-level bootstrap; 4 geometries of
    # each parent resampled together, never treat their electrons as independent.
    rng=np.random.default_rng(20261104)
    boots=[];point=np.arange(len(source))
    for _ in range(reps):
        sample=np.concatenate([rng.choice(z,size=len(z),replace=True)
                               for _,z in sorted(strata.items())])
        boots.append(sample)
    def value(g,metric,ii,stat,source_weighted=False):
        z=raw[g][metric][ii]
        ok=np.isfinite(z)
        z=z[ok]
        if len(z)==0:return float('nan')
        if not source_weighted:
            return float(z.mean() if stat=='mean' else z.std(ddof=0))
        w=source_weights[ii][ok]
        w=w/w.sum()
        mean=float(np.dot(w,z))
        return mean if stat=='mean' else float(np.sqrt(np.dot(w,(z-mean)**2)))
    result={'schema':'NA60_SPACETIME_PARENT_PRECISION_REPORT_V1',
            'population':'reserved_holdout' if final_holdout else 'development',
            'n_parent_events':len(source),'n_true_Garfield_avalanche_realizations':
             sum(load_json(chunk_dir(i)/'complete.json')['n_full_Garfield_avalanches']
                 for i in range(begin//GROUP,stop//GROUP)),
            'source':'12 equal-sized station x signed-angle strata, with optional Geant4-eligible-stratum weighted PROXY moments; not an exactly flux-random source sample',
            'source_proxy_weighting':'Original eligible Geant4 source crossing counts for 12 station/signed-angle bins; exact parent-level unstratified cross-section requires another sampling design',
            'source_proxy_effective_sample_size':round(1.0/float(np.sum(source_weights**2)),1),
            'bootstrap':'95pct percentile stratified parent resampling, four geometries paired',
            'n_bootstrap_replicates':reps,'per_geometry':{},'FW_minus_noFW':{},
            'thresholds':{'rms_half_width':rms_target,
                          'mean_gain_relative_half_width':gain_target,
                          'paired_sector_RMS_half_width':paired_target}}
    passes=[]
    for g in GEOMS:
        data={}
        for metric,stat in (('gain','mean'),('hemi','rms'),('TB','rms')):
            z=value(g,metric,point,'mean' if stat=='mean' else 'rms')
            vals=np.array([value(g,metric,b,'mean' if stat=='mean' else 'rms') for b in boots])
            lo,hi=quantile_interval(vals)
            hw=(hi-lo)/2
            threshold=gain_target*abs(z) if metric=='gain' else rms_target
            ready=math.isfinite(hw) and hw<=threshold
            passes.append(ready)
            data[metric]={'statistic':stat,'estimate':round(z,7),
                          'CI95':[lo,hi],'halfwidth':round(hw,7),
                          'threshold':round(threshold,7), 'criterion_met':ready,
                          'defined_count':int(np.isfinite(raw[g][metric]).sum())}
        data['pooled_electron_weighted_gain_balanced_sample']=collected[g]/primaries if primaries else None
        data['source_proxy_moments']={}
        for metric,stat in (('gain','mean'),('hemi','rms'),('TB','rms')):
            z=value(g,metric,point,'mean' if stat=='mean' else 'rms',True)
            vals=np.array([value(g,metric,b,'mean' if stat=='mean' else 'rms',True)
                           for b in boots])
            ci=quantile_interval(vals);hw=(ci[1]-ci[0])/2
            threshold=gain_target*abs(z) if metric=='gain' else rms_target
            passes.append(hw<=threshold)
            data['source_proxy_moments'][metric]={'estimate':round(z,7),
                'CI95':ci,'halfwidth':round(hw,7),
                'criterion_met':bool(hw<=threshold)}
        # Tail quantiles are exploratory diagnostics, not precision-stopping
        # metrics. In particular 99th-percentile gain needs many more parents.
        valid=raw[g]['gain'][np.isfinite(raw[g]['gain'])]
        data['gain_tail_quantiles_balanced']={f'q{q}':round(float(np.quantile(valid,int(q)/100)),5)
              for q in ('50','90','95','99')} if len(valid) else {}
        result['per_geometry'][g]=data
    for p in (3,4):
        no=f'p{p}_noFW';fw=f'p{p}_FW75';key=f'p{p}_FW75_minus_noFW_TB_RMS'
        observed=value(fw,'TB',point,'rms')-value(no,'TB',point,'rms')
        vals=[value(fw,'TB',b,'rms')-value(no,'TB',b,'rms') for b in boots]
        ci=quantile_interval(vals);hw=(ci[1]-ci[0])/2
        ready=hw<=paired_target
        passes.append(ready)
        result['FW_minus_noFW'][key]={'estimate':round(observed,7),
             'CI95':ci,'halfwidth':round(hw,7),'criterion_met':ready}
        proxy=value(fw,'TB',point,'rms',True)-value(no,'TB',point,'rms',True)
        proxvals=[value(fw,'TB',b,'rms',True)-value(no,'TB',b,'rms',True)
                  for b in boots]
        ci2=quantile_interval(proxvals);hw2=(ci2[1]-ci2[0])/2
        passes.append(hw2<=paired_target)
        result['FW_minus_noFW'][key]['source_proxy']={
            'estimate':round(proxy,7),'CI95':ci2,'halfwidth':round(hw2,7),
            'criterion_met':bool(hw2<=paired_target)}
    result['all_prespecified_precision_criteria_met']=bool(all(passes))
    result['scientific_readiness_note']=(
       'Only mean collected gain and ion birth sector asymmetry are checked. '
       'This does NOT validate near-wire rare tails, continuous 3D+time ion clouds, '
       'electron response, cathode Shockley-Ramo strip charge, or source generator systematics. '
       'Never tune LUT on reserved holdout if using it as independent validation.')
    REPORTS.mkdir(parents=True,exist_ok=True)
    filename=REPORTS/(('FINAL_reserved_' if final_holdout else 'development_')+f'{len(source):04d}_parents.json')
    if filename.exists():
        print('Existing versioned report preserved (not overwritten):',filename)
        print(json.dumps(result,indent=2))
        return
    write_json(filename,result)
    print(json.dumps(result,indent=2))
    print('PASS immutable independent-parent precision report:',filename)


def selftest():
    aa=[{'source':'Jpsi','chunk':2,'event':j,'station':s,
         'theta_u_in_deg':(-1 if q==-1 else 1)}
        for j,(s,q) in enumerate(((s,q) for s in range(6) for q in (-1,1)))]
    assert len({parent_key(z) for z in aa})==12
    assert all(k in GEOMS for k in ('p3_noFW','p3_FW75'))
    assert source_split(list(range(30)))==[list(range(10)),list(range(10,20)),list(range(20,30))]
    pool={(s,q):[(0,('x',),{'source':'Jpsi','chunk':'2','event':str(s*2+(q+1)//2)})]
          for s in range(6) for q in (-1,1)}
    selected,n=select(pool,1)
    assert len(selected)==12 and len({parent_key(x) for x in selected})==12
    print('PASS independent parent selection / 12 strata / chunking synthetic selftest; NO physics')


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('mode',choices=('selftest','prepare','preflight','run','status','report'))
    ap.add_argument('--n',type=int,default=2400,help='Freeze n=2400 upfront; run a few blocks at a time')
    ap.add_argument('--max-chunks',type=int,default=1,help='New batches per call, 10 parents each')
    ap.add_argument('--hours',type=float,default=8)
    ap.add_argument('--reps',type=int,default=2000)
    ap.add_argument('--rms-halfwidth',type=float,default=.03)
    ap.add_argument('--gain-rel-halfwidth',type=float,default=.10)
    ap.add_argument('--paired-halfwidth',type=float,default=.03)
    ap.add_argument('--final-holdout',action='store_true')
    a=ap.parse_args()
    if a.mode=='selftest':selftest()
    elif a.mode=='prepare':prepare(a.n)
    elif a.mode=='preflight':preflight()
    elif a.mode=='run':run(a.max_chunks,a.hours)
    elif a.mode=='status':status()
    else:report(a.reps,a.rms_halfwidth,a.gain_rel_halfwidth,
                a.paired_halfwidth,a.final_holdout)

if __name__=='__main__':main()
