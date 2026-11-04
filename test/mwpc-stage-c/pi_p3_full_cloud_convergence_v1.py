#!/usr/bin/env python3
"""NA60+ p3 ONLY: SHA-frozen full Garfield ion-birth cloud pilot and compression study.

Use original SHA-verified D500 *initial* electron positions, NEVER a previous
anode endpoint or a muon surrogate. A20 and A100 old sample indices are excluded.
Two physically distinct 3-mm geometries, separate analyses and decisions.

Modes: selftest | prepare | preflight | run | status | analyze.
No cathode induced signal, ion drift, Geant4 muon rerun, or old-file mutation.
A full-ion cloud here means all microscopic ion-IONIZATION-BIRTH callbacks,
not every subsequent ion trajectory.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import itertools
import json
import math
import os
import random
import shutil
import subprocess
import time
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

from pi_audit_avalanche_lut_reuse import BASE, CONFIGS, EXPORT_FIELDS, locations, sha
from pi_run_lut500_overnight import offset_for, STRIDE, source_maps
from pi_garfield_runtime_env import garfield_runtime_env

REPO=Path(__file__).resolve().parents[2]
ROOT=BASE/'pi_p3_full_cloud_convergence_v1'
INPUT=ROOT/'input'
CHUNKS=ROOT/'full_cloud_chunks'
CPP=REPO/'test/mwpc-avalanche/phaseB_native_muon_joint_cloud_lut_v1.cxx'
EXE=REPO/'test/mwpc-avalanche/build/mwpc_phase_b_joint_cloud_lut_v1'
GEOMS=('p3_noFW','p3_FW75')
GROUP=8
FIELDS=tuple(EXPORT_FIELDS)
SUFFIX=('_trials.csv','_birth_reservoir.csv','_timing.csv',
        '_moments.csv','_electron_endpoint_sample.csv')
CAPS=(32,128,512)
RNG_TAG='NA60-P3-CLOUD-UNSEEN-D500-20261104-V1'


def must(c,msg):
    if not c:raise ValueError(msg)
def digest(s):return hashlib.sha256(s.encode('utf-8')).digest()
def stable(*items):return int.from_bytes(digest('|'.join(map(str,items)))[:8],'big')
def write_json(path,item):
    path.write_text(json.dumps(item,sort_keys=True,indent=2,allow_nan=False)+'\n')
def read_json(path):return json.loads(path.read_text())
def csv_records(path,fields=None):
    with path.open(newline='') as f:
        rd=csv.DictReader(f)
        if fields is not None:must(tuple(rd.fieldnames or ())==tuple(fields),
                                  'Wrong CSV columns: '+str(path))
        yield from rd

def write_csv(path,rows,fields=FIELDS):
    with path.open('x',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields)
        w.writeheader();w.writerows(rows)

def geometry_distance(u,v,g):
    # p3 local cell [-1.5,1.5], nearest anode at local u=0.
    da=math.hypot(u,v)
    if g.endswith('noFW'):
        return da,math.inf
    return da,min(math.hypot(u-1.5,v),math.hypot(u+1.5,v))

def classify(u,v,g):
    da,df=geometry_distance(u,v,g)
    if da<=.015 or df<=.0375:
        return None # 30-um anode diameter, 75-um FW diameter: in metal
    if da<.35:return 'near_anode'
    if df<.35:return 'near_field'
    if min(da,df)<1.:return 'transition'
    return 'bulk'

def stratum(u,v,g):
    c=classify(u,v,g)
    return None if c is None else (c,'positive' if v>=0 else 'negative')

def sample_for(g,box):
    # Exactly the frozen D500 permutation. Pass >=100 excludes *all* original
    # A20 passes [0,20) and selected A100 extra passes [20,100).
    passage=100+stable(RNG_TAG,g,box)%400
    return (offset_for(g,box)+STRIDE*passage)%500

def selection_from_sources(per_geom):
    must(per_geom>0 and per_geom<=96 and per_geom%24==0,
         'Use 24, 48, 72 or 96 selected starts per geometry')
    source,idxsha=source_maps() # verifies original immutable D500 provenance
    groups={}
    counts={}
    sha_orig={}
    for g,pitch,fw in CONFIGS:
        if g not in GEOMS:continue
        _,raw=locations(500,pitch,fw)
        must(sha(raw)==source[g]['lut_csv_sha256'], 'D500 source SHA changed '+g)
        sha_orig[g]=sha(raw)
        pool=defaultdict(list);found=set()
        targets={b:sample_for(g,b) for b in range(900)}
        with raw.open(newline='') as f:
            reader=csv.DictReader(f)
            must(set(('box_id','sample','iu','iv','u_start_mm','v_start_mm',
                      'w_start_mm','random_seed','drift_ok','collected',
                      'field_wire_collected','wire_index','phi_hit_deg',
                      'drift_time_ns','last_step_energy_eV')).issubset(reader.fieldnames or []),
                 'Missing frozen D500 column '+g)
            for r in reader:
                b=int(r['box_id']);s=int(r['sample'])
                if b not in targets or s!=targets[b]:continue
                must(b not in found,'Duplicate D500 selected candidate '+g+'/'+str(b))
                found.add(b)
                iu,iv=int(r['iu']),int(r['iv'])
                u,v,w=float(r['u_start_mm']),float(r['v_start_mm']),float(r['w_start_mm'])
                must(b==iu*30+iv and abs(w)<1.e-8 and
                     -1.5+iu*.1-1.e-6<=u<=-1.5+(iu+1)*.1+1.e-6 and
                     -3.+iv*.2-1.e-6<=v<=-3.+(iv+1)*.2+1.e-6,
                     'D500 cell/position corrupted '+g+'/'+str(b))
                region=stratum(u,v,g)
                if region is None:continue
                item={
                    'lut_trials':500,'geometry':g,
                    'original_lut_sha256':source[g]['lut_csv_sha256'],
                    'box_id':b,'iu':iu,'iv':iv,'sample':s,
                    'random_seed':int(r['random_seed']),
                    'u_start_mm':r['u_start_mm'],'v_start_mm':r['v_start_mm'],
                    'w_start_mm':r['w_start_mm'],
                    'old_drift_ok':r['drift_ok'], 'old_collected':r['collected'],
                    'old_field_wire_collected':r['field_wire_collected'],
                    'old_wire_index':r['wire_index'],
                    'old_phi_hit_deg':r['phi_hit_deg'],
                    'old_drift_time_ns':r['drift_time_ns'],
                    'old_last_step_energy_eV':r['last_step_energy_eV']}
                must(set(item)==set(FIELDS), 'LUT exact source columns changed')
                pool[region].append(item)
        must(len(found)==900, 'Not 900 distinct D500 available source points '+g)
        regions=('near_anode','transition','bulk') if g.endswith('noFW') else (
                 'near_anode','near_field','transition','bulk')
        # Equal strata for each geometry; no mixing across g or hemisphere.
        labels=[(a,side) for a in regions for side in ('positive','negative')]
        each=per_geom//len(labels)
        must(each*len(labels)==per_geom, 'Uneven space/time sector allocation')
        selected={}
        counts[g]={}
        for label in labels:
            ordered=sorted(pool[label], key=lambda r:(
                stable(RNG_TAG,'region-order',g,label,r['box_id']),r['box_id']))
            counts[g][str(label)]={'available':len(ordered),'chosen':each}
            must(len(ordered)>=each,
                 'Not enough true D500 cells for '+str((g,label,each,len(ordered))))
            selected[label]=ordered[:each]
        # Round-robin sectors: early checkpoint already spans near-wire + bulk.
        ordered=[]
        for i in range(each):
            for label in labels:ordered.append((selected[label][i],label))
        must(len(ordered)==per_geom,'Selection allocation error '+g)
        groups[g]=ordered
    all_rows=[];context=[]
    for i in range(per_geom):
        for g in GEOMS:
            row,label=groups[g][i];all_rows.append(row)
            context.append({'geometry':g,'box_id':row['box_id'],
                            'sample':row['sample'],'regime':label[0],
                            'v_sign':label[1]})
    must(len(all_rows)==per_geom*2,'Wrong total selected D500 starts')
    must(len({(r['geometry'],r['box_id'],r['sample']) for r in all_rows})==len(all_rows),
         'Duplicate pilot start')
    return all_rows,context,sha_orig,idxsha,counts

def check_plan():
    p=INPUT/'exact_D500_p3_unseen.csv'
    m=read_json(INPUT/'selection_manifest.json')
    rows=list(csv_records(p,FIELDS))
    tags=read_json(INPUT/'strata.json')
    must(m.get('schema')=='NA60_P3_CLOUD_INPUT_V1' and
         m.get('input_sha256')==sha(p) and
         m.get('strata_sha256')==sha(INPUT/'strata.json') and
         m.get('n_selected')==len(rows)==len(tags) and
         m.get('n_selected')==2*m['per_geometry'] and
         m['per_geometry']%24==0 and len(rows)%GROUP==0,
         'Pilot selection input identity differs')
    keys=[(r['geometry'],int(r['box_id']),int(r['sample'])) for r in rows]
    must(len(keys)==len(set(keys)) and
         all(a in GEOMS and sample_for(a,b)==s for a,b,s in keys),
         'Pilot original D500 trial identity/independence changed')
    must(all((r['geometry'],int(r['box_id']),int(r['sample']))==
             (t['geometry'],int(t['box_id']),int(t['sample']))
             for r,t in zip(rows,tags)), 'Pilot stratum metadata identity changed')
    must(all(r['geometry'] in GEOMS for r in rows),'A p4 geometry entered p3 pilot')
    return rows,tags,m

def prepare(per_geometry):
    must(not ROOT.exists(),'Refuse to change existing frozen pilot: '+str(ROOT))
    selected,tags,shas,indexsha,counts=selection_from_sources(per_geometry)
    INPUT.mkdir(parents=True)
    path=INPUT/'exact_D500_p3_unseen.csv'
    write_csv(path,selected)
    write_json(INPUT/'strata.json',tags)
    write_json(INPUT/'selection_manifest.json',{
        'schema':'NA60_P3_CLOUD_INPUT_V1','per_geometry':per_geometry,
        'n_selected':len(selected),'original_frozen_D500_sha256':shas,
        'original_index_manifest_sha256':indexsha,
        'source_selection':'One distinct archived original D500 start per spatial cell, stable pass>=100; excludes all A20/A100 passes; independent physical microscopic collisions; NOT muon flux sample',
        'input_sha256':sha(path),'strata_sha256':sha(INPUT/'strata.json'),
        'region_counts':counts,'group_size':GROUP,'geometries':GEOMS,
        'B_T':0,'anode_V':1800,'cathode_gap_mm':3,
        'gas':'Ar70CO2_30_293.15K_760T'})
    check_plan()
    print('PASS immutable p3-only full-cloud pilot input',len(selected),
          'exact D500 starts; '+str(per_geometry)+' per geometry',flush=True)
    for g in GEOMS:print('STRATA',g,counts[g],flush=True)

def check_cpp():
    must(CPP.is_file() and EXE.is_file(),'Build new isolated mwpc_phase_b_joint_cloud_lut_v1 C++ target')
    source=CPP.read_text()
    must('--lut-sidecar-full-births' in source and
         'std::numeric_limits<std::size_t>::max()' in source and
         '(*exactMomentsOut)' in source,
         'Full-birth target source does not match documented capabilities')
    return {'cpp_sha256':sha(CPP),'binary_sha256':sha(EXE)}

def preflight():
    rows,_,m=check_plan();binary=check_cpp()
    first=INPUT/'preflight8.csv'
    subset=rows[:GROUP]
    if not first.exists():write_csv(first,subset)
    must(list(csv_records(first,FIELDS))==subset,'Preflight input changed')
    proc=subprocess.run([str(EXE),'--validate-lut-input',str(first)],
                        capture_output=True,text=True,env=garfield_runtime_env())
    must(proc.returncode==0 and 'PASS parsed 8 ORIGINAL frozen LUT20/500 exact-start rows' in proc.stdout,
         'Native C++ strict D500 parser failed: '+proc.stdout+'\n'+proc.stderr)
    meta={'schema':'NA60_P3_CLOUD_PREFLIGHT_V1','binary':binary,
          'input_sha256':m['input_sha256'],'native_parser':'PASS 8 exact D500 sites',
          'n_selected':len(rows)}
    target=INPUT/'preflight.json'
    if target.exists():must(read_json(target)==meta,'Preflight physics identity changed')
    else:write_json(target,meta)
    print('PASS real C++ parser 8 unique unused D500 starts and frozen binary',flush=True)

def verified_binary():
    pre=read_json(INPUT/'preflight.json')
    rows,_,m=check_plan();binary=check_cpp()
    must(pre['binary']==binary and pre['input_sha256']==m['input_sha256'] and
         pre['n_selected']==len(rows),'Physics binary changed since preflight')
    return binary

def chunk_path(i):return CHUNKS/f'chunk_{i:04d}'
def prefix(d):return d/'cloud'
def batch(rows):return [rows[i:i+GROUP] for i in range(0,len(rows),GROUP)]

def check_outputs(d, expected, fullscan=True):
    # No acceptance on old deterministic drift flags: real avalanche may differ.
    p=prefix(d)
    trials=list(csv_records(Path(str(p)+SUFFIX[0])))
    moments=list(csv_records(Path(str(p)+SUFFIX[3])))
    times=list(csv_records(Path(str(p)+SUFFIX[2])))
    must(len(trials)==len(moments)==len(times)==len(expected),
         'Missing trial/timing/exact moment in '+str(d))
    expect={(g,int(r['box_id']),int(r['sample'])):r for r in expected
            for g in [r['geometry']]}
    seen=set();count_by_seed={}
    for x in trials:
        g=x['configuration']; b=int(x['box_id']);s=int(x['sample'])
        k=(g,b,s)
        must(k in expect and k not in seen and x['original_geometry']==g and
             x['original_lut_sha256']==expect[k]['original_lut_sha256'] and
             int(x['original_random_seed'])==int(expect[k]['random_seed']) and
             x['lut_trials']=='500','Mismatched original D500 site in Garfield trial')
        must(int(x['capped'])==0 and x['successful']=='1',
             'Incomplete/capped microscopic avalanche; preserve raw chunk: '+str(k))
        must(float(x['seed_w_mm'])==0. and
             abs(float(x['seed_u_mm'])-float(expect[k]['u_start_mm']))<1e-6 and
             abs(float(x['seed_v_mm'])-float(expect[k]['v_start_mm']))<1e-6,
             'Initial electron coordinate corrupted '+str(k))
        count=int(x['ion_births'])
        must(count>=0 and int(x['garfield_ni'])>=0,
             'Invalid microscopic ion count '+str(k))
        seen.add(k)
        # C++ sidecar has primary_index=0, crossing_id is unique LUT identity.
        sk=(x['crossing_id'],g,int(x['primary_index']))
        must(int(x['primary_index'])==0 and sk not in count_by_seed,
             'Malformed avalanche trial label')
        count_by_seed[sk]=count
    must(len(seen)==len(expected),'Missing physical trial')
    gotmom=set()
    for r in moments:
        key=(r['crossing_id'],r['configuration'],int(r['primary_index']))
        must(key in count_by_seed and key not in gotmom and
             int(r['ion_births'])==count_by_seed[key],
             'Exact birth-moment label/count mismatch')
        gotmom.add(key)
    must(len(gotmom)==len(count_by_seed),'Incomplete exact birth moments')
    must(len({(r['crossing_id'],r['configuration'],int(r['primary_index']))
              for r in times})==len(expected),'Timing identities do not match')
    samples=Counter();total=0
    if fullscan:
        for r in csv_records(Path(str(p)+SUFFIX[1])):
            key=(r['crossing_id'],r['configuration'],int(r['primary_index']))
            must(key in count_by_seed,'Unknown ion birth source '+str(key))
            must(int(r['seed_ion_births_total'])==count_by_seed[key] and
                 abs(float(r['reservoir_weight'])-1.)<1e-9,
                 'Full cloud must have weight 1 and complete ion count')
            must(all(math.isfinite(float(r[c])) for c in ('u_mm','v_mm','w_mm','t_ns')),
                 'Nonfinite ion coordinate')
            samples[key]+=1;total+=1
        must(all(samples[key]==n for key,n in count_by_seed.items()),
             'Full microscopic ion cloud truncated or incomplete')
    return {'n_avalanches':len(expected),'n_full_ion_births':total,
            'n_zero_avalanches':sum(n==0 for n in count_by_seed.values()),
            'n_nonzero_avalanches':sum(n>0 for n in count_by_seed.values())}

def audit_chunk(i,expected,binary):
    d=chunk_path(i);seal=d/'complete.json'
    must(seal.is_file(),'Unsealed/incomplete chunk must be inspected: '+str(d))
    rec=read_json(seal)
    must(rec.get('schema')=='NA60_P3_FULL_CLOUD_SHA_CHUNK_V1' and
         rec.get('chunk')==i and rec.get('n_avalanches')==len(expected) and
         rec.get('binary')==binary,'Chunk identity differs')
    for name,digest in rec['sha256_files'].items():
        must(sha(d/name)==digest,'Frozen chunk checksum differs: '+str(d/name))
    must(list(csv_records(d/'exact_D500_chunk.csv',FIELDS))==expected,
         'Completed chunk original D500 source changed')
    return rec

def run(max_chunks,hours):
    rows,_,m=check_plan(); binary=verified_binary()
    batches=batch(rows)
    must(max_chunks>=1 and hours>0,'Positive chunk/time budgets required')
    CHUNKS.mkdir(exist_ok=True)
    existing=0;attempt=0;start=time.monotonic()
    for i,subset in enumerate(batches):
        d=chunk_path(i)
        if d.exists():
            audit_chunk(i,subset,binary)
            existing+=1
            continue
        if attempt>=max_chunks or time.monotonic()-start>=hours*3600:
            print('SAFE STOP',existing+attempt,'/',len(batches),flush=True)
            return
        # Full clouds may produce millions of rows; reserve disk space.
        space=shutil.disk_usage(ROOT).free
        must(space>2*1024**3,'Less than 2 GiB remaining: full birth cloud requires room')
        d.mkdir() # NEVER reuse an interrupted folder automatically.
        write_csv(d/'exact_D500_chunk.csv',subset)
        logfile=d/'physics.log'
        print('FULL MICROSCOPIC GARFIELD p3 ONLY chunk',i,
              'avalanches',len(subset),flush=True)
        start_chunk=time.monotonic()
        with logfile.open('x') as f:
            cp=subprocess.run([str(EXE),'--lut-sidecar-full-births',
                               str(d/'exact_D500_chunk.csv'),str(prefix(d))],
                               cwd=d,stdout=f,stderr=subprocess.STDOUT,
                               env=garfield_runtime_env())
        elapsed=time.monotonic()-start_chunk
        must(cp.returncode==0 and
             'PASS: independent full-start microscopic avalanche sidecar' in
             logfile.read_text(errors='replace'),
             'Native Garfield failed/partial data preserved for inspection '+str(d))
        record=check_outputs(d,subset,fullscan=True)
        paths=['exact_D500_chunk.csv','physics.log']+[prefix(d).name+s for s in SUFFIX]
        seal={'schema':'NA60_P3_FULL_CLOUD_SHA_CHUNK_V1','chunk':i,
              'n_avalanches':len(subset),'binary':binary,
              'n_ion_births':record['n_full_ion_births'],
              'zero_avalanches':record['n_zero_avalanches'],
              'elapsed_wall_seconds':elapsed,
              'sha256_files':{name:sha(d/name) for name in paths}}
        write_json(d/'complete.json',seal)
        audit_chunk(i,subset,binary)
        attempt+=1
        print('PASS SHA FULL BIRTHS',i,'ions',record['n_full_ion_births'],
              'zero',record['n_zero_avalanches'],'wall_s',round(elapsed,1),flush=True)
    print('PASS ALL COMPLETE',len(rows),'independent exact-start p3 avalanches',flush=True)

def completed():
    rows,tags,meta=check_plan()
    binary=verified_binary()
    batches=batch(rows)
    done=[];missing=False
    for i,x in enumerate(batches):
        d=chunk_path(i)
        if d.exists():
            must(not missing,'Gap in SHA chunks')
            done.append(audit_chunk(i,x,binary))
        else:missing=True
    return rows,tags,meta,done

def status():
    rows,tags,meta,done=completed()
    print('STATUS p3 only, selected',len(rows),'parents (D500 electrons, NOT muons)')
    print('VERIFIED CHUNKS',len(done),'/',len(rows)//GROUP)
    print('FULL MICROSCOPIC AVA',sum(x['n_avalanches'] for x in done))
    print('EXACT BIRTH POINTS',sum(x['n_ion_births'] for x in done))
    print('p3_noFW / p3_FW75 only; no cathode signals')


def wrapped_u(a):
    # Match std::round away from zero in C++ near nearest anode.
    q=a/3.
    nearest=np.where(q>=0,np.floor(q+.5),np.ceil(q-.5))
    return a-3*nearest

def cloud_stats(points):
    # Physical phase of every birth relative to nearest REAL anode.
    u,v,w,t=points.T
    du=wrapped_u(u);radius=np.hypot(du,v)
    phi=np.arctan2(v,du)
    above=np.count_nonzero((phi>0)&(phi<math.pi))
    below=np.count_nonzero((phi<0)&(phi>-math.pi))
    top=np.count_nonzero((phi>=math.pi/4)&(phi<3*math.pi/4))
    bottom=np.count_nonzero((phi>=-3*math.pi/4)&(phi<-math.pi/4))
    n=len(points)
    return {'up':above,'down':below,'top':top,'bottom':bottom,
            'top_frac':top/n,'bottom_frac':bottom/n,
            'A_TB':(top-bottom)/(top+bottom) if top+bottom else math.nan,
            'A_hemi':(above-below)/(above+below) if above+below else math.nan,
            'mean_radius_mm':float(np.mean(radius)),
            'sd_radius_mm':float(np.std(radius)),
            'mean_t_ns':float(np.mean(t)),
            'sd_t_ns':float(np.std(t)),
            'q90_t_ns':float(np.quantile(t,.9)),
            'q90_radius_mm':float(np.quantile(radius,.9)),
            'mean_w_mm':float(np.mean(w)),
            'corr_radius_time':float(np.corrcoef(radius,t)[0,1]) if
             n>2 and np.std(radius)>1.e-12 and np.std(t)>1.e-12 else math.nan,
            'phi':phi,'t':t}

def joint_hist(data,min_t,max_t):
    # 12 angle x 6 time bins. Keep edges fixed *from exact cloud*.
    span=max_t-min_t
    norm=(data['t']-min_t)/span if span>1.e-12 else np.zeros(len(data['t']))
    v=np.histogram2d(data['phi'],np.clip(norm,0,1),
                     bins=(np.linspace(-math.pi,math.pi,13),
                           np.linspace(0,1+1.e-9,7)))[0]
    must(v.sum()>0,'Empty normalized joint distribution')
    return v/v.sum()

def compression_errors(full, sample, f=None, f_joint=None):
    if f is None:f=cloud_stats(full)
    s=cloud_stats(sample)
    time_min,time_max=float(np.min(f['t'])),float(np.max(f['t']))
    if f_joint is None:f_joint=joint_hist(f,time_min,time_max)
    s_joint=joint_hist(s,time_min,time_max)
    fields=('mean_radius_mm','sd_radius_mm','mean_t_ns','sd_t_ns',
            'q90_t_ns','q90_radius_mm','mean_w_mm','top_frac','bottom_frac',
            'A_TB','A_hemi','corr_radius_time')
    out={f'abs_{a}':abs(f[a]-s[a]) if math.isfinite(f[a]) and
         math.isfinite(s[a]) else None for a in fields}
    out['joint_phi_time_TV']=float(.5*np.sum(np.abs(f_joint-s_joint)))
    out['mean_radius_rel_abs']=(abs(f['mean_radius_mm']-s['mean_radius_mm'])/
        f['mean_radius_mm'] if f['mean_radius_mm']>1.e-12 else None)
    out['mean_t_abs_ns']=abs(f['mean_t_ns']-s['mean_t_ns'])
    return out

def verify_moments(full,exact,trial):
    n=len(full)
    must(n==int(exact['ion_births'])==int(trial['ion_births']),
         'Exact cloud/trial/moments count mismatch')
    if not n:return
    v=cloud_stats(full)
    required=(('mean_birth_radius_mm','mean_radius_mm'),
              ('sd_birth_radius_mm','sd_radius_mm'),
              ('mean_birth_time_ns','mean_t_ns'),
              ('sd_birth_time_ns','sd_t_ns'),
              ('mean_birth_w_mm','mean_w_mm'))
    for col,name in required:
        a,b=float(exact[col]),v[name]
        must(math.isclose(a,b,rel_tol=1.e-5,abs_tol=1.e-5),
             'Full cloud does not reproduce exact '+col)
    for col,metric in (('ion_up','up'),('ion_down','down'),
                       ('ion_top','top'),('ion_bottom','bottom')):
        must(int(trial[col])==v[metric],
             'Full cloud does not reproduce exact all-ion angular sector '+col)

def trial_key(row):return (row['crossing_id'],row['configuration'],int(row['primary_index']))
def full_cloud_iterator(path):
    # C++ RunLutSidecar writes *one whole avalanche at a time*. Stream
    # groupby, avoiding loading all (potentially >million) full-ion points.
    with path.open(newline='') as f:
        reader=csv.DictReader(f)
        for key,rs in itertools.groupby(reader,trial_key):
            data=[]
            for r in rs:
                must(float(r['reservoir_weight'])==1.,'Not a genuine complete cloud')
                data.append([float(r[col]) for col in ('u_mm','v_mm','w_mm','t_ns')])
            yield key,np.asarray(data,dtype=np.float64)

def quantile(v,p):return float(np.quantile(v,p)) if len(v) else None

def analyze(reps,partial=False):
    must(8<=reps<=500,'Require 8..500 subsampling reps')
    rows,labels,m,done=completed()
    must(done,'No completed native full-cloud physics chunks')
    must(partial or len(done)*GROUP==len(rows),
         'Incomplete; analyze with --partial for diagnostic only, or finish physics')
    n=len(done)*GROUP
    target=ROOT/('analysis_partial_%03d'%len(done) if partial else 'analysis_complete')
    must(not target.exists(),'Preserve existing analysis directory: '+str(target))
    truth=[];seedplan={}
    for r,t in zip(rows[:n],labels[:n]):
        k=(r['geometry'],int(r['box_id']),int(r['sample']))
        seedplan[k]=t
    for i in range(len(done)):
        d=chunk_path(i);p=prefix(d)
        trials=list(csv_records(Path(str(p)+SUFFIX[0])))
        moments=list(csv_records(Path(str(p)+SUFFIX[3])))
        must(len(trials)==len(moments)==GROUP,'Full cloud native row mismatch')
        bytrial={trial_key(r):r for r in trials}
        bymom={trial_key(r):r for r in moments}
        expected=set(bytrial)
        must(expected==set(bymom),'Different trial/moment identities')
        processed=set()
        for key,full in full_cloud_iterator(Path(str(p)+SUFFIX[1])):
            must(key in bytrial and key not in processed,'Unrecognized/repeated ion source')
            trial=bytrial[key];exact=bymom[key]
            if int(trial['ion_births'])==0:
                must(False,'A zero ion trial unexpectedly appears in nonempty CSV')
            verify_moments(full,exact,trial)
            processed.add(key)
            g=trial['configuration'];b=int(trial['box_id']);s=int(trial['sample'])
            label=seedplan[g,b,s]
            fullstats=cloud_stats(full)
            fulljoint=joint_hist(fullstats,float(np.min(fullstats['t'])),
                                 float(np.max(fullstats['t'])))
            entry={'geometry':g,'box_id':b,'sample':s,
                   'regime':label['regime'],'v_sign':label['v_sign'],
                   'n_ion_births':len(full),'exact_top_frac':fullstats['top_frac'],
                   'exact_A_TB':fullstats['A_TB'],
                   'cap_study':{}}
            rg=random.Random(stable(RNG_TAG,'subsample',g,b,s))
            for cap in CAPS:
                errors=[]
                for rep in range(reps):
                    k=min(cap,len(full))
                    indices=rg.sample(range(len(full)),k)
                    sampled=full[indices,:]
                    errors.append(compression_errors(full,sampled,fullstats,fulljoint))
                names=errors[0]
                entry['cap_study'][str(cap)]={name:{
                    'mean':float(np.mean(vals)),'q90':quantile(vals,.9),
                    'maximum':float(max(vals))} if vals else None
                    for name in names for vals in
                      [[e[name] for e in errors if e[name] is not None]]}
            truth.append(entry)
        for key in expected-processed:
            trial=bytrial[key]
            must(int(trial['ion_births'])==0,
                 'Full cloud missing for nonzero avalanche '+str(key))
            g,b,s=trial['configuration'],int(trial['box_id']),int(trial['sample'])
            label=seedplan[g,b,s]
            truth.append({'geometry':g,'box_id':b,'sample':s,
                'regime':label['regime'],'v_sign':label['v_sign'],
                'n_ion_births':0,'exact_top_frac':None,'exact_A_TB':None,
                'cap_study':{}})
    must(len(truth)==n and len({(x['geometry'],x['box_id'],x['sample']) for x in truth})==n,
         'Final analysis missing full microscopic avalanches')
    summaries={}
    for g in GEOMS:
        allrows=[r for r in truth if r['geometry']==g]
        must(len(allrows)==n//2,'Independent geometry cardinality')
        s={'n':len(allrows),'n_nonzero':sum(r['n_ion_births']>0 for r in allrows),
           'n_zero':sum(r['n_ion_births']==0 for r in allrows),
           'total_exact_ion_births':sum(r['n_ion_births'] for r in allrows),
           'by_regime':{}}
        for regime in ('near_anode','near_field','transition','bulk'):
            subset=[r for r in allrows if r['regime']==regime]
            if not subset:continue
            region={'n':len(subset),'n_nonzero':sum(r['n_ion_births']>0 for r in subset)}
            for cap in CAPS:
                val=[r['cap_study'][str(cap)] for r in subset if r['n_ion_births']>0]
                region[str(cap)]={
                    key:{'mean_across_avalanches':float(np.mean([a[key]['mean'] for a in val
                                                 if a[key] is not None])),
                         'q90_across_avalanches':float(np.quantile([a[key]['mean'] for a in val
                                                 if a[key] is not None],.9))}
                    for key in ('joint_phi_time_TV','abs_A_TB','abs_top_frac',
                                'mean_radius_rel_abs','abs_q90_t_ns')
                    if any(a[key] is not None for a in val)
                }
            s['by_regime'][regime]=region
        summaries[g]=s
    target.mkdir(parents=True)
    outputs={'summary.json':{
                'schema':'NA60_P3_FULL_BIRTH_COMPRESS_ANALYSIS_V1',
                'selected_starts':n,'completed_chunks':len(done),
                'reps_per_cap_per_avalanche':reps,
                'warning':'These are deliberately stratified D500 starts, NOT a muon sample and NOT a cathode signal validation.',
                'no_automatic_pass':'All benchmark thresholds must be fixed before downstream signal analysis.',
                'summary_by_geometry':summaries},
             'per_avalanche.json':truth}
    for name,obj in outputs.items():write_json(target/name,obj)
    write_json(target/'manifest.json',{
        'schema':'NA60_P3_FULL_CLOUD_ANALYSIS_SHA_V1',
        'selection_sha256':m['input_sha256'],
        'input_chunk_sha256':{f'chunk_{i:04d}/complete.json':sha(chunk_path(i)/'complete.json') for i in range(len(done))},
        'output_sha256':{name:sha(target/name) for name in outputs}})
    for g in GEOMS:
        print('EXACT CLOUD',g,'full births',summaries[g]['total_exact_ion_births'],
              'nonzero avalanches',summaries[g]['n_nonzero'],
              'zero',summaries[g]['n_zero'])
        for regime,s in summaries[g]['by_regime'].items():
            key='joint_phi_time_TV'; k='128'
            val=s[k].get(key)
            print('REGIME',g,regime,'n',s['n'],
                  '128-point joint φ/time TV mean',round(val['mean_across_avalanches'],4) if val else 'n/a')
    print('PASS SHA exact full birth clouds, no Garfield response simulation,',target,flush=True)


def selftest():
    must(classify(0.,0.,'p3_noFW') is None,'Anode metal not recognized')
    must(classify(1.5,0.,'p3_FW75') is None,'FW metal not recognized')
    must(classify(.2,.01,'p3_noFW')=='near_anode','Near anode classification')
    must(classify(1.47,.2,'p3_FW75')=='near_field','Field wire classification')
    must(classify(1.,2.0,'p3_noFW')=='bulk','Bulk classification')
    for g in GEOMS:
        for b in (0,1,10,899):
            old={(offset_for(g,b)+STRIDE*k)%500 for k in range(100)}
            must(sample_for(g,b) not in old,'Pilot reused A20/A100 sample')
    toy=np.asarray([[0.,.10,0.,0.],[.05,.15,.01,.2],
                    [-.10,-.20,.02,.4],[.12,.09,.03,.9]],float)
    stats=cloud_stats(toy)
    must(stats['up']==3 and stats['down']==1 and stats['top']==2,
         'Toy azimuth wrong')
    vals=compression_errors(toy,toy)
    must(max(abs(x) for x in vals.values() if x is not None)<1.e-10,
         'Self-consistent cloud compression nonzero')
    must(len(list(itertools.groupby([1,1,2,2,3])))==3,'Streaming test')
    print('PASS p3 2-geometry input strata, D500 >=100 disjointness, exact-cloud measures; NO Garfield',flush=True)


def main():
    a=argparse.ArgumentParser(description=__doc__)
    a.add_argument('mode',choices=('selftest','prepare','preflight','run','status','analyze'))
    a.add_argument('--per-geometry',type=int,default=24)
    a.add_argument('--max-chunks',type=int,default=1)
    a.add_argument('--hours',type=float,default=6.)
    a.add_argument('--reps',type=int,default=64)
    a.add_argument('--partial',action='store_true')
    z=a.parse_args()
    if z.mode=='selftest':selftest()
    elif z.mode=='prepare':prepare(z.per_geometry)
    elif z.mode=='preflight':preflight()
    elif z.mode=='run':run(z.max_chunks,z.hours)
    elif z.mode=='status':status()
    else:analyze(z.reps,z.partial)
if __name__=='__main__':main()
