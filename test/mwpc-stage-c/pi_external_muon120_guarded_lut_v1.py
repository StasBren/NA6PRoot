#!/usr/bin/env python3
"""Freeze and run a NEW 120-parent Geant4 archive holdout, then blind-compare
A20/A100 vs fixed guarded-neighbor A20/A100 with shared-cell bootstrap.

No new Geant4 is generated: 120 UNSEEN original J/psi chunk-2 parent events
are selected, disjoint from the 120 original chunk-2 muon parent events.
TrackHeed and full Garfield AvalancheMicroscopic WILL be calculated anew
on the fixed, cached original Geant4 crossings. Explicit phases, SHA seals,
10-parent checkpoints, no overwriting original LUT/holdout, no implicit upload.
"""
from __future__ import annotations
import argparse
import csv
import json
import math
import random
import time
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pi_holdout120_A20_from_D500 as legacy
from pi_2m_native_six_geometry_fold import crossing_provenance, SOURCE
from pi_prepare_muon_avalanche_b1 import FIELDS, validate, sha, seed
from pi_spatial_lut_muon120_crn_v1 import GEOMS, load_native as load_lut_native, donors_for, compile_neighbors
from pi_A100_crn_hierarchical_check import bin_for
from pi_spatial_lut_muon120_hier_v1 import simulate, metrics, report, MODES

REPO = Path(__file__).resolve().parents[2]
BASE = REPO/'test_runs/mwpc_stage_c'
ROOT = BASE/'pi_external_new120_guarded_lut_blind_v1'
INPUT = ROOT/'input_unseen_chunk2_120'
REF = ROOT/'new_full_garfield_120'
EVAL = ROOT/'blind_evaluation_fixed_guard_v1'
OLD_INPUT = legacy.INPUT
OLD_REF = legacy.REF
SOURCE_CHUNK = 2
N_PARENT = 120
LUT_CODE = REPO/'test/mwpc-stage-c/pi_spatial_lut_muon120_crn_v1.py'
HIER_CODE = REPO/'test/mwpc-stage-c/pi_spatial_lut_muon120_hier_v1.py'
OLD_RUNNER_CODE = REPO/'test/mwpc-stage-c/pi_holdout120_A20_from_D500.py'
PARAMS = {'radius_bins':2,'neighbor_probability':0.75,'local_probability':0.25,
          'neighbor_kernel':'exp(-(du**2+dv**2)/2)',
          'guard':'same v hemisphere; near-axis u sign; no donor/receiver within 0.35mm of wire centre',
          'training_bootstrap':'resample 20 original A20 per spatial cell ONCE globally/replicate, plus original 80 new A100 cells',
          'analysis_protocol':'90pct predictive interval of ion births; MAE of predictive mean; four correlated geometries per physical parent; no parameter fitting on this new sample'}

# Delegate the expensive original, rigorously checked physical C++ pipeline to
# the frozen 927-line legacy runner, replacing ONLY its module-level PATHS.
# No calls to legacy.prepare() or legacy.evaluate() which belong to the old 120.
legacy.ROOT=ROOT
legacy.INPUT=INPUT
legacy.REF=REF
legacy.COMP=ROOT/'unused_legacy_A20_comparison'


def old_parents():
    p=OLD_INPUT/'native_muon_crossings_input.csv'
    m=OLD_INPUT/'input_manifest.json'
    if not p.is_file() or not m.is_file():
        raise FileNotFoundError('Previously used original 120 source input missing')
    original=json.loads(m.read_text())
    if original.get('status')!='complete' or original.get('n_muons')!=120 or original.get('sha256_input')!=sha(p):
        raise ValueError('Old 120 input provenance mismatch')
    rr=legacy.read_csv(p,FIELDS)
    if len(rr)!=120:raise ValueError('Original 120 count changed')
    return {(x['source'],str(x['chunk']),str(x['event'])) for x in rr}, sha(p)


def fingerprint():
    return {'fixed_guarded_lut_source_sha256':sha(LUT_CODE),
            'fixed_hier_bootstrap_source_sha256':sha(HIER_CODE),
            'original_complete_full_garfield_driver_sha256':sha(OLD_RUNNER_CODE)}


def prepare():
    if ROOT.exists():raise FileExistsError('Separate new external holdout already exists: '+str(ROOT))
    forbidden, old_sha=old_parents()
    path,proof=crossing_provenance('Jpsi',SOURCE_CHUNK)
    source=SOURCE['Jpsi']
    if not proof.get('csv_sha256') or sha(path)!=proof['csv_sha256']:
        raise ValueError('Original Geant4 chunk-2 SHA failed')
    # Select without looking at Garfield, LUT predictions, ion births, or truth.
    # Use 160 best candidate rows per station/signed-angle stratum, excluding
    # ALL previous 120 parent events (not merely previously used track IDs).
    pools=defaultdict(list);counts=Counter();rejected=0
    with path.open(newline='') as f:
        rd=csv.DictReader(f)
        required=set(FIELDS)-{'source','chunk','momentum_GeV','heed_seed'}
        if not required.issubset(set(rd.fieldnames or ())):
            raise ValueError('Frozen Geant4 crossing source columns missing')
        for r in rd:
            if (source,str(SOURCE_CHUNK),str(r['event'])) in forbidden:
                rejected+=1
                continue
            checked=validate(r)
            if checked is None:continue
            mom,vals=checked
            if abs(vals['u_out_mm']-vals['u_in_mm'])>8 or abs(vals['w_out_mm']-vals['w_in_mm'])>40:
                continue
            station=int(r['station'])
            if not 0<=station<6:continue
            sign= -1 if float(r['theta_u_in_deg'])<0 else 1
            stratum=(station,sign);counts[stratum]+=1
            src_id=(source,SOURCE_CHUNK,r['event'],r['hit_index'],r['track_id'],station,r['det_id'])
            # New independent fixed SHA ranking, different salt from old120.
            ranking=seed('EXTERNAL-UNSEEN-MUON120-FIXED-GUARD-20261031',*src_id)
            rec={k:r[k] for k in FIELDS if k not in ('source','chunk','momentum_GeV','heed_seed')}
            rec['source']=source;rec['chunk']=SOURCE_CHUNK
            rec['momentum_GeV']=f'{mom:.12g}'
            rec['heed_seed']=seed('EXTERNAL-HEED-NEW-MUON-20261031',*src_id)
            arr=pools[stratum];arr.append((ranking,tuple(map(str,src_id)),rec))
            arr.sort(key=lambda x:(x[0],x[1]))
            if len(arr)>160:arr.pop()
    selected,by_stratum=legacy.choose_holdout_from_pool(pools,per_stratum=10)
    original_parent={(x['source'],str(x['chunk']),str(x['event'])) for x in selected}
    if len(original_parent)!=N_PARENT or original_parent&forbidden:
        raise ValueError('New parents overlap with old 120 Geant4 parent events')
    if len({legacy.id_of(x) for x in selected})!=120:
        raise ValueError('New parent track identity non-unique')
    INPUT.mkdir(parents=True)
    inp=INPUT/'native_muon_crossings_input.csv'
    legacy.write_csv(inp,FIELDS,selected)
    manifest={'schema':'PI_INDEPENDENT_GEANT4_MUON_HOLDOUT120_V1',
              'new_campaign_schema':'NA60_EXTERNAL_120_OLD_CHUNK2_UNSEEN_PARENT_FIXED_GUARD_V1',
              'status':'complete','n_muons':120,'n_distinct_parent_events':120,
              'source_chunk':SOURCE_CHUNK,'source_id':source,
              'sha256_input':sha(inp),
              'original_old120_input_sha256':old_sha,
              'frozen_source_sha256':proof['csv_sha256'],
              'original_root_sha256':proof['immutable_input']['input_root_sha256'],
              'selection':'new SHA ranking, 10 per each of 12 signed-angle/station strata; 160 candidate rows; exclude ALL original 120 parent event IDs before selection; original 9-wire valid-domain cuts',
              'original_lut_A20_A100_immutable':True,
              'models_fixed_before_seeing_new_full_Garfield':PARAMS,
              'source_code_seals_before_new_reference':fingerprint(),
              'old_parent_excluded_count':120,
              'candidate_rows_rejected_from_old_parents':rejected,
              'n_candidates_by_stratum':{str(k):counts[k] for k in sorted(counts)},
              'original_selected_parent_keys_sha256':sha(inp)}
    legacy.write_json(INPUT/'input_manifest.json',manifest)
    check()
    print('PASS PREPARED 120 NEW previously unused original Geant4 parent events;',
          'disjoint from old 120, 10 per station/signed-angle stratum.',flush=True)
    print('Frozen protocol/params/checksums:',INPUT/'input_manifest.json',flush=True)


def check():
    old,old_sha=old_parents()
    rows,manifest=legacy.payload_rows()  # old validator, redirected to NEW input only
    if manifest.get('new_campaign_schema')!='NA60_EXTERNAL_120_OLD_CHUNK2_UNSEEN_PARENT_FIXED_GUARD_V1':
        raise ValueError('This is not the frozen NEW external campaign')
    if manifest.get('original_old120_input_sha256')!=old_sha or manifest.get('models_fixed_before_seeing_new_full_Garfield')!=PARAMS:
        raise ValueError('Old120 provenance or frozen hyperparameters changed')
    if manifest.get('source_code_seals_before_new_reference')!=fingerprint():
        raise ValueError('Frozen analysis or physics source Python changed; halt before new reference use')
    original_path,proof=crossing_provenance('Jpsi',SOURCE_CHUNK)
    if (manifest.get('frozen_source_sha256')!=proof['csv_sha256'] or
        manifest.get('original_root_sha256')!=proof['immutable_input']['input_root_sha256']):
        raise ValueError('Original Geant4 source changed')
    now={(r['source'],str(r['chunk']),str(r['event'])) for r in rows}
    if len(now)!=120 or now&old or len({legacy.id_of(r) for r in rows})!=120:
        raise ValueError('Overlap with original parents or duplicate new muons')
    return rows,manifest


def status():
    print('NEW ROOT',ROOT)
    for p in (INPUT/'input_manifest.json',REF/'run_config.json',REF/'manifest.json',EVAL/'manifest.json'):
        print(p.relative_to(ROOT),p.exists())
    if (INPUT/'input_manifest.json').is_file():
        check();print('PASS new parent disjointness/source/protocol SHA')
    if REF.is_dir():
        legacy.status()


def build_external_records(tables,targets,rows):
    from pi_run_muon_lut_validation import PRIMARY_FIELDS, SEED_FIELDS
    from pi_run_native_muon_b1 import CROSSING_COLUMNS
    primaries=legacy.read_csv(REF/'muon_reference_heed_primaries.csv',PRIMARY_FIELDS)
    seeds=legacy.read_csv(REF/'muon_reference_seed.csv',SEED_FIELDS)
    cross=legacy.read_csv(REF/'muon_reference_crossing.csv',CROSSING_COLUMNS)
    grouped=defaultdict(list)
    for e in primaries:grouped[e['crossing_id']].append(e)
    seeddict={}
    for z in seeds:
        k=z['crossing_id'],z['configuration'],int(z['primary_index'])
        if k in seeddict:raise ValueError('Duplicate full Garfield primary index')
        seeddict[k]=z
    # First field of CROSSING_COLUMNS is the C++ composite source muon identity.
    truth={(r[CROSSING_COLUMNS[0]],r['configuration']):r for r in cross}
    if len(truth)!=480 or len(grouped)!=120:
        raise ValueError('Expected 120 x 4 independent full Garfield parent truth and shared Heed populations')
    rec=[]
    for mu in rows:
        parent=legacy.id_of(mu)
        e=sorted(grouped[parent],key=lambda q:int(q['primary_index']))
        if not e:raise ValueError('No TrackHeed primaries for new parent')
        if len({int(z['primary_index']) for z in e})!=len(e):raise ValueError('Duplicate primary')
        for g in GEOMS:
            bcnt=Counter()
            pitch=3. if g.startswith('p3_') else 4.
            for z in e:
                bin_=bin_for(g,mu,z)
                bcnt[bin_]+=1
                source=seeddict[parent,g,int(z['primary_index'])]
                true_u=float(source['seed_u_mm'])
                local=true_u-legacy.round_away(true_u/pitch)*pitch
                # A direct guard against a mismatched physical TrackHeed cell.
                # The physical phase must match the frozen LUT bin implementation.
                testbox=legacy.spatial_cell(g,int(mu['station']),float(mu['u_in_mm']),
                                           float(z['u_rel_cm']),float(z['v_abs_cm']))[0]
                if bin_!=testbox:raise ValueError('Disagreement between independently implemented LUT bin assignments')
                _,expected_local_u,_=legacy.spatial_cell(g,int(mu['station']),float(mu['u_in_mm']),
                                                           float(z['u_rel_cm']),float(z['v_abs_cm']))
                if (abs(local-expected_local_u)>1e-3 and
                    abs(abs(local-expected_local_u)-pitch)>1e-3):
                    raise ValueError('New full Garfield primary wire phase disagrees with LUT primary')
                if abs(float(source['seed_v_mm'])-float(z['v_abs_cm'])*10)>1e-4:
                    raise ValueError('Microscopic Garfield / common TrackHeed electron v mismatch')
                if int(source['primary_index'])!=int(z['primary_index']):
                    raise ValueError('Original primary-to-Garfield reference identity mismatch')
            t=truth[parent,g]
            if int(t['heed_primaries'])!=len(e):
                raise ValueError('Full Garfield reference primary count mismatch')
            rec.append({'parent':parent,'geom':g,'counts':tuple(sorted(bcnt.items())),
                        'n_primary':len(e),
                        'n_target':sum(n for b,n in bcnt.items() if (g,b) in targets),
                        'ref_gain':int(t['ion_births']),
                        'ref_up':int(t['ion_up']),'ref_down':int(t['ion_down']),
                        'ref_top':int(t['ion_top']),'ref_bottom':int(t['ion_bottom'])})
    if len(rec)!=480 or sum(r['n_primary'] for r in rec)%4:
        raise ValueError('Wrong external record geometry multiplicity')
    for i in range(120):
        if tuple(rec[4*i+j]['geom'] for j in range(4))!=GEOMS or len({rec[4*i+j]['parent'] for j in range(4)})!=1:
            raise ValueError('New physical parent geometry alignment changed')
    return rec,len(primaries)


def prepare_generic_training(tables,rec,total_primaries):
    """Exact frozen global bootstrap training-pool construction, generalized
    ONLY from 5881 old Heed electrons to the counted new Heed population.
    The physics kernel/weights/settings are unchanged.
    """
    neighbors=compile_neighbors(tables,True)
    result={};parents=None
    for g in GEOMS:
        rr=[r for r in rec if r['geom']==g]
        if len(rr)!=N_PARENT:raise ValueError('Expected new 120 events')
        ids=[];lengths=[]
        for r in rr:
            ids.extend(b for b,n in r['counts'] for _ in range(n))
            lengths.append(r['n_primary'])
        if len(ids)!=sum(lengths) or len(ids)!=total_primaries:
            raise ValueError('Lost/duplicated new TrackHeed primaries')
        ids=np.asarray(ids,np.intp)
        old,extra=tables[g];ncell=len(old)
        target=np.zeros(ncell,bool)
        extraarray=np.zeros((ncell,80,5),np.int64)
        for b,v in extra.items():target[b]=True;extraarray[b]=v
        donor=np.zeros((ncell,24),np.intp)
        cdf=np.ones((ncell,24),float)
        has=np.zeros(ncell,bool)
        for b,(c,w) in enumerate(neighbors[g]):
            if len(c)>24:raise ValueError('Neighbor kernel outside fixed 5x5 grid')
            if len(c):donor[b,:len(c)]=c;cdf[b,:len(c)]=w;has[b]=True
        pp=[r['parent'] for r in rr]
        if parents is None:parents=pp
        elif parents!=pp:raise ValueError('Lost one-to-one parent alignment')
        result[g]={'r':rr,'old':old,'extra':extraarray,'tgt':target,'ids':ids,
                   'offs':np.cumsum([0]+lengths[:-1]),
                   'donor':donor,'cdf':cdf,'has':has}
    return result


def evaluate(reps,seedval):
    if not 500<=reps<=10000:raise ValueError('Use 500..10000 bootstrap repetitions')
    if EVAL.exists():raise FileExistsError('Preserve original blind evaluation directory: '+str(EVAL))
    rows,meta=check()
    _,refstat,refsha=legacy.verify_full()  # true independent Geant4/Heed/Garfield references are mandatory
    tables,targets,*rest=load_lut_native()  # all prior 84k+16k archive source SHA checked
    rec,n_primary=build_external_records(tables,targets,rows)
    if not 0<n_primary<100000:raise ValueError('Unexpected count of new Heed primaries')
    data=prepare_generic_training(tables,rec,n_primary)
    t0=time.monotonic()
    # Same original physical bootstrap implementation and frozen guarded 5x5
    # kernel as old120; no optimization on external truth allowed.
    sampled=simulate(data,reps,seedval)
    per=metrics(data,sampled)
    src={'new_external_parent_input_sha256':meta['sha256_input'],
         'new_full_Garfield_manifest_sha256':refsha,
         'original_unseen_geant4_chunk2_source_sha256':meta['frozen_source_sha256'],
         'old_120_used_parent_input_sha256':meta['original_old120_input_sha256'],
         'analysis_scripts_sha256':fingerprint(),
         'source_conditioning':'New Geant4 parent EVENT identities excluded from previous 120 before microscopic run; original 84k A20+16k A100 training frozen'}
    result=report(per,reps,src)
    result['external_holdout']=True
    result['old_source_population_n_Heed_electrons']=5881
    result['new_holdout_Heed_electrons']=n_primary
    result['runtime_seconds']=round(time.monotonic()-t0,3)
    result['sampling_design']='10 per station/signed-angle stratum, NOT unbiased actual muon flux distribution'
    result['comparison_limits']='Exact ion-birth and 5-sector response only; no cathode induced U/V current or signal, no magnetic field; one true microscopic realization per new parent and configuration.'
    EVAL.mkdir(parents=True)
    perfile=EVAL/'per_external_muon_4_models.csv'
    with perfile.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(per[0]));w.writeheader();w.writerows(per)
    sf=EVAL/'summary.json';sf.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    (EVAL/'manifest.json').write_text(json.dumps({'schema':'NA60_EXTERNAL_NEW120_FIXED_GUARDED_LUT_EVAL_V1',
          'status':'complete','files_sha256':{p.name:sha(p) for p in (perfile,sf)}},indent=2)+'\n')
    print(json.dumps({'external':True,'n_new_Heed_primaries':n_primary,
          'overall':result['overall'],
          'parent_clustered_coverage_deltas':result['parent_clustered_coverage_deltas']},indent=2),flush=True)
    print('PASS NEW external 120 Geant4 parents, cached new TrackHeed and full Garfield vs FROZEN A20/A100/guarded LUT; original data unchanged',flush=True)


def selftest():
    # Synthetic checks only; do NOT load original Geant4 or run physics.
    from pi_spatial_lut_muon120_hier_v1 import selftest as older_test
    older_test()
    old=[{'source':'Jpsi','chunk':'2','event':str(i),'track_id':str(i)} for i in range(120)]
    new=[{'source':'Jpsi','chunk':'2','event':str(i+120),'track_id':str(i+120)} for i in range(120)]
    assert not {(x['source'],x['chunk'],x['event']) for x in old}&{(x['source'],x['chunk'],x['event']) for x in new}
    sample=np.zeros((120,4,4,5),np.int64)
    assert sample.shape==(120,4,4,5)
    # Regression of new-source variable Heed cardinality (not old fixed 5881).
    fake={}; rec=[]
    for g in GEOMS:
        cells=900 if g.startswith('p3_') else 1200
        fake[g]=(np.zeros((cells,20,5),dtype=np.int64),{})
    for i in range(120):
        for g in GEOMS:
            rec.append({'parent':str(i),'geom':g,'counts':((47,1),),
                        'n_primary':1,'n_target':0})
    prepared=prepare_generic_training(fake,rec,120)
    assert all(prepared[g]['ids'].size==120 and
               prepared[g]['offs'].shape==(120,) for g in GEOMS)
    print('PASS synthetic NEW Heed cardinality and fixed guarded donors in separate external holdout',flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('mode',choices=('selftest','prepare','check','preflight','run','finish','evaluate','status'))
    p.add_argument('--hours',type=float,default=7.75)
    p.add_argument('--max-chunks',type=int,default=12)
    p.add_argument('--reps',type=int,default=3000)
    p.add_argument('--seed',type=int,default=20261101)
    args=p.parse_args()
    if args.mode=='selftest':selftest()
    elif args.mode=='prepare':prepare()
    elif args.mode=='check':check();print('PASS new external parent IDs and model pre-freeze',flush=True)
    elif args.mode=='preflight':check();legacy.preflight()
    elif args.mode=='run':check();legacy.run(args.hours,args.max_chunks)
    elif args.mode=='finish':check();legacy.finish()
    elif args.mode=='evaluate':evaluate(args.reps,args.seed)
    else:status()


if __name__=='__main__':main()