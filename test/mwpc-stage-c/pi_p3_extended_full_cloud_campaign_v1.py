#!/usr/bin/env python3
"""NA60+ 3mm-only independent full Garfield ion-birth cloud EXTENDED pilot.

Separately assess no-field-wire and FW75 geometries. Source coordinates are
unaltered SHA-frozen original D500 *initial electrons*, not drift endpoints.
Uses 40 independent D500 passes 100..139 PER SOURCE CELL; original A20/A100
use passes 0..99. Also explicitly excludes every previous 48-cloud pilot start.

This module reuses the already-tested physical full-birth native C++ and
strict per-chunk SHA validation from pi_p3_full_cloud_convergence_v1, but
redirects its Python campaign paths to a NEW isolated directory. Nothing in
old files is changed, and no induced cathode signal is computed here.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

import pi_p3_full_cloud_convergence_v1 as core
from pi_audit_avalanche_lut_reuse import BASE, CONFIGS, EXPORT_FIELDS, locations, sha
from pi_run_lut500_overnight import offset_for, STRIDE, source_maps

REPO = Path(__file__).resolve().parents[2]
ROOT = BASE/'pi_p3_extended_full_cloud_v1'
INPUT = ROOT/'input'
CHUNKS = ROOT/'full_cloud_chunks'
OLD_PILOT = BASE/'pi_p3_full_cloud_convergence_v1'
GEOMS = ('p3_noFW','p3_FW75')
GROUP = 8
PASS_START, PASS_END = 100, 140  # no A20/A100 and enough near-wire statistics
CAPS = (32,128,256,512,1024)
SCHEMA = 'NA60_P3_EXTENDED_FULL_CLOUD_INPUT_V1'
SELECTION_TAG = 'NA60-P3-EXTENDED-FULL-CLOUD-20261104-V1'
FIELDS = tuple(EXPORT_FIELDS)


def must(x,why):
    if not x: raise ValueError(why)


def load(p):return json.loads(p.read_text())
def save(p,obj):core.write_json(p,obj)
def score(*items):
    b='|'.join(map(str,items)).encode('utf-8')
    return int.from_bytes(hashlib.sha256(b).digest()[:8],'big')


def sample_for(geom,box,pass_index):
    """Exact frozen 500-sample D500 permutation, no new simulated inputs."""
    must(geom in GEOMS and 0<=box<900 and
         PASS_START<=pass_index<PASS_END, 'Invalid reserved D500 trial pass')
    return (offset_for(geom,box)+STRIDE*pass_index)%500


def old_pilot_keys():
    """Check original 48-site pilot provenance, without changing its data."""
    man=load(OLD_PILOT/'input/selection_manifest.json')
    path=OLD_PILOT/'input/exact_D500_p3_unseen.csv'
    must(man.get('schema')=='NA60_P3_CLOUD_INPUT_V1' and
         man.get('n_selected')==48 and man.get('per_geometry')==24 and
         sha(path)==man.get('input_sha256'),
         'Previous 48-site full-cloud source not SHA-valid')
    rows=list(core.csv_records(path,FIELDS))
    keys={(r['geometry'],int(r['box_id']),int(r['sample'])) for r in rows}
    must(len(rows)==len(keys)==48 and
         {r['geometry'] for r in rows}==set(GEOMS),
         'Previous 48-site exclusion set changed')
    return keys,man['input_sha256']


def choose_balanced(bins,geom,label,count):
    """Balance over DISTINCT D500 spatial cells, then extra trials/cell.

    Important around metal: only a few D500 cells exist near the anode or
    FW, so reusing independent starts in a cell is not pseudo-replication of
    the same avalanche. Each D500 (cell,sample) has its own collision seed.
    Inference/bootstrap must, however, be clustered by spatial cell when
    testing the spatially conditional model, not naive independent-ion CI.
    """
    groups=bins.get(label,{})
    cells=sorted(groups,key=lambda b:(score(SELECTION_TAG,'box',geom,label,b),b))
    ordered={b:sorted(groups[b],key=lambda rr:(
        score(SELECTION_TAG,'start',geom,label,b,rr[0]['sample']),rr[0]['sample']))
             for b in cells}
    chosen=[];seen=set()
    for j in range(PASS_END-PASS_START):
        for b in cells:
            if j>=len(ordered[b]):continue
            row,pass_index=ordered[b][j]
            key=(geom,int(row['box_id']),int(row['sample']))
            must(key not in seen,'Duplicate candidate after round-robin')
            seen.add(key)
            chosen.append((row,pass_index))
            if len(chosen)==count:return chosen,len(cells),sum(map(len,groups.values()))
    raise ValueError(f'Not enough unique gas-phase D500 starts for {geom}/{label}: '
                     f'need {count}, found {len(chosen)}; never invent more starts')


def selection_from_sources(per_geometry):
    must(per_geometry in (480,960), 'Choose exactly 480 or 960 starts PER geometry')
    oldkeys,oldsha=old_pilot_keys()
    source,indexsha=source_maps()  # original D500 SHA manifest
    rows_by_geom={};strata_by_geom={};count_by_geom={};source_sha={}
    for geom,pitch,field_diameter in CONFIGS:
        if geom not in GEOMS:continue
        _,raw=locations(500,pitch,field_diameter)
        must(sha(raw)==source[geom]['lut_csv_sha256'],
             'Original immutable D500 source checksum changed '+geom)
        source_sha[geom]=sha(raw)
        candidates=defaultdict(lambda:defaultdict(list))
        passage_for={b:{sample_for(geom,b,k):k
                        for k in range(PASS_START,PASS_END)}
                     for b in range(900)}
        must(all(len(t)==PASS_END-PASS_START for t in passage_for.values()),
             'D500 pass mapping is not one-to-one')
        hit=0;exclude=0
        with raw.open(newline='') as fh:
            reader=csv.DictReader(fh)
            columns=set(reader.fieldnames or [])
            need={'box_id','sample','iu','iv','u_start_mm','v_start_mm',
                  'w_start_mm','random_seed','drift_ok','collected',
                  'field_wire_collected','wire_index','phi_hit_deg',
                  'drift_time_ns','last_step_energy_eV'}
            must(need.issubset(columns),'Frozen D500 missing column: '+geom)
            for source_row in reader:
                box=int(source_row['box_id']);sample=int(source_row['sample'])
                must(0<=box<900,'Bad D500 spatial cell')
                pass_index=passage_for[box].get(sample)
                if pass_index is None:continue
                hit+=1
                ident=(geom,box,sample)
                if ident in oldkeys:
                    exclude+=1;continue
                iu,iv=int(source_row['iu']),int(source_row['iv'])
                u,v,w=[float(source_row[k]) for k in
                       ('u_start_mm','v_start_mm','w_start_mm')]
                must(box==iu*30+iv and math.isfinite(u) and
                     math.isfinite(v) and math.isfinite(w) and abs(w)<1.e-8 and
                     -1.5+iu*.1-1.e-6<=u<=-1.5+(iu+1)*.1+1.e-6 and
                     -3+iv*.2-1.e-6<=v<=-3+(iv+1)*.2+1.e-6,
                     'Original D500 spatial site invalid '+str(ident))
                label=core.stratum(u,v,geom)
                if label is None:continue # physically inside anode/FW metal
                record={'lut_trials':500,'geometry':geom,
                    'original_lut_sha256':source[geom]['lut_csv_sha256'],
                    'box_id':box,'iu':iu,'iv':iv,'sample':sample,
                    'random_seed':int(source_row['random_seed']),
                    'u_start_mm':source_row['u_start_mm'],
                    'v_start_mm':source_row['v_start_mm'],
                    'w_start_mm':source_row['w_start_mm'],
                    'old_drift_ok':source_row['drift_ok'],
                    'old_collected':source_row['collected'],
                    'old_field_wire_collected':source_row['field_wire_collected'],
                    'old_wire_index':source_row['wire_index'],
                    'old_phi_hit_deg':source_row['phi_hit_deg'],
                    'old_drift_time_ns':source_row['drift_time_ns'],
                    'old_last_step_energy_eV':source_row['last_step_energy_eV']}
                must(set(record)==set(FIELDS),'Original source record changed')
                candidates[label][box].append((record,pass_index))
        must(hit==900*(PASS_END-PASS_START),
             'Missing/duplicated D500 reserved trial positions '+geom)
        regions=('near_anode','transition','bulk') if geom=='p3_noFW' else (
                 'near_anode','near_field','transition','bulk')
        labels=[(regime,side) for regime in regions
                for side in ('positive','negative')]
        per_group=per_geometry//len(labels)
        must(per_group*len(labels)==per_geometry,
             'Each geometry must be evenly balanced over its OWN strata')
        counts={};chosen={}
        for label in labels:
            selected,nbox,available=choose_balanced(candidates,geom,label,per_group)
            chosen[label]=selected
            counts[str(label)]={'eligible_sites':available,'eligible_boxes':nbox,
                                'chosen':per_group,
                                'chosen_distinct_boxes':len({int(row['box_id']) for row,k in selected})}
        # Balanced round robin within geometry; early SHA checkpoint spans regimes.
        ordered=[]
        for i in range(per_group):
            for label in labels:ordered.append((chosen[label][i],label))
        must(len(ordered)==per_geometry,'Nonuniform strata allocation '+geom)
        rows_by_geom[geom]=[pair[0][0] for pair in ordered]
        strata_by_geom[geom]=[{'geometry':geom,'box_id':int(pair['box_id']),
            'sample':int(pair['sample']),'pass_index':pass_index,
            'regime':label[0],'v_sign':label[1]}
            for (pair,pass_index),label in ordered]
        count_by_geom[geom]={'strata':counts,'excluded_pilot_sites':exclude,
                             'native_candidates_reserved':hit}
    allrows=[];tags=[]
    for i in range(per_geometry):
        for geom in GEOMS:
            allrows.append(rows_by_geom[geom][i])
            tags.append(strata_by_geom[geom][i])
    keys={(r['geometry'],int(r['box_id']),int(r['sample'])) for r in allrows}
    must(len(keys)==len(allrows)==2*per_geometry and not keys.intersection(oldkeys),
         'Duplicate or old-pilot D500 initial electron starts')
    return allrows,tags,source_sha,indexsha,count_by_geom,oldsha


def check_plan():
    source_path=INPUT/'exact_D500_p3_unseen.csv'
    manifest=load(INPUT/'selection_manifest.json')
    rows=list(core.csv_records(source_path,FIELDS))
    tags=load(INPUT/'strata.json')
    must(manifest.get('schema')==SCHEMA and
         manifest.get('per_geometry') in (480,960) and
         manifest.get('n_selected')==2*manifest['per_geometry']==len(rows)==len(tags) and
         manifest.get('group_size')==GROUP and
         tuple(manifest.get('geometries',[]))==GEOMS and
         sha(source_path)==manifest.get('input_sha256') and
         sha(INPUT/'strata.json')==manifest.get('strata_sha256'),
         'Extended p3 source/manifests changed')
    oldkeys,oldsha=old_pilot_keys()
    must(oldsha==manifest['original_48pilot_input_sha256'],
         'Immutable old48 selection changed')
    observed=set();strata=Counter();last_geom=None
    for i,(r,t) in enumerate(zip(rows,tags)):
        geom=r['geometry'];box=int(r['box_id']);s=int(r['sample'])
        k=(geom,box,s)
        must(k not in observed and k not in oldkeys and
             geom==GEOMS[i%2] and geom==t['geometry'] and
             box==int(t['box_id']) and s==int(t['sample']) and
             int(t['pass_index']) in range(PASS_START,PASS_END) and
             sample_for(geom,box,int(t['pass_index']))==s and
             core.stratum(float(r['u_start_mm']),float(r['v_start_mm']),geom)==
             (t['regime'],t['v_sign']),
             'Reused, metal-embedded or shifted D500 selection '+str(k))
        observed.add(k);strata[geom,t['regime'],t['v_sign']]+=1
    for geom in GEOMS:
        labels=[('near_anode','positive'),('near_anode','negative'),
                ('transition','positive'),('transition','negative'),
                ('bulk','positive'),('bulk','negative')]
        if geom=='p3_FW75':
            labels[2:2]=[('near_field','positive'),('near_field','negative')]
        each=manifest['per_geometry']//len(labels)
        for reg,side in labels:
            must(strata[geom,reg,side]==each,
                 'Spatial near-wire/bulk strata no longer balanced '+str((geom,reg,side)))
    return rows,tags,manifest


def prepare(n_per_geom):
    must(not ROOT.exists(),'Refuse to overwrite any previous extended pilot: '+str(ROOT))
    rows,tags,source_shas,indexsha,counts,oldsha=selection_from_sources(n_per_geom)
    INPUT.mkdir(parents=True)
    source_path=INPUT/'exact_D500_p3_unseen.csv'
    core.write_csv(source_path,rows)
    save(INPUT/'strata.json',tags)
    save(INPUT/'selection_manifest.json',{
        'schema':SCHEMA,'n_selected':len(rows),'per_geometry':n_per_geom,
        'group_size':GROUP,'geometries':GEOMS,'caps':CAPS,
        'pass_range_reserved':[PASS_START,PASS_END],
        'original_48pilot_input_sha256':oldsha,
        'original_D500_source_SHA256':source_shas,
        'original_D500_index_manifest_SHA256':indexsha,
        'selection':'Balanced stratum+hemisphere original D500 initial gas electrons; '
            '40 independent random passes 100..139 per source cell; '
            'exclude all A20/A100 and old48 original starts; no new TrackHeed/Geant4',
        'counts_by_geometry':counts,
        'anode_voltage_V':1800,'field_wire_diameter_um':75,
        'gas':'Ar/CO2 70:30, 293.15K, 760 Torr',
        'input_sha256':sha(source_path),'strata_sha256':sha(INPUT/'strata.json')})
    check_plan()
    print('PASS PREPARED extended full-ion cloud starts',len(rows),
          'exact D500 distinct (geometry,box,sample),',n_per_geom,'PER p3 geometry',flush=True)
    for g in GEOMS:
        print('STRATA',g,json.dumps(counts[g],sort_keys=True),flush=True)


def configure_core():
    # Reuse strict production-tested v1 physics reader/writer/checkpoint code,
    # but redirect ALL the module-level paths in this process to an isolated
    # new directory. Never patch its on-disk source or its C++ executable.
    core.ROOT=ROOT
    core.INPUT=INPUT
    core.CHUNKS=CHUNKS
    core.CAPS=CAPS
    core.RNG_TAG=SELECTION_TAG
    core.check_plan=check_plan


def selftest():
    must(GEOMS==core.GEOMS and GROUP==core.GROUP,
         'Unexpected prior validated physical runner interface')
    must(CAPS==(32,128,256,512,1024), 'Expected cloud cap benchmark')
    must(core.classify(0.,0.,'p3_noFW') is None,
         'Anode metal must not be a valid D500 start')
    must(core.classify(1.5,0.,'p3_FW75') is None,
         'Field wire metal must not be a valid D500 start')
    for g in GEOMS:
        for b in (0,1,14,450,899):
            prev={(offset_for(g,b)+STRIDE*k)%500 for k in range(100)}
            new=[sample_for(g,b,k) for k in range(PASS_START,PASS_END)]
            must(len(set(new))==len(new)==40 and not prev.intersection(new),
                 'New exact-start source reuses original 100 random passes')
    # A single spatial cell can contribute more than one *independent* trial,
    # but no exact (geometry,box,sample) identity can appear more than once.
    fake={0:{i:[({'box_id':i,'sample':100+j},100+j)
                for j in range(4)] for i in range(3)}}
    x,nb,total=choose_balanced(fake,'p3_noFW',0,8)
    must(len(x)==8 and nb==3 and total==12 and
         len({(r['box_id'],r['sample']) for r,k in x})==8 and
         len({r['box_id'] for r,k in x})==3,
         'Synthetic balanced 3-cell independent-source selection incorrect')
    toy=__import__('numpy').array([[0.,.1,0.,0.],[.01,.02,0.,.2],
                                   [-.03,-.1,0.,.8],[.04,.05,0.,1.]],float)
    for cap in CAPS:
        v=core.compression_errors(toy,toy)
        must(all(x is None or abs(x)<1e-10 for x in v.values()),
             'Full ion cloud should self-match at all reservoir caps')
    # A live status reader MUST tolerate exactly one unsealed (in-flight)
    # physical chunk, while rejecting gaps and inconsistent future seals.
    import tempfile
    with tempfile.TemporaryDirectory() as work:
        paths=[Path(work)/f'chunk_{i:04d}' for i in range(4)]
        must(scan_chunk_state(paths)==([],None),
             'Empty campaign must be reportable')
        paths[0].mkdir()
        must(scan_chunk_state(paths)==([],0),
             'Active first chunk incorrectly treated as corrupt')
        (paths[0]/'complete.json').write_text('{}')
        paths[1].mkdir()
        must(scan_chunk_state(paths)==([0],1),
             'Sealed first + active second chunk must be reportable')
        (paths[1]/'complete.json').write_text('{}')
        must(scan_chunk_state(paths)==([0,1],None),
             'Two sealed chunks incorrectly reported')
        paths[3].mkdir()
        try:
            scan_chunk_state(paths)
            raise AssertionError('A gap between physical chunks was accepted')
        except ValueError as exc:
            must('Gap in physical checkpoint' in str(exc),
                 'Wrong error for a skipped physical checkpoint')
        paths[2].mkdir()
        must(scan_chunk_state(paths)==([0,1],2),
             'Two consecutive pending physical chunks accepted')
        (paths[2]/'complete.json').write_text('{}')
        must(scan_chunk_state(paths)==([0,1,2],3),
             'Final unfinished chunk not reported as pending')
    print('PASS extended p3 SHA source selection, disjoint 40-pass D500 starts, '
          'stratum/cell balancing, 5-cap cloud metrics and live-status checkpoint '
          'safety; NO new Garfield',flush=True)


def scan_chunk_state(paths):
    """Read-only snapshot of checkpoint directories.

    A writer creates a chunk directory and physics.log BEFORE it writes the
    final SHA seal. Therefore the newest *unsealed* chunk is expected while
    the background Garfield++ subprocess is alive. Never count it as done.
    Reject gaps, multiple unfinished chunks, or a sealed chunk after a
    pending chunk; those indicate an inconsistent physical campaign.
    """
    verified_ids=[]
    pending=None
    missing_seen=False
    for i,d in enumerate(paths):
        if not d.exists():
            missing_seen=True
            continue
        must(d.is_dir(), 'Checkpoint path is not a directory: '+str(d))
        must(not missing_seen,
             'Gap in physical checkpoint directories before '+str(d))
        seal=d/'complete.json'
        if seal.is_file():
            must(pending is None,
                 'A later sealed chunk exists beyond an unsealed chunk: '+str(d))
            verified_ids.append(i)
        else:
            must(not seal.exists(),
                 'Malformed physical checkpoint seal path: '+str(seal))
            must(pending is None,
                 'Multiple unsealed physical chunks found: '+str(d))
            pending=i
    return verified_ids,pending


def status():
    # Deliberately NOT core.completed(): the strict completed() API is used
    # by final analysis and must fail on unsealed chunks. Status, however,
    # must be safe to query WHILE the one permitted background writer runs.
    rows,_,manifest=check_plan()
    binary=core.verified_binary()
    blocks=core.batch(rows)
    ids,pending=scan_chunk_state([core.chunk_path(i)
                                  for i in range(len(blocks))])
    done=[core.audit_chunk(i,blocks[i],binary) for i in ids]
    n=sum(x['n_avalanches'] for x in done)
    births=sum(x['n_ion_births'] for x in done)
    must(n==len(ids)*GROUP and len(ids)==(pending if pending is not None
                                         else len(ids)),
         'Non-contiguous/invalid SHA-verified completed chunk census')
    print('PASS SHA extended cloud COMPLETED',n,'of',len(rows),
          'exact D500 electrons (sealed only)')
    print('VERIFIED CHUNKS',len(done),'/',len(blocks))
    print('FULL MICROSCOPIC ION BIRTHS (verified)',births)
    print('p3_noFW',n//2,'p3_FW75',n//2,
          '; no ion drift or cathode signal')
    if pending is not None:
        d=core.chunk_path(pending)
        print('UNSEALED CHUNK',pending,'NOT COUNTED; may be ACTIVE or INTERRUPTED')
        print('Check the writer PID and latest log. Do NOT launch a second run')
        print('or remove',str(d),'while the original writer may be running.')
    elif len(done)<len(blocks):
        print('NEXT NOT YET STARTED CHUNK',len(done))
    else:
        print('PASS ALL SHA-SEALED 120 / 120 FULL CLOUD CHUNKS')


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('mode',choices=('selftest','prepare','preflight','run','status','analyze'))
    ap.add_argument('--per-geometry',type=int,default=480)
    ap.add_argument('--max-chunks',type=int,default=1)
    ap.add_argument('--hours',type=float,default=6.)
    ap.add_argument('--reps',type=int,default=32)
    ap.add_argument('--partial',action='store_true')
    args=ap.parse_args()
    configure_core()
    if args.mode=='selftest':selftest()
    elif args.mode=='prepare':prepare(args.per_geometry)
    elif args.mode=='preflight':core.preflight()
    elif args.mode=='run':core.run(args.max_chunks,args.hours)
    elif args.mode=='status':status()
    else:core.analyze(args.reps,args.partial)

if __name__=='__main__':main()
