#!/usr/bin/env python3
"""Strict read-only audit of SHA-frozen external new120 Geant4/Heed/Garfield.

Fixes ONLY a new analysis path: microscopic Garfield unsuccessful transport
from a primary electron physically INSIDE 30-um anode metal has zero ion
births, though an old endpoint-radius heuristic counts the electron as
collected. No old physics, original validator, SHA files, or source modified.

Unknown unsuccessful starts, cap, missing primaries, altered CSV/checkpoints,
and incongruent births fail CLOSED.
"""
from __future__ import annotations
import argparse
import json
import math
from collections import Counter, defaultdict

GEOMS=('p3_noFW','p3_FW75','p4_noFW','p4_FW75')
ANODE_RADIUS_MM=.015
FIELD_RADIUS_MM=.0375
METAL_TOL_MM=.000001


def must(ok,msg):
    if not ok:raise ValueError('EXTERNAL METAL PROVENANCE: '+msg)


def nearest_anode(row):
    u=float(row['seed_u_mm']);v=float(row['seed_v_mm'])
    must(math.isfinite(u) and math.isfinite(v),'non-finite electron seed')
    geom=row['configuration']
    must(geom in GEOMS,'unknown geometry '+geom)
    pitch=3. if geom.startswith('p3_') else 4.
    # Match the original C++ nearest-wire phase; 9 anodes in patch [-4,4].
    n=int(math.floor(abs(u/pitch)+.5))*(1 if u>=0 else -1)
    return float(n*pitch),math.hypot(u-n*pitch,v),n


def valid_anode_metal_zero(row):
    u0,d,n=nearest_anode(row)
    if not (-4<=n<=4 and d < ANODE_RADIUS_MM-METAL_TOL_MM):
        return False
    cols=('garfield_ne','garfield_ni','ion_births','collected_e','attached',
          'other','ion_up','ion_down','ion_top','ion_bottom')
    if any(int(row[k])<0 for k in cols):return False
    return (int(row['successful'])==0 and int(row['capped'])==0 and
            int(row['garfield_ne'])==1 and int(row['garfield_ni'])==0 and
            int(row['ion_births'])==0 and
            int(row['attached'])==0 and
            int(row['collected_e'])<=1 and int(row['other'])<=1 and
            int(row['collected_e'])+int(row['other'])<=1 and
            all(int(row[k])==0 for k in
                ('ion_up','ion_down','ion_top','ion_bottom')))


def audit():
    import pi_external_muon120_guarded_lut_v1 as ext
    from pi_run_muon_lut_validation import (SEED_FIELDS,PRIMARY_FIELDS,
          CLUSTER_FIELDS,field_wire_blocked_zero_start)
    from pi_run_native_muon_b1 import CROSSING_COLUMNS

    parents,selected_manifest=ext.check()
    legacy=ext.legacy
    root=ext.REF
    manifestpath=root/'manifest.json'
    must(manifestpath.is_file(),'missing finished full physics reference')
    meta=json.loads(manifestpath.read_text())
    must(meta.get('status')=='complete' and
         meta.get('n_physically_distinct_muons')==120 and
         meta.get('n_completed_chunks')==12,
         'reference must be complete 120 real parents x four geometries')
    must(meta.get('original_geant_input_sha256')==selected_manifest['sha256_input'],
         'Geant4 source manifest differs')
    # Original binary provenance and all 12 original microscopic chunk SHA seals.
    binary=legacy.check_executable()
    must(meta.get('cpp_executable_sha256')==binary,'old native compiled binary differs')
    for i,subset in enumerate(legacy.groups(parents)):
        rec=legacy.read_seal(legacy.chunk_dir(i),i,subset,binary)
        must(legacy.sha(legacy.chunk_dir(i)/'complete.json')==
             meta['original_chunks_sha256'][f'chunk_{i:02d}'],
             f'completed original chunk sha changed: {i}')
    must(len(meta.get('outputs_sha256',{}))>=5,'incomplete aggregated physics SHA list')
    for name,dig in meta['outputs_sha256'].items():
        must(legacy.sha(root/name)==dig,'altered physical archive file '+name)

    from pi_prepare_muon_avalanche_b1 import FIELDS
    from pi_run_muon_lut_validation import id_of
    expected_parents={id_of(x) for x in parents}
    must(len(expected_parents)==120,'not 120 independent parent IDs')
    primary=legacy.read_csv(root/'muon_reference_heed_primaries.csv',PRIMARY_FIELDS)
    clusters=legacy.read_csv(root/'muon_reference_heed_clusters.csv',CLUSTER_FIELDS)
    seeds=legacy.read_csv(root/'muon_reference_seed.csv',SEED_FIELDS)
    cross=legacy.read_csv(root/'muon_reference_crossing.csv',CROSSING_COLUMNS)
    must(len(primary)>0 and len(seeds)==4*len(primary),
         'lost/duplicated Garfield primaries')
    pidx=defaultdict(set)
    pcluster=Counter();cidx=defaultdict(set)
    for row in primary:
        x=row['crossing_id'];i=int(row['primary_index']);ci=int(row['cluster_index'])
        must(x in expected_parents and i>=0 and ci>=0 and i not in pidx[x],
             'invalid primary index or parent')
        pidx[x].add(i);pcluster[x,ci]+=1
        must(all(math.isfinite(float(row[k])) for k in
                 ('u_rel_cm','v_abs_cm','w_rel_cm','t_ns')),
             'non-finite TrackHeed electron')
    must(set(pidx)==expected_parents and
         all(ids==set(range(len(ids))) for ids in pidx.values()),
         'missing TrackHeed parent or non-contiguous electron indices')
    for row in clusters:
        x=row['crossing_id'];i=int(row['cluster_index'])
        must(x in expected_parents and i>=0 and i not in cidx[x] and
             int(row['n_accepted'])==pcluster[x,i] and
             int(row['n_original'])>=int(row['n_accepted']),
             'TrackHeed cluster membership invalid')
        cidx[x].add(i)
    must(set(cidx)==expected_parents and
         all(ids==set(range(len(ids))) for ids in cidx.values()),
         'missing TrackHeed clusters')

    crossing={(r[CROSSING_COLUMNS[0]],r['configuration']):r for r in cross}
    must(len(crossing)==480 and all((x,g) in crossing
          for x in expected_parents for g in GEOMS),
         'missing or duplicate four-geometry full reference crossing')
    for (x,g),r in crossing.items():
        must(int(r['heed_primaries'])==len(pidx[x]) and
             int(r['heed_clusters'])==len(cidx[x]),
             'crossing Heed primary/cluster cardinality disagreement')
    seen=set();ions=Counter();sector=defaultdict(Counter)
    failures=Counter();affected=defaultdict(set);details=[]
    for row in seeds:
        x=row['crossing_id'];g=row['configuration'];idx=int(row['primary_index'])
        k=(x,g,idx)
        must(x in expected_parents and g in GEOMS and
             idx in pidx[x] and k not in seen,'duplicate/orphan Garfield seed')
        seen.add(k)
        must(row['capped']=='0','full microscopic avalanche cap encountered '+str(k))
        birth=int(row['ion_births']);ne=int(row['garfield_ne']);ni=int(row['garfield_ni'])
        must(birth>=0 and ne>=0 and ni>=0,'negative physical microscopic counts '+str(k))
        for s in ('ion_up','ion_down','ion_top','ion_bottom'):
            must(int(row[s])>=0,'negative ion sector '+str(k))
        must(int(row['ion_up'])+int(row['ion_down'])<=birth and
             int(row['ion_top'])+int(row['ion_bottom'])<=birth,
             'ion sector larger than all births '+str(k))
        ions[x,g]+=birth
        sector[x,g]['up']+=int(row['ion_up'])
        sector[x,g]['down']+=int(row['ion_down'])
        sector[x,g]['top']+=int(row['ion_top'])
        sector[x,g]['bottom']+=int(row['ion_bottom'])
        if row['successful']=='1':continue
        # A physically impossible gas-origin electron inside 30um anode metal
        # does not initiate drift. It has zero births. The old approximate
        # collection-radius classifier MAY nevertheless set collected_e=1.
        if valid_anode_metal_zero(row):tag='anode_metal_zero'
        elif field_wire_blocked_zero_start(row):tag='field_wire_metal_zero'
        else:raise RuntimeError('UNEXPLAINED unsuccessful Garfield seed; '
             'original SHA archive preserved: '+str(k)+' '+str(row))
        failures[g,tag]+=1;affected[g].add(x)
        center,d,n=nearest_anode(row)
        details.append({'parent':x,'geometry':g,'primary_index':idx,
                        'classification':tag,'seed_u_mm':float(row['seed_u_mm']),
                        'seed_v_mm':float(row['seed_v_mm']),
                        'nearest_anode_axis_mm':center,
                        'radial_distance_from_anode_axis_mm':d,
                        'raw_collected_e':int(row['collected_e']),
                        'raw_successful':0,'ion_births':0})
    must(len(seen)==4*len(primary),'missing Garfield avalanche seed')
    for (x,g),cr in crossing.items():
        must(ions[x,g]==int(cr['ion_births']),
             'sum of microscopic ion births differs from muon truth '+x+' '+g)
        for s,k in (('up','ion_up'),('down','ion_down'),
                    ('top','ion_top'),('bottom','ion_bottom')):
            must(sector[x,g][s]==int(cr[k]),
                 'muon sector does not agree with primary ion census '+x+' '+g)
    info={
      'schema':'NA60_EXTERNAL_NEW120_METAL_START_READONLY_V1',
      'status':'PASS all original SHA seals and physical parent/ion-sector census',
      'reference_manifest_sha256':legacy.sha(manifestpath),
      'original_source_sha256':selected_manifest['sha256_input'],
      'n_parents':120,'n_original_trackheed_electrons':len(primary),
      'n_garfield_electron_avalanches':len(seeds),
      'by_geometry':{g:{'n_failed_in_anode_metal':failures[g,'anode_metal_zero'],
                        'n_failed_in_field_wire_metal':failures[g,'field_wire_metal_zero'],
                        'n_affected_parents':len(affected[g]),
                        'affected_parent_ids':sorted(affected[g])}
                     for g in GEOMS},
      'exceptions':details,
      'warning':('Zero ion births for in-metal starts are retained in the count '
                  'truth. Any raw collected_e=1 for these starts is a loose '
                  'endpoint-radius classification, NOT a valid electron '
                  'transport or physical collected charge. Never use it as '
                  'an input to future electron induced-signal physics.')}
    return info


def selftest():
    base={'configuration':'p3_noFW','seed_u_mm':'0.00237593484574722',
          'seed_v_mm':'0.0107657602328344','successful':'0','capped':'0',
          'garfield_ne':'1','garfield_ni':'0','ion_births':'0',
          'collected_e':'1','attached':'0','other':'0','ion_up':'0',
          'ion_down':'0','ion_top':'0','ion_bottom':'0'}
    must(valid_anode_metal_zero(base),'known exact archived 30um anode start rejected')
    _,distance,_=nearest_anode(base)
    must(abs(distance-.011025)<.0001,'wrong radial distance')
    for name,value in (('seed_v_mm','.020'),('ion_births','1'),
                       ('garfield_ne','2'),('capped','1'),('successful','1'),
                       ('other','1')):
        q=dict(base);q[name]=value
        must(not valid_anode_metal_zero(q),'unsound acceptance '+name)
    q=dict(base);q['configuration']='p3_FW75'
    must(valid_anode_metal_zero(q),'anode metal classification varies by FW incorrectly')
    print('PASS exact observed anode-metal geometry (r=11um < 15um), strict zeros and fail-closed checks; NO Garfield')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('mode',choices=('selftest','audit'))
    args=p.parse_args()
    if args.mode=='selftest':selftest()
    else:
        result=audit()
        print(json.dumps(result,indent=2,sort_keys=True))
        print('PASS external independent new120 SHA-sealed full Garfield: '
              'only explicitly proven in-metal zero-start exceptions')
if __name__=='__main__':main()
