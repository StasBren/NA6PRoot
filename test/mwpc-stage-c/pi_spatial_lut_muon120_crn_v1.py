#!/usr/bin/env python3
"""Read-only source-conditioned evaluation of spatially pooled Garfield A20/A100
avalanche libraries on 120 actual Geant4+TrackHeed parent muons (480 correlated
muon-by-geometry comparisons). NO new detector or Garfield simulations.

Archived mode: --published-dir with all four trials_*, trials_extra_*,
native_muon_crossings_input.csv, muon_reference_heed_primaries.csv,
muon_reference_crossing.csv, muon_A20_vs_targeted_A100.csv,
target_cells.csv and old manifest.json. Each source SHA-verified.
Native mode: --native uses existing project's SHA-verified native input loader.

The 120 Heed source positions were used previously to select A100 cells and
spatial pooling hyperparameters were tuned on high-variance A100-selected
cells; this is therefore NOT blind external model validation.
"""
from __future__ import annotations
import argparse, collections, csv, hashlib, json, math, statistics
from pathlib import Path
import numpy as np

GEOMS = ('p3_noFW','p3_FW75','p4_noFW','p4_FW75')
COLUMNS = ('ion_births','ion_up','ion_down','ion_top','ion_bottom')
MODES = ('A20','A100','guarded_A20','guarded_A100','unguarded_A20','unguarded_A100')
EXPECTED = {
 'native_muon_crossings_input.csv':'64338d9af0cff9b96e50b6d44b8f8b71fc5f29929aa89d71568c4d37384be9ba',
 'muon_reference_heed_primaries.csv':'66e0db1366f4995034fe8f32996d0afc0a9ea44ad8b0af1287bdc8280b9130a9',
 'muon_reference_crossing.csv':'45dc019079830dff083c8c148e273c36153fd1e4b1291ff20e287bcf8d5adc35',
 'trials_extra_p3_noFW.csv':'4922e5930dfa3a1ffa7f6dfe76f1408e61e76c008fa6c040fe14d64ec0e74f67',
 'trials_extra_p3_FW75.csv':'a9b7c2cfafdc5d910b0b5b779468ecfdd166ec4005331855a8c7a793a02c4ab4',
 'trials_extra_p4_noFW.csv':'c24f411f28e96197924d51bbb6906de6723a9eb6a93372648e24e4b3fd1576ac',
 'trials_extra_p4_FW75.csv':'b7e2a8d7d184509f0e7b98739067c81a20b64f74549fc2a6a68bdaa1912a4e22',
 'target_cells.csv':'1f3320e0b336266eec5aa1fe41acacaf25baf3d421b59dbdde5788f7df31b8b6',
}
A20_SHA='8ba639e3f082192d716d328775e13743ca72f9035ebfe703b1f722a202fb006c'
A100_SHA='2a8520c61dc38dc2b32aa31319c3665e36d9ed86f9119e8f3737f32849d7da82'

def sha(path):
 h=hashlib.sha256()
 with open(path,'rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()

def rows(path):
 with open(path,newline='') as f:return list(csv.DictReader(f))

def load_archive(root):
 root=Path(root);m=json.loads((root/'manifest.json').read_text())
 if sha(root/'manifest.json')!=A20_SHA or m.get('n_completed_avalanches')!=84000:raise ValueError('Wrong A20 manifest')
 for n,h in EXPECTED.items():
  if sha(root/n)!=h:raise ValueError('Source SHA mismatch: '+n)
 targets={(x['geometry'],int(x['box_id'])) for x in rows(root/'target_cells.csv')}
 if len(targets)!=200:raise ValueError('Wrong A100 target set')
 for g in GEOMS:
  if sum(k[0]==g for k in targets)!=50:raise ValueError('Target count invalid')
 tables={}
 for g in GEOMS:
  ncell=900 if g.startswith('p3_') else 1200
  train=[[] for _ in range(ncell)];extra=[[] for _ in range(ncell)]
  for p,buf,per in [(root/f'trials_{g}.csv',train,20),(root/f'trials_extra_{g}.csv',extra,80)]:
   if p.name.startswith('trials_') and not p.name.startswith('trials_extra_'):
    if sha(p)!=m['files_sha256'][p.name]:raise ValueError('Old20 source SHA mismatch: '+p.name)
   seen=set()
   with open(p,newline='') as f:
    for r in csv.DictReader(f):
     b=int(r['box_id']);sample=int(r['sample']);key=(b,sample)
     if key in seen or int(r['capped'])!=0 or int(r['successful'])!=1:raise ValueError('Invalid microscopic sample')
     if r['configuration']!=g or r['original_geometry']!=g:raise ValueError('Geometric source identity mismatch')
     seen.add(key)
     buf[b].append((sample,tuple(int(r[k]) for k in COLUMNS)))
   if len(seen)!=(ncell if per==20 else 50)*per:raise ValueError('Incorrect trial cardinality')
  for box in range(ncell):
   if len(train[box])!=20 or len(extra[box])!=(80 if (g,box) in targets else 0):raise ValueError('Wrong training cells')
   if {x[0] for x in train[box]} & {x[0] for x in extra[box]}:raise ValueError('A20/A100 not independent samples')
  a=np.asarray([[x[1] for x in sorted(c)] for c in train],dtype=np.int64)
  b={i:np.asarray([x[1] for x in sorted(extra[i])],dtype=np.int64) for i in range(ncell) if extra[i]}
  if np.any(a[...,1]+a[...,2]>a[...,0]) or np.any(a[...,3]+a[...,4]>a[...,0]):raise ValueError('Ion sector counts not physical')
  tables[g]=(a,b)
 return tables,targets,rows(root/'native_muon_crossings_input.csv'),rows(root/'muon_reference_heed_primaries.csv'),rows(root/'muon_reference_crossing.csv'),rows(root/'muon_A20_vs_targeted_A100.csv'),{'input':'SHA-pinned 84k old20 + SHA-pinned added16k full Garfield and 120 real muons','A20_sha':A20_SHA,'A100_sha':A100_SHA}

def load_native():
 """Read the native WSL SHA-sealed source archive via existing verified code."""
 from pi_A100_crn_hierarchical_check import load_native as source_loader, assert_pool
 original,target,mu,prim,cross,old,provenance,source_root=source_loader()
 assert_pool(original,target)
 tables={}
 for g in GEOMS:
  n=900 if g.startswith('p3_') else 1200
  a=np.asarray([original[g][b][:20] for b in range(n)],dtype=np.int64)
  extra={b:np.asarray(original[g][b][20:],dtype=np.int64) for b in range(n) if len(original[g][b])==100}
  if a.shape!=(n,20,5) or len(extra)!=50 or any(x.shape!=(80,5) for x in extra.values()):
   raise ValueError('Native source dimension or old20 plus new80 physics corrupted')
  tables[g]=(a,extra)
 return tables,target,mu,prim,cross,old,provenance,source_root


def make_records(tables,targets,mu,prim,cross,original):
 from collections import Counter,defaultdict
 pr=defaultdict(list)
 for p in prim:pr[p['crossing_id']].append(p)
 ref={(r['source|chunk|event|hit_index|track_id|station|det_id'],r['configuration']):r for r in cross}
 cmp={(r['crossing_id'],r['configuration']):r for r in original}
 if len(mu)!=120 or len(prim)!=5881 or len(ref)!=480 or len(cmp)!=480:raise ValueError('Wrong source cardinality')
 out=[]
 for m in mu:
  mid='|'.join(str(m[x]) for x in ('source','chunk','event','hit_index','track_id','station','det_id'))
  if mid not in pr:raise ValueError('Missing original TrackHeed primaries')
  p=pr[mid]
  for g in GEOMS:
   pitch=3. if g.startswith('p3_') else 4.
   nw=math.floor((307.2 if int(m['station'])==0 else 512.)/pitch+0.5)
   wirephase=pitch*0.5 if nw%2==0 else 0.
   def rnd(x):return math.floor(x+.5) if x>=0 else math.ceil(x-.5)
   wire=m['u_in_mm'];wire=float(wire)
   shift=wire-(rnd((wire-wirephase)/pitch)*pitch+wirephase)
   cells=[]
   for e in p:
    u=10.*float(e['u_rel_cm'])+shift
    u=u-rnd(u/pitch)*pitch
    v=10.*float(e['v_abs_cm'])
    if not (-3.-1e-8<=v<=3.+1e-8):raise ValueError('Outside gas')
    iu=max(0,min(30 if g.startswith('p3_') else 40,math.floor((u+pitch/2.)/.1)))
    iu=min(iu,(30 if g.startswith('p3_') else 40)-1)
    iv=max(0,min(29,math.floor((v+3.)/.2)))
    cells.append(int(iu*30+iv))
   counts=tuple(sorted(Counter(cells).items()))
   t=ref[mid,g]
   o=cmp[mid,g]
   birth=int(t['ion_births'])
   if birth!=int(o['full_reference_ion_births']) or len(p)!=int(o['n_real_Heed_electrons']):raise ValueError('Full Garfield truth not matched')
   selected=sum(n for box,n in counts if (g,box) in targets)
   if selected!=int(o['n_electrons_in_selected_A100_cells']):raise ValueError('Wrong Heed spatial mapping')
   arr=tables[g][0]
   if any(not 0<=box<len(arr) for box,n in counts):raise ValueError('Invalid Heed cell')
   out.append(dict(parent=mid,geom=g,counts=counts,n_primary=len(p),n_target=selected,
                   ref_gain=birth,ref_up=int(t['ion_up']),ref_down=int(t['ion_down']),
                   ref_top=int(t['ion_top']),ref_bottom=int(t['ion_bottom'])))
 if len(out)!=480 or sum(x['n_target']==0 for x in out)!=158:raise ValueError('120 parent records/geometry counts changed')
 return out

def donors_for(g,cell,guard=True,radius=2):
 """A20 neighbor donors; conservative same-ion-hemisphere and wire safeguards.

 Prevent pooling across anode upper/lower half (v=0) for ion-sector response;
 near central anode or optional field wire (~0.35 mm center-to-metal), do
 not pool across bins. For near-axis |v|<=0.6mm also do not mix u signs.
 These geometric guard rules are predeclared, not tuned to muon truth.
 """
 nu=30 if g.startswith('p3_') else 40;pitch=3. if g.startswith('p3_') else 4.
 i,j=divmod(cell,30)
 u=-pitch/2+(i+.5)*.1;v=-3+(j+.5)*.2
 fw=('FW75' in g)
 def wire_dist(x,y):
  d=math.hypot(x,y)
  if fw:d=min(d,math.hypot(x-pitch/2,y),math.hypot(x+pitch/2,y))
  return d
 near=wire_dist(u,v)<0.35
 result=[]
 if guard and near:return (np.array([],dtype=np.int32),np.array([],dtype=float))
 for du in range(-radius,radius+1):
  for dv in range(-radius,radius+1):
   if du==0 and dv==0:continue
   x,y=i+du,j+dv
   if not(0<=x<nu and 0<=y<30):continue
   un=-pitch/2+(x+.5)*.1;vn=-3+(y+.5)*.2
   if guard:
    if (v>0)!=(vn>0):continue
    if (abs(v)<=.6 or abs(vn)<=.6) and (u>0)!=(un>0):continue
    if wire_dist(un,vn)<.35:continue
   result.append((x*30+y,math.exp(-.5*(du*du+dv*dv))))
 if not result:return np.array([],dtype=np.int32),np.array([],dtype=float)
 ids,weights=zip(*result)
 return np.asarray(ids,dtype=np.int32),np.cumsum(np.asarray(weights,dtype=float)/sum(weights))

def compile_neighbors(tables,guard):
 return {g:[donors_for(g,i,guard) for i in range(len(tables[g][0]))] for g in GEOMS}

def evaluate_one(r,table,neighbors,replicas,seed):
 """Paired common-random-number per-electron sampling with exact A20 vs A100
 identity in all unmodified cells. Independent draws across primaries, shared
 random streams across predictors. A100 target uses 20% original20 + 80% new80.
 """
 a,extras=table
 rng=np.random.default_rng(int.from_bytes(hashlib.sha256(f'MU120_POOL_V1|{seed}|{r["parent"]}|{r["geom"]}'.encode()).digest()[:8],'big'))
 v=np.zeros((len(MODES),replicas,5),np.int64)
 n_guard=0;n_unguard=0;n_smooth=0
 for box,n in r['counts']:
  draw=rng.random((5,replicas,n))
  idx20=np.minimum(19,(draw[0]*20).astype(np.intp))
  o=a[box,idx20]
  if box in extras:
   ix80=np.minimum(79,(draw[0]*80).astype(np.intp))
   new=np.where((draw[1]<.2)[...,None],o,extras[box][ix80])
  else:new=o
  v[0]+=o.sum(axis=1);v[1]+=new.sum(axis=1)
  for x,mix in ((2,'guarded'),(4,'unguarded')):
   ids,cdf=neighbors[mix][r['geom']][box]
   if not len(ids):
    v[x]+=o.sum(axis=1);v[x+1]+=new.sum(axis=1)
    continue
   other=(draw[2]>=.25)
   nb=ids[np.searchsorted(cdf,draw[3],side='right').clip(0,len(ids)-1)]
   donated=a[nb, np.minimum(19,(draw[4]*20).astype(np.intp))]
   v[x]+=np.where(other[...,None],donated,o).sum(axis=1)
   v[x+1]+=np.where(other[...,None],donated,new).sum(axis=1)
   if mix=='guarded':n_guard+=n;n_smooth+=n
   else:n_unguard+=n
 metrics={}
 for i,mode in enumerate(MODES):
  g=v[i,:,0].astype(float)
  ref=r['ref_gain'];lo,hi=np.quantile(g,[.05,.95]);mean=np.mean(g)
  metrics.update({f'{mode}_mean':float(mean),f'{mode}_p05':float(lo),f'{mode}_p95':float(hi),
                  f'{mode}_covered':int(lo<=ref<=hi),f'{mode}_abs_mean_error':float(abs(mean-ref)),
                  f'{mode}_width':float(hi-lo)})
  for label,a_idx,b_idx in [('A_hemi',1,2),('A_TB',3,4)]:
   den=v[i,:,a_idx]+v[i,:,b_idx]
   vals=(v[i,:,a_idx]-v[i,:,b_idx])[den>0]/den[den>0]
   ar=r['ref_up'] if label=='A_hemi' else r['ref_top']
   br=r['ref_down'] if label=='A_hemi' else r['ref_bottom']
   if ar+br and len(vals):
    l,h=np.quantile(vals,[.05,.95]);metrics[f'{mode}_{label}_covered']=int(l<=(ar-br)/(ar+br)<=h)
   else: metrics[f'{mode}_{label}_covered']=''
 return dict(parent=r['parent'],geometry=r['geom'],n_primary=r['n_primary'],n_target=r['n_target'],
             true_ion_births=ref,n_guard_poolable=n_guard,n_unguarded_poolable=n_unguard,**metrics)

def pct_bootstrap_delta(per,nboot=3000,seed=45):
 rng=np.random.default_rng(seed)
 original=np.array([[a[f'{m}_covered'] for m in MODES] for a in per],dtype=float)
 # Each physical parent contributes exactly four geometry rows, one per GEOMS.
 if len(per)!=480 or not all(per[4*i+j]['geometry']==GEOMS[j] for i in range(120) for j in range(4)):
  raise ValueError('Parent/geometry grouping lost')
 differences=(original[:,2:]-np.tile(original[:,0:2],(1,2))) # paired [guard A20 vs A20, guard A100 vs A100, unguard A20 vs A20, unguard A100 vs A100]
 differences=differences.reshape(120,4,4).mean(axis=1)
 boot=differences[rng.integers(120,size=(nboot,120))].mean(axis=1)
 return {name: {'delta_pp':float(100*differences[:,i].mean()),'parent_bootstrap95_pp':np.quantile(100*boot[:,i],[.025,.975]).tolist()} for i,name in enumerate(('guarded_A20_minus_A20','guarded_A100_minus_A100','unguarded_A20_minus_A20','unguarded_A100_minus_A100'))}

def main():
 ap=argparse.ArgumentParser(description=__doc__)
 ap.add_argument('mode',choices=('run','selftest'))
 ap.add_argument('--published-dir',type=Path,default=Path('/mnt/data'))
 ap.add_argument('--native',action='store_true',help='Use SHA-sealed original WSL native chunks and archived 120 Geant4 reference')
 ap.add_argument('--output',type=Path,default=None,help='New directory only; cannot overwrite old results')
 ap.add_argument('--reps',type=int,default=1500)
 ap.add_argument('--seed',type=int,default=20261031)
 args=ap.parse_args()
 if args.mode=='selftest':
  for g in GEOMS:
   a,b=donors_for(g,(15 if g.startswith('p3_') else 20)*30+15,True)
   assert not len(a)  # protected wire center
   assert all((k%30)>=15 for k in donors_for(g,15*30+24,True)[0])
  v=np.zeros((20,5),int);v[:,0]=np.arange(20)
  t={'p3_noFW':(np.array([v]),{})}
  demo={'parent':'demo','geom':'p3_noFW','counts':((0,3),),'n_primary':3,'n_target':0,'ref_gain':24,'ref_up':0,'ref_down':0,'ref_top':0,'ref_bottom':0}
  neigh={s:{'p3_noFW':[(np.array([],int),np.array([],float))]} for s in ('guarded','unguarded')}
  out=evaluate_one(demo,t['p3_noFW'],neigh,500,29)
  assert all(out[f'{m}_mean']==out['A20_mean'] and out[f'{m}_covered']==out['A20_covered'] for m in MODES)
  d=pct_bootstrap_delta([dict(geometry=GEOMS[i%4],**{f'{m}_covered': 1 for m in MODES}) for i in range(480)],nboot=10)
  assert all(t['delta_pp']==0 for t in d.values())
  print('PASS synthetic protected-cell A20/A100/guarded CRN identity; no detector physics')
  return
 if not 300<=args.reps<=10000:raise ValueError('Use 300..10000 conditional draws per real muon geometry')
 if args.native:
  table,targets,mu,pri,cross,compare,provenance,source_root=load_native()
  output=args.output or source_root/'compare_frozen_120_real_muons_guarded_spatial_pooling_v1'
 else:
  output=args.output
  if output is None:ap.error('Published archived mode requires explicit --output to avoid overwriting original archives')
  table,targets,mu,pri,cross,compare,provenance=load_archive(args.published_dir)
 if output.exists():raise FileExistsError(f'Cannot overwrite output: {output}')
 rec=make_records(table,targets,mu,pri,cross,compare)
 neighbors={'guarded':compile_neighbors(table,True),'unguarded':compile_neighbors(table,False)}
 result=[]
 for k,r in enumerate(rec):
  z=evaluate_one(r,table[r['geom']],neighbors,args.reps,args.seed)
  if r['n_target']==0:
   for name in ('A20','guarded_A20','unguarded_A20'):
    corresponding={'A20':'A100','guarded_A20':'guarded_A100','unguarded_A20':'unguarded_A100'}[name]
    for fld in ('mean','p05','p95','covered'):
     if z[f'{name}_{fld}']!=z[f'{corresponding}_{fld}']:
      raise AssertionError('Common random numbers changed unchanged A20/A100 response '+str((r['parent'],r['geom'],name,fld)))
  result.append(z)
  if (k+1)%40==0:print(f'PASS paired pooled LUT {k+1}/480, {args.reps} replicas',flush=True)
 cover={g:{m:sum(v['geometry']==g and v[f'{m}_covered'] for v in result)/120 for m in MODES} for g in GEOMS}
 overall={m:sum(v[f'{m}_covered'] for v in result)/480 for m in MODES}
 mae={g:{m:float(np.mean([v[f'{m}_abs_mean_error'] for v in result if v['geometry']==g])) for m in MODES} for g in GEOMS}
 stat={'schema':'NA60_MWPC_CRN_GUARDED_A20_A100_SPATIAL_POOLING_MU120_V1','status':'complete',
       'n_parent_muons':120,'n_correlated_geometry_tests':480,'reps_per_test':args.reps,
       'conditional_gain_90PI_coverage_by_geometry':cover,'conditional_gain_90PI_overall':overall,
       'mean_abs_error_of_predictive_mean_by_geometry':mae,
       'paired_parent_bootstrap90PI_coverage_deltas':pct_bootstrap_delta(result),
       'n_Heed_electron_starts_guard_poolable':sum(x['n_guard_poolable'] for x in result),
       'n_Heed_electron_starts_unguarded_poolable':sum(x['n_unguarded_poolable'] for x in result),
       'configuration':{'radius_bins':2,'local_probability':.25,'neighbor_probability':.75,
        'kernel':'exp[-(du^2+dv^2)/2] on original A20 neighbor cells only; local A20 or A100 independent Garfield samples',
        'guard':'no crossing v=0; preserve u sign near wire axis; suppress pooling within 0.35 mm of central anode or optional FW wire (measured at bin centers); never wrap edges'},
       'caveats':('Only 120 previously used parent muons, 4 correlated geometries; source-conditioned A100 cells and pooling hyperparameters selected from same high-variance cells, not unseen external validation. ' 
        'Each true full Garfield muon observation is a single stochastic realization. This is conditional empirical resampling only, not cell-training bootstrap, and does not resolve independent finite-library uncertainty. ' 
        'Within-cell exact start-position mismatch and ion-birth spatial charge clouds remain unresolved. No ion drift / U-V induced strip signals.'),
       'provenance':provenance}
 output.mkdir(parents=True)
 per=output/'per_muon_6models.csv';summary=output/'summary.json'
 with per.open('w',newline='') as f:
  writer=csv.DictWriter(f,fieldnames=list(result[0]));writer.writeheader();writer.writerows(result)
 summary.write_text(json.dumps(stat,indent=2,sort_keys=True)+'\n')
 (output/'manifest.json').write_text(json.dumps({'status':'complete','source_A20':A20_SHA,'source_A100':A100_SHA,'files_sha256':{per.name:sha(per),summary.name:sha(summary)}},indent=2)+'\n')
 print(json.dumps({k:stat[k] for k in ['conditional_gain_90PI_overall','paired_parent_bootstrap90PI_coverage_deltas','n_Heed_electron_starts_guard_poolable']},indent=2),flush=True)
 print('PASS source-conditioned matched 120-muon ×4 geometry spatial-pooling LUT evaluation; original physics unchanged',output,flush=True)

if __name__=='__main__':main()