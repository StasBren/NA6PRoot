#!/usr/bin/env python3
"""Read-only within-cell MWPC avalanche position-dependence diagnostic.

Checks SHA of original A20 (84k) and targeted new A80 (16k); analyses only
new independent 80 trials per selected cell (chosen using ORIGINAL A20).
Fourfold cross-validation compares constant versus linear/quadratic dependence
on within-cell position for log(1+ion births) and for zero-birth probability.

No Garfield/Heed/Geant4 invocation. The test does NOT calibrate muon response;
it tests whether subcell coordinate information improves held-out predictions.
Usage:
  python pi_A100_subcell_spatial_diagnostic.py \\
    --old-dir test_runs/mwpc_stage_c/pi_lut500_avalanche_overnight_v1/run_20percell_firstpass_v1 \\
    --new-dir test_runs/mwpc_stage_c/pi_A100_targeted_50cells_per_geometry_v1 \\
    --out-dir test_runs/mwpc_stage_c/pi_A100_targeted_50cells_per_geometry_v1/subcell_spatial_diagnostic_v1
"""
import argparse, csv, hashlib, json, math, statistics
from collections import defaultdict
from pathlib import Path
import numpy as np

GEOMS=('p3_noFW','p3_FW75','p4_noFW','p4_FW75')
SHA_A100_FILES={
 'trials_extra_p3_noFW.csv':'4922e5930dfa3a1ffa7f6dfe76f1408e61e76c008fa6c040fe14d64ec0e74f67',
 'trials_extra_p3_FW75.csv':'a9b7c2cfafdc5d910b0b5b779468ecfdd166ec4005331855a8c7a793a02c4ab4',
 'trials_extra_p4_noFW.csv':'c24f411f28e96197924d51bbb6906de6723a9eb6a93372648e24e4b3fd1576ac',
 'trials_extra_p4_FW75.csv':'b7e2a8d7d184509f0e7b98739067c81a20b64f74549fc2a6a68bdaa1912a4e22',
 'target_cells.csv':'1f3320e0b336266eec5aa1fe41acacaf25baf3d421b59dbdde5788f7df31b8b6',
}
EXPECTED_NEW_MANIFEST='2a8520c61dc38dc2b32aa31319c3665e36d9ed86f9119e8f3737f32849d7da82'
EXPECTED_OLD_MANIFEST='8ba639e3f082192d716d328775e13743ca72f9035ebfe703b1f722a202fb006c'

def digest(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()

def must_sha(p,expected):
 assert p.is_file(),f'Missing {p}'
 seen=digest(p)
 if seen!=expected:raise ValueError(f'SHA mismatch {p}: {seen} expected {expected}')

def rows(p):
 with open(p,newline='') as f: return list(csv.DictReader(f))

def load(old,new):
 oldm=old/'manifest.json'
 must_sha(oldm,EXPECTED_OLD_MANIFEST)
 m=json.loads(oldm.read_text())
 assert m['status']=='complete' and m['n_completed_avalanches']==84000
 # Optional physical manifest present on WSL, omitted from flat published archive.
 newm=new/'manifest.json'
 if newm.is_file() and digest(newm)!=EXPECTED_OLD_MANIFEST:
  must_sha(newm,EXPECTED_NEW_MANIFEST)
  j=json.loads(newm.read_text()); assert j['status']=='complete' and j['new_full_Garfield_avalanches']==16000
  for name in SHA_A100_FILES: assert j['files_sha256'][name]==SHA_A100_FILES[name]
 for name,expected in SHA_A100_FILES.items():must_sha(new/name,expected)
 tr=rows(new/'target_cells.csv')
 target={(r['geometry'],int(r['box_id'])) for r in tr}
 assert len(target)==200 and {g:sum(a==g for a,b in target) for g in GEOMS}==dict.fromkeys(GEOMS,50)
 selected={}
 old_selected={}
 for g in GEOMS:
  p=old/f'trials_{g}.csv';must_sha(p,m['files_sha256'][p.name]);original=rows(p)
  n=18000 if g.startswith('p3') else 24000
  assert len(original)==n
  old_subset=defaultdict(list)
  for r in original:
   b=int(r['box_id'])
   if (g,b) in target:old_subset[b].append(r)
  assert len(old_subset)==50 and all(len(a)==20 for a in old_subset.values())
  rs=rows(new/f'trials_extra_{g}.csv')
  assert len(rs)==4000
  gnew=defaultdict(list)
  for r in rs:
   b=int(r['box_id'])
   if (g,b) not in target:raise ValueError('Unexpected new avalanche cell')
   gnew[b].append(r)
  assert len(gnew)==50 and all(len(a)==80 for a in gnew.values())
  for b in gnew:
   a=old_subset[b];bb=gnew[b]
   if {int(t['sample']) for t in a}&{int(t['sample']) for t in bb}:raise ValueError('Duplicate D500 original start')
   for q in a+bb:
    assert int(q['box_id'])==b and q['configuration']==g
    assert not int(q['capped']) and int(q['successful'])==1
    iu=int(q['iu']);iv=int(q['iv']);assert iu*30+iv==b
    # Check geometry and center within exact 0.1×0.2 mm bin.
    pitch=3 if g.startswith('p3') else 4
    uc=-pitch/2+iu*.1+.05
    vc=-3+iv*.2+.1
    du=float(q['seed_u_mm'])-uc
    dv=float(q['seed_v_mm'])-vc
    if not (-.0500001 <= du <= .0500001 and -.1000001<=dv<=.1000001):
     raise ValueError('Original drift-start coordinate leaves assigned cell')
  selected[g]=gnew
  old_selected[g]=old_subset
 return selected,old_selected

def design(rs,features):
 du=[];dv=[]
 for r in rs:
  g=r['configuration'];pitch=3 if g.startswith('p3') else 4
  iu=int(r['iu']);iv=int(r['iv'])
  du.append((float(r['seed_u_mm'])-(-pitch/2+iu*.1+.05))/.05)
  dv.append((float(r['seed_v_mm'])-(-3+iv*.2+.1))/.1)
 u=np.asarray(du);v=np.asarray(dv)
 if features==1: return np.column_stack([u,v])
 if features==2: return np.column_stack([u,v,u*u-1/3,v*v-1/3,u*v])
 raise ValueError

def fit(X,y,alpha):
 # Ridge coefficient shrinkage, intercept NOT regularized.
 m=X.mean(axis=0);my=y.mean()
 Z=X-m
 beta=np.linalg.solve(Z.T@Z+alpha*np.eye(X.shape[1]),Z.T@(y-my))
 return my-m@beta,beta

def validate_cv(rs,seed=0,folds=4):
 n=len(rs);assert n==80
 # The 80 new responses are independent of target cell selection. Keep old20 out of training/test.
 rng=np.random.default_rng(seed)
 ids=rng.permutation(n);fold_id=np.empty(n,dtype=int)
 fold_id[ids]=np.arange(n)%folds
 a=design(rs,1);b=design(rs,2)
 gain=np.array([int(q['ion_births']) for q in rs],dtype=float)
 l=np.log1p(gain);z=(gain==0).astype(float)
 results={k:0. for k in ('base_log','lin_log','quad_log','base_zero','lin_zero','quad_zero','base_abs_gain','lin_abs_gain')}
 for f in range(folds):
  te=fold_id==f;tr=~te
  # log transform handles heavy right tail, zero events remain explicit.
  for targ,prefix in [(l,'log'),(z,'zero')]:
   base=np.repeat(targ[tr].mean(),te.sum())
   if prefix=='zero':
    # Regression is deliberately clipped to [0,1] probability range for Brier score.
    base=np.clip(base,0,1)
   results['base_'+prefix]+=np.sum((targ[te]-base)**2)
   for name,X,lam in [('lin',a,5.0),('quad',b,10.0)]:
    intercept,beta=fit(X[tr],targ[tr],lam)
    pred=intercept+X[te]@beta
    if prefix=='zero':pred=np.clip(pred,0,1)
    results[name+'_'+prefix]+=np.sum((targ[te]-pred)**2)
  # Physical ion-yield baseline and linear predictor, direct squared loss (tails sensitive).
  raw=np.log1p(gain[tr]).mean() # intentionally not used to predict gain mean
  constant=gain[tr].mean(); intercept,beta=fit(a[tr],gain[tr],5.0)
  results['base_abs_gain']+=np.sum(np.abs(gain[te]-constant))
  results['lin_abs_gain']+=np.sum(np.abs(gain[te]-np.maximum(intercept+a[te]@beta,0)))
 return {k:float(v/n) for k,v in results.items()}

def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--old-dir',default='/mnt/data')
 p.add_argument('--new-dir',default='/mnt/data')
 p.add_argument('--out-dir',default='/mnt/data/MWPC_A100_SUBCELL_20261031')
 args=p.parse_args()
 old=Path(args.old_dir);new=Path(args.new_dir);out=Path(args.out_dir)
 if out.exists():raise FileExistsError(f'Output exists, protect it: {out}')
 newer,old20=load(old,new)
 result=[]
 for geom in GEOMS:
  for b,rs in sorted(newer[geom].items()):
   seed=int.from_bytes(hashlib.sha256(f'NA60_SUBCELL_{geom}_{b}'.encode()).digest()[:8],'big')
   cv=validate_cv(rs,seed)
   lg=np.array([int(q['ion_births']) for q in rs])
   ol=np.array([int(q['ion_births']) for q in old20[geom][b]])
   result.append(dict(geometry=geom,box_id=b,n_old=20,n_new=80,
     old_mean_gain=float(ol.mean()),new_mean_gain=float(lg.mean()),
     old_p_zero=float(np.mean(ol==0)),new_p_zero=float(np.mean(lg==0)),**cv))
 groups={}
 for geom in GEOMS:
  xs=[r for r in result if r['geometry']==geom]
  x={k:float(sum(r[k] for r in xs)) for k in ['base_log','lin_log','quad_log','base_zero','lin_zero','quad_zero','base_abs_gain','lin_abs_gain']}
  groups[geom]={
    'n_target_cells':len(xs),'n_unseen_from_selection_avalanches':sum(r['n_new'] for r in xs),
    'log_gain_linear_cv_mse_reduction_pct':100*(x['base_log']-x['lin_log'])/x['base_log'],
    'log_gain_quadratic_cv_mse_reduction_pct':100*(x['base_log']-x['quad_log'])/x['base_log'],
    'zero_birth_linear_cv_brier_reduction_pct':100*(x['base_zero']-x['lin_zero'])/x['base_zero'],
    'zero_birth_quadratic_cv_brier_reduction_pct':100*(x['base_zero']-x['quad_zero'])/x['base_zero'],
    'raw_gain_linear_cv_abs_error_reduction_pct':100*(x['base_abs_gain']-x['lin_abs_gain'])/x['base_abs_gain'],
    'cells_log_gain_linear_improved':sum(r['lin_log']<r['base_log'] for r in xs),
    'cells_log_gain_quadratic_improved':sum(r['quad_log']<r['base_log'] for r in xs),
    'cells_zero_linear_improved':sum(r['lin_zero']<r['base_zero'] for r in xs),
  }
 total={k:float(sum(r[k] for r in result)) for k in ['base_log','lin_log','quad_log','base_zero','lin_zero','quad_zero','base_abs_gain','lin_abs_gain']}
 all_summ={
  'n_target_cells':len(result),'n_new_independent_avalanches':len(result)*80,
  'cv_log_gain_linear_reduction_pct':100*(total['base_log']-total['lin_log'])/total['base_log'],
  'cv_log_gain_quad_reduction_pct':100*(total['base_log']-total['quad_log'])/total['base_log'],
  'cv_zero_linear_reduction_pct':100*(total['base_zero']-total['lin_zero'])/total['base_zero'],
  'cv_zero_quad_reduction_pct':100*(total['base_zero']-total['quad_zero'])/total['base_zero'],
  'cv_raw_gain_linear_absolute_reduction_pct':100*(total['base_abs_gain']-total['lin_abs_gain'])/total['base_abs_gain'],
 }
 # Cluster bootstrap at unit of independent selected cell. Heed source selection of cells is fixed.
 rng=np.random.default_rng(20261031);boots={}
 for metric,base,cand in [('log_linear','base_log','lin_log'),('log_quad','base_log','quad_log'),('zero_linear','base_zero','lin_zero'),('zero_quad','base_zero','quad_zero')]:
  arr=np.array([(r[base],r[cand]) for r in result]);n=len(arr)
  deltas=[]
  for _ in range(4000):
   ids=rng.integers(n,size=n);a,b=arr[ids].sum(axis=0)
   deltas.append(100*(a-b)/a)
  boots[metric]=[float(x) for x in np.quantile(deltas,[.025,.975])]
 summary={
  'schema':'NA60_MWPC_A100_SUBCELL_CV_V1','method':'Four-fold held-out CV on independent 80 new A100 avalanches per selected cell; old20 untouched and omitted from CV',
  'old_original_A20_manifest_sha256':EXPECTED_OLD_MANIFEST,
  'new_original_A100_manifest_sha256':EXPECTED_NEW_MANIFEST,
  'status':'complete','geometries':groups,'overall':all_summ,'bootstrap_95pct_across_target_cells':boots,
  'caveats':'Diagnostics only targeted 200 cells (selected from noisy original A20); independent new80 avoids direct ranking leakage. Position features within (0.1 x 0.2) mm bins. Cross-validation tests predictive point loss for log-gain and zero probability; does NOT validate full gain distribution/90% coverage, and does not establish interpolation validity for Heed position distributions. B=0, no ion drift/cathode signal.',
 }
 out.mkdir(parents=True)
 with (out/'per_cell.csv').open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(result[0]));w.writeheader();w.writerows(result)
 (out/'summary.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False)+'\n')
 (out/'manifest.json').write_text(json.dumps({'status':'complete','inputs_sha256':{'A20manifest':EXPECTED_OLD_MANIFEST,'A100manifest':EXPECTED_NEW_MANIFEST},'files_sha256':{p.name:digest(p) for p in out.iterdir() if p.is_file()}},indent=2)+'\n')
 print(json.dumps(summary,indent=2))
 print('PASS subcell 200 selected cells x 80 truly independent new avalanche starts, SHA verified; no new Garfield physics')

if __name__=='__main__': main()