#!/usr/bin/env python3
"""Read-only candidate Polya/hurdle-Gamma 90%-interval holdout on targeted 200 cells.

Train original 20 empirical Garfield avalanches per selected cell; independently
validate using 80 new microscopic avalanches per cell. Exactly 16,000 test events.
For simplicity and interpretable baseline, compare:
 - nonparametric empirical 5th/95th percentile from 20 trials (NumPy linear)
 - a two-stage probability model: atom at zero (Beta(1/2,1/2) smoothing),
   positive gain drawn from Gamma (Polya continuous-limit), MLE shape & scale.

These are SINGLE-ELECTRON response interval tests, not predictions for full muon
ion yield or electrical readout; not a validated alternative detector simulator.
"""
import argparse,csv,hashlib,json,math
from collections import defaultdict
from pathlib import Path
import numpy as np
from scipy.stats import gamma

GEOMS=('p3_noFW','p3_FW75','p4_noFW','p4_FW75')
SHA_A100={
 'trials_extra_p3_noFW.csv':'4922e5930dfa3a1ffa7f6dfe76f1408e61e76c008fa6c040fe14d64ec0e74f67',
 'trials_extra_p3_FW75.csv':'a9b7c2cfafdc5d910b0b5b779468ecfdd166ec4005331855a8c7a793a02c4ab4',
 'trials_extra_p4_noFW.csv':'c24f411f28e96197924d51bbb6906de6723a9eb6a93372648e24e4b3fd1576ac',
 'trials_extra_p4_FW75.csv':'b7e2a8d7d184509f0e7b98739067c81a20b64f74549fc2a6a68bdaa1912a4e22',
 'target_cells.csv':'1f3320e0b336266eec5aa1fe41acacaf25baf3d421b59dbdde5788f7df31b8b6',
}
A20_MANIFEST='8ba639e3f082192d716d328775e13743ca72f9035ebfe703b1f722a202fb006c'
A100_MANIFEST='2a8520c61dc38dc2b32aa31319c3665e36d9ed86f9119e8f3737f32849d7da82'

def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()

def verified(p,s):
 if not p.exists() or sha(p)!=s:raise ValueError(f'SHA mismatch/missing {p}')

def read(p):
 with open(p,newline='') as f:return list(csv.DictReader(f))

def fit_gamma(x):
 x=x[x>0]
 if len(x)<3:raise ValueError('Insufficient positive avalanche training events')
 a,loc,scale=gamma.fit(x,floc=0)
 if not (math.isfinite(a) and math.isfinite(scale) and a>0 and scale>0 and loc==0):
  raise ValueError('Invalid fitted conditional positive gamma')
 return a,scale

def bounds_gamma(train,zero_prior=.5):
 n=len(train);zero=int(np.sum(train==0))
 p0=(zero+zero_prior)/(n+2*zero_prior)
 a,scale=fit_gamma(train)
 def q(p):
  if p<=p0:return 0.
  return float(gamma.ppf((p-p0)/(1-p0),a,loc=0,scale=scale))
 return p0,a,scale,q(.05),q(.95)

def main():
 ap=argparse.ArgumentParser(description=__doc__)
 ap.add_argument('--old-dir',default='/mnt/data')
 ap.add_argument('--new-dir',default='/mnt/data')
 ap.add_argument('--out-dir',default='/mnt/data/MWPC_A100_POLYA_20261031')
 x=ap.parse_args();old=Path(x.old_dir);new=Path(x.new_dir);out=Path(x.out_dir)
 if out.exists():raise FileExistsError('Existing immutable output: '+str(out))
 verified(old/'manifest.json',A20_MANIFEST)
 original=json.loads((old/'manifest.json').read_text())
 assert original['status']=='complete' and original['n_completed_avalanches']==84000
 if (new/'manifest.json').is_file() and sha(new/'manifest.json')!=A20_MANIFEST:
  verified(new/'manifest.json',A100_MANIFEST)
  m=json.loads((new/'manifest.json').read_text());assert m['status']=='complete' and m['new_full_Garfield_avalanches']==16000
 for name,d in SHA_A100.items():verified(new/name,d)
 target={(r['geometry'],int(r['box_id'])) for r in read(new/'target_cells.csv')}
 assert len(target)==200
 res=[]
 for g in GEOMS:
  p=old/f'trials_{g}.csv';verified(p,original['files_sha256'][p.name])
  oldg=defaultdict(list);newg=defaultdict(list);s_old=defaultdict(set);s_new=defaultdict(set)
  rs=read(p);assert len(rs)==(18000 if g.startswith('p3_') else 24000)
  for r in rs:
   b=int(r['box_id'])
   if (g,b) in target:
    oldg[b].append(int(r['ion_births']));s_old[b].add(int(r['sample']))
  for r in read(new/f'trials_extra_{g}.csv'):
   b=int(r['box_id']);assert (g,b) in target
   newg[b].append(int(r['ion_births']));s_new[b].add(int(r['sample']))
  assert len(oldg)==len(newg)==50
  for b in oldg:
   train=np.asarray(oldg[b],dtype=float);test=np.asarray(newg[b],dtype=float)
   assert len(train)==20 and len(test)==80 and not(s_old[b]&s_new[b])
   emp_lo,emp_hi=np.quantile(train,[.05,.95])
   p0,a,scale,pol_lo,pol_hi=bounds_gamma(train)
   res.append(dict(geometry=g,box_id=b,training_n=20,test_n=80,
       p0_estimate=p0,polya_shape=a,polya_scale=scale,
       old20_mean_gain=float(train.mean()),new80_mean_gain=float(test.mean()),
       old20_zero_probability=float((train==0).mean()),new80_zero_probability=float((test==0).mean()),
       empirical_lo=emp_lo,empirical_hi=emp_hi,polya_lo=pol_lo,polya_hi=pol_hi,
       empirical_test_below=int(np.sum(test<emp_lo)),empirical_test_above=int(np.sum(test>emp_hi)),
       empirical_test_covered=int(np.sum((test>=emp_lo)&(test<=emp_hi))),
       polya_test_below=int(np.sum(test<pol_lo)),polya_test_above=int(np.sum(test>pol_hi)),
       polya_test_covered=int(np.sum((test>=pol_lo)&(test<=pol_hi))),
       empirical_width=float(emp_hi-emp_lo),polya_width=float(pol_hi-pol_lo)))
 gsummary={}
 for g in GEOMS:
  r=[x for x in res if x['geometry']==g]
  gsummary[g]={
   'n_selected_cells':50,'n_new_independent_tests':4000,
   'empirical90_coverage':sum(a['empirical_test_covered'] for a in r)/4000,
   'polya90_coverage':sum(a['polya_test_covered'] for a in r)/4000,
   'empirical_tail_below':sum(a['empirical_test_below'] for a in r)/4000,
   'empirical_tail_above':sum(a['empirical_test_above'] for a in r)/4000,
   'polya_tail_below':sum(a['polya_test_below'] for a in r)/4000,
   'polya_tail_above':sum(a['polya_test_above'] for a in r)/4000,
   'empirical_avg_width':sum(a['empirical_width'] for a in r)/50,
   'polya_avg_width':sum(a['polya_width'] for a in r)/50,
   'median_polya_shape':float(np.median([a['polya_shape'] for a in r])),
  }
 rnd=np.random.default_rng(20261031)
 d=np.array([[(r['polya_test_covered']-r['empirical_test_covered'])/80,
               r['polya_test_covered']/80,r['empirical_test_covered']/80] for r in res])
 ci=np.quantile([d[rnd.integers(200,size=200),:].mean(axis=0) for _ in range(4000)], [.025,.975], axis=0)
 report={'schema':'NA60_MWPC_A20_POLYA_NEW80_INDEPENDENT_HOLDOUT_V1','status':'complete','total_new_physics_holdout':16000,
 'method':'Original 20 full Garfield avalanches per selected 200 cells train empirical quantiles or zero-atom plus fitted continuous-Gamma/Polya, new 80 are independent test. 90% nominal intervals single-electron ion births.',
 'by_geometry':gsummary,
 'overall_empirical_coverage':sum(a['empirical_test_covered'] for a in res)/16000,
 'overall_polya_coverage':sum(a['polya_test_covered'] for a in res)/16000,
 'overall_polya_minus_empirical_coverage_cluster95':ci[:,0].tolist(),
 'overall_polya_coverage_cluster95':ci[:,1].tolist(),
 'overall_empirical_coverage_cluster95':ci[:,2].tolist(),
 'caveats':'Cells preselected using old20 high uncertainty, so results apply only to these cells. New80 independent of selection ensures fair heldout training comparison within selected cells. Pólya/Gamma continuous positive tail + separate Bernoulli zero mass; no explicit field or transport; fit to 20 avalanches fragile. Single-electron ion births coverage is NOT full muon event coverage. No ion charge weighting, ion drift, B-field or strip signal.'}
 out.mkdir(parents=True)
 with (out/'per_cell.csv').open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(res[0]));w.writeheader();w.writerows(res)
 (out/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
 (out/'manifest.json').write_text(json.dumps({'status':'complete','input_sha256':{'old_original_A20':A20_MANIFEST,'new_original_A100':A100_MANIFEST},'files_sha256':{p.name:sha(p) for p in out.iterdir() if p.is_file()}},indent=2)+'\n')
 print(json.dumps(report,indent=2))
 print('PASS Pólya/gamma vs A20 empirical independent new80 holdout; 16,000 full Garfield physically computed trials REUSED; no new physics')

if __name__=='__main__':main()