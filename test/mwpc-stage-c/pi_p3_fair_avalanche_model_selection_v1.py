#!/usr/bin/env python3
"""NA60+ p=3 mm ONLY: fair A20-trained model comparison on unseen new80.

Competing models use the SAME original A20 source (20 microscopic trials/cell):
  local20                 = empirical 20 outcomes in requested cell
  guarded20               = frozen protected 25% local/75% 5x5 neighbor kernel
  geometry_hurdle_gamma   = zero+Gamma fitted only to ORIGINAL A20 of SAME geometry
No RKF or cathode calculation. All source checksums checked. No old files changed.
A100 new80 (50 selected cells/geometry) is the trial-level test, NOT training.
Selected cells are source-conditioned: this is development evaluation, not final
independent representative muon validation.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,math,statistics
from collections import defaultdict
from pathlib import Path
import numpy as np
from scipy.special import gammaln
from scipy.stats import gamma

HERE=Path(__file__).resolve().parents[2]
BASE=HERE/'test_runs/mwpc_stage_c'
OLD=BASE/'pi_lut500_avalanche_overnight_v1/run_20percell_firstpass_v1'
EXTRA=BASE/'pi_A100_targeted_50cells_per_geometry_v1'
OUTPUT=BASE/'pi_p3_fair_avalanche_model_selection_v1'
GEOMS=('p3_noFW','p3_FW75')
# Published independent A100 SHA values; sourced from frozen targeted50x4 manifest.
EXTRA_SHA={
 'trials_extra_p3_noFW.csv':'4922e5930dfa3a1ffa7f6dfe76f1408e61e76c008fa6c040fe14d64ec0e74f67',
 'trials_extra_p3_FW75.csv':'a9b7c2cfafdc5d910b0b5b779468ecfdd166ec4005331855a8c7a793a02c4ab4',
}
# Read exact manifest values first. These two strings may be wrong in legacy notes;
# actual WSL complete A100 manifest, when present, always governs extra data.
TARGET_SHA='1f3320e0b336266eec5aa1fe41acacaf25baf3d421b59dbdde5788f7df31b8b6'
MODES=('local20','guarded20','geometry_hurdle_gamma')
METRICS=('CRPS','q95_pinball','q99_pinball','mean_abs_error','one_draw_MAE',
         'zero_Brier','central90_coverage','above95','above99','interval_width')
COLS=('ion_births','collected_e','ion_up','ion_down','ion_top','ion_bottom')

def require(ok,msg):
    if not ok:raise ValueError(msg)
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for s in iter(lambda:f.read(1024*1024),b''):h.update(s)
    return h.hexdigest()
def rows(p):
    with p.open(newline='') as f:return list(csv.DictReader(f))
def write_json(p,x):
    p.write_text(json.dumps(x,indent=2,sort_keys=True,allow_nan=False)+'\n')

def load(old_root,extra_root):
    m=json.loads((old_root/'manifest.json').read_text())
    require(m.get('status')=='complete' and m.get('n_completed_avalanches')==84000,
            'Original A20 manuscript / number does not match')
    extra_manifest=extra_root/'manifest.json'
    em=json.loads(extra_manifest.read_text()) if extra_manifest.is_file() else None
    if em is not None and em.get('new_full_Garfield_avalanches')!=16000:
        # Folder with both original A20 and extra A100 CSV for local development.
        em=None
    train={};test={}
    for g in GEOMS:
        a20=old_root/f'trials_{g}.csv'
        require(sha(a20)==m['files_sha256'][a20.name], 'Original A20 SHA: '+g)
        cells=defaultdict(list);ids=defaultdict(set)
        for r in rows(a20):
            require(r['configuration']==g and r['original_geometry']==g,'A20 geometry')
            b=int(r['box_id']);s=int(r['sample']);ids[b].add(s)
            cells[b].append([int(r[c]) for c in COLS])
        require(len(cells)==900 and all(len(x)==20 and len(ids[b])==20
                 for b,x in cells.items()), 'A20 900 x20 incomplete: '+g)
        train[g]={b:np.array(a,dtype=float) for b,a in cells.items()}
        p=extra_root/f'trials_extra_{g}.csv'
        if em is not None:
            require(sha(p)==em['files_sha256'][p.name], 'A100 SHA: '+g)
        else:
            # Archived local export can merge A20 and A100 in one folder.
            require(sha(p)==EXTRA_SHA[p.name], 'Unknown test A100 SHA '+g)
        cells=defaultdict(list);samples=defaultdict(set)
        for r in rows(p):
            require(r['configuration']==g and r['original_geometry']==g,'A100 geometry')
            b=int(r['box_id']);s=int(r['sample'])
            require(s not in ids[b] and s not in samples[b], 'Old/new identity overlap')
            samples[b].add(s);cells[b].append([int(r[c]) for c in COLS])
        require(len(cells)==50 and all(len(x)==80 and len(samples[b])==80
                for b,x in cells.items()), 'A100 new80 selected cells incomplete '+g)
        test[g]={b:np.asarray(a,dtype=float) for b,a in cells.items()}
    target=extra_root/'target_cells.csv'
    if em is not None:require(sha(target)==em['files_sha256'][target.name], 'A100 target sha')
    else:
        require(sha(target)==TARGET_SHA, 'A100 target fixture sha')
    targets={g:set() for g in GEOMS}
    for r in rows(target):
        if r['geometry'] in GEOMS:targets[r['geometry']].add(int(r['box_id']))
    for g in GEOMS:require(targets[g]==set(test[g]), 'A100 target cell IDs mismatch '+g)
    return train,test,{'old_sha':sha(old_root/'manifest.json'),
         'extra_sha':sha(extra_manifest) if em is not None else 'fixture_by_source_sha',
         'targets_sha':sha(target)}

def neighbor_ids(g,b):
    # Exact frozen guarded A20 center-wire heuristic, no v=0 crossing.
    i,j=divmod(b,30)
    u=-1.5+(i+.5)*.1;v=-3+(j+.5)*.2
    wire=lambda a,c:min([math.hypot(a,c)]+([
          math.hypot(a-1.5,c), math.hypot(a+1.5,c)] if g.endswith('FW75') else []))
    if wire(u,v)<.35:return []
    donors=[]
    for di in range(-2,3):
        for dj in range(-2,3):
            if di==dj==0:continue
            ii,jj=i+di,j+dj
            if not(0<=ii<30 and 0<=jj<30):continue
            un=-1.5+(ii+.5)*.1;vn=-3+(jj+.5)*.2
            if (vn>0)!=(v>0):continue
            if (abs(v)<=.6 or abs(vn)<=.6) and (un>0)!=(u>0):continue
            if wire(un,vn)<.35:continue
            donors.append((ii*30+jj,math.exp(-.5*(di*di+dj*dj))))
    z=sum(w for _,w in donors)
    return [(b,w/z) for b,w in donors] if z else []

def mixture(g,b,train):
    local=train[g][b]
    other=neighbor_ids(g,b)
    if not other:return local,np.full(len(local),1./len(local)),False
    data=[local];weights=[np.full(len(local),.25/len(local))]
    for ix,w in other:
        v=train[g][ix];data.append(v);weights.append(np.full(len(v),.75*w/len(v)))
    return np.concatenate(data,axis=0),np.concatenate(weights),True

def quantile_emp(x,w,p):
    idx=np.argsort(x,kind='stable');y=x[idx];v=w[idx]
    return float(y[np.searchsorted(np.cumsum(v),p,side='left').clip(0,len(y)-1)])

def weighted_crps(x,w,truth):
    """E|X-y|-.5 E|X-X'|, sorted O(n log n), exact weighted empirical."""
    order=np.argsort(x);x=x[order];w=w[order]
    cw=np.cumsum(w);cx=np.cumsum(w*x)
    pair=np.sum(w*(x*np.r_[0,cw[:-1]]-np.r_[0,cx[:-1]]))
    return np.sum(w[:,None]*np.abs(x[:,None]-truth[None,:]),axis=0)-pair

def fit_gamma(data):
    # All 18k originals, geometry-specific. No target new80 truth used.
    x=np.concatenate([y[:,0] for y in data.values()]);pzero=(np.sum(x==0)+.5)/(len(x)+1.)
    positive=x[x>0]
    k,loc,theta=gamma.fit(positive,floc=0)
    require(loc==0 and k>0 and theta>0 and math.isfinite(k),'Gamma fit failed')
    return dict(pzero=float(pzero),shape=float(k),scale=float(theta),
                n_train=int(len(x)),n_train_zero=int(np.sum(x==0)))

def gamma_crps(y,p):
    """Analytic CRPS of point mass at zero + Gamma(k,theta), y>=0."""
    z=p['pzero'];k=p['shape'];s=p['scale'];mu=k*s
    F=gamma.cdf(y,k,scale=s)
    Fplus=gamma.cdf(y,k+1,scale=s)
    exabs=y*(2*F-1)+mu*(1-2*Fplus)
    einter=2*mu*math.exp(gammaln(k+.5)-gammaln(k+1)) / math.sqrt(math.pi)
    return z*y+(1-z)*exabs-z*(1-z)*mu-.5*(1-z)**2*einter

def get_gq(p,q):
    return 0. if q<=p['pzero'] else float(gamma.ppf((q-p['pzero'])/(1-p['pzero']),
                                                   p['shape'],scale=p['scale']))

def scored_emp(truth,donors,w):
    x=donors[:,0];y=truth[:,0];zero=float(np.sum(w*(x==0)))
    q05=quantile_emp(x,w,.05);q95=quantile_emp(x,w,.95);q99=quantile_emp(x,w,.99)
    mean=float(np.sum(w*x));crps=weighted_crps(x,w,y)
    return score(y,mean,zero,crps,q05,q95,q99)

def scored_gamma(truth,p):
    y=truth[:,0];q05,q95,q99=(get_gq(p,z) for z in (.05,.95,.99))
    return score(y,(1-p['pzero'])*p['shape']*p['scale'],p['pzero'],
                 gamma_crps(y,p),q05,q95,q99)

def score(y,mean,pzero,crps,q05,q95,q99):
    pin=lambda a,z:float(np.mean(np.maximum(a*(y-z),(a-1)*(y-z))))
    return dict(CRPS=float(np.mean(crps)),q95_pinball=pin(.95,q95),
                q99_pinball=pin(.99,q99),mean_abs_error=float(abs(np.mean(y)-mean)),
                one_draw_MAE=float(np.mean(abs(y-mean))),
                zero_Brier=float(np.mean(((y==0).astype(float)-pzero)**2)),
                central90_coverage=float(np.mean((y>=q05)&(y<=q95))),
                above95=float(np.mean(y>q95)),above99=float(np.mean(y>q99)),
                interval_width=float(q95-q05))

def angular_summary(truth,donors,w):
    # Per ION sector fraction, physical exact ion counts (counts sum to births).
    res={};totaltruth=float(truth[:,0].sum());totalpred=float(np.sum(w*donors[:,0]))
    for name,j in (('up',2),('down',3),('top',4),('bottom',5)):
        if totaltruth<=0 or totalpred<=0:res[name]=None
        else:res[name]=dict(truth=float(truth[:,j].sum()/totaltruth),
             predicted=float(np.sum(w*donors[:,j])/totalpred),
             absolute_error=float(abs(truth[:,j].sum()/totaltruth-
                                     np.sum(w*donors[:,j])/totalpred)))
    return res

def cell_class(g,b):
    i,j=divmod(b,30);u=-1.5+(i+.5)*.1;v=-3+(j+.5)*.2
    dist=min([math.hypot(u,v)]+([math.hypot(u-1.5,v),math.hypot(u+1.5,v)]
               if g.endswith('FW75') else []))
    return 'near_wire' if dist<.35 else ('near_midplane' if abs(v)<.6 else 'bulk')

def ci(x,reps,rng):
    a=np.array(x,float);require(np.all(np.isfinite(a)),'Nonfinite score')
    vals=a[rng.integers(len(a),size=(reps,len(a)))].mean(axis=1)
    return dict(mean=float(a.mean()),CI95=[float(k) for k in np.quantile(vals,[.025,.975])],n=len(a))

def analyze(train,test,source,reps,seed,out):
    require(not out.exists(),'Immutable output exists: '+str(out))
    rng=np.random.default_rng(seed)
    records=[];fits={}
    for g in GEOMS:
        p=fit_gamma(train[g]);fits[g]=p
        for b,truth in sorted(test[g].items()):
            local=train[g][b];wlocal=np.full(20,1/20)
            pooled,wpool,pooledflag=mixture(g,b,train)
            d={ 'local20':scored_emp(truth,local,wlocal),
                'guarded20':scored_emp(truth,pooled,wpool),
                'geometry_hurdle_gamma':scored_gamma(truth,p)}
            sector={ 'local20':angular_summary(truth,local,wlocal),
                     'guarded20':angular_summary(truth,pooled,wpool)}
            records.append(dict(geometry=g,box_id=b,region=cell_class(g,b),
               n_new80=80,has_eligible_neighbors=pooledflag,
               n_guarded_donors=len(pooled),metrics=d,sector=sector))
    summary={};decision={}
    for g in GEOMS:
        a=[r for r in records if r['geometry']==g]
        require(len(a)==50,'Missing 50 original selected test cells '+g)
        summary[g]={'all_cells':{},'by_region':{},'paired_delta':{},'sector':{}}
        for mode in MODES:
            summary[g]['all_cells'][mode]={m:ci([r['metrics'][mode][m] for r in a],reps,rng)
               for m in METRICS}
        for mode in ('guarded20','geometry_hurdle_gamma'):
            summary[g]['paired_delta'][mode+'_minus_local20']={m:ci(
               [r['metrics'][mode][m]-r['metrics']['local20'][m] for r in a],reps,rng)
               for m in ('CRPS','mean_abs_error','q95_pinball','q99_pinball',
                         'above95','above99','zero_Brier')}
        for region in ('near_wire','near_midplane','bulk'):
            rr=[r for r in a if r['region']==region]
            summary[g]['by_region'][region]={'n_cells':len(rr)}
            if rr:
                for mode in MODES:
                    summary[g]['by_region'][region][mode]={m:ci(
                        [r['metrics'][mode][m] for r in rr],min(reps,3000),rng)
                        for m in ('CRPS','mean_abs_error','above99')}
        for mode in ('local20','guarded20'):
            summary[g]['sector'][mode]={name:ci([r['sector'][mode][name]['absolute_error']
              for r in a if r['sector'][mode][name] is not None],reps,rng)
              for name in ('up','down','top','bottom')}
        baseline=summary[g]['all_cells']['local20']['CRPS']['mean']
        guard=summary[g]['all_cells']['guarded20']['CRPS']['mean']
        gam=summary[g]['all_cells']['geometry_hurdle_gamma']['CRPS']['mean']
        decision[g]={'best_count_CRPS_on_selected_cells':min(
            MODES,key=lambda name:summary[g]['all_cells'][name]['CRPS']['mean']),
            'guarded_CRPS_relative_change_vs_local20_pct':100*(guard/baseline-1),
            'geometry_gamma_CRPS_relative_change_vs_local20_pct':100*(gam/baseline-1),
            'joint_cloud_eligible':'only local20/guarded20 preserve exact shared ion-sector histories; all A20 samples have at most 5 space/time points',
            'holdout_status':'selected-cell old20-new80 development; NOT final blind muon'}
    payload=dict(schema='NA60_P3_FAIR_20TRAIN_80TEST_V1',status='complete',
      model_per_geometry=GEOMS,train='18000 old A20 independent physical trials per geometry',
      validation='4000 different new A100 microscopic outcomes in 50 selected cells per geometry',
      warning='Selection of test cells source-conditioned by old A20 and old muons; not representative final holdout',
      warning2='Full ion birth cloud, drift times and future cathode response NOT tested',
      details='Gamma count-only challenger trained on SAME A20 physical source; local/guarded are joint sector empirical responses',
      source=source,method={'guarded':'original guarded25local/75neighbors; radius 2; Gaussian center-cell weights; no cross v=0/near-wire','gamma':'global pzero + MLE positive Gamma per geometry (not an RKF hybrid)'},
      gamma_parameters=fits,model_rank_dev=decision,by_geometry=summary)
    out.mkdir(parents=True)
    write_json(out/'summary.json',payload)
    write_json(out/'cell_metrics.json',records)
    seal={'schema':'NA60_P3_FAIR_SHA_V1','sha256':{
        p.name:sha(p) for p in (out/'summary.json',out/'cell_metrics.json')},
        'sources':source,'seed':seed,'bootstrap_reps':reps}
    write_json(out/'manifest.json',seal)
    for g in GEOMS:
        print('\n'+g+' : FAIR OLD20 TRAIN -> NEW80 SELECTED-CELL TEST')
        for model in MODES:
            s=summary[g]['all_cells'][model]
            print('%-22s CRPS %9.3f   MAE(cell mean) %9.3f   zero Brier %.5f   above95 %.3f%%   above99 %.3f%%'%(model,s['CRPS']['mean'],s['mean_abs_error']['mean'],s['zero_Brier']['mean'],100*s['above95']['mean'],100*s['above99']['mean']))
        for cand in ('guarded20','geometry_hurdle_gamma'):
            d=summary[g]['paired_delta'][cand+'_minus_local20']['CRPS']
            print('  ΔCRPS %s - local20: %+.3f 95%% cell-CI [%.3f,%.3f]'%(cand,d['mean'],*d['CI95']))
        print('  near_wire selected cells:',summary[g]['by_region']['near_wire']['n_cells'])
    print('PASS SHA-verified p3-only fair selected-cell test; no Garfield run:',out)

def selftest():
    x=np.array([0.,1.,2.,10.]);w=np.full(4,.25);y=np.array([0.,3.,8.])
    crps=weighted_crps(x,w,y)
    brute=np.array([np.sum(w*abs(x-z))-.5*np.sum(np.outer(w,w)*abs(x[:,None]-x[None,:])) for z in y])
    require(np.allclose(crps,brute),'Empirical CRPS identity')
    q={'pzero':.17,'shape':1.7,'scale':200.}
    # numerical integration: CRPS(F,y) = integral (F(z)-1_{z>=y})^2 dz
    from scipy.integrate import quad
    for obs in (0.,100.,450.):
        F=lambda z:q['pzero']+(1-q['pzero'])*gamma.cdf(z,q['shape'],scale=q['scale'])
        numeric=quad(lambda z:F(z)**2,0,obs,epsabs=1.e-7)[0]+quad(lambda z:(1-F(z))**2,obs,4000,epsabs=1.e-7)[0]
        require(abs(gamma_crps(np.array([obs]),q)[0]-numeric)<.05,'Gamma+zero CRPS identity')
    for g in GEOMS:
        for b in (0,200,449,899):
            nei=neighbor_ids(g,b)
            require(all(0<=j<900 and w>0 for j,w in nei),'Invalid physical neighbor')
            if nei:require(abs(sum(w for _,w in nei)-1)<1.e-12,'Neighbor weights')
    print('PASS synthetic weighted CRPS / hurdle-Gamma analytic CRPS / physical guards; NO physics')

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('mode',choices=('selftest','audit'))
    p.add_argument('--old-dir',type=Path,default=OLD)
    p.add_argument('--extra-dir',type=Path,default=EXTRA)
    p.add_argument('--out',type=Path,default=OUTPUT)
    p.add_argument('--bootstrap',type=int,default=4000)
    p.add_argument('--seed',type=int,default=20261104)
    a=p.parse_args()
    if a.mode=='selftest':selftest()
    else:
        require(200<=a.bootstrap<=10000,'Use 200..10000 cell-bootstrap replicas')
        d,t,s=load(a.old_dir,a.extra_dir)
        analyze(d,t,s,a.bootstrap,a.seed,a.out)
if __name__=='__main__':main()