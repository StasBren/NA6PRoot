#!/usr/bin/env python3
"""Global finite-training bootstrap of guarded neighbor-pooled A20/A100 for
120 actual Geant4+TrackHeed parent muons and 4 B=0 MWPC geometries.
READ ONLY: no Garfield/Heed/Geant4/Magboltz simulation or original file edits.
"""
import argparse,csv,json,math,time
from pathlib import Path
import numpy as np
from pi_spatial_lut_muon120_crn_v1 import (GEOMS,sha,load_archive,load_native,make_records,compile_neighbors)
MODES=('A20','A100','guarded_A20','guarded_A100')

def prepare(tables,rec):
    neigh=compile_neighbors(tables,True);data={};parents=None
    for g in GEOMS:
        rr=[r for r in rec if r['geom']==g]
        if len(rr)!=120:raise ValueError('120 parent muons required')
        ids=[];lens=[]
        for r in rr:
            cells=[b for b,n in r['counts'] for _ in range(n)]
            if len(cells)!=r['n_primary']:raise ValueError('Heed count mismatch')
            ids.extend(cells);lens.append(len(cells))
        ids=np.array(ids,np.intp);old,extra=tables[g];ncell=len(old)
        if len(ids)!=5881:raise ValueError('Wrong original Heed population')
        tgt=np.zeros(ncell,bool);more=np.zeros((ncell,80,5),np.int64)
        for b,v in extra.items():tgt[b]=True;more[b]=v
        donor=np.zeros((ncell,24),np.intp);cdf=np.ones((ncell,24),float);has=np.zeros(ncell,bool)
        for b,(c,w) in enumerate(neigh[g]):
            if len(c)>24:raise ValueError('Neighbor kernel invalid')
            if len(c):donor[b,:len(c)]=c;cdf[b,:len(c)]=w;has[b]=True
        pp=[r['parent'] for r in rr]
        if parents is None:parents=pp
        elif parents!=pp:raise ValueError('Unpaired parent ordering')
        data[g]={'r':rr,'old':old,'extra':more,'tgt':tgt,'ids':ids,
            'offs':np.cumsum([0]+lens[:-1]),'donor':donor,'cdf':cdf,'has':has}
    return data

def simulate(data,reps,seed):
    rng=np.random.default_rng(seed)
    out=np.zeros((reps,120,4,4,5),np.int64)
    for t in range(reps):
        for gi,g in enumerate(GEOMS):
            q=data[g];old=q['old'];extra=q['extra'];b=q['ids'];ncell=len(old);n=len(b)
            # One single bootstrapped cell library per global iteration, shared
            # among all physical parents and primary electrons in that cell.
            boot20=rng.integers(20,size=(ncell,20),dtype=np.intp)
            boot80=rng.integers(80,size=(ncell,80),dtype=np.intp)
            orig=old[b,boot20[b,rng.integers(20,size=n,dtype=np.intp)]]
            extended=orig.copy()
            choose=q['tgt'][b] & (rng.random(n)>=.2)
            if np.any(choose):
                ix=np.flatnonzero(choose);cell=b[ix]
                extended[ix]=extra[cell,boot80[cell,rng.integers(80,size=len(ix),dtype=np.intp)]]
            new_neigh=q['has'][b] & (rng.random(n)>=.25)
            po=orig.copy();pe=extended.copy()
            if np.any(new_neigh):
                ix=np.flatnonzero(new_neigh);cell=b[ix];u=rng.random(len(ix))
                ranks=np.sum(u[:,None]>q['cdf'][cell],axis=1).clip(0,23)
                dest=q['donor'][cell,ranks]
                vals=old[dest,boot20[dest,rng.integers(20,size=len(ix),dtype=np.intp)]]
                po[ix]=vals;pe[ix]=vals
            for mi,v in enumerate((orig,extended,po,pe)):
                out[t,:,gi,mi,:]=np.add.reduceat(v,q['offs'],axis=0)
            if t==0:
                for pi,r in enumerate(q['r']):
                    if not r['n_target'] and not np.array_equal(out[t,pi,gi,0],out[t,pi,gi,1]):
                        raise ValueError('CRN changed unselected A20/A100 physics')
        if (t+1)%max(50,reps//10)==0:
            print(f'PASS global shared-cell bootstrap {t+1}/{reps} × 120 parents × 4 geometries',flush=True)
    return out

def metrics(data,x):
    refs=np.array([[data[g]['r'][i]['ref_gain'] for g in GEOMS] for i in range(120)])
    pred=x[:,:,:,:,0].astype(float)
    lo,hi=np.quantile(pred,[.05,.95],axis=0);avg=np.mean(pred,axis=0)
    if pred.shape[1:]!=(120,4,4):raise ValueError('Unexpected prediction dimensions')
    covered=(lo<=refs[:,:,None])&(refs[:,:,None]<=hi)
    out=[]
    for i in range(120):
        for gi,g in enumerate(GEOMS):
            r=data[g]['r'][i];row={'parent':r['parent'],'geometry':g,'n_primaries':r['n_primary'],
              'n_target':r['n_target'],'full_Garfield_ions':r['ref_gain']}
            for mi,m in enumerate(MODES):
                row[m+'_covered']=int(covered[i,gi,mi]);row[m+'_mean']=float(avg[i,gi,mi]);
                row[m+'_p05']=float(lo[i,gi,mi]);row[m+'_p95']=float(hi[i,gi,mi]);
                row[m+'_MAE']=float(abs(avg[i,gi,mi]-refs[i,gi]));row[m+'_width']=float(hi[i,gi,mi]-lo[i,gi,mi])
                for label,k,j in [('hemi',1,2),('topbottom',3,4)]:
                    a=x[:,i,gi,mi,k];b=x[:,i,gi,mi,j];mask=(a+b)>0
                    truth=data[g]['r'][i];ra=truth['ref_up'] if k==1 else truth['ref_top'];rb=truth['ref_down'] if k==1 else truth['ref_bottom']
                    if ra+rb>0 and np.count_nonzero(mask)>=100:
                        q=np.quantile((a[mask]-b[mask])/(a[mask]+b[mask]),[.05,.95]);v=(ra-rb)/(ra+rb)
                        row[m+'_'+label+'_covered']=int(q[0]<=v<=q[1])
                    else:row[m+'_'+label+'_covered']=''
            out.append(row)
    return out

def report(rows,reps,source):
    if len(rows)!=480 or not all(rows[4*i+j]['geometry']==g for i in range(120) for j,g in enumerate(GEOMS)):
        raise ValueError('Parent and geometry alignment broken')
    pergeo={};overall={}
    for geo in GEOMS:
        rr=[r for r in rows if r['geometry']==geo]
        pergeo[geo]={m:one(rr,m) for m in MODES}
    for m in MODES:overall[m]=one(rows,m)
    z=np.array([[r[m+'_covered'] for m in MODES] for r in rows],float).reshape(120,4,4)
    rng=np.random.default_rng(20261101);ix=rng.integers(120,size=(12000,120))
    deltas={}
    for name,k,j in [('guarded_A20_minus_A20',2,0),('guarded_A100_minus_A100',3,1)]:
        d=(z[:,:,k]-z[:,:,j]).mean(axis=1)
        deltas[name]={'delta_pp':100*float(np.mean(d)),'parent_cluster_bootstrap95_pp':np.quantile(100*np.mean(d[ix],axis=1),[.025,.975]).tolist()}
    return {'schema':'NA60_MWPC_GUARDED_LUT_GLOBAL_CELL_TRAINING_BOOTSTRAP_MU120_V1',
      'status':'complete','n_parent_muons':120,'n_correlated_geometry_tests':480,
      'global_cell_training_bootstrap_reps':reps,'overall':overall,'by_geometry':pergeo,
      'parent_clustered_coverage_deltas':deltas,'pooling_parameters':{'radius_bins':2,'neighbor_probability':.75,'local_probability':.25,'guard':'same v hemisphere, protected wire centers and near-axis sign'},
      'critical_interpretation':'Source-conditioned diagnostic: parent TrackHeed positions helped choose A100 cells and guided pooling hyperparameter selection. Shared cell bootstrap approximates finite Garfield training uncertainty but not an exact posterior; physical 120 full Garfield true outcomes have one random history each. No independent external validation and no induced cathode strips.',
      'source':source}

def one(rr,m):
    def sector(name):
        y=[r[m+'_'+name+'_covered'] for r in rr if r[m+'_'+name+'_covered']!='']
        return {'n_defined':len(y),'coverage':float(np.mean(y)) if y else None}
    return {'n':len(rr),'covered':sum(r[m+'_covered'] for r in rr),'coverage':float(np.mean([r[m+'_covered'] for r in rr])),
      'mean_abs_error':float(np.mean([r[m+'_MAE'] for r in rr])),
      'mean_90PI_width':float(np.mean([r[m+'_width'] for r in rr])),
      'hemisphere':sector('hemi'),'top_bottom':sector('topbottom')}

def selftest():
    rng=np.random.default_rng(101)
    a=np.zeros((24,20,5),np.int64);a[:,:,0]=10;a[:,:,1]=4;a[:,:,2]=6
    ix=rng.integers(20,size=(24,20));points=np.array([4,5,4,11]);draw=a[points,ix[points,np.array([1,5,7,13])]]
    assert draw.shape==(4,5) and np.all(draw[:,0]==10)
    assert all(np.sum(draw[:,1:3],axis=1)==draw[:,0])
    # Full synthetic 120×4 hierarchy regression: all parents reuse one cell,
    # no target and no pooling; old/new/guarded results exactly identical.
    fake={}
    for g in GEOMS:
        base=np.zeros((24,20,5),np.int64)
        base[:,10:,0]=100; base[:,10:,1]=40; base[:,10:,2]=60
        fake[g]={ 'r':[{'parent':str(i),'n_target':0} for i in range(120)],
          'old':base,'extra':np.zeros((24,80,5),np.int64),
          'tgt':np.zeros(24,bool),'ids':np.zeros(120,np.intp),
          'offs':np.arange(120,dtype=np.intp),
          'has':np.zeros(24,bool),'donor':np.zeros((24,24),np.intp),
          'cdf':np.ones((24,24),float)}
    sampled=simulate(fake,100,17)
    assert np.array_equal(sampled[:,:,:,0,:],sampled[:,:,:,1,:])
    assert np.array_equal(sampled[:,:,:,0,:],sampled[:,:,:,2,:])
    assert np.array_equal(sampled[:,:,:,0,:],sampled[:,:,:,3,:])
    assert np.all(sampled[...,0]==sampled[...,1]+sampled[...,2])
    print('PASS full synthetic 120-muon global shared-cell bootstrap CRN/sector identity')
    print('PASS offline shared training-pool indexing and correlated ion-sector selftest; no physics',flush=True)

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('mode',choices=('selftest','run'))
    p.add_argument('--published-dir',type=Path);p.add_argument('--native',action='store_true');p.add_argument('--output',type=Path)
    p.add_argument('--reps',type=int,default=1000);p.add_argument('--seed',type=int,default=20261031)
    args=p.parse_args()
    if args.mode=='selftest':selftest();return
    if args.native==bool(args.published_dir):p.error('Specify exactly one of --native and --published-dir')
    if args.native:
        tables,targets,mu,prim,cross,old,src,root=load_native()
        output=args.output or root/'diagnose_guarded_spatial_pool_global_shared_training_bootstrap_v1'
    else:
        tables,targets,mu,prim,cross,old,src=load_archive(args.published_dir)
        output=args.output
        if output is None:p.error('Published archive mode requires a new --output directory')
    if output.exists():raise FileExistsError('Refusing to overwrite '+str(output))
    rec=make_records(tables,targets,mu,prim,cross,old)
    data=prepare(tables,rec);start=time.monotonic()
    pred=simulate(data,args.reps,args.seed)
    rr=metrics(data,pred);summ=report(rr,args.reps,src)
    summ['runtime_seconds']=round(time.monotonic()-start,3)
    output.mkdir(parents=True)
    with (output/'per_muon_4models.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rr[0]));writer.writeheader();writer.writerows(rr)
    (output/'summary.json').write_text(json.dumps(summ,sort_keys=True,indent=2)+'\n')
    files=[output/'per_muon_4models.csv',output/'summary.json']
    (output/'manifest.json').write_text(json.dumps({'status':'complete','files_sha256':{q.name:sha(q) for q in files}},indent=2)+'\n')
    print(json.dumps({'overall':summ['overall'],'parent_clustered_coverage_deltas':summ['parent_clustered_coverage_deltas'],'runtime_seconds':summ['runtime_seconds']},indent=2),flush=True)
    print('PASS 120 parents × 4 geoms, shared global finite-cell-training bootstrap; NO Garfield physics rerun',output,flush=True)

if __name__=='__main__':main()