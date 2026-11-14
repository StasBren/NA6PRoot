#!/usr/bin/env python3
"""Resume-safe 1000-real-Geant4-muon Penning057 MWPC collecting-wire truth run.
One FULL AvalancheMicroscopic call per fresh TrackHeed electron, using an
isolated clone of the 30k native engine. This script NEVER runs Garfield here;
run it under WSL where Garfield++ is installed.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,math,os,subprocess,sys,time,zipfile
from collections import defaultdict,Counter
from pathlib import Path

SCRIPT=Path(__file__).resolve()
REPO=SCRIPT.parents[2]
SELECTION=REPO/'test/mwpc-stage-c/p3_wiretruth_1000_selection.csv'
SEL_META=REPO/'test/mwpc-stage-c/p3_wiretruth_1000_selection.json'
OUT=REPO/'test_runs/mwpc_stage_c/p3_noFW_penning057_wiretruth_1000_v1'
BINARY=REPO/'test/mwpc-avalanche/build/mwpc_p3_noFW_penning057_wiretruth_v1'
FILES=['crossing','seed','ion_birth_sample','heed_primaries','heed_clusters',
       'timing','muon_timing','wire_seed','wire_muon']

def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for chunk in iter(lambda:f.read(4*1024*1024),b''):
            h.update(chunk)
    return h.hexdigest()
def readcsv(path):
    with open(path,newline='') as f:return list(csv.DictReader(f))
def writejson(path,obj):
    path.write_text(json.dumps(obj,indent=2,sort_keys=True)+'\n')
def selection():
    meta=json.loads(SEL_META.read_text());h=sha(SELECTION)
    if h!=meta['selected_csv_sha256']:raise RuntimeError('Frozen 1000 input SHA mismatch')
    rows=readcsv(SELECTION)
    if len(rows)!=1000 or len({tuple(r[k] for k in ('source','chunk','event','hit_index','track_id','station','det_id')) for r in rows})!=1000:
        raise RuntimeError('Frozen selection not 1000 distinct Geant4 crossings')
    if Counter(r['source'].lower().split('40')[0].capitalize() for r in rows)!= {'Jpsi':334,'Phi':333,'Omega':333}:
        raise RuntimeError('Frozen channel allocation mismatch')
    return rows,meta

def native_id(r):
    return '|'.join(r[k] for k in ('source','chunk','event','hit_index','track_id','station','det_id'))

def chunk_input(i):return OUT/'inputs'/f'chunk_{i:04d}.csv'
def chunk_output(i):return OUT/'chunks'/f'chunk_{i:04d}'

def prepare():
    rows,meta=selection()
    OUT.mkdir(parents=True,exist_ok=True)
    config=OUT/'locked_selection.json'
    lock={'source_sha256':sha(SELECTION),'selected_rows':1000,
          'input_chunks':50,'physics':'p3_noFW Penning Ar 0.57 3mm pitch 3+3mm 1800V B0',
          'native_binary':'mwpc_p3_noFW_penning057_wiretruth_v1',
          'note':'New Heed seeds; not bitwise replay of 30k old TrackHeed primaries.'}
    if config.exists() and json.loads(config.read_text())!=lock:
        raise RuntimeError('Output path already locked to different selection. Refusing overwrite')
    if not config.exists():writejson(config,lock)
    (OUT/'inputs').mkdir(exist_ok=True);(OUT/'chunks').mkdir(exist_ok=True)
    hdr=list(rows[0].keys())
    assert all(r.keys()==rows[0].keys() for r in rows)
    for i in range(50):
        p=chunk_input(i)
        if p.exists():
            existing=readcsv(p)
            if existing!=rows[20*i:20*i+20]:
                raise RuntimeError('Existing input chunk differs: '+str(p))
            continue
        with p.open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=hdr,lineterminator='\n');w.writeheader();w.writerows(rows[20*i:20*i+20])
    print('PREPARED: 50 immutable input chunks x 20 Geant4 muons; 334/333/333 channels; source SHA',sha(SELECTION),flush=True)

def validate_chunk(d,inputpath):
    pref=d/'microscopic'
    paths={n:Path(str(pref)+'_'+n+'.csv') for n in FILES}
    if any(not f.exists() for f in paths.values()):raise RuntimeError('Missing one or more native outputs')
    crossings=readcsv(paths['crossing']);seeds=readcsv(paths['seed'])
    wb=readcsv(paths['wire_seed']);wm=readcsv(paths['wire_muon'])
    if len(crossings)!=20:raise RuntimeError('Not exactly 20 true muon results')
    crossing_h=list(crossings[0]);idcol=crossing_h[0]
    ids={native_id(r) for r in readcsv(inputpath)}
    if set(r[idcol] for r in crossings)!=ids:raise RuntimeError('Wrong muon/crossing IDs')
    seedidx={(r['crossing_id'],int(r['primary_index'])):r for r in seeds}
    if len(seedidx)!=len(seeds):raise RuntimeError('Duplicate Heed primary rows')
    ws=defaultdict(lambda:[0,0])
    for r in wb:
        k=(r['crossing_id'],int(r['primary_index']))
        if k not in seedidx:raise RuntimeError('Wire truth has no original seed')
        ws[k][0]+=int(r['collected_e']);ws[k][1]+=int(r['ion_births'])
    for k,r in seedidx.items():
        q,i=ws[k]
        if q!=int(r['collected_e']) or i!=int(r['ion_births']):
            raise RuntimeError('SEED wire-sum != original Garfield '+str(k))
    wmby=defaultdict(lambda:[0,0]);active=defaultdict(set)
    for r in wm:
        cid=r['crossing_id'];q=int(r['collected_e']);i=int(r['ion_births'])
        wmby[cid][0]+=q;wmby[cid][1]+=i
        if q>0:active[cid].add(int(r['wire_index_global']))
    seedby=defaultdict(lambda:[0,0,0])
    capcount=0;failed=0
    for r in seeds:
        cid=r['crossing_id'];seedby[cid][0]+=int(r['collected_e'])
        seedby[cid][1]+=int(r['ion_births']);seedby[cid][2]+=1
        capcount+=int(r['capped']);failed+=int(r['successful'])!=1
    for r in crossings:
        cid=r[idcol];q=int(r['collected_avalanche_electrons']);nb=int(r['ion_births'])
        if wmby[cid]!=[q,nb] or seedby[cid][:2]!=[q,nb]:
            raise RuntimeError('MUON charge/ion sum != original Garfield '+cid)
        if seedby[cid][2]!=int(r['heed_primaries']):
            raise RuntimeError('MUON actual seeds != Geant4/Heed record '+cid)
    # Critical: reject a censored spectrum instead of silently treating it as reference.
    if capcount:raise RuntimeError('Garfield had capped avalanches: '+str(capcount))
    # The failed AvalancheElectron() flag is recorded; not all flags must be 1.
    return dict(muons=len(crossings),primary_electrons=len(seeds),
                ions=sum(int(r['ion_births']) for r in crossings),
                total_collected_e=sum(int(r['collected_avalanche_electrons']) for r in crossings),
                two_or_more_wires=sum(len(v)>=2 for v in active.values()),
                capped=capcount,unsuccessful=failed,
                files={n:{'sha256':sha(p),'bytes':p.stat().st_size} for n,p in paths.items()})

def check_binary():
    if not BINARY.exists():raise RuntimeError('Missing Garfield binary. Build CMake target first: '+str(BINARY))

def run(max_chunks):
    prepare();check_binary()
    done=0;start=time.monotonic()
    for i in range(50):
        d=chunk_output(i);receipt=d/'complete.json';inputpath=chunk_input(i)
        if receipt.exists():
            r=json.loads(receipt.read_text())
            if r['input_sha256']!=sha(inputpath):raise RuntimeError('Input has changed for completed chunk')
            continue
        if d.exists():
            raise RuntimeError(f'Unsealed partial output found at {d}; inspect physics.log before any manual quarantine')
        if done>=max_chunks:break
        d.mkdir(parents=True)
        pref=d/'microscopic';log=d/'physics.log'
        argv=[str(BINARY),'--wiretruth-p3-noFW',str(inputpath),str(pref)]
        print('RUN native full Garfield++',i,'x20 events,',flush=True)
        ts=time.monotonic()
        with log.open('w') as f:
            f.write(' '.join(argv)+'\n');f.flush()
            p=subprocess.run(argv,stdout=f,stderr=subprocess.STDOUT,cwd=d)
        if p.returncode!=0:
            raise RuntimeError(f'Garfield chunk {i:04d} failed; inspect {log}; no complete.json written')
        stats=validate_chunk(d,inputpath)
        stats.update(input_sha256=sha(inputpath),native_binary_sha256=sha(BINARY),
                     chunk=i,wall_seconds=time.monotonic()-ts,
                     physics_log_sha256=sha(log))
        writejson(receipt,stats)
        print('PASS chunk',i,':',stats['muons'],'muons',stats['primary_electrons'],
              'Heed electrons','wire truth saved','wall_sec',round(stats['wall_seconds'],1),flush=True)
        done+=1
    print('DONE currently sealed',status(quiet=True)['sealed_chunks'],'of 50. New chunks this invocation',done,
          'elapsed seconds',round(time.monotonic()-start,1),flush=True)

def status(quiet=False):
    n=0;prim=0;ions=0;mu=0;unsealed=[]
    for i in range(50):
        d=chunk_output(i);receipt=d/'complete.json'
        if receipt.exists():
            r=json.loads(receipt.read_text());n+=1;prim+=r['primary_electrons'];mu+=r['muons'];ions+=r['ions']
        elif d.exists():unsealed.append(i)
    x={'sealed_chunks':n,'muons':mu,'primary_electrons':prim,'ion_births':ions,'unsealed_chunks':unsealed,'output':str(OUT)}
    if not quiet:print(json.dumps(x,indent=2))
    return x

def audit():
    rows,_=selection();done=status(quiet=True)
    for i in range(50):
        d=chunk_output(i);p=d/'complete.json'
        if not p.exists():continue
        r=json.loads(p.read_text());ip=chunk_input(i)
        if r['input_sha256']!=sha(ip):raise RuntimeError('Input changed: '+str(ip))
        for n,entry in r['files'].items():
            f=d/('microscopic_'+n+'.csv')
            if sha(f)!=entry['sha256']:raise RuntimeError('Tampered native output: '+str(f))
        val=validate_chunk(d,ip)
        for k in ('muons','primary_electrons','ions','total_collected_e','two_or_more_wires'):
            if val[k]!=r[k]:raise RuntimeError('Physics receipt mismatch '+str(i)+' '+k)
    print('AUDIT PASS',done['sealed_chunks'],'sealed chunks /',done['muons'],'actual muons')

def analyse():
    stat=status(quiet=True)
    angles=[0,5,10,15,20,30,60,180]
    x=defaultdict(lambda:{'n':0,'n_multi':0,'sum_q2':0.,'n_active':0})
    bychan=Counter();mu=0
    for i in range(50):
        d=chunk_output(i)
        if not (d/'complete.json').exists():continue
        wire=defaultdict(list)
        for r in readcsv(d/'microscopic_wire_muon.csv'):
            wire[r['crossing_id']].append(int(r['collected_e']))
        cr=readcsv(d/'microscopic_crossing.csv');cid_name=list(cr[0])[0]
        for row in cr:
            cid=row[cid_name];source=cid.split('|',1)[0].lower()
            channel='Jpsi' if source.startswith('jpsi') else ('Phi' if source.startswith('phi') else 'Omega')
            theta=abs(float(row['theta_u_in_deg']))
            bin_=next((f'{angles[j]}-{angles[j+1]}' for j in range(len(angles)-1) if angles[j]<=theta<angles[j+1]),'>=180')
            charges=sorted(wire[cid],reverse=True);total=sum(charges)
            q2=(charges[1]/total if total>0 and len(charges)>1 else 0.)
            # Treat cases of zero collected charge separately: no meaningful q2.
            for group in ('All',channel):
                t=x[(group,bin_)];t['n']+=1;t['n_multi']+=int(q2>0.05)
                t['sum_q2']+=q2;t['n_active']+=int(total>0)
            mu+=1;bychan[channel]+=1
    dest=OUT/'reference_wire_sharing_by_angle.csv'
    with dest.open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['channel','abs_theta_u_bin_deg','muons','active_charge_muons','n_second_wire_gt_5pct','fraction_all','mean_Q2_fraction'])
        for (group,b),a in sorted(x.items()):
            w.writerow([group,b,a['n'],a['n_active'],a['n_multi'],a['n_multi']/a['n'],a['sum_q2']/a['n']])
    md=['# First 1000 full Garfield++ reference — collecting-anode truth','',
       f'- Completed native muons: **{mu}**','- Real Geant4 crossings and freshly simulated TrackHeed/Garfield avalanches.',
       '- `Q2` is the second-largest collected-electron count divided by the total collected count **for the same muon**.',
       '- Zero-charge events remain in the denominator; efficiency requires separate reporting.','',
       '## Current counts','',f'Jpsi {bychan["Jpsi"]}, Phi {bychan["Phi"]}, Omega {bychan["Omega"]}.',
       '',f'Detailed table: `{dest.name}`.',
       '', '**Not yet a Fast4/Fast5 model-vs-truth comparison**: replay both fast models using the saved per-primary Heed inputs for these exact 1000 events.']
    (OUT/'reference_wire_summary.md').write_text('\n'.join(md)+'\n')
    print('WROTE',dest,OUT/'reference_wire_summary.md')

def export_compact():
    """Upload-friendly archive: wire truth + Heed primaries, but NOT the
    potentially huge 128-point ion birth reservoir CSV; it stays on WSL.
    """
    n=status(quiet=True)['muons']
    if n==0:raise RuntimeError('No completed full-Garfield chunks to export')
    out=OUT/f'wiretruth_compact_{n}_muons.zip'
    if out.exists():
        print('Compact archive already present:',out);return
    tmp=out.with_name(out.name+'.part')
    selected=['wire_muon','wire_seed','heed_primaries','heed_clusters',
              'crossing','seed','timing','muon_timing']
    with zipfile.ZipFile(tmp,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        z.write(SELECTION,'selection_1000.csv')
        z.write(SEL_META,'selection_1000.json')
        for f in ('locked_selection.json','reference_wire_sharing_by_angle.csv','reference_wire_summary.md'):
            p=OUT/f
            if p.exists():z.write(p,f)
        for i in range(50):
            d=chunk_output(i);receipt=d/'complete.json'
            if not receipt.exists():continue
            base=f'chunks/chunk_{i:04d}/'
            z.write(chunk_input(i),base+'input.csv')
            z.write(receipt,base+'complete.json')
            for name in selected:
                path=d/('microscopic_'+name+'.csv')
                z.write(path,base+name+'.csv')
    tmp.rename(out)
    print('COMPACT WIRE TRUTH ARCHIVE:',out)

def selftest():
    rows,_=selection();assert len(rows)==1000
    assert sum(1 for r in rows if abs(float(r['theta_u_in_deg']))>=10)>30
    assert set(int(r['station']) for r in rows)==set(range(6))
    print('SELFTEST PASS: 1000 frozen genuine crossings, all stations and elevated-angle events included')

def main():
    p=argparse.ArgumentParser();p.add_argument('command',choices=['selftest','prepare','run','status','audit','analyse','export']);p.add_argument('--max-chunks',type=int,default=50)
    a=p.parse_args()
    if a.command=='selftest':selftest()
    elif a.command=='prepare':prepare()
    elif a.command=='run':run(a.max_chunks)
    elif a.command=='status':status()
    elif a.command=='audit':audit()
    elif a.command=='analyse':analyse()
    elif a.command=='export':export_compact()
if __name__=='__main__':main()
