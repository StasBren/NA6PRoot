#!/usr/bin/env python3
"""Reproducible 100k-crossing/channel/station LUT500-vs-midline comparison.

Use inside NA6PRoot feature/mwpc-stage-c checkout with frozen native exports.
This fast variant uses the project's *Heed-calibrated marked-Poisson* ionization
sampler along REAL Geant4 chords, not a new true TrackHeed replay.

For matched absolute anode charge, provide the J/psi gain-donor NPZ. Both
assignments share the same per-primary avalanche gains, and LUT500 is sampled
only conditional on collection to avoid double counting zero-gain avalanches.
"""
from __future__ import annotations
import argparse
import csv
import gzip
import hashlib
import heapq
import json
import os
import sys
import time
import zipfile
from collections import defaultdict
from pathlib import Path

import numpy as np
from wire_core import (stable_hash, rng_for, build_lut_kernel,
                       observables, compare_one, STATION_WIDTHS, THRESHOLDS)

HERE=Path(__file__).resolve().parent
CHANNELS=("Jpsi","Phi","Omega")
STATIONS=(0,3,5)
NATIVE_REQUIRED={'event','hit_index','track_id','station','det_id','pdg',
   'path_length_mm','u_in_mm','v_in_mm','w_in_mm',
   'u_out_mm','v_out_mm','w_out_mm',
   'theta_u_in_deg','theta_w_in_deg',
   'pu_in_GeV','pv_in_GeV','pw_in_GeV'}
FIELDS=('source','chunk','event','hit_index','track_id','station','det_id','pdg',
        'path_length_mm','u_in_mm','v_in_mm','w_in_mm',
        'u_out_mm','v_out_mm','w_out_mm',
        'theta_u_in_deg','theta_w_in_deg',
        'pu_in_GeV','pv_in_GeV','pw_in_GeV')


def sha(path):
    digest=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(4*1024*1024), b''):
            digest.update(block)
    return digest.hexdigest()


def native_interface(repo):
    stage=Path(repo)/'test/mwpc-stage-c'
    if not stage.is_dir(): raise FileNotFoundError(stage)
    sys.path.insert(0,str(stage))
    from pi_2m_native_six_geometry_fold import crossing_provenance, SOURCE
    from pi_six_geometry_native_timing import birth_sample,HEED_TRACKS,HEED_CLUSTERS
    from run_pi_per_chamber_response import read_heed
    return crossing_provenance,SOURCE,birth_sample,HEED_TRACKS,HEED_CLUSTERS,read_heed


def validate(row, station):
    """Preflight physical muon gas chord; do not mix with non-muons."""
    if int(row['station']) != station or abs(int(row['pdg']))!=13: return False
    vals=[float(row[x]) for x in ('path_length_mm','u_in_mm','v_in_mm',
          'w_in_mm','u_out_mm','v_out_mm','w_out_mm',
          'theta_u_in_deg','theta_w_in_deg',
          'pu_in_GeV','pv_in_GeV','pw_in_GeV')]
    if not np.all(np.isfinite(vals)): return False
    if not 0.05<vals[0]<100: return False
    if max(abs(vals[2]),abs(vals[5]))>3.05: return False
    # Stable station geometry; keep near-edge muons and explicitly report
    # collection outside finite wire bounds rather than censoring them.
    if max(abs(vals[1]),abs(vals[4]))>STATION_WIDTHS[station]/2+0.05: return False
    d=np.array([vals[4]-vals[1],vals[5]-vals[2],vals[6]-vals[3]])
    length=np.linalg.norm(d)
    if abs(length-vals[0])>max(0.03,0.01*vals[0]): return False
    moment=np.array(vals[9:12]);p=np.linalg.norm(moment)
    if not 0.11 < p < 1.e5: return False
    # A very slightly bent native step is permitted, but incompatible
    # momentum-vs-chord direction is not reinterpreted as a physical muon.
    if float(np.dot(d,moment)/(length*p))<.995: return False
    return True


def select(args):
    crossing_provenance,SOURCE,*_=native_interface(args.repo)
    out=Path(args.out);out.mkdir(parents=True,exist_ok=True)
    counters={}
    for channel in CHANNELS:
        # Per-channel memory limit: keep only N candidates for each station.
        heaps={st:[] for st in STATIONS}
        stats={st:dict(unique_muon_station=0,valid=0,rejected=0,
                       redundant_module_crossings=0) for st in STATIONS}
        for chunk in range(1,args.chunks+1):
            csvfile,original=crossing_provenance(channel,chunk)  # SHA verified
            best={st:{} for st in STATIONS}
            with csvfile.open(newline='',encoding='utf-8') as f:
                reader=csv.DictReader(f)
                if not NATIVE_REQUIRED.issubset(reader.fieldnames or []):
                    raise ValueError(f'Native CSV missing physical fields: {csvfile}')
                for r in reader:
                    st=int(r['station'])
                    if st not in STATIONS or abs(int(r['pdg']))!=13:continue
                    if not validate(r,st):
                        stats[st]['rejected']+=1
                        continue
                    stats[st]['valid']+=1
                    key=(r['event'],r['track_id'])
                    present=best[st].get(key)
                    # One selected gas chord for each (physical muon, station),
                    # even in overlap regions. Longest representative crossing.
                    if present is not None:
                        stats[st]['redundant_module_crossings']+=1
                    if present is None or (float(r['path_length_mm']),-int(r['hit_index']),-int(r['det_id'])) > (
                        float(present['path_length_mm']),-int(present['hit_index']),-int(present['det_id'])):
                        best[st][key]=r
            for st in STATIONS:
                for (event,trackid),row in best[st].items():
                    stats[st]['unique_muon_station']+=1
                    score=stable_hash('selected-100k-uniform-crossings-v1',
                                      channel,chunk,event,trackid,st)
                    record={x:row[x] for x in FIELDS if x not in ('source','chunk')}
                    record['source']=SOURCE[channel];record['chunk']=str(chunk)
                    heap=heaps[st]
                    item=(-score,tuple(str(record[f]) for f in FIELDS))
                    if len(heap)<args.n:
                        heapq.heappush(heap,item)
                    elif item>heap[0]:
                        heapq.heapreplace(heap,item)
            if chunk%5==0 or chunk==args.chunks:
                print(f'{channel}: scanned native chunk {chunk}/{args.chunks}',flush=True)
        for st in STATIONS:
            rows=sorted(heaps[st],key=lambda x:-x[0]); nfound=stats[st]['unique_muon_station']
            if len(rows)<args.n:
                raise RuntimeError(f'Insufficient MS{st} {channel} muons: '+
                                   f'{len(rows):,} eligible of {args.n:,} requested')
            path=out/f'{channel}_MS{st}_selected_{args.n}.csv.gz'
            if path.exists() and not args.overwrite:
                raise FileExistsError(f'Existing selection preserved: {path}')
            with gzip.open(path,'wt',newline='',encoding='utf-8',compresslevel=4) as f:
                w=csv.writer(f);w.writerow(FIELDS)
                for _,data in rows:w.writerow(data)
            stats[st]['n_selected']=len(rows)
            stats[st]['selected_sha256']=sha(path)
            stats[st]['selected_path']=str(path)
            stats[st]['source']='frozen 2M parent Geant4 native gas-hit exports'
            counters[f'{channel}_MS{st}']=stats[st]
        print(f'Finished selecting {args.n:,}/station for {channel}',flush=True)
    manifest={
       'protocol':'NA60_LUT500_VS_MIDLINE_UNIFORM_HASH_RANKED_V1',
       'generated_parent_events_per_channel':args.chunks*50000,
       'candidate_unit':'one physical muon track per station, longest valid module gas chord',
       'all_station_statistics':'one crossing per selected muon in station; not parent-aggregated',
       'selection':'smallest deterministic 64-bit BLAKE2b scores among all eligible tracks',
       'stations':list(STATIONS), 'channels':list(CHANNELS),'n_per_sample':args.n,
       'selection_seed_semantics':'fixed ID hash, pilot is first 10k of 100k selection',
       'samples':counters,
       'limitations':['Angle distribution inherits frozen Geant4 acceptance',
                       'Overlap modules reduced to longest gas chord per muon/station',
                       'Reject invalid/non-muon gas chords only; near-edge responses audited separately']
    }
    m=out/'selection_manifest.json'
    m.write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
    print('PASS selected:',m)


def read_lut_zip(path):
    with zipfile.ZipFile(path) as z:
        matches=[x for x in z.namelist() if x.endswith('/garfield_electrons.csv')]
        if len(matches)!=1:raise ValueError('LUT500 ZIP does not contain exactly one response CSV')
        with z.open(matches[0]) as f:
            txt=(line.decode('utf-8') for line in f)
            reader=csv.DictReader(txt)
            required={'iu','iv','collected','wire_index'}
            if not required.issubset(reader.fieldnames or []): raise ValueError('LUT columns missing')
            ii=[];vv=[];yes=[];delta=[]
            for row in reader:
                ii.append(int(row['iu']));vv.append(int(row['iv']))
                yes.append(int(row['collected']))
                delta.append(int(row['wire_index']) if int(row['collected']) else 0)
    if len(ii)!=450000:
        raise ValueError(f'Expected 450,000 original LUT500 trials; got {len(ii):,}')
    return build_lut_kernel(ii,vv,yes,delta)


class DonorGains:
    def __init__(self,filename):
        with np.load(filename,allow_pickle=False) as z:
            if not {'cell','gain'}.issubset(z.files):
                raise ValueError('Expected J/psi calibration keys: cell, gain')
            cells=z['cell'].astype('int32');self.gains=z['gain'].astype('float64')
        if len(cells)!=len(self.gains) or np.any((cells<0)|(cells>=900)) or \
                np.any(self.gains<0) or not np.all(np.isfinite(self.gains)):
            raise ValueError('Invalid calibration donor library')
        self.order=np.argsort(cells,kind='stable')
        self.counts=np.bincount(cells,minlength=900)
        if np.any(self.counts==0):raise ValueError('No gain donors for some LUT cells')
        self.start=np.r_[0,np.cumsum(self.counts[:-1])]

    def draw(self,cells,rng):
        cells=np.asarray(cells,dtype=np.int32)
        offsets=rng.integers(0,self.counts[cells])
        return self.gains[self.order[self.start[cells]+offsets]]


def write_events(path,rows):
    if not rows:return
    with gzip.open(path,'wt',newline='',encoding='utf-8',compresslevel=4) as f:
        w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)


def analyze_sample(path,output,st,args,*,birth_sample,lam,marks,kernel,donors):
    from wire_core import nearest_wire,cell_ids
    out=[]; nrequested=args.limit
    with gzip.open(path,'rt',newline='',encoding='utf-8') as f:
        reader=csv.DictReader(f)
        for i,row in enumerate(reader):
            if i>=nrequested:break
            channel=path.name.split('_MS')[0]
            idparts=(channel,row['chunk'],row['event'],row['track_id'],
                     row['station'],row['det_id'],row['hit_index'])
            pos,_=birth_sample(row,lam,marks,stable_hash('same-Heed-v1',*idparts))
            u=pos[:,0];v=pos[:,1];k,local=nearest_wire(u,st)
            cells=cell_ids(local,v)
            if donors is None:
                weights=np.ones(len(u),dtype=np.float64)
            else:
                weights=donors.draw(cells,rng_for('shared-Jpsi-donors-v1',*idparts))
            geo,lut,ge,le=compare_one(u,v,weights,st,kernel,
                                     rng=rng_for('LUT500-drift-v1',*idparts))
            g=observables(geo);l=observables(lut)
            d={k:row[k] for k in ('source','chunk','event','track_id','station','det_id',
                                  'theta_u_in_deg','theta_w_in_deg','path_length_mm')}
            d.update(n_primaries=len(u),n_source_cells=len(np.unique(k)),
                     total_gain_drawn=float(weights.sum()),
                     geometry_edge_charge=ge,lut_edge_charge=le)
            for prefix,obs in (('geometry',g),('lut',l)):
                d.update({f'{prefix}_{name}':value for name,value in obs.items()})
            out.append(d)
            if (i+1)%10000==0:print(f'{channel} MS{st}: {i+1:,}/{nrequested:,}',flush=True)
    if len(out)!=nrequested: raise RuntimeError(f'{path} has only {len(out):,} rows')
    write_events(output,out)
    return out


def estimate_cluster_se(records,geo_key,lut_key):
    """Parent-clustered sandwich approximation for paired probability diff."""
    groups=defaultdict(lambda:[0.,0])
    n=len(records)
    if n==0:return None
    diff=[]
    for r in records:
        d=int(r[lut_key])-int(r[geo_key]);diff.append(d)
        cluster=(r['chunk'],r['event'])
        g=groups[cluster];g[0]+=d;g[1]+=1
    mu=float(np.mean(diff));m=len(groups)
    if m<=1:return {'diff':mu,'se':None,'ci95':None,'n_parent_clusters':m}
    resid=np.array([s-mu*c for s,c in groups.values()],dtype=float)
    se=float(np.sqrt(m/(m-1)*np.sum(resid**2))/n)
    return {'diff':mu,'se':se,'ci95':[mu-1.96*se,mu+1.96*se],
            'n_parent_clusters':m}


def summarize(records,channel,station):
    n=len(records)
    out={'channel':channel,'station':station,'n_crossings':n,
         'source':'frozen native Geant4 plus Heed-calibrated mark sampler',
         'denominator':'all selected muon gas crossings including zero response',
         'threshold_comparison':{},'multiplicity':{},'angle_bins':[]}
    for threshold in (1,5,10):
        key=f'n_gt{threshold}pct'
        gr=np.array([int(r[f'geometry_{key}']) for r in records]);
        lu=np.array([int(r[f'lut_{key}']) for r in records]);
        for mult in (2,3):
            name=f'ge{mult}_anodes_gt{threshold}pct'
            gmean=float(np.mean(gr>=mult));lmean=float(np.mean(lu>=mult))
            paired=[dict(chunk=r['chunk'],event=r['event'],
                         geom=int(gr[i]>=mult),lut=int(lu[i]>=mult))
                    for i,r in enumerate(records)]
            out['threshold_comparison'][name]={
               'geometry_rate':gmean,'lut_rate':lmean,
               'lut_minus_geometry':estimate_cluster_se(paired,'geom','lut')}
    for method in ('geometry','lut'):
        multiplicities=np.array([int(r[f'{method}_n_gt5pct']) for r in records]);
        out['multiplicity'][method]={str(m):int(np.count_nonzero(multiplicities==m)) for m in range(5)}
        out['multiplicity'][method]['5+']=int(np.count_nonzero(multiplicities>=5))
        out['multiplicity'][method]['max_n_positive']=int(max(int(r[f'{method}_n_positive']) for r in records))
        out['multiplicity'][method]['rate_three_or_more_5pct']=float(np.mean(multiplicities>=3))
        out['multiplicity'][method]['mean_q2_fraction']=float(np.mean([r[f'{method}_q2_fraction'] for r in records]))
    th=np.array([abs(float(r['theta_u_in_deg'])) for r in records])
    for lo,hi in ((0,5),(5,10),(10,15),(15,20),(20,30),(30,180)):
        pick=np.where((th>=lo)&(th<hi))[0]
        if not len(pick):continue
        row=dict(theta_u_abs_deg=[lo,hi],n=int(len(pick)))
        for method in ('geometry','lut'):
            row[f'{method}_two_5pct']=float(np.mean([
                    int(records[i][f'{method}_n_gt5pct'])>=2 for i in pick]))
            row[f'{method}_three_5pct']=float(np.mean([
                    int(records[i][f'{method}_n_gt5pct'])>=3 for i in pick]))
        out['angle_bins'].append(row)
    out['total_primary_electrons']=int(sum(int(r['n_primaries']) for r in records))
    out['events_with_zero_drawn_gain']=int(sum(float(r['total_gain_drawn'])==0 for r in records))
    out['events_with_geometric_edge_loss']=int(sum(float(r['geometry_edge_charge'])>0 for r in records))
    out['events_with_lut_edge_loss']=int(sum(float(r['lut_edge_charge'])>0 for r in records))
    return out


def plot_summary(summary,path):
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        return
    rows=summary['angle_bins']; x=[(r['theta_u_abs_deg'][0]+r['theta_u_abs_deg'][1])/2 for r in rows]
    # use upper edge for rightmost open-ended bin only for plotting: left at 30.
    if x:x[-1]=30 if rows[-1]['theta_u_abs_deg'][1]==180 else x[-1]
    fig,ax=plt.subplots(figsize=(8,4.9))
    for method,style in (('geometry','--'),('lut','-')):
        ax.plot(x,[100*r[f'{method}_two_5pct'] for r in rows],
                linestyle=style,marker='o',label=f'{method}: >=2 anodes')
        ax.plot(x,[100*r[f'{method}_three_5pct'] for r in rows],
                linestyle=style,marker='s',label=f'{method}: >=3 anodes')
    ax.set(xlabel='Absolute local muon incidence angle |theta_u| [degrees]',
           ylabel='Muon crossings meeting 5% charge threshold [%]',
           title=f"{summary['channel']} / MS{summary['station']} — geometric vs LUT500")
    ax.grid(alpha=0.25);ax.legend(fontsize=8)
    fig.tight_layout();fig.savefig(path,dpi=150);plt.close(fig)


def compare(args):
    out=Path(args.out);sel=json.loads((out/'selection_manifest.json').read_text())
    if sel.get('n_per_sample')!=args.selected_n:
        raise ValueError('Mismatch between selection and requested selected size')
    crossing_provenance,SOURCE,birth_sample,HT,HC,read_heed=native_interface(args.repo)
    if not args.lut or not args.lut.is_file():
        raise FileNotFoundError('Original NA60_LUT500.zip required (--lut)')
    if args.weights=='donor' and (not args.calibration or not args.calibration.is_file()):
        raise FileNotFoundError('J/psi response NPZ required for --weights donor (--calibration)')
    lam,marks=read_heed(HT,HC)
    kernel=read_lut_zip(args.lut)
    donors=DonorGains(args.calibration) if args.weights=='donor' else None
    dst=out/f'comparison_{args.limit}_{args.weights}'
    if dst.exists() and not args.overwrite: raise FileExistsError(dst)
    dst.mkdir(parents=True,exist_ok=True)
    aggregate=[]
    for channel in CHANNELS:
        for st in STATIONS:
            sample=out/f'{channel}_MS{st}_selected_{args.selected_n}.csv.gz'
            original=sel['samples'][f'{channel}_MS{st}']
            if sha(sample)!=original['selected_sha256']:raise ValueError('Selection SHA changed: '+str(sample))
            eventsfile=dst/f'events_{channel}_MS{st}.csv.gz'
            records=analyze_sample(sample,eventsfile,st,args,birth_sample=birth_sample,
                                   lam=lam,marks=marks,kernel=kernel,donors=donors)
            summary=summarize(records,channel,st)
            (dst/f'summary_{channel}_MS{st}.json').write_text(json.dumps(summary,indent=2)+'\n')
            plot_summary(summary,dst/f'angle_{channel}_MS{st}.png')
            aggregate.append(summary)
            print(f'PASS: {channel} MS{st}, 2 anodes >5%: '
                  f"geometry {100*summary['threshold_comparison']['ge2_anodes_gt5pct']['geometry_rate']:.2f}%, "
                  f"LUT {100*summary['threshold_comparison']['ge2_anodes_gt5pct']['lut_rate']:.2f}%",flush=True)
    manifest={'method':'paired drift assignment with fixed common Heed starts/charge draws',
              'version':'NA60_P3_NOFW_MIDLINE_VS_LUT500_20261118',
              'execution':'Geant4 frozen 2M per channel; HEED-CALIBRATED MARKED-POISSON (not native TrackHeed per crossing)',
              'weight_semantics':'same Jpsi donor gain for geometry and LUT' if donors else 'unit primary charge (no gain)',
              'lut_conditional_on_collection':True,'notes':[
                 'Geometry and LUT share same electron starts and weights; only the physical anode assignment differs',
                 'One crossing per physical muon per station; multi-anode counts unrestricted at event level',
                 'LUT offsets relative per electron are -1,0,+1 only',
                 'Finite physical chamber edge losses tracked and not silently reassigned',
                 'No Shockley-Ramo cathode response or VMM thresholds',
                 'True TrackHeed primary coordinates would require separate native Heed-only export/replay'],
              'lut_sha256':sha(args.lut),'heed_tracks_sha256':sha(HT),
              'heed_clusters_sha256':sha(HC),'calibration_sha256':sha(args.calibration) if donors else None,
              'selection_manifest_sha256':sha(out/'selection_manifest.json'),
              'pilot_or_full':'pilot' if args.limit<args.selected_n else 'full',
              'n_crossings_per_channel_station':args.limit,
              'results':[{'channel':x['channel'],'station':x['station'],
                          'two_5pct':x['threshold_comparison']['ge2_anodes_gt5pct'],
                          'three_5pct':x['threshold_comparison']['ge3_anodes_gt5pct']}
                         for x in aggregate]}
    (dst/'run_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    table=[]
    for smry in aggregate:
        s2=smry['threshold_comparison']['ge2_anodes_gt5pct']
        s3=smry['threshold_comparison']['ge3_anodes_gt5pct']
        ci=s2['lut_minus_geometry']['ci95']
        table.append({'channel':smry['channel'],'station':f"MS{smry['station']}",
            'n_crossings':smry['n_crossings'],
            'geometry_ge2_5pct_rate':s2['geometry_rate'],
            'lut_ge2_5pct_rate':s2['lut_rate'],
            'diff_ge2_5pct_rate':s2['lut_minus_geometry']['diff'],
            'diff_ge2_5pct_ci95_low':ci[0] if ci else None,
            'diff_ge2_5pct_ci95_high':ci[1] if ci else None,
            'geometry_ge3_5pct_rate':s3['geometry_rate'],
            'lut_ge3_5pct_rate':s3['lut_rate'],
            'diff_ge3_5pct_rate':s3['lut_minus_geometry']['diff'],
            'max_n_positive_geometry':smry['multiplicity']['geometry']['max_n_positive'],
            'max_n_positive_lut':smry['multiplicity']['lut']['max_n_positive'],
            'n_lut_edge_losses':smry['events_with_lut_edge_loss']})
    with (dst/'comparison_table.csv').open('w',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=list(table[0]));writer.writeheader();writer.writerows(table)
    lines=['# LUT500 vs geometric midline — paired 5% threshold','',
           '| Channel | Station | N | Geometry ≥2 | LUT ≥2 | LUT−geometry (pp) | Geometry ≥3 | LUT ≥3 |',
           '|---|---|---:|---:|---:|---:|---:|---:|']
    for r in table:
        lines.append(f"| {r['channel']} | {r['station']} | {r['n_crossings']:,} | "
            f"{r['geometry_ge2_5pct_rate']*100:.3f}% | {r['lut_ge2_5pct_rate']*100:.3f}% | "
            f"{r['diff_ge2_5pct_rate']*100:+.3f} | "
            f"{r['geometry_ge3_5pct_rate']*100:.3f}% | {r['lut_ge3_5pct_rate']*100:.3f}% |")
    lines.extend(['','All rates use all selected muon crossings; zero-gain events are retained.',
       'The charged-primary model is Heed-calibrated marked Poisson, not exact per-track TrackHeed.',
       'Station-specific final MuID segmentation and induced cathode signals are not modeled.'])
    (dst/'comparison_table.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    print('PASS compared:',dst)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('mode',choices=('select','compare'))
    p.add_argument('--repo',type=Path,default=Path.home()/'na6p/src/NA6PRoot-mwpc')
    p.add_argument('--out',type=Path,default=Path.home()/'na6p/studies/lut500_vs_geometry_100k')
    p.add_argument('--chunks',type=int,default=40,help='50k-parent native Geant4 chunks per source')
    p.add_argument('--n',type=int,default=100000,help='sample size per channel-station for select')
    p.add_argument('--selected-n',type=int,default=100000)
    p.add_argument('--limit',type=int,default=10000,help='first N of selected sample to analyze')
    p.add_argument('--lut',type=Path,help='original NA60_LUT500.zip')
    p.add_argument('--calibration',type=Path,help='J/psi donor NPZ with cell,gain arrays')
    p.add_argument('--weights',choices=('donor','unit'),default='donor')
    p.add_argument('--overwrite',action='store_true',help='Explicitly replace existing output')
    args=p.parse_args()
    if not 1<=args.chunks<=40 or args.n<=0 or not 1<=args.limit<=args.selected_n:
        p.error('Expected chunks 1..40, n>0, and 1 <= limit <= selected-n')
    if args.mode=='select':select(args)
    else:compare(args)

if __name__=='__main__':main()