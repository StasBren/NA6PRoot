#!/usr/bin/env python3
"""Read-only development-only plots from SHA-audited full Garfield muons.

Plots (a) mean collected G with parent CI, (b) 41-bin whole-muon ion A
for hemisphere/TB, (c) ion-BIRTH phi using saved 128-point weighted samples.
No sampled birth distribution is ever labelled exact cathode strip signal.
Reserved final validation parents are never opened.
"""
from __future__ import annotations
import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import pi_adaptive_muon_spacetime_campaign_v1 as camp

HERE=Path(__file__).resolve()
LABELS={'p3_noFW':'3 mm / no FW','p3_FW75':'3 mm / FW75',
        'p4_noFW':'4 mm / no FW','p4_FW75':'4 mm / FW75'}

def build():
    source,meta,done=camp.read_completed()
    n=min(len(done)*camp.GROUP,meta['dev_parents'])
    if n<120 or n%camp.GROUP:
        raise RuntimeError('Need at least 120 SHA-valid development muons')
    rep=camp.REPORTS/f'development_{n:04d}_parents.json'
    if not rep.is_file():
        raise FileNotFoundError('Run the campaign report before plotting: '+str(rep))
    summary=json.loads(rep.read_text())
    if summary['n_parent_events']!=n:
        raise ValueError('Report parent cardinality mismatch')
    out=camp.ROOT/f'plots_development_{n:04d}'
    if out.exists():raise FileExistsError('Existing plot set preserved: '+str(out))
    # Deliberately perform expensive SHA source verification BEFORE output dir.
    files=[]
    for i in range(n//camp.GROUP):
        d=camp.chunk_dir(i)
        for name in ('muon_reference_crossing.csv','muon_reference_ion_birth_sample.csv'):
            files.append((d/name,camp.sha(d/name)))
    out.mkdir()
    allmuon={g:{'hemi':[],'TB':[],'G':[]} for g in camp.GEOMS}
    # Source is one muon gas crossing per independent parent; 4 geometries
    # from same parent are matched, but each parent counts once per curve.
    for i in range(n//camp.GROUP):
        with (camp.chunk_dir(i)/'muon_reference_crossing.csv').open(newline='') as f:
            for r in csv.DictReader(f):
                g=r['configuration']
                for col,key in (('A_ion_hemi','hemi'),('A_ion_TB','TB'),
                                ('G_eff_per_Heed_primary','G')):
                    if r[col]:allmuon[g][key].append(float(r[col]))
    # One single figure: gain means + 95% parent-bootstrap interval.
    fig,ax=plt.subplots(figsize=(8.2,4.8))
    mu=[];elo=[];ehi=[]
    for g in camp.GEOMS:
        d=summary['per_geometry'][g]['gain'];v=d['estimate']
        mu.append(v);elo.append(max(0.,v-d['CI95'][0]));ehi.append(max(0.,d['CI95'][1]-v))
    xs=np.arange(len(mu))
    ax.errorbar(xs,mu,yerr=[elo,ehi],fmt='o',capsize=4)
    ax.set_xticks(xs,labels=[LABELS[g] for g in camp.GEOMS],rotation=15)
    ax.set_yscale('log');ax.set_ylabel('Mean collected effective gain / Heed electron')
    ax.set_title(f'True microscopic muons: {n} independent parents; 95% bootstrap CI')
    ax.grid(axis='y',alpha=.2)
    fig.tight_layout()
    fig.savefig(out/'01_mean_effective_gain_CI.png',dpi=190)
    plt.close(fig)
    # Crossing-count-normalized probability P(A), NOT ion weighted.
    fig,axes=plt.subplots(2,2,figsize=(11,8.5),sharex=True,sharey=True)
    edges=np.linspace(-1,1,42)
    for i,p in enumerate((3,4)):
        for j,key in enumerate(('hemi','TB')):
            ax=axes[i,j]
            for geom,style in ((f'p{p}_noFW','-'),(f'p{p}_FW75','--')):
                x=np.asarray(allmuon[geom][key],float)
                if len(x):
                    count,_=np.histogram(x,bins=edges)
                    ax.step((edges[1:]+edges[:-1])/2,count/max(1,count.sum()),
                            where='mid',linestyle=style,label=LABELS[geom])
            ax.set_title(f'{p} mm: '+('hemisphere' if key=='hemi' else 'cathode-facing ±45°'))
            ax.set_xlim(-1,1);ax.grid(alpha=.15)
            ax.legend(fontsize=8)
            if j==0:ax.set_ylabel('Probability per bin / physical muon')
            if i==1:ax.set_xlabel('Whole-muon ion-birth asymmetry A')
    fig.suptitle(f'Full Garfield ion-birth sector asymmetry, {n} true parents (undefined excluded)')
    fig.tight_layout()
    fig.savefig(out/'02_true_muon_ion_asymmetry_distributions.png',dpi=190)
    plt.close(fig)
    # Angular density of each ionisation birth, reconstructed from the
    # weighted UNIFORM per-avalanche ion position reservoir, not complete cloud.
    nbins=72
    edges=np.linspace(-180,180,nbins+1)
    hist={g:np.zeros(nbins,float) for g in camp.GEOMS}
    for i in range(n//camp.GROUP):
        with (camp.chunk_dir(i)/'muon_reference_ion_birth_sample.csv').open(newline='') as f:
            for r in csv.DictReader(f):
                g=r['configuration'];pitch=3. if g.startswith('p3_') else 4.
                u=float(r['u_mm']);v=float(r['v_mm'])
                # Source C++ uses std::round: ties away from zero.
                if not all(map(math.isfinite,(u,v))):raise ValueError('Nonfinite birth position')
                z=u/pitch
                nearest=math.floor(z+.5) if z>=0 else math.ceil(z-.5)
                du=u-nearest*pitch
                phi=math.degrees(math.atan2(v,du))
                weight=float(r['reservoir_weight'])
                if not math.isfinite(weight) or weight<=0:raise ValueError('Invalid true ion sample weight')
                k=min(nbins-1,max(0,int((phi+180)*nbins/360)))
                hist[g][k]+=weight
    fig,axes=plt.subplots(2,1,figsize=(10,8),sharex=True)
    for i,p in enumerate((3,4)):
        ax=axes[i]
        for g,style in ((f'p{p}_noFW','-'),(f'p{p}_FW75','--')):
            y=hist[g]
            if y.sum()>0:ax.step((edges[1:]+edges[:-1])/2,y/y.sum(),where='mid',
                                 linestyle=style,label=LABELS[g])
        for a,b in ((45,135),(-135,-45)):
            ax.axvspan(a,b,alpha=.09)
        ax.set_ylabel('Ion-weighted probability per 5° bin')
        ax.set_title(f'{p} mm: true microscopic ion birth azimuth about nearest anode')
        ax.legend();ax.grid(alpha=.15)
    axes[-1].set_xticks(np.arange(-180,181,45))
    axes[-1].set_xlabel('Ion birth azimuth φ [degrees]; shaded = cathode-facing sectors')
    fig.suptitle('128-point/avalanche weighted sample, NOT exact space-time cloud or strip signal')
    fig.tight_layout()
    fig.savefig(out/'03_weighted_ion_birth_azimuth.png',dpi=190)
    plt.close(fig)
    info={'schema':'NA60_PI_SASHA_TRUE_MICRO_DEV_FIGURES_V1',
          'status':'development_only','n_independent_parents':n,
          'origin_sha256':[{'path':str(p),'sha256':s} for p,s in files],
          'report_sha256':camp.sha(rep), 'plotter_sha256':camp.sha(HERE),
          'plots':{p.name:camp.sha(p) for p in sorted(out.glob('*.png'))},
          'plot3_note':'Source C++ reservoir is capped at 128 points per seed, weighted by N_ion_births / n_saved. The exact ionisation-birth sector census is in crossing/seed outputs; the sampled azimuth histogram is approximate.'}
    (out/'manifest.json').write_text(json.dumps(info,indent=2)+'\n')
    print('PASS three SHA-source-audited development-only figures:',out)

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('mode',choices=('build',))
    args=ap.parse_args()
    build()
