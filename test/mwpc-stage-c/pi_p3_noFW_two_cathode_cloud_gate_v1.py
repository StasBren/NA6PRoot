#!/usr/bin/env python3
"""SHA-audited small complete real Garfield ion cloud for 2 whole-cathode tests.

Select one full p3_noFW cloud (512..1024 real ion birth points), preserve all
positions/times, and generate nested unbiased macroion subsamples (32/128/512).
Run separate C++ RKF with identical actual p3_noFW wire and two entire cathode
weighting fields; compare FULL-source true-ion drift to compressed macroions.
Original D500, Garfield cloud chunks, ion births are READ ONLY.
No U/V strips; Ar+ in Ar mobility used as approximate gas model.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,random
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

REPO=Path(__file__).resolve().parents[2]
BASE=REPO/'test_runs/mwpc_stage_c'
SOURCE=BASE/'pi_p3_fullcloud_learning5_v1'/'full_cloud_chunks'
OUT=BASE/'pi_p3_noFW_two_cathode_cloud_gate_v1'
SIZES=(32,128,512)
KEY=('crossing_id','configuration','primary_index')
BIRTH=('u_mm','v_mm','w_mm','t_ns')
def must(ok,msg):
    if not ok:raise RuntimeError(msg)
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()
def read(path):
    with path.open(newline='') as f:
        yield from csv.DictReader(f)
def identity(row):
    return tuple(row[x] for x in KEY)
def find_cloud(lo,hi,search_chunks):
    for j in range(search_chunks):
        d=SOURCE/f'chunk_{j:04d}'
        seal_path=d/'complete.json'
        must(seal_path.is_file(),f'Missing sealed full-cloud chunk {j}: '+str(d))
        seal=json.loads(seal_path.read_text())
        must(seal['chunk']==j and seal['schema']=='NA60_P3_FULL_CLOUD_SHA_CHUNK_V1',
             'Invalid full-cloud provenance')
        trial=d/'cloud_trials.csv'
        birth=d/'cloud_birth_reservoir.csv'
        must(sha(trial)==seal['sha256_files'][trial.name]
             and sha(birth)==seal['sha256_files'][birth.name],
             'SHA mismatch in original exact Full-Cloud Garfield archives')
        candidates={}
        for row in read(trial):
            if row['configuration']=='p3_noFW':
                n=int(row['ion_births'])
                if lo<=n<=hi:
                    candidates[identity(row)]=n
        if not candidates:continue
        # Stream file; each record is an original microscopic ion birth,
        # NOT a 128-position reservoir from 11,000 real muons.
        found={}
        for row in read(birth):
            key=identity(row)
            if key in candidates:
                must(float(row['reservoir_weight'])==1.,
                     'FullCloud archive is unexpectedly sampled/weighted')
                found.setdefault(key,[]).append([float(row[x]) for x in BIRTH])
        for key,n in candidates.items():
            if key in found:
                must(len(found[key])==n,'FullCloud ion birth census incomplete')
                return np.array(found[key]),dict(chunk=j,full_ion_births=n,
                    Garfield_crossing=key[0],geometry=key[1],primary_index=key[2],
                    source_chunk_seal_sha256=sha(seal_path),
                    cloud_trials_sha256=sha(trial),
                    cloud_births_sha256=sha(birth))
    raise RuntimeError(f'No original genuine full Garfield p3_noFW cloud in [{lo},{hi}] births among first {search_chunks} chunks')
def csvwrite(path,points,weight):
    with path.open('x',newline='') as f:
        w=csv.writer(f)
        w.writerow(('u_mm','v_mm','w_mm','t_ns','weight'))
        for row in points:w.writerow([*(f'{x:.15g}' for x in row),f'{weight:.15g}'])
def prepare(lo,hi,chunks):
    must(not OUT.exists(),'Refusing to overwrite real-cloud reference '+str(OUT))
    data,info=find_cloud(lo,hi,chunks)
    must(len(data)<=hi and len(data)>=lo,'Incorrect exact full cloud')
    N=len(data)
    # Same RNG permutation gives nested subsets, reducing extra noise
    # when comparing 32/128/512 macroion response.
    chosen=random.Random(20261108).sample(range(N),N)
    OUT.mkdir(parents=True,exist_ok=False)
    for k in SIZES:
        must(k<=N,'Full reference must contain at least 512 original ions')
        csvwrite(OUT/f'ions_{k:04d}.csv',data[chosen[:k]],N/k)
    csvwrite(OUT/'ions_full.csv',data,1.)
    info.update(schema='P3_NOFW_FULL_CLOUD_TWO_CATHODE_ION_GATE_V1',
                nested_sample_seed=20261108, sizes=SIZES,
                model='Ion RKF drift and whole-cathode Ramo signal (C++); no electron current',
                ion_mobility='Ar+ in Ar provisional; mixture time calibration not established',
                input_file_sha256={p.name:sha(p) for p in OUT.iterdir()})
    (OUT/'source_receipt.json').write_text(json.dumps(info,indent=2)+'\n')
    print('PASS frozen TRUE full Garfield cloud extracted (no Garfield rerun):',
          N,'ions, chunk',info['chunk'])
    print('PASS nested unbiased 32/128/512 samples plus ALL original full ion points')
    print('OUTPUT',OUT)
def compare():
    info=json.loads((OUT/'source_receipt.json').read_text())
    full=OUT/'full_waveforms.csv'
    must(full.is_file(),'Run the full two-cathode Garfield drift benchmark first')
    def waveform(k):
        p=OUT/f'{k}_waveforms.csv'
        with p.open(newline='') as f:
            r=list(csv.DictReader(f))
        a=np.asarray([[float(x['time_ns']),float(x['plus_ion_fC_per_ns']),
                       float(x['minus_ion_fC_per_ns'])] for x in r])
        return a
    ref=waveform('full')
    report={'schema':'P3_NOFW_ION_MACROION_CONVERGENCE_V1',
            'source_receipt':info,'error_metrics':{},'source_sha256':sha(OUT/'source_receipt.json')}
    fig,axs=plt.subplots(2,1,figsize=(10,7),sharex=True)
    for j,side in enumerate(('plus','minus'),start=1):
        axs[j-1].plot(ref[:,0]/1000.,ref[:,j],label='ALL real Garfield ion starts',lw=1.7)
    for k in SIZES:
        x=waveform(f'{k:04d}')
        must(x.shape==ref.shape and np.allclose(x[:,0],ref[:,0],atol=1e-9),
             'Different ion signal time grids; cannot compare')
        report['error_metrics'][str(k)]={}
        for j,side in enumerate(('plus','minus'),start=1):
            actual=ref[:,j];pred=x[:,j]
            # Integrated signed charge alone can conceal incorrect time shape.
            l1=float(np.sum(np.abs(pred-actual))/max(1e-40,np.sum(np.abs(actual))))
            dt=ref[1,0]-ref[0,0]
            Q0=float(np.sum(actual)*dt);Qk=float(np.sum(pred)*dt)
            qp=float(abs(Qk-Q0)/max(1e-40,abs(Q0)))
            report['error_metrics'][str(k)][side]=dict(L1_relative_waveform_error=l1,
                relative_integrated_signed_charge_error=qp,
                Q_full_fC=Q0,Q_macro_fC=Qk)
            axs[j-1].plot(x[:,0]/1000.,pred,alpha=.6,lw=1,
                          label=f'{k} weighted macroions')
    for j,side in enumerate(('Upper plane','Lower plane')):
        axs[j].set_ylabel(f'{side} current [fC/ns]')
        axs[j].legend(fontsize=8);axs[j].grid(alpha=.12)
    axs[-1].set_xlabel('Time [microseconds]')
    fig.suptitle('Two complete cathodes: exact ion starts vs spatial subsampling')
    fig.tight_layout()
    fn=OUT/'two_cathode_waveform_compression.png'
    must(not fn.exists(),'Refusing to overwrite previous graph')
    fig.savefig(fn,dpi=180);plt.close(fig)
    report['output_plot_sha256']=sha(fn)
    p=OUT/'compression_comparison.json'
    must(not p.exists(),'Refusing to overwrite previous comparison')
    p.write_text(json.dumps(report,indent=2)+'\n')
    print('PASS full-vs-macroion ion-current comparison',report['error_metrics'])
    print('OUTPUT',fn)
def selftest():
    rng=random.Random(3)
    z=rng.sample(range(50),50)
    assert len(set(z[:32]))==32 and set(z[:32]).issubset(z[:40])
    assert identity(dict(zip(KEY,('LUT3','p3_noFW','0')))==('LUT3','p3_noFW','0')
    print('PASS exact-original-full-ion selection, nested subsampling semantics')
def main():
    a=argparse.ArgumentParser(description=__doc__)
    a.add_argument('mode',choices=('selftest','prepare','compare'))
    a.add_argument('--min-births',type=int,default=512)
    a.add_argument('--max-births',type=int,default=1024)
    a.add_argument('--scan-chunks',type=int,default=80)
    z=a.parse_args()
    if z.mode=='selftest':selftest()
    elif z.mode=='prepare':prepare(z.min_births,z.max_births,z.scan_chunks)
    else:compare()
if __name__=='__main__':main()
