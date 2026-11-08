#!/usr/bin/env python3
"""Read-only live performance monitor for 3x10k original-Geant4 Penning PBS runs.

Never launches/modifies physics. Reads only per-task verified reports, timestamps
and exit codes, so it is lightweight even for massive birth CSV archives.
"""
import argparse
import json
import time
from pathlib import Path

CHANNELS=("Jpsi","Phi","Omega")


def summary(root,channel):
    camp=root/f"original2m_{channel.lower()}_10000_v1"
    manifest=camp/"manifest.json"
    if not manifest.is_file():
        return None
    meta=json.loads(manifest.read_text())
    assert meta.get("channel")==channel and meta.get("n_muons")==10000
    done=0; nruns=0; failed=[]; active=[]; notstarted=0
    muons=0; nseed=0; caps=0; ionbirths=0; collected=0
    sumwall=0.; starts=[]; ends=[]
    for i in range(500):
        d=camp/"results"/f"chunk_{i:04d}"
        if not d.is_dir():
            notstarted+=1
            continue
        sp=d/"started_epoch.txt"
        if sp.is_file():
            starts.append(float(sp.read_text().strip()))
        ef=d/"exit_status.txt"
        if not ef.is_file():
            active.append(i)
            continue
        if ef.read_text().strip()!="0" or not (d/"verified_report.json").is_file():
            failed.append(i)
            continue
        rr=json.loads((d/"verified_report.json").read_text())
        assert rr["chunk_index"]==i and rr["n_muons"]==20
        assert rr["schema"]==meta["schema"]
        wall=d/"wall_seconds.txt"
        if not wall.is_file():
            failed.append(i);continue
        sec=float(wall.read_text().strip())
        assert 0<=sec<10000000
        if sp.is_file():
            ends.append(float(sp.read_text().strip())+sec)
        nruns+=1;muons+=20;sumwall+=sec
        nseed+=rr["n_true_heed_primaries"]
        caps+=rr["n_capped_avalanches"]
        ionbirths+=rr["ion_births"]
        collected+=rr["collected_electrons"]
    return dict(channel=channel,done=nruns,muons=muons,wall=sumwall,
                primaries=nseed,capped=caps,births=ionbirths,collected=collected,
                active=active,failed=failed,notstarted=notstarted,
                starts=starts,ends=ends)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--base",type=Path,
                   default=Path("/storage/cgrp/stanis/mwpc_penning057"))
    p.add_argument("--worker-slots",type=int,default=10,
                   help="Maximum desired concurrently running worker slots, for conditional ETA")
    a=p.parse_args()
    rows=[summary(a.base,ch) for ch in CHANNELS]
    known=[v for v in rows if v is not None]
    print("PENNING r=0.57 | original Geant4 muon crossings | all three channels")
    for channel,r in zip(CHANNELS,rows):
        if r is None:
            print(f"{channel:<5}: NOT PREPARED")
            continue
        avg=r["wall"]/r["muons"] if r["muons"] else None
        print(f"{r['channel']:<5}: {r['muons']:>5}/10000 muons "
              f"{r['done']:>3}/500 verified jobs, active {len(r['active'])}, "
              f"failed {len(r['failed'])}, unstarted {r['notstarted']}")
        if avg is not None:
            print(f"       {avg:.2f} worker seconds/muon, "
                  f"{r['primaries']} full Garfield avalanches, "
                  f"capped {r['capped']}, ion births {r['births']:,}")
        if r["failed"]:print("       failed indices:",r["failed"][:15])
    total_mu=sum(r["muons"] for r in known)
    total_work=sum(r["wall"] for r in known)
    starts=[v for r in known for v in r["starts"]]
    ends=[v for r in known for v in r["ends"]]
    if starts:
        origin=min(starts)
        last=max(ends) if ends else origin
        if any(r["active"] for r in known):last=max(last,time.time())
        elapsed=(last-origin)/3600
        print(f"Elapsed since first worker start: {elapsed:.2f} h "
              f"(completed span / ongoing run), completed muons: {total_mu}")
        if elapsed>0 and total_mu:
            print(f"Observed aggregate completed throughput: {total_mu/elapsed:.1f} muons/h "
                  f"(includes queue/PBS ramp-up)")
    if total_mu:
        avg=total_work/total_mu
        print(f"Mean cumulative worker time per muon: {avg:.2f} s "
              f"({total_work/3600:.2f} CPU-worker hours already completed)")
        if a.worker_slots>0:
            remain=30000-total_mu
            print(f"Conditional remaining compute at {a.worker_slots} fully occupied "
                  f"worker slots: ~{remain*avg/a.worker_slots/3600:.1f} h "
                  f"if future muons match observed mean (NOT a scheduler ETA).")
    print("No failed jobs are retried automatically; never overwrite a partial task directory.")


if __name__=="__main__":
    main()
