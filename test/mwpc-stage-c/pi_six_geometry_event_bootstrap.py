#!/usr/bin/env python3
"""First-chunk six-geometry MWPC *PARENT-EVENT* statistical uncertainties.

Re-folds the existing 3000 native Geant4 generated parents PER SOURCE
using the EXACT shared-Heed birth seed and already frozen Garfield LUTs.
This repeats Python LOOKUPS, not Geant4/Garfield physics simulations.
Parent event IDs (including generated events with NO chamber hit) are the
independent resampling unit. All chambers/hits/electrons arising from a
parent are weighted together in the bootstrap. Each source is kept
separate; no fake per-electron independence.

Stage scope: station-integrated ten-ns timing PDFs, conditional mean,
centered RMS, collection fraction, approach-angle and impact-angle
histograms. Per-anode and entrance-angle-conditioned parent-bootstrap
uncertainties require a subsequent more selective event-indexed summary.

Input must include the compact station reference, independently produced
by pi_six_geometry_native_dashboard_pages.py from the first full fold.
Fail closed if bootstrap base totals do not match the original full fold.
"""
from __future__ import annotations
import argparse
import csv
import json
import math
import os
import time
from pathlib import Path

import numpy as np

from pi_six_geometry_native_timing import (
    OUT, EXPORT, MAPS, HEED_TRACKS, HEED_CLUSTERS,
    SOURCES, FIRST_CHUNK_EXPECTED, GEOMETRY, PITCHES,
    stable_seed, birth_sample, response_for_births, load_map, sha, tag,
)
from run_pi_per_chamber_response import read_heed

N_PARENT=3000
NCFG=6
NST=6
T0=5
T1=T0+31
AN=T1
AS=AN+1
AQ=AS+1
AH0=AQ+1
AH1=AH0+18
PN=AH1
PH0=PN+1
NF=PH0+36

# Features: 0=Heed primaries; 1=anode collected; 2=known drift time;
# 3=sum time; 4=sum squared time; 5:36=31 time bins including overflow.
# 36:39=angle finite N,sum,sum2; 39:57=18 approach-angle bins.
# 57=phi finite N; 58:94=36 signed impact-angle bins.


def finite_moments(n,s,q):
    if n<=0:
        return None,None
    m=float(s/n)
    return m,float(math.sqrt(max(0.,q/n-m*m)))


def scalar_metrics(vec):
    p,col,known=vec[:3]
    mean,rms=finite_moments(known,vec[3],vec[4])
    amean,arms=finite_moments(vec[AN],vec[AS],vec[AQ])
    return dict(
        n_primaries=int(round(p)),
        n_anode_collected=int(round(col)),
        n_known_t=int(round(known)),
        anode_collection_fraction=float(col/p) if p>0 else None,
        mean_drift_ns=mean,
        centered_rms_drift_ns=rms,
        time_hist=[int(round(v)) for v in vec[T0:T1]],
        time_pdf_fraction=(vec[T0:T1]/known).tolist() if known>0 else None,
        n_known_approach_angle=int(round(vec[AN])),
        mean_approach_angle_deg=amean,
        centered_rms_approach_angle_deg=arms,
        angle_hist=[int(round(v)) for v in vec[AH0:AH1]],
        approach_angle_pdf_fraction=(
            vec[AH0:AH1]/vec[AN]).tolist() if vec[AN]>0 else None,
        n_known_impact_phi=int(round(vec[PN])),
        phi_hist=[int(round(v)) for v in vec[PH0:NF]],
        impact_phi_pdf_fraction=(
            vec[PH0:NF]/vec[PN]).tolist() if vec[PN]>0 else None,
    )


def ensemble_errors(boot):
    """Return bootstrap SD for ratios, arithmetic mean, and centered RMS."""
    p=boot[:,0]
    c=boot[:,1]
    n=boot[:,2]
    t=np.divide(boot[:,3],n,out=np.full_like(n,np.nan),
                where=n>0)
    v=np.divide(boot[:,4],n,out=np.full_like(n,np.nan),
                where=n>0)-t*t
    rms=np.sqrt(np.maximum(0.,v))
    angle_n=boot[:,AN]
    angmean=np.divide(boot[:,AS],angle_n,
        out=np.full_like(n,np.nan),where=angle_n>0)
    angv=np.divide(boot[:,AQ],angle_n,
        out=np.full_like(n,np.nan),where=angle_n>0)-angmean**2
    angrams=np.sqrt(np.maximum(0.,angv))
    coll=np.divide(c,p,out=np.full_like(n,np.nan),where=p>0)
    def sd(x):
        z=np.asarray(x)[np.isfinite(x)]
        return float(z.std(ddof=1)) if len(z)>=2 else None
    def bins(numer,den):
        fractions=np.divide(numer,den[:,None],
             out=np.full_like(numer,np.nan),where=den[:,None]>0)
        return [sd(fractions[:,j]) for j in range(numer.shape[1])]
    return dict(
        anode_collection_fraction=sd(coll),
        mean_drift_ns=sd(t),
        centered_rms_drift_ns=sd(rms),
        mean_approach_angle_deg=sd(angmean),
        centered_rms_approach_angle_deg=sd(angrms),
        time_pdf_fraction=bins(boot[:,T0:T1],n),
        approach_angle_pdf_fraction=bins(boot[:,AH0:AH1],angle_n),
        impact_phi_pdf_fraction=bins(boot[:,PH0:NF],boot[:,PN]),
    )


def preflight_maps(root):
    maps={}
    provenance={}
    for pitch in PITCHES:
        for fw in (0.,75.):
            key=tag(pitch,fw)
            pool,meta=load_map(pitch,fw,root)
            maps[key]=(pitch,fw,pool)
            provenance[key]=meta
    return maps,provenance


def event_cube(source,raw,maps,lam,marks):
    """Dense first-chunk source-event x station x config x sufficient stats."""
    cube=np.zeros((N_PARENT,NST,NCFG,NF),dtype=np.float64)
    seen=0
    map_entries=list(maps.items())
    with raw.open(newline="",encoding="utf-8") as file:
        reader=csv.DictReader(file)
        required={"event","hit_index","track_id","pdg","station","det_id",
                  "path_length_mm","u_in_mm","u_out_mm","v_in_mm",
                  "v_out_mm","w_in_mm","w_out_mm"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError("Missing original native segment fields")
        for row in reader:
            if abs(int(row["pdg"]))!=13:
                continue
            event=int(row["event"]);st=int(row["station"])
            det=int(row["det_id"])
            if not (0<=event<N_PARENT and 0<=st<NST and det in GEOMETRY
                    and GEOMETRY[det].station==st):
                raise ValueError("Invalid parent/installed TGeo module")
            seed=stable_seed("native-six-pitch-shared-Heed-v1",source,
                             "chunk0001",row["event"],row["hit_index"],
                             row["track_id"])
            xyz,q=birth_sample(row,lam,marks,seed)
            u,v=xyz[:,0],xyz[:,1]
            seen+=1
            for ci,(_, (pitch,fw,pool)) in enumerate(map_entries):
                finite,rawcol,wires,response=response_for_births(
                    u,v,q,pitch,st,pool)
                col=response[finite]
                stat=cube[event,st,ci]
                stat[0]+=len(q)
                stat[1]+=len(col)
                if not len(col):
                    continue
                t=col[:,2]
                good=np.isfinite(t)&(t>=0)
                t=t[good]
                stat[2]+=len(t)
                stat[3]+=float(t.sum())
                stat[4]+=float(np.square(t).sum())
                stat[T0:T1]+=np.bincount(
                    np.minimum((t//10).astype(np.int64),30),minlength=31)
                angle=col[:,3]
                angle=angle[np.isfinite(angle)&(angle>=0)&(angle<=180)]
                stat[AN]+=len(angle)
                stat[AS]+=float(angle.sum())
                stat[AQ]+=float(np.square(angle).sum())
                stat[AH0:AH1]+=np.bincount(
                    np.minimum((angle//10).astype(np.int64),17),
                    minlength=18)
                phi=col[:,4]
                phi=phi[np.isfinite(phi)&(phi>=-180)&(phi<=180)]
                stat[PN]+=len(phi)
                stat[PH0:NF]+=np.bincount(
                    np.clip(np.floor((phi+180)/10).astype(np.int64),0,35),
                    minlength=36)
    if seen!=FIRST_CHUNK_EXPECTED[source]:
        raise ValueError(f"Wrong native first chunk for {source}: "
                         f"{seen}/{FIRST_CHUNK_EXPECTED[source]}")
    # Independent sums over parent event reproduce exact collected, known,
    # and 10ns/angle histogram population.
    totals=cube.sum(axis=0)
    if (np.any(totals[:,:,1]>totals[:,:,0]) or
        not np.array_equal(totals[:,:,T0:T1].sum(axis=-1),
                           totals[:,:,2]) or
        not np.array_equal(totals[:,:,AH0:AH1].sum(axis=-1),
                           totals[:,:,AN]) or
        not np.array_equal(totals[:,:,PH0:NF].sum(axis=-1),
                           totals[:,:,PN])):
        raise AssertionError("Event-level Heed/Garfield sufficient stats broken")
    return cube,seen


def reference_check(source,totals,reference):
    """Check every station/config against the earlier FULL native 496k fold."""
    lookup={(g["source"],g["station"],g["config"]):g
            for g in reference["groups"]}
    configs=list(reference["meta"]["garfield_luts"])
    for st in range(NST):
        for ci,config in enumerate(configs):
            stored=lookup[(source,st,config)]
            vec=totals[st,ci]
            if (stored["primaries"]!=int(vec[0]) or
                stored["collected"]!=int(vec[1]) or
                stored["known"]!=int(vec[2]) or
                stored["hist"]!=[int(x) for x in vec[T0:T1]] or
                stored["angle_n"]!=int(vec[AN]) or
                stored["angle_hist"]!=[int(x) for x in vec[AH0:AH1]] or
                stored["phi_n"]!=int(vec[PN]) or
                stored["phi_hist"]!=[int(x) for x in vec[PH0:NF]]):
                raise AssertionError(
                    f"Full native fold differs from event-indexed fold: "
                    f"{source} MS{st} {config}")
            mean,rms=finite_moments(vec[2],vec[3],vec[4])
            for k,pred in (("mean_ns",mean),("rms_ns",rms)):
                actual=stored[k]
                if (actual is None)!=(pred is None) or (
                    actual is not None and abs(actual-pred)>1.e-6):
                    raise AssertionError(
                        f"Native fold mean/RMS mismatch: {source} MS{st}, "
                        f"{config}: {k} {actual} vs {pred}")
    return True


def source_bootstrap(source,cube,configs,replicas,seed):
    rng=np.random.default_rng(seed)
    # One multinomial count for each generated parent, shared over all
    # 6 stations and all 6 detector configurations in each bootstrap draw.
    weights=rng.multinomial(N_PARENT,np.full(N_PARENT,1./N_PARENT),
                            size=replicas).astype(np.float64)
    weighted=(weights @ cube.reshape(N_PARENT,-1)).reshape(
        replicas,NST,NCFG,NF)
    baseline=cube.sum(axis=0)
    result=[]
    for station in range(NST):
        for ci,config in enumerate(configs):
            row=dict(source=source,station=station,config=config,
                n_generated_parents=N_PARENT,
                statistic=scalar_metrics(baseline[station,ci]),
                bootstrap_sd=ensemble_errors(weighted[:,station,ci,:]),
            )
            result.append(row)
    return result


def run(args):
    t0=time.perf_counter()
    reference=json.loads(args.reference.read_text(encoding="utf-8"))
    maps,provenance=preflight_maps(args.maps)
    configs=list(maps)
    if list(reference["meta"]["garfield_luts"])!=configs:
        raise ValueError("Garfield geometry order inconsistent with reference")
    for key,v in provenance.items():
        if v["csv_sha256"]!=reference["meta"]["garfield_luts"][key]["csv_sha256"]:
            raise ValueError("LUT SHA differs from full native reference")
    if (sha(args.heed_tracks)!=reference["meta"]["heed_tracks_sha"] or
        sha(args.heed_clusters)!=reference["meta"]["heed_clusters_sha"]):
        raise ValueError("Heed library differs from original native fold")
    lam,marks=read_heed(args.heed_tracks,args.heed_clusters)
    args.out.mkdir(parents=True,exist_ok=True)
    rows=[]
    for source,channel in SOURCES.items():
        raw=args.export_dir/f"{channel}_geant_gas_crossings.csv"
        if sha(raw)!=reference["meta"]["native_geant4_export_sha256"][source]:
            raise ValueError("Native Geant4 source file has changed")
        began=time.perf_counter()
        cube,count=event_cube(source,raw,maps,lam,marks)
        reference_check(source,cube.sum(axis=0),reference)
        results=source_bootstrap(source,cube,configs,args.replicas,
                                 args.seed+list(SOURCES).index(source))
        rows.extend(results)
        print(f"PASS {source}: {count:,} native gas hits, "
              f"{N_PARENT:,} parent resampling units; "
              f"all 36 station/config distributions match full fold; "
              f"wall {time.perf_counter()-began:.1f}s",flush=True)
        del cube
    report={
        "model":"six frozen p=3,3.5,4 mm, +/-75um field wires; "
                "shared-Heed first-chunk Geant4 fold",
        "uncertainty":"BOOTSTRAP SD over 3000 GENERATED parent events per "
              "source, including parents without gas hits. Weights shared "
              "across all correlated gas segments, electrons, stations and "
              "six Garfield responses in one source. These are STATISTICAL "
              "errors only; not Garfield LUT finite-trial/model uncertainty.",
        "replicas":args.replicas,"bootstrap_seed":args.seed,
        "time_bins_ns":"[0,10),...,[290,300),[300,infinity)",
        "electron_angle_bins_deg":"[0,10),...,[170,180]",
        "impact_phi_bins_deg":"[-180,-170),...,[170,180]",
        "scope":"FIRST 3000 GENERATED PARENT DIMUON EVENTS PER CHANNEL; "
                 "SOURCE×STATION×CONFIG INTEGRATED only. "
                 "No chamber/anode/incident-muon-angle bootstrap yet.",
        "full_native_station_reference_sha256":sha(args.reference),
        "garfield_luts":provenance,
        "wall_seconds":time.perf_counter()-t0,
        "groups":rows,
    }
    target=args.out/"pi_six_geometry_parent_bootstrap.json"
    tmp=target.with_suffix(target.suffix+".tmp")
    tmp.write_text(json.dumps(report,ensure_ascii=False,indent=2,
                              allow_nan=False)+"\n",encoding="utf-8")
    os.replace(tmp,target)
    table=args.out/"pi_six_geometry_parent_bootstrap.csv"
    with table.open("w",newline="",encoding="utf-8") as f:
        w=csv.writer(f)
        w.writerow(["source","station","config","n_generated_parents",
                    "n_primary","n_collected","n_known_time",
                    "mean_ns","mean_error_ns","rms_ns","rms_error_ns",
                    "anode_collection_pct","collection_error_pp"])
        for r in rows:
            m,e=r["statistic"],r["bootstrap_sd"]
            def pct(v):return None if v is None else 100*v
            w.writerow([r["source"],r["station"],r["config"],N_PARENT,
                m["n_primaries"],m["n_anode_collected"],m["n_known_t"],
                m["mean_drift_ns"],e["mean_drift_ns"],
                m["centered_rms_drift_ns"],
                e["centered_rms_drift_ns"],
                pct(m["anode_collection_fraction"]),
                pct(e["anode_collection_fraction"])])
    print("PASS: native source/station six-geometry PARENT bootstrap; "
          "all base histograms verified against original first-chunk fold")
    print("Bootstrap JSON:",target)
    print("Compact PI summary CSV:",table)
    return report


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--reference",type=Path,
                    default=OUT/"pages/pi_six_geometry_station_reference.json")
    ap.add_argument("--maps",type=Path,default=MAPS)
    ap.add_argument("--export-dir",type=Path,default=EXPORT)
    ap.add_argument("--heed-tracks",type=Path,default=HEED_TRACKS)
    ap.add_argument("--heed-clusters",type=Path,default=HEED_CLUSTERS)
    ap.add_argument("--out",type=Path,default=OUT/"parent_bootstrap")
    ap.add_argument("--replicas",type=int,default=160)
    ap.add_argument("--seed",type=int,default=261021)
    args=ap.parse_args()
    if args.replicas<25:
        ap.error("At least 25 parent bootstrap replicas required")
    if not args.reference.is_file():
        ap.error("Run pi_six_geometry_native_dashboard_pages.py first "
                 "to make independent verified station reference")
    run(args)


if __name__=="__main__":
    main()
