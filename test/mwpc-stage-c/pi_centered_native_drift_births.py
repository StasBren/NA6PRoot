#!/usr/bin/env python3
"""Sample native Geant4/Heed electron BIRTH coordinates by muon angle.

STRATIFIED pilot across three sources, six stations, and four |theta_u|
strata. Reuses the same deterministic first-chunk Heed primary ionisation
draw from centered pilot, but simulates new TRUE microscopic Garfield drift
lines, rather than extrapolating between birth and endpoint.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import heapq
import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np
from mwpc_checkerboard_projection import build_chambers
from pi_native_centered_timing_pilot import BASE,SOURCES,EXPORT
from pi_centered_electron_sampler import stable_seed
from run_pi_per_chamber_response import read_heed

OUT=BASE/"pi_centered_native_drift_trajectories"
GEOM={c.det_id:c for group in build_chambers() for c in group}
CATEGORIES=((0.,2.,'0–2°'),(2.,5.,'2–5°'),
            (5.,10.,'5–10°'),(10.,180.01,'≥10°'))


def angle_group(theta):
    t=abs(float(theta))
    if not math.isfinite(t) or t>180.:
        raise ValueError("Invalid muon entrance angle")
    for i,(lo,hi,_) in enumerate(CATEGORIES):
        if lo<=t<hi:
            return i
    raise AssertionError("Unbinned physical entrance angle")


def stratified_selection(export_dir,limit):
    """Keep lowest stable hash per group, not first K events in input order."""
    heaps=defaultdict(list)
    exposure=defaultdict(int)
    for source,channel in SOURCES.items():
        path=export_dir/f"{channel}_geant_gas_crossings.csv"
        if not path.is_file():
            raise FileNotFoundError(path)
        with path.open(newline='',encoding='utf-8') as stream:
            reader=csv.DictReader(stream)
            fields={'event','hit_index','track_id','station','det_id','pdg',
                    'u_in_mm','u_out_mm','v_in_mm','v_out_mm',
                    'w_in_mm','w_out_mm','theta_u_in_deg','path_length_mm'}
            if not fields.issubset(reader.fieldnames or []):
                raise ValueError(f"Native export missing TGeo/angle columns: {path}")
            for raw in reader:
                if abs(int(raw['pdg']))!=13:
                    continue
                st=int(raw['station']);det=int(raw['det_id'])
                if det not in GEOM or GEOM[det].station!=st:
                    raise ValueError("Native detector/station mismatch")
                theta=float(raw['theta_u_in_deg'])
                angle=angle_group(theta)
                key=(source,st,angle)
                exposure[key]+=1
                score=stable_seed("PI_CENTERED_DRIFT_1D_ANGLE_STRATIFIED_V1",
                         source,st,raw['event'],raw['hit_index'],raw['track_id'])
                # source×station×angle gets its own equal-length curated sample.
                item=(-score,(int(raw['event']),int(raw['hit_index'])),raw)
                heap=heaps[key]
                if len(heap)<limit:
                    heapq.heappush(heap,item)
                elif score<-heap[0][0]:
                    heapq.heapreplace(heap,item)
    return {k:[entry[2] for entry in sorted(heap,reverse=True)]
            for k,heap in heaps.items()},dict(exposure)


def sample(args):
    lam,marks=read_heed(args.heed_tracks,args.heed_clusters)
    selection,exposure=stratified_selection(args.export_dir,args.hits_per_group)
    args.out.mkdir(parents=True,exist_ok=True)
    start_path=args.out/"native_heed_starts.csv"
    label_path=args.out/"native_heed_start_metadata.csv"
    start_cols=('electron_id','u_cell_mm','v_birth_mm','w_start_mm','seed')
    label_cols=('electron_id','source','station','det_id','event','hit_index',
                'track_id','theta_u_deg','theta_w_deg','angle_group',
                'u_physical_birth_mm','v_physical_birth_mm','w_physical_birth_mm',
                'wire_center_u_mm')
    counts=defaultdict(int)
    rejections=defaultdict(int)
    electron_id=0
    with start_path.open('w',newline='',encoding='utf-8') as fs, \
         label_path.open('w',newline='',encoding='utf-8') as fm:
        ws=csv.DictWriter(fs,fieldnames=start_cols);ws.writeheader()
        wm=csv.DictWriter(fm,fieldnames=label_cols);wm.writeheader()
        for (source,st,group),rows in sorted(selection.items()):
            for row in rows:
                d=float(row['path_length_mm'])
                if not math.isfinite(d) or d<=0:
                    rejections['invalid_path']+=1;continue
                r=np.random.default_rng(stable_seed(
                    "centered_native_pilot_v1",source,"chunk0001",
                    int(row['event']),int(row['hit_index']),int(row['track_id'])))
                nc=int(r.poisson(lam*d))
                cluster=marks[r.integers(len(marks),size=nc)]
                fractions=np.repeat(r.random(nc),cluster)
                xyz0=np.array([float(row[f"{x}_in_mm"]) for x in 'uvw'])
                xyz1=np.array([float(row[f"{x}_out_mm"]) for x in 'uvw'])
                if not np.all(np.isfinite(xyz0)) or not np.all(np.isfinite(xyz1)):
                    raise ValueError("Invalid Geant4 gas entrance/exit positions")
                if max(abs(xyz0[1]),abs(xyz1[1]))>3.05:
                    raise ValueError("Wrong gas half-gap in native Geant4 input")
                nwire=77 if st==0 else 128
                # Up to K physical electron origins per muon; no conditioning
                # on collection and no use of old Garfield endpoint statuses.
                for j,fraction in enumerate(fractions[:args.electrons_per_hit]):
                    u,v,w=xyz0+fraction*(xyz1-xyz0)
                    nearest=int(math.floor(u/4+(nwire-1)/2+.5))
                    if not 0<=nearest<nwire:
                        rejections['outside_finite_grid']+=1;continue
                    center=(nearest-(nwire-1)/2)*4.
                    relative_u=u-center
                    if not (-2.0001<=relative_u<=2.0001):
                        raise AssertionError("Physical nearest-wire u phase")
                    if not -3<v<3:
                        rejections['on_cathode_boundary']+=1;continue
                    if math.hypot(relative_u,v)<.0151:
                        rejections['inside_anode_metal']+=1;continue
                    seed=int(stable_seed("PI_ACTUAL_MICRO_DRIFT_TRACK_V1",
                        source,st,int(row['event']),int(row['hit_index']),
                        int(row['track_id']),j)%(2**31-1))+1
                    ws.writerow(dict(electron_id=electron_id,
                        u_cell_mm=f"{relative_u:.12g}",v_birth_mm=f"{v:.12g}",
                        w_start_mm=0,seed=seed))
                    wm.writerow(dict(electron_id=electron_id,source=source,
                        station=st,det_id=int(row['det_id']),
                        event=int(row['event']),hit_index=int(row['hit_index']),
                        track_id=int(row['track_id']),theta_u_deg=row['theta_u_in_deg'],
                        theta_w_deg=row.get('theta_w_in_deg',''),
                        angle_group=group,
                        u_physical_birth_mm=f"{u:.12g}",
                        v_physical_birth_mm=f"{v:.12g}",
                        w_physical_birth_mm=f"{w:.12g}",
                        wire_center_u_mm=f"{center:.12g}"))
                    electron_id+=1
                    counts[f"{source}|MS{st}|{CATEGORIES[group][2]}"]+=1
    if electron_id==0:
        raise ValueError("No real native Heed electron starting points selected")
    for path in (start_path,label_path):
        if not path.is_file():
            raise FileNotFoundError(path)
    manifest={
       'geometry_id':'MWPC_GARFIELD_CENTERED_3p0_3p0_B0_V1',
       'pilot':'STRATIFIED SAMPLE of true native Geant4 first-chunk muon angles',
       'scope':'First 3000 parent events/source, 3 sources, 6 stations, '
               'up to K muon gas-hit segments per angle stratum; '
               'up to M Heed electrons from each selected hit',
       'not_unbiased_full_campaign':True,
       'selection':'Lowest stable SHA256 ranks by source×station×|theta_u|; '
                   'per-hit cap on Heed primary electrons',
       'nominal_anode_pitch_mm':4.,'cathode_half_gap_mm':3.,
       'B_T':0.,'hits_per_group':args.hits_per_group,
       'electrons_per_hit':args.electrons_per_hit,
       'native_muon_hit_exposures':{
           f'{s}|MS{st}|{CATEGORIES[ang][2]}':n
           for (s,st,ang),n in exposure.items()},
       'electron_starts_per_group':dict(counts),
       'skipped_start_counts':dict(rejections),
       'generated_native_Heed_electron_starts':electron_id,
       'native_heed_starts_sha256':hashlib.sha256(start_path.read_bytes()).hexdigest(),
       'metadata_sha256':hashlib.sha256(label_path.read_bytes()).hexdigest(),
    }
    (args.out/"native_heed_trace_manifest.json").write_text(
         json.dumps(manifest,indent=2,ensure_ascii=False)+"\n",
         encoding='utf-8')
    print(f"PASS: {electron_id:,} true native-angle Heed electron birth "
          "positions for NEW full microscopic Garfield trajectories")
    print("Heed/Geant4 starts:",start_path)
    return manifest


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--export-dir',type=Path,default=EXPORT)
    p.add_argument('--heed-tracks',type=Path,
      default=BASE/"c_heed_library_2p4/heed_2p4_tracks.csv")
    p.add_argument('--heed-clusters',type=Path,
      default=BASE/"c_heed_library_2p4/heed_2p4_clusters.csv")
    p.add_argument('--hits-per-group',type=int,default=8)
    p.add_argument('--electrons-per-hit',type=int,default=12)
    p.add_argument('--out',type=Path,default=OUT)
    args=p.parse_args()
    if args.hits_per_group<1 or args.electrons_per_hit<1:
        p.error("Pilot sample sizes must be positive")
    sample(args)


if __name__=='__main__':
    main()
