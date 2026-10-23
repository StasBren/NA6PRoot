#!/usr/bin/env python3
"""Chunk-aware, station-sharded six-Garfield native electron-response fold.

The input is 50,000 real Geant4 dimuon parents/source/chunk, not 3,000.
This new pipeline REUSES the original validated Heed microbirth sampler
and all six existing Garfield response maps. NO new Geant4/Garfield.

Process one (source, Geant4 chunk, station) at a time to bound RAM. Write
gzip-compressed sparse group counters with unrounded electron moments;
merge COMPLETE 2x50k/source station shards to source-station dashboards.
Keep source+chunk+local event number in stable Heed seed to prevent
correlations or collisions when native event numbers reset per chunk.
No per-electron statistical errors claimed: later parent-level bootstrap.
"""
from __future__ import annotations
import argparse
import csv
import gzip
import hashlib
import json
import os
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
from pi_six_geometry_native_timing import (
    MAPS, HEED_TRACKS, HEED_CLUSTERS, PITCHES, SOURCES, ANGLES,
    Count, aggregate, validate_conservation, load_map, tag, sha, REQUIRED_HIT,
)
from pi_2m_export_completed_chunks import DEST as EXPORT, SOURCE
from pi_six_geometry_layout_v3 import render as render_v3, make_portal
from pi_six_geometry_native_timing_html import build as render_detail
from run_pi_per_chamber_response import read_heed
from run_pi_geant4_campaign import REPO

OUTPUT=REPO/"test_runs/mwpc_stage_c/pi_six_geometry_2m_timing_v1"
CHANNELS=("Jpsi","Phi","Omega")
VERSION="PI_NATIVE_CHUNK_AWARE_2M_SIX_GEOM_V1"
FIELDS=("segments","primaries","collected","known","unknown","edge")
MOMENTS={
    "summed":"sum_drift_ns","squared":"sum_drift_ns2",
    "incidence_sum":"sum_incidence_deg",
    "incidence_squared":"sum_incidence_deg2",
    "phi_sum":"sum_phi_deg","phi_squared":"sum_phi_deg2",
}
COUNTS={"incidence_n":"angle_n","phi_n":"phi_n"}


def sha_file(path):
    # Constant RAM SHA256 even when CSV input is >300 MiB.
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda:f.read(4*1024*1024),b""):
            h.update(block)
    return h.hexdigest()


def config_pools():
    pools={}
    provenance={}
    for p in PITCHES:
        for d in (0.,75.):
            k=tag(p,d)
            pool,details=load_map(p,d,MAPS)
            pools[k]=(p,d,pool)
            provenance[k]=details
    if len(pools)!=6:
        raise AssertionError("Six independent Garfield LUTs required")
    return pools,provenance


def crossing_provenance(channel,index):
    folder=EXPORT/channel.lower()/f"chunk{index:04d}"
    infofile=folder/"crossings_manifest.json"
    csvfile=folder/"geant_muon_crossings.csv"
    if not infofile.is_file() or not csvfile.is_file():
        raise FileNotFoundError(f"Missing validated real Geant4 chunk export: {folder}")
    doc=json.loads(infofile.read_text(encoding="utf-8"))
    frozen=doc["immutable_input"]
    if (frozen["channel"]!=channel or frozen["source"]!=SOURCE[channel] or
        frozen["chunk"]!=index or frozen["n_generated_parent_events"]!=50000 or
        frozen["global_parent_index_interval"]!=[(index-1)*50000,index*50000]):
        raise ValueError(f"Wrong source/chunk/generated-event provenance: {infofile}")
    if sha_file(csvfile)!=doc["csv_sha256"]:
        raise ValueError(f"Native CSV SHA changed after verification: {csvfile}")
    return csvfile,doc


def snapshot(bucket):
    row=bucket.done()
    for attribute,name in MOMENTS.items():
        row[name]=float(getattr(bucket,attribute))
    return row


def restore(bucket,row):
    for k in FIELDS:
        setattr(bucket,k,getattr(bucket,k)+int(row[k]))
    bucket.hist+=np.asarray(row["hist"],dtype=np.int64)
    bucket.angle_hist+=np.asarray(row["angle_hist"],dtype=np.int64)
    bucket.phi_hist+=np.asarray(row["phi_hist"],dtype=np.int64)
    for attribute,name in MOMENTS.items():
        setattr(bucket,attribute,getattr(bucket,attribute)+float(row[name]))
    bucket.incidence_n+=int(row["angle_n"])
    bucket.phi_n+=int(row["phi_n"])


def as_row(key,bucket):
    src,cfg,st,det,wire,axis,anglebin=key
    return dict(source=src,config=cfg,station=st,det_id=det,
                anode=wire,axis=axis,anglebin=anglebin,**snapshot(bucket))


def key(row):
    return (row["source"],row["config"],row["station"],
            row["det_id"],row["anode"],row["axis"],row["anglebin"])


def frozen_fold_inputs(csvmeta,provenance,channel,index,station,limit):
    return dict(
        model=VERSION,source=SOURCE[channel],chunk=index,station=station,
        generated_parent_events_in_input_chunk=50000,
        pilot_muon_gas_segment_limit=limit,
        native_csv_sha256=csvmeta["csv_sha256"],
        native_input_root_sha256=csvmeta["immutable_input"]["input_root_sha256"],
        gas_crossings=csvmeta["n_gas_crossings"],
        muon_gas_crossings=csvmeta["n_muon_gas_crossings"],
        heed_tracks_sha256=sha_file(HEED_TRACKS),
        heed_clusters_sha256=sha_file(HEED_CLUSTERS),
        garfield_map_provenance=provenance,
        fold_source_sha256=sha_file(
            Path(__file__).parent/"pi_six_geometry_native_timing.py"),
        fold_driver_sha256=sha_file(Path(__file__)),
        rng_semantics="stable_seed('native-six-pitch-shared-Heed-v1',"
                      "source,chunk,local_event,hit_index,track_id);"
                      " shared births and q across six pitch fields",
    )


def shard_path(out,channel,index,station):
    return out/channel.lower()/f"MS{station}"/f"chunk{index:04d}.json.gz"


def fold_one(channel,index,station,limit,pools,provenance,lam,marks,out):
    path=shard_path(out,channel,index,station)
    manifest=path.with_suffix(".manifest.json")
    csvpath,csvmeta=crossing_provenance(channel,index)
    immutable=frozen_fold_inputs(csvmeta,provenance,channel,index,station,limit)
    if manifest.exists():
        saved=json.loads(manifest.read_text(encoding="utf-8"))
        if saved.get("immutable_inputs")!=immutable or not path.is_file():
            raise ValueError(f"Incompatible existing fold checkpoint: {manifest}")
        if sha_file(path)!=saved.get("gzip_sha256"):
            raise ValueError(f"Corrupted folded native station data: {path}")
        print(f"SKIP native Heed+6 LUT {channel} chunk{index:04d} MS{station}: "
              f"{saved['n_muon_gas_hits']:,} muon gas segments",flush=True)
        return
    if path.exists():
        raise ValueError(f"Unmanifested fold shard: {path}")
    def rows():
        with csvpath.open(newline="",encoding="utf-8") as f:
            reader=csv.DictReader(f)
            if not REQUIRED_HIT.issubset(reader.fieldnames or []):
                raise ValueError(f"Unexpected native TGeo CSV schema: {csvpath}")
            for row in reader:
                if int(row["station"])!=station or abs(int(row["pdg"]))!=13:
                    continue
                event=int(row["event"])
                if not 0<=event<50000:
                    raise ValueError(f"Unphysical local parent event ID: {event}")
                yield SOURCE[channel],row
    t0=time.perf_counter()
    groups,processed=aggregate(rows(),pools,lam,marks,limit=limit,
             chunk_key=f"chunk{index:04d}",source_overview=False,
             station_filter=station)
    validate_conservation(groups)
    n=processed.get(SOURCE[channel],0)
    if n<1:
        raise ValueError(f"No muon hits for station MS{station} in {csvpath}")
    rows_serial=[as_row(k,b) for k,b in sorted(groups.items())]
    path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_name(path.name+".incomplete")
    with gzip.open(tmp,"wt",encoding="utf-8",compresslevel=4) as f:
        json.dump({"meta":immutable,"group_count":len(rows_serial),
                   "groups":rows_serial},f,separators=(",",":"),
                  ensure_ascii=False,allow_nan=False)
    os.replace(tmp,path)
    meta=dict(immutable_inputs=immutable,gzip_sha256=sha_file(path),
              n_muon_gas_hits=n,n_groups=len(rows_serial),
              n_generated_parent_events_in_chunk=50000,
              measured_fold_seconds=time.perf_counter()-t0,
              description="Heed/LUT microscopic sampled electron statistics;"
                          " do not treat individual electrons as parent events")
    tmp_m=manifest.with_name(manifest.name+".tmp")
    tmp_m.write_text(json.dumps(meta,indent=2,ensure_ascii=False)+"\n",
                     encoding="utf-8")
    os.replace(tmp_m,manifest)
    print(f"PASS folded {channel} chunk{index:04d} MS{station}: "
          f"{n:,} actual muon gas hits, {len(rows_serial):,} groups, "
          f"{path.stat().st_size/1048576:.1f} MiB gzip, "
          f"{meta['measured_fold_seconds']:.1f}s",flush=True)


def merge(out,limit,number):
    if number!=2:
        raise ValueError("PI benchmark merge requires exactly 2 full Geant4 chunks/source")
    if limit:
        raise ValueError("PILOT shard cannot be labelled a 100k-parent timing report")
    portal=out/"pages"
    portal.mkdir(parents=True,exist_ok=True)
    all_refs=[]
    for channel in CHANNELS:
        for station in range(6):
            combined=defaultdict(Count)
            provenance=[]
            hits=0
            for index in range(1,number+1):
                path=shard_path(out,channel,index,station)
                manifest=path.with_suffix(".manifest.json")
                if not manifest.is_file() or not path.is_file():
                    raise FileNotFoundError(f"Missing complete native shard: {path}")
                m=json.loads(manifest.read_text(encoding="utf-8"))
                if m["immutable_inputs"]["pilot_muon_gas_segment_limit"]!=0 or \
                   m["immutable_inputs"]["source"]!=SOURCE[channel] or \
                   m["immutable_inputs"]["station"]!=station or \
                   m["immutable_inputs"]["chunk"]!=index or \
                   m["n_generated_parent_events_in_chunk"]!=50000 or \
                   sha_file(path)!=m["gzip_sha256"]:
                    raise ValueError(f"Inconsistent merge inputs: {path}")
                with gzip.open(path,"rt",encoding="utf-8") as f:
                    shard=json.load(f)
                if shard["group_count"]!=m["n_groups"]:
                    raise ValueError(f"Truncated native Heed/LUT fold: {path}")
                for row in shard["groups"]:
                    if row["source"]!=SOURCE[channel] or row["station"]!=station:
                        raise ValueError("Shard contains other source/station")
                    restore(combined[key(row)],row)
                hits+=m["n_muon_gas_hits"]
                provenance.append(m)
            validate_conservation(combined)
            if hits<=0:
                raise ValueError("Zero total muon gas hit source exposure")
            meta={
                "model_version":VERSION,
                "data_scope":"SIX-GEOMETRY NATIVE Geant4/Heed/Garfield, "
                             "100,000 GENERATED PARENTS PER SOURCE",
                "n_generated_parent_events_per_source":number*50000,
                "gas_hit_counts":{src:(hits if src==SOURCE[channel] else 0)
                      for src in SOURCES},
                "garfield_luts":provenance[0]["immutable_inputs"][
                    "garfield_map_provenance"],
                "angle_bins_deg":ANGLES.tolist(),
                "time_bin_width_ns":10,
                "geometries":6,
                "source_station":f"{SOURCE[channel]} MS{station}",
                "n_muon_gas_crossings_station":hits,
                "chunk_manifest_sha256":[sha_file(shard_path(
                   out,channel,i,station).with_suffix(".manifest.json"))
                   for i in range(1,number+1)],
                "statistical_errors":"PARENT-EVENT BOOTSTRAP NOT YET COMPUTED; "
                    "single-electron uncertainties MUST NOT be used instead",
            }
            payload={"meta":meta,"group_count":len(combined),
                     "groups":[dict(as_row(k,b),**{}) for k,b
                               in sorted(combined.items())]}
            path=out/"merged"/channel.lower()/f"MS{station}.json"
            path.parent.mkdir(parents=True,exist_ok=True)
            temp=path.with_suffix(".json.incomplete")
            temp.write_text(json.dumps(payload,ensure_ascii=False,
                          separators=(",",":"),allow_nan=False),
                            encoding="utf-8")
            os.replace(temp,path)
            dest=portal/SOURCE[channel]/f"MS{station}_layout_v3.html"
            render_v3(path,dest,SOURCE[channel],station)
            render_detail(path,portal/SOURCE[channel]/f"MS{station}.html")
            all_refs.append((SOURCE[channel],station))
            print(f"PASS merged source {SOURCE[channel]} MS{station}: "
                  f"{hits:,} real native muon gas hits, "
                  f"{len(combined):,} groups",flush=True)
    launcher=make_portal(portal,all_refs)
    print("PASS: complete 100k-generated-parent/channel first two ROOT "
          "chunks refolded through six frozen Garfield LUTs.")
    print("New independent high-stat dashboard:",launcher)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("mode",choices=("smoke","fold","merge"))
    p.add_argument("--channel",choices=CHANNELS+("all",),default="Jpsi")
    p.add_argument("--station",type=int,default=0,
                   help="-1 = all 6 MS; 0..5 = one station")
    p.add_argument("--chunks",type=int,default=1,
                   help="Count of consecutive completed 50k parent chunks")
    p.add_argument("--max-muon-hits-per-station",type=int,default=0,
                   help="smoke only; e.g. 200 true muon gas hits/MS")
    args=p.parse_args()
    if args.station not in (-1,0,1,2,3,4,5):
        p.error("station must be -1 (all) or MS0..MS5")
    if not 1<=args.chunks<=40:
        p.error("chunks must be in 1..40")
    if args.mode=="smoke":
        if not 1<=args.max_muon_hits_per_station<=1000:
            p.error("smoke needs --max-muon-hits-per-station in [1,1000]")
        if args.chunks!=1:
            p.error("smoke permits exactly one 50k parent input chunk")
    elif args.max_muon_hits_per_station:
        p.error("Only smoke mode permits a muon gas segment limit")
    if args.mode=="merge":
        merge(OUTPUT/"full",0,args.chunks)
        return
    limit=args.max_muon_hits_per_station
    out=OUTPUT/("smoke" if limit else "full")
    pools,provenance=config_pools()
    lam,marks=read_heed(HEED_TRACKS,HEED_CLUSTERS)
    for ch in (CHANNELS if args.channel=="all" else (args.channel,)):
        for i in range(1,args.chunks+1):
            for st in (range(6) if args.station==-1 else (args.station,)):
                fold_one(ch,i,st,limit,pools,provenance,lam,marks,out)


if __name__=="__main__":
    main()
