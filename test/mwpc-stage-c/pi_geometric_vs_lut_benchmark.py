#!/usr/bin/env python3
"""Matched real-native-muon geometry/Poisson-Heed vs Garfield LUT experiment.

Same source+chunk+event Geant4 gas hits and SAME Heed cluster/multiplicity
realization, all six pitches/field wires. Cluster birth points follow the
straight native gas segment uniformly in length. The geometric model
assigns EACH WHOLE CLUSTER to nearest *finite* physical anode (binomial
conditional on Nclusters and segment-length fraction for a two-cell path).
No detector charge loss or analog discriminator in that baseline.

LUT model propagates the same primary electrons through the SIX frozen
microscopic Garfield response maps (per-electron collection, wire shift,
time); finite wire loss and transport across midpoint are permitted.

Timing uses perf_counter_ns around ONLY method computation, shared Heed
sampling separately, CSV scan and LUT loads separately. This is a speed
comparison of GEOMETRIC POSTPROCESSING vs LOOKUP POSTPROCESSING, NOT
a timing measurement of microscopic Garfield++ LUT GENERATION.
NO Geant4 or Garfield++ executable runs here.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
import os
import time
from collections import defaultdict
from pathlib import Path

import numpy as np

from pi_six_geometry_native_timing import (
    MAPS,HEED_TRACKS,HEED_CLUSTERS,PITCHES,
    birth_sample,load_map,tag,n_anodes,response_for_births,stable_seed,
)
from pi_2m_export_completed_chunks import (
    DEST as EXPORT,SOURCE,PARENTS_PER_CHUNK,SUITE,
)
from pi_pitch_fieldwire_pilot import load_column_timing
from pi_2m_native_six_geometry_fold import crossing_provenance,sha_file
from run_pi_per_chamber_response import read_heed
from run_pi_geant4_campaign import REPO

OUT=REPO/"test_runs/mwpc_stage_c/pi_geometric_vs_lut_benchmark_v1"
VERSION="PI_NATIVE_MATCHED_GEOMETRY_LUT_BENCH_V1"
THRESHOLDS=(1,2,3,5)
CHANNELS=("Jpsi","Phi","Omega")


def shared_heed_births(row,lam,marks,seed):
    """Mirror original RNG sequence; additionally retain whole clusters.

    Each cluster i has exactly one path coordinate and mark m_i>=1.
    Quantiles q are one per electron, as in birth_sample; checking exact
    identity against the original sampler is done in the synthetic test.
    """
    rng=np.random.default_rng(seed)
    length=float(row["path_length_mm"])
    if not 0<length<1000 or not math.isfinite(length):
        raise ValueError("Invalid Geant4 gas path")
    uvw=np.array([[float(row[f"{axis}_in_mm"]) for axis in "uvw"],
                  [float(row[f"{axis}_out_mm"]) for axis in "uvw"]],
                 dtype=np.float64)
    if not np.isfinite(uvw).all() or abs(uvw[:,1]).max()>3.05:
        raise ValueError("Invalid TGeo Garfield chamber-local position")
    ncluster=int(rng.poisson(lam*length))
    marks_drawn=marks[rng.integers(0,len(marks),size=ncluster)] if ncluster \
        else np.empty(0,dtype=np.int64)
    if np.any(marks_drawn<1):
        raise AssertionError("Empirical Heed positive cluster mark required")
    tcluster=rng.random(ncluster)
    telectron=np.repeat(tcluster,marks_drawn)
    positions=(uvw[0]+telectron[:,None]*(uvw[1]-uvw[0])) if \
        len(telectron) else np.empty((0,3),dtype=np.float64)
    q=rng.random(len(telectron))
    return positions,q,tcluster,marks_drawn,length,uvw


def geometric_wire_indices(u0,u1,tcluster,pitch,station):
    """Nearest existing physical anode at each shared cluster birthplace.

    Midpoint split has exact crossing fractions along the *3D segment*
    because straight-path arc length is proportional to t in [0,1].
    """
    n=n_anodes(station,pitch)
    centers=(np.arange(n,dtype=np.float64)-(n-1)/2)*pitch
    boundaries=(centers[1:]+centers[:-1])/2
    ucluster=u0+(u1-u0)*tcluster
    wire=np.searchsorted(boundaries,ucluster,side="right")
    # For near-normal u1=u0, unique full path in one Voronoi cell.
    if u1==u0:
        lengths=np.array([1.],dtype=np.float64)
    else:
        tcut=(boundaries-u0)/(u1-u0)
        cut=np.sort(tcut[(tcut>0)&(tcut<1)])
        if len(cut):
            cut=np.unique(cut)
        lengths=np.diff(np.concatenate(([0.],cut,[1.])))
    if not math.isclose(float(lengths.sum()),1.0,abs_tol=1.e-11):
        raise AssertionError("Unphysical physical-wire path partition")
    return wire,lengths


def two_cell_probabilities(lengths,ncluster,mu):
    """Return conditional binomial and unconditional Poisson-thinning P2.

    Only defined for a track intersecting EXACTLY two wire Voronoi cells.
    For given N, K ~ Binomial(N,p); mark draws do NOT split cluster.
    """
    if len(lengths)!=2:
        return None
    p=float(lengths[0])
    conditional=(1-p**ncluster-(1-p)**ncluster) if ncluster else 0.
    unconditional=(-math.expm1(-mu*p))*(-math.expm1(-mu*(1-p)))
    return conditional,unconditional,min(p,1-p)


def select_muon_segments(csvpath,station,nwanted,seed):
    """Deterministic uniform reservoir sample over all selected gas hits.

    Sampling unit is gas crossing, NOT generated parent event. Report as
    method-matched physics benchmark, not parent-statistics confidence.
    nwanted=0 loads EVERY selected station gas segment (can need lots RAM).
    """
    rng=np.random.default_rng(seed)
    rows=[]
    n_seen=0
    t0=time.perf_counter_ns()
    with csvpath.open(newline="",encoding="utf-8") as f:
        reader=csv.DictReader(f)
        necessary={"event","hit_index","track_id","pdg","station",
                   "path_length_mm","u_in_mm","u_out_mm",
                   "v_in_mm","v_out_mm","w_in_mm","w_out_mm"}
        if not necessary.issubset(reader.fieldnames or []):
            raise ValueError("Native Geant4 CSV missing expected gas hit columns")
        for row in reader:
            if int(row["station"])!=station or abs(int(row["pdg"]))!=13:
                continue
            event=int(row["event"])
            if not 0<=event<PARENTS_PER_CHUNK:
                raise ValueError("Local Geant4 event outside native chunk")
            n_seen+=1
            if not nwanted or len(rows)<nwanted:
                rows.append(row)
            else:
                k=int(rng.integers(0,n_seen))
                if k<nwanted: rows[k]=row
    if n_seen==0:
        raise ValueError("No real muon gas hits in this station")
    return rows,n_seen,(time.perf_counter_ns()-t0)*1.e-9


def compute(rows,channel,chunk,station,lam,marks,pools):
    measures={cfg:defaultdict(int) for cfg in pools}
    timings={cfg:{"geometry_ns":0,"lut_ns":0} for cfg in pools}
    shared_heed_ns=0
    two_cell={cfg:dict(n=0,p_cond_sum=0.,p_poisson_sum=0.,
                       p_minor_sum=0.,geo_two=0,lut_two=0) for cfg in pools}
    stats={cfg:{"clusters":0,"electrons":0,"collected":0,
                "migrated":0,"n":0} for cfg in pools}
    for row in rows:
        seed=stable_seed("native-six-pitch-shared-Heed-v1",SOURCE[channel],
                         f"chunk{chunk:04d}",row["event"],row["hit_index"],
                         row["track_id"])
        t0=time.perf_counter_ns()
        positions,q,tcluster,mult,length,uvwin=shared_heed_births(
            row,lam,marks,seed)
        shared_heed_ns+=time.perf_counter_ns()-t0
        u=positions[:,0];v=positions[:,1]
        ne=len(q)
        for cfg,(pitch,fieldwire,pool) in pools.items():
            t0=time.perf_counter_ns()
            geomcluster,segments=geometric_wire_indices(
                uvwin[0,0],uvwin[1,0],tcluster,pitch,station)
            geomelectron=np.repeat(geomcluster,mult)
            nwire=n_anodes(station,pitch)
            geometrichits=np.bincount(geomelectron,minlength=nwire)
            gactive={threshold:int(np.count_nonzero(
                geometrichits>=threshold)) for threshold in THRESHOLDS}
            timings[cfg]["geometry_ns"]+=time.perf_counter_ns()-t0
            t0=time.perf_counter_ns()
            finite,rawcol,wire,samples=response_for_births(
                u,v,q,pitch,station,pool)
            collected=wire[finite]
            luthits=np.bincount(collected,minlength=nwire)
            lactive={threshold:int(np.count_nonzero(
                luthits>=threshold)) for threshold in THRESHOLDS}
            migrated=int(np.count_nonzero(
                collected!=geomelectron[finite]))
            timings[cfg]["lut_ns"]+=time.perf_counter_ns()-t0
            m=measures[cfg]
            m["gas_hits"]+=1
            m["clusters"]+=len(mult)
            m["primary_electrons"]+=ne
            m["lut_collected_electrons"]+=len(collected)
            m["lut_lost_electrons"]+=ne-len(collected)
            m["lut_migrated_collected_electrons"]+=migrated
            for threshold in THRESHOLDS:
                m[f"geo_multiwire_atleast_{threshold}e"]+=int(
                    gactive[threshold]>=2)
                m[f"lut_multiwire_atleast_{threshold}e"]+=int(
                    lactive[threshold]>=2)
            predicted=two_cell_probabilities(segments,len(mult),lam*length)
            if predicted is not None:
                c=two_cell[cfg]
                c["n"]+=1
                c["p_cond_sum"]+=predicted[0]
                c["p_poisson_sum"]+=predicted[1]
                c["p_minor_sum"]+=predicted[2]
                c["geo_two"]+=int(gactive[1]>=2)
                c["lut_two"]+=int(lactive[1]>=2)
    return measures,timings,shared_heed_ns,two_cell


def original_microscopic_garfield_times(pools):
    """Recover ACTUAL old microscopic Garfield runtime, if fully recorded.

    IMPORTANT: this is upfront six-LUT map-generation time, measured when
    the physical Garfield electrons were originally drifted. The 6×
    column timing includes per-column process startup (not lookup costs).
    If any sidecar is missing, report unknown rather than fabricate.
    """
    output={}
    for cfg,(pitch,fieldwire,pool) in pools.items():
        directory=MAPS/cfg
        config=json.loads((directory/"frozen_config.json").read_text())
        ncols=int(config["n_u"])
        percol=[]
        all_recorded=True
        for index in range(ncols):
            filename=directory/"columns"/f"iu_{index:03d}.csv"
            if not filename.is_file():
                all_recorded=False
                break
            elapsed=load_column_timing(filename)
            if elapsed is None:
                all_recorded=False
                break
            percol.append(elapsed)
        output[cfg]={
            "physical_garfield_lut_creation_measured_seconds":(
                sum(percol) if all_recorded and len(percol)==ncols else None),
            "columns_measured":len(percol),
            "columns_total":ncols,
            "n_garfield_microscopic_electron_traces":(
                int(config["n_u"])*int(config["n_v"])*
                int(config["samples_per_box"])),
            "measurement_includes_repeated_column_process_startup":True,
            "measurement_source":"Garfield original per-column timing.json "
                                  "checkpoints, CRC-validated against CSV",
        }
    return output


def measured_native_geant4_chunk(channel,chunk):
    path=SUITE/channel.lower()/f"chunk{chunk:04d}"/"completed.json"
    if not path.is_file():
        return None
    doc=json.loads(path.read_text())
    if (doc.get("events")!=50000 or doc.get("channel")!=channel
        or doc.get("chunk")!=chunk):
        raise ValueError(f"Not a matching completed Geant4 chunk: {path}")
    return {
        "n_generated_parent_events":doc["events"],
        "geant4_native_transport_seconds":doc.get("wall_seconds",{}).get("geant4"),
        "native_root_validation_seconds":doc.get("wall_seconds",{}).get("root_validation"),
        "checkpoint":"completed.json with measured monotonic elapsed times",
    }


def summarize(rows,seen,read_sec,channel,chunk,station,lam,marks,
              pools,load_sec,out):
    total_start=time.perf_counter_ns()
    counts,times,common_ns,two=compute(rows,channel,chunk,station,
                                      lam,marks,pools)
    compute_s=(time.perf_counter_ns()-total_start)*1.e-9
    summary=[]
    for cfg in pools:
        c=counts[cfg];n=c["gas_hits"]
        t=times[cfg]
        model=two[cfg]
        outrow={
            "configuration":cfg,
            "pitch_mm":pools[cfg][0],
            "field_wire_diameter_um":pools[cfg][1],
            "n_muon_gas_segments_sampled":n,
            "n_heed_clusters":c["clusters"],
            "n_heed_primary_electrons":c["primary_electrons"],
            "lut_collection_fraction":(
                c["lut_collected_electrons"]/c["primary_electrons"]
                if c["primary_electrons"] else None),
            "lut_migration_fraction_among_collected":(
                c["lut_migrated_collected_electrons"]/
                    c["lut_collected_electrons"]
                if c["lut_collected_electrons"] else None),
            "geometry_method_seconds":t["geometry_ns"]*1.e-9,
            "lut_method_seconds":t["lut_ns"]*1.e-9,
            "lut_over_geometry_processing_ratio":(
                t["lut_ns"]/t["geometry_ns"]
                if t["geometry_ns"]>0 else None),
            "geometry_plus_shared_heed_seconds":(
                (common_ns+t["geometry_ns"])*1.e-9),
            "lut_plus_shared_heed_seconds":(
                (common_ns+t["lut_ns"])*1.e-9),
            "two_geometric_cell_gas_segments":model["n"],
            "two_cell_expected_binomial_conditional_fraction":(
                model["p_cond_sum"]/model["n"] if model["n"] else None),
            "two_cell_expected_poisson_unconditional_fraction":(
                model["p_poisson_sum"]/model["n"] if model["n"] else None),
            "two_cell_observed_geo_multiwire_fraction":(
                model["geo_two"]/model["n"] if model["n"] else None),
            "two_cell_observed_lut_multiwire_fraction":(
                model["lut_two"]/model["n"] if model["n"] else None),
        }
        for threshold in THRESHOLDS:
            outrow[f"P_geo_at_least_two_wires_ge_{threshold}_electrons"]=(
                c[f"geo_multiwire_atleast_{threshold}e"]/n)
            outrow[f"P_lut_at_least_two_wires_ge_{threshold}_electrons"]=(
                c[f"lut_multiwire_atleast_{threshold}e"]/n)
        summary.append(outrow)
    microscopic=original_microscopic_garfield_times(pools)
    for item in summary:
        item.update(microscopic[item["configuration"]])
    result={
        "version":VERSION,
        "real_input":{
            "source":SOURCE[channel],"chunk":chunk,"station":station,
            "generated_parent_events_in_source_chunk":50000,
            "sampling_unit":"MUON GAS SEGMENT, not generated parent event",
            "n_selected_station_muon_segments_available":seen,
            "n_sampled_muon_segments":len(rows),
            "selection":"uniform seeded reservoir over station muon gas segments",
            "sampling_seed":20261023,
            "source_chunk_local_event_key":"source+chunk+event",
        },
        "shared_heed":{
            "cluster_density_per_mm":lam,
            "heed_marks_sha256":sha_file(HEED_CLUSTERS),
            "heed_tracks_sha256":sha_file(HEED_TRACKS),
            "heedsampling_seconds":common_ns*1.e-9,
            "interpretation":"ONE identical Poisson cluster/empirical marks "
                  "birth population for the geometry and LUT methods. "
                  "All members of a cluster belong to one geometric wire.",
        },
        "timing":{
            "clock":"time.perf_counter_ns, sequential Python on this WSL host",
            "native_upstream_geant4_50k_generated_parents":
                measured_native_geant4_chunk(channel,chunk),
            "microscopic_lut_generation_source":
                "Stored timing.json sidecars from original Garfield physics "
                "runs, separately timed from this postprocessing benchmark",
            "source_csv_read_and_sampling_seconds":read_sec,
            "six_garfield_lut_load_seconds":load_sec,
            "matched_six_methods_compute_wall_seconds":compute_s,
            "per_config_scoped_model_times_exclude_Heed_and_disk_io":True,
            "direct_microscopic_Garfield_simulation_per_electron":"NOT run",
            "historic_Garfield_LUT_creation_wall_time":
                "Read ONLY from original validated Garfield per-column "
                "timing sidecars when all columns present; no new traces.",
            "native_Geant4_generation_wall_time":"separately measured",
            "note":"Two models are postprocessing approaches, not "
                   "independent full Geant4+ionization+Garfield simulations. "
                   "Shared Heed and lookup setup costs are identified.",
        },
        "comparison":summary,
    }
    out.parent.mkdir(parents=True,exist_ok=True)
    if out.exists():
        raise ValueError(f"Refusing to overwrite prior physics benchmark: {out}")
    temp=out.with_suffix(out.suffix+".incomplete")
    temp.write_text(json.dumps(result,indent=2,ensure_ascii=False,
                   allow_nan=False)+"\n",encoding="utf-8")
    os.replace(temp,out)
    print("BENCHMARK: matched native Heed geometry/binomial vs frozen Garfield LUT")
    print(f"Source: {SOURCE[channel]}; chunk{chunk:04d}; MS{station}. "
          f"Random sample: {len(rows):,}/{seen:,} actual MUON GAS crossings")
    print(f"Measured shared Heed sampling: {common_ns*1e-9:.3f} s")
    print(f"Measured input read/reservoir: {read_sec:.3f} s")
    print(f"Measured six LUT loads: {load_sec:.3f} s")
    print(f"Measured six-pitch model loop including comparisons: {compute_s:.3f} s")
    print("Config                       Geometry(s)   LUT(s)   LUT/geometry"
          "    P>=2 geom   P>=2 LUT")
    for x in summary:
        print(f"{x['configuration']:<27}"
              f" {x['geometry_method_seconds']:>10.3f}"
              f" {x['lut_method_seconds']:>8.3f}"
              f" {x['lut_over_geometry_processing_ratio']:>11.2f}"
              f" {x['P_geo_at_least_two_wires_ge_1_electrons']:>11.4f}"
              f" {x['P_lut_at_least_two_wires_ge_1_electrons']:>11.4f}")
    print("Original full microscopic Garfield LUT creation wall times "
          "(SEPARATE, per 6 configurations, if all timing sidecars exist):")
    for cfg in pools:
        value=microscopic[cfg]["physical_garfield_lut_creation_measured_seconds"]
        print(f"  {cfg}: " + (f"{value:.1f} s measured" if value is not None
            else "NOT RECORDED for all grid columns"))
    print("PASS: matched first-order geometry vs six frozen Garfield LUT "
          "time and multiwire comparison")
    print("Frozen result:",out)
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--channel",choices=CHANNELS,default="Jpsi")
    p.add_argument("--chunk",type=int,default=1)
    p.add_argument("--station",type=int,default=0)
    p.add_argument("--sample-muon-hits",type=int,default=2000)
    p.add_argument("--out",type=Path,default=None)
    args=p.parse_args()
    if not 1<=args.chunk<=40 or not 0<=args.station<6 or (
        not 1<=args.sample_muon_hits<=20000):
        p.error("chunk 1..40, station 0..5, sample hits 1..20k required")
    csvpath,manifest=crossing_provenance(args.channel,args.chunk)
    time_start=time.perf_counter_ns()
    pools={tag(p,d):(p,d,load_map(p,d,MAPS)[0]) for p in PITCHES
           for d in (0.,75.)}
    load_sec=(time.perf_counter_ns()-time_start)*1.e-9
    lam,marks=read_heed(HEED_TRACKS,HEED_CLUSTERS)
    rows,seen,read_sec=select_muon_segments(
        csvpath,args.station,args.sample_muon_hits,20261023)
    path=args.out or OUT/(
        f"{args.channel.lower()}_chunk{args.chunk:04d}_MS{args.station}_"
        f"sample{len(rows)}.json")
    summarize(rows,seen,read_sec,args.channel,args.chunk,args.station,
              lam,marks,pools,load_sec,path)


if __name__=="__main__":
    main()
