#!/usr/bin/env python3
"""Measured 6M-parent run stage times + paired Geom/binomial vs six LUT clock.

No Geant4, Heed executable or Garfield++ is invoked by this script.
Paired method benchmark invokes the already verified Py-Heed sampler and
six FROZEN Gaussian/Garfield response maps for the SAME 2000 real native
gas crossings per source/station in late chunk0040. Exact wall times are
measured with perf_counter_ns in pi_geometric_vs_lut_benchmark.py.

Definitions:
- Geant4 elapsed: sum completed.json native simulator wall seconds (each
  chunk's program internally uses multiple worker processes).
- ROOT native export elapsed: sum original CSV manifest measured times.
- Six-Garfield-Heed fold elapsed: sum unique source/chunk/station shard
  measured wall seconds INCLUDING shared Heed, response, grouping and gzip.
- Geometry vs frozen-LUT time: method-exclusive Python wall seconds, same
  source/chunk/station/muon gas crossings, same Heed sampled births.
- Actual microscopic Garfield LUT generation time: original per-column
  physical runs, one-time, if original independently verified timing exists.
Summed stage wall seconds are sequential stage sum, not exact workstation
makespan when e.g. parallel worker wall time or external I/O differs.
"""
from __future__ import annotations
import argparse
import csv
import json
import os
import subprocess
import sys
from pathlib import Path
from statistics import median

from pi_2m_export_completed_chunks import SUITE,DEST as EXPORT,CHANNELS
from pi_2m_native_six_geometry_fold import OUTPUT
from pi_geometric_vs_lut_benchmark import CHANNELS as BENCH_CHANNELS
from pi_six_geometry_native_timing import PITCHES
from pi_pitch_fieldwire_pilot import tag
from run_pi_geant4_campaign import REPO

HERE=Path(__file__).resolve().parent
OUT=REPO/"test_runs/mwpc_stage_c/pi_weekend_2m_operations_v1"
BENCH=OUT/"paired_runtime_benchmark_v1"
N_CHUNKS=40
SAMPLE=2000


def path_benchmark(ch,st):
    return BENCH/f"{ch.lower()}_chunk0040_MS{st}_sample{SAMPLE}.json"


def benchmark():
    BENCH.mkdir(parents=True,exist_ok=True)
    for ch in CHANNELS:
        for st in range(6):
            outfile=path_benchmark(ch,st)
            if outfile.exists():
                saved=json.loads(outfile.read_text(encoding="utf-8"))
                obs=saved["real_input"]
                if (obs["chunk"]!=40 or obs["station"]!=st
                    or obs["n_sampled_muon_segments"]!=SAMPLE
                    or obs["source"]!=dict(Jpsi="jpsi40-bore",Phi="phi40-bore",
                                            Omega="omega40-bore")[ch]
                    or len(saved["comparison"])!=6):
                    raise ValueError(f"Invalid completed paired runtime benchmark: {outfile}")
                print(f"SKIP verified matched runtime: {ch} MS{st}",flush=True)
                continue
            cmd=[sys.executable,str(HERE/"pi_geometric_vs_lut_benchmark.py"),
                 "--channel",ch,"--chunk","40","--station",str(st),
                 "--sample-muon-hits",str(SAMPLE),
                 "--out",str(outfile)]
            print(f"RUN matched Poisson/binomial vs Garfield LUT clocks: {ch} MS{st}",flush=True)
            subprocess.run(cmd,check=True,cwd=REPO)
    print("PASS: 18 matched 2,000-muon-gas-hit benchmark physics+method clocks")


def report():
    rows=[]
    counts=dict(geant4=0,export=0,lut_fold=0)
    clocks=dict(geant4=0.,geant4_validation=0.,root_export=0.,six_lut_fold=0.)
    for ch in CHANNELS:
        for i in range(1,N_CHUNKS+1):
            cp=SUITE/ch.lower()/f"chunk{i:04d}"/"completed.json"
            if cp.exists():
                item=json.loads(cp.read_text(encoding="utf-8"))
                if item["channel"]!=ch or item["chunk"]!=i or item["events"]!=50000:
                    raise ValueError(f"Invalid Geant4 elapsed provenance: {cp}")
                clocks["geant4"]+=float(item["wall_seconds"]["geant4"])
                clocks["geant4_validation"]+=float(item["wall_seconds"]["root_validation"])
                counts["geant4"]+=1
            exp=EXPORT/ch.lower()/f"chunk{i:04d}"/"crossings_manifest.json"
            if exp.exists():
                item=json.loads(exp.read_text(encoding="utf-8"))
                if item["immutable_input"]["channel"]!=ch or item["immutable_input"]["chunk"]!=i:
                    raise ValueError(f"Invalid ROOT export provenance: {exp}")
                clocks["root_export"]+=float(item["measured_export_wall_seconds"])
                counts["export"]+=1
            for st in range(6):
                fold=OUTPUT/"full"/ch.lower()/f"MS{st}"/f"chunk{i:04d}.manifest.json"
                if not fold.exists():continue
                item=json.loads(fold.read_text(encoding="utf-8"))
                f=item["immutable_inputs"]
                if f["chunk"]!=i or f["station"]!=st:
                    raise ValueError(f"Wrong response shard: {fold}")
                clocks["six_lut_fold"]+=float(item["measured_fold_seconds"])
                counts["lut_fold"]+=1
    paired=[]
    for ch in CHANNELS:
        for st in range(6):
            file=path_benchmark(ch,st)
            if file.exists():
                doc=json.loads(file.read_text(encoding="utf-8"))
                if len(doc["comparison"])!=6 or doc["real_input"]["chunk"]!=40:
                    raise ValueError(f"Wrong paired reference: {file}")
                paired.append(doc)
    configs=[tag(p,d) for p in (3.,3.5,4.) for d in (0.,75.)]
    paired_metrics=[]
    for cfg in configs:
        geometry=0.;lut=0.;heed=0.;n=0;deltas=[]
        phys=[]
        for doc in paired:
            datum=next((x for x in doc["comparison"]
                        if x["configuration"]==cfg),None)
            if datum is None:raise ValueError(f"Config missing in paired benchmark: {cfg}")
            geometry+=datum["geometry_method_seconds"]
            lut+=datum["lut_method_seconds"]
            heed+=doc["shared_heed"]["heedsampling_seconds"]
            n+=datum["n_muon_gas_segments_sampled"]
            deltas.append(datum["lut_over_geometry_processing_ratio"])
            value=datum.get("physical_garfield_lut_creation_measured_seconds")
            if value is not None:phys.append(value)
        paired_metrics.append({
            "config":cfg,"n_matched_muon_gas_crossings":n,
            "pure_geometry_seconds":geometry,
            "frozen_lut_lookup_seconds":lut,
            "LUT_over_geometry_ratio":lut/geometry if geometry else None,
            "median_per_station_LUT_over_geometry_ratio":median(deltas) if deltas else None,
            "shared_heed_seconds_across_station_samples":heed,
            "geometry_plus_shared_heed_seconds":geometry+heed,
            "lut_plus_shared_heed_seconds":lut+heed,
            "original_microscopic_Garfield_LUT_creation_seconds":(
                phys[0] if len(phys)==len(paired) and phys else None),
        })
    summary={
        "version":"PI_2M_PER_SOURCE_CLOCK_REPORT_V1",
        "generated_parent_events_per_source_target":2000000,
        "parent_events_total_target":6000000,
        "simulation_stage_completed_chunks":counts,
        "simulation_stage_expected_chunks":{"geant4":120,"export":120,"lut_fold":720},
        "measured_stage_elapsed_seconds":clocks,
        "matched_first_2000_actual_muon_gas_crossings_per_source_station_chunk40":
            len(paired),
        "matched_expected_18_source_station_samples":18,
        "paired_clock_comparison":paired_metrics,
        "interpretation":{
            "stage_wall_clock":"sum of measured completed chunk elapsed times; excludes idle, user pauses",
            "paired_comparison":"same native muon gas hits, same original Heed cluster births, same primary electron random quantiles",
            "geometry":"nearest finite-anode assignment of WHOLE Poisson-Heed clusters using straight path",
            "lut":"actual six precomputed Garfield++ microscopic electron transport lookup",
            "microscopic":"map production once, not repeated for every Geant4 electron; includes per-column startup",
            "statistical":"clock is hardware/run dependent, not a physics observable; no Monte Carlo uncertainty inferred",
        }
    }
    OUT.mkdir(parents=True,exist_ok=True)
    tmp=OUT/"runtime_comparison.json.incomplete"
    tmp.write_text(json.dumps(summary,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    os.replace(tmp,OUT/"runtime_comparison.json")
    lines=[
        "NA60+ six-geometry RUN TIME COMPARISON (actual measured clocks)",
        "Generated parents: target 2,000,000 per source (6,000,000 total).",
        f"Native Geant4 chunks: {counts['geant4']}/120; ROOT exports: {counts['export']}/120;",
        f"Six-LUT source/chunk/station fold shards: {counts['lut_fold']}/720;",
        f"Paired 2k-hit Jpsi/Phi/Omega x MS0..MS5 on chunk0040: {len(paired)}/18.",
        "Completed-chunk sequential elapsed-time SUMS (not wall-clock makespan):",
    ]
    for key,value in clocks.items():
        lines.append(f"  {key:24s} {value:10.1f} s | {value/3600:8.3f} h")
    lines.append("Paired method-exclusive physics postprocess times, same Heed births:")
    for x in paired_metrics:
        if x["n_matched_muon_gas_crossings"]:
            lines.append(f" {x['config']:24s} {x['n_matched_muon_gas_crossings']:>6,d}"
                         f" gas hits; geometric {x['pure_geometry_seconds']:6.3f}s;"
                         f" LUT {x['frozen_lut_lookup_seconds']:6.3f}s;"
                         f" LUT/geom {x['LUT_over_geometry_ratio']:.2f}x;"
                         f" old microscopic Garfield map build "
                         f"{x['original_microscopic_Garfield_LUT_creation_seconds']}s")
    lines+=["","IMPORTANT: Geant4 and Heed are COMMON for both methods.",
            "The geometry model does not predict time-of-drift, physical "
            "electron loss, Garfield field propagation or electronics.",
            "The time ratio refers to two CPU postprocessing methods, "
            "not an end-to-end electronics/physics equivalence."]
    text="\n".join(lines)+"\n"
    (OUT/"runtime_comparison.txt").write_text(text,encoding="utf-8")
    print(text,flush=True)
    return summary


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("mode",choices=("benchmark","report"))
    args=p.parse_args()
    if args.mode=="benchmark":benchmark()
    report()


if __name__=="__main__":
    main()
