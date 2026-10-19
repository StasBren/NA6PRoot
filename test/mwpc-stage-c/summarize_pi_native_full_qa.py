#!/usr/bin/env python3
"""Summarize and CONSERVATION-CHECK the FINISHED native full Geant4 PI QA.

Standard-library only. Does not read or rewrite any Geant4 ROOT files;
uses the four already generated full-QA products. In particular,
per-generated-event GAS HIT rate is not the probability per incoming muon.
"""
from __future__ import annotations
import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path

from mwpc_checkerboard_projection import SOURCES, build_chambers

HERE=Path(__file__).resolve().parent
REPO=HERE.parents[1]
DEFAULT=REPO/"test_runs/mwpc_stage_c/pi_native_full_quality"
GEOMETRY={c.det_id:c for group in build_chambers() for c in group}
SOURCE_LABEL={"jpsi40-bore":"J/psi","phi40-bore":"phi","omega40-bore":"omega"}
VARIABLES=(
    "edep_eV","path_mm","theta_u_in_deg","theta_u_out_deg",
    "theta_w_in_deg","theta_w_out_deg","delta_theta_u_deg",
    "delta_theta_w_deg",
)


def load(path):
    if not path.is_file():
        raise FileNotFoundError(path)
    with path.open(newline="",encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def number(value):
    return float(value) if value not in ("",None) else math.nan


def require(cond,message):
    if not cond:
        raise ValueError("PI_QA_CONSERVATION_FAIL: "+message)


def checks_and_readout(root):
    report=json.loads((root/"pi_native_all_QA_report.json").read_text(
        encoding="utf-8"))
    require(report["full_campaign_analyzed"] is True,
            "Quality report covers a subset, not full 201 chunks")
    require(int(report["processed_chunks_per_source"])==
            int(report["total_chunks_per_source"])==67,
            "Wrong completed chunk total")
    events=int(report["generated_events_per_source"])
    require(events==200000,"Expected 200,000 generated events per source")

    chambers=load(root/"pi_native_all_chamber_quality.csv")
    stations=load(root/"pi_native_all_station_quality.csv")
    histograms=load(root/"pi_native_all_histograms.csv")
    require(len(chambers)==3*298,"Expected exactly 894 module rows")
    require(len(stations)==3*6,"Expected exactly 18 station rows")

    mod={}
    stmap={}
    for r in chambers:
        src,det=r["source"],int(r["det_id"])
        key=src,det
        require(src in SOURCES and det in GEOMETRY,
                "Unexpected chamber "+str(key))
        require(key not in mod,"Duplicate chamber "+str(key))
        require(int(r["station"])==GEOMETRY[det].station,
                "Wrong station assignment for DET"+str(det))
        require(int(r["gas_hit_rows"])>=0,"Negative gas-hit count")
        mod[key]=r
    for r in stations:
        src,st=r["source"],int(r["station"])
        key=src,st
        require(src in SOURCES and 0<=st<6 and key not in stmap,
                "Bad/duplicate station row "+str(key))
        stmap[key]=r
    for source in SOURCES:
        for det in GEOMETRY:
            require((source,det) in mod,"Missing chamber "+str((source,det)))
        for st in range(6):
            require((source,st) in stmap,"Missing station "+str((source,st)))

    # Histograms are conditional per NATIVE gas-hit SEGMENT; zero-exposure
    # chambers contain bins with zero hit count, never NaN probabilities.
    hist=defaultdict(lambda:{"count":0,"n":None,"seen":set(),"first":None,
                              "last":None,"bad":0,"far":0})
    for r in histograms:
        src,det,var=r["source"],int(r["det_id"]),r["variable"]
        require((src,det) in mod and var in VARIABLES,
                "Unexpected histogram "+str((src,det,var)))
        lo,hi=float(r["bin_lower"]),float(r["bin_upper"])
        n=int(r["gas_hit_exposure"])
        k=int(r["hit_count"])
        p=number(r["conditional_bin_fraction"])
        expect=int(mod[src,det]["gas_hit_rows"])
        require(hi>lo and k>=0 and n==expect,
                "Bad histogram bounds/exposure for "+str((src,det,var)))
        if n:
            require(math.isfinite(p) and abs(p-k/n)<1.e-8,
                    "Wrong conditional histogram fraction "+str((src,det,var)))
        else:
            require(k==0 and math.isnan(p),
                    "Zero-exposure histogram cannot have a fraction")
        key=src,det,var
        acc=hist[key]
        bounds=lo,hi
        require(bounds not in acc["seen"],"Duplicate histogram bin "+str(key))
        acc["seen"].add(bounds)
        acc["count"]+=k
        require(acc["n"] is None or acc["n"]==n,
                "Inconsistent exposure "+str(key))
        acc["n"]=n
        if var.startswith("theta_") and (hi<=-10 or lo>=10):
            acc["far"]+=k
    require(len(hist)==3*298*len(VARIABLES),
            "Expected 8 complete histograms per detector")
    for (src,det,var),acc in hist.items():
        expect=int(mod[src,det]["gas_hit_rows"])
        require(acc["count"]==expect,
                "Histogram counts do not sum to gas hits: "+
                str((src,det,var,acc["count"],expect)))
        bounds=sorted(acc["seen"])
        require(all(math.isclose(bounds[i][1],bounds[i+1][0],
                 abs_tol=1.e-9,rel_tol=1.e-12)
                 for i in range(len(bounds)-1)),
                "Histogram has gaps/overlaps: "+str((src,det,var)))

    # True overlap is "same muon crosses multiple physical chambers";
    # repeated segments in one physical chamber are DIFFERENT. These
    # counters are not station incoming-muon denominators.
    station_checks=0
    for src in SOURCES:
        src_gas=0
        for st in range(6):
            s=stmap[src,st]
            n=int(s["gas_hit_rows"])
            sum_modules=sum(int(mod[src,c.det_id]["gas_hit_rows"])
                            for c in GEOMETRY.values() if c.station==st)
            require(sum_modules==n,
                    "Sum of detector gas hits != station gas hits for "+
                    str((src,st,sum_modules,n)))
            uniq=int(s["unique_gas_crossing_muon_tracks"])
            overlap=int(s["overlap_muon_tracks"])
            extra=int(s["extra_chamber_hits_per_overlaps"])
            repeated=int(s["repeated_gas_segments_in_same_module"])
            require(0<=overlap<=uniq and extra>=overlap and repeated>=0,
                    "Impossible station overlap counters "+str((src,st)))
            # Gas hit rows may exclude flagged zero-length hits. Only
            # enforce this inequality, not an equality with *all* raw hits.
            require(n<=uniq+extra+repeated,
                    "Detector gas hits exceed distinct+overlap+reentry "+
                    str((src,st)))
            station_checks+=1
            src_gas+=n
        metrics=report["per_source"][src]
        require(int(metrics["total_detector_hits"])==
                int(metrics["muon_hits"])+int(metrics["nonmuon_hits"]),
                "ALL particle counters inconsistent: "+src)
        require(src_gas==int(metrics["muon_hits"])-
                int(metrics["invalid_muon_hits"]),
                "Full-source gas-hit conservation failed: "+src)

    out=[]
    out.append("NA60+ / DI-CE MWPC — NATIVE Geant4 full-campaign PI readout")
    out.append("="*72)
    out.append("3 channels × 200,000 generated events, 67 chunks/channel.")
    out.append("FULL CAMPAIGN: true MCKine-tagged Geant4 gas hits, 298 modules.")
    out.append("QA conservation: PASS — 894 chambers, 18 stations, 7152")
    out.append("module histogram spectra; all exposure counts and 201 chunk")
    out.append("completion checks are consistent.")
    out.append("")
    for source in SOURCES:
        metrics=report["per_source"][source]
        label=SOURCE_LABEL[source]
        total=int(metrics["muon_hits"])
        other=int(metrics["nonmuon_hits"])
        invalid=int(metrics["invalid_muon_hits"])
        unique=int(metrics["unique_tracks_with_gas"])
        through=int(metrics["tracks_crossing_gas_in_all_six"])
        n_zero=sum(int(mod[source,d]["gas_hit_rows"])==0 for d in GEOMETRY)
        n_low=sum(0<int(mod[source,d]["gas_hit_rows"])<100 for d in GEOMETRY)
        out.append(f"SOURCE {label}")
        out.append(f"  Generated events     : {events:,}")
        out.append(f"  Native muon gas hits : {total:,}  ({total/events:.4f} hits/generated event)")
        out.append(f"  Other-particle hits  : {other:,}")
        out.append(f"  Invalid muon segments: {invalid:,}")
        out.append(f"  Unique muon IDs in gas: {unique:,}")
        out.append(f"  All-six-gas muon IDs : {through:,}")
        out.append(f"  Empty/low-exposure modules: {n_zero} / {n_low} (low means 1..99 hits)")
        out.append("  MS | Gas hits | Unique gas tracks | Overlap tracks | "
                   "Reentries | <Edep> keV | <path> mm")
        for st in range(6):
            row=stmap[source,st]
            out.append(
                f"  {st:>2} | {int(row['gas_hit_rows']):>8,} | "
                f"{int(row['unique_gas_crossing_muon_tracks']):>17,} | "
                f"{int(row['overlap_muon_tracks']):>14,} | "
                f"{int(row['repeated_gas_segments_in_same_module']):>9,} | "
                f"{number(row['mean_edep_eV'])/1000:>10.4f} | "
                f"{number(row['mean_path_mm']):>9.3f}"
            )
        out.append("")
    flagged=[]
    for source in SOURCES:
        for det,c in GEOMETRY.items():
            r=mod[source,det]
            n=int(r["gas_hit_rows"])
            if n==0:
                continue
            long=int(r["long_path_over_50mm"])
            angle=int(r["large_direction_change"])
            if long or angle:
                flagged.append((long+angle,source,det,c.station,n,long,angle))
    flagged.sort(reverse=True)
    out.append("CHAMBER QA FLAGS (not automatically detector problems)")
    out.append(f"  Modules with ≥1 recorded long-path/direction flag: "
               f"{len(flagged)} / 894 source×modules")
    for total,source,det,st,n,long,angle in flagged[:15]:
        out.append(f"  {SOURCE_LABEL[source]:>5} MS{st} DET{det:03d}: "
                   f"{n:,} gas hits, path>50mm: {long}, "
                   f"|delta theta|>5deg: {angle}")
    if not flagged:
        out.append("  No flags at the QA thresholds.")
    out.extend([
        "",
        "INTERPRETATION NOTES",
        "1. Gas hit counts per GENERATED event are NOT probabilities per incoming",
        "   station muon. Reference-plane denominators remain to be evaluated.",
        "2. Overlapping physical chambers can both record the SAME muon.",
        "3. Counts and Edep are Geant4 truth; no anode-gain/discriminator",
        "   or Heed-Garfield probability is included in these spectra.",
        "4. Histogram-count conservation does not, by itself, establish",
        "   expected Landau-like Edep shapes or all angular correlations.",
        "5. Regions with small n require large statistical uncertainties.",
    ])
    out_file=root/"pi_native_QA_readout.txt"
    out_file.write_text("\n".join(out)+"\n",encoding="utf-8")
    for line in out:
        print(line)
    print(f"\nPASS: full campaign QA readout: {out_file}")
    return out_file


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input",type=Path,default=DEFAULT)
    args=parser.parse_args()
    checks_and_readout(args.input)


if __name__=="__main__":
    main()
