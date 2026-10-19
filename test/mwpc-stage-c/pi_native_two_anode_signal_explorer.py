#!/usr/bin/env python3
"""Interactive TWO-ANODE native Geant4 muon event signals, centered 3+3mm.

Selects real FIRST-CHUNK native muon gas-hit segments with EXACTLY TWO
distinct finite physical anodes collecting primary electrons (N_all=2)
from the already frozen, validated centered Heed/Garfield timing pilot.
Replays EXACT SAME ionization/collection via that pilot's per-hit SHA256
seed; then digitizes each anode independently using explicit UNCALIBRATED
anode-only gain/current/shaping/noise surrogate. The illustrative sum of
the TWO INDEPENDENT channel outputs is NOT a physical summed readout.
No new ROOT, Geant4 or microscopic Garfield transport. No cathode strips.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path

import numpy as np

from mwpc_checkerboard_projection import build_chambers
from pi_centered_electron_sampler import (
    GEOMETRY_ID, ionize_and_drift, load_centered_map, multiplicities, stable_seed
)
from pi_anode_digitizer import AnodeConfig, digitize_hit, empirical_gains
from run_pi_per_chamber_response import read_heed
from pi_native_centered_timing_pilot import (
    BASE, EXPORT, CENTERED, OUT, SOURCES, GEOMETRY, CUTOFFS
)

HERE=Path(__file__).resolve().parent
VIEW=BASE/"pi_fixed_wire_maps/pi_two_anode_signal_explorer.html"
TEMPLATE=HERE/"pi_two_anode_signal_explorer_template.html"
INPUT=OUT/"pi_centered_native_hit_multiplicities.csv"
N_SAMPLES=7
# Display deliberately selected cases, NOT a randomly weighted
# detector-efficiency or pitch comparison sample.


def choose_examples(rows,max_examples):
    """Spread available signed muon angles for EACH (source,physical DET).

    Return at most max_examples per chamber-source, including both
    early two-anode collection and intentionally late second-anode
    timing outcomes if available. Sampling affects VIEW ONLY; raw
    selection/response pilot remains untouched.
    """
    bins=defaultdict(list)
    for row in rows:
        if int(row["N_all"])!=2:
            continue
        src=row["source"];det=int(row["det_id"])
        if src not in SOURCES or det not in GEOMETRY:
            raise ValueError("Unknown source/chamber in frozen native pilot")
        if int(row["station"])!=GEOMETRY[det].station:
            raise ValueError("Frozen pilot station mismatch")
        bins[src,det].append(row)
    result={}
    for key,vals in bins.items():
        vals.sort(key=lambda r:(
            float(r["theta_u_in_deg"]),
            int(r["event"]),int(r["hit_index"])
        ))
        chosen={}
        def take(r):
            chosen[(int(r["event"]),int(r["hit_index"]),int(r["track_id"]))]=r
        # Sample quantiles of signed angle (not absolute); inclusive
        # angle extrema to make oblique cases visible in the event menu.
        if len(vals)<=max_examples:
            for r in vals:take(r)
        else:
            for k in range(max_examples):
                ix=round(k*(len(vals)-1)/(max_examples-1))
                take(vals[ix])
            # Ensure one example with late secondary wire, provided this
            # changes the first-100ns multiplicity. Keep total fixed.
            late=[r for r in vals if int(r["N_t100"])<2]
            if late and len(chosen)==max_examples:
                pick=late[len(late)//2]
                existing=sorted(chosen,key=lambda identity:
                    abs(float(chosen[identity]["theta_u_in_deg"])))
                # Replace one interior muon-angle case (not both extremes).
                if len(existing)>2:
                    chosen.pop(existing[len(existing)//2])
                    take(pick)
        result[key]=list(chosen.values())[:max_examples]
    return result


def build(args):
    if args.max_per_chamber<2:
        raise ValueError("Need at least two examples per chamber")
    for p in (args.pilot,args.garfield,args.manifest,args.heed_tracks,
              args.heed_clusters,args.template,args.export_dir):
        if not p.exists():
            raise FileNotFoundError("Missing existing native centered pilot/model: "+str(p))

    pool,geometry_meta=load_centered_map(args.garfield,args.manifest)
    lam,marks=read_heed(args.heed_tracks,args.heed_clusters)
    gain_bank=empirical_gains(args.gain_csv) if args.gain_csv else None
    cfg=AnodeConfig()
    cfg.check()

    with args.pilot.open(newline="",encoding="utf-8") as f:
        reader=csv.DictReader(f)
        needed={"source","event","hit_index","track_id","station",
                "det_id","N_all","N_t100","theta_u_in_deg",
                "theta_w_in_deg","edep_eV","gas_path_mm"}
        if not needed.issubset(reader.fieldnames or []):
            raise ValueError("Incomplete first-chunk centered pilot CSV")
        candidates=choose_examples(reader,args.max_per_chamber)
    if not candidates:
        raise ValueError("No two-anode native Geant4 gas hits in pilot CSV")

    by_source=defaultdict(dict)
    for (source,det),records in candidates.items():
        for pilot in records:
            identity=(int(pilot["event"]),int(pilot["hit_index"]),
                      int(pilot["track_id"]))
            if identity in by_source[source]:
                raise ValueError("Repeated native event/hit/track identity")
            by_source[source][identity]=(det,pilot)

    examples=defaultdict(list)
    processed=0
    for source,channel in SOURCES.items():
        input_path=args.export_dir/f"{channel}_geant_gas_crossings.csv"
        if not input_path.is_file():
            raise FileNotFoundError("Missing native first-chunk ROOT export "+str(input_path))
        found=set()
        with input_path.open(newline="",encoding="utf-8") as f:
            reader=csv.DictReader(f)
            for raw in reader:
                identity=(int(raw["event"]),int(raw["hit_index"]),
                          int(raw["track_id"]))
                if identity not in by_source[source]:
                    continue
                if identity in found:
                    raise ValueError("Nonunique native source/event/hit/track")
                found.add(identity)
                det,pilot=by_source[source][identity]
                st=int(raw["station"])
                if (abs(int(raw["pdg"]))!=13 or int(raw["det_id"])!=det or
                    st!=GEOMETRY[det].station):
                    raise ValueError("Native PDG or chamber ID mismatch")
                if (abs(float(raw["theta_u_in_deg"])-
                        float(pilot["theta_u_in_deg"]))>1.e-8 or
                    abs(float(raw["theta_w_in_deg"])-
                        float(pilot["theta_w_in_deg"]))>1.e-8):
                    raise ValueError("Native vs pilot entrance angles disagree")
                rng=np.random.default_rng(stable_seed(
                    "centered_native_pilot_v1",source,"chunk0001",
                    *identity))
                outcome=ionize_and_drift(raw,rng,lam,marks,pool)
                mult=multiplicities(outcome,CUTOFFS)
                if mult[-1]!=2 or mult[1]!=int(pilot["N_t100"]):
                    raise ValueError("Replayed Heed/Garfield seed does not match "
                                     "frozen PI native pilot anode multiplicity")
                wires=sorted(set(int(x) for x in outcome["anodes"]))
                if len(wires)!=2:
                    raise AssertionError("Pilot required exactly two collected anodes")
                rng2=np.random.default_rng(stable_seed(
                    "PI_TWO_ANODE_UNCALIBRATED_DIGITIZER_V1",
                    source,"chunk0001",*identity))
                endpoints=zip(outcome["anodes"].tolist(),
                              outcome["times"].tolist())
                wf=digitize_hit(endpoints,rng2,cfg,gain_bank)
                a,b=wf[str(wires[0])],wf[str(wires[1])]
                if len(a["samples_mV"])!=len(b["samples_mV"]):
                    raise AssertionError("Independent anode sample axes differ")
                # Illustrative arithmetic sum of two independent output
                # voltages. No actual analog hardware summing simulated.
                combined=[round(x+y,3) for x,y in
                          zip(a["samples_mV"],b["samples_mV"])]
                early=[]
                for wire in wires:
                    ts=outcome["times"][outcome["anodes"]==wire]
                    valid=ts[np.isfinite(ts)&(ts>=0)]
                    early.append(round(float(np.min(valid)),3)
                                 if valid.size else None)
                delta=(round(early[1]-early[0],3)
                       if None not in early else None)
                c=GEOMETRY[det]
                source_key=f"{source}|{det}"
                example={
                    "src":source,"st":st,"det":det,
                    "layer":c.layer,
                    "event":identity[0],"track":identity[2],
                    "hitIndex":identity[1],
                    "pdg":int(raw["pdg"]),
                    "thetaU":round(float(raw["theta_u_in_deg"]),4),
                    "thetaW":round(float(raw["theta_w_in_deg"]),4),
                    "edep_eV":round(float(raw["edep_eV"]),3),
                    "path_mm":round(float(raw["path_length_mm"]),4),
                    "nPrimary":outcome["nprimary"],
                    "nCollected":outcome["ncollected"],
                    "N_t100":mult[1],
                    "wires":wires,
                    "anodes":[
                        {**v,"wire":wire,
                         "firstElectron_ns":early[i],
                         "arrivalTimes_ns":sorted(
                            round(float(t),3) for t in
                            outcome["times"][(outcome["anodes"]==wire)&
                            np.isfinite(outcome["times"])&
                            (outcome["times"]>=0)]
                          )}
                        for i,(v,wire) in enumerate(zip((a,b),wires))
                    ],
                    "sumSamples_mV":combined,
                    "deltaFirstElectron_ns":delta,
                }
                examples[source_key].append(example)
                processed+=1
        expected=set(by_source[source])
        if found!=expected:
            raise ValueError(f"Missing selected native Geant4 gas hits "
                             f"in {source}: {len(expected-found)} of {len(expected)}")
        print(f"{source}: replayed {len(found):,} actual Geant4 two-anode "
              f"muon segments (fixed native Heed seed), "
              "analog outputs distinct and summed for illustration",flush=True)

    for arr in examples.values():
        arr.sort(key=lambda v:(v["thetaU"],v["event"],v["hitIndex"]))
    # Small PI event catalog complements the interactive multi-megabyte
    # standalone HTML, permitting quick choice of angle and lag examples.
    event_catalog=args.html.parent/"pi_two_anode_event_catalog.csv"
    event_catalog.parent.mkdir(parents=True,exist_ok=True)
    fields=("source","station","det_id","event","track_id","hit_index",
            "theta_u_deg","theta_w_deg","wire_a","wire_b",
            "electrons_a","electrons_b","first_a_ns","first_b_ns",
            "delta_first_b_minus_a_ns","both_channels_over_5mV",
            "n_anodes_with_primary_electrons_by_100ns")
    with event_catalog.open("w",newline="",encoding="utf-8") as stream:
        writer=csv.DictWriter(stream,fieldnames=fields);writer.writeheader()
        for events in examples.values():
            for v in events:
                a,b=v["anodes"]
                writer.writerow({
                    "source":v["src"],"station":v["st"],
                    "det_id":v["det"],"event":v["event"],
                    "track_id":v["track"],"hit_index":v["hitIndex"],
                    "theta_u_deg":v["thetaU"],"theta_w_deg":v["thetaW"],
                    "wire_a":v["wires"][0],"wire_b":v["wires"][1],
                    "electrons_a":a["nPrimaryCollected"],
                    "electrons_b":b["nPrimaryCollected"],
                    "first_a_ns":v["anodes"][0]["firstElectron_ns"],
                    "first_b_ns":v["anodes"][1]["firstElectron_ns"],
                    "delta_first_b_minus_a_ns":v["deltaFirstElectron_ns"],
                    "both_channels_over_5mV":int(
                        any(x>=5 for x in a["samples_mV"]) and
                        any(x>=5 for x in b["samples_mV"])),
                    "n_anodes_with_primary_electrons_by_100ns":v["N_t100"],
                })
    layout=[{
        "st":c.station,"det":c.det_id,
        "x":c.x,"y":c.y,"hx":c.half_x,"hy":c.half_y,
        "row":c.row,"col":c.col,"layer":c.layer,
    } for group in build_chambers() for c in group]
    payload={
        "meta":{
            "geometry_id":GEOMETRY_ID,
            "sample":"Selected 2-wire true native Geant4 gas-hit segments from "
                     "first 3000 events per source. CURATED, not unbiased.",
            "field_sha256":geometry_meta["response_csv_sha256"],
            "gainSource":("unconditioned Stage A" if gain_bank is not None
                          else "uncalibrated Polya surrogate"),
            "digitizer":{k:v for k,v in cfg.__dict__.items()},
            "waveModel":"uncalibrated surrogate induced electron/ion charge "
                        "+ nominal CR-RC shaper, independent per wire; "
                        "NOT measured anode weighting current/ASIC",
            "hypotheticalSum":"Arithmetic sum of two independent channel "
                              "output waveforms; NOT a measured/common "
                              "summed-anode hardware channel",
            "timeOrigin":"t=0 at local Heed primary electron birth "
                         "approximation; detailed cluster creation time offsets "
                         "were not simulated",
            "units":{"time":"ns","voltage":"mV",
                     "muonAngles":"degrees"},
        },
        "layout":layout,"sources":list(SOURCES),
        "examples":dict(examples),
    }
    html=args.template.read_text(encoding="utf-8")
    marker="__PI_TWO_ANODE_EVENT_DATA__"
    if html.count(marker)!=1:
        raise ValueError("Expected exactly one two-anode data marker")
    args.html.parent.mkdir(parents=True,exist_ok=True)
    args.html.write_text(
        html.replace(marker,json.dumps(payload,ensure_ascii=False,
            separators=(",",":")).replace("</","<\\/")),
        encoding="utf-8")
    report={
        "geometry_id":GEOMETRY_ID,
        "source_pi_pilot":str(args.pilot),
        "pilot_sha256":hashlib.sha256(args.pilot.read_bytes()).hexdigest(),
        "centered_garfield_sha256":geometry_meta["response_csv_sha256"],
        "source_counts":{src:sum(len(v) for k,v in examples.items()
                             if k.startswith(src+"|")) for src in SOURCES},
        "two_wire_example_count":processed,
        "conservation_checks":"native event/hit/track IDs, fixed Heed seed "
                              "replay N_all==2, N_t100 matching",
        "sum_is_illustrative_not_real_hardware":True,
        "wave_model":"UNVALIDATED parametric gain/electron+ion+CR-RC waveform",
        "no_cathode_strips":True,
    }
    (args.html.parent/"pi_two_anode_signal_explorer_manifest.json").write_text(
        json.dumps(report,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(f"PASS: {processed} interactive native two-anode examples "
          "with independent waveforms and hypothetical summed output")
    print("Output:",args.html)
    print("Sortable PI two-wire angular/timing event catalog:",event_catalog)
    return report


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pilot",type=Path,default=INPUT)
    ap.add_argument("--export-dir",type=Path,default=EXPORT)
    ap.add_argument("--garfield",type=Path,
                    default=CENTERED/"c_transport_box_electrons.csv")
    ap.add_argument("--manifest",type=Path,
                    default=CENTERED/"pi_centered_transport_manifest.json")
    ap.add_argument("--heed-tracks",type=Path,
                    default=BASE/"c_heed_library_2p4/heed_2p4_tracks.csv")
    ap.add_argument("--heed-clusters",type=Path,
                    default=BASE/"c_heed_library_2p4/heed_2p4_clusters.csv")
    ap.add_argument("--gain-csv",type=Path,default=None)
    ap.add_argument("--max-per-chamber",type=int,default=N_SAMPLES)
    ap.add_argument("--template",type=Path,default=TEMPLATE)
    ap.add_argument("--html",type=Path,default=VIEW)
    args=ap.parse_args()
    build(args)


if __name__=="__main__":
    main()
