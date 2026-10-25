#!/usr/bin/env python3
"""Diagnose the ZERO field-wire endpoint count and reduced anode collection.

Read-only, offline side-by-side comparison of all original B=0, 20-sample
Garfield LUTs (or 50 sample LUTs via --lut-root). Pair trials exactly by
(pitch, box_id, sample) and retain Garfield endpoint status, initial electron
positions, and whether each electron reaches an anode/field wire.

This diagnostic does not rerun Garfield/Heed/Geant4, modify any LUT, or alter
the running 6M production pipeline. Count differences are microscopic
trial-level, not independent detector-event-level statistical uncertainties.

Example:
 python3 test/mwpc-stage-c/pi_field_wire_collection_loss_diagnostic.py
 python3 test/mwpc-stage-c/pi_field_wire_collection_loss_diagnostic.py \
    --lut-root test_runs/mwpc_stage_c/pi_pitch_fieldwire_50_v1/pilot \
    --out test_runs/mwpc_stage_c/pi_fw_loss_diagnostic_50_v1

Garfield documented status: -7 is attachment, -1 left drift area,
-5 left drift medium, -16 transport cut, -17 time window.
Noncollected status 0 can also be a geometric non-anode endpoint.
"""
from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
import json
import math
from pathlib import Path

from pi_pitch_fieldwire_pilot import PITCHES, tag, DU, NV, ROOT
from run_pi_geant4_campaign import REPO

CLASS={-7:"attachment (-7)",-1:"left drift area (-1)",
       -5:"outside drift medium (-5)",-16:"transport cut (-16)",
       -17:"time window (-17)",-3:"calculation error (-3)"}


def get_state(row):
    anode=int(row["collected"])
    field=int(row["field_wire_collected"])
    if anode+field>1 or anode not in (0,1) or field not in (0,1):
        raise ValueError("Invalid mutually exclusive electrode endpoints")
    if anode:return "anode"
    if field:return "field_wire"
    status=int(row["status"])
    u=float(row["u_end_mm"])
    v=float(row["v_end_mm"])
    if status==0 and not (math.isfinite(u) and math.isfinite(v)):
        return "missing endpoint (status 0)"
    return CLASS.get(status,f"other (status {status})")


def read_lut(path,pitch,samples):
    if not path.is_file():raise FileNotFoundError(path)
    by_key={}
    with path.open(newline="",encoding="utf-8") as f:
        for row in csv.DictReader(f):
            key=(int(row["box_id"]),int(row["sample"]))
            if key in by_key:raise ValueError(f"Repeated LUT trial {path} {key}")
            if (not 0<=key[0]<round(pitch/DU)*NV or
                not 0<=key[1]<samples or
                int(row["iu"])*NV+int(row["iv"])!=key[0]):
                raise ValueError(f"Invalid pitch/box/sample indexing: {path} {key}")
            by_key[key]=row
    expected=round(pitch/DU)*NV*samples
    if len(by_key)!=expected:raise ValueError(
        f"Incomplete grid: {path}: {len(by_key)}/{expected}")
    return by_key


def analyze_pair(pitch,root):
    p0=root/tag(pitch,0.)
    pf=root/tag(pitch,75.)
    meta0=json.loads((p0/"frozen_config.json").read_text(encoding="utf-8"))
    metaf=json.loads((pf/"frozen_config.json").read_text(encoding="utf-8"))
    for k in ("pitch_mm","samples_per_box","gas","magnetic_field_t",
              "gap_minus_mm","gap_plus_mm","anode_voltage_v","anode_diameter_um",
              "du_mm","dv_mm","n_u","n_v","mode","seed",
              "cpp_sha256","binary_sha256"):
        if meta0.get(k)!=metaf.get(k):
            raise ValueError(f"Incompatible independently frozen LUTs: {k}")
    if (meta0.get("pitch_mm")!=pitch or
        meta0.get("field_wire_diam_um")!=0. or
        metaf.get("field_wire_diam_um")!=75. or
        meta0.get("magnetic_field_t")!=0.):
        raise ValueError("Unexpected pitch/FW/B physics")
    samples=int(meta0["samples_per_box"])
    baseline=read_lut(p0/"garfield_electrons.csv",pitch,samples)
    with_fw=read_lut(pf/"garfield_electrons.csv",pitch,samples)
    if baseline.keys()!=with_fw.keys():
        raise ValueError("The two Garfield samples do not share (box,sample)")
    outcomes0=Counter()
    outcomesf=Counter()
    transitions=Counter()
    changed_starts=0
    net_cell=defaultdict(lambda:{"n":0,"anode_without":0,"anode_with":0,
                                  "attachment_without":0,"attachment_with":0})
    for k,r0 in baseline.items():
        rf=with_fw[k]
        state0=get_state(r0)
        statef=get_state(rf)
        outcomes0[state0]+=1
        outcomesf[statef]+=1
        transitions[(state0,statef)]+=1
        if (r0["random_seed"]!=rf["random_seed"]):
            raise ValueError(f"Paired electron RNG seed mismatch: {k}")
        if (abs(float(r0["u_start_mm"])-float(rf["u_start_mm"]))>1e-9 or
            abs(float(r0["v_start_mm"])-float(rf["v_start_mm"]))>1e-9):
            changed_starts+=1
        i,iv=divmod(k[0],NV)
        d=net_cell[(i,iv)]
        d["n"]+=1
        d["anode_without"]+=state0=="anode"
        d["anode_with"]+=statef=="anode"
        d["attachment_without"]+=state0=="attachment (-7)"
        d["attachment_with"]+=statef=="attachment (-7)"
    n=len(baseline)
    a0=outcomes0["anode"]
    af=outcomesf["anode"]
    lost=sum(q for (a,b),q in transitions.items() if a=="anode" and b!="anode")
    gained=sum(q for (a,b),q in transitions.items() if a!="anode" and b=="anode")
    if a0-af!=lost-gained:
        raise AssertionError("Paired conservation failed")
    delta=(af-a0)/n
    return {
        "pitch_mm":pitch,"samples_per_cell":samples,"n_independent_garfield_trials":n,
        "fraction_anode_without_fw":a0/n,
        "fraction_anode_with_fw":af/n,
        "delta_anode_collection_percentage_points":100*delta,
        "field_wire_collection_with_fw":outcomesf["field_wire"],
        "lost_from_anode_on_adding_fw":lost,
        "gained_to_anode_on_adding_fw":gained,
        "net_anode_loss":lost-gained,
        "pairs_with_different_start_positions":changed_starts,
        "states_without_fw":dict(sorted(outcomes0.items())),
        "states_with_fw":dict(sorted(outcomesf.items())),
        "state_deltas":{k:outcomesf[k]-outcomes0[k]
                        for k in sorted(outcomes0.keys()|outcomesf.keys())},
        "transitions":[{"without_fw":a,"with_fw":b,"count":number}
             for (a,b),number in transitions.most_common()],
        "spatial_grid":[{
            "iu":i,"iv":iv,
            "u_center_mm":-pitch/2+(i+.5)*DU,
            "v_center_mm":-3.+(iv+.5)*.2,
            "n":d["n"],
            "fraction_anode_without":d["anode_without"]/d["n"],
            "fraction_anode_with":d["anode_with"]/d["n"],
            "delta_anode":(d["anode_with"]-d["anode_without"])/d["n"],
            "attachment_without":d["attachment_without"]/d["n"],
            "attachment_with":d["attachment_with"]/d["n"],
        } for (i,iv),d in sorted(net_cell.items())],
    }


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lut-root",type=Path,default=ROOT/"pilot")
    parser.add_argument("--out",type=Path,
                        default=REPO/"test_runs/mwpc_stage_c/pi_fw_loss_diagnostic_20_v1")
    args=parser.parse_args()
    args.out.mkdir(parents=True,exist_ok=True)
    results=[analyze_pair(p,args.lut_root) for p in PITCHES]
    summary={
        "schema":"NA60_STAGEC_PAIRED_FIELDDWIRE_LOSS_DIAGNOSTIC_V1",
        "source_lut_root":str(args.lut_root),
        "physical_note":"Direct microscopic Garfield 20- or 50-trial per box "
           "uniform spatial grid, NOT the weighted J/psi native Heed population; "
           "status -7 is electron attachment; zero field-wire collection "
           "is zero within this specific endpoint classification",
        "results":results,
    }
    tmp=args.out/"paired_field_wire_loss.json.incomplete"
    tmp.write_text(json.dumps(summary,indent=2,allow_nan=False)+"\n",
                   encoding="utf-8")
    tmp.replace(args.out/"paired_field_wire_loss.json")
    for row in results:
        p=row["pitch_mm"]
        print(f"\n{p:.1f} mm · {row['n_independent_garfield_trials']:,} "
              "microscopic trials")
        print(f"  anode: {100*row['fraction_anode_without_fw']:.3f}% → "
              f"{100*row['fraction_anode_with_fw']:.3f}% "
              f"(Δ={row['delta_anode_collection_percentage_points']:+.3f} pp)")
        print(f"  lost anode: {row['lost_from_anode_on_adding_fw']:,} "
              f"gained anode: {row['gained_to_anode_on_adding_fw']:,}; "
              f"field wire: {row['field_wire_collection_with_fw']:,}")
        print(f"  paired start-position changes: "
              f"{row['pairs_with_different_start_positions']:,}")
        print("  Garfield physical endpoint classes (w/o FW → with FW):")
        for name in sorted(row["states_without_fw"].keys()|
                           row["states_with_fw"].keys()):
            print(f"    {name:30s}: "
                  f"{row['states_without_fw'].get(name,0):6d} → "
                  f"{row['states_with_fw'].get(name,0):6d} "
                  f"({row['state_deltas'][name]:+d})")
        for name in sorted(set(r["with_fw"] for r in row["transitions"]
                              if r["without_fw"]=="anode" and
                                 r["with_fw"]!="anode")):
            n=sum(r["count"] for r in row["transitions"]
                  if r["without_fw"]=="anode" and r["with_fw"]==name)
            print(f"    ANODE LOST TO {name:30s}: {n:,}")
    print("\nPASS: pure read-only 6-LUT Garfield endpoint-status "
          "and paired-(box,sample) diagnostic")
    print("Compact physical and spatial grid report:",
          args.out/"paired_field_wire_loss.json")


if __name__=="__main__":
    main()
