#!/usr/bin/env python3
"""Read-only, one-pass physical audit of ALL 120 full native Geant4 CSVs.

Finds every muon gas crossing rejected by the frozen Heed birth sampler
BEFORE expensive six-LUT folding. Only accepts EXACT zero-length same-point
in-gas chords with Edep <1 eV (Poisson(lambda*0)=0), records each such row.
All other abnormalities fail closed. Input CSV/checkpoint Sha256 validated
by crossing_provenance(). Original Geant4/ROOT/CSV never modified.
"""
from __future__ import annotations
import argparse
import csv
import json
import os
from collections import Counter
from pathlib import Path

from diagnose_pi_2m_native_gas_paths import classify
from pi_2m_native_six_geometry_fold import crossing_provenance
from pi_2m_zero_path_compatible_resume import coincident_zero_chord,VERSION,hashed
from run_pi_geant4_campaign import REPO

OUT=REPO/"test_runs/mwpc_stage_c/pi_weekend_2m_operations_v1/native_gas_path_audit.json"
CHANNELS=("Jpsi","Phi","Omega")


def audit():
    per_csv=[]
    exceptions=[]
    for channel in CHANNELS:
        for chunk in range(1,41):
            path,meta=crossing_provenance(channel,chunk)
            counts=Counter()
            muons=0
            for_all=0
            with path.open(newline="",encoding="utf-8") as f:
                reader=csv.DictReader(f)
                for row in reader:
                    for_all+=1
                    if abs(int(row["pdg"]))!=13:
                        continue
                    muons+=1
                    class_name,value,distance=classify(row)
                    if class_name=="valid":
                        counts["valid"]+=1
                        continue
                    if (class_name=="degenerate_zero_length_identical_endpoints" and
                        coincident_zero_chord(row)):
                        counts["accepted_zero_chord_sub_eV"]+=1
                        exceptions.append({
                            "source":meta["immutable_input"]["source"],
                            "channel":channel,"chunk":chunk,
                            "event":int(row["event"]),
                            "hit_index":int(row["hit_index"]),
                            "track_id":int(row["track_id"]),
                            "station":int(row["station"]),
                            "det_id":int(row["det_id"]),
                            "edep_eV":float(row["edep_eV"]),
                            "length_mm":0,
                            "native_CSV_SHA256":meta["csv_sha256"],
                        })
                        continue
                    bad={
                        "source":channel,"chunk":chunk,"station":row["station"],
                        "event":row["event"],"hit_index":row["hit_index"],
                        "track_id":row["track_id"],"class":class_name,
                        "edep_eV":row["edep_eV"],"path_length_mm":row["path_length_mm"],
                        "endpoint_distance_mm":distance,
                    }
                    raise RuntimeError(
                        "FAIL CLOSED: unexpected native Geant4/Heed crossing "
                        +json.dumps(bad,ensure_ascii=False))
            if muons!=meta["n_muon_gas_crossings"] or (
                for_all!=meta["n_gas_crossings"]):
                raise AssertionError(
                    f"Native CSV count changed: {channel} chunk{chunk:04d}")
            per_csv.append({
                "channel":channel,"chunk":chunk,
                "csv_sha256":meta["csv_sha256"],
                "n_gas_crossings":for_all,
                "n_muon_crossings":muons,
                "physics_categories":dict(counts),
            })
            print(f"AUDIT {channel} chunk{chunk:04d}: "
                  f"{muons:,} muon hits, "
                  f"{counts['accepted_zero_chord_sub_eV']} "
                  "strict zero-length exceptions",flush=True)
    return {
        "version":"PI_NATIVE_120_CSV_PHYSICAL_PATH_AUDIT_V1",
        "zero_path_policy":VERSION,
        "zero_path_wrapper_sha256":hashed(Path(__file__).parent/
                                         "pi_2m_zero_path_compatible_resume.py"),
        "complete_source_chunks":len(per_csv),
        "n_valid_muon_crossings":sum(r["physics_categories"].get("valid",0)
                                     for r in per_csv),
        "n_accepted_zero_chords":len(exceptions),
        "exceptions":exceptions,"chunks":per_csv,
        "interpretation":"All native gas crossing rows retained. Accepted "
           "0 mm chords with identical endpoints and Edep<1eV generate "
           "Poisson(0) Heed clusters; all other paths retain original model.",
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out",type=Path,default=OUT)
    args=p.parse_args()
    result=audit()
    args.out.parent.mkdir(parents=True,exist_ok=True)
    if args.out.is_file():
        existing=json.loads(args.out.read_text(encoding="utf-8"))
        if existing!=result:
            raise ValueError("Previous 120 CSV physics-path audit differs "
                             "from current frozen data; refusing overwrite")
    else:
        tmp=args.out.with_name(args.out.name+".incomplete")
        tmp.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",
                       encoding="utf-8")
        os.replace(tmp,args.out)
    print(f"PASS: validated {result['complete_source_chunks']}/120 "
          "original native exports and every muon gas path. "
          f"Strict 0 mm, Edep<1eV exceptions: {len(result['exceptions'])}."
          f" Ledger: {args.out}",flush=True)

if __name__=="__main__":
    main()
