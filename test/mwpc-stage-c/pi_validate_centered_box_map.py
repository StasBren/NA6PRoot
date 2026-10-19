#!/usr/bin/env python3
"""Validate and stamp a FRESH centered 3+3mm Garfield transport box map.

Checks actual response CSV grid, exact per-box trials, local geometry indices
and run provenance. Never accepts old 2+4mm map as a centered map. No Garfield
or Geant4 calls are made by this validator.
"""
from __future__ import annotations
import argparse
from collections import Counter
import csv
import hashlib
import json
import math
from pathlib import Path
import re

N_U=20
N_V=30
U_MIN=-2.0
V_MIN=-3.0
DU=DV=0.2
REQUIRED={
    "iu","iv","box_id","u_center_mm","v_center_mm","sample",
    "collected","drift_time_ns","wire_index",
}


def validate_map(csv_file:Path, n_per_box:int, progress:Path|None=None):
    if n_per_box<=0:
        raise ValueError("Expected strictly positive electrons per box")
    if "gm3p0_gp3p0" not in str(csv_file.parent):
        raise ValueError("Centered map must be in a distinct gm3p0_gp3p0 directory")
    if not csv_file.is_file():
        raise FileNotFoundError(csv_file)
    seen=set()
    boxes=Counter()
    total=collected=unknown=0
    with csv_file.open(newline="",encoding="utf-8") as source:
        reader=csv.DictReader(source)
        missing=REQUIRED-set(reader.fieldnames or ())
        if missing:
            raise ValueError(f"Transport map missing {sorted(missing)}")
        for row in reader:
            iu,iv,bid,sample=map(int,(row["iu"],row["iv"],row["box_id"],row["sample"]))
            if not (0<=iu<N_U and 0<=iv<N_V and bid==iu*N_V+iv):
                raise ValueError(f"Invalid 3+3 cell indices: {iu=} {iv=} {bid=}")
            if not 0<=sample<n_per_box:
                raise ValueError(f"Invalid sample number in cell {bid}: {sample}")
            key=(bid,sample)
            if key in seen:
                raise ValueError(f"Repeated electron sample: {key}")
            seen.add(key)
            uc=float(row["u_center_mm"])
            vc=float(row["v_center_mm"])
            if (not math.isclose(uc,U_MIN+(iu+0.5)*DU,abs_tol=1e-7) or
                not math.isclose(vc,V_MIN+(iv+0.5)*DV,abs_tol=1e-7)):
                raise ValueError(f"Wrong centered local (u,v) for box {bid}")
            is_collected=int(row["collected"])
            if is_collected not in (0,1):
                raise ValueError("Collected outcome must be 0 or 1")
            if is_collected:
                # 9 explicit wire electrodes: -4 .. +4. Relative index
                # is still NOT a finite MS0/MS1 physical anode number.
                wire=int(row["wire_index"])
                if not -4<=wire<=4:
                    raise ValueError(f"Unphysical Garfield collection wire {wire}")
                try:
                    dt=float(row["drift_time_ns"])
                except (TypeError,ValueError):
                    dt=math.nan
                if math.isfinite(dt):
                    if dt<0:
                        raise ValueError("Negative Garfield collected-electron drift time")
                else:
                    unknown+=1
                collected+=1
            boxes[bid]+=1
            total+=1
    if len(boxes)!=N_U*N_V:
        raise ValueError(f"Map has {len(boxes)}/600 Garfield start boxes")
    inconsistent=[(bid,n) for bid,n in boxes.items() if n!=n_per_box]
    if inconsistent:
        raise ValueError(f"Map incomplete: {inconsistent[:8]}")
    if total!=N_U*N_V*n_per_box or len(seen)!=total:
        raise ValueError(f"Map total samples {total}; expected {600*n_per_box}")

    if progress is not None:
        if not progress.is_file():
            raise FileNotFoundError(f"Missing Garfield run provenance log: {progress}")
        log=progress.read_text(encoding="utf-8",errors="replace")
        gap=re.findall(r"gap \(-,\+\) \[mm\]\s*:\s*([0-9.]+)\s*,\s*([0-9.]+)",log)
        fields=re.findall(r"\bB \[T\]\s*:\s*([-+0-9.eE]+)",log)
        if not gap or not fields:
            raise ValueError("Provenance log lacks gap/B settings")
        if any(not(math.isclose(float(gm),3.) and math.isclose(float(gp),3.))
               for gm,gp in gap):
            raise ValueError("Run log records non-centered gap settings")
        if any(not math.isclose(float(b),0.,abs_tol=1e-12) for b in fields):
            raise ValueError("B!=0: not the declared centered B=0 reference")
        if "COMPLETE" not in log:
            raise ValueError("No complete Garfield map run in progress log")
    return {
        "geometry_id":"MWPC_GARFIELD_CENTERED_3p0_3p0_B0_V1",
        "gap_minus_mm":3.,"gap_plus_mm":3.,
        "wire_pitch_mm":4.,"wire_diameter_um":30.,
        "anode_voltage_V":1800.,
        "gas":"Ar/CO2 70:30, 293.15 K, 760 Torr",
        "B_local_u_T":0.,
        "start_energy_eV":0.1,
        "u_bin_min_mm":U_MIN,"v_bin_min_mm":V_MIN,
        "du_mm":DU,"dv_mm":DV,
        "n_u":N_U,"n_v":N_V,"electrons_per_box":n_per_box,
        "total_response_rows":total,
        "collected_response_rows":collected,
        "collected_unknown_drift_times":unknown,
        "response_csv_sha256":hashlib.sha256(csv_file.read_bytes()).hexdigest(),
        "provenance_log_checked":progress is not None,
    }


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input",type=Path,required=True)
    ap.add_argument("--per-box",type=int,default=50)
    ap.add_argument("--progress",type=Path,required=True)
    ap.add_argument("--manifest",type=Path,required=True)
    args=ap.parse_args()
    summary=validate_map(args.input,args.per_box,args.progress)
    args.manifest.parent.mkdir(parents=True,exist_ok=True)
    args.manifest.write_text(json.dumps(summary,indent=2)+"\n",encoding="utf-8")
    print(f"PASS: validated {summary['total_response_rows']} independent centered "
          f"3+3mm Garfield single-electron transport responses across 600 boxes")
    print(f"COLLECTED {summary['collected_response_rows']}; missing drift time "
          f"{summary['collected_unknown_drift_times']}")
    print("CENTERED GEOMETRY ID:",summary["geometry_id"])
    print("Manifest:",args.manifest)


if __name__=="__main__":
    main()
