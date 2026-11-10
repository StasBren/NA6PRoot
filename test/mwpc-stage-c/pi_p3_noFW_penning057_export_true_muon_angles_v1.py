#!/usr/bin/env python3
"""Read-only compact true Geant4 1000-parent muon-direction export.

From the frozen same-Penning-ON native 1000 Geant4 crossing input,
recover *true* longitudinal incidence direction and theta_w_in_deg
to compare with longitudinal Gaussian ion source width. No Garfield.
"""
from __future__ import annotations
import csv
import hashlib
import json
import math
import sys
from pathlib import Path

BASE=Path(__file__).resolve().parents[2]
ROOT=BASE/"test_runs/mwpc_stage_c/pi_p3_noFW_jpsi1000_penning057_v1"
INPUT=ROOT/"input/native_muon_crossings_input.csv"
MANIFEST=ROOT/"input/input_manifest.json"
OUT=ROOT/"true_muon_incidence_v1"
CSV=OUT/"NA60_PENNING_ON_REAL_1000_MUON_DIRECTIONS.csv"
INFO=OUT/"provenance.json"
REQUIRED=(
    "source","chunk","event","hit_index","track_id","station",
    "det_id","pdg","u_in_mm","v_in_mm","w_in_mm",
    "u_out_mm","v_out_mm","w_out_mm",
    "theta_u_in_deg","theta_w_in_deg","path_length_mm",
    "edep_eV","pu_in_GeV","pv_in_GeV","pw_in_GeV",
    "momentum_GeV","heed_seed")
COLUMNS=(
    "parent_muon","crossing_id","station","det_id","pdg",
    "u_in_mm","v_in_mm","w_in_mm","u_out_mm","v_out_mm","w_out_mm",
    "delta_u_mm","delta_v_mm","delta_w_mm","chord_length_mm",
    "geant_path_length_mm","geant_theta_u_in_deg",
    "geant_theta_w_in_deg","true_unit_d_u",
    "true_unit_d_v","true_unit_d_w",
    "abs_true_unit_d_w",
    "angle_to_positive_w_axis_deg",
    "acute_angle_to_w_axis_deg",
    "inclination_to_plane_perpendicular_to_w_deg",
    "muon_momentum_GeV")

def need(ok,message):
    if not ok:raise RuntimeError(message)

def sha(file):
    digest=hashlib.sha256()
    with Path(file).open("rb") as f:
        for b in iter(lambda:f.read(1<<20),b""):
            digest.update(b)
    return digest.hexdigest()

def compute_row(i,r):
    a=[float(r[c+"_in_mm"]) for c in ("u","v","w")]
    b=[float(r[c+"_out_mm"]) for c in ("u","v","w")]
    ds=[b[j]-a[j] for j in range(3)]
    chord=math.sqrt(sum(x*x for x in ds))
    path=float(r["path_length_mm"])
    need(chord>.05 and abs(path-chord)<.005*chord+.02,
         "True original muon gas chord not geometrically consistent")
    dd=[x/chord for x in ds]
    mom=[float(r[c]) for c in ("pu_in_GeV","pv_in_GeV","pw_in_GeV")]
    pm=math.sqrt(sum(x*x for x in mom))
    need(pm>0 and sum(dd[j]*mom[j] for j in range(3))/pm>.999,
         "Geant4 muon chord and momentum have incompatible directions")
    w=max(-1.,min(1.,dd[2]))
    identifier="|".join([r["source"]]+[
        str(int(r[k])) for k in
        ("chunk","event","hit_index","track_id","station","det_id")
    ])
    values=dict(
        parent_muon=i,crossing_id=identifier,
        station=int(r["station"]),det_id=int(r["det_id"]),
        pdg=int(r["pdg"]),
        u_in_mm=a[0],v_in_mm=a[1],w_in_mm=a[2],
        u_out_mm=b[0],v_out_mm=b[1],w_out_mm=b[2],
        delta_u_mm=ds[0],delta_v_mm=ds[1],delta_w_mm=ds[2],
        chord_length_mm=chord,geant_path_length_mm=path,
        geant_theta_u_in_deg=float(r["theta_u_in_deg"]),
        geant_theta_w_in_deg=float(r["theta_w_in_deg"]),
        true_unit_d_u=dd[0],true_unit_d_v=dd[1],true_unit_d_w=w,
        abs_true_unit_d_w=abs(w),
        angle_to_positive_w_axis_deg=math.degrees(math.acos(w)),
        acute_angle_to_w_axis_deg=math.degrees(math.acos(abs(w))),
        inclination_to_plane_perpendicular_to_w_deg=
            math.degrees(math.asin(abs(w))),
        muon_momentum_GeV=float(r["momentum_GeV"])
    )
    need(set(values)==set(COLUMNS),"Internal true angle CSV schema changed")
    return values

def export():
    need(INPUT.is_file() and MANIFEST.is_file(),
         "Original selected real Geant4 1000 muons are not present")
    manifest=json.loads(MANIFEST.read_text())
    need(manifest["n_parents"]==1000 and
         manifest["geometry"]=="p3_noFW" and
         sha(INPUT)==manifest["selected_input_sha256"],
         "Real original Geant4 muon SHA differs from frozen Penning-ON input")
    need(not OUT.exists(),"Original true-angle export already exists; no overwrite")
    with INPUT.open(newline="") as f:
        reader=csv.DictReader(f)
        need(reader.fieldnames is not None and
             set(REQUIRED).issubset(reader.fieldnames) and
             len(set(reader.fieldnames))==len(reader.fieldnames),
             "Original real muon input has unexpected coordinate schema")
        data=[compute_row(i,r) for i,r in enumerate(reader)]
    need(len(data)==1000 and
         len(set(row["crossing_id"] for row in data))==1000,
         "Original true 1000 Geant4 parent crossing identities changed")
    OUT.mkdir(parents=False,exist_ok=False)
    with CSV.open("x",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=COLUMNS)
        writer.writeheader()
        writer.writerows(data)
    INFO.write_text(json.dumps({
        "schema":"NA60_PENNING_ON_REAL_1000_TRUE_DIRECTION_EXPORT_V1",
        "original_muon_source_sha256":sha(INPUT),
        "selected_manifest_sha256":sha(MANIFEST),
        "new_csv_sha256":sha(CSV),
        "sample":"1000 real Geant4 J/psi muon gas crossings, Penning ON selection",
        "physics_note":"true_unit_d_w=(w_out-w_in)/3D gas chord; "
                       "acute angle to wire axis = acos(abs(true_unit_d_w)); "
                       "this is NOT an avalanche-electron arrival angle",
    },indent=2)+"\n")
    print("PASS 1000 ORIGINAL true Geant4 muon directions; same verified SHA")
    print("CSV",CSV,"bytes",CSV.stat().st_size)
    print("No Garfield/Heed rerun, no modification to 1000 or 30k data")

def selftest():
    q={
        "source":"test","chunk":"3","event":"17","hit_index":"4",
        "track_id":"2","station":"0","det_id":"5","pdg":"13",
        "u_in_mm":"0","v_in_mm":"-3","w_in_mm":"0",
        "u_out_mm":"0","v_out_mm":"3","w_out_mm":"6",
        "path_length_mm":str(math.sqrt(72)),
        "pu_in_GeV":"0","pv_in_GeV":"1","pw_in_GeV":"1",
        "theta_u_in_deg":"0","theta_w_in_deg":"45",
        "momentum_GeV":str(math.sqrt(2)),
    }
    row=compute_row(0,q)
    need(abs(row["abs_true_unit_d_w"]-math.sqrt(.5))<1e-12 and
         abs(row["acute_angle_to_w_axis_deg"]-45)<1e-12 and
         row["crossing_id"]=="test|3|17|4|2|0|5",
         "True Geant4 longitudinal incidence calculation failed")
    print("PASS true direction / angle / crossing ID selftest; no file writes")

if __name__=="__main__":
    args=sys.argv[1:]
    need(len(args)==1 and args[0] in ("selftest","export"),
         "Usage: python SCRIPT selftest | export")
    if args[0]=="selftest":selftest()
    else:export()
