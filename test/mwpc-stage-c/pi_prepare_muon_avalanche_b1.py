#!/usr/bin/env python3
"""Versioned, no-Geant4 selection of REAL native J/psi muon gas crossings for B1.

Reads the existing frozen Geant4 CSV (including already-exported 3D momentum).
Selects deterministic, angle-stratified physical crossings. Does NOT generate
independent random angles or claim the selected pilot follows the full angle PDF.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path

from pi_2m_native_six_geometry_fold import crossing_provenance, SOURCE
from pi_centered_electron_sampler import REQUIRED_HIT

HERE=Path(__file__).resolve().parent
REPO=HERE.parents[1]
OUT=REPO/"test_runs/mwpc_stage_c/pi_muon_avalanche_b1_v1"
EXTRA=("pu_in_GeV","pv_in_GeV","pw_in_GeV")
FIELDS=("source","chunk","event","hit_index","track_id","station","det_id","pdg",
        "u_in_mm","v_in_mm","w_in_mm","u_out_mm","v_out_mm","w_out_mm",
        "theta_u_in_deg","theta_w_in_deg","path_length_mm","edep_eV",
        "pu_in_GeV","pv_in_GeV","pw_in_GeV","momentum_GeV","heed_seed")
def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for b in iter(lambda:f.read(1024*1024),b""):h.update(b)
    return h.hexdigest()
def seed(*words):
    h=hashlib.sha256("\x1f".join(map(str,words)).encode()).digest()
    return 1+int.from_bytes(h[:4],"little")%2_000_000_000
def validate(row):
    if abs(int(row["pdg"]))!=13 or int(row["station"]) not in range(6):return None
    if not set(REQUIRED_HIT|set(EXTRA)).issubset(row):raise ValueError("missing native schema")
    vals={}
    for field in ("u_in_mm","v_in_mm","w_in_mm","u_out_mm","v_out_mm","w_out_mm",
                  "theta_u_in_deg","theta_w_in_deg","path_length_mm","edep_eV",*EXTRA):
        vals[field]=float(row[field])
    if not all(math.isfinite(t) for t in vals.values()):raise ValueError("nonfinite Geant4 physics")
    length=math.sqrt(sum((vals[f"{a}_out_mm"]-vals[f"{a}_in_mm"])**2 for a in "uvw"))
    if not .05<length<100 or abs(length-vals["path_length_mm"])>max(.02,length*.005):
        return None
    if max(abs(vals["v_in_mm"]),abs(vals["v_out_mm"]))>3.05:return None
    # 9-wire periodic avalanche model is not reliable next to finite
    # chamber edges. Require 4 pitches of clearance on either side.
    half_width=153.6 if int(row["station"])==0 else 256.0
    if max(abs(vals["u_in_mm"]),abs(vals["u_out_mm"]))>=half_width-16.0:
        return None
    p=math.sqrt(sum(vals[f"p{a}_in_GeV"]**2 for a in "uvw"))
    if p<.11 or p>100000:return None
    du=[(vals[f"{a}_out_mm"]-vals[f"{a}_in_mm"])/length for a in "uvw"]
    pdir=[vals[f"p{a}_in_GeV"]/p for a in "uvw"]
    dot=sum(a*b for a,b in zip(du,pdir))
    if dot<.999:return None # near-straight native gas chord must match incoming muon direction
    if max(abs(vals["theta_u_in_deg"]),abs(vals["theta_w_in_deg"]))>180:return None
    return p,vals
def make(max_events=8,chunk=1,output=None):
    if not 2<=max_events<=24:raise ValueError("2..24 real muon gas hits for small pilot")
    path,meta=crossing_provenance("Jpsi",chunk)
    output=Path(output or OUT/f"smoke{max_events}_chunk{chunk:04d}")
    if output.exists():raise FileExistsError(f"Preserving existing input pilot: {output}")
    candidates={}
    counts={}
    with path.open(newline="",encoding="utf-8") as f:
        rd=csv.DictReader(f)
        if not set(REQUIRED_HIT|set(EXTRA)).issubset(rd.fieldnames or []):
            raise ValueError("Original geant_muon_crossings.csv lacks expected 3D momenta")
        for row in rd:
            res=validate(row)
            if res is None:continue
            p,vals=res
            st=int(row["station"]);side=-1 if vals["theta_u_in_deg"]<0 else 1
            key=(st,side)
            counts[key]=counts.get(key,0)+1
            # Deterministic hash-ranked candidate = uniform within every stratum.
            keyid=(SOURCE["Jpsi"],chunk,row["event"],row["hit_index"],row["track_id"],st,row["det_id"])
            score=seed("muon-aval-b1-selection-v1",*keyid)
            # Keep the top two independent ORIGINAL Geant4 crossings per
            # signed-angle/station stratum. The previous 2/8 B1 smoke picks
            # remain byte-for-byte identical: the first ranked candidate in
            # each stratum is untouched.
            rec={k:row[k] for k in FIELDS if k not in ("source","chunk","momentum_GeV","heed_seed")}
            rec["source"]=SOURCE["Jpsi"];rec["chunk"]=chunk
            rec["momentum_GeV"]=f"{p:.12g}"
            rec["heed_seed"]=seed("real-Heed-once-v1",*keyid)
            pool=candidates.setdefault(key,[])
            pool.append((score,rec))
            # Stable sort on score alone: tie-breaking by input order is
            # exactly what the old if(score<best) smoke selector did.
            pool.sort(key=lambda x:x[0])
            if len(pool)>2:pool.pop()
    # Sweep all 12 strata before using any second-crossing candidate.
    # This keeps B1's first 2/8 selected physical crossings unchanged.
    ordered=[]
    for rank in (0,1):
        for mirror in (False,True):
            for st in range(6):
                side=(1 if st%2==0 else -1)*(-1 if mirror else 1)
                if (st,side) in candidates and len(candidates[st,side])>rank:
                    ordered.append(((st,side),candidates[st,side][rank][1]))
    selected=ordered[:max_events]
    if len(selected)<max_events:raise ValueError("Insufficient signed-angle x station samples")
    seen=set();output_rows=[]
    for key,rec in selected:
        ident=(rec["chunk"],rec["event"],rec["hit_index"],rec["track_id"],rec["station"],rec["det_id"])
        if ident in seen:raise ValueError("Duplicated physical gas crossing")
        seen.add(ident);output_rows.append(rec)
    output.mkdir(parents=True)
    filename=output/"native_muon_crossings_input.csv"
    with filename.open("w",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=FIELDS)
        writer.writeheader();writer.writerows(output_rows)
    manifest={
        "schema":"PI_TRUE_GEANT4_MUON_B1_INPUT_V1",
        "status":"complete","selection":"deterministic hash-ranked top 2 per station/signed-angle stratum; cycle 1 covers all 12 strata before cycle 2; alternating signs across MS0..MS5 then mirrored signs; physical u-edge 16mm excluded; NOT incidence-angle PDF representative",
        "channel":"Jpsi","source":SOURCE["Jpsi"],"chunk":chunk,"n_crossings":len(output_rows),
        "original_native_csv":str(path),"original_native_sha256":meta["csv_sha256"],
        "original_root_sha256":meta["immutable_input"]["input_root_sha256"],
        "original_merged_generated_parents":meta["immutable_input"]["n_generated_parent_events"],
        "angle_stratum_valid_candidates":{f"MS{s}_{'neg' if side<0 else 'pos'}":n
                                         for (s,side),n in sorted(counts.items())},
        "selected_ids":[[x[k] for k in ("chunk","event","hit_index","track_id","station","det_id")]
                        for x in output_rows],
        "selected_momenta_GeV":[float(x["momentum_GeV"]) for x in output_rows],
        "code_sha256":sha(__file__),"input_sha256":sha(filename),
        "primary_ionization":"TRUE TrackHeed per physical native Geant4 crossing; not marked Poisson fold",
        "geometry":"symmetric 3+3mm, four p3/p4 FW variants downstream",
    }
    (output/"input_manifest.json").write_text(json.dumps(manifest,indent=2)+"\n")
    print(f"PASS: {len(output_rows)} verified native J/psi muon crossings: {filename}")
    for k,r in selected:
        print(f"  MS{k[0]} theta_u={'negative' if k[1]<0 else 'nonnegative'}: "
              f"theta_u={float(r['theta_u_in_deg']):+.3f}, "
              f"theta_w={float(r['theta_w_in_deg']):+.3f}, "
              f"p={float(r['momentum_GeV']):.3f}GeV/c, L={float(r['path_length_mm']):.3f}mm")
    return output

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--events",type=int,default=8)
    ap.add_argument("--chunk",type=int,default=1)
    ap.add_argument("--out",type=Path)
    a=ap.parse_args()
    if a.chunk<1 or a.chunk>2:ap.error("Pilot only: frozen 100k J/psi parent chunks 1 or 2")
    make(a.events,a.chunk,a.out)

if __name__=="__main__":main()
