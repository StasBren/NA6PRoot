#!/usr/bin/env python3
"""Read-only physical audit of native TGeo gas segment lengths before Heed.

Never modifies ROOT, CSV, Heed/Garfield maps, existing fold shards or
frozen provenance.  Explicitly classifies:
 - zero path with entry == exit (physically zero ionisation length),
 - negative, nonfinite, overlong, or inconsistent path lengths,
 - out-of-gas endpoints, plus signed-angle nonfinite values.
CSV is the exact immutable ROOT export; use only source/chunk manifests.
"""
from __future__ import annotations
import argparse
import csv
import json
import math
from collections import Counter
from pathlib import Path

from pi_2m_native_six_geometry_fold import crossing_provenance
from pi_2m_export_completed_chunks import SOURCE

FIELDS=("u_in_mm","v_in_mm","w_in_mm","u_out_mm","v_out_mm","w_out_mm")
EXAMPLES_PER_CLASS=6


def classify(row):
    try:
        path=float(row["path_length_mm"])
        coords=[float(row[key]) for key in FIELDS]
        theta=[float(row["theta_u_in_deg"]),
               float(row["theta_w_in_deg"])]
    except (ValueError,TypeError,KeyError):
        return "unparseable_length_or_coordinates",None,None
    if not math.isfinite(path):
        return "nonfinite_length",path,None
    if not all(math.isfinite(x) for x in coords):
        return "nonfinite_endpoints",path,None
    u0,v0,w0,u1,v1,w1=coords
    geometric=math.sqrt((u1-u0)**2+(v1-v0)**2+(w1-w0)**2)
    if path==0. and geometric==0. and -3.05<=v0<=3.05:
        return "degenerate_zero_length_identical_endpoints",path,geometric
    if path<0.:
        return "negative_length",path,geometric
    if path==0.:
        return "zero_length_nonzero_displacement",path,geometric
    if path>=1000.:
        return "over_1000mm_length",path,geometric
    if max(abs(v0),abs(v1))>3.05:
        return "gas_normal_endpoint_outside_nominal_gap",path,geometric
    if not all(math.isfinite(t) and -180<=t<=180 for t in theta):
        return "nonfinite_or_invalid_signed_angles",path,geometric
    # Roundoff from source C++ 12-digit significant floating text can
    # be ~1e-9 mm absolute for 3D coordinates at high distances.
    # This flags gross mistakes, not natural decimal formatting error.
    if abs(path-geometric)>max(1.e-6,1.e-6*path):
        return "csv_length_endpoint_disagreement",path,geometric
    return "valid",path,geometric


def scan(channel,chunk,station,*,examples=EXAMPLES_PER_CLASS):
    csvfile,metadata=crossing_provenance(channel,chunk)
    counts=Counter()
    first={}
    n_station=0
    with csvfile.open(newline="",encoding="utf-8") as f:
        reader=csv.DictReader(f)
        for row in reader:
            if int(row["station"])!=station or abs(int(row["pdg"]))!=13:
                continue
            n_station+=1
            typ,value,geometric=classify(row)
            counts[typ]+=1
            if typ!="valid" and len(first.setdefault(typ,[]))<examples:
                first[typ].append({
                    "event":row["event"],
                    "hit_index":row["hit_index"],
                    "track_id":row["track_id"],
                    "det_id":row["det_id"],
                    "pdg":row["pdg"],
                    "length_mm":row["path_length_mm"],
                    "calculated_endpoint_distance_mm":geometric,
                    "endpoints_mm":{x:row.get(x) for x in FIELDS},
                    "edep_eV":row.get("edep_eV"),
                    "time_ns":row.get("time_ns"),
                    "status_start":row.get("status_start"),
                    "status_end":row.get("status_end"),
                })
    if not n_station:
        raise ValueError(f"No native muon crossings in {channel} {chunk} MS{station}")
    return {
        "source":SOURCE[channel],"channel":channel,
        "chunk":chunk,"station":station,
        "root_export_sha256":metadata["csv_sha256"],
        "n_muon_crossings":n_station,
        "categories":dict(counts),
        "offending_rows":first,
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--channel",choices=("Jpsi","Phi","Omega"),default="Jpsi")
    p.add_argument("--chunk",type=int,default=9)
    p.add_argument("--station",type=int,default=3)
    p.add_argument("--all-stations",action="store_true")
    p.add_argument("--out",type=Path,default=None)
    args=p.parse_args()
    if not 1<=args.chunk<=40 or args.station not in range(6):
        p.error("chunk 1..40; station 0..5")
    stations=range(6) if args.all_stations else (args.station,)
    reports=[scan(args.channel,args.chunk,st) for st in stations]
    for r in reports:
        print(f"PHYSICAL NATIVE GAS AUDIT: {r['channel']} "
              f"chunk{r['chunk']:04d} MS{r['station']}; "
              f"{r['n_muon_crossings']:,} muon gas crossings")
        for k,n in sorted(r["categories"].items()):
            print(f"  {k}: {n:,}")
            for row in r["offending_rows"].get(k,[]):
                print("    event={event} hit={hit_index} track={track_id} "
                      "det={det_id} L={length_mm} mm; "
                      "endpoint_d={calculated_endpoint_distance_mm} mm; "
                      "Edep={edep_eV} eV".format(**row))
    if args.out:
        args.out.parent.mkdir(parents=True,exist_ok=True)
        if args.out.exists():
            raise FileExistsError(f"Refuse overwrite diagnostic: {args.out}")
        args.out.write_text(json.dumps({"reports":reports},indent=2)+"\n")
        print("Saved immutable read-only diagnostic:",args.out)
    print("PASS: original ROOT/CSV untouched. "
          "Stop to interpret any non-valid category before changing Heed.")
if __name__=="__main__":
    main()
