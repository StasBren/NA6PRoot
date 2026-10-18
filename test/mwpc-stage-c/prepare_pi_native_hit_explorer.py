#!/usr/bin/env python3
"""Convert VERIFIED native Stage-C Geant4 gas-hit MCKine-joined CSV into
the same selected-track input schema as the historical event viewer.

Critically uses real TGeo chamber local u/v and Geant4 gas entry/exit
segments. Geant4 Edep is preserved INDEPENDENTLY of resampled Heed
primary ionization, avoiding double-counting as if deposited energy
were an exact electron-cluster measurement.

Only input chunks requested by the launcher are copied. No full-run
merging or new Geant4 or new Garfield execution.
"""
from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

from mwpc_checkerboard_projection import SOURCES, build_chambers

HERE=Path(__file__).resolve().parent
REPO=HERE.parents[1]
C={c.det_id:c for station in build_chambers() for c in station}
SOURCE_CHANNEL=dict(zip(SOURCES,("jpsi","phi","omega")))
FIELDS=(
    "source","job","event","track","pdg","station","det_id","layer",
    "x_in_cm","y_in_cm","z_in_cm","x_out_cm","y_out_cm","z_out_cm",
    "chamber_center_y_cm","gas_center_z_cm",
    "u_in_mm","v_in_mm","u_out_mm","v_out_mm",
    "gas_path_mm","theta_u_deg","theta_w_deg","p_GeV",
    "native_hit","edep_eV","time_ns","hit_index",
)


def convert(source_path,out_path,source,job):
    if not source_path.is_file():
        raise FileNotFoundError(source_path)
    count=0
    edep=0.
    skipped=0
    with source_path.open(newline="",encoding="utf-8") as inp, out_path.open(
        "w",newline="",encoding="utf-8"
    ) as target:
        reader=csv.DictReader(inp)
        required={
            "event","hit_index","track_id","pdg","station","det_id",
            "x_in_cm","y_in_cm","z_in_cm","x_out_cm","y_out_cm","z_out_cm",
            "u_in_mm","v_in_mm","u_out_mm","v_out_mm",
            "path_length_mm","edep_eV","time_ns",
            "px_in_GeV","py_in_GeV","pz_in_GeV",
            "theta_u_in_deg","theta_w_in_deg",
        }
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(f"{source_path}: missing {required-set(reader.fieldnames or [])}")
        writer=csv.DictWriter(target,fieldnames=FIELDS)
        writer.writeheader()
        for raw in reader:
            pdg=int(raw["pdg"])
            if abs(pdg)!=13:
                skipped+=1
                continue
            det,st=int(raw["det_id"]),int(raw["station"])
            if det not in C or C[det].station!=st:
                raise ValueError(f"Bad native station/det ID: MS{st}/DET{det}")
            xyz=tuple(float(raw[k]) for k in (
                "x_in_cm","y_in_cm","z_in_cm",
                "x_out_cm","y_out_cm","z_out_cm"))
            u0,u1=float(raw["u_in_mm"]),float(raw["u_out_mm"])
            v0,v1=float(raw["v_in_mm"]),float(raw["v_out_mm"])
            length=float(raw["path_length_mm"])
            deposit=float(raw["edep_eV"])
            xyzp=tuple(float(raw[k]) for k in (
                "px_in_GeV","py_in_GeV","pz_in_GeV"))
            p=math.sqrt(sum(x*x for x in xyzp))
            if (not all(math.isfinite(x) for x in
                        (*xyz,u0,u1,v0,v1,length,deposit,*xyzp)) or
                length<=0 or deposit<0 or p<=0):
                raise ValueError("Nonfinite/invalid native gas hit")
            # Global nominal MWPC orientation maps local u to Y and v to Z.
            # These centre values are needed only by old compatibility
            # columns: native Heed and wire folding USE TRUE LOCAL u/v.
            yc=0.5*(xyz[1]+xyz[4])-0.05*(u0+u1)
            zc=0.5*(xyz[2]+xyz[5])-0.05*(v0+v1)
            writer.writerow({
                "source":source,"job":job,"event":raw["event"],
                "track":raw["track_id"],"pdg":pdg,
                "station":st,"det_id":det,"layer":C[det].layer,
                **{k:raw[k] for k in (
                    "x_in_cm","y_in_cm","z_in_cm",
                    "x_out_cm","y_out_cm","z_out_cm",
                    "u_in_mm","v_in_mm","u_out_mm","v_out_mm",
                )},
                "chamber_center_y_cm":yc,"gas_center_z_cm":zc,
                "gas_path_mm":length,
                "theta_u_deg":raw["theta_u_in_deg"],
                "theta_w_deg":raw["theta_w_in_deg"],
                "p_GeV":p,"native_hit":1,"edep_eV":deposit,
                "time_ns":raw["time_ns"],"hit_index":raw["hit_index"],
            })
            count+=1; edep+=deposit
    print(f"{source}: {count:,} real gas muon hits exported; "
          f"{skipped:,} non-muon detector hits excluded; "
          f"Geant4 deposited {edep/1e6:.3f} MeV",flush=True)
    return count


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--export-dir",type=Path,required=True)
    ap.add_argument("--out",type=Path,required=True)
    ap.add_argument("--chunk",default="chunk0001")
    args=ap.parse_args()
    args.out.mkdir(parents=True,exist_ok=True)
    for source,ch in SOURCE_CHANNEL.items():
        path=args.export_dir/f"{ch}_geant_gas_crossings.csv"
        target=args.out/f"crossings_{source}_ms2_530.csv"
        count=convert(path,target,source,f"native_{ch}_{args.chunk}")
        if count==0:
            raise RuntimeError(f"Zero true muon gas hits in {source}")


if __name__=="__main__":
    main()
