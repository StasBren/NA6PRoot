#!/usr/bin/env python3
"""Read-only digest of sealed 16k selected-cell RKF+Polya crossfit.

Does not change any physics, rerun the crossfit, or touch old archives.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path

BASE = Path(__file__).resolve().parents[2] / "test_runs/mwpc_stage_c"
DEFAULT = BASE / "pi_RKF_Polya_cell_crossfit_v1"
GEOMS = ("p3_noFW","p3_FW75","p4_noFW","p4_FW75")
METRICS = ("cell_mean_absolute_error","one_draw_MAE","PI90_coverage",
           "upper95_exceedance","upper99_exceedance","q95_pinball",
           "q99_pinball","zero_Brier")
DELTA = ("cell_mean_absolute_error","PI90_coverage","q95_pinball",
         "q99_pinball","zero_Brier")

def sha(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for block in iter(lambda:f.read(1024*1024),b""):h.update(block)
    return h.hexdigest()

def check(ok,explanation):
    if not ok:raise ValueError(explanation)

def read(p):return json.loads(p.read_text())

def fmt(s,percentage=False):
    z=100 if percentage else 1
    t=f'{s["value"]*z:.4f}'
    a,b=(x*z for x in s["CI95_cell_bootstrap"])
    return f'{t} [95% CI {a:.4f}, {b:.4f}]'

def report(directory):
    root=Path(directory)
    manifest=read(root/"manifest.json")
    check(manifest.get("schema")=="NA60_RKF_POLYA_CELLS_SHA_V1",
          "Wrong or missing SHA receipt; abort")
    files=("summary.json","fold_parameters.json","per_cell_results.json")
    check(set(manifest.get("output_SHA256",{}))==set(files),
          "Incomplete output SHA file list")
    for name in files:
        check(sha(root/name)==manifest["output_SHA256"][name],
              "Output changed since sealed analysis: "+name)
    data=read(root/"summary.json")
    cells=read(root/"per_cell_results.json")
    folds=read(root/"fold_parameters.json")
    check(data["schema"]=="NA60_RKF_POLYA_CELLS_V1" and
          data["status"]=="complete" and
          data["total_selected_cells"]==200 and
          data["total_validation_avalanches"]==16000,
          "Wrong physical study/cardinality")
    check(len(cells)==200 and len(folds)==8,"Incomplete files")
    ids=set()
    for r in cells:
        k=(r["geometry"],int(r["box_id"]))
        check(k[0] in GEOMS and k not in ids and r["n"]==80,
              "Bad validation cell identity/count")
        ids.add(k)
    check(all(sum(g==r["geometry"] for r in cells)==50 for g in GEOMS),
          "Missing geometry")
    for g in GEOMS:
        s=data["summary"][g]
        check(s["n_validation_cells"]==50 and
              s["n_validation_avalanches"]==4000,
              "Invalid summary "+g)
        print("\n"+"="*78)
        print(g,"(50 withheld-cell predictions; 80 full Garfield observations/cell)")
        print("model order: RKF+Polya | original A20 | geometry-only Polya")
        for m in METRICS:
            x=[s[tag][m]["value"] for tag in
               ("RKF_Polya","A20","geom_Polya")]
            unit="%" if m in ("PI90_coverage","upper95_exceedance",
                              "upper99_exceedance") else ""
            a=[round(v*100 if unit else v,4) for v in x]
            print(f"  {m:27s} {a[0]:>11} | {a[1]:>11} | {a[2]:>11} {unit}")
        print("  Positive Δ means RKF+Polya has a LARGER metric.")
        print("  For MAE/pinball/Brier: Δ<0 is better. Coverage Δ≈0 preferred.")
        for control in ("geom_Polya","A20"):
            print(f"  RKF+Polya minus {control} (paired-cell bootstrap):")
            for m in DELTA:
                result=s["RKF_minus_"+control][m]
                print(f"    {m:25s}: {fmt(result,percentage=(m=='PI90_coverage'))}")
    print("\nCAVEATS: 200 cells were previously selected by A20/old muons;")
    print("16k selected-cell validation is NOT a new-geometry or blind muon test.")
    print("RKF and geometry-only share p_zero by construction.")
    print("Nominal PI90 coverage near 95% may follow a ~17% atom at zero.")
    print("Signal, ion-cloud coordinates and arrival-time are UNTESTED.")
    print("PASS SHA-verified read-only summary; no Garfield and no new outputs.")

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--dir",type=Path,default=DEFAULT)
    a=p.parse_args()
    report(a.dir)
if __name__=="__main__":main()
