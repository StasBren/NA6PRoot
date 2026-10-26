#!/usr/bin/env python3
"""Build a ONE-SLIDE ROOT-style PI presentation from real 20-vs-500 CSV.

Reads the previously computed, identical-exposure J/psi MS0
phi_comparison.csv; launches a ROOT C++ macro making two vertically
stacked TH1D histograms (20 above, 500 below), matched x/y scales and
individual 1-sigma Garfield-LUT errors. Saves ROOT, PDF and PNG plus a
single-slide PPTX when python-pptx is installed.

No Garfield, Geant4, Heed or frozen 720-fold jobs are run or changed.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import shutil
import subprocess
from pathlib import Path

HERE=Path(__file__).resolve().parent
REPO=HERE.parents[1]
INPUT=REPO/"test_runs/mwpc_stage_c/pi_pitch_fieldwire_20_vs_500_expected_v1/phi_comparison.csv"
OUTPUT=REPO/"test_runs/mwpc_stage_c/pi_lut_root_slide_v1"
MACRO=HERE/"pi_lut_root_slide.C"
GEOMS=("p3p0_field75p0um","p3p0_anode_only")


def validate_csv(path,geometry):
    with path.open(newline="",encoding="utf-8") as stream:
        r=csv.DictReader(stream)
        expected={"geometry","phi_lower_deg","phi_upper_deg",
                  "p20","se20","p500","se500"}
        if not expected.issubset(r.fieldnames or []):
            raise ValueError("Missing phi_comparison.csv fields")
        selected=[row for row in r if row["geometry"]==geometry]
    if len(selected)!=36:
        raise ValueError(f"Expected 36 angular bins, got {len(selected)}")
    sums={"p20":0.,"p500":0.}
    for i,row in enumerate(sorted(selected,key=lambda z:float(z["phi_lower_deg"]))):
        lo=float(row["phi_lower_deg"])
        hi=float(row["phi_upper_deg"])
        if not math.isclose(lo,-180.+i*10.,abs_tol=1e-6) or not \
           math.isclose(hi,lo+10.,abs_tol=1e-6):
            raise ValueError(f"Invalid angular bin {i}")
        for key in ("p20","se20","p500","se500"):
            v=float(row[key])
            if not math.isfinite(v) or v<0:
                raise ValueError(f"Invalid {key} in bin {i}")
        sums["p20"]+=float(row["p20"])
        sums["p500"]+=float(row["p500"])
    if any(abs(sums[name]-1.)>1e-4 for name in sums):
        raise ValueError(f"PDFs do not normalize: {sums}")
    return sums


def save_pptx(png,out):
    try:
        from pptx import Presentation
        from pptx.util import Inches
    except ImportError:
        print("NOTE: python-pptx not installed; PDF, PNG and ROOT will "
              "still be available. Optional: python3 -m pip install python-pptx",
              flush=True)
        return False
    ppt=Presentation()
    ppt.slide_width=Inches(13.333333)
    ppt.slide_height=Inches(7.5)
    slide=ppt.slides.add_slide(ppt.slide_layouts[6])
    slide.shapes.add_picture(str(png),0,0,width=ppt.slide_width,
                             height=ppt.slide_height)
    ppt.save(str(out))
    return True


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csv",type=Path,default=INPUT,
                        help="Original generated phi_comparison.csv")
    parser.add_argument("--geometry",choices=GEOMS,
                        default="p3p0_field75p0um")
    parser.add_argument("--out",type=Path,default=OUTPUT)
    parser.add_argument("--validate-only",action="store_true",
                        help="Check measured input without requiring ROOT")
    a=parser.parse_args()
    if not a.csv.is_file():
        raise FileNotFoundError(
            f"Missing real angular data: {a.csv}\n"
            "First run pi_lut_20_vs_500_expected.py --station 0 --chunks 2")
    normalizations=validate_csv(a.csv,a.geometry)
    print(f"PASS: measured CSV, 36 normalized angular bins, "
          f"geometry={a.geometry}; {normalizations}",flush=True)
    if a.validate_only:
        return
    if not MACRO.is_file():
        raise FileNotFoundError(f"Missing ROOT plotting macro: {MACRO}")
    root=shutil.which("root")
    if not root:
        raise EnvironmentError("ROOT executable 'root' not on PATH")
    a.out.mkdir(parents=True,exist_ok=True)
    prefix=a.geometry+"_slide"
    expected=[a.out/(prefix+"."+ext) for ext in ("png","pdf","root")]
    # Never use stale successful results to mask ROOT macro errors.
    for target in expected:
        if target.exists():
            target.unlink()
    # Standard WSL POSIX paths; no arbitrary shell invocation.
    def quote_root(p):
        s=str(p.resolve())
        if any(char in s for char in ('"',"\n","\r")):
            raise ValueError("Unexpected filename characters")
        return '"'+s+'"'
    command=(str(MACRO.resolve())+"("+
             quote_root(a.csv)+","+quote_root(a.out)+
             ',"'+a.geometry+'")')
    process=subprocess.run(
        [root,"-l","-b","-q",command],cwd=REPO,
        stdout=subprocess.PIPE,stderr=subprocess.STDOUT,
        text=True,check=False)
    print(process.stdout,flush=True)
    if process.returncode or "PASS: ROOT TH1D" not in process.stdout or \
       any(not p.is_file() or p.stat().st_size==0 for p in expected):
        raise RuntimeError("ROOT slide rendering failed; see preceding log")
    pptx=a.out/(prefix+".pptx")
    made=save_pptx(expected[0],pptx)
    report={
        "analysis":"J/psi MS0, p=3mm, original 20 vs independent 500",
        "geometry":a.geometry,"angular_bins":36,
        "data_sha256":hashlib.sha256(a.csv.read_bytes()).hexdigest(),
        "presentation_statement":"Expected per-cell independent Garfield "
            "statistical uncertainty scales as 1/sqrt(K), nominal 5x "
            "reduction from 20 to 500",
        "probability_source":"Fixed Geant4 muon crossing exposure; "
            "analytically expected Heed spatial occupancy; "
            "not the individual parent-event Heed realization",
        "error_bar_meaning":"Individual 1-sigma microscopic LUT uncertainty",
        "ROOT_pdf":str(expected[1]),
        "ROOT_file":str(expected[2]),
        "PNG":str(expected[0]),
        "PPTX":str(pptx) if made else None,
    }
    (a.out/(prefix+"_provenance.json")).write_text(
        json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print("PASS: one PI slide, actual ROOT histograms and 1-sigma error bars")
    for file in expected+([pptx] if made else []):
        print("OUTPUT",file,flush=True)


if __name__=="__main__":
    main()
