#!/usr/bin/env python3
"""Build a ONE-SLIDE ROOT-style PI presentation from real 20-vs-500 CSV.

Reads the previously computed, identical-exposure J/psi MS0
phi_comparison.csv; draws two vertically stacked ROOT-style histograms
(20 above, 500 below) with matched axes and 1-sigma Garfield-LUT errors.
Matplotlib renders PNG, inserted into one-slide PPTX. No ROOT executable,
ROOT file or PDF is needed.

No Garfield, Geant4, Heed or frozen 720-fold jobs are run or changed.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

HERE=Path(__file__).resolve().parent
REPO=HERE.parents[1]
INPUT=REPO/"test_runs/mwpc_stage_c/pi_pitch_fieldwire_20_vs_500_expected_v1/phi_comparison.csv"
OUTPUT=REPO/"test_runs/mwpc_stage_c/pi_lut_root_slide_v1"
GEOMS=("p3p0_field75p0um","p3p0_anode_only")


def validate_csv(path: Path, geometry: str):
    """Return ordered, normalized 36-bin P(phi) and finite-LUT errors."""
    with path.open(newline="", encoding="utf-8") as handle:
        reader=csv.DictReader(handle)
        needed={"geometry","phi_lower_deg","phi_upper_deg",
                "p20","se20","p500","se500"}
        if not needed.issubset(reader.fieldnames or []):
            raise ValueError("Missing required angular CSV columns")
        rows=[row for row in reader if row["geometry"]==geometry]
    if len(rows)!=36:
        raise ValueError(f"Expected 36 angular bins for {geometry}, got {len(rows)}")
    rows.sort(key=lambda r:float(r["phi_lower_deg"]))
    values={k:[] for k in ("p20","se20","p500","se500")}
    for i,row in enumerate(rows):
        lo=float(row["phi_lower_deg"])
        hi=float(row["phi_upper_deg"])
        if not (math.isfinite(lo) and math.isfinite(hi)
                and math.isclose(lo,-180+10*i,abs_tol=1e-6)
                and math.isclose(hi,lo+10,abs_tol=1e-6)):
            raise ValueError(f"Wrong or duplicate phi bin {i}")
        for key in values:
            v=float(row[key])
            if not math.isfinite(v) or v<0:
                raise ValueError(f"Invalid {key} in bin {i}")
            values[key].append(v)
    for key in ("p20","p500"):
        if not math.isclose(sum(values[key]),1.,abs_tol=1e-4,rel_tol=0.):
            raise ValueError("Angular distribution not normalized: "+key)
    return values


def study_statistics(csv_path: Path, geometry: str, data: dict):
    """Read matched native Geant4/Heed exposure from the companion report.

    Never label a new slide with outdated or mismatched exposure counts.
    Collection is geometry-dependent; use LUT500 to state the approximate
    anode-collected exposure shown in the new (lower) histogram.
    """
    path=csv_path.with_name("comparison.json")
    if not path.is_file():
        raise FileNotFoundError(
            "Missing companion study metadata: "+str(path)+
            ". Re-run the isolated pi_lut_20_vs_500_expected.py comparison.")
    report=json.loads(path.read_text(encoding="utf-8"))
    scope=report["scope"]
    if scope.get("source")!="jpsi40-bore" or int(scope["station"])!=0:
        raise ValueError("Slide expects the J/psi MS0 native Geant4 study")
    if scope.get("heed_occupancy_mode")!="analytic_expectation":
        raise ValueError("Slide expects the analytically expected Heed exposure")
    if int(scope["selected_muon_gas_crossings"])!=int(
            scope["all_muon_gas_crossings_seen"]):
        raise ValueError("Slide statistics require complete selected chunk exposure")
    chunks=int(scope["chunks"])
    if not 1<=chunks<=40:
        raise ValueError("Unexpected number of J/psi Geant4 chunks")
    results=report["results"][geometry]
    for key,version in (("20","LUT20"),("500","LUT500")):
        expected=results[version]["phi_probabilities"]
        measured=data["p"+key]
        if len(expected)!=36 or any(
                not math.isclose(float(x),float(y),rel_tol=0.,abs_tol=1e-10)
                for x,y in zip(expected,measured)):
            raise ValueError("JSON and CSV angular distributions do not match")
    primaries=float(scope["heed_primaries"])
    collected=float(results["LUT500"]["expected_collected_heed_electrons"])
    if not(0.<collected<=primaries):
        raise ValueError("Unphysical expected electron collection")
    return {
        "generated_parents":chunks*50000,
        "ms0_crossings":int(scope["selected_muon_gas_crossings"]),
        "expected_primaries":primaries,
        "expected_collected_lut500":collected,
    }


def require_pptx():
    try:
        from pptx import Presentation
        from pptx.util import Inches
    except ImportError as exc:
        raise RuntimeError(
            "python-pptx is required for the PPTX output. "
            "On Ubuntu 24.04, use a virtual environment rather than "
            "sudo apt install python3-pptx: "
            "python3 -m venv --system-site-packages "
            "test_runs/mwpc_stage_c/.venv_pi_slides; "
            "test_runs/mwpc_stage_c/.venv_pi_slides/bin/python "
            "-m pip install python-pptx; "
            "then run the script with that venv's Python."
        ) from exc
    return Presentation, Inches


def save_pptx(png,out,presentation,inches):
    ppt=presentation()
    ppt.slide_width=inches(13.333333)
    ppt.slide_height=inches(7.5)
    slide=ppt.slides.add_slide(ppt.slide_layouts[6])
    slide.shapes.add_picture(str(png),0,0,width=ppt.slide_width,
                             height=ppt.slide_height)
    ppt.save(str(out))
    if not out.is_file() or out.stat().st_size==0:
        raise RuntimeError("PowerPoint slide was not saved")



def render_matplotlib(data, png, geometry, stats):
    """Reproduce scientific ROOT-style histogram presentation without ROOT."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
        from matplotlib.ticker import MultipleLocator, MaxNLocator
    except ImportError as exc:
        raise RuntimeError(
            "matplotlib is missing in the slide Python environment. "
            "Install via .venv_pi_slides/bin/python -m pip install matplotlib"
        ) from exc
    edges=np.arange(-180.,181.,10.)
    centers=(edges[:-1]+edges[1:])/2
    max_y=1.22*max(
        max(p+e for p,e in zip(data["p20"],data["se20"])),
        max(p+e for p,e in zip(data["p500"],data["se500"])))
    fig,axes=plt.subplots(2,1,figsize=(13.3333,7.5),dpi=160)
    fig.patch.set_facecolor("white")
    fig.subplots_adjust(left=.10,right=.955,top=.80,bottom=.268,hspace=.48)
    fig.text(.10,.961,
        "Electron-hit azimuth: convergence of the microscopic Garfield LUT",
        fontsize=19,weight="bold",ha="left",va="top")
    label="field wire 75 μm" if geometry=="p3p0_field75p0um" else "no field wire"
    fig.text(.10,.913,
        f"J/ψ  |  MS0  |  pitch 3.0 mm  |  {label}  |  B = 0",
        fontsize=12,color="#444444",ha="left",va="top")
    fig.text(.10,.867,
        f"Statistics: {stats['generated_parents']:,} generated J/ψ events  |  "
        f"{stats['ms0_crossings']:,} MS0 muon gas crossings  |  "
        f"{stats['expected_primaries']/1e6:.2f}M expected Heed electrons  |  "
        f"{stats['expected_collected_lut500']/1e6:.2f}M anode-collected (LUT500)",
        fontsize=10.4,color="#333333",ha="left",va="top")
    for ax,key,count,color,label in (
        (axes[0],"20",20,"#343434","Original LUT"),
        (axes[1],"500",500,"#136999","Updated LUT")):
        p=np.asarray(data["p"+key],dtype=float)
        se=np.asarray(data["se"+key],dtype=float)
        ax.stairs(p,edges,color=color,linewidth=1.65)
        ax.errorbar(centers,p,yerr=se,fmt="o",markersize=2.7,
            color=color,ecolor=color,elinewidth=.9,capsize=1.8,capthick=.8)
        ax.set_xlim(-180,180)
        ax.set_ylim(0,max_y)
        ax.set_xticks(np.arange(-180,181,30))
        ax.xaxis.set_minor_locator(MultipleLocator(10))
        ax.yaxis.set_major_locator(MaxNLocator(5))
        ax.tick_params(axis="both",labelsize=10,direction="in",top=True,right=True)
        ax.tick_params(axis="both",which="minor",direction="in",top=True,right=True)
        for spine in ax.spines.values():
            spine.set_linewidth(.8)
        ax.set_xlabel(r"$\phi_{\mathrm{hit}}$ [degrees]",fontsize=12,labelpad=3)
        ax.set_ylabel("Probability / 10° bin",fontsize=11)
        ax.text(.012,.965,f"{label}: {count} trajectories / spatial cell",
            transform=ax.transAxes,fontsize=12,ha="left",va="top")
    fig.text(.10,.142,
        r"$\sigma_{\mathrm{LUT}}\propto 1/\sqrt{K}$"
        r"$\qquad \sqrt{500/20}=5$"
        "  —  fivefold reduction in expected statistical uncertainty",
        fontsize=14,ha="left",va="center",color="#111111")
    fig.text(.10,.098,
        r"Per spatial cell: $\mathrm{SE}(p_j)\approx"
        r"\sqrt{p_j(1-p_j)/K}$"
        r"; $K$ = independent microscopic Garfield drifts.",
        fontsize=10.5,ha="left",va="center",color="#333333")
    fig.text(.10,.052,
        "Millions of Heed electrons reuse the finite LUT. "
        "Error bars: ±1σ from the weighted Garfield samples, not Heed counting noise.",
        fontsize=9,ha="left",va="center",color="#555555")
    fig.savefig(png,dpi=160,facecolor="white")
    plt.close(fig)
    if not png.is_file() or png.stat().st_size==0:
        raise RuntimeError("PNG plotting failed: "+str(png))

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
    data=validate_csv(a.csv,a.geometry)
    stats=study_statistics(a.csv,a.geometry,data)
    normalizations={"p20":sum(data["p20"]),"p500":sum(data["p500"])}
    print(f"PASS: measured CSV, 36 normalized angular bins, "
          f"geometry={a.geometry}; {normalizations}",flush=True)
    print(f"PASS: exposure {stats['generated_parents']:,} J/ψ parents; "
          f"{stats['ms0_crossings']:,} MS0 crossings; "
          f"{stats['expected_primaries']/1e6:.2f}M expected Heed, "
          f"{stats['expected_collected_lut500']/1e6:.2f}M collected (LUT500)",
          flush=True)
    if a.validate_only:
        return

    presentation,inches=require_pptx()
    a.out.mkdir(parents=True,exist_ok=True)
    prefix=a.geometry+"_slide"
    png=a.out/(prefix+".png")
    pptx=a.out/(prefix+".pptx")
    render_matplotlib(data,png,a.geometry,stats)
    save_pptx(png,pptx,presentation,inches)
    print("PASS: PNG and PPTX slide; no ROOT executable required",flush=True)
    print("OUTPUT",png,flush=True)
    print("OUTPUT",pptx,flush=True)


if __name__=="__main__":
    main()
