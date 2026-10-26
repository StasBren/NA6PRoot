#!/usr/bin/env python3
"""Read-only report for Magboltz eta(E) along the actual MWPC field map.

Does not simulate drift, and deliberately does NOT infer total electron
survival from pointwise eta: survival depends on the integrated drift path.
"""
from __future__ import annotations
import argparse,csv,json,math,statistics
from pathlib import Path
from collections import defaultdict

SCHEMA={"kind","geometry","u_mm","v_mm","E_V_cm","transport_ok",
        "eta_per_cm","alpha_per_cm","drift_cm_per_ns",
        "diff_long_sqrt_cm","diff_trans_sqrt_cm"}
E_MIN,E_MAX=100.,500000.

def read(path):
    groups=defaultdict(list)
    with path.open(newline="",encoding="utf-8") as f:
        reader=csv.DictReader(f)
        if not SCHEMA.issubset(reader.fieldnames or []):
            raise ValueError("Missing Garfield attachment CSV columns")
        for r in reader:
            if r["kind"] not in ("gas_sweep","chamber"):
                raise ValueError("Unexpected kind "+r["kind"])
            if int(r["transport_ok"]) not in (0,1):
                raise ValueError("Bad transport status")
            row={k:(float(r[k]) if k in (
                "u_mm","v_mm","E_V_cm","eta_per_cm","alpha_per_cm",
                "drift_cm_per_ns","diff_long_sqrt_cm",
                "diff_trans_sqrt_cm") else r[k]) for k in SCHEMA}
            row["transport_ok"]=int(r["transport_ok"])
            if not math.isfinite(row["E_V_cm"]) or row["E_V_cm"]<0:
                raise ValueError("Invalid field magnitude")
            groups[(r["kind"],r["geometry"])].append(row)
    if len(groups[("gas_sweep","n_a_n")])<20:
        raise ValueError("Incomplete eta(E) sweep")
    for name in ("p3_no_fw","p3_field75um"):
        if len(groups[("chamber",name)])<10:
            raise ValueError("Too few sampled chamber field locations: "+name)
    return groups

def summarise(path,out):
    groups=read(path)
    result={"source":str(path),"units":"eta,alpha: 1/cm; E: V/cm",
            "gas_table_E_range_V_cm":[E_MIN,E_MAX],
            "samples":{}}
    for (kind,geometry),rows in sorted(groups.items()):
        valid=[r for r in rows if r["transport_ok"]]
        e=[r["E_V_cm"] for r in valid]
        eta=[r["eta_per_cm"] for r in valid]
        result["samples"][kind+"/"+geometry]={
            "total":len(rows),"valid":len(valid),
            "E_min":min(e) if e else None,"E_max":max(e) if e else None,
            "eta_min":min(eta) if eta else None,
            "eta_max":max(eta) if eta else None,
            "eta_negative":sum(x<0 for x in eta),
            "eta_gt_1_per_cm":sum(x>1 for x in eta),
            "E_below_grid":sum(x<E_MIN for x in e),
            "E_above_grid":sum(x>E_MAX for x in e)}
    out.mkdir(parents=True,exist_ok=True)
    target=out/"attachment_field_summary.json"
    target.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    for k,v in result["samples"].items():
        print(k,v)
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("No matplotlib; numerical JSON complete")
        return
    fig,axes=plt.subplots(2,1,figsize=(9.8,7.8))
    sweep=[r for r in groups[("gas_sweep","n_a_n")] if r["transport_ok"]]
    axes[0].plot([r["E_V_cm"] for r in sweep],
                 [r["eta_per_cm"] for r in sweep],color="#343434")
    axes[0].set_xscale("log")
    axes[0].set_yscale("symlog",linthresh=1e-3)
    axes[0].set_ylabel("Attachment coefficient eta [1/cm]")
    axes[0].set_title("Loaded Magboltz attachment coefficient vs electric field")
    axes[0].grid(alpha=.2)
    for geom,color in (("p3_no_fw","#1d77a1"),
                       ("p3_field75um","#bf784b")):
        rows=[r for r in groups[("chamber",geom)] if r["transport_ok"]]
        axes[1].scatter([r["E_V_cm"] for r in rows],
                        [r["eta_per_cm"] for r in rows],
                        label=geom,color=color,s=18,alpha=.7)
    axes[1].set_xscale("log")
    axes[1].set_yscale("symlog",linthresh=1e-3)
    axes[1].set_xlabel("Actual MWPC electric field [V/cm]")
    axes[1].set_ylabel("Tabulated eta at sampled gas locations [1/cm]")
    axes[1].legend(frameon=False)
    axes[1].grid(alpha=.2)
    fig.text(.05,.012,"Local eta is not survival probability; survival depends on integral of eta along trajectories.",fontsize=9)
    fig.tight_layout(rect=(0,.03,1,1))
    fig.savefig(out/"attachment_eta_vs_field.png",dpi=150)
    plt.close(fig)
    print("PASS",target,out/"attachment_eta_vs_field.png")

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("csv",type=Path)
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    summarise(args.csv,args.out)

if __name__=="__main__":
    main()
