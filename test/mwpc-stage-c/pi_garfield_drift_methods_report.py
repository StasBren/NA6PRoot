#!/usr/bin/env python3
"""Read-only summary of small matched-start Garfield three-method pilot."""
import argparse,csv,json,math,statistics
from pathlib import Path
from collections import defaultdict

METHODS=("microscopic","mc","rkf")
GEOMS=("p3_no_fw","p3_field75um")

def read(path):
    data=defaultdict(list)
    with path.open(newline="",encoding="utf-8") as f:
        reader=csv.DictReader(f)
        required={"geometry","method","iu","iv","trial","random_seed",
                  "u_start_mm","v_start_mm","drift_wall_us",
                  "anode_collected","field_wire_collected","other_end",
                  "phi_hit_deg","drift_time_ns"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError("Missing benchmark CSV columns")
        for row in reader:
            if row["geometry"] not in GEOMS or row["method"] not in METHODS:
                raise ValueError("Unexpected method or geometry")
            if sum(int(row[k]) for k in
                ("anode_collected","field_wire_collected","other_end"))!=1:
                raise ValueError("Electron classification is not exclusive")
            data[(row["geometry"],row["method"])].append(row)
    for geom in GEOMS:
        seen=None
        for method in METHODS:
            rows=data[(geom,method)]
            if not rows:raise ValueError("No rows "+str((geom,method)))
            starts={(int(x["iu"]),int(x["iv"]),int(x["trial"]),
                     int(x["random_seed"]),float(x["u_start_mm"]),
                     float(x["v_start_mm"])) for x in rows}
            if len(starts)!=len(rows):
                raise ValueError("Repeated starts for "+method)
            if seen is not None and seen!=starts:
                raise ValueError("Unmatched start coordinates")
            seen=starts
    return data

def stats(rows):
    picked=[r for r in rows if int(r["anode_collected"])]
    wall=[float(r["drift_wall_us"]) for r in rows]
    if not all(math.isfinite(t) and t>=0 for t in wall):
        raise ValueError("Invalid drift timing")
    drift=[float(r["drift_time_ns"]) for r in picked]
    if not all(math.isfinite(t) and t>=0 for t in drift):
        raise ValueError("Bad anode drift time")
    phi=[0]*12
    for row in picked:
        deg=float(row["phi_hit_deg"])
        if not math.isfinite(deg) or not -180<=deg<=180:
            raise ValueError("Bad collected impact azimuth")
        phi[max(0,min(11,int((deg+180.)/30.)))]+=1
    attachment_status=sum(int(r.get("status",0))==-7 for r in rows)
    return {"starts":len(rows),"anode":len(picked),
        "status_minus7_attachment_fraction":attachment_status/len(rows),
        "mc_attachment_setting":next(iter({r.get("mc_attachment","unrecorded") for r in rows}))
            if len({r.get("mc_attachment","unrecorded") for r in rows})==1 else "mixed",
        "anode_fraction":len(picked)/len(rows),
        "fieldwire_fraction":sum(int(r["field_wire_collected"]) for r in rows)/len(rows),
        "mean_drift_ns":statistics.mean(drift) if drift else None,
        "rms_drift_ns":statistics.pstdev(drift) if drift else None,
        "transport_wall_us_total":sum(wall),
        "mean_us_per_start":statistics.mean(wall),
        "median_us_per_start":statistics.median(wall),
        "phi_probability_per_30deg":[x/max(1,len(picked)) for x in phi]}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("csv",type=Path)
    p.add_argument("--out",type=Path,
      default=Path("test_runs/mwpc_stage_c/pi_garfield_drift_methods_v1"))
    a=p.parse_args()
    data=read(a.csv)
    result={}
    for geom in GEOMS:
        vals={m:stats(data[(geom,m)]) for m in METHODS}
        ref=vals["microscopic"]
        for m in ("mc","rkf"):
            x=vals[m]
            x["speedup_transport_only"]=(
                ref["transport_wall_us_total"]/x["transport_wall_us_total"]
                if x["transport_wall_us_total"]>0 else None)
            x["phi_TV_vs_micro_30deg"]=.5*sum(abs(v-w) for v,w in zip(
                x["phi_probability_per_30deg"],ref["phi_probability_per_30deg"]))
        result[geom]=vals
        print("\n"+geom+"; N="+str(ref["starts"])+" matched starts")
        for m,v in vals.items():
            print(f"  {m:12s} anode={v['anode_fraction']:.3f} "+
                f"attachment(-7)={v['status_minus7_attachment_fraction']:.3f} "+
                f"drift mean={v['mean_drift_ns']} ns "+
                f"CPU={v['mean_us_per_start']:.1f} us/e "+
                f"speedup={v.get('speedup_transport_only',1.)}")
    a.out.mkdir(parents=True,exist_ok=True)
    (a.out/"method_comparison.json").write_text(
        json.dumps({"source":str(a.csv),"geometries":result,
          "warning":"Only an exploratory 25-box sparse matched-start pilot. "+
          "MC and RKF require precomputed Magboltz coefficients. "+
          "30-degree histograms are low-statistical-precision. "+
          "Speed ratios exclude gas-table preparation and executable startup. "+
          "Test MC step convergence and gas-table E-grid before large LUT."},
          indent=2)+"\n",encoding="utf-8")
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("NOTE: matplotlib unavailable; JSON summary is complete")
    else:
        fig,axes=plt.subplots(1,2,figsize=(12.4,4.4),sharey=True)
        phi_edges=list(range(-180,181,30))
        for ax,geom in zip(axes,GEOMS):
            vals=result[geom]
            for method,color in (("microscopic","#404040"),
                                 ("mc","#127ca4"),("rkf","#b9783b")):
                ax.stairs(vals[method]["phi_probability_per_30deg"],
                          phi_edges,color=color,linewidth=1.5,label=method)
            ax.set_title(geom)
            ax.set_xlabel("Endpoint azimuth [degrees]")
            ax.legend(frameon=False)
        axes[0].set_ylabel("Probability / 30 degree bin")
        fig.suptitle("Garfield drift-method pilot: identical start sites (sparse)")
        fig.tight_layout()
        fig.savefig(a.out/"method_phi_distributions.png",dpi=150)
        plt.close(fig)
    print("PASS:",a.out/"method_comparison.json")
if __name__=="__main__":
    main()
