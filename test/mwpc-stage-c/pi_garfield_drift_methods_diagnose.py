#!/usr/bin/env python3
"""Diagnose WHY the small Garfield AvalancheMC run misses anode collection.

READ-ONLY on matched_5um.csv / matched_2um.csv. Unlike conventional
collection reports, this breaks down ALL rejected endpoints by Garfield
termination status and the distance to each physical conductor, to
distinguish actual losses from method-dependent endpoint conventions.

Endpoint distances to wires are diagnostics, not proof of impact:
status=-5 and status=-1 can mean different geometric endings across methods.
Never retroactively relabel an electron as collected just from proximity.
"""
from __future__ import annotations

import argparse
import json
import math
from collections import Counter
from pathlib import Path

from pi_garfield_drift_methods_report import read,GEOMS,METHODS

ANODE_RADIUS_MM=0.015
FW_RADIUS_MM=0.0375
PITCH_MM=3.0
HALF_WIRES=4


def nearest_distances(row,field):
    """Distance to nearest finite conductor (mm), endpoint may be absent."""
    u=float(row["u_end_mm"])
    v=float(row["v_end_mm"])
    if not math.isfinite(u) or not math.isfinite(v):
        return None,None,None
    anode=min(math.hypot(u-i*PITCH_MM,v)
              for i in range(-HALF_WIRES,HALF_WIRES+1))
    fw=(min(math.hypot(u-(i+0.5)*PITCH_MM,v)
            for i in range(-HALF_WIRES,HALF_WIRES))
        if field else None)
    cath=min(abs(v+3.0),abs(v-3.0))
    return anode,fw,cath


def diagnose(csvpath):
    groups=read(csvpath)
    outcome={}
    for geometry in GEOMS:
        result={}
        field=geometry=="p3_field75um"
        micro={(
            r["iu"],r["iv"],r["trial"],r["random_seed"]
        ):r for r in groups[(geometry,"microscopic")]}
        for method in METHODS:
            rows=groups[(geometry,method)]
            status_all=Counter()
            status_other=Counter()
            status_lost_micro_collected=Counter()
            distances=Counter()
            by_flag=Counter()
            examples=[]
            for row in rows:
                status=int(row["status"])
                status_all[status]+=1
                a=int(row["anode_collected"])
                f=int(row["field_wire_collected"])
                o=int(row["other_end"])
                by_flag["anode"]+=a
                by_flag["field_wire"]+=f
                by_flag["other"]+=o
                near_anode,near_fw,near_cath=nearest_distances(row,field)
                if o:
                    status_other[status]+=1
                    if near_anode is None:
                        distances["missing/nonfinite endpoint"]+=1
                    elif near_anode<=1.5*ANODE_RADIUS_MM:
                        distances["within 1.5 anode radii (but rejected)"]+=1
                    elif near_anode<=3*ANODE_RADIUS_MM:
                        distances["anode 1.5-3 radii"]+=1
                    elif near_anode<=6*ANODE_RADIUS_MM:
                        distances["anode 3-6 radii"]+=1
                    elif field and near_fw is not None and near_fw<=3*FW_RADIUS_MM:
                        distances["within 3 field-wire radii"]+=1
                    elif near_cath is not None and near_cath<=0.05:
                        distances["within 50 um of cathode"]+=1
                    else:
                        distances["elsewhere in gas/outside"]+=1
                    k=(row["iu"],row["iv"],row["trial"],row["random_seed"])
                    if method!="microscopic" and int(micro[k]["anode_collected"]):
                        status_lost_micro_collected[status]+=1
                        if len(examples)<5:
                            examples.append(dict(
                                iu=int(row["iu"]),iv=int(row["iv"]),
                                trial=int(row["trial"]),
                                status=status,
                                anode_distance_um=(near_anode*1000
                                    if near_anode is not None else None),
                                field_wire_distance_um=(near_fw*1000
                                    if near_fw is not None else None),
                                cathode_distance_um=(near_cath*1000
                                    if near_cath is not None else None),
                                u_end_mm=float(row["u_end_mm"]),
                                v_end_mm=float(row["v_end_mm"])))
            result[method]=dict(
                n=len(rows),
                collection=dict(by_flag),
                status_all={str(k):v for k,v in sorted(status_all.items())},
                status_rejected={str(k):v for k,v in sorted(status_other.items())},
                location_rejected=dict(sorted(distances.items())),
                rejected_but_microscopic_collected_by_status={
                    str(k):v for k,v in sorted(status_lost_micro_collected.items())
                },
                five_example_rejections=examples,
            )
            print("\n"+geometry+" / "+method+": N="+str(len(rows)))
            print("  collection:",dict(by_flag))
            print("  statuses (all):",dict(status_all))
            print("  statuses (other):",dict(status_other))
            print("  rejected endpoint location:",dict(distances))
            if method!="microscopic":
                print("  failures where SAME birth site microscopic collected:",
                      dict(status_lost_micro_collected))
        outcome[geometry]=result
    return {"source":str(csvpath),"geometries":outcome,
      "warning":"Garfield status is not alone a hit label: inspect endpoint "
      "and conductor geometry together; this diagnostic does not "
      "reclassify any simulated electron."}


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("csv",type=Path)
    ap.add_argument("--out",type=Path,default=None)
    args=ap.parse_args()
    result=diagnose(args.csv)
    if args.out is not None:
        if args.out.exists():
            raise FileExistsError("Refusing to overwrite old diagnosis: "+str(args.out))
        args.out.parent.mkdir(parents=True,exist_ok=True)
        args.out.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
        print("PASS diagnostic file:",args.out)


if __name__=="__main__":
    main()
