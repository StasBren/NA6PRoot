#!/usr/bin/env python3
"""Prove parent-level grouping, detector conservation, and electron totals."""
from __future__ import annotations
import csv
import gzip
import tempfile
from pathlib import Path
from pi_parent_asymmetry_geometry_dashboard import aggregate_lut20,RAW_FIELDS,GEOMS
from pi_pitch_fieldwire_pilot import tag


def test_parent_vs_crossing():
    with tempfile.TemporaryDirectory(prefix="parent_geometry_fixture_") as tmp:
        csvfile=Path(tmp)/"tiny.csv.gz"
        old={"configs":[{"key":tag(3.,0.)},{"key":tag(3.,75.)}],
             "stations":{}}
        summary={"stations":{}}
        for st in range(6):
            d0=str(st*10+1);d1=str(st*10+2)
            old["stations"][str(st)]={
                "gas_crossings":2,"heed_primaries":6,
                "layout":[{"id":int(d0),"row":0,"col":0,"layer":0,
                          "x":0,"y":0,"hx":1,"hy":1},
                         {"id":int(d1),"row":0,"col":1,"layer":0,
                          "x":2,"y":0,"hx":1,"hy":1}],
                "n_physical_chambers":2,
                "station_parent":{"configurations":{
                    tag(3.,0.):{"parents_with_gas_crossing":1},
                    tag(3.,75.):{"parents_with_gas_crossing":1}}},
                "chambers":{d0:{},d1:{}}}
            summary["stations"][str(st)]={"gas_crossings":2,"configs":{
                cfg:{"n_upper":3,"n_lower":3} for _,_,cfg in GEOMS}}
        with gzip.open(csvfile,"wt",encoding="utf-8",newline="") as f:
            wr=csv.DictWriter(f,fieldnames=RAW_FIELDS)
            wr.writeheader()
            for st in range(6):
                for crossing in range(2):
                    for pitch,fw,cfg in GEOMS:
                        up,down=((3,0) if crossing==0 else (0,3))
                        wr.writerow({
                            "source":"jpsi40-bore","chunk":1,"event":0,
                            "hit_index":st*2+crossing,"track_id":1+crossing,
                            "station":st,"det_id":st*10+1+crossing,
                            "pitch_mm":pitch,"field_wire_um":fw,
                            "lut_trials_per_cell":20,
                            "n_primary":3,"n_anode_finite":3,
                            "up":up,"down":down,"hemisphere_boundary":0})
        result=aggregate_lut20(csvfile,old,summary)
        for pitch in (3.,4.):
            for st in range(6):
                s=result[pitch]["stations"][str(st)]
                assert s["gas_crossings"]==2
                for fw in (0.,75.):
                    cfg=tag(pitch,fw)
                    parent=s["station_parent"]["configurations"][cfg]
                    assert parent["parents_with_gas_crossing"]==1
                    assert parent["valid_A"]==1
                    assert parent["mean_A"]==0.
                    assert parent["total_up"]==3 and parent["total_down"]==3
                    assert parent["hist"][20]==1
                    # DETs are individually fully asymmetric, NOT additive A.
                    a=s["chambers"][str(st*10+1)]["configurations"][cfg]
                    b=s["chambers"][str(st*10+2)]["configurations"][cfg]
                    assert a["mean_A"]==1 and b["mean_A"]==-1
                    assert a["total_up"]+b["total_up"]==parent["total_up"]
                    assert a["total_down"]+b["total_down"]==parent["total_down"]
    print("PASS parent-level 3/4mm LUT20, same-parent crossing aggregation, "
          "nonlinear A, physical DET electron conservation")


if __name__=="__main__":test_parent_vs_crossing()
