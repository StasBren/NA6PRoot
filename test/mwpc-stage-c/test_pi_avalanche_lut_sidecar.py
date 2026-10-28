#!/usr/bin/env python3
"""Offline exact-start avalanche sidecar CSV contract and conservation tests.

No Geant4, Garfield, frozen LUT reads or physical output are modified.
"""
from __future__ import annotations
import csv
import tempfile
from pathlib import Path

from pi_run_avalanche_lut_sidecar import (
    COLUMNS, BIRTH_COLUMNS, TIME_COLUMNS, validate_outputs, choose_rows,
)
from pi_audit_avalanche_lut_reuse import CONFIGS


def write_csv(path,columns,rows):
    with path.open("w",newline="") as f:
        wr=csv.DictWriter(f,fieldnames=columns)
        wr.writeheader()
        wr.writerows(rows)


def main():
    # Generate the same 3 x 6 x 3 synthetic per-geometry selected keys.
    sources=[]
    for k in (20,500):
        for geom,pitch,fw in CONFIGS:
            for iu in (0,round(pitch/.1)//2,round(pitch/.1)-1):
                for iv in (0,7,14,15,22,29):
                    for j in (0,k//2,k-1):
                        sources.append({"lut_trials":str(k),"geometry":geom,
                            "iu":str(iu),"iv":str(iv),"sample":str(j),
                            "box_id":str(iu*30+iv)})
    assert len(sources)==432
    sample8=choose_rows(sources,"smoke8")
    assert len(sample8)==8 and len({(x["lut_trials"],x["geometry"]) for x in sample8})==8
    assert choose_rows(sources,"full432")==sources

    origin={"lut_trials":"20","geometry":"p3_noFW","box_id":"465",
            "sample":"10","original_lut_sha256":"a"*64,
            "random_seed":"9023","u_start_mm":"0.012345678901",
            "v_start_mm":"0.123456789012","w_start_mm":"0",
            "old_collected":"1","old_field_wire_collected":"0"}
    crossing="LUT20|p3_noFW|465|10"
    trial={column:"" for column in COLUMNS}
    trial.update(lut_trials="20",original_geometry="p3_noFW",
        original_lut_sha256="a"*64,box_id="465",iu="15",iv="15",
        sample="10",original_random_seed="9023",
        original_u_start_mm=origin["u_start_mm"],
        original_v_start_mm=origin["v_start_mm"],
        original_w_start_mm="0",old_collected="1",
        old_field_wire_collected="0",crossing_id=crossing,
        configuration="p3_noFW",primary_index="0",
        garfield_seed="299223",seed_u_mm=origin["u_start_mm"],
        seed_v_mm=origin["v_start_mm"],seed_w_mm="0",
        garfield_ne="5",garfield_ni="3",ion_births="3",
        collected_e="5",successful="1",capped="0",ion_up="2",
        ion_down="1",ion_top="1",ion_bottom="1")
    timer={column:"" for column in TIME_COLUMNS}
    timer.update(crossing_id=crossing,configuration="p3_noFW",
       wall_seconds="0.09",garfield_seed="299223")
    births=[]
    for _ in range(3):
        b={column:"" for column in BIRTH_COLUMNS}
        b.update(crossing_id=crossing,seed_ion_births_total="3",
                 reservoir_weight="1.0")
        births.append(b)

    with tempfile.TemporaryDirectory() as tmp:
        root=Path(tmp)
        def write():
            write_csv(root/"lut_avalanche_trials.csv",COLUMNS,[trial])
            write_csv(root/"lut_avalanche_birth_reservoir.csv",BIRTH_COLUMNS,births)
            write_csv(root/"lut_avalanche_timing.csv",TIME_COLUMNS,[timer])
        write()
        report=validate_outputs(root,[origin])
        assert report["n_trials"]==report["n_timed"]==1
        assert report["n_reservoir_ion_positions"]==3
        assert report["total_ion_births"]==3
        assert abs(report["mean_seed_wall_seconds"]-.09)<1.e-10

        trial["seed_u_mm"]="0.7"
        write()
        try:validate_outputs(root,[origin])
        except ValueError:pass
        else:raise AssertionError("Original exact-start position shift not caught")
        trial["seed_u_mm"]=origin["u_start_mm"]

        births.pop()
        write()
        try:validate_outputs(root,[origin])
        except ValueError:pass
        else:raise AssertionError("Incomplete ion birth reservoir not caught")
    print("PASS avalanche exact-start sidecar: 8 geometry/LUT selections, original coordinate, timing and ion-reservoir conservation")


if __name__=="__main__":main()
