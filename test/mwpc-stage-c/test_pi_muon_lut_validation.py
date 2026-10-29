#!/usr/bin/env python3
"""Offline real-muon and conditional-LUT identity mathematics tests.

Does NOT execute Geant4, Garfield++, rclone or any real detector physics;
the actual B1 build and check20 are mandatory after this smoke.
"""
from __future__ import annotations

import csv
import math
import tempfile
from pathlib import Path

from pi_run_muon_lut_validation import (
    GEOMS,N_U,PRIMARY_FIELDS,CLUSTER_FIELDS,SEED_FIELDS,LUT_FIELDS,
    spatial_cell,entry_phase_mm,round_away,stable_seed,percent,
    asym,summarize_samples,read_csv,field_wire_blocked_zero_start,
)


def main():
    assert tuple(GEOMS)==tuple(N_U)
    assert len(LUT_FIELDS)==39
    assert round_away(.5)==1 and round_away(-.5)==-1
    assert entry_phase_mm(0,0,3.)==1.5
    assert entry_phase_mm(0,0,4.)==0.
    assert entry_phase_mm(0,1,3.)==0.
    assert entry_phase_mm(0,1,4.)==2.
    assert entry_phase_mm(0,1,3.)==0.
    assert spatial_cell("p4_noFW",0,0.,0.,0.)==(20*30+15,0.,0.)
    assert spatial_cell("p3_FW75",1,0.,0.,0.)==(15*30+15,0.,0.)
    assert spatial_cell("p3_noFW",0,0.,0.,0.)==(15,-1.5,0.)
    assert spatial_cell("p4_FW75",1,0.,0.,0.)==(15,-2.,0.)
    assert spatial_cell("p3_noFW",1,0.,-.015,.01)[0]==13*30+15
    assert spatial_cell("p3_noFW",1,0.,+.015,-.01)[0]==16*30+14
    assert spatial_cell("p4_noFW",0,0.,+.015,+.01)[0]==21*30+15
    assert stable_seed("muon","p3",0)==stable_seed("muon","p3",0)
    assert stable_seed("muon","p3",0)!=stable_seed("muon","p4",0)
    assert asym(0,0) is None
    assert asym(3,1)==.5
    assert math.isclose(percent([0,10,20],.5),10.)
    assert math.isclose(percent([0,10],.5),5.)
    a=summarize_samples([0,100,200,300])
    assert a["mean"]==150 and a["p50"]==150 and a["sd"]>0
    # Pure input parser fails on extra or shuffled columns, because the
    # cached native Heed electrons MUST have exact source provenance.
    with tempfile.TemporaryDirectory() as tmp:
        p=Path(tmp)/"primaries.csv"
        with p.open("w",newline="") as f:
            w=csv.DictWriter(f,fieldnames=PRIMARY_FIELDS)
            w.writeheader()
            row={k:"0" for k in PRIMARY_FIELDS}
            row.update(crossing_id="real|1|4|2|1|0|17",
                       primary_index="0",cluster_index="0")
            w.writerow(row)
        data=read_csv(p,PRIMARY_FIELDS)
        assert len(data)==1
        assert data[0]["cluster_index"]=="0"
        try:read_csv(p,CLUSTER_FIELDS)
        except ValueError:pass
        else:raise AssertionError("Different Heed cluster/electron schema accepted")
    # A real Heed conduction electron embedded in the field-wire metal
    # is a PHYSICAL zero-ion outcome, despite Garfield successful=0.
    # Field-wire AddWire diameter=75um -> metal radius=0.0375mm.
    # Existing original reference must remain unchanged; classifier only.
    wire={k:"0" for k in SEED_FIELDS}
    wire.update(configuration="p4_FW75",successful="0",capped="0",
                garfield_ne="1",garfield_ni="0",ion_births="0",
                seed_u_mm="-2.00535685975886",
                seed_v_mm="-0.00391288042673045")
    assert field_wire_blocked_zero_start(wire)
    # Exact C++ 'other' endpoint semantics: a failed initial electron
    # inside metal can produce one non-collected/non-attached endpoint.
    wire["other"]="1"
    assert field_wire_blocked_zero_start(wire)
    wire["other"]="2"
    assert not field_wire_blocked_zero_start(wire)
    wire["other"]="1"
    for x,y in [(-2.00920622158948,.0229186340901303),
                (-2.00612397367103,.0239032733895584)]:
        wire.update(seed_u_mm=str(x),seed_v_mm=str(y))
        assert field_wire_blocked_zero_start(wire)
    wire.update(seed_u_mm="-2.06",seed_v_mm="0")
    assert not field_wire_blocked_zero_start(wire)
    wire.update(seed_u_mm="-2.005",seed_v_mm="0")
    for geom in ("p4_noFW","p3_noFW"):
        wire["configuration"]=geom
        assert not field_wire_blocked_zero_start(wire)
    wire["configuration"]="p4_FW75"
    for k,v in (("ion_births","1"),("capped","1"),("garfield_ne","2"),
                ("garfield_ni","1"),("successful","1"),
                ("collected_e","1"),("attached","1")):
        old=wire[k];wire[k]=v
        assert not field_wire_blocked_zero_start(wire)
        wire[k]=old
    wire["configuration"]="p3_FW75"
    wire.update(seed_u_mm="-1.505",seed_v_mm="0")
    assert field_wire_blocked_zero_start(wire)
    wire["configuration"]="p4_FW75"
    wire.update(seed_u_mm="-2.005",seed_v_mm=".0375")
    assert not field_wire_blocked_zero_start(wire)
    print("PASS 20 native muon×4 geometry LUT-grid phase, spatial lookup, "
          "strict shared-Heed identity schema, deterministic independent "
          "replicates, and undefined asymmetries (OFFLINE only)")


if __name__=="__main__":
    main()
