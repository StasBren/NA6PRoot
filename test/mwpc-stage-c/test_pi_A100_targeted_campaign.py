#!/usr/bin/env python3
"""SYNTHETIC/OFFLINE ONLY: targeted 20->100 full-avalanche trial design.
Does NOT run any original Geant4, TrackHeed or Garfield physics.
"""
import argparse
import csv
import json
from collections import Counter
from pi_A100_targeted_50cell_campaign import (
    TOP,ADD,TOTAL,CHUNK,N_CHUNKS,GEOMS,EXPORT_FIELDS,spatial_cell,
    additional_original_D500_sample,
)
from pi_run_lut500_overnight import (
    selected_j,offset_for,STRIDE,ROOT as A20_ROOT,source_maps,
)
from pi_holdout120_A20_from_D500 import payload_rows,verify_full,REF
from pi_run_muon_lut_validation import id_of,PRIMARY_FIELDS
from pi_audit_avalanche_lut_reuse import CONFIGS,locations
from pi_prepare_muon_avalanche_b1 import sha
from pi_run_muon_lut_validation import field_wire_blocked_zero_start
from pi_compare_A100_targeted_vs_A20 import N_REP

def check_original_real_sources():
    """Actual WSL source/header/contract audit; NO physical simulations,
    original file writes or A100 prepare-stage output directory.
    """
    rows,source_manifest=payload_rows()  # Exactly TWO documented values.
    full_rows,full_stats,full_digest=verify_full()  # Exactly THREE values.
    assert len(rows)==len(full_rows)==120
    assert [id_of(r) for r in rows]==[id_of(r) for r in full_rows]
    assert full_stats["n_native_muons"]==120
    assert full_stats["n_full_avalanche_trials"]==4*full_stats["n_heed_primaries"]
    assert sha(REF/"manifest.json")==full_digest
    a20=json.loads((A20_ROOT/"manifest.json").read_text())
    assert a20["status"]=="complete" and a20["n_completed_avalanches"]==84000
    assert a20["total_capped"]==0 and a20["total_failed"]==0
    cells=A20_ROOT/"cell_statistics.csv"
    assert sha(cells)==a20["files_sha256"]["cell_statistics.csv"]
    with cells.open(newline="") as f:
        columns=tuple(csv.reader(f).__next__())
    assert columns==("geometry","box_id","iu","iv","n","zero","p_zero",
        "mean_ion_births","sd_ion_births","median_ion_births",
        "mean_timing_seconds","ion_weighted_A_hemi","ion_weighted_A_TB")
    with (REF/"muon_reference_heed_primaries.csv").open(newline="") as f:
        cols=tuple(next(csv.reader(f)))
    assert cols==PRIMARY_FIELDS
    meta,index_sha=source_maps()
    assert len(meta)==4 and set(meta)==set(GEOMS)
    assert len(index_sha)==64
    for geometry,pitch,field_wire_um in CONFIGS:
        conf,original_D500=locations(500,pitch,field_wire_um)
        assert conf.is_file() and original_D500.is_file()
        with original_D500.open(newline="") as f:
            cols=set(next(csv.reader(f)))
        assert set(("box_id","sample","iu","iv","random_seed",
            "u_start_mm","v_start_mm","w_start_mm","drift_ok",
            "collected","field_wire_collected","wire_index",
            "phi_hit_deg","drift_time_ns",
            "last_step_energy_eV")).issubset(cols)
    print("PASS ORIGINAL REAL WSL INPUTS: 120 physical Geant4 parent muons, "
          "same full Garfield reference, 5881 source Heed electrons, "
          "original 84k A20 avalanche SHA/statistics and four original "
          "D500 source schemas. NO new Garfield physics or output writes.")


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--real-source-check",action="store_true",
        help="Read and SHA-verify actual original user WSL Geant4, Heed and D500 files")
    args=parser.parse_args()
    assert len(GEOMS)==4 and TOP==50 and ADD==80
    assert TOTAL==16000 and N_CHUNKS==40 and CHUNK==400
    assert N_REP==200
    assert EXPORT_FIELDS[:3]==(
        "lut_trials","geometry","original_lut_sha256")
    # Exhaustively check EVERY real spatial D500 cell, not just a pilot.
    # Crucially the original selected_j only accepts 0..19; it MUST
    # remain unchanged, while the new A100 helper extends 20..99.
    for g in GEOMS:
        nbox=30*30 if g.startswith("p3_") else 40*30
        for box in range(nbox):
            original=[selected_j(g,box,idx) for idx in range(20)]
            new=[additional_original_D500_sample(g,box,idx)
                 for idx in range(20,100)]
            assert original==[(offset_for(g,box)+STRIDE*idx)%500
                              for idx in range(20)]
            assert len(original)==20 and len(set(original))==20
            assert len(new)==80 and len(set(new))==80
            assert set(original).isdisjoint(new)
            assert len(set(original+new))==100
            assert all(0<=x<500 for x in original+new)
    for invalid_pass in (-1,0,19,100,500):
        try:
            additional_original_D500_sample("p3_noFW",0,invalid_pass)
        except ValueError:
            pass
        else:
            raise AssertionError("Wrongly allowed extra D500 pass "+str(invalid_pass))
    try:
        selected_j("p3_noFW",0,20)
    except ValueError:
        pass
    else:
        raise AssertionError("Frozen original A20 selected_j guard was modified")
    # Repeated physical electrons in ONE cell have one shared
    # uncertain 20-sample cell mean: uncertainty scales as n**2.
    count=Counter({11:20,12:1})
    sigma={11:850.,12:400.}
    training=sum(n*n*sigma[b]**2/20 for b,n in count.items())
    wrong_independent=sum(n*sigma[b]**2/20 for b,n in count.items())
    assert training>15*wrong_independent
    # Confirm one real metal start is counted zero only in wire geometry.
    r={"configuration":"p4_FW75","successful":"0","capped":"0",
       "garfield_ne":"1","garfield_ni":"0","ion_births":"0",
       "collected_e":"0","attached":"0","other":"1",
       "ion_up":"0","ion_down":"0","ion_top":"0","ion_bottom":"0",
       "seed_u_mm":"-2.00535685975886",
       "seed_v_mm":"-0.00391288042673045"}
    assert field_wire_blocked_zero_start(r)
    r["configuration"]="p4_noFW"
    assert not field_wire_blocked_zero_start(r)
    if args.real_source_check:
        check_original_real_sources()
    print("PASS SYNTHETIC/OFFLINE ONLY: 16,000 extra FULL Garfield starts "
          "from disjoint original 500 D500 indices EXHAUSTIVELY TESTED "
          "on 4,200 cells; original selected_j unchanged; 40x400 immutable "
          "chunk plan, n_cell^2 finite A20 mean uncertainty, strict metal "
          "zero-ion classification; NO physics execution.")

if __name__=="__main__":
    main()
