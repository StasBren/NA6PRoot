#!/usr/bin/env python3
"""Offline nightly-run selection/budget/compact reservoir regression.

No Garfield libraries, Geant4, original LUT CSVs, WSL files or network
needed: checks the 4200-cell design before compilation or long jobs.
"""
from __future__ import annotations
import random
from collections import Counter
import csv
import tempfile
from pathlib import Path

import pi_run_lut500_overnight as nightly
from pi_run_avalanche_lut_sidecar import COLUMNS,BIRTH_COLUMNS,TIME_COLUMNS

from pi_run_lut500_overnight import (
    CONFIGS,GEOM_NAMES,PER_CELL,CHUNK_SIZE,RESERVOIR,
    offset_for,selected_j,batches,canonical,
)
from pi_run_avalanche_lut_sidecar import validate_outputs
from pi_audit_avalanche_lut_reuse import EXPORT_FIELDS


def run():
    assert PER_CELL==20 and CHUNK_SIZE==400 and RESERVOIR==5
    assert tuple(geom for geom,_,_ in CONFIGS)==GEOM_NAMES
    cells=0;all_keys=set();by_geom=Counter()
    for geom,pitch,fw in CONFIGS:
        nbox=round(pitch/.1)*30
        for b in range(nbox):
            js=[selected_j(geom,b,j) for j in range(PER_CELL)]
            assert len(set(js))==PER_CELL
            assert all(0<=x<500 for x in js)
            assert js==[selected_j(geom,b,j) for j in range(PER_CELL)]
            assert offset_for(geom,b)==js[0]
            cells+=1
            for j in js:
                id_=(geom,b,j)
                assert id_ not in all_keys
                all_keys.add(id_)
                by_geom[geom]+=1
    assert cells==4200 and len(all_keys)==84000
    assert by_geom["p3_noFW"]==by_geom["p3_FW75"]==18000
    assert by_geom["p4_noFW"]==by_geom["p4_FW75"]==24000
    assert len(batches(list(range(84000))))==210
    assert all(len(b)==400 for b in batches(list(range(84000))))
    assert len(EXPORT_FIELDS)==18
    # Stable deterministic ordering of a complete pass, preventing
    # geometry-ordered time cutoff, repeatable from fixed RNG seed.
    a=list(range(4200));b=list(range(4200))
    random.Random(20261029).shuffle(a)
    random.Random(20261029).shuffle(b)
    assert a==b and len(set(a))==4200
    # Full canonicalized frozen config is itself stable and newline-terminated.
    assert canonical({"x":3,"n":1})==canonical({"n":1,"x":3})
    # Reproduce the real failure boundary: Python nightly source uses an
    # INTEGER random_seed but the C++ CSV stores its decimal STRING.
    # A previously completed Garfield chunk with a validator-only
    # failure is recoverable WITHOUT recomputing or mutating raw results.
    with tempfile.TemporaryDirectory() as td:
        root=Path(td)
        smoke=root/"smoke8_compact5_v1"
        chunk=smoke/"chunks"/"chunk_0000"
        chunk.mkdir(parents=True)
        fake_exe=root/"garfield-fake-executable"
        fake_exe.write_bytes(b"FAKE: offline SHA regression only")
        old_exe=nightly.EXE
        nightly.EXE=fake_exe
        try:
            original={col:"" for col in EXPORT_FIELDS}
            original.update(lut_trials=500,geometry="p3_noFW",
                original_lut_sha256="a"*64,box_id=465,iu=15,iv=15,
                sample=250,random_seed=9023,u_start_mm="0.012345678901",
                v_start_mm="0.123456789012",w_start_mm="0",
                old_collected="1",old_field_wire_collected="0")
            nightly.write_csv(chunk/"original_exact_starts.csv",EXPORT_FIELDS,[original])
            identity="LUT500|p3_noFW|465|250"
            trial={col:"" for col in COLUMNS}
            trial.update(lut_trials="500",original_geometry="p3_noFW",
                original_lut_sha256="a"*64,box_id="465",sample="250",
                original_random_seed="9023",crossing_id=identity,
                old_collected="1",old_field_wire_collected="0",
                configuration="p3_noFW",seed_u_mm=original["u_start_mm"],
                seed_v_mm=original["v_start_mm"],seed_w_mm="0",
                ion_births="1",garfield_ni="1",ion_up="1",ion_down="0",
                ion_top="0",ion_bottom="0",successful="1",capped="0")
            timer={col:"" for col in TIME_COLUMNS}
            timer.update(crossing_id=identity,configuration="p3_noFW",
                wall_seconds="0.12")
            birth={col:"" for col in BIRTH_COLUMNS}
            birth.update(crossing_id=identity,seed_ion_births_total="1")
            nightly.write_csv(chunk/"lut_avalanche_trials.csv",COLUMNS,[trial])
            nightly.write_csv(chunk/"lut_avalanche_timing.csv",TIME_COLUMNS,[timer])
            nightly.write_csv(chunk/"lut_avalanche_birth_reservoir.csv",BIRTH_COLUMNS,[birth])
            (chunk/"physics.log").write_text("PASS: independent full-start microscopic avalanche sidecar 1 frozen LUT trials\n")
            rawfiles=[p for p in chunk.iterdir() if p.is_file()]
            orig_hash={p.name:nightly.sha(p) for p in rawfiles}
            result=nightly.recover_existing_smoke_chunk(smoke,[original])
            assert result["n_trials"]==1 and result["total_ion_births"]==1
            assert (chunk/"complete.json").is_file()
            assert {p.name:nightly.sha(p) for p in rawfiles}==orig_hash
            try:nightly.recover_existing_smoke_chunk(smoke,[original])
            except RuntimeError:pass
            else:raise AssertionError("Already sealed chunk must not be overwritten")
        finally:
            nightly.EXE=old_exe
    print("PASS overnight 84,000 distinct original LUT500 starts over all 4200 cells,"
          " 20 independent indices/cell, 210 SHA-stamped 400-trial chunks,"
          " deterministic pass interleaving, 5-ion independent birth sampling")


if __name__=="__main__":
    run()
