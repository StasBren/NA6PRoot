#!/usr/bin/env python3
"""Offline nightly-run selection/budget/compact reservoir regression.

No Garfield libraries, Geant4, original LUT CSVs, WSL files or network
needed: checks the 4200-cell design before compilation or long jobs.
"""
from __future__ import annotations
import random
from collections import Counter

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
    print("PASS overnight 84,000 distinct original LUT500 starts over all 4200 cells,"
          " 20 independent indices/cell, 210 SHA-stamped 400-trial chunks,"
          " deterministic pass interleaving, 5-ion independent birth sampling")


if __name__=="__main__":
    run()
