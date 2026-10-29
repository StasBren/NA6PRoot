#!/usr/bin/env python3
"""No-physics regression for the independent 120 original Geant4-muon holdout.

ALL inputs below are synthetic tests of selection, SHA sealing and physics
input identity. This DOES NOT claim to run or validate Garfield++ or TrackHeed.
"""
from __future__ import annotations

import json
import tempfile
from collections import Counter
from pathlib import Path

import pi_holdout120_A20_from_D500 as h

def fake_pool():
    # Synthetic Geant4 identities; 12 strata × 16 distinct fake tracks.
    pool={}
    for st in range(6):
        for sign in (-1,1):
            rows=[]
            for j in range(16):
                row={
                    "source":"FAKE-NOT-GEANT4",
                    "chunk":"2","event":f"{st}_{sign}_{j}",
                    "track_id":"13","station":str(st),
                }
                rows.append((j,("fake",str(st),str(sign),str(j)),row))
            pool[st,sign]=rows
    return pool

def run():
    assert h.N_MUONS==120 and h.N_CHUNKS==12 and h.GROUP==10
    assert h.REPLICATES==200 and h.SOURCE_CHUNK==2
    arr=list(range(120))
    chunks=h.groups(arr)
    assert len(chunks)==12
    assert all(len(c)==10 for c in chunks)
    assert [x for sub in chunks for x in sub]==arr
    pool=fake_pool()
    selection,stat=h.choose_holdout_from_pool(pool,10)
    assert len(selection)==120 and len(stat)==12
    assert all(v==10 for v in stat.values())
    assert len({(r["source"],r["chunk"],r["event"],r["track_id"])
                for r in selection})==120
    # Two distinct chamber gas hits with the SAME original physical
    # Geant4 particle MUST NOT count as two independent muons.
    pool=fake_pool()
    pool[0,-1][0][2]["event"]=pool[0,1][0][2]["event"]
    selection,stat=h.choose_holdout_from_pool(pool,10)
    assert len(selection)==120 and len({
        (r["source"],r["chunk"],r["event"],r["track_id"])
        for r in selection
    })==120
    pool=fake_pool()
    pool[0,1]=pool[0,1][:4]
    try:h.choose_holdout_from_pool(pool,10)
    except ValueError:pass
    else:raise AssertionError("Accepted a 10-muon stratum with only four candidates")
    # An existing original Garfield checkpoint is immutable and must
    # reject changed input or even a changed mock physics log.
    with tempfile.TemporaryDirectory() as td:
        c=Path(td)/"chunk_00";c.mkdir()
        original=[]
        for i in range(10):
            r={k:"" for k in h.FIELDS}
            r.update(source="SYNTHETIC-UNIT-TEST-ONLY",chunk="2",event=str(i),
                     hit_index="1",track_id="13",station="0",det_id="0")
            original.append(r)
        inp=c/"native_muon_crossings_input.csv"
        h.write_csv(inp,h.FIELDS,original)
        log=c/"physics.log"
        log.write_text("FAKE SHA CHECKPOINT UNIT TEST; NO GARFIELD\n")
        fake_binary="FAKE-NOT-REAL-CPP-SHA"
        seal={
            "schema":"PI_10_REAL_GEANT4_MUON_GARFIELD_SHA_CHECKPOINT_V1",
            "status":"complete","chunk":0,"n_muons":10,
            "cpp_executable_sha256":fake_binary,
            "input_sha256":h.sha(inp),
            "outputs_sha256":{inp.name:h.sha(inp),log.name:h.sha(log)},
        }
        h.write_json(c/"complete.json",seal)
        assert h.read_seal(c,0,original,fake_binary)["status"]=="complete"
        try:h.read_seal(c,0,list(reversed(original)),fake_binary)
        except ValueError:pass
        else:raise AssertionError("Wrong original physical muon order was accepted")
        try:h.read_seal(c,0,original,"DIFFERENT-BINARY")
        except ValueError:pass
        else:raise AssertionError("Wrong compiled C++ executable was accepted")
        log.write_text("TAMPERED")
        try:h.read_seal(c,0,original,fake_binary)
        except ValueError:pass
        else:raise AssertionError("Physics log SHA tampering was accepted")
    print("PASS OFFLINE ONLY: A20@D500 naming, 120 unique source Geant4 track "
          "selection in 12 balanced signed-angle/station strata, 12x10 chunk "
          "identity, resume seal and tamper/wrong-binary refusal. "
          "NO Geant4/Heed/Garfield physics executed.")

if __name__=="__main__":
    run()
