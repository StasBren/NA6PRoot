#!/usr/bin/env python3
"""Simulation-free checkpoint SHA, source/chunk parent-ID validation tests."""
from __future__ import annotations
import json
import tempfile
from pathlib import Path
from pi_2m_export_completed_chunks import (
    checkpoint,export_one,SOURCE,PARENTS_PER_CHUNK,
)
from run_pi_geant4_campaign import BASE_SEED,sha256


def fixture(root,ch,idx):
    loc=root/ch.lower()/f"chunk{idx:04d}"
    loc.mkdir(parents=True)
    frozen={}
    for name in ("HitsMuonSpec.root","geometry.root","MCKine.root"):
        p=loc/name
        p.write_bytes(f"synthetic {ch} {idx} {name}".encode())
        frozen[name]=sha256(p)
    (loc/"completed.json").write_text(json.dumps({
        "channel":ch,"chunk":idx,"events":PARENTS_PER_CHUNK,
        "seed":BASE_SEED[ch]+(idx-1)*251,
        "verified_gas_step_tree":True,
        "verified_virtual_muon_plane_tree":True,
        "output_sha256":frozen,
    }),encoding="utf-8")
    return loc


def must_fail(f,msg):
    try:
        f()
    except (ValueError,FileNotFoundError):
        return
    raise AssertionError("Expected strict checkpoint failure: "+msg)


def main():
    with tempfile.TemporaryDirectory(prefix="pi_2m_export_synthetic_") as t:
        suite=Path(t)/"campaign"
        out=Path(t)/"out"
        for ch in ("Jpsi","Phi","Omega"):
            assert SOURCE[ch].endswith("40-bore")
            for i in (1,2):
                fixture(suite,ch,i)
                loc,meta=checkpoint(suite,ch,i)
                assert loc.name==f"chunk{i:04d}"
                assert meta["seed"]==BASE_SEED[ch]+(i-1)*251
                assert export_one(None,suite,out,ch,i,dry=True) is None
        assert not out.exists()  # all dry-run, no exported files
        hits=suite/"jpsi"/"chunk0001"/"HitsMuonSpec.root"
        hits.write_bytes(b"modified after validation")
        must_fail(lambda:checkpoint(suite,"Jpsi",1),
                  "tampered original muon crossings ROOT")
        done=suite/"phi"/"chunk0002"/"completed.json"
        info=json.loads(done.read_text());info["seed"]+=1
        done.write_text(json.dumps(info))
        must_fail(lambda:checkpoint(suite,"Phi",2),
                  "wrong random seed for local chunk")
        missing=checkpoint(suite,"Jpsi",3)
        assert missing is None
    print("PASS: native 2M Geant4 ROOT exporter 6 dry chunk preflights, "
          "all source/chunk/seed pairs distinct, SHA tamper detection, "
          "NO ROOT or Garfield execution")


if __name__=="__main__":
    main()
