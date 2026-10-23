#!/usr/bin/env python3
"""Archive exact bytes of frozen 6M-parent Geant4 production inputs.

This avoids a subtle irreversible-data loss: campaign.json retained only
the SHA of field maps / simcuts and some installed binaries, not their
actual bytes. Do not modify frozen inputs or run Geant4 if anything has
changed. Idempotent and write-once snapshots; only archives code/field/
installed-executable inputs, NOT ROOT event data.
"""
from __future__ import annotations
import hashlib
import json
import os
import shutil
from pathlib import Path

from pi_2m_export_completed_chunks import SUITE
from run_pi_geant4_campaign import REPO

OUT=REPO/"test_runs/mwpc_stage_c/pi_weekend_2m_operations_v1"/"frozen_inputs"
MANIFEST=OUT/"snapshot_manifest.json"


def sha(path:Path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda:f.read(4*1024*1024),b""):
            h.update(chunk)
    return h.hexdigest()


def source_files(frozen):
    specs={
        "generator_original.C":(frozen["generator"],frozen["generator_sha256"]),
        "mwpc_hook_original.C":(frozen["hook"],frozen["hook_sha256"]),
        "original_validator.C":(frozen["validator"],frozen["validator_sha256"]
            if "validator_sha256" in frozen else None),
        "na6psim_original":(frozen["simulator"],frozen["simulator_sha256"]),
        "na6psim_parallel_original":(
            frozen["parallel_runner"],frozen["parallel_runner_sha256"]),
        "libsimLib_original.so":(
            frozen["loaded_libsimLib"],frozen["loaded_libsimLib_sha256"]),
    }
    for original,expected in frozen.get("installed_field_and_cut_sha256",{}).items():
        key="installed_field_or_cut_"+Path(original).name
        if key in specs:raise ValueError("Duplicate field/cut basename")
        specs[key]=(original,expected)
    return specs


def main():
    campaign=SUITE/"campaign.json"
    if not campaign.is_file():raise FileNotFoundError(campaign)
    manifest=json.loads(campaign.read_text(encoding="utf-8"))
    frozen=manifest["frozen_inputs"]
    if (frozen.get("events_per_source")!=2_000_000 or
        frozen.get("chunk_size")!=50000 or
        frozen.get("mode")!="full"):
        raise RuntimeError("Refusing snapshot of a different production campaign")
    specs=source_files(frozen)
    OUT.mkdir(parents=True,exist_ok=True)
    receipts={}
    for alias,(original,expected) in specs.items():
        path=Path(original).resolve()
        if not path.is_file():
            raise FileNotFoundError(f"Original frozen production input missing: {path}")
        actual=sha(path)
        if expected is not None and actual!=expected:
            raise ValueError(
                f"Original frozen physics or binary changed: {path}\n"
                f"expected {expected}, actual {actual}; do not start Geant4")
        target=OUT/alias
        if target.exists():
            if sha(target)!=actual:
                raise ValueError(f"Existing immutable snapshot disagrees: {target}")
        else:
            temp=OUT/(alias+".incomplete")
            if temp.exists():temp.unlink()
            shutil.copy2(path,temp)
            if sha(temp)!=actual:
                raise RuntimeError(f"Snapshot copy failed checksum: {path}")
            os.replace(temp,target)
        receipts[alias]={"original_path":str(path),"sha256":actual,
                         "original_frozen_sha256":expected,
                         "size_bytes":target.stat().st_size}
    payload={
        "schema":"PI_2M_IMMUTABLE_PRODUCTION_INPUT_ASSETS_V1",
        "campaign":str(campaign),
        "campaign_sha256":sha(campaign),
        "assets":receipts,
        "limitations":"ROOT/Geant4 external runtime libraries not embedded; "
                      "preserve WSL runtime and software versions separately",
    }
    if MANIFEST.is_file():
        existing=json.loads(MANIFEST.read_text(encoding="utf-8"))
        if existing!=payload:
            raise ValueError(f"Existing immutable snapshot manifest mismatch: {MANIFEST}")
    else:
        tmp=MANIFEST.with_name(MANIFEST.name+".incomplete")
        tmp.write_text(json.dumps(payload,indent=2)+"\n",encoding="utf-8")
        os.replace(tmp,MANIFEST)
    print(f"FROZEN INPUT SNAPSHOT PASS: {len(receipts)} physics/source/binary "
          f"assets copied and checksum verified. Archive: {OUT}",flush=True)


if __name__=="__main__":
    main()
