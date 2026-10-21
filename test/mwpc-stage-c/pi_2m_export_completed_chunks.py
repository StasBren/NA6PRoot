#!/usr/bin/env python3
"""Export VERIFIED full-geometry Geant4 parent chunks into native TGeo gas hits.

Read-only on Geant4 ROOT checkpoint files and frozen Garfield/Heed LUTs.
Each ROOT chunk's local event index 0..49999 is NEVER interpreted as a
global parent ID. The pair (source, chunk, event) uniquely identifies the
generated parent. Same ROOT macro used by first-chunk validated pipeline.
An immutable per-chunk SHA manifest prevents mixing reordered/changed
inputs and allows restart without exporting completed chunks again.

No Geant4 physics, no Garfield physics and no Heed transport in this tool.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import os
import subprocess
import sys
import tempfile
import time
from pathlib import Path

from run_pi_geant4_campaign import CHANNELS, BASE_SEED, REPO, sha256

HERE=Path(__file__).resolve().parent
SUITE=REPO/"test_runs/mwpc_stage_c/pi_geant4_2m_parent_campaign_v1"
DEST=REPO/"test_runs/mwpc_stage_c/pi_six_geometry_2m_native_exports_v1"
EXPORT_MACRO=HERE/"export_c_geant_muon_crossings.C"
PARENTS_PER_CHUNK=50_000
SOURCE=dict(Jpsi="jpsi40-bore",Phi="phi40-bore",Omega="omega40-bore")
REQUIRED_FIELDS=(
    "event","hit_index","track_id","pdg","det_id","station",
    "u_in_mm","v_in_mm","w_in_mm","u_out_mm","v_out_mm","w_out_mm",
    "path_length_mm","theta_u_in_deg","theta_w_in_deg",
)


def checkpoint(suite,channel,i):
    chunk=suite/channel.lower()/f"chunk{i:04d}"
    marker=chunk/"completed.json"
    if not marker.is_file():
        return None
    info=json.loads(marker.read_text(encoding="utf-8"))
    if (info.get("channel")!=channel or info.get("chunk")!=i or
        info.get("events")!=PARENTS_PER_CHUNK or
        info.get("seed")!=BASE_SEED[channel]+(i-1)*251 or
        info.get("verified_gas_step_tree") is not True or
        info.get("verified_virtual_muon_plane_tree") is not True):
        raise ValueError(f"Invalid full-geometry completed checkpoint: {marker}")
    needed=("HitsMuonSpec.root","geometry.root","MCKine.root")
    frozen=info.get("output_sha256",{})
    for name in needed:
        path=chunk/name
        if not path.is_file() or not path.stat().st_size:
            raise ValueError(f"Missing physical checkpoint ROOT file: {path}")
        if frozen.get(name)!=sha256(path):
            raise ValueError(f"Native ROOT file changed after checkpoint: {path}")
    return chunk,info


def export_one(root_exe,suite,dest,channel,i,dry=False):
    found=checkpoint(suite,channel,i)
    if found is None:
        raise ValueError(f"No validated {channel} chunk{i:04d}")
    chunk,info=found
    macro_sha=sha256(EXPORT_MACRO)
    work=dest/channel.lower()/f"chunk{i:04d}"
    out=work/"geant_muon_crossings.csv"
    manifest=work/"crossings_manifest.json"
    root_hashes={name:info["output_sha256"][name]
                 for name in ("HitsMuonSpec.root","geometry.root","MCKine.root")}
    immutable={
        "source":SOURCE[channel],"channel":channel,"chunk":i,
        "n_generated_parent_events":PARENTS_PER_CHUNK,
        "global_parent_index_interval":[(i-1)*PARENTS_PER_CHUNK,
                                        i*PARENTS_PER_CHUNK],
        "seed":info["seed"],
        "input_root_sha256":root_hashes,
        "exporter_macro_sha256":macro_sha,
        "source_protocol":"verified full 298-module native Geant4",
        "native_event_ID_semantics":"local to one ROOT chunk; pair with chunk",
    }
    if manifest.exists():
        meta=json.loads(manifest.read_text(encoding="utf-8"))
        if meta.get("immutable_input")!=immutable:
            raise ValueError(f"Existing native CSV does not match frozen ROOT inputs: {manifest}")
        if not out.exists() or sha256(out)!=meta.get("csv_sha256"):
            raise ValueError(f"Existing native CSV corrupted: {out}")
        print(f"SKIP verified {channel} chunk{i:04d}: "
              f"{meta['n_gas_crossings']:,} native rows",flush=True)
        return meta
    if out.exists():
        raise ValueError(f"Unmanifested CSV exists; refusing overwrite: {out}")
    if dry:
        print(f"READY to export {channel} chunk{i:04d} "
              f"from verified real Geant4 ({PARENTS_PER_CHUNK:,} parents)")
        return None
    if not root_exe:
        raise ValueError("ROOT executable required for real export")
    work.mkdir(parents=True,exist_ok=True)
    tmp=work/"geant_muon_crossings.csv.incomplete"
    if tmp.exists():
        tmp.unlink()  # This file is never a completed/manifested result.
    # ROOT ACLiC reads THIS checkout's existing export macro and writes
    # a TEMPORARY CSV which is only published after strict validation.
    invocation=(f'{EXPORT_MACRO}+("{chunk/"HitsMuonSpec.root"}",'
                f'"{chunk/"geometry.root"}","{tmp}",'
                f'"{chunk/"MCKine.root"}")')
    log=work/"export.log"
    t0=time.perf_counter()
    with log.open("w",encoding="utf-8") as handle:
        run=subprocess.run(
            [root_exe,"-l","-b","-q",invocation],
            cwd=REPO,stdout=handle,stderr=subprocess.STDOUT,
            check=False)
    if run.returncode or "GEANT MWPC crossing export COMPLETE" not in (
        log.read_text(encoding="utf-8",errors="replace")):
        raise RuntimeError(
            f"Failed ROOT native 298-det TGeo export: {log}\n"+
            "\n".join(log.read_text(errors="replace").splitlines()[-30:]))
    if not tmp.is_file() or not tmp.stat().st_size:
        raise RuntimeError(f"Missing complete output CSV: {tmp}")
    count=0;muons=0;invalid=0
    with tmp.open(newline="",encoding="utf-8") as handle:
        reader=csv.DictReader(handle)
        if not set(REQUIRED_FIELDS).issubset(reader.fieldnames or []):
            raise ValueError("Native CSV field schema does not match original exporter")
        for row in reader:
            event=int(row["event"])
            if not 0<=event<PARENTS_PER_CHUNK:
                raise ValueError(f"ROOT export leaked event index outside chunk: {event}")
            count+=1
            if abs(int(row["pdg"]))==13:
                muons+=1
            elif int(row["pdg"])==0:
                invalid+=1
    if invalid:
        raise ValueError(f"Unjoined muon PDG truth in source {channel}: {invalid}")
    meta={
        "immutable_input":immutable,
        "csv_sha256":sha256(tmp),
        "n_gas_crossings":count,"n_muon_gas_crossings":muons,
        "csv_size_bytes":tmp.stat().st_size,
        "measured_export_wall_seconds":time.perf_counter()-t0,
    }
    # Publish the data atomically BEFORE publishing the completion manifest.
    os.replace(tmp,out)
    mtemp=manifest.with_name(manifest.name+".tmp")
    mtemp.write_text(json.dumps(meta,indent=2,ensure_ascii=False)+"\n",
                     encoding="utf-8")
    os.replace(mtemp,manifest)
    print(f"PASS: native export {channel} chunk{i:04d}: "
          f"{count:,} crossings, {muons:,} muon crossings, "
          f"{out.stat().st_size/1048576:.1f} MiB, "
          f"{meta['measured_export_wall_seconds']:.1f}s",flush=True)
    return meta


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--suite",type=Path,default=SUITE)
    p.add_argument("--out",type=Path,default=DEST)
    p.add_argument("--root",default="root",
                   help="Full absolute path of the ROOT executable")
    p.add_argument("--max-chunks-per-source",type=int,default=2)
    p.add_argument("--dry-run",action="store_true")
    args=p.parse_args()
    if not 1<=args.max_chunks_per_source<=40:
        p.error("expected 1..40 50k native chunks per source")
    suite=(args.suite).resolve()
    frozen=json.loads((suite/"campaign.json").read_text(encoding="utf-8"))[
        "frozen_inputs"]
    if (frozen.get("events_per_source")!=2_000_000 or
        frozen.get("chunk_size")!=PARENTS_PER_CHUNK or
        frozen.get("channels")!=list(CHANNELS)):
        raise ValueError("Refusing incompatible Geant4 native campaign")
    if not EXPORT_MACRO.is_file():
        raise FileNotFoundError(EXPORT_MACRO)
    selected=[]
    for channel in CHANNELS:
        for i in range(1,args.max_chunks_per_source+1):
            item=checkpoint(suite,channel,i)
            if item is None:
                raise ValueError(f"Benchmark incomplete: {channel} chunk{i:04d}")
            selected.append((channel,i))
    if not args.dry_run:
        exe=Path(args.root).resolve()
        if not exe.is_file():
            raise FileNotFoundError(f"ROOT executable missing: {exe}")
    else:
        exe=None
    for channel,i in selected:
        export_one(str(exe) if exe else None,suite,args.out,
                   channel,i,dry=args.dry_run)
    print(f"PASS: {len(selected)} verified Geant4 chunk exports "
          "for native Garfield/Heed folding; no new Geant4.")


if __name__=="__main__":
    main()
