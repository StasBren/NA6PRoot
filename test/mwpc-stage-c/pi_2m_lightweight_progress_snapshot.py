#!/usr/bin/env python3
"""Compact (<100 kB) resumable 6M parent Heed/6xGarfield status snapshot.

Read-only in source physics data; never touch original ROOT, native CSV,
six Garfield LUTs or main frozen Heed manifest SHAs. Can run WHILE folds
are being produced: require BOTH gzip+manifest presence and atomically
publish a descriptive status JSON for quick Google Drive reading.
The counts show existing checkpoint manifests, not cryptographic rechecks
of the gzip payload (done by the native resume process).
"""
from __future__ import annotations
import json
import os
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

REPO=Path(__file__).resolve().parents[2]
BASE=REPO/"test_runs/mwpc_stage_c"
FOLDS=BASE/"pi_six_geometry_2m_timing_v1/full"
OPS=BASE/"pi_weekend_2m_operations_v1"
DEST=OPS/"lightweight_progress_snapshot.json"
CHANNELS=("jpsi","phi","omega")
CHUNKS=40
STATIONS=6


def build():
    summaries=[]
    totals=defaultdict(int)
    for source in CHANNELS:
        for st in range(STATIONS):
            row=dict(source=source,station=st,expected_chunks=CHUNKS,
                     completed_chunks=0,muon_gas_segments=0,
                     measured_fold_seconds=0.0,checkpoint_bytes=0,
                     zero_path_ledgers=0,zero_path_exceptions=0,
                     completed_chunk_indices=[])
            for chunk in range(1,CHUNKS+1):
                p=FOLDS/source/f"MS{st}"/f"chunk{chunk:04d}.json.gz"
                manifest=p.with_suffix(".manifest.json")
                if not (p.is_file() and manifest.is_file()):
                    continue
                try:
                    doc=json.loads(manifest.read_text(encoding="utf-8"))
                    inputs=doc["immutable_inputs"]
                    if (int(inputs["chunk"])!=chunk or
                        int(inputs["station"])!=st):
                        raise ValueError("misidentified frozen shard manifest")
                    row["completed_chunks"]+=1
                    row["completed_chunk_indices"].append(chunk)
                    row["muon_gas_segments"]+=int(doc["n_muon_gas_hits"])
                    row["measured_fold_seconds"]+=float(
                        doc.get("measured_fold_seconds",0))
                    row["checkpoint_bytes"]+=int(p.stat().st_size)
                    ledger=p.with_suffix(".zero_path_policy.json")
                    if ledger.is_file():
                        exc=json.loads(ledger.read_text(encoding="utf-8"))
                        if exc.get("folded_gzip_sha256")!=doc.get("gzip_sha256"):
                            raise ValueError("separate exception ledger SHA mismatch")
                        row["zero_path_ledgers"]+=1
                        row["zero_path_exceptions"]+=int(
                            exc["degenerate_crossing_count"])
                except (KeyError,ValueError,TypeError) as e:
                    raise RuntimeError(
                        f"Cannot trust existing completed fold manifest {manifest}: {e}"
                    ) from e
            summaries.append(row)
            totals["fold_shards_completed"]+=row["completed_chunks"]
            totals["zero_path_exceptions_in_completed_shards"]+=row["zero_path_exceptions"]
            totals["collected_checkpoint_bytes"]+=row["checkpoint_bytes"]
    audit=OPS/"native_gas_path_audit.json"
    aud=None
    if audit.is_file():
        d=json.loads(audit.read_text(encoding="utf-8"))
        aud={"validated_export_chunks":d.get("complete_source_chunks"),
             "native_muon_gas_crossings_without_zero_chords":
                d.get("n_valid_muon_crossings"),
             "zero_path_exceptions_over_all_exports":
                d.get("n_accepted_zero_chords")}
    return {
        "schema":"NA60PLUS_6M_PARENT_LIVE_STATUS_V1",
        "created_utc":datetime.now(timezone.utc).isoformat(),
        "physics_caveat":"Observes existing manifest+gzip pairs only; "
           "not a full ROOT, gzip SHA, or final 40-chunk physics merge.",
        "generated_parent_events_expected":6000000,
        "original_native_Geant4_chunks_expected":120,
        "fold_shards_expected":len(CHANNELS)*STATIONS*CHUNKS,
        "fold_shards_completed":totals["fold_shards_completed"],
        "fold_shards_remaining":720-totals["fold_shards_completed"],
        "checkpoint_gzip_bytes_present":totals["collected_checkpoint_bytes"],
        "zero_path_exceptions_in_completed_shards":
            totals["zero_path_exceptions_in_completed_shards"],
        "native_full_120_chunk_audit":aud,
        "source_station_details":summaries,
    }


def main():
    doc=build()
    DEST.parent.mkdir(parents=True,exist_ok=True)
    tmp=DEST.with_name(DEST.name+".incomplete")
    tmp.write_text(json.dumps(doc,indent=2,ensure_ascii=False)+"\n",
                   encoding="utf-8")
    os.replace(tmp,DEST)
    print(f"PASS: {doc['fold_shards_completed']}/720 complete frozen fold shards;"
          f" {doc['fold_shards_remaining']} remaining;"
          f" status snapshot {DEST}",flush=True)


if __name__=="__main__":
    main()
