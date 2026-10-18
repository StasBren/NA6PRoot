#!/usr/bin/env python3
"""Read-only audit of the finished 3-channel native Geant4 MWPC campaign.

Use only the frozen campaign manifest and per-chunk validated checkpoints;
never assume CAMPAIGN_COMPLETE merely because a shell stopped running.
No Geant4 rerun, no merging, no ROOT imports and no file modifications.
"""
from __future__ import annotations
import argparse
import json
import math
import re
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
DEFAULT = REPO / "test_runs/mwpc_stage_c/pi_geant4_hit_validation_v1"
FILES = ("MCKine.root", "HitsMuonSpec.root", "geometry.root",
         "MuonReferencePlanes.root", "na6pLayout.ini")
CHANNELS = ("Jpsi", "Phi", "Omega")


def audit(suite: Path, verbose: bool = False) -> dict:
    manifest_path = suite / "campaign.json"
    if not manifest_path.is_file():
        raise FileNotFoundError(f"No frozen campaign manifest: {manifest_path}")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    frozen = manifest["frozen_inputs"]
    if frozen["mode"] != "full" or list(frozen["channels"]) != list(CHANNELS):
        raise ValueError("This is not the frozen full three-channel campaign")
    total = int(frozen["events_per_source"])
    size = int(frozen["chunk_size"])
    nworkers = int(frozen["workers"])
    nchunks = math.ceil(total / size)
    if not (total > 0 and size > 0 and 1 <= nworkers <= 32):
        raise ValueError("Invalid frozen event count or chunk size")
    result = {
        "status": "UNKNOWN", "campaign": str(suite),
        "expected_events_per_source": total,
        "expected_chunks_per_source": nchunks,
        "channels": {}, "problems": [],
    }
    for channel in CHANNELS:
        completed_events = 0
        complete = 0
        steps = 0
        reference_states = 0
        validation_data = 0
        files_bytes = 0
        for index in range(nchunks):
            chunk = suite / channel.lower() / f"chunk{index+1:04d}"
            done = chunk / "completed.json"
            n = min(size, total - index*size)
            seed = {"Jpsi":202610181, "Phi":202710181,
                    "Omega":202810181}[channel] + index*251
            if not done.is_file():
                result["problems"].append(f"{channel}: chunk{index+1:04d} not complete")
                continue
            try:
                entry = json.loads(done.read_text(encoding="utf-8"))
                if (entry.get("events") != n or entry.get("seed") != seed or
                    entry.get("channel") != channel or
                    entry.get("chunk") != index + 1):
                    raise ValueError("wrong event count, seed or source")
                if not entry.get("verified_gas_step_tree"):
                    raise ValueError("gas step check absent")
                if not entry.get("verified_virtual_muon_plane_tree"):
                    raise ValueError("reference-plane check absent")
                for f in FILES:
                    p = chunk / f
                    if not p.is_file() or p.stat().st_size == 0:
                        raise ValueError(f"missing {f}")
                    files_bytes += p.stat().st_size
                original = chunk / "raw_workers"
                if not original.is_dir():
                    raise ValueError("missing original raw_workers")
                for iw in range(min(nworkers,n)):
                    worker = original / f"worker{iw:03d}"
                    if not worker.is_dir():
                        raise ValueError(f"missing original worker {iw}")
                log = chunk / "validation.log"
                if not log.is_file():
                    raise ValueError("missing native ROOT validation log")
                entries = [
                    line for line in log.read_text(errors="replace").splitlines()
                    if "PI_GEANT_CHUNK_PASS" in line
                ]
                if len(entries) != 1:
                    raise ValueError("expected exactly one native ROOT PASS marker")
                pattern = r"events=(\d+)\s+virtualPlaneStates=(\d+)\s+gasSteps=(\d+)"
                m = re.search(pattern, entries[0])
                if not m or int(m[1]) != n:
                    raise ValueError("malformed ROOT PASS/count")
                steps += int(m[3])
                reference_states += int(m[2])
                validation_data += 1
            except (ValueError, KeyError, OSError) as exc:
                result["problems"].append(f"{channel}: chunk{index+1:04d}: {exc}")
                continue
            completed_events += n
            complete += 1
        result["channels"][channel] = {
            "completed_chunks":complete,"required_chunks":nchunks,
            "verified_events":completed_events,"expected_events":total,
            "raw_gas_steps":steps,
            "interpolated_muon_reference_plane_rows":reference_states,
            "bytes_five_mandatory_outputs":files_bytes,
        }
        print(f"{channel:<6} {complete:3d}/{nchunks:<3d} chunks | "
              f"events {completed_events:,}/{total:,} | "
              f"Geant4 Edep steps {steps:,} | "
              f"virtual-plane states {reference_states:,}",flush=True)
    # Require BOTH the separate campaign completion certificate and every
    # completed chunk. Merely reaching the elapsed-time budget is not enough.
    marker = suite / "campaign_complete.json"
    if not marker.is_file():
        result["problems"].append("campaign_complete.json missing")
    else:
        try:
            cert = json.loads(marker.read_text(encoding="utf-8"))
            if (cert.get("events_per_source") != total or
                cert.get("n_chunks") != nchunks*len(CHANNELS) or
                cert.get("channels") != list(CHANNELS)):
                result["problems"].append("campaign completion certificate disagrees with manifest")
        except (ValueError, KeyError):
            result["problems"].append("invalid campaign completion certificate")
    result["status"] = "COMPLETE" if not result["problems"] else "INCOMPLETE"
    print("FINAL CAMPAIGN STATUS:",result["status"],flush=True)
    if result["problems"]:
        print("Incompleteness details (first 15):",flush=True)
        for problem in result["problems"][:15]:
            print(" -",problem,flush=True)
    if verbose:
        print(json.dumps(result,ensure_ascii=False,indent=2))
    return result


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--campaign",type=Path,default=DEFAULT)
    ap.add_argument("--verbose",action="store_true")
    args=ap.parse_args()
    try:
        result=audit(args.campaign,args.verbose)
    except (OSError,ValueError,KeyError) as exc:
        print(f"CAMPAIGN_AUDIT_ERROR: {exc}",flush=True)
        raise SystemExit(2) from exc
    if result["status"] != "COMPLETE":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
