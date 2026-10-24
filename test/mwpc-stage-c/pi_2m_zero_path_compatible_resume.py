#!/usr/bin/env python3
"""Resume the exact existing frozen six-LUT fold with a documented, strict
zero-length gas crossing exception; preserve ALL previous checkpoint SHAs.

Critical: DO NOT modify pi_six_geometry_native_timing.py or
pi_2m_native_six_geometry_fold.py. Their source SHAs are stored in all
already-completed manifests (including the initial 100k/sample); modifying
either would invalidate the whole ongoing campaign.

This wrapper replaces only the in-memory call to birth_sample with a
narrow exception for numerically EXACTLY zero chord length AND
coincident finite in-gas endpoints and Edep < 1 eV. The Heed path-length
Poisson intensity is zero, hence 0 generated clusters/electrons. The native
MuonSpec gas crossing itself is NOT discarded: in aggregate() its
"segments" counter is incremented for all 6 geometries and the normal
signed-angle/in-chamber provenance is retained. All non-degenerate paths
use the original sampler byte-for-byte with the same seed.

An explicit exception ledger is written for EVERY newly folded shard,
with its already-verified SHA and separate policy source SHA. Original
per-shard immutable checkpoint keys remain IDENTICAL, allowing the old
100k chunks to coexist. Original ROOT, CSV and LUT are never modified.
"""
from __future__ import annotations
import csv
import hashlib
import json
import math
import os
import sys
from pathlib import Path
import numpy as np

import pi_six_geometry_native_timing as timing
import pi_2m_native_six_geometry_fold as fold

ORIGINAL=timing.birth_sample
VERSION="PI_NATIVE_ZERO_CHORD_STRICT_COMPATIBILITY_V1"
HERE=Path(__file__).resolve()
NORMALIZED_EDEP_LIMIT_EV=1.0


def hashed(p):
    h=hashlib.sha256()
    with Path(p).open("rb") as stream:
        for chunk in iter(lambda:stream.read(4*1024*1024),b""):
            h.update(chunk)
    return h.hexdigest()


def coincident_zero_chord(row):
    try:
        length=float(row["path_length_mm"])
        edep=float(row["edep_eV"])
        points=np.asarray(
            [[float(row[f"{a}_in_mm"]) for a in "uvw"],
             [float(row[f"{a}_out_mm"]) for a in "uvw"]],dtype=np.float64)
    except (KeyError,ValueError,TypeError) as err:
        raise ValueError("Bad native crossing schema, not a physical zero chord") from err
    if (length != 0. or not np.isfinite(points).all() or
        not np.array_equal(points[0],points[1]) or
        abs(points[:,1]).max()>3.05 or not math.isfinite(edep) or
        not 0<=edep<NORMALIZED_EDEP_LIMIT_EV):
        return False
    return True


seen=[]
def guarded_birth_sample(row,lam,marks,seed):
    # Do not change the original physics or its RNG for valid entries.
    if not coincident_zero_chord(row):
        return ORIGINAL(row,lam,marks,seed)
    seen.append({
        "event":int(row["event"]),"hit_index":int(row["hit_index"]),
        "track_id":int(row["track_id"]),"pdg":int(row["pdg"]),
        "det_id":int(row["det_id"]),"station":int(row["station"]),
        "length_mm":0.,"edep_eV":float(row["edep_eV"]),
        "native_coords_mm":{a:float(row[f"{a}_in_mm"]) for a in "uvw"},
    })
    # Native segment still counted by unchanged aggregate for 6 models.
    # Zero path => Poisson(lambda*0)=Poisson(0): zero primary electrons.
    return np.empty((0,3),dtype=np.float64),np.empty(0,dtype=np.float64)


def fold_one_recorded(channel,index,station,limit,pools,provenance,lam,marks,out):
    path=fold.shard_path(out,channel,index,station)
    original_manifest=path.with_suffix(".manifest.json")
    existed=original_manifest.exists()
    seen.clear()
    RAW_FOLD_ONE(channel,index,station,limit,pools,provenance,lam,marks,out)
    if existed:
        return
    # The main native completed manifest remains identical in schema and
    # source hashes to all previously completed 100k+ chunks.
    if not path.is_file() or not original_manifest.is_file():
        raise AssertionError(f"Expected validated new folded shard: {path}")
    completed=json.loads(original_manifest.read_text(encoding="utf-8"))
    ledger=path.with_suffix(".zero_path_policy.json")
    if ledger.exists():
        raise FileExistsError(f"Refuse overwrite physics exception record: {ledger}")
    content={
        "policy":VERSION,
        "policy_sha256":hashed(HERE),
        "frozen_original_aggregate_source_sha256":hashed(Path(timing.__file__)),
        "frozen_original_fold_source_sha256":hashed(Path(fold.__file__)),
        "source":fold.SOURCE[channel],"chunk":index,"station":station,
        "native_input_csv_sha256":
            completed["immutable_inputs"]["native_csv_sha256"],
        "folded_gzip_sha256":completed["gzip_sha256"],
        "original_completed_manifest_sha256":hashed(original_manifest),
        "interpretation":"zero chord and identical in-gas endpoints, "
                         "Edep<1eV: keep native gas crossing; 0 primary "
                         "Heed electrons; all non-zero paths untouched.",
        "degenerate_crossing_count":len(seen),
        "degenerate_crossings":list(seen),
    }
    temp=ledger.with_name(ledger.name+".incomplete")
    temp.write_text(json.dumps(content,indent=2)+"\n",encoding="utf-8")
    os.replace(temp,ledger)
    if seen:
        print(f"PHYSICS EXCEPTION: {channel} chunk{index:04d} MS{station}: "
              f"{len(seen)} true degenerate zero-length native crossing(s); "
              "zero Heed primaries, preserved in native gas exposure. "
              f"Audit: {ledger}",flush=True)
    else:
        print(f"PASS exact original Heed physics on "
              f"{channel} chunk{index:04d} MS{station}; "
              "zero-length exceptions: 0",flush=True)


RAW_FOLD_ONE=fold.fold_one


def run():
    if len(sys.argv)<2 or sys.argv[1] not in ("fold","smoke"):
        raise SystemExit("Usage: python3 pi_2m_zero_path_compatible_resume.py "
                         "fold --channel all --station -1 --chunks 40")
    timing.birth_sample=guarded_birth_sample
    fold.fold_one=fold_one_recorded
    fold.main()

if __name__=="__main__":
    run()
