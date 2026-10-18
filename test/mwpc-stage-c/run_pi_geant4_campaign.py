#!/usr/bin/env python3
"""Resumable modern Geant4 MWPC validation campaign with frozen provenance.

This is a HIGH-FIDELITY benchmark in the CURRENT 298-chamber material layout.
It keeps all native ROOT truth, actual gas entry/exit and energy deposits,
the optional individual positive-Edep Geant4 gas steps, and the TGeo geometry.
It does NOT claim that those gas interactions can be re-used unchanged after
moving, replacing, or materially modifying chambers.

For geometry variations, use the *separate* historical reference-plane
LocalTrajectoryAudit.root trajectory library for inexpensive projections.
For a new material geometry rerun this campaign (can share seeds and generator
configuration to reduce statistical variation, not guarantee identical tracks).

Example:
  python3 test/mwpc-stage-c/run_pi_geant4_campaign.py --mode smoke
  python3 test/mwpc-stage-c/run_pi_geant4_campaign.py --mode full \
      --events-per-source 30000 --chunk-size 3000 --workers 4 --hours 2.5

A whole chunk is the atomic checkpoint. Interrupted chunks are preserved under
an _interrupted directory and are never mistaken for completed events.
"""
from __future__ import annotations

import argparse
import datetime as dt
import fcntl
import hashlib
import json
import math
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

CHANNELS = ("Jpsi", "Phi", "Omega")
BASE_SEED = {"Jpsi": 202610181, "Phi": 202710181, "Omega": 202810181}
ESSENTIAL = ("MCKine.root", "HitsMuonSpec.root", "geometry.root",
             "MuonReferencePlanes.root", "na6pLayout.ini")
STUDY = "pi_geant4_hit_validation_v1"
HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json_atomic(path: Path, data: object):
    tmp = path.with_name(path.name + ".tmp")
    with tmp.open("w", encoding="utf-8") as stream:
        json.dump(data, stream, ensure_ascii=False, indent=2)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(tmp, path)


def require_program(name):
    program = shutil.which(name)
    if not program:
        raise RuntimeError(
            f"Required {name} missing from PATH. "
            "Source the existing NA6PRoot ROOT/Geant4 environment first."
        )
    return str(Path(program).resolve())


def setup(args):
    if args.mode == "smoke":
        if args.events_per_source != 3 or args.chunk_size != 3:
            raise RuntimeError(
                "Smoke mode requires its fixed 3-event test settings"
            )
        suite = (REPO / "test_runs/mwpc_stage_c" /
                 (STUDY + "_smoke")).resolve()
        args.workers = 1
        args.hours = 0.
    else:
        suite = (REPO / "test_runs/mwpc_stage_c" /
                 STUDY).resolve()
    if args.suite:
        suite = args.suite.expanduser().resolve()
    if (not 1 <= args.events_per_source <= 10**8 or
        not 1 <= args.chunk_size <= 50000 or
        not 1 <= args.workers <= 32 or
        not math.isfinite(args.hours) or args.hours < 0):
        raise ValueError("Invalid event, chunk, worker count, or hour budget")
    sim = require_program("na6psim")
    parallel = require_program("na6psim_parallel")
    require_program("hadd")
    root = require_program("root")

    # Freeze the actual dynamically loaded detector library too; changing
    # libsimLib.so without changing na6psim would otherwise silently mix
    # ROOT hits from different module construction/stepping code versions.
    try:
        dynamic = subprocess.check_output(
            ["ldd", sim], text=True, stderr=subprocess.STDOUT,
            timeout=15,
        )
    except Exception as exc:
        raise RuntimeError("Could not inspect native simulator library bindings") from exc
    simlib = None
    for row in dynamic.splitlines():
        if "libsimLib.so" in row and "=>" in row:
            candidate = row.split("=>", 1)[1].strip().split()[0]
            if candidate == "not" or not Path(candidate).is_file():
                raise RuntimeError("libsimLib.so missing in runtime loader path")
            simlib = str(Path(candidate).resolve())
            break
    if simlib is None:
        raise RuntimeError("Native executable does not link to expected libsimLib.so")

    generator = REPO / "test/genDimuonBgEvent.C"
    hook = REPO / "test/mwpc-dimuon-validation/mwpcDimuonHooks.C"
    validator = HERE / "validate_pi_geant4_chunk.C"
    for p in (generator, hook, validator):
        if not p.is_file():
            raise FileNotFoundError(p)

    inputs = {
        "protocol": STUDY,
        "mode": args.mode,
        "events_per_source": args.events_per_source,
        "chunk_size": args.chunk_size,
        "workers": args.workers,
        "channels": list(CHANNELS),
        "base_seeds": BASE_SEED,
        "physics_settings": "beam.energyPerNucleon=40",
        "digitization": False,
        "gas_step_diagnostics": True,
        "portable_reference_plane_audit": True,
        "reference_plane_schedule": "280..900cm/10cm + 365,430,495cm",
        "generator": str(generator),
        "generator_sha256": sha256(generator),
        "hook": str(hook),
        "hook_sha256": sha256(hook),
        "validator": str(validator),
        "simulator": sim,
        "simulator_sha256": sha256(Path(sim)),
        "loaded_libsimLib": simlib,
        "loaded_libsimLib_sha256": sha256(Path(simlib)),
        "parallel_runner": parallel,
        "parallel_runner_sha256": sha256(Path(parallel)),
        "root": root,
    }
    suite.mkdir(parents=True, exist_ok=True)
    manifest_path = suite / "campaign.json"
    if manifest_path.is_file():
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        if manifest.get("frozen_inputs") != inputs:
            raise RuntimeError(
                "Found different frozen simulation settings/executable. "
                f"Refusing to mix runs under {suite}. "
                "Use a different --suite or restore original tools/settings."
            )
    else:
        try:
            commit = subprocess.check_output(
                ["git", "-C", str(REPO), "rev-parse", "HEAD"],
                text=True, timeout=6,
            ).strip()
        except Exception:
            commit = "unknown"
        manifest = {
            "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
            "git_commit_at_creation": commit,
            "frozen_inputs": inputs,
            "physics_scope": (
                "Modern Geant4-VMC native full apparatus with current "
                "298-module MWPC material geometry, FTFP_BERT, installed "
                "field maps, and forced signal -> mu+mu-"
            ),
            "output_semantics": {
                "MCKine.root": "Generated parent and signal daughters / ancestry",
                "HitsMuonSpec.root": (
                    "Actual per-track per-module gas entrance/exit "
                    "x,p, accumulated Edep and time; gasStepsMuonSpec tree "
                    "stores positive-Edep step positions and momenta"
                ),
                "MuonReferencePlanes.root": (
                    "Geometry-independent 10-cm Z-grid muon transport snapshots "
                    "(linearly interpolated Geant4 steps; unremapped track IDs)"
                ),
                "geometry.root": "Geometry of THIS run, needed for correct local axes",
                "na6pLayout.ini": "Exact NA6PRoot configurable layout snapshot",
            },
            "limitations": [
                "Large Geant4 gas Edep depends on the material layout; "
                "do not reuse as truth for altered chamber volume placements",
                "The gas step tree stores the ORIGINAL track ID while "
                "MuonSpec hit vectors may store REMAPPED track IDs. "
                "Never join without handling the remapping.",
                "Saved gas-step diagnostics only include positive-Edep steps",
                "Forced signal dimuons are not absolute physical rates",
                "Reusing identical RNG seed in changed geometry does not "
                "guarantee identical individual particle trajectories",
                "Virtual Z-plane positions and momenta are linearly "
                "interpolated WITHIN actual Geant4 steps; they are not new "
                "sensitive detector hits or exact physics crossings",
                "Existing historical LocalTrajectoryAudit.root reference-plane "
                "data are a distinct independent dataset, not event-matched",
            ],
        }
        write_json_atomic(manifest_path, manifest)
    print("CAMPAIGN:", suite, flush=True)
    print(
        "Edep steps AND passive muon reference planes enabled. Every "
        "completed chunk retains truth, gas hits, spatial steps, Z-plane "
        "snapshots, geometry, and layout INI.",
        flush=True,
    )
    return suite, inputs


def simulate(args, suite, inputs):
    started = time.monotonic()
    end_time = started + 3600 * args.hours if args.hours > 0 else math.inf
    manifest = suite / "campaign.json"
    sim = inputs["parallel_runner"]
    root = inputs["root"]
    chunks_per_source = math.ceil(args.events_per_source / args.chunk_size)
    count_new = 0
    for ich, channel in enumerate(CHANNELS):
        for part in range(chunks_per_source):
            n = min(args.chunk_size,
                    args.events_per_source-part*args.chunk_size)
            seed = BASE_SEED[channel] + part*251
            out = suite / channel.lower() / f"chunk{part+1:04d}"
            done = out / "completed.json"
            if done.is_file():
                info = json.loads(done.read_text(encoding="utf-8"))
                if info.get("events") != n or info.get("seed") != seed:
                    raise RuntimeError(f"Checkpoint count/seed mismatch: {out}")
                if not all((out / fn).is_file() and
                           (out / fn).stat().st_size > 0
                           for fn in ESSENTIAL):
                    raise RuntimeError(
                        f"Completed checkpoint missing original outputs: {out}"
                    )
                print(f"SKIP completed {channel} chunk{part+1:04d}",flush=True)
                continue
            if time.monotonic() >= end_time:
                print(
                    "CAMPAIGN_PAUSED wall-clock limit reached between chunks. "
                    "Repeat the identical command to resume.",flush=True
                )
                return
            out.parent.mkdir(parents=True, exist_ok=True)
            if out.exists():
                stamp = dt.datetime.now().strftime("%Y%m%d_%H%M%S_%f")
                old = out.with_name(out.name+"_interrupted_"+stamp)
                out.rename(old)
                print(f"Archived interrupted chunk: {old}", flush=True)
            out.mkdir()
            command = [
                sim, "--workers", str(args.workers),
                "-n", str(n), "-r", str(seed),
                "-g", (f'{inputs["generator"]}+'
                       f'(1,"{channel}",0.,5.,0.,6.,false)'),
                "-u", inputs["hook"]+"+",
                "--doDigitization", "false",
                "--configKeyValues",
                f"beam.energyPerNucleon=40;keyval.output_dir={out}",
            ]
            write_json_atomic(out/"command.json",{
                "argv": command,
                "seed":seed,
                "events":n,
                "channel":channel,
                "original_event_range_in_campaign": [
                    part*args.chunk_size, part*args.chunk_size+n,
                ],
                "env": {
                    "NA6P_MWPC_STORE_GAS_STEPS": "1",
                    "NA6P_MWPC_STORE_REFERENCE_PLANES": "1",
                },
                "notes":"Each chunk has distinct random seed; event indices "
                        "inside ROOT are local to that chunk.",
            })
            env = os.environ.copy()
            env["NA6P_MWPC_STORE_GAS_STEPS"] = "1"
            env["NA6P_MWPC_STORE_REFERENCE_PLANES"] = "1"
            log = out / "simulation.log"
            print(
                f"RUN {channel} chunk{part+1:04d}: "
                f"{n:,} events, {args.workers} workers, seed={seed}",
                flush=True,
            )
            with log.open("w",encoding="utf-8") as stream:
                result = subprocess.run(
                    command, cwd=REPO, env=env, stdout=stream,
                    stderr=subprocess.STDOUT, check=False,
                )
            if result.returncode != 0:
                raise RuntimeError(
                    f"Geant4 failed rc={result.returncode} ({log}).\n"
                    + "\n".join(log.read_text(errors="replace").splitlines()[-30:])
                )
            missing = [
                f for f in ESSENTIAL
                if not (out/f).is_file() or (out/f).stat().st_size==0
            ]
            if missing:
                raise RuntimeError(
                    f"Geant4 did not produce mandatory outputs {missing}. "
                    f"See {log}"
                )
            # ROOT checks event counts, object branches and the gas-step tree.
            validator = f'{inputs["validator"]}+("{out}",{n})'
            verify_log = out / "validation.log"
            with verify_log.open("w", encoding="utf-8") as stream:
                check = subprocess.run(
                    [root,"-l","-b","-q",validator],
                    cwd=REPO,env=env,stdout=stream,
                    stderr=subprocess.STDOUT,check=False,
                )
            valid_text = verify_log.read_text(errors="replace")
            if check.returncode != 0 or "PI_GEANT_CHUNK_PASS" not in valid_text:
                raise RuntimeError(
                    f"Geant4 output ROOT validation FAILED: {verify_log}\n"
                    + "\n".join(valid_text.splitlines()[-35:])
                )
            checksums = {name: sha256(out/name) for name in ESSENTIAL}
            write_json_atomic(done,{
                "channel":channel,
                "chunk":part+1,
                "events":n,
                "seed":seed,
                "completed_utc":dt.datetime.now(dt.timezone.utc).isoformat(),
                "output_sha256":checksums,
                "verified_gas_step_tree":True,
                "verified_virtual_muon_plane_tree":True,
            })
            count_new += 1
            print(
                f"PASS {channel} chunk{part+1:04d}; "
                f"new completed chunks {count_new}; "
                f"elapsed {(time.monotonic()-started)/60:.1f} min",
                flush=True,
            )
    print("CAMPAIGN_COMPLETE all three channels and every frozen chunk.",flush=True)
    write_json_atomic(suite / "campaign_complete.json",{
        "utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "channels":list(CHANNELS),
        "events_per_source":args.events_per_source,
        "n_chunks":chunks_per_source*len(CHANNELS),
    })


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--mode",choices=["smoke","full"],required=True)
    ap.add_argument("--events-per-source",type=int,default=None)
    ap.add_argument("--chunk-size",type=int,default=None)
    ap.add_argument("--workers",type=int,default=4)
    ap.add_argument("--hours",type=float,default=2.5,
                    help="Budget checked BETWEEN whole chunks, 0 = unlimited")
    ap.add_argument("--suite",type=Path,default=None,
                    help="Fresh custom output directory for changed settings")
    args = ap.parse_args()
    if args.events_per_source is None:
        args.events_per_source = 3 if args.mode=="smoke" else 30000
    if args.chunk_size is None:
        args.chunk_size = 3 if args.mode=="smoke" else 3000
    suite,inputs = setup(args)
    lock=(suite/".campaign.lock").open("a")
    try:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    except BlockingIOError as e:
        raise RuntimeError("Another campaign process is already running here") from e
    simulate(args,suite,inputs)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"CAMPAIGN_ERROR: {exc}", file=sys.stderr, flush=True)
        raise SystemExit(1)
