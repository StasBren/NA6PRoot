#!/usr/bin/env python3
"""SHA-protected 1000 REAL J/psi Geant4 muon crossings -> TrackHeed -> full Garfield.

EXACTLY p3_noFW: 3-mm anode pitch, NO field wires, 30-um anodes,
Ar:CO2 70:30, cathode +/-3mm, +1800V, B=0. One real Geant4 crossing
per UNIQUE generated J/psi parent event. No Poisson surrogate or LUT.

Saves full ion count/angle census and 128 unbiased ion-birth (u,v,w,t)
reservoir points PER microscopic avalanche (not full ion clouds or
cathode signals). Per-muon and per-seed wall-times, gas init, file I/O
and SHA-sealed 20-parent chunks are tracked separately.

Modes: selftest | prepare | preflight | run | status | audit
Never edit previous 4-geometry Garfield, A20/A100, D500 or full-cloud archives.
Source defaults to J/psi chunk 3, distinct from previous chunk-2 holdouts.
On interruption, never overwrite/remove unsealed last chunk.
"""
from __future__ import annotations

import argparse
import csv
import fcntl
import hashlib
import json
import math
import os
import statistics
import subprocess
import time
from collections import Counter
from pathlib import Path

from pi_2m_native_six_geometry_fold import crossing_provenance, SOURCE
from pi_prepare_muon_avalanche_b1 import FIELDS, seed, sha, validate
from pi_garfield_runtime_env import garfield_runtime_env
from pi_run_native_muon_b1 import CROSSING_COLUMNS
from pi_run_muon_lut_validation import id_of

REPO = Path(__file__).resolve().parents[2]
ROOT = REPO / "test_runs/mwpc_stage_c/pi_p3_noFW_jpsi1000_native_v1"
INPUT = ROOT / "input"
CHUNKS = ROOT / "native_chunks"
N = 1000
GROUP = 20
G = "p3_noFW"
RESERVOIR = 128
SOURCE_CHUNK_DEFAULT = 3
TAG = "NA60-JPSI-1000-ONE-GEANT-PARENT-ONE-CHORD-P3NOFW-20261105-V1"
CXX = REPO / "test/mwpc-avalanche/phaseB_native_muon_p3_noFW_jpsi_runtime_v1.cxx"
EXE = REPO / "test/mwpc-avalanche/build/mwpc_phase_b_native_p3_noFW_jpsi_bench_v1"
FILES = (
    "native_muon_crossings_input.csv", "physics.log",
    "muon_reference_crossing.csv", "muon_reference_seed.csv",
    "muon_reference_ion_birth_sample.csv",
    "muon_reference_heed_primaries.csv", "muon_reference_heed_clusters.csv",
    "muon_reference_timing.csv", "muon_reference_muon_timing.csv",
)


def must(condition, message):
    if not condition:
        raise ValueError(message)


def jread(path):
    return json.loads(Path(path).read_text())


def jwrite(path, obj):
    Path(path).write_text(
        json.dumps(obj, indent=2, sort_keys=True, allow_nan=False) + "\n"
    )


def csvread(path, required=None):
    with Path(path).open(newline="") as f:
        rd = csv.DictReader(f)
        must(rd.fieldnames is not None, "Empty CSV " + str(path))
        if required:
            must(set(required).issubset(rd.fieldnames),
                 "Missing CSV columns in " + str(path))
        for r in rd:
            yield r


def csvwrite(path, rows):
    with Path(path).open("x", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)


def check_binary():
    must(CXX.is_file() and EXE.is_file(),
         "Build ONLY the new one-geometry p3_noFW target first: " + str(EXE))
    t = CXX.read_text()
    must("RunGeometry(gas,GEOMS.front()" in t and
         '--muon-cache-physics-p3-noFW' in t and
         'p3_noFW ONLY' in t,
         "One-geometry C++ selection not found; refusing four-geometry benchmark")
    return {"cpp_sha256": sha(CXX), "exe_sha256": sha(EXE)}


def parent_key(r):
    return str(r["source"]), str(r["chunk"]), str(r["event"])


def source(nchunk):
    p, meta = crossing_provenance("Jpsi", nchunk)
    must(sha(p) == meta["csv_sha256"], "Original J/psi Geant4 CSV SHA changed")
    return p, meta


def build_selection(nchunk):
    path, meta = source(nchunk)
    must(int(nchunk) >= 3, "Use a separate J/psi chunk >=3; protected previous holdouts used chunk 2")
    chosen = {}
    candidates = 0
    for r in csvread(path, set(FIELDS) - {"source", "chunk", "momentum_GeV", "heed_seed"}):
        checked = validate(r)
        if checked is None:
            continue
        p, v = checked
        # Stay strictly in the existing nine-wire analytic Garfield domain.
        if abs(v["u_out_mm"] - v["u_in_mm"]) > 8. or \
           abs(v["w_out_mm"] - v["w_in_mm"]) > 40.:
            continue
        candidates += 1
        orig = (SOURCE["Jpsi"], int(nchunk), r["event"],
                r["hit_index"], r["track_id"], r["station"], r["det_id"])
        rec = {k: r[k] for k in FIELDS
               if k not in ("source", "chunk", "momentum_GeV", "heed_seed")}
        rec.update({
            "source": SOURCE["Jpsi"], "chunk": int(nchunk),
            "momentum_GeV": f"{p:.12g}",
            "heed_seed": seed(TAG, "heed", *orig)
        })
        pk = parent_key(rec)
        # If one generated J/psi event has several crossings, take exactly
        # one via a frozen SHA-random score independent of Garfield physics.
        rank = (hashlib.sha256((TAG + "|crossing|" +
                                "|".join(map(str, orig))).encode()).hexdigest(),
                tuple(map(str, orig)))
        if pk not in chosen or rank < chosen[pk][0]:
            chosen[pk] = (rank, rec)
    must(len(chosen) >= N,
         f"Only {len(chosen)} eligible independent J/psi parents in original chunk {nchunk}; "
         "no files created. Try a different frozen chunk >=3.")
    # SHA-ranked sample across all unique eligible parents.
    eligible = []
    for pk, (key, rec) in chosen.items():
        h = hashlib.sha256((TAG + "|parent|" + "|".join(pk)).encode()).hexdigest()
        eligible.append((h, pk, rec))
    eligible.sort()
    selected = [r for _, _, r in eligible[:N]]
    must(len(selected) == N and
         len({parent_key(r) for r in selected}) == N and
         len({id_of(r) for r in selected}) == N, "Nonindependent selected muons")
    by_station = Counter(int(r["station"]) for r in selected)
    by_sign = Counter("neg" if float(r["theta_u_in_deg"]) < 0 else "pos"
                      for r in selected)
    return selected, {
        "n_original_eligible_gas_crossings": candidates,
        "n_independent_eligible_parent_events": len(chosen),
        "selected_per_station": {str(s): by_station[s] for s in range(6)},
        "selected_per_incidence_sign": dict(by_sign),
        "selection_note": "One SHA-random eligible crossing per unique J/psi parent, then "
                          "1000 SHA-random parents from frozen native chunk. "
                          "Not a fully flux-unbiased chamber-crossing distribution."
    }, meta


def prepare(nchunk):
    must(not ROOT.exists(), "Existing benchmark root; do not overwrite: " + str(ROOT))
    selected, stats, source_meta = build_selection(nchunk)
    binary = check_binary()
    ROOT.mkdir(parents=True)
    INPUT.mkdir()
    out = INPUT / "native_muon_crossings_input.csv"
    csvwrite(out, selected)
    jwrite(INPUT / "input_manifest.json", {
        "schema": "NA60_JPSI_1000_P3NOFW_GARFIELD_BENCH_V1",
        "n_parents": N, "group_size": GROUP,
        "geometry": G, "reservoir_per_avalanche": RESERVOIR,
        "source": SOURCE["Jpsi"], "source_chunk": nchunk,
        "original_geant_csv_sha256": source_meta["csv_sha256"],
        "original_geant_root_sha256":
            source_meta["immutable_input"]["input_root_sha256"],
        "selected_input_sha256": sha(out),
        "executable_provenance_prepared": binary,
        "selection": stats,
    })
    check_input()
    print("PASS prepared exactly 1000 independent true Geant4 J/psi parent crossings", flush=True)
    print("GEOMETRY p3_noFW ONLY; 20 parents/chunk, 50 chunks total")
    print("SOURCE frozen chunk", nchunk, "eligible parent events",
          stats["n_independent_eligible_parent_events"])
    print("STATIONS", stats["selected_per_station"])
    print("ANGLE SIGNS", stats["selected_per_incidence_sign"])
    print("NO Heed or Garfield physics by prepare")


def check_input():
    meta = jread(INPUT / "input_manifest.json")
    inp = INPUT / "native_muon_crossings_input.csv"
    must(meta["schema"] == "NA60_JPSI_1000_P3NOFW_GARFIELD_BENCH_V1"
         and meta["n_parents"] == N and meta["group_size"] == GROUP
         and meta["geometry"] == G
         and sha(inp) == meta["selected_input_sha256"],
         "Frozen benchmark source changed")
    original, source_meta = source(meta["source_chunk"])
    must(source_meta["csv_sha256"] == meta["original_geant_csv_sha256"],
         "Frozen native original J/psi source changed")
    rr = list(csvread(inp, FIELDS))
    must(len(rr) == N and
         len({parent_key(r) for r in rr}) == N and
         len({id_of(r) for r in rr}) == N,
         "Wrong 1000 independent native muon identities")
    must(all(r["source"] == SOURCE["Jpsi"]
             and int(r["chunk"]) == meta["source_chunk"]
             and validate(r) is not None for r in rr),
         "Invalid selected native J/psi crossing")
    return rr, meta


def preflight():
    rr, meta = check_input()
    physics = check_binary()
    p = INPUT / "first20_native_preflight.csv"
    if not p.exists():
        csvwrite(p, rr[:GROUP])
    else:
        must([id_of(r) for r in csvread(p, FIELDS)] ==
             [id_of(r) for r in rr[:GROUP]],
             "Preflight source crossing selection changed")
    env = garfield_runtime_env()
    for command in ("--validate-input", "--validate-heed"):
        cp = subprocess.run([str(EXE), command, str(p)],
                            capture_output=True, text=True, env=env,
                            cwd=INPUT)
        must(cp.returncode == 0 and "PASS:" in cp.stdout,
             "NATIVE preflight " + command + " failed:\n" +
             (cp.stdout + "\n" + cp.stderr)[-7000:])
        print("PASS", command, cp.stdout.strip().splitlines()[-1], flush=True)
    receipt = {
        "schema": "NA60_JPSI_P3NOFW_PREFLIGHT_V1",
        "input_sha256": meta["selected_input_sha256"],
        "preflight_sha256": sha(p), "binary": physics
    }
    target = INPUT / "preflight.json"
    if target.exists():
        must(jread(target) == receipt, "Frozen native preflight has changed")
    else:
        jwrite(target, receipt)
    print("PASS 20 REAL native Heed parents checked; NO avalanches by preflight")


def checked_physics():
    r = jread(INPUT / "preflight.json")
    physics = check_binary()
    must(r["input_sha256"] == jread(INPUT / "input_manifest.json")["selected_input_sha256"]
         and r["binary"] == physics,
         "Native C++/EXE changed after preflight: do not mix physics")
    return physics


def chunks(rr):
    return [rr[i:i+GROUP] for i in range(0, N, GROUP)]


def folder(i):
    return CHUNKS / f"chunk_{i:04d}"


def validate_native(d, selected):
    """Strict one-geometry event/seed/Heed/timing conservation, source disjoint."""
    prefix = d / "muon_reference"
    base = {id_of(r): r for r in selected}
    observed = {}
    for r in csvread(d / "muon_reference_crossing.csv", CROSSING_COLUMNS):
        key = r[CROSSING_COLUMNS[0]]
        must(key in base and key not in observed and r["configuration"] == G,
             "Unexpected native crossing/geometry identity")
        src = base[key]
        must(int(r["station"]) == int(src["station"]) and
             int(r["det_id"]) == int(src["det_id"]),
             "Native station/detector mismatch")
        for k, source_k in (("path_mm", "path_length_mm"),
                            ("theta_u_in_deg", "theta_u_in_deg"),
                            ("theta_w_in_deg", "theta_w_in_deg"),
                            ("momentum_GeV", "momentum_GeV")):
            must(abs(float(r[k]) - float(src[source_k])) < 1.e-5,
                 "Native original Geant4 track changed: " + k)
        n = int(r["heed_primaries"])
        total = int(r["ion_births"])
        must(0 <= n <= 1000 and total >= 0
             and int(r["heed_failed"]) == 0
             and int(r["ion_up"]) + int(r["ion_down"]) +
                 int(r["ion_boundary"]) == total
             and int(r["ion_top"]) + int(r["ion_bottom"]) +
                 int(r["ion_side"]) == total,
             "Native Heed/ion sector conservation failed")
        observed[key] = r
    must(len(observed) == len(selected), "Missing parent/geometry native response")

    prim = Counter()
    for r in csvread(d / "muon_reference_heed_primaries.csv",
                     ("crossing_id", "primary_index", "cluster_index")):
        k = r["crossing_id"]
        must(k in base and int(r["primary_index"]) == prim[k],
             "Missing/duplicate original Heed primary index")
        prim[k] += 1
    for k, r in observed.items():
        must(prim[k] == int(r["heed_primaries"]),
             "Heed stored primary source not complete")

    seeds = {}
    seed_total = Counter()
    seed_collected = Counter()
    seed_zeros = Counter()
    seed_caps = Counter()
    seed_failed = Counter()
    for r in csvread(d / "muon_reference_seed.csv",
                     ("crossing_id", "configuration", "primary_index",
                      "ion_births", "collected_e", "successful", "capped",
                      "ion_up", "ion_down", "ion_top", "ion_bottom")):
        k = r["crossing_id"]
        idx = int(r["primary_index"])
        key = k, idx
        must(k in base and idx in range(prim[k]) and key not in seeds
             and r["configuration"] == G, "Unexpected/missing native seed")
        n = int(r["ion_births"])
        must(n >= 0 and
             int(r["ion_up"]) + int(r["ion_down"]) <= n and
             int(r["ion_top"]) + int(r["ion_bottom"]) <= n,
             "Garfield birth sectors invalid")
        seeds[key] = n
        seed_total[k] += n
        seed_collected[k] += int(r["collected_e"])
        seed_zeros[k] += int(r["collected_e"]) == 0
        seed_caps[k] += int(r["capped"])
        seed_failed[k] += int(r["successful"]) == 0
    must(len(seeds) == sum(prim.values()), "Missing independent microscopic avalanches")
    for k, r in observed.items():
        must(seed_total[k] == int(r["ion_births"]) and
             seed_collected[k] == int(r["collected_avalanche_electrons"]) and
             seed_zeros[k] == int(r["zero_seeds"]) and
             seed_caps[k] == int(r["capped_seeds"]) and
             seed_failed[k] == int(r["failed_seeds"]),
             "Native muon totals do not match all per-seed Garfield outcomes")

    timing = {}
    seed_wall_sum = 0.
    for r in csvread(d / "muon_reference_timing.csv",
                     ("crossing_id", "configuration", "primary_index",
                      "wall_seconds", "ion_births", "capped")):
        key = r["crossing_id"], int(r["primary_index"])
        must(key in seeds and key not in timing
             and r["configuration"] == G and
             int(r["ion_births"]) == seeds[key] and
             math.isfinite(float(r["wall_seconds"])) and float(r["wall_seconds"]) >= 0,
             "Per-seed native Garfield timing/identity changed")
        timing[key] = float(r["wall_seconds"])
        seed_wall_sum += timing[key]
    must(set(timing) == set(seeds), "Not all microscopic seeds have wall times")

    # Reservoir is spatially subsampled, but all ion births are counted
    # exactly. Verify 128 unbiased points and weights per nonzero avalanche.
    counts = Counter()
    for r in csvread(d / "muon_reference_ion_birth_sample.csv",
                     ("crossing_id", "configuration", "primary_index",
                      "u_mm", "v_mm", "w_mm", "t_ns",
                      "seed_ion_births_total", "reservoir_weight")):
        k = (r["crossing_id"], int(r["primary_index"]))
        must(k in seeds and r["configuration"] == G and seeds[k] > 0 and
             int(r["seed_ion_births_total"]) == seeds[k],
             "Reservoir ion source identity/weight changed")
        for c in ("u_mm", "v_mm", "w_mm", "t_ns"):
            must(math.isfinite(float(r[c])), "Invalid Garfield birth space-time")
        expected = seeds[k] / min(RESERVOIR, seeds[k])
        must(abs(float(r["reservoir_weight"]) - expected) <
             1.e-7 * max(1., expected), "Incorrect ion-birth reservoir weight")
        counts[k] += 1
        must(counts[k] <= RESERVOIR, "Ion-birth reservoir cap exceeded")
    must(all(counts[key] == min(RESERVOIR, n) for key, n in seeds.items()),
         "Missing original ion birth reservoir points")

    muon_time = {}
    muon_total_wall = 0.
    for r in csvread(d / "muon_reference_muon_timing.csv",
                     ("crossing_id", "n_Heed_primaries", "muon_total_wall_seconds",
                      "heed_wall_seconds", "garfield_and_output_wall_seconds")):
        k = r["crossing_id"]
        must(k in base and k not in muon_time and
             int(r["n_Heed_primaries"]) == prim[k],
             "Muon timing source mismatch")
        values = [float(r[f]) for f in ("muon_total_wall_seconds",
                   "heed_wall_seconds", "garfield_and_output_wall_seconds")]
        # The gap between Heed and Garfield timers also includes writing
        # the Heed clusters/primary cache. It must NOT be forced to zero:
        # full per-muon wall time legitimately exceeds these two components.
        must(all(math.isfinite(x) and x >= 0. for x in values)
             and values[0] + 0.001 >= values[1] + values[2],
             "Native per-muon wall-time decomposition invalid")
        muon_time[k] = values[0]
        muon_total_wall += values[0]
    must(len(muon_time) == len(selected), "Missing per-muon physical timings")

    return {
        "n_parents": len(selected),
        "n_trackheed_electrons": sum(prim.values()),
        "n_full_garfield_avalanches": len(seeds),
        "n_ion_births": sum(seed_total.values()),
        "n_zero_collected_seeds": sum(seed_zeros.values()),
        "n_capped_seeds": sum(seed_caps.values()),
        "n_unsuccessful_seeds": sum(seed_failed.values()),
        "n_reservoir_ion_points": sum(counts.values()),
        "sum_seed_garfield_wall_seconds": seed_wall_sum,
        "sum_muon_wall_seconds": muon_total_wall,
        "muon_time_s_min": min(muon_time.values()),
        "muon_time_s_max": max(muon_time.values())
    }


def audit_chunk(i, batch, phys, verify_files=True):
    d = folder(i)
    seal = d / "complete.json"
    must(seal.is_file(), "Unsealed or missing physics checkpoint " + str(d))
    m = jread(seal)
    must(m["schema"] == "NA60_JPSI_20PARENT_P3NOFW_SHA_V1"
         and m["chunk"] == i and m["n_parents"] == len(batch)
         and m["physics"] == phys,
         "Native checkpoint physics/source metadata changed")
    must(m["source_sha256"] == sha(d / FILES[0]) and
         [id_of(x) for x in csvread(d / FILES[0], FIELDS)] ==
         [id_of(x) for x in batch],
         "Native checkpoint source identities changed")
    if verify_files:
        must(set(m["file_sha256"]) == set(FILES),
             "Wrong set of SHA-protected native files")
        for name, digest in m["file_sha256"].items():
            must(sha(d / name) == digest, "Full-Garfield output changed: " + str(d/name))
    return m


def run(max_chunks, hours):
    rr, _ = check_input()
    physics = checked_physics()
    CHUNKS.mkdir(exist_ok=True, parents=True)
    lock = CHUNKS / ".writer.lock"
    with lock.open("a+") as fh:
        try:
            fcntl.flock(fh, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise RuntimeError("Another J/psi 1000 p3_noFW Garfield writer is ACTIVE")
        batches = chunks(rr)
        started = time.monotonic()
        completed = 0
        added = 0
        env = garfield_runtime_env()
        for i, batch in enumerate(batches):
            d = folder(i)
            if d.exists():
                audit_chunk(i, batch, physics)
                completed += 1
                continue
            if added >= max_chunks or time.monotonic() - started >= hours*3600:
                print(f"SAFE STOP SHA-sealed {completed}/{len(batches)} chunks, "
                      f"{completed*GROUP}/{N} REAL muon parents", flush=True)
                return
            d.mkdir()  # never retry/overwrite any unsealed source data
            csvwrite(d / "native_muon_crossings_input.csv", batch)
            log = d / "physics.log"
            print("FULL GARFIELD p3_noFW only chunk", i, "of", len(batches),
                  "20 real J/psi parents", flush=True)
            t0 = time.monotonic()
            with log.open("x") as f:
                cp = subprocess.run(
                    [str(EXE), "--muon-cache-physics-p3-noFW",
                     str(d / "native_muon_crossings_input.csv"),
                     str(d / "muon_reference")],
                    stdout=f, stderr=subprocess.STDOUT,
                    cwd=d, env=env)
            wall = time.monotonic() - t0
            must(cp.returncode == 0 and
                 "PASS: immutable same-Heed primary/cluster cache and full "
                 "Garfield microscopic reference emitted; p3_noFW ONLY" in
                 log.read_text(errors="replace"),
                 "Native full Garfield incomplete; KEEP unsealed chunk: " + str(d))
            stats = validate_native(d, batch)
            # Physics exact counts are audited before final SHA sealing.
            receipt = {
                "schema": "NA60_JPSI_20PARENT_P3NOFW_SHA_V1",
                "chunk": i, "n_parents": len(batch),
                "physics": physics, "source_sha256": sha(d / FILES[0]),
                "actual_native_wall_seconds": wall,
                "stats": stats,
                "file_sha256": {name: sha(d / name) for name in FILES}
            }
            # Atomic seal: status readers never observe partial JSON while
            # the physical subprocess finishes its last I/O and validation.
            receipt_tmp = d / "complete.json.tmp"
            jwrite(receipt_tmp, receipt)
            receipt_tmp.replace(d / "complete.json")
            audit_chunk(i, batch, physics)
            completed += 1
            added += 1
            print("PASS SHA chunk", i, "parents", stats["n_parents"],
                  "full avalanches", stats["n_full_garfield_avalanches"],
                  "births", stats["n_ion_births"],
                  "native wall s", round(wall, 2),
                  "mean s/muon", round(wall / GROUP, 3),
                  "capped", stats["n_capped_seeds"],
                  "failed", stats["n_unsuccessful_seeds"], flush=True)
        print("PASS ALL COMPLETE 1000 independent native J/psi parents p3_noFW ONLY", flush=True)


def scanned_status(full_audit=False):
    rr, _ = check_input()
    physics = checked_physics()
    batches = chunks(rr)
    results = []
    pending = None
    gap = False
    for i, batch in enumerate(batches):
        d = folder(i)
        if d.exists():
            must(not gap, "Physics chunks have non-contiguous checkpoint directories")
            if (d/"complete.json").is_file():
                must(pending is None, "Sealed chunk exists beyond an unsealed checkpoint")
                results.append(audit_chunk(i, batch, physics, verify_files=full_audit))
            else:
                must(pending is None, "More than one unsealed native chunk")
                pending = i
        else:
            gap = True
    return results, pending


def status(audit=False):
    results, pending = scanned_status(full_audit=audit)
    parents = sum(r["stats"]["n_parents"] for r in results)
    avalanches = sum(r["stats"]["n_full_garfield_avalanches"] for r in results)
    physics_wall = sum(r["actual_native_wall_seconds"] for r in results)
    n_heed = sum(r["stats"]["n_trackheed_electrons"] for r in results)
    n_caps = sum(r["stats"]["n_capped_seeds"] for r in results)
    n_failed = sum(r["stats"]["n_unsuccessful_seeds"] for r in results)
    print(f"SHA VERIFIED P3 NO-FW JPSI MUONS {parents} / {N}")
    print(f"SEALED CHUNKS {len(results)} / {N//GROUP}")
    print(f"ACTUAL MICROSCOPIC GARFIELD AVALANCHES {avalanches}")
    print(f"TRACKHEED ELECTRONS {n_heed}")
    print(f"ACTUAL NATIVE C++ PHYSICS+IO WALL s {physics_wall:.3f}")
    print(f"MEASURED s / MUON (including per-chunk gas init and IO) "
          f"{physics_wall/parents:.3f}" if parents else "NO PHYSICS YET")
    print(f"MICROSCOPIC SEEDS WITH CAPS {n_caps}, UNSUCCESSFUL {n_failed}")
    if parents:
        # Extrapolation assumes SOURCE-RANDOM selected parents are similar
        # and a 20-parent chunk is used on the same machine.
        forecast = N * physics_wall / parents
        remain = (N - parents) * physics_wall / parents
        print(f"PROVISIONAL 1000-MUON SEQUENTIAL WALL h {forecast/3600:.3f}")
        print(f"PROVISIONAL REMAINING WALL h {remain/3600:.3f}")
        print(f"TOTAL avg true Heed electrons / muon {n_heed/parents:.3f}")
        print(f"TOTAL avg measured seconds / avalanche "
              f"{sum(r['stats']['sum_seed_garfield_wall_seconds'] for r in results)/max(1,avalanches):.4f}")
    if pending is not None:
        print("UNSEALED CHUNK", pending, ": may be ACTIVE or INTERRUPTED; NOT COUNTED")
        print("Do not restart/delete without checking process and physics.log")
    if len(results) == len(chunks(check_input()[0])):
        print("PASS ALL 1000 P3 NOFW PARENTS SHA SEALED")
    elif pending is None:
        print("NEXT UNUSED CHUNK", len(results))
    if audit:
        print("PASS all completed checkpoint native data SHA re-audited")


def selftest():
    must(G == "p3_noFW" and N == 1000 and GROUP == 20 and
         N%GROUP == 0 and RESERVOIR == 128,
         "Wrong one-geometry benchmark configuration")
    fake = [{"source": "Jpsi", "chunk": 3, "event": str(k),
             "hit_index": str(k), "track_id": "1",
             "station": str(k%6), "det_id": "0"} for k in range(N)]
    must(len({parent_key(r) for r in fake}) == N and
         len({id_of(r) for r in fake}) == N and
         len(chunks(fake)) == 50,
         "Parent independence/chunking selftest failed")
    must(not ROOT.exists() or ROOT.is_dir(), "Invalid existing output path")
    print("PASS 1000 independent-parent identity, 50 x 20 chunk plan, one geometry")
    print("NO Geant4, TrackHeed or Garfield physics by selftest")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("mode", choices=("selftest", "prepare", "preflight",
                                    "run", "status", "audit"))
    p.add_argument("--source-chunk", type=int, default=SOURCE_CHUNK_DEFAULT)
    p.add_argument("--max-chunks", type=int, default=1)
    p.add_argument("--hours", type=float, default=6.)
    args = p.parse_args()
    must(args.max_chunks > 0 and args.hours > 0, "Invalid run budget")
    if args.mode == "selftest":
        selftest()
    elif args.mode == "prepare":
        prepare(args.source_chunk)
    elif args.mode == "preflight":
        preflight()
    elif args.mode == "run":
        run(args.max_chunks, args.hours)
    elif args.mode == "status":
        status()
    else:
        status(audit=True)


if __name__ == "__main__":
    main()
