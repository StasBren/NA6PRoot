#!/usr/bin/env python3
"""Non-overwriting 10,000-real-muon Penning-ON PBS array campaign.

Only standard-library Python. Source is a NEW verified selection from the
original 2M-parent/channel Geant4 production (Jpsi, Phi or Omega).
gas-crossing CSV: ONE real crossing per distinct generated parent.
500 independent 20-muon chunks, one full TrackHeed + Garfield microscopic
avalanche for each genuine Heed primary, Penning r=0.57 enabled in C++.
No LUT and no Penning-OFF calculations. All runtime products in /storage.
"""
import argparse
import csv
import hashlib
import json
import shutil
import subprocess
from collections import defaultdict
from pathlib import Path

N, GROUP = 10000, 20
SCHEMA = "NA60_P3NOFW_ORIGINAL2M_10K_PENNING057_PBS_V1"
CHANNEL_SOURCE = {"Jpsi": "jpsi40-bore", "Phi": "phi40-bore", "Omega": "omega40-bore"}
HEADER = (
    "source", "chunk", "event", "hit_index", "track_id", "station", "det_id",
    "pdg", "u_in_mm", "v_in_mm", "w_in_mm", "u_out_mm", "v_out_mm", "w_out_mm",
    "theta_u_in_deg", "theta_w_in_deg", "path_length_mm", "edep_eV",
    "pu_in_GeV", "pv_in_GeV", "pw_in_GeV", "momentum_GeV", "heed_seed"
)
OUTPUTS = (
    "penning_crossing.csv", "penning_seed.csv", "penning_ion_birth_sample.csv",
    "penning_heed_primaries.csv", "penning_heed_clusters.csv",
    "penning_timing.csv", "penning_muon_timing.csv"
)


def require(value, message):
    if not value:
        raise RuntimeError(message)


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for buf in iter(lambda: f.read(1 << 20), b""):
            h.update(buf)
    return h.hexdigest()


def readcsv(path):
    with Path(path).open(newline="") as f:
        reader = csv.DictReader(f)
        require(reader.fieldnames is not None, "Empty CSV: " + str(path))
        require(len(reader.fieldnames) == len(set(reader.fieldnames)),
                "Duplicate CSV column: " + str(path))
        return reader.fieldnames, list(reader)


def parent_key(row):
    return (row["source"], int(row["chunk"]), int(row["event"]))


def crossing_id(row):
    return "|".join([row["source"]] + [str(int(row[k])) for k in (
        "chunk", "event", "hit_index", "track_id", "station", "det_id"
    )])


def source_rows(path, channel):
    columns, rows = readcsv(path)
    require(tuple(columns) == HEADER, "Input columns differ from native C++ format")
    require(len(rows) == N, f"Expected {N} real crossings; got {len(rows)}")
    parents, crossings = set(), set()
    for row in rows:
        p, c = parent_key(row), crossing_id(row)
        require(p not in parents and c not in crossings,
                "Duplicate J/psi parent or crossing: " + str(p))
        parents.add(p)
        crossings.add(c)
        require(row["source"] == CHANNEL_SOURCE[channel],
                "Wrong original-Geant4 source in selected crossing CSV")
        require(abs(int(row["pdg"])) == 13 and int(row["heed_seed"]) > 0,
                "Invalid muon or Heed seed")
        require(0 <= int(row["station"]) <= 5,
                "Invalid station")
    return rows


def metadata(campaign):
    path = campaign / "manifest.json"
    require(path.is_file(), "Campaign missing: " + str(path))
    data = json.loads(path.read_text())
    require(data["schema"] == SCHEMA and data["n_muons"] == N and
            data["muons_per_chunk"] == GROUP, "Wrong campaign manifest")
    return data


def prepare(source, campaign, binary, channel, selection_manifest):
    source = source.resolve(strict=True)
    binary = binary.resolve(strict=True)
    require(binary.is_file() and binary.stat().st_mode & 0o111,
            "Penning high-Heed executable not found or not executable")
    require(not campaign.exists(), "Campaign already exists; refusing overwrite")
    require(channel in CHANNEL_SOURCE, "Source must be Jpsi, Phi or Omega")
    require(selection_manifest.is_file(), "Original 2M source selection manifest missing")
    selection = json.loads(selection_manifest.read_text())
    require(selection["schema"]=="NA60_PENNING057_ORIGINAL_2M_SOURCE_SELECTION_V1" and
            selection["channel"]==channel and
            selection["source"]==CHANNEL_SOURCE[channel] and
            selection["n_selected"]==N and
            selection["n_source_chunks"]==40 and
            selection["n_per_chunk"]==250,
            "Not the correct provenance-checked original 2M event selection")
    rows = source_rows(source,channel)  # validate before creating output.
    source_sha, binary_sha = digest(source), digest(binary)
    require(source_sha == selection["input_sha256"],
            "Selected CSV does not match original 2M parent source manifest")
    campaign.mkdir(parents=True, exist_ok=False)
    (campaign / "chunks").mkdir()
    (campaign / "results").mkdir()
    (campaign / "logs").mkdir()
    saved = campaign / "frozen_original2m_10000_input.csv"
    shutil.copyfile(source, saved)
    require(digest(saved) == source_sha, "Source copy changed")
    shutil.copyfile(selection_manifest, campaign / "selection_manifest.json")
    require(digest(campaign / "selection_manifest.json") == digest(selection_manifest),
            "Original Geant4 parent selection manifest changed")
    chunkmeta = {}
    for i in range(N // GROUP):
        items = rows[GROUP * i:GROUP * (i + 1)]
        target = campaign / "chunks" / f"chunk_{i:04d}.csv"
        with target.open("x", newline="") as f:
            w = csv.DictWriter(f, fieldnames=HEADER, lineterminator="\n")
            w.writeheader()
            w.writerows(items)
        chunkmeta[f"{i:04d}"] = {
            "sha256": digest(target),
            "crossing_ids": [crossing_id(x) for x in items],
            "parent_keys": [list(parent_key(x)) for x in items],
        }
    manifest = {
        "schema": SCHEMA,
        "channel": channel,
        "original_2m_selection_manifest_sha256": digest(selection_manifest),
        "n_muons": N,
        "muons_per_chunk": GROUP,
        "n_chunks": N // GROUP,
        "original_source_path": str(source),
        "frozen_source_sha256": source_sha,
        "binary_path": str(binary),
        "binary_sha256": binary_sha,
        "physics": "p3_noFW Ar/CO2(70/30) 1800V r_Penning=0.57; full TrackHeed and Garfield++",
        "chunks": chunkmeta,
    }
    (campaign / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    (campaign / "binary.sha256").write_text(f"{binary_sha}  {binary}\n")
    print(f"PASS: 10,000 original {channel} parents, 40x250, in 500 distinct 20-muon inputs")
    print("Source SHA256:", source_sha)
    print("Penning executable SHA256:", binary_sha)
    print("Campaign:", campaign)


def check_input(campaign, index, binary=None):
    meta = metadata(campaign)
    require(0 <= index < N // GROUP, "Index must be 0..499")
    key = f"{index:04d}"
    file = campaign / "chunks" / f"chunk_{key}.csv"
    require(digest(file) == meta["chunks"][key]["sha256"],
            "Frozen chunk input SHA mismatch")
    if binary is not None:
        require(digest(binary) == meta["binary_sha256"],
                "Executable SHA changed since campaign preparation")
    columns, rows = readcsv(file)
    require(tuple(columns) == HEADER and len(rows) == GROUP,
            "Frozen chunk row/header mismatch")
    require([crossing_id(r) for r in rows] ==
            meta["chunks"][key]["crossing_ids"],
            "Frozen crossing identities changed")
    return meta, rows


def verify_chunk(campaign, index, workdir, binary):
    meta, rows = check_input(campaign, index, binary)
    expected = [crossing_id(r) for r in rows]
    require(all((workdir / p).is_file() for p in OUTPUTS),
            "Missing native CSV output(s)")
    _, crossings = readcsv(workdir / "penning_crossing.csv")
    _, times = readcsv(workdir / "penning_muon_timing.csv")
    _, seeds = readcsv(workdir / "penning_seed.csv")
    _, primaries = readcsv(workdir / "penning_heed_primaries.csv")
    _, timing_rows = readcsv(workdir / "penning_timing.csv")
    require(len(crossings) == len(times) == GROUP,
            "Native C++ did not complete exactly 20 real muons")
    require([r["source|chunk|event|hit_index|track_id|station|det_id"]
             for r in crossings] == expected, "Crossing order/identity differs")
    require([r["crossing_id"] for r in times] == expected,
            "Muon timing identities differ")
    require(len(primaries) == len(timing_rows) == len(seeds) and len(seeds) > 0,
            "TrackHeed primaries, Garfield seeds and per-seed timings differ")
    s_collect, s_birth = defaultdict(int), defaultdict(int)
    seen = set()
    capped = collected = births = zeros = 0
    for row in seeds:
        k = (row["crossing_id"], int(row["primary_index"]))
        require(k[0] in expected and k not in seen,
                "Duplicate/unknown Garfield primary electron")
        seen.add(k)
        c = int(row["collected_e"])
        b = int(row["ion_births"])
        require(c >= 0 and b >= 0, "Negative physical electron or ion count")
        s_collect[k[0]] += c
        s_birth[k[0]] += b
        collected += c
        births += b
        zeros += (c == 0)
        capped += int(row["capped"])
    for row in crossings:
        k = row["source|chunk|event|hit_index|track_id|station|det_id"]
        require(int(row["collected_avalanche_electrons"]) == s_collect[k] and
                int(row["ion_births"]) == s_birth[k],
                "Per-muon electron/ion conservation mismatch")
    report = {
        "schema": SCHEMA,
        "chunk_index": index,
        "input_sha256": meta["chunks"][f"{index:04d}"]["sha256"],
        "binary_sha256": meta["binary_sha256"],
        "n_muons": GROUP,
        "n_true_heed_primaries": len(seeds),
        "n_capped_avalanches": capped,
        "n_zero_collected": zeros,
        "collected_electrons": collected,
        "ion_births": births,
        "physics_note": "Capped avalanches, if any, are censored, not physical gain values",
    }
    (workdir / "verified_report.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(f"PASS verified chunk {index:04d}: 20 real J/psi muons, "
          f"{len(seeds)} microscopic avalanches, "
          f"{collected} anode electrons, {births} ion births, {capped} capped",
          flush=True)


def status(campaign, audit=False):
    meta = metadata(campaign)
    require(digest(campaign / "frozen_original2m_10000_input.csv") ==
            meta["frozen_source_sha256"], "Frozen 10k source SHA mismatch")
    done, failed, active, missing = [], [], [], 0
    for i in range(N // GROUP):
        key = f"{i:04d}"
        out = campaign / "results" / f"chunk_{key}"
        if not out.exists():
            missing += 1
            continue
        exitfile = out / "exit_status.txt"
        if not exitfile.is_file():
            active.append(i)
            continue
        if exitfile.read_text().strip() != "0" or not (out / "verified_report.json").is_file():
            failed.append(i)
            continue
        report = json.loads((out / "verified_report.json").read_text())
        require(report["chunk_index"] == i and report["schema"] == SCHEMA and
                report["binary_sha256"] == meta["binary_sha256"] and
                report["input_sha256"] == meta["chunks"][key]["sha256"] and
                report["n_muons"] == GROUP,
                "Completed task report provenance mismatch: " + key)
        if audit:
            require(digest(campaign / "chunks" / f"chunk_{key}.csv") ==
                    meta["chunks"][key]["sha256"], "Changed chunk input")
            cp = subprocess.run(["sha256sum", "-c", "outputs.sha256"], cwd=out,
                                capture_output=True, text=True)
            require(cp.returncode == 0, "Output SHA audit FAILED: " + key +
                    "\n" + cp.stdout[-1000:] + cp.stderr[-1000:])
        wallfile = out / "wall_seconds.txt"
        require(wallfile.is_file(), "Missing wall seconds: " + key)
        done.append((report, float(wallfile.read_text().strip())))
    n_mu = sum(r["n_muons"] for r, _ in done)
    n_seed = sum(r["n_true_heed_primaries"] for r, _ in done)
    n_cap = sum(r["n_capped_avalanches"] for r, _ in done)
    n_birth = sum(r["ion_births"] for r, _ in done)
    n_coll = sum(r["collected_electrons"] for r, _ in done)
    core_wall = sum(sec for _, sec in done)
    print(f"{meta['channel']}: COMPLETED {len(done)}/500 jobs = {n_mu}/10000 real muons")
    print(f"IN PROGRESS / INCOMPLETE: {len(active)}; FAILED: {len(failed)}; NOT STARTED: {missing}")
    print(f"True Garfield avalanches {n_seed}; capped {n_cap}; "
          f"ion births {n_birth}; collected electrons {n_coll}")
    if done:
        print(f"SUM worker wall seconds {core_wall:.0f}; "
              f"mean {core_wall/n_mu:.2f} seconds per real muon (includes startup/I/O)")
        print("Note: sum of per-job wall time is NOT elapsed time of concurrent PBS array.")
    if failed:
        print("Failed indices (inspect output, DO NOT overwrite):", failed[:25])
    if active:
        print("Incomplete/running indices:", active[:25])
    if audit:
        require(len(done) == 500 and not failed and not active and n_cap == 0,
                "Complete scientific audit requires 500 successful chunks and zero caps")
        print("PASS ALL 10000 REAL PENNING-ON MUONS SHA-AUDITED, ZERO CAPS")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("command", choices=("prepare", "check-input", "verify-chunk",
                                         "status", "audit"))
    ap.add_argument("--campaign", required=True, type=Path)
    ap.add_argument("--source", type=Path)
    ap.add_argument("--selection-manifest", type=Path)
    ap.add_argument("--channel", choices=sorted(CHANNEL_SOURCE))
    ap.add_argument("--binary", type=Path)
    ap.add_argument("--index", type=int)
    ap.add_argument("--workdir", type=Path)
    a = ap.parse_args()
    campaign = a.campaign.resolve()
    if a.command == "prepare":
        require(a.source and a.binary and a.channel and a.selection_manifest,
                "prepare needs --source --selection-manifest --channel --binary")
        prepare(a.source,campaign,a.binary,a.channel,a.selection_manifest)
    elif a.command in ("check-input", "verify-chunk"):
        require(a.index is not None and a.binary,
                "Need --index and --binary")
        if a.command == "check-input":
            check_input(campaign, a.index, a.binary)
            print(f"PASS frozen input and binary for chunk {a.index}", flush=True)
        else:
            require(a.workdir, "verify-chunk needs --workdir")
            verify_chunk(campaign, a.index, a.workdir, a.binary)
    else:
        status(campaign, a.command == "audit")


if __name__ == "__main__":
    main()
