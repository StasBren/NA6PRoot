#!/usr/bin/env python3
"""Read-only, full SHA-256 audit of the 3 x 10k Penning-ON original-2M campaigns.

Accept a separate recovery only when the original result is not successful.
Verify every original Geant4 input chunk, all seven native Garfield++ CSV
outputs, the native physics report, frozen source/selection manifest, and
binary fingerprint. All 500 distinct 20-muon chunks must pass, with zero
capped avalanches. Never modifies an existing campaign.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
import time
from pathlib import Path

from pi_penning057_original2m_cluster10k_v1 import (
    GROUP, N, OUTPUTS, SCHEMA, check_input, digest, metadata, require,
)

CHANNELS = ("Jpsi", "Phi", "Omega")
EXPECTED_SEAL_NAMES = set(OUTPUTS) | {"verified_report.json"}


def sealed_files(output: Path) -> int:
    seal_path = output / "outputs.sha256"
    require(seal_path.is_file(), f"Missing output checksum manifest: {seal_path}")
    lines = seal_path.read_text().splitlines()
    hashes = {}
    for line in lines:
        require(len(line) >= 67 and line[64:66] == "  " and
                re.fullmatch(r"[0-9a-f]{64}", line[:64]) is not None,
                f"Invalid checksum record in {seal_path}: {line[:90]!r}")
        name = line[66:]
        require(name in EXPECTED_SEAL_NAMES and name not in hashes,
                f"Unexpected or repeated sealed filename {name!r} in {seal_path}")
        hashes[name] = line[:64]
    require(set(hashes) == EXPECTED_SEAL_NAMES,
            f"Missing CSV or report files from SHA seal: {seal_path}")
    count = 0
    for name, sha in hashes.items():
        actual = output / name
        require(actual.is_file(), f"Missing sealed output: {actual}")
        require(digest(actual) == sha, f"SHA-256 mismatch: {actual}")
        count += 1
    return count


def audit_channel(root: Path, channel: str, report_every: int) -> dict:
    campaign = root / f"original2m_{channel.lower()}_10000_v1"
    meta = metadata(campaign)
    require(meta["channel"] == channel and meta["n_chunks"] == 500,
            f"Wrong campaign channel/chunk count: {channel}")
    require(meta["schema"] == SCHEMA and meta["n_muons"] == N and
            meta["muons_per_chunk"] == GROUP,
            f"Unexpected physics campaign schema for {channel}")
    require(digest(campaign / "frozen_original2m_10000_input.csv") ==
            meta["frozen_source_sha256"], f"Frozen source changed: {channel}")
    require(digest(campaign / "selection_manifest.json") ==
            meta["original_2m_selection_manifest_sha256"],
            f"Original Geant4 selection provenance changed: {channel}")
    binary = Path(meta["binary_path"])
    require(binary.is_file() and digest(binary) == meta["binary_sha256"],
            f"Penning-ON executable no longer has frozen SHA: {binary}")
    recovered, failures_preserved, sealed_count = [], [], 0
    muons = avalanches = caps = births = collected = 0
    all_wall = 0
    for index in range(500):
        key = f"{index:04d}"
        check_input(campaign, index)  # SHA + distinct native Geant4 crossing IDs
        main = campaign / "results" / f"chunk_{key}"
        recovery = campaign / "recovery_runtime_v1" / f"chunk_{key}"

        def successful(directory: Path) -> bool:
            return (directory.is_dir() and
                    (directory / "exit_status.txt").is_file() and
                    (directory / "exit_status.txt").read_text().strip() == "0" and
                    (directory / "verified_report.json").is_file())

        main_ok = successful(main)
        recovery_ok = successful(recovery)
        require(not (main_ok and recovery_ok),
                f"Duplicate successful original/recovery: {channel}:{key}")
        if recovery.is_dir():
            require(recovery_ok and not main_ok,
                    f"Invalid/unfinished recovery: {channel}:{key}")
            output = recovery
            recovered.append(index)
            if main.is_dir():
                failures_preserved.append(index)
        else:
            require(main_ok, f"Unfinished, failed or missing chunk: {channel}:{key}")
            output = main

        report = json.loads((output / "verified_report.json").read_text())
        require(report.get("schema") == SCHEMA and
                report.get("chunk_index") == index and
                report.get("n_muons") == GROUP and
                report.get("input_sha256") == meta["chunks"][key]["sha256"] and
                report.get("binary_sha256") == meta["binary_sha256"],
                f"Physics/input provenance verification failed: {channel}:{key}")
        require(isinstance(report.get("n_true_heed_primaries"), int)
                and report["n_true_heed_primaries"] > 0 and
                isinstance(report.get("n_capped_avalanches"), int) and
                report["n_capped_avalanches"] == 0,
                f"Missing Heed primaries or capped avalanche: {channel}:{key}")

        sealed_count += sealed_files(output)
        wall = output / "wall_seconds.txt"
        require(wall.is_file(), f"Missing worker runtime: {channel}:{key}")
        all_wall += float(wall.read_text().strip())
        muons += GROUP
        avalanches += report["n_true_heed_primaries"]
        caps += report["n_capped_avalanches"]
        births += report["ion_births"]
        collected += report["collected_electrons"]
        if (index + 1) % report_every == 0:
            print(f"{channel}: SHA checked {index+1}/500 chunks "
                  f"({sealed_count} sealed files)", flush=True)

    require(muons == N and caps == 0 and sealed_count == 500 * len(EXPECTED_SEAL_NAMES),
            f"Not all chunks passed strict audit: {channel}")
    result = {
        "channel": channel,
        "muons": muons,
        "avalanches": avalanches,
        "capped": caps,
        "ion_births": births,
        "collected_electrons": collected,
        "worker_seconds": all_wall,
        "recovered_indices": recovered,
        "preserved_original_failures": failures_preserved,
        "sha_sealed_files": sealed_count,
    }
    print(f"PASS {channel}: {muons} distinct real Geant4 muons, "
          f"{avalanches} microscopic Garfield avalanches, 0 caps; "
          f"{len(recovered)} recoveries {recovered}; all {sealed_count} "
          f"native/report files SHA checked", flush=True)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--root", type=Path,
        default=Path("/storage/cgrp/stanis/mwpc_penning057"))
    parser.add_argument("--channel", choices=(*CHANNELS, "All"), default="All")
    parser.add_argument("--report-every", type=int, default=25)
    args = parser.parse_args()
    require(args.report_every >= 1, "--report-every must be >=1")
    channels = CHANNELS if args.channel == "All" else (args.channel,)
    start = time.monotonic()
    results = [audit_channel(args.root, c, args.report_every) for c in channels]
    if len(results) == 3:
        require(sum(r["muons"] for r in results) == 30000 and
                sum(r["capped"] for r in results) == 0,
                "Complete campaign not verified")
        print("PASS ALL THREE CHANNELS: 30,000 real muons; "
              "1,500 distinct verified chunks, zero caps; "
              "all native data and recovery results full SHA-256 audited.")
        print("Total microscopic avalanches:",
              sum(r["avalanches"] for r in results))
        print("Total ion births:", sum(r["ion_births"] for r in results))
        print("Total recovery chunks:", sum(len(r["recovered_indices"]) for r in results))
    print(f"Audit elapsed wall time: {time.monotonic()-start:.1f} seconds")


if __name__ == "__main__":
    try:
        main()
    except (RuntimeError, AssertionError, ValueError, OSError) as exc:
        print(f"FAIL STRICT SCIENTIFIC AUDIT: {exc}", file=sys.stderr)
        sys.exit(1)
