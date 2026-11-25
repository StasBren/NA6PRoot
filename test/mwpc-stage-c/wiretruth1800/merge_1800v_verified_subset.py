#!/usr/bin/env python3
"""Merge precisely the 48 sealed 1800 V ROOT shards (5620 muons) in a NEW area.

This does not run Garfield++, does not alter the original 1800 V campaign and
does not modify its SHA-locked executables. Statistical operations and plots
are performed by the pre-existing, checksum-locked C++/ROOT merger.
"""
from __future__ import annotations
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile

from full_wiretruth1800_campaign import (
    SCHEMA, sha, require, load, save_json, copyhash,
    write_checksums, copytree_verify, timestamp,
)

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
FROZEN_REPO = Path.home() / "NA6PRoot-wiretruth1800-analysis"
FROZEN = FROZEN_REPO / "test_runs/wiretruth1800_FULL_20261123T145716Z"
MERGER = FROZEN_REPO / "test/mwpc-stage-c/build-wiretruth1800-analysis/wiretruth1800_merge_full"
FILES = ("wiretruth1800_8k.root", "channel_station_summary.csv",
         "asymmetry_denominators.csv")
EXPECTED_MISSING = (0, 5, 7, 16, 17, 27, 32)
EXPECTED_SHARDS = 48
EXPECTED_MUONS = 5620
SUBSET_SCHEMA = "NA60_1800V_WIRETRUTH_VERIFIED_SUBSET_5620_V1"


def file_sha_table(directory: Path) -> dict[str, str]:
    file = directory / "SHA256SUMS"
    require(file.is_file(), "missing source SHA256SUMS: " + str(directory))
    result = {}
    for line in file.read_text(encoding="utf-8").splitlines():
        value, sep, name = line.partition("  ")
        require(sep and len(value) == 64 and name and "/" not in name
                and name not in result and not name.startswith("."),
                "malformed source checksum entry")
        result[name] = value
    for name in (*FILES, "sealed.json"):
        require(name in result and (directory/name).is_file(),
                "incomplete sealed source shard: " + str(directory/name))
    require(sha(directory/"sealed.json") == result["sealed.json"],
            "source sealed metadata SHA mismatch")
    return result


def frozen_plan():
    require(FROZEN.is_dir(), "frozen source campaign not available")
    original = load(FROZEN/"locked/campaign.json")
    require(original["schema"] == SCHEMA and original["n_muons"] == 6460
            and original["n_shards"] == 55,
            "frozen original campaign does not match expected 6460/55")
    require(MERGER.is_file() and
            sha(MERGER) == original["merge_binary_sha256"],
            "frozen original C++/ROOT merger binary SHA changed")
    return original


def prepare(out: Path):
    out = out.resolve()
    require(str(out).startswith(str((REPO/"test_runs").resolve()) + "/"),
            "output must be a new path inside the dashboard worktree test_runs")
    require(not out.exists(), "output directory already exists; do not overwrite")
    original = frozen_plan()
    good, missing = [], []
    for item in original["shards"]:
        idx = int(item["id"])
        source = FROZEN/"shards"/f"shard_{idx:04d}"
        if not (source/"sealed.json").is_file():
            missing.append(idx)
            continue
        record = load(source/"sealed.json")
        require(record["schema"] == SCHEMA
                and record["shard"] == idx
                and record["manifest_sha256"] == item["manifest_sha256"]
                and record["n_muons"] == item["n_muons"]
                and record["analysis_binary_sha256"] == original["analysis_binary_sha256"],
                f"shard {idx} provenance is not compatible with frozen campaign")
        hashes = file_sha_table(source)
        good.append({
            "id": idx, "channel": item["channel"],
            "n_muons": item["n_muons"],
            "manifest_sha256": item["manifest_sha256"],
            "original_dir": str(source),
            "hashes": {name: hashes[name] for name in (*FILES, "sealed.json")},
        })
    require(tuple(missing) == EXPECTED_MISSING,
            "missing shards changed; review source before freezing: " + str(missing))
    muons = sum(s["n_muons"] for s in good)
    require(len(good) == EXPECTED_SHARDS and muons == EXPECTED_MUONS,
            f"unexpected subset {len(good)} shards / {muons} muons")
    subset = {
        "schema": SUBSET_SCHEMA, "created_utc": timestamp(),
        "source_campaign": str(FROZEN),
        "source_campaign_sha256": sha(FROZEN/"locked/campaign.json"),
        "source_analysis_git": original["analysis_git"],
        "source_analysis_binary_sha256": original["analysis_binary_sha256"],
        "merger_binary": str(MERGER),
        "merger_binary_sha256": original["merge_binary_sha256"],
        "excluded_missing_ids": missing, "n_shards": len(good),
        "n_muons": muons, "shards": good,
        "scope": "verified completed subset, not all 6460 muons",
        "caution": "Seven failed shards excluded; no missing events filled with zeros",
    }
    (out/"locked").mkdir(parents=True)
    (out/"pbs-logs").mkdir()
    (out/"attempts").mkdir()
    save_json(out/"locked/subset.json", subset)
    print(f"PREPARE_PASS isolated verified subset: {len(good)} shards "
          f"{muons} muons; omitted={missing}")
    print(f"OUTPUT={out}")
    print("NEXT: submit one separate PBS ROOT merge job. No Garfield jobs.")


def run(out: Path):
    require(os.environ.get("PBS_JOBID"), "run requires PBS compute node")
    scratch = Path(os.environ.get("TMPDIR", ""))
    require(str(scratch) not in ("", ".") and scratch.is_dir(),
            "PBS TMPDIR missing")
    subset = load(out/"locked/subset.json")
    original = frozen_plan()
    require(subset["schema"] == SUBSET_SCHEMA
            and subset["source_campaign_sha256"] == sha(FROZEN/"locked/campaign.json")
            and subset["merger_binary_sha256"] == original["merge_binary_sha256"],
            "partial snapshot provenance drift")
    require(not (out/"merged").exists(), "partial merged output already exists")
    with tempfile.TemporaryDirectory(prefix="w1800_subset_", dir=scratch) as tmp:
        base = Path(tmp)
        listing = base/"shard_dirs.txt"
        with listing.open("x", encoding="utf-8") as f:
            for item in subset["shards"]:
                idx = item["id"]
                source = Path(item["original_dir"])
                locked = original["shards"][idx]
                record = load(source/"sealed.json")
                require(record["manifest_sha256"] == locked["manifest_sha256"]
                        and record["n_muons"] == item["n_muons"]
                        and record["analysis_binary_sha256"] ==
                            original["analysis_binary_sha256"],
                        f"shard {idx} changed after selection")
                observed = file_sha_table(source)
                for name, wanted in item["hashes"].items():
                    require(observed[name] == wanted, f"shard {idx} checksum list changed")
                target = base/"inputs"/f"shard_{idx:04d}"
                target.mkdir(parents=True)
                for name in FILES:
                    copyhash(source/name, target/name, item["hashes"][name])
                f.write(str(target)+"\n")
                print(f"SUBSET_SOURCE_SHA_PASS id={idx} muons={item['n_muons']}",
                      flush=True)
        output = base/"output"
        output.mkdir()
        print(f"SUBSET_ROOT_MERGE_BEGIN shards={subset['n_shards']} "
              f"muons={subset['n_muons']}", flush=True)
        subprocess.run([str(MERGER), str(listing), str(output),
                        str(subset["n_muons"])], check=True)
        shutil.copy2(FROZEN/"locked/campaign.json", output/"frozen_campaign.json")
        shutil.copy2(out/"locked/subset.json", output/"verified_subset_manifest.json")
        save_json(output/"sealed.json", {
            "schema": SUBSET_SCHEMA, "n_shards": subset["n_shards"],
            "n_muons": subset["n_muons"], "source_campaign_sha256":
                subset["source_campaign_sha256"],
            "merger_binary_sha256": subset["merger_binary_sha256"],
            "pbs_job_id": os.environ["PBS_JOBID"], "sealed_utc": timestamp()
        })
        write_checksums(output)
        target = out/"attempts"/f"merged_staging_{os.getpid()}"
        copytree_verify(output, target)
        require(not (out/"merged").exists(),
                "another merge output appeared")
        target.rename(out/"merged")
        print(f"SUBSET_MERGE_PASS muons={subset['n_muons']} "
              f"output={out/'merged'}", flush=True)


def status(out: Path):
    subset = load(out/"locked/subset.json")
    print(f"FROZEN_SUBSET={subset['n_muons']} muons "
          f"{subset['n_shards']} shards missing={subset['excluded_missing_ids']}")
    merged = out/"merged"
    if (merged/"SHA256SUMS").is_file():
        print(f"MERGED_READY={merged}")
        print("ROOT, CSV, PNG and provenance checksums exist; run sha256sum -c "
              "for independent final audit")
    else:
        print("MERGED_NOT_YET_PRESENT; check the dedicated PBS job logs")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("mode", choices=("prepare", "run", "status"))
    p.add_argument("--output", required=True, type=Path)
    args = p.parse_args()
    if args.mode == "prepare":
        prepare(args.output)
    elif args.mode == "run":
        run(args.output)
    else:
        status(args.output)


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        raise SystemExit("SUBSET_FAIL: " + str(e))
