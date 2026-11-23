#!/usr/bin/env python3
"""Read-only snapshot, verified PBS sharding, status and exact ROOT merge orchestration.
NO numerical physics/statistical fits in Python: all physics and merged
statistics are produced by compiled C++/ROOT applications.
No Garfield++ jobs are run; original source and archives never modified.
"""
from __future__ import annotations
import argparse
import csv
import datetime as dt
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time

SCHEMA = "NA60_1800V_FULL_WIRE_ASYM_SHARDED_V1"
SOURCE_SCHEMA = "NA60_WIRETRUTH30K_PENNING057_PRODUCTION_V1"
ORIGINAL_GIT = "599e29522edb6fcec9dfcca2904308cb3186edce"
SOURCE_ROOT = Path("/storage/cgrp/stanis/mwpc_penning057_wiretruth30k_v1")
ORIGINAL_ROOT = Path("/storage/cgrp/stanis/mwpc_penning057")
CHANNELS = ("Jpsi", "Phi", "Omega")
MANIFEST_HEADER = [
    "channel", "chunk", "original_input", "wire_muon",
    "wire_seed", "crossing", "seed",
]
INPUT_FILES = {
    "wire_muon": "microscopic_wire_muon.csv",
    "wire_seed": "microscopic_wire_seed.csv",
    "crossing": "microscopic_crossing.csv",
    "seed": "microscopic_seed.csv",
}
REPO = Path(__file__).resolve().parents[3]
ANALYZER = REPO / "test/mwpc-stage-c/build-wiretruth1800-analysis/wiretruth1800_distributions"
MERGER = REPO / "test/mwpc-stage-c/build-wiretruth1800-analysis/wiretruth1800_merge_full"
SHARD_OUTPUT = ("wiretruth1800_8k.root", "channel_station_summary.csv",
                "asymmetry_denominators.csv", "analysis_readme.txt")
MERGED_OUTPUT = ("wiretruth1800_full.root", "channel_station_summary.csv",
                 "asymmetry_denominators.csv", "merged_root_qa.txt",
                 "frozen_campaign.json")


def require(ok: bool, message: str) -> None:
    if not ok:
        raise RuntimeError("1800V FULL CAMPAIGN: " + message)


def timestamp() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")


def stage(message: str) -> None:
    print(f"W1800_FULL_STAGE {timestamp()} {message}", flush=True)


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as inp:
        for chunk in iter(lambda: inp.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def copyhash(source: Path, dest: Path, wanted: str, size: int | None = None) -> None:
    require(source.is_file(), f"missing original source {source}")
    dest.parent.mkdir(parents=True, exist_ok=True)
    h = hashlib.sha256()
    n = 0
    with source.open("rb") as inp, dest.open("xb") as out:
        while True:
            b = inp.read(1024 * 1024)
            if not b:
                break
            h.update(b)
            n += len(b)
            out.write(b)
    require(h.hexdigest() == wanted,
            f"sha256 mismatch for immutable original {source}")
    if size is not None:
        require(n == size, f"file size mismatch for original {source}")


def load(path: Path):
    with path.open(encoding="utf-8") as inp:
        return json.load(inp)


def save_json(path: Path, obj) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8") as out:
        json.dump(obj, out, indent=2, sort_keys=True)
        out.write("\n")
        out.flush()
        os.fsync(out.fileno())


def git() -> str:
    return subprocess.check_output(
        ["git", "-C", str(REPO), "rev-parse", "HEAD"], text=True
    ).strip()


def original_locks():
    selpath = SOURCE_ROOT / "locked/selection.json"
    binpath = SOURCE_ROOT / "locked/binary.json"
    selection, binary = load(selpath), load(binpath)
    require(selection["schema"] == SOURCE_SCHEMA and binary["schema"] == SOURCE_SCHEMA,
            "wrong 1800V production schema")
    require(selection["git_commit"] == binary["git_commit"] == ORIGINAL_GIT,
            "wrong original physics Git lock")
    require(selection["cxx_sha256"] == binary["cxx_sha256"],
            "original physics source lock changed")
    require(selection["physics"]["anode_voltage_V"] == 1800 and
            selection["physics"]["pitch_mm"] == 3 and
            selection["physics"]["penning_r"] == .57,
            "original 1800V field geometry differs from locked baseline")
    require(selection["input_strategy"] ==
            "direct_read_original2m_source_chunk_lazy_sha256_verify",
            "different Geant4 source selection strategy")
    return selection, binary, sha(selpath), sha(binpath)


def source_path(selection, channel, index):
    name = channel.lower()
    source = Path(selection["source_campaigns"][channel]["original_chunk_dir"])
    expected_dir = ORIGINAL_ROOT / f"original2m_{name}_10000_v1/chunks"
    require(source == expected_dir, "unknown original frozen Geant4 source directory")
    return source / f"chunk_{index:04d}.csv"


def completed_receipt(selection, binary, channel, index):
    d = SOURCE_ROOT / "results" / channel.lower() / f"chunk_{index:04d}"
    p = d / "verified_report.json"
    if not p.is_file():
        return None
    x = load(p)
    require(x["schema"] == SOURCE_SCHEMA and x["channel"] == channel and
            x["chunk"] == index, f"bad channel/index/schema in completed receipt {p}")
    require(x["git_commit"] == ORIGINAL_GIT and
            x["binary_sha256"] == binary["binary_sha256"] and
            x["cxx_sha256"] == selection["cxx_sha256"],
            f"source physics provenance mismatch {p}")
    require(x["input_sha256"] ==
            selection["input_sha256"][f"{channel}:{index:04d}"],
            f"Geant4 source mismatch in completed receipt {p}")
    require(x["counts"]["muons"] == 20 and
            x["counts"]["capped"] == 0 and
            x["counts"]["failed_return"] == 0,
            f"invalid incomplete 20-muon physical sample in receipt {p}")
    for key in INPUT_FILES:
        rec = x["files"][key]
        require(len(rec["sha256"]) == 64 and int(rec["bytes"]) > 0,
                f"invalid SHA or byte count in {p}: {key}")
    return x, p, d


def job_list(root: Path):
    plan = load(root / "locked/campaign.json")
    require(plan["schema"] == SCHEMA, "wrong campaign manifest schema")
    require(plan["analysis_git"] == git(),
            "analysis source HEAD differs from frozen campaign; do not mix code versions")
    require(plan["analysis_binary_sha256"] == sha(ANALYZER),
            "C++ executable was modified after snapshot preparation")
    require(plan["merge_binary_sha256"] == sha(MERGER),
            "C++ ROOT merger was modified after snapshot preparation")
    selection, binary, sel_sha, bin_sha = original_locks()
    require(plan["source_selection_sha256"] == sel_sha and
            plan["source_binary_lock_sha256"] == bin_sha,
            "original frozen source locks changed")
    shards = plan["shards"]
    require(len(shards) == plan["n_shards"] and
            plan["n_muons"] == 20 * plan["n_chunks"],
            "bad campaign total counts")
    seen = set()
    for i, shard in enumerate(shards):
        require(shard["id"] == i and shard["channel"] in CHANNELS and
                len(shard["indices"]) >= 1, "malformed shard plan")
        manifest = root / "locked/shards" / f"shard_{i:04d}.tsv"
        require(sha(manifest) == shard["manifest_sha256"],
                "shard manifest changed after preparation")
        for index in shard["indices"]:
            require((shard["channel"], index) not in seen,
                    "duplicate real crossing input chunk in full plan")
            seen.add((shard["channel"], index))
    require(len(seen) == plan["n_chunks"], "incomplete or duplicated full plan")
    return plan, selection, binary


def make_manifest(path: Path, rows) -> None:
    with path.open("x", encoding="utf-8", newline="") as f:
        wr = csv.writer(f, delimiter="\t", lineterminator="\n")
        wr.writerow(MANIFEST_HEADER)
        wr.writerows(rows)
        f.flush()
        os.fsync(f.fileno())


def manifest_rows(path: Path):
    with path.open(encoding="utf-8", newline="") as inp:
        r = csv.DictReader(inp, delimiter="\t")
        require(r.fieldnames == MANIFEST_HEADER, "invalid 1800V input manifest schema")
        rows = list(r)
    require(len(rows) > 0, "empty shard input manifest")
    return rows


def prepare(root: Path, per_shard: int) -> None:
    require(1 <= per_shard <= 20, "per-shard chunk count should be 1..20")
    require(not root.exists(), "campaign output already exists: " + str(root))
    require(ANALYZER.is_file() and MERGER.is_file(),
            "compile both ROOT/C++ executables using full1800.sh build")
    require(root.is_relative_to(REPO / "test_runs"),
            "full campaign must live under independent worktree test_runs")
    stage("PREPARE_BEGIN enumerate SHA-receipt-backed completed original 1800V chunks")
    selection, binary, sel_sha, bin_sha = original_locks()
    eligible = {c: [] for c in CHANNELS}
    for channel in CHANNELS:
        for idx in range(500):
            item = completed_receipt(selection, binary, channel, idx)
            if item is not None:
                x, receipt_path, result = item
                eligible[channel].append((idx, x, receipt_path, result))
        stage(f"RECEIPT_ENUMERATED channel={channel} verified20_chunks={len(eligible[channel])}")
    n_chunks = sum(len(a) for a in eligible.values())
    require(n_chunks > 0, "no complete original 1800V chunks at snapshot time")
    root.mkdir(parents=True, exist_ok=False)
    (root / "locked/shards").mkdir(parents=True)
    (root / "shards").mkdir()
    (root / "attempts").mkdir()
    (root / "pbs-logs").mkdir()
    jobs = []
    for channel in CHANNELS:
        items = eligible[channel]
        for start in range(0, len(items), per_shard):
            subset = items[start:start+per_shard]
            i = len(jobs)
            rows = []
            receipt_shas = []
            indices = []
            for index, receipt, receipt_path, d in subset:
                indices.append(index)
                expected_original = source_path(selection, channel, index)
                rows.append([channel.lower(), index, expected_original,
                             *(d / INPUT_FILES[key] for key in INPUT_FILES)])
                receipt_shas.append(sha(receipt_path))
            path = root / "locked/shards" / f"shard_{i:04d}.tsv"
            make_manifest(path, rows)
            jobs.append({
                "id": i, "channel": channel, "indices": indices,
                "receipt_sha256": receipt_shas,
                "manifest_sha256": sha(path), "n_muons": 20*len(indices),
            })
    plan = {
        "schema": SCHEMA,
        "analysis_git": git(),
        "analysis_binary_sha256": sha(ANALYZER),
        "merge_binary_sha256": sha(MERGER),
        "created_utc": timestamp(),
        "source_root": str(SOURCE_ROOT),
        "source_selection_sha256": sel_sha,
        "source_binary_lock_sha256": bin_sha,
        "physics_git": ORIGINAL_GIT,
        "original_physics_binary_sha256": binary["binary_sha256"],
        "pitch_mm": 3.0, "voltage_V": 1800, "per_shard_chunks": per_shard,
        "n_chunks": n_chunks, "n_muons": n_chunks * 20,
        "n_shards": len(jobs), "shards": jobs,
        "channel_muons": {c: 20*len(v) for c,v in eligible.items()}
    }
    save_json(root / "locked/campaign.json", plan)
    stage("PREPARE_PASS original physics only; zero Garfield jobs "
          f"n_chunks={n_chunks} n_muons={20*n_chunks} "
          f"n_shards={len(jobs)} channel_muons={plan['channel_muons']}")
    print(f"FROZEN_CAMPAIGN={root}\nSHARDS=0-{len(jobs)-1}", flush=True)


def audit_files(folder: Path) -> None:
    control = folder / "SHA256SUMS"
    require(control.is_file(), f"missing SHA256SUMS in {folder}")
    n = 0
    for raw in control.read_text(encoding="utf-8").splitlines():
        digest, sep, name = raw.partition("  ")
        require(sep and len(digest)==64 and name and "/" not in name and
                name not in (".", "..") and not name.startswith("."),
                "unsafe/malformed checksum path " + raw)
        f = folder/name
        require(f.is_file() and sha(f)==digest,
                f"output file mutated or incomplete {f}")
        n += 1
    kind = MERGED_OUTPUT if (folder/"wiretruth1800_full.root").exists() else SHARD_OUTPUT
    require(n >= len(kind), "incomplete sealed ROOT output")
    for name in kind:
        require((folder/name).is_file(),
                f"expected {('merged' if kind is MERGED_OUTPUT else 'shard')} output absent: {name}")
    require((folder / "sealed.json").is_file(), "missing sealed provenance")


def write_checksums(folder: Path) -> None:
    out = folder / "SHA256SUMS"
    require(not out.exists(), "refusing overwrite of SHA256SUMS")
    files = sorted(p for p in folder.iterdir() if p.is_file())
    with out.open("x", encoding="utf-8") as f:
        for p in files:
            f.write(f"{sha(p)}  {p.name}\n")
    audit_files(folder)


def copytree_verify(source: Path, dest: Path) -> None:
    require(not dest.exists(), "output or staging path already exists")
    shutil.copytree(source, dest)
    audit_files(dest)


def run_shard(root: Path, index: int) -> None:
    plan, selection, binary = job_list(root)
    require(0 <= index < plan["n_shards"], "invalid array shard id")
    shard = plan["shards"][index]
    final = root / "shards" / f"shard_{index:04d}"
    if final.exists():
        audit_files(final)
        require(load(final/"sealed.json")["manifest_sha256"] ==
                shard["manifest_sha256"], "existing shard has different manifest")
        stage(f"SKIP already-verified completed shard={index}")
        return
    scratch = Path(os.environ.get("TMPDIR", ""))
    require(str(scratch) not in ("", ".") and scratch.is_dir() and
            os.access(scratch, os.W_OK),
            "PBS TMPDIR is required for local staged immutable input/hash/ROOT processing")
    require(os.environ.get("PBS_JOBID", ""),
            "per-shard full analysis must run on PBS compute node")
    manifest = root / "locked/shards" / f"shard_{index:04d}.tsv"
    rows = manifest_rows(manifest)
    require(len(rows) == len(shard["indices"]), "shard list length changed")
    claim = root / "attempts" / f"claim_{index:04d}"
    claim.mkdir(exist_ok=False)
    (claim/"owner_job_id").write_text(os.environ["PBS_JOBID"]+"\n",encoding="utf-8")
    # On abnormal termination keep claim for operator diagnosis; no mutation
    # to original physics products, no hidden duplicate PBS reruns.
    with tempfile.TemporaryDirectory(prefix=f"w1800_{index:04d}_",dir=scratch) as tmp:
        tmp = Path(tmp)
        safe_rows = []
        stage(f"SHARD_BEGIN={index} n_chunks={len(rows)} channel={shard['channel']} "
              f"scratch={tmp}")
        for j, row in enumerate(rows):
            channel = shard["channel"]
            chunkid = int(row["chunk"])
            require(channel.lower() == row["channel"] and
                    chunkid == shard["indices"][j],
                    "shard original source index differs from frozen manifest")
            receipt_item = completed_receipt(selection,binary,channel,chunkid)
            require(receipt_item is not None, "source verified receipt disappeared")
            receipt, p, _ = receipt_item
            require(sha(p) == shard["receipt_sha256"][j],
                    "original completed receipt changed after snapshot")
            local = tmp / "inputs" / f"{channel.lower()}_{chunkid:04d}"
            source = source_path(selection,channel,chunkid)
            require(str(source) == row["original_input"], "original path mismatch")
            target = local/"original_input.csv"
            copyhash(source,target,receipt["input_sha256"])
            staged = {"original_input":target}
            for key,filename in INPUT_FILES.items():
                origin = Path(row[key])
                require(origin == p.parent/filename,
                        "physical archive path changed after frozen manifest")
                t = local/filename
                expected = receipt["files"][key]
                copyhash(origin,t,expected["sha256"],int(expected["bytes"]))
                staged[key] = t
            safe_rows.append([channel.lower(),chunkid,
                              *(staged[key] for key in MANIFEST_HEADER[2:])])
            stage(f"SHARD_SOURCE_SHA_PASS={index} {channel}:{chunkid:04d}")
        checked = tmp/"checked_manifest.tsv"
        make_manifest(checked,safe_rows)
        stage(f"SHARD_PHYSICS_BEGIN={index}")
        output = tmp/"output"
        subprocess.run([str(ANALYZER),str(checked),str(output),"--no-plots"],check=True)
        stage(f"SHARD_PHYSICS_PASS={index}")
        for filename in SHARD_OUTPUT:
            require((output/filename).is_file(),
                    f"missing post-physics ROOT/CSV output {filename}")
        with (output/"channel_station_summary.csv").open(newline="") as f:
            counts = {r["scope"]:r for r in csv.DictReader(f)}
        require(int(counts["all"]["n_muons"])==shard["n_muons"],
                "C++ output independent Geant4 muon count differs from shard plan")
        shutil.copy2(manifest,output/"frozen_original_source_manifest.tsv")
        save_json(output/"sealed.json",{
            "schema":SCHEMA,"shard":index,"channel":shard["channel"],
            "indices":shard["indices"],"n_muons":shard["n_muons"],
            "manifest_sha256":shard["manifest_sha256"],
            "analysis_git":plan["analysis_git"],
            "analysis_binary_sha256":plan["analysis_binary_sha256"],
            "source_physics_binary_sha256":plan["original_physics_binary_sha256"],
            "pbs_job_id":os.environ["PBS_JOBID"],"sealed_utc":timestamp()
        })
        write_checksums(output)
        temp = root/"attempts"/f"staging_{index:04d}_{os.getpid()}"
        stage(f"SHARD_COPY_BEGIN={index}")
        copytree_verify(output,temp)
        require(not final.exists(), "another worker created same original shard")
        temp.rename(final)
        stage(f"SHARD_PERSIST_PASS={index} muons={shard['n_muons']}")
    (claim/"owner_job_id").unlink()
    claim.rmdir()


def status(root: Path) -> None:
    plan = load(root/"locked/campaign.json")
    successes=[];missing=[];bad=[]
    for shard in plan["shards"]:
        index=shard["id"]
        f=root/"shards"/f"shard_{index:04d}"
        if (f/"sealed.json").is_file():
            try:
                a=load(f/"sealed.json")
                require(a["manifest_sha256"]==shard["manifest_sha256"] and
                        a["n_muons"]==shard["n_muons"],
                        "completed job manifest differs from frozen plan")
                successes.append(index)
            except Exception as e:
                bad.append((index,str(e)))
        else:
            missing.append(index)
    done=sum(plan["shards"][i]["n_muons"] for i in successes)
    print(f"ORIGINAL_1800V_FROZEN muons={plan['n_muons']} "
          f"chunks={plan['n_chunks']} shards={plan['n_shards']} "
          f"analysis_git={plan['analysis_git']}")
    print(f"COMPLETE_SHARDS={len(successes)} COMPLETE_MUONS={done} "
          f"REMAINING_MUONS={plan['n_muons']-done} "
          f"MISSING_SHARDS={len(missing)} BAD_SHARDS={len(bad)}")
    print("MISSING_IDS="+",".join(map(str,missing)))
    if bad:
        print("INVALID_SHARDS="+str(bad))
    if (root/"merged/SHA256SUMS").exists():
        print("MERGED_OUTPUT_EXISTS="+str(root/"merged"))


def merge(root: Path) -> None:
    plan, _, _ = job_list(root)
    final = root/"merged"
    require(not final.exists(), "merged final result already exists; never overwrite")
    require(os.environ.get("PBS_JOBID",""),
            "combine all ROOT distributions only inside a PBS compute job")
    scratch=Path(os.environ.get("TMPDIR",""))
    require(str(scratch) not in ("",".") and scratch.is_dir() and
            os.access(scratch,os.W_OK),"PBS TMPDIR required for merger")
    with tempfile.TemporaryDirectory(prefix="w1800merge_",dir=scratch) as work:
        work=Path(work)
        listfile=work/"shard_dirs.txt"
        total=0
        with listfile.open("x",encoding="utf-8") as lf:
            for shard in plan["shards"]:
                n=shard["id"]
                original=root/"shards"/f"shard_{n:04d}"
                require(original.is_dir(),f"missing completed shard {n}")
                audit_files(original)
                meta=load(original/"sealed.json")
                require(meta["manifest_sha256"]==shard["manifest_sha256"] and
                        meta["analysis_binary_sha256"]==plan["analysis_binary_sha256"] and
                        meta["n_muons"]==shard["n_muons"],
                        f"provenance changed for original shard {n}")
                local=work/"inputs"/f"shard_{n:04d}"
                local.mkdir(parents=True)
                for filename in ("wiretruth1800_8k.root",
                                 "channel_station_summary.csv",
                                 "asymmetry_denominators.csv"):
                    expected=dict((s.split("  ",1)[1],s.split("  ",1)[0])
                                  for s in (original/"SHA256SUMS").read_text().splitlines())[filename]
                    copyhash(original/filename,local/filename,expected)
                lf.write(str(local)+"\n")
                total+=shard["n_muons"]
                stage(f"MERGE_SOURCE_SHA_PASS={n} total_selected_muons={total}")
        require(total==plan["n_muons"], "muon total differs from campaign plan")
        output=work/"output"
        output.mkdir()
        stage(f"MERGE_ROOT_BEGIN shards={plan['n_shards']} expected_muons={total}")
        subprocess.run([str(MERGER),str(listfile),str(output),str(total)],check=True)
        stage("MERGE_ROOT_PASS preparing provenance and final SHA")
        shutil.copy2(root/"locked/campaign.json",output/"frozen_campaign.json")
        save_json(output/"sealed.json",{
            "schema":SCHEMA, "n_muons":total,"n_chunks":plan["n_chunks"],
            "n_shards":plan["n_shards"],"analysis_git":plan["analysis_git"],
            "analysis_binary_sha256":plan["analysis_binary_sha256"],
            "merge_binary_sha256":plan["merge_binary_sha256"],
            "pbs_job_id":os.environ["PBS_JOBID"],
            "sealed_utc":timestamp()
        })
        write_checksums(output)
        temp=root/"attempts"/f"merged_staging_{os.getpid()}"
        stage("MERGE_COPY_BEGIN")
        copytree_verify(output,temp)
        require(not final.exists(),"merger result appeared simultaneously")
        temp.rename(final)
        stage(f"MERGE_PASS final_muons={total} combined_output={final}")


def main() -> None:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("mode",choices=("prepare","run-shard","status","merge"))
    p.add_argument("--campaign",type=Path,required=True)
    p.add_argument("--per-shard",type=int,default=6)
    p.add_argument("--shard",type=int)
    a=p.parse_args()
    root=a.campaign.resolve()
    if a.mode=="prepare":prepare(root,a.per_shard)
    elif a.mode=="run-shard":
        require(a.shard is not None,"run-shard requires --shard")
        run_shard(root,a.shard)
    elif a.mode=="status":status(root)
    else:merge(root)


if __name__=="__main__":
    try:
        main()
    except Exception as e:
        print(f"FAILED 1800V FULL: {type(e).__name__}: {e}",file=sys.stderr,flush=True)
        sys.exit(2)
