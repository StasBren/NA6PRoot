#!/usr/bin/env python3
"""One-file, exact, read-only 1800 V 5620-muon Heed+Garfield ROOT export.

Run on PBS. Adds Heed cluster / primary-start trees to a SHA-checked COPY of
the previously merged 5620-muon ROOT. No physics reruns or source mutation.
"""
from __future__ import annotations
import argparse
import csv
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile

from full_wiretruth1800_campaign import (
    require, sha, load, save_json, timestamp, copyhash, manifest_rows,
    original_locks, completed_receipt, source_path, INPUT_FILES
)
from merge_1800v_verified_subset import (
    frozen_plan, FROZEN, FILES, file_sha_table, SUBSET_SCHEMA
)

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
BINARY = REPO / "test/mwpc-stage-c/build-wiretruth1800-heed-enrichment/enrich_heed_root"
SCHEMA = "NA60_1800V_5620_HEED_ENRICHED_V1"
OUTPUT = "wiretruth1800_enriched_5620.root"
NEEDED = ("heed_clusters", "heed_primaries", "seed", "crossing")
CLUSTER_COLUMNS = ("crossing_id", "cluster_index", "u_rel_cm",
                   "v_abs_cm", "w_rel_cm", "t_ns", "n_original", "n_accepted")
PRIMARY_COLUMNS = (
    "crossing_id", "primary_index", "cluster_index",
    "u_rel_cm", "v_abs_cm", "w_rel_cm", "t_ns",
    "raw_heed_energy_eV", "raw_dx", "raw_dy", "raw_dz",
    "seed_u_mm", "seed_v_mm", "seed_w_mm", "seed_energy_eV",
    "garfield_seed", "garfield_ne", "garfield_ni", "ion_births",
    "collected_e", "attached", "other", "successful", "capped"
)


def read_csv(path: Path, required: tuple[str, ...]):
    with path.open(encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f)
        require(reader.fieldnames is not None
                and len(reader.fieldnames) == len(set(reader.fieldnames))
                and set(required).issubset(reader.fieldnames),
                "invalid input CSV header: " + str(path))
        yield from reader


def checksums(directory: Path) -> dict[str, str]:
    control = directory / "SHA256SUMS"
    require(control.is_file(), "missing SHA256SUMS " + str(directory))
    result = {}
    for line in control.read_text(encoding="utf-8").splitlines():
        digest, sep, name = line.partition("  ")
        require(sep and len(digest) == 64 and name not in result
                and name and "/" not in name and not name.startswith("."),
                "malformed checksum list")
        result[name] = digest
    return result


def check_files(directory: Path):
    result = checksums(directory)
    require(OUTPUT in result and "sealed.json" in result,
            "incomplete enriched output")
    for name, want in result.items():
        require((directory/name).is_file() and sha(directory/name) == want,
                "output checksum mismatch: " + name)


def write_files(directory: Path):
    target = directory/"SHA256SUMS"
    require(not target.exists(), "refusing duplicate SHA256SUMS")
    with target.open("x", encoding="utf-8") as f:
        for p in sorted(directory.iterdir()):
            if p.is_file():
                f.write(f"{sha(p)}  {p.name}\n")
    check_files(directory)


def plan_and_subset(subset_root: Path):
    original = frozen_plan()
    root = subset_root/"locked/subset.json"
    part = load(root)
    require(part["schema"] == SUBSET_SCHEMA
            and part["n_muons"] == 5620 and part["n_shards"] == 48,
            "unexpected dashboard subset schema/size")
    require(part["source_campaign_sha256"] == sha(FROZEN/"locked/campaign.json"),
            "source campaign snapshot changed")
    require(part["source_campaign"] == str(FROZEN),
            "unrecognized frozen campaign")
    require(tuple(part["excluded_missing_ids"]) == (0,5,7,16,17,27,32),
            "missing shard selection changed")
    require(len({x["id"] for x in part["shards"]}) == 48,
            "duplicate subset shard")
    merged = subset_root/"merged"
    source_hashes = checksums(merged)
    root_file = merged/"wiretruth1800_full.root"
    require(source_hashes.get("wiretruth1800_full.root") == sha(root_file),
            "existing dashboard ROOT file changed")
    return original, part, root_file, source_hashes["wiretruth1800_full.root"]


def prepare(out: Path, subset_root: Path):
    out = out.resolve()
    subset_root = subset_root.resolve()
    require(str(out).startswith(str((REPO/"test_runs").resolve()) + "/"),
            "enriched export must be a new dashboard worktree test_runs path")
    require(not out.exists(), "output exists; never overwrite")
    require(BINARY.is_file(), "build ROOT importer before preparation")
    source, part, root_file, root_sha = plan_and_subset(subset_root)
    manifest = {
        "schema": SCHEMA, "created_utc": timestamp(),
        "source_subset": str(subset_root),
        "source_subset_sha256": sha(subset_root/"locked/subset.json"),
        "source_root": str(root_file),
        "source_root_sha256": root_sha,
        "source_campaign_sha256": sha(FROZEN/"locked/campaign.json"),
        "frozen_cxx_merger_sha256": source["merge_binary_sha256"],
        "new_converter": str(BINARY), "new_converter_sha256": sha(BINARY),
        "n_muons": 5620, "n_shards": 48,
        "missing_ids": part["excluded_missing_ids"],
        "trees_preserved": ["muons","wires","seeds","seed_wires"],
        "trees_added": ["heed_clusters","heed_primaries"],
        "physical_caveat": "No full drift paths or individual electron endpoints are archived",
    }
    (out/"locked").mkdir(parents=True)
    (out/"pbs-logs").mkdir()
    (out/"attempts").mkdir()
    save_json(out/"locked/export.json", manifest)
    print(f"HEED_PREPARE_PASS 5620 muons 48 sealed shards input_root_sha={root_sha}")
    print("HEED_EXPORT_DIRECTORY=" + str(out))


def verify_source_files(path: Path, receipt, names=NEEDED):
    source = {}
    for key in names:
        data = receipt["files"][key]
        p = path/f"microscopic_{key}.csv"
        require(p.is_file() and p.stat().st_size == int(data["bytes"])
                and sha(p) == data["sha256"],
                "original Garfield CSV mutated: " + str(p))
        source[key] = p
    return source


def append_chunk(source: dict[str, Path], cf, pf):
    crossing_rows = list(read_csv(source["crossing"], (
        "source|chunk|event|hit_index|track_id|station|det_id",
        "heed_clusters", "heed_primaries")))
    expect = {}
    for x in crossing_rows:
        id_ = x["source|chunk|event|hit_index|track_id|station|det_id"]
        require(id_ and id_ not in expect, "duplicate crossing in original chunk")
        expect[id_] = (int(x["heed_clusters"]),int(x["heed_primaries"]))
    require(len(expect) == 20, "expected 20 independent Geant4 muons per chunk")
    clusters: dict[tuple[str,int], int] = {}
    counts_cl = {id_: 0 for id_ in expect}
    accepted_cl = {id_: 0 for id_ in expect}
    nclusters = 0
    for x in read_csv(source["heed_clusters"], CLUSTER_COLUMNS):
        id_=x["crossing_id"]; j=int(x["cluster_index"])
        orig=int(x["n_original"]); accepted=int(x["n_accepted"])
        require(id_ in expect and j >= 0 and orig >= accepted >= 0
                and (id_,j) not in clusters, "invalid original Heed cluster")
        clusters[(id_,j)] = accepted
        counts_cl[id_] += 1
        accepted_cl[id_] += accepted
        cf.writerow([x[k] for k in CLUSTER_COLUMNS])
        nclusters += 1
    for id_, (nc,np) in expect.items():
        require(counts_cl[id_] == nc and accepted_cl[id_] == np
                and sorted(j for (c,j) in clusters if c == id_) == list(range(nc)),
                "Heed cluster count/index mismatch")
    seeds = {}
    seed_keys = ("crossing_id","primary_index","seed_u_mm","seed_v_mm",
                 "seed_w_mm","seed_energy_eV","garfield_seed","garfield_ne",
                 "garfield_ni","ion_births","collected_e","attached",
                 "other","successful","capped")
    for x in read_csv(source["seed"], seed_keys):
        key=(x["crossing_id"],int(x["primary_index"]))
        require(key[0] in expect and key not in seeds and key[1]>=0
                and int(x["capped"])==0, "invalid/duplicate native seed")
        seeds[key]=x
    prim_counts = {id_:0 for id_ in expect}
    cluster_prim_counts: dict[tuple[str,int],int] = {}
    used = set()
    nprim=0
    for x in read_csv(source["heed_primaries"], (
        "crossing_id","primary_index","cluster_index","u_rel_cm",
        "v_abs_cm","w_rel_cm","t_ns","raw_heed_energy_eV",
        "raw_dx","raw_dy","raw_dz")):
        key=(x["crossing_id"],int(x["primary_index"]))
        ckey=(key[0],int(x["cluster_index"]))
        require(key in seeds and key not in used and ckey in clusters,
                "primary missing Garfield seed or parent Heed cluster")
        used.add(key);prim_counts[key[0]] += 1
        cluster_prim_counts[ckey] = cluster_prim_counts.get(ckey,0)+1
        s=seeds[key]
        pf.writerow([x[k] for k in ("crossing_id","primary_index",
                                    "cluster_index","u_rel_cm","v_abs_cm",
                                    "w_rel_cm","t_ns","raw_heed_energy_eV",
                                    "raw_dx","raw_dy","raw_dz")] +
                   [s[k] for k in seed_keys[2:]])
        nprim+=1
    require(len(used)==len(seeds), "Heed/Garfield primary population mismatch")
    for id_, (_,np) in expect.items():
        require(prim_counts[id_]==np
                and sorted(i for (c,i) in used if c==id_) == list(range(np)),
                "wrong primary indices/count per Geant4 crossing")
    for k,n in clusters.items():
        require(cluster_prim_counts.get(k,0)==n,
                "Heed electron children differ from stored n_accepted")
    return len(expect),nclusters,nprim


def run(out: Path):
    require(os.environ.get("PBS_JOBID"), "submit via PBS compute node")
    scratch=Path(os.environ.get("TMPDIR",""))
    require(str(scratch) not in ("",".") and scratch.is_dir()
            and os.access(scratch,os.W_OK), "PBS TMPDIR missing")
    meta=load(out/"locked/export.json")
    require(meta["schema"]==SCHEMA and meta["n_muons"]==5620,
            "wrong export metadata")
    subset_root=Path(meta["source_subset"])
    original, subset, root_file, root_sha=plan_and_subset(subset_root)
    require(meta["source_subset_sha256"]==sha(subset_root/"locked/subset.json")
            and meta["source_root_sha256"]==root_sha
            and meta["source_campaign_sha256"]==sha(FROZEN/"locked/campaign.json")
            and meta["new_converter_sha256"]==sha(BINARY),
            "ROOT conversion provenance drift")
    require(not (out/"enriched").exists(), "output already exists")
    selection,binary,selsha,binsha=original_locks()
    require(selsha==original["source_selection_sha256"]
            and binsha==original["source_binary_lock_sha256"],
            "original 1800 V physics locks changed")
    with tempfile.TemporaryDirectory(prefix="w1800_heed_",dir=scratch) as td:
        tmp=Path(td)
        stages=tmp/"import"
        stages.mkdir()
        clusters_tsv=stages/"clusters.tsv"
        primaries_tsv=stages/"primaries.tsv"
        tally={"muons":0,"clusters":0,"primaries":0,"chunks":0}
        with clusters_tsv.open("x",encoding="utf-8",newline="") as fc, \
             primaries_tsv.open("x",encoding="utf-8",newline="") as fp:
            cw=csv.writer(fc,delimiter="\t",lineterminator="\n")
            pw=csv.writer(fp,delimiter="\t",lineterminator="\n")
            cw.writerow(CLUSTER_COLUMNS)
            pw.writerow(PRIMARY_COLUMNS)
            for item in subset["shards"]:
                idx=item["id"]
                locked=original["shards"][idx]
                require(locked["id"]==idx and
                        locked["manifest_sha256"]==item["manifest_sha256"] and
                        locked["n_muons"]==item["n_muons"],
                        "subset shard differs from original lock")
                folder=FROZEN/"locked/shards"/f"shard_{idx:04d}.tsv"
                require(sha(folder)==locked["manifest_sha256"],
                        "original locked shard manifest changed")
                rows=manifest_rows(folder)
                require(len(rows)==len(locked["indices"]),
                        "original shard chunk count changed")
                for row,chunkid in zip(rows,locked["indices"]):
                    channel=locked["channel"]
                    require(row["channel"]==channel.lower()
                            and int(row["chunk"])==chunkid
                            and row["original_input"]==str(source_path(selection,channel,chunkid)),
                            "chunk identity changed")
                    got=completed_receipt(selection,binary,channel,chunkid)
                    require(got is not None, "completed original Garfield receipt missing")
                    receipt,rpath,source_dir=got
                    position=locked["indices"].index(chunkid)
                    require(sha(rpath)==locked["receipt_sha256"][position],
                            "frozen receipt hash changed")
                    checked=verify_source_files(source_dir,receipt)
                    nm,nc,np=append_chunk(checked,cw,pw)
                    tally["muons"]+=nm
                    tally["clusters"]+=nc
                    tally["primaries"]+=np
                    tally["chunks"]+=1
                print(f"HEED_SOURCE_SHA_PASS shard={idx} "
                      f"total_muons={tally['muons']} total_primaries={tally['primaries']}",
                      flush=True)
        require(tally["muons"]==5620 and tally["chunks"]==281,
                "subset expansion does not reproduce original 5620 muons")
        rootcopy=tmp/OUTPUT
        copyhash(root_file,rootcopy,root_sha)
        subprocess.run([str(BINARY),str(rootcopy),str(clusters_tsv),
                        str(primaries_tsv),"5620"],check=True)
        output=tmp/"completed"
        output.mkdir()
        rootcopy.rename(output/OUTPUT)
        (output/"README_DATA.txt").write_text(
            "NA60+ 1800 V independent archived 5620 Geant4 muons.\n"
            "Trees: muons, wires, seeds, seed_wires, heed_clusters, heed_primaries.\n"
            "Join heed_primaries -> heed_clusters on crossing_id,cluster_index.\n"
            "Join heed_primaries -> seed_wires on crossing_id,primary_index.\n"
            "t_ns is Heed initial electron time; mean_t_collected_ns from "
            "seed_wires includes avalanche development, NOT pure drift time.\n"
            "Heed u_rel_cm/w_rel_cm are local; primary seed_u_mm refers to the "
            "Garfield periodic wire patch. Do not identify it with global wire u.\n"
            "raw_heed_energy_eV/raw_d* are NOT meaningful initial kinematics "
            "when TrackHeed delta electron transport is enabled.\n"
            "No full electron trajectories, individual endpoint coordinates or "
            "full ion-birth coordinate census was archived. Original ion birth "
            "reservoir CSV remains on WIPP, NOT in this ROOT.\n"
            "Excluded original 7 failed analysis shards. Do not relabel as 6460.\n",
            encoding="utf-8"
        )
        save_json(output/"sealed.json",{
            "schema":SCHEMA, "n_muons":5620,
            "tally":tally, "source_root_sha256":root_sha,
            "source_subset_sha256":meta["source_subset_sha256"],
            "converter_sha256":meta["new_converter_sha256"],
            "pbs_job_id":os.environ["PBS_JOBID"],"sealed_utc":timestamp()
        })
        write_files(output)
        target=out/"attempts"/f"enriched_staging_{os.getpid()}"
        require(not target.exists(), "duplicate staging destination")
        shutil.copytree(output,target)
        check_files(target)
        require(not (out/"enriched").exists(), "another export appeared")
        target.rename(out/"enriched")
        print(f"HEED_EXPORT_PASS {tally} FILE={out/'enriched'/OUTPUT}",
              flush=True)


def status(out: Path):
    meta=load(out/"locked/export.json")
    print(f"HEED_EXPORT_SOURCE=5620 muons, source={meta['source_root']}")
    target=out/"enriched"
    if (target/"SHA256SUMS").exists():
        check_files(target)
        print(f"HEED_EXPORT_READY={target/OUTPUT}")
        print(f"ROOT_SIZE_BYTES={(target/OUTPUT).stat().st_size}")
    else:
        print("HEED_EXPORT_PENDING; inspect dedicated PBS job logs")


def main():
    p=argparse.ArgumentParser()
    p.add_argument("mode",choices=("prepare","run","status"))
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--subset",type=Path)
    a=p.parse_args()
    out=a.output.resolve()
    if a.mode=="prepare":
        require(a.subset is not None,"--subset required for preparation")
        prepare(out,a.subset)
    elif a.mode=="run":run(out)
    else:status(out)


if __name__=="__main__":
    try:main()
    except Exception as exc:
        raise SystemExit("HEED_EXPORT_FAIL: "+str(exc))
