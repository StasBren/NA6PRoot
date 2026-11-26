#!/usr/bin/env python3
"""Resumable 48-shard Heed+Garfield enrichment of the verified 5620-muon ROOT.

No Garfield++ reruns. Each PBS worker verifies original physics receipts,
streams four original CSVs once to its local scratch with SHA-256, materializes
two normalized tables, and atomically seals an immutable shard export.

A separate PBS reducer checks all 48 seals, joins their TSVs and calls the
unchanged ROOT C++ enrichment executable against a copy of verified ROOT.
"""
from __future__ import annotations

import argparse
import csv
import os
from pathlib import Path
import shutil
import subprocess
import tempfile

from full_wiretruth1800_campaign import (
    require, sha, load, save_json, timestamp, copyhash,
    manifest_rows, original_locks, completed_receipt, source_path,
)
from merge_1800v_verified_subset import FROZEN, frozen_plan
from export_1800v_heed_root import (
    BINARY, OUTPUT, CLUSTER_COLUMNS, PRIMARY_COLUMNS, SCHEMA,
    append_chunk, checksums, plan_and_subset,
)

STAGE_SCHEMA = "NA60_1800V_HEED_SHARDED_5620_STAGE_V1"
OUTFILES = ("clusters.tsv", "primaries.tsv")
N_EXPECTED = 48
CHUNKS_EXPECTED = 281
MUONS_EXPECTED = 5620


def stage_dir(out: Path, sid: int):
    return out/"shard_exports"/f"shard_{sid:04d}"


def sha_manifest(d: Path) -> dict[str,str]:
    hashes = checksums(d)
    require(set(hashes) == {"clusters.tsv","primaries.tsv","sealed.json"},
            f"bad sealed shard checksum set: {d}")
    for name, wanted in hashes.items():
        p=d/name
        require(p.is_file() and sha(p)==wanted,
                f"shard output failed SHA-256 validation: {p}")
    return hashes


def write_shard_checksums(directory: Path):
    """Checksum exactly three immutable payloads (never SHA256SUMS itself)."""
    expected=("clusters.tsv","primaries.tsv","sealed.json")
    require({p.name for p in directory.iterdir()}==set(expected),
            "unexpected file in sealed shard staging area")
    path=directory/"SHA256SUMS"
    require(not path.exists(),"refusing existing SHA256SUMS")
    with path.open("x",encoding="utf-8") as f:
        for name in sorted(expected):
            f.write(f"{sha(directory/name)}  {name}\n")
    return sha_manifest(directory)


def selftest():
    """Fast local sealing QA without touching original physics outputs."""
    with tempfile.TemporaryDirectory(prefix="heed_shard_selftest_") as td:
        d=Path(td)
        for name,data in (
            ("clusters.tsv","crossing_id\tcluster_index\n"),
            ("primaries.tsv","crossing_id\tprimary_index\n"),
            ("sealed.json",'{"schema":"selftest"}\n')):
            (d/name).write_text(data,encoding="utf-8")
        hashes=write_shard_checksums(d)
        require(set(hashes)=={"clusters.tsv","primaries.tsv","sealed.json"},
                "selftest returned wrong checksum set")
        with (d/"primaries.tsv").open("a",encoding="utf-8") as f:
            f.write("tampered\n")
        tamper_rejected=False
        try:
            sha_manifest(d)
        except Exception:
            tamper_rejected=True
        require(tamper_rejected,"tampered shard was accepted")
    print("HEED_SHARD_SELFTEST_PASS checksum manifest and tamper rejection")


def source_state(out: Path):
    meta=load(out/"locked/export.json")
    require(meta["schema"]==SCHEMA and meta["n_shards"]==N_EXPECTED and
            meta["n_muons"]==MUONS_EXPECTED and
            meta["new_converter_sha256"]==sha(BINARY),
            "export snapshot/converter mismatch")
    subsetroot=Path(meta["source_subset"])
    original,subset,source_root,root_hash=plan_and_subset(subsetroot)
    require(sha(subsetroot/"locked/subset.json")==meta["source_subset_sha256"]
            and sha(FROZEN/"locked/campaign.json")==meta["source_campaign_sha256"]
            and root_hash==meta["source_root_sha256"],
            "source ROOT/selection drift since export snapshot")
    require(len(subset["shards"])==N_EXPECTED,
            "unexpected number of verified shards")
    return meta,original,subset,source_root,root_hash


def stage(out: Path, array_index: int):
    require(os.environ.get("PBS_JOBID"),"stage must run under PBS")
    scratch=Path(os.environ.get("TMPDIR",""))
    require(str(scratch) not in ("",".") and scratch.is_dir() and
            os.access(scratch,os.W_OK), "PBS scratch missing")
    meta,original,subset,source_root,root_sha=source_state(out)
    require(0<=array_index<len(subset["shards"]),"bad PBS array index")
    item=subset["shards"][array_index]
    sid=item["id"]
    target=stage_dir(out,sid)
    target.parent.mkdir(parents=True,exist_ok=True)
    if target.is_dir():
        sha_manifest(target)
        seal=load(target/"sealed.json")
        require(seal["schema"]==STAGE_SCHEMA and
                seal["shard"]==sid and
                seal["manifest_sha256"]==item["manifest_sha256"],
                "existing shard doesn't match frozen source selection")
        print(f"HEED_SHARD_SKIP_ALREADY_VERIFIED={sid}",flush=True)
        return
    print(f"HEED_STAGE_ENTER utc={timestamp()} array={array_index} "
          f"shard={sid} expected_muons={item['n_muons']}",flush=True)
    print(f"HEED_LOCKS_START utc={timestamp()} shard={sid}",flush=True)
    selection,binary,selsha,binsha=original_locks()
    require(selsha==original["source_selection_sha256"] and
            binsha==original["source_binary_lock_sha256"],
            "source physics locks changed")
    print(f"HEED_LOCKS_PASS utc={timestamp()} shard={sid}",flush=True)
    locked=original["shards"][sid]
    require(locked["id"]==sid and
            locked["manifest_sha256"]==item["manifest_sha256"] and
            locked["n_muons"]==item["n_muons"],"shard freeze mismatch")
    src_manifest=FROZEN/"locked/shards"/f"shard_{sid:04d}.tsv"
    require(sha(src_manifest)==item["manifest_sha256"],"original TSV changed")
    rows=manifest_rows(src_manifest)
    require(len(rows)==len(locked["indices"])==
            len(locked["receipt_sha256"]),"chunk/receipt assignment changed")
    with tempfile.TemporaryDirectory(prefix=f"heed_{sid:04d}_",dir=scratch) as td:
        tmp=Path(td)
        output=tmp/"output"
        output.mkdir()
        tally={"muons":0,"clusters":0,"primaries":0,"chunks":0}
        with (output/"clusters.tsv").open("x",encoding="utf-8",newline="") as fc,\
             (output/"primaries.tsv").open("x",encoding="utf-8",newline="") as fp:
            cw=csv.writer(fc,delimiter="\t",lineterminator="\n")
            pw=csv.writer(fp,delimiter="\t",lineterminator="\n")
            cw.writerow(CLUSTER_COLUMNS);pw.writerow(PRIMARY_COLUMNS)
            for row,chunkid,receiptsha in zip(
                    rows,locked["indices"],locked["receipt_sha256"]):
                channel=locked["channel"]
                require(row["channel"]==channel.lower() and
                        int(row["chunk"])==chunkid and
                        row["original_input"]==str(source_path(selection,channel,chunkid)),
                        "chunk source mismatch")
                print(f"HEED_CHUNK_START utc={timestamp()} shard={sid} "
                      f"channel={channel} chunk={chunkid}",flush=True)
                got=completed_receipt(selection,binary,channel,chunkid)
                require(got is not None,"missing verified physics receipt")
                receipt,receipt_path,source_dir=got
                require(sha(receipt_path)==receiptsha,
                        "receipt SHA differs from frozen campaign")
                staged={}
                one=tmp/"input"/f"chunk_{chunkid:04d}"
                one.mkdir(parents=True,exist_ok=True)
                for k in ("crossing","seed","heed_clusters","heed_primaries"):
                    orig=source_dir/f"microscopic_{k}.csv"
                    r=receipt["files"][k]
                    require(orig.is_file() and orig.stat().st_size==int(r["bytes"]),
                            "missing or wrong-sized source CSV: "+str(orig))
                    local=one/orig.name
                    print(f"HEED_SHA_COPY_START utc={timestamp()} shard={sid} "
                          f"chunk={chunkid} key={k} bytes={r['bytes']}",flush=True)
                    copyhash(orig,local,r["sha256"],int(r["bytes"]))
                    staged[k]=local
                    print(f"HEED_SHA_COPY_PASS utc={timestamp()} shard={sid} "
                          f"chunk={chunkid} key={k}",flush=True)
                nm,nc,np=append_chunk(staged,cw,pw)
                for k,value in (("muons",nm),("clusters",nc),
                                ("primaries",np),("chunks",1)):
                    tally[k]+=value
                print(f"HEED_CHUNK_PASS utc={timestamp()} shard={sid} "
                      f"chunk={chunkid} muons={tally['muons']} "
                      f"primaries={tally['primaries']}",flush=True)
                shutil.rmtree(one)
        require(tally["muons"]==item["n_muons"] and
                tally["chunks"]==len(rows),"local shard counts do not conserve muons")
        save_json(output/"sealed.json",{
            "schema":STAGE_SCHEMA,"shard":sid,"array_index":array_index,
            "n_muons":tally["muons"],"tally":tally,
            "manifest_sha256":item["manifest_sha256"],
            "source_subset_sha256":meta["source_subset_sha256"],
            "physics_selection_sha256":selsha,
            "physics_binary_lock_sha256":binsha,
            "pbs_job_id":os.environ["PBS_JOBID"],
            "sealed_utc":timestamp(),
        })
        write_shard_checksums(output)
        tmpdest=out/"attempts"/f"staging_shard_{sid:04d}_{os.getpid()}"
        require(not tmpdest.exists() and not target.exists(),
                "another shard export is already present")
        shutil.copytree(output,tmpdest)
        sha_manifest(tmpdest)
        tmpdest.rename(target)
        print(f"HEED_SHARD_PASS shard={sid} "
              f"muons={tally['muons']} clusters={tally['clusters']} "
              f"primaries={tally['primaries']}",flush=True)


def inventory(out: Path, strict=False):
    _,_,subset,_,_=source_state(out)
    passed=[];bad=[];missing=[]
    for item in subset["shards"]:
        sid=item["id"]
        d=stage_dir(out,sid)
        if not d.exists():
            missing.append(sid)
            continue
        try:
            sha_manifest(d)
            seal=load(d/"sealed.json")
            meta=load(out/"locked/export.json")
            require(seal["schema"]==STAGE_SCHEMA and
                    seal["shard"]==sid and
                    seal["manifest_sha256"]==item["manifest_sha256"] and
                    seal["source_subset_sha256"]==meta["source_subset_sha256"],
                    "invalid shard seal")
            # Check exact manifest correspondence and population; avoid
            # relying on array order when grouping per channel.
            require(seal["manifest_sha256"]==item["manifest_sha256"] and
                    seal["n_muons"]==item["n_muons"],
                    "bad shard selection/count")
            passed.append((item,d,seal))
        except Exception as exc:
            bad.append((sid,str(exc)))
    print(f"HEED_STAGE_STATUS complete={len(passed)}/{N_EXPECTED} "
          f"missing={missing} invalid={bad}",flush=True)
    if strict:
        require(not missing and not bad and len(passed)==N_EXPECTED,
                "cannot merge incomplete/invalid shard exports")
    return passed


def combine(out: Path):
    require(os.environ.get("PBS_JOBID"),"combine must run on PBS")
    scratch=Path(os.environ.get("TMPDIR",""))
    require(str(scratch) not in ("",".") and scratch.is_dir(),
            "PBS scratch required")
    meta,original,subset,root_file,root_sha=source_state(out)
    require(not (out/"enriched").exists(),"final output already exists")
    selected=inventory(out,strict=True)
    totals={k:sum(item[2]["tally"][k] for item in selected)
            for k in ("muons","clusters","primaries","chunks")}
    require(totals["muons"]==MUONS_EXPECTED and
            totals["chunks"]==CHUNKS_EXPECTED,
            "sealed shard totals do not reproduce expected original crossings")
    with tempfile.TemporaryDirectory(prefix="heed_reduce_",dir=scratch) as td:
        tmp=Path(td)
        clusters=tmp/"clusters.tsv"
        primaries=tmp/"primaries.tsv"
        for filename,header,dest in (
                ("clusters.tsv",CLUSTER_COLUMNS,clusters),
                ("primaries.tsv",PRIMARY_COLUMNS,primaries)):
            with dest.open("xb") as result:
                result.write(("\t".join(header)+"\n").encode())
                for item,folder,seal in selected:
                    p=folder/filename
                    expected=checksums(folder)[filename]
                    local=tmp/f"stage_{item['id']:04d}_{filename}"
                    print(f"HEED_REDUCE_COPY_START utc={timestamp()} "
                          f"shard={item['id']} file={filename}",flush=True)
                    copyhash(p,local,expected)
                    with local.open("rb") as source:
                        line=source.readline().decode().strip()
                        require(line=="\t".join(header),
                                "wrong normalized TSV header")
                        shutil.copyfileobj(source,result,1024*1024)
                    local.unlink()
                    print(f"HEED_REDUCE_COPY_PASS utc={timestamp()} "
                          f"shard={item['id']} file={filename}",flush=True)
        target=tmp/OUTPUT
        print(f"HEED_REDUCE_ROOT_COPY_START utc={timestamp()}",flush=True)
        copyhash(root_file,target,root_sha)
        print(f"HEED_REDUCE_ROOT_COPY_PASS utc={timestamp()}",flush=True)
        print(f"HEED_REDUCE_ROOT_ENRICH_START utc={timestamp()}",flush=True)
        subprocess.run([str(BINARY),str(target),
                        str(clusters),str(primaries),"5620"],check=True)
        print(f"HEED_REDUCE_ROOT_ENRICH_PASS utc={timestamp()}",flush=True)
        result=tmp/"result"
        result.mkdir()
        target.rename(result/OUTPUT)
        (result/"README_DATA.txt").write_text(
            "NA60+ 1800 V, 5620 original independent Geant4 muons.\n"
            "ROOT trees: muons, wires, seeds, seed_wires, heed_clusters, heed_primaries.\n"
            "join cluster<->primary on crossing_id + cluster_index; "
            "primary<->seed_wires on crossing_id + primary_index.\n"
            "Units: Heed local u_rel_cm,w_rel_cm,v_abs_cm in cm; "
            "seed_u_mm,seed_v_mm,seed_w_mm in mm.\n"
            "mean_t_collected_ns includes avalanche development, not pure drift.\n"
            "No complete electron drift trajectories or individual endpoint positions.\n"
            "No complete ion-birth coordinate census in this ROOT. Original "
            "weighted ion birth reservoir archived separately on WIPP.\n"
            "Exactly 48 sealed successful shards; 7 excluded original analysis shards.\n",
            encoding="utf-8")
        save_json(result/"sealed.json",{
            "schema":SCHEMA,"n_muons":MUONS_EXPECTED,"tally":totals,
            "source_root_sha256":root_sha,
            "source_subset_sha256":meta["source_subset_sha256"],
            "converter_sha256":meta["new_converter_sha256"],
            "n_source_shards":N_EXPECTED,"pbs_job_id":os.environ["PBS_JOBID"],
            "sealed_utc":timestamp()
        })
        final_names=(OUTPUT,"README_DATA.txt","sealed.json")
        require({p.name for p in result.iterdir()}==set(final_names),
                "unexpected file in final ROOT result before sealing")
        with (result/"SHA256SUMS").open("x",encoding="utf-8") as f:
            for name in sorted(final_names):
                f.write(f"{sha(result/name)}  {name}\n")
        finalstage=out/"attempts"/f"final_staging_{os.getpid()}"
        require(not finalstage.exists(),"duplicate result staging")
        shutil.copytree(result,finalstage)
        finalhash=checksums(finalstage)
        for name,wanted in finalhash.items():
            require(sha(finalstage/name)==wanted,"final output SHA mismatch")
        finalstage.rename(out/"enriched")
        print(f"HEED_EXPORT_PASS file={out/'enriched'/OUTPUT} "
              f"size={(out/'enriched'/OUTPUT).stat().st_size} "
              f"total_primaries={totals['primaries']}",flush=True)


def main():
    p=argparse.ArgumentParser()
    p.add_argument("mode",choices=("selftest","stage","gate-one","status","gate","combine"))
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--array-index",type=int)
    a=p.parse_args()
    out=a.output.resolve()
    if a.mode=="selftest":
        selftest()
    elif a.mode=="stage":
        require(a.array_index is not None,"--array-index required")
        stage(out,a.array_index)
    elif a.mode=="gate-one":
        require(a.array_index is not None, "pilot index is required")
        _,_,subset,_,_=source_state(out)
        item=subset["shards"][a.array_index]
        folder=stage_dir(out,item["id"])
        require(folder.is_dir(),"pilot output does not exist")
        sha_manifest(folder)
        seal=load(folder/"sealed.json")
        require(seal["schema"]==STAGE_SCHEMA and
                seal["manifest_sha256"]==item["manifest_sha256"] and
                seal["n_muons"]==item["n_muons"],
                "pilot output source identity/count mismatch")
        print(f"HEED_PILOT_GATE_PASS array={a.array_index} "
              f"shard={item['id']} muons={seal['n_muons']}")
    elif a.mode=="status":
        inventory(out)
        if (out/"enriched"/OUTPUT).exists():
            print("HEED_ENRICHED_READY="+str(out/"enriched"/OUTPUT))
    elif a.mode=="gate":
        inventory(out,strict=True)
        print("HEED_SHARD_GATE_PASS all 48 sealed source inputs")
    else:combine(out)


if __name__=="__main__":
    try:main()
    except Exception as exc:raise SystemExit("HEED_SHARDED_FAIL: "+str(exc))
