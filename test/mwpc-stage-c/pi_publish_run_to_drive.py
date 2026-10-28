#!/usr/bin/env python3
"""Safely publish immutable summaries/figures to NA60+ Google Drive with rclone.

This script runs on the user's WSL host. It uses a separately configured
rclone Drive remote; ChatGPT's Drive connector does NOT authenticate rclone.
Default is dry-run; pass --upload for actual transfer. Never uses 'sync',
never deletes remote files and never overwrites a different existing file.

Remote root is anchored to the exact already-created avalanche_3p3_runs
folder ID (so path lookup does not depend on My Drive/shared-drive roots).
"""
from __future__ import annotations
import argparse
import hashlib
import json
import re
import shutil
import subprocess
from datetime import datetime,timezone
from pathlib import Path

HERE=Path(__file__).resolve().parent
REPO=HERE.parents[1]
BUNDLE_ROOT=REPO/"test_runs/mwpc_stage_c/drive_publications_v1"
DRIVE_FOLDER_ID="1Xbh45V86BhVu-WLYSM5qtosqPEG363M7"
ALLOWED_EXT={".png",".jpg",".jpeg",".svg",".pdf",".csv",".md",".txt",".log"}
ALLOWED_NAMES={"summary.json","manifest.json","run_config.json",
               "readme.md","README.md","provenance.json"}
MAX_FILE_BYTES=20*1024*1024
MAX_TOTAL_BYTES=100*1024*1024

def sha(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for piece in iter(lambda:f.read(1<<20),b""):h.update(piece)
    return h.hexdigest()

def prepare(source,run_id):
    source=Path(source).resolve()
    if not source.is_dir() or not (source/"summary.json").is_file():
        raise FileNotFoundError("Expected completed run directory with summary.json: "+str(source))
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]{2,95}",run_id) or ".." in run_id:
        raise ValueError("Run ID must be a unique, safe identifier (letters/numbers/_.-)")
    dest=BUNDLE_ROOT/run_id
    if dest.exists():raise FileExistsError("Preserving existing publication bundle: "+str(dest))
    if source==dest or BUNDLE_ROOT in source.parents:
        raise ValueError("Refusing to recursively publish an existing publication bundle")
    files=[]
    nbytes=0
    for item in sorted(source.iterdir()):
        if not item.is_file() or item.is_symlink():continue
        if item.name not in ALLOWED_NAMES and item.suffix.lower() not in ALLOWED_EXT:
            continue
        if item.name.lower().endswith(".gz") or item.name.lower().endswith(".zip"):
            continue
        size=item.stat().st_size
        if size>MAX_FILE_BYTES:continue
        if nbytes+size>MAX_TOTAL_BYTES:
            print("SKIP due to 100MiB publication cap:",item.name)
            continue
        files.append(item);nbytes+=size
    if source/"manifest.json" not in files and (source/"manifest.json").is_file():
        raise ValueError("Original manifest is too large or was unexpectedly excluded")
    if source/"summary.json" not in files:
        raise ValueError("Scientific summary too large for automated publication")
    dest.mkdir(parents=True)
    inventory={}
    for p in files:
        target=dest/p.name
        shutil.copy2(p,target)
        if sha(p)!=sha(target):raise IOError("Copy checksum mismatch: "+p.name)
        inventory[p.name]={"sha256":sha(target),"size_bytes":target.stat().st_size}
    def git_commit():
        return subprocess.check_output(["git","-C",str(REPO),"rev-parse","HEAD"],
                                       text=True).strip()
    manifest={
        "schema":"NA60PLUS_MWPC_DRIVE_PUBLICATION_V1",
        "run_id":run_id,"created_utc":datetime.now(timezone.utc).isoformat(),
        "source_local_path":str(source),
        "source_sha256_files":inventory,"total_bytes":nbytes,
        "git_commit":git_commit(),
        "target_google_drive_folder_id":DRIVE_FOLDER_ID,
        "physics_note":"Original physics files unmodified; not an independent reanalysis",
        "excluded_large_or_raw_files":True,
    }
    (dest/"publication_manifest.json").write_text(
        json.dumps(manifest,indent=2)+"\n",encoding="utf-8")
    print("Bundle:",dest)
    print("Included:",", ".join(inventory))
    print("Excluded by design: raw crossings.csv.gz, large HTML/ZIP, large CSV.")
    return dest

def publish(bundle,remote,upload=False):
    if not re.fullmatch(r"[A-Za-z0-9_-]+",remote):
        raise ValueError("Remote must be a configured rclone remote NAME (not a URL or path)")
    binary=shutil.which("rclone")
    if binary is None:
        raise RuntimeError("rclone not found: install and run 'rclone config' in WSL first")
    remotes=subprocess.check_output([binary,"listremotes"],text=True).splitlines()
    if remote+":" not in remotes:
        raise RuntimeError(f"No rclone remote '{remote}:' configured; available: {remotes}")
    bundle=Path(bundle)
    if not (bundle/"publication_manifest.json").is_file():
        raise ValueError("Publication manifest missing")
    m=json.loads((bundle/"publication_manifest.json").read_text())
    run_id=m["run_id"]
    for filename,props in m["source_sha256_files"].items():
        if sha(bundle/filename)!=props["sha256"]:
            raise ValueError("Publication bundle changed since preparation: "+filename)
    dest=f"{remote}:{run_id}"
    args=[binary,"copy",str(bundle),dest,
          "--drive-root-folder-id",DRIVE_FOLDER_ID,
          "--immutable","--checksum","--stats-one-line"]
    if not upload:
        args.append("--dry-run")
    print("rclone mode:", "UPLOAD" if upload else "DRY RUN")
    print("Remote folder: "+f"https://drive.google.com/drive/folders/{DRIVE_FOLDER_ID}")
    print("Immutable run: "+dest)
    subprocess.run(args,check=True)
    if upload:
        # This is an independent read of files after writing to Google Drive,
        # not an assertion that the user's UI already shows them.
        check=[binary,"check",str(bundle),dest,"--drive-root-folder-id",DRIVE_FOLDER_ID,
               "--one-way"]
        subprocess.run(check,check=True)
        print("PASS: uploaded immutable run and rclone readback check completed")
    else:
        print("DRY RUN only; add --upload to transfer and verify.")
def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--source",type=Path,required=True)
    ap.add_argument("--run-id",required=True)
    ap.add_argument("--remote",required=True,
                    help="Name from 'rclone listremotes', e.g. gdrive (NOT credentials)")
    ap.add_argument("--upload",action="store_true",
                    help="Without this, only copy to local staging and rclone dry-run")
    a=ap.parse_args()
    bundle=prepare(a.source,a.run_id)
    publish(bundle,a.remote,a.upload)

if __name__=="__main__":main()
