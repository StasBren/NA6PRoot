#!/usr/bin/env bash
# ROOT/C++ physical per-wire analysis of the 1800 V ~8k completed crossings.
# Only Python use here: SHA-256 verification of original input/4 audited CSVs
# and creation of immutable path list; ALL physics/statistics is C++/ROOT.
# Usage: bash run_wiretruth1800_analysis.sh build|run|submit|status
set -Eeo pipefail
export GARFIELD_INSTALL=${GARFIELD_INSTALL:-/cvmfs/sft-nightlies.cern.ch/lcg/nightlies/dev3/Mon/Garfield++/HEAD/x86_64-el9-gcc13-opt}
if [[ -f "$GARFIELD_INSTALL/Garfield++-env.sh" ]];then
  source "$GARFIELD_INSTALL/Garfield++-env.sh"
fi
set -u
SCRIPT="$(readlink -f "${BASH_SOURCE[0]}")"
HERE="$(dirname "$SCRIPT")"
REPO="$(cd "$HERE/../../.." && pwd)"
SRC="$HERE/wiretruth1800_distributions.cc"
BUILD="$REPO/test/mwpc-stage-c/build-wiretruth1800-analysis"
BIN="$BUILD/wiretruth1800_distributions"
ROOT="${MWPC_1800_ROOT:-/storage/cgrp/stanis/mwpc_penning057_wiretruth30k_v1}"
OLD="${MWPC_1800_OLD_ROOT:-/storage/cgrp/stanis/mwpc_penning057}"
OUT="${MWPC_1800_ANALYSIS_OUT:-$REPO/test_runs/wiretruth1800_8k_snapshot_$(date -u +%Y%m%dT%H%M%SZ)}"
MODE="${1:-status}"
# WIPP 'N' was confirmed by successful earlier 1800V production.
# Unlike the inaccessible normE queue, this queue may be used by the
# original MWPC project account; user can override after checking qstat -Qf.
PBS_QUEUE="${MWPC_1800_PBS_QUEUE:-N}"
[[ "$PBS_QUEUE" =~ ^[A-Za-z][A-Za-z0-9_-]*$ ]] || {
  echo "Invalid PBS queue name: $PBS_QUEUE" >&2;exit 2;
}
CHUNK_LIMIT="${MWPC_1800_CHUNKS_PER_CHANNEL:-0}"
[[ "$CHUNK_LIMIT" =~ ^[0-9]+$ ]] || { echo "Invalid MWPC_1800_CHUNKS_PER_CHANNEL" >&2;exit 2; }
fail(){ echo "1800 V analysis ERROR: $*" >&2;exit 2;}
check(){
  [[ "$ROOT" == "/storage/cgrp/stanis/mwpc_penning057_wiretruth30k_v1" ]] ||
    fail "Unknown source production root. Do not mix scientific campaigns."
  [[ -s "$ROOT/locked/selection.json" && -s "$ROOT/locked/binary.json" ]] ||
    fail "Missing locked original 1800V source/physics identity"
  [[ ! -e "$OUT" ]] || fail "Analysis OUT exists: $OUT. Choose a fresh folder, never overwrite."
}
build(){
  command -v root-config >/dev/null || fail "ROOT not loaded. Use original CERN CVMFS environment."
  command -v c++ >/dev/null || fail "C++ compiler missing."
  mkdir -p "$BUILD"
  echo "C++ compiler=$(c++ -dumpfullversion) ROOT=$(root-config --version)"
  c++ -std=c++17 -O2 -Wall -Wextra \
    $(root-config --cflags) "$SRC" \
    $(root-config --libs) -o "$BIN"
  [[ -x "$BIN" ]] || fail "Compilation did not produce ROOT/C++ executable."
  echo "PASS compiled $BIN; no Geant4/Heed/Garfield jobs launched."
}
run(){
  check
  [[ -x "$BIN" ]] || fail "First run build mode (no compilation inside PBS)."
  mkdir -p "$(dirname "$OUT")"
  SNAP="$OUT.frozen_chunks.tsv"
  [[ ! -e "$SNAP" ]] || fail "Snapshot already exists: $SNAP. Choose another output prefix."
  python3 - "$ROOT" "$OLD" "$SNAP" "$CHUNK_LIMIT" <<'PY'
import sys, os, json, hashlib
from pathlib import Path
root,old,dst=map(Path,sys.argv[1:4])
chunk_limit=int(sys.argv[4])
assert chunk_limit>=0
selection=json.loads((root/"locked/selection.json").read_text())
binary=json.loads((root/"locked/binary.json").read_text())
schema="NA60_WIRETRUTH30K_PENNING057_PRODUCTION_V1"
expected_git="599e29522edb6fcec9dfcca2904308cb3186edce"
assert selection["schema"]==binary["schema"]==schema
assert selection["git_commit"]==binary["git_commit"]==expected_git
assert selection["cxx_sha256"]==binary["cxx_sha256"]
expected_bin=binary["binary_sha256"]
assert len(expected_bin)==64 and all(c in "0123456789abcdef" for c in expected_bin)
print("LOCKED original 1800V binary sha256",expected_bin,flush=True)
phys=selection["physics"]
assert phys["anode_voltage_V"]==1800.0
assert phys["pitch_mm"]==3.0 and phys["penning_r"]==0.57
assert selection["input_strategy"]=="direct_read_original2m_source_chunk_lazy_sha256_verify"
assert root.resolve()!=old.resolve()
def hash_file(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda:f.read(1024*1024),b""):
            h.update(chunk)
    return h.hexdigest()
parts={"wire_muon":"microscopic_wire_muon.csv",
       "wire_seed":"microscopic_wire_seed.csv",
       "crossing":"microscopic_crossing.csv",
       "seed":"microscopic_seed.csv"}
counts={x:0 for x in ("jpsi","phi","omega")}
partial=0
with dst.open("x") as out:
    out.write("channel\tchunk\toriginal_input\twire_muon\twire_seed\tcrossing\tseed\n")
    for channel in ("Jpsi","Phi","Omega"):
        ch=channel.lower()
        srcdir=Path(selection["source_campaigns"][channel]["original_chunk_dir"])
        assert srcdir.resolve()==(old/f"original2m_{ch}_10000_v1/chunks").resolve()
        for i in range(500):
            if chunk_limit and counts[ch]>=chunk_limit:
                break
            d=root/"results"/ch/f"chunk_{i:04d}"
            receipt=d/"verified_report.json"
            if not receipt.exists():
                partial+=d.is_dir()
                continue
            x=json.loads(receipt.read_text())
            assert (x["schema"],x["channel"],x["chunk"])==(schema,channel,i),receipt
            assert x["git_commit"]==expected_git and x["binary_sha256"]==expected_bin,receipt
            assert x["cxx_sha256"]==selection["cxx_sha256"],receipt
            assert x["counts"]["muons"]==20 and x["counts"]["capped"]==0,receipt
            assert x["counts"]["failed_return"]==0,receipt
            src=srcdir/f"chunk_{i:04d}.csv"
            assert x["input_sha256"]==selection["input_sha256"][f"{channel}:{i:04d}"],src
            assert hash_file(src)==x["input_sha256"],src
            paths=[]
            for key,name in parts.items():
                f=d/name
                record=x["files"][key]
                assert f.is_file() and f.stat().st_size==record["bytes"],f
                assert hash_file(f)==record["sha256"],f
                paths.append(str(f))
            out.write("\t".join([ch,str(i),str(src),*paths])+"\n")
            counts[ch]+=1
            if counts[ch]%25==0:
                print(f"SHA-VERIFIED {channel}: {counts[ch]} complete 20-muon chunks",flush=True)
    out.flush()
    os.fsync(out.fileno())
total=sum(counts.values())
assert total>=1,"No verified source chunks"
print("FROZEN SNAPSHOT",counts,"chunks=",total,"muons=",20*total,
      "per_channel_limit=",chunk_limit,"partial dirs=",partial,
      "manifest=",dst,flush=True)
PY
  echo "ROOT/C++ analysis of original 1800V physical per-wire distributions"
  "$BIN" "$SNAP" "$OUT"
  cp "$SNAP" "$OUT/frozen_chunk_manifest.tsv"
  {
    echo "analysis_git_commit=$(git -C "$REPO" rev-parse HEAD)"
    echo "per_channel_chunk_limit=$CHUNK_LIMIT"
    echo "original_physics_git=599e29522edb6fcec9dfcca2904308cb3186edce"
    echo "original_physics_binary_sha256=$(python3 -c 'import json,sys;print(json.load(open(sys.argv[1]))["binary_sha256"])' "$ROOT/locked/binary.json")"
    echo "pbs_queue=$PBS_QUEUE"
    echo "root_version=$(root-config --version)"
    echo "gcc_version=$(c++ -dumpfullversion)"
  } > "$OUT/provenance.txt"
  (
    cd "$OUT"
    sha256sum *.root *.csv *.png *.txt *.tsv > SHA256SUMS
    sha256sum -c SHA256SUMS >/dev/null
  )
  echo "PASS sealed $OUT"
}
submit(){
  check
  [[ -x "$BIN" ]] || fail "build before submitting (source-only ROOT compile on login)."
  command -v /opt/pbs/bin/qsub >/dev/null || fail "OpenPBS qsub missing."
  mkdir -p "$REPO/test_runs/pbs-logs"
  [[ "$OUT" == "$REPO"/* ]] || fail "Choose unique OUT inside new analysis worktree."
  [[ ! -e "$OUT.frozen_chunks.tsv" ]] || fail "Snapshot exists; choose a fresh output."
  local_jobname=w1800ana
  local_walltime=04:00:00
  if [[ "$CHUNK_LIMIT" == "1" ]];then
    local_jobname=w1800pil
    local_walltime=01:00:00
  fi
  echo "Using PBS queue=$PBS_QUEUE, name=$local_jobname, mem=4gb, walltime=$local_walltime, io=0.1"
  JOB=$(/opt/pbs/bin/qsub -q "$PBS_QUEUE" -N "$local_jobname" \
    -l select=1:ncpus=1:mem=4gb -l "walltime=$local_walltime" -l io=0.1 \
    -v "MWPC_1800_ROOT=$ROOT,MWPC_1800_OLD_ROOT=$OLD,MWPC_1800_ANALYSIS_OUT=$OUT,MWPC_1800_ANALYSIS_REPO=$REPO,MWPC_1800_CHUNKS_PER_CHANNEL=$CHUNK_LIMIT,MWPC_1800_PBS_QUEUE=$PBS_QUEUE" \
    -o "$REPO/test_runs/pbs-logs" -e "$REPO/test_runs/pbs-logs" \
    "$HERE/pbs_wiretruth1800_analysis.sh")
  echo "PBS_SUBMITTED $JOB"
  echo "OUTPUT_AFTER_SUCCESS $OUT"
  printf '%s\t%s\t%s\n' "$(date -u +%FT%TZ)" "$JOB" "$OUT" \
    >> "$REPO/test_runs/pbs-logs/submissions.tsv"
}
case "$MODE" in
  build) build;;
  run) run;;
  submit) submit;;
  pilot)
    [[ "$CHUNK_LIMIT" == "0" ]] || fail "Set no external chunk limit when using pilot"
    CHUNK_LIMIT=1
    submit
    ;;
  status)
    echo "Original 1800V production root: $ROOT"
    echo "Output requested: $OUT"
    if [[ -e "$OUT/SHA256SUMS" ]];then
      (cd "$OUT" && sha256sum -c SHA256SUMS)
    else echo "No sealed analysis output at that path (yet).";fi
    ;;
  *) fail "Usage: bash run_wiretruth1800_analysis.sh {build|run|pilot|submit|status}";;
esac
