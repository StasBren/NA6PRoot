#!/usr/bin/env bash
# Frozen ALL-completed 1800V physical Garfield wire/asymmetry root shards.
# One-time preparation, bounded OpenPBS N/io0.1 array, strict combined merge.
# Never launches new Garfield, never modifies original physical archives.
set -Eeo pipefail
export GARFIELD_INSTALL="${GARFIELD_INSTALL:-/cvmfs/sft-nightlies.cern.ch/lcg/nightlies/dev3/Mon/Garfield++/HEAD/x86_64-el9-gcc13-opt}"
if [[ -f "$GARFIELD_INSTALL/Garfield++-env.sh" ]];then
  # WIPP CVMFS initialization may reference unset ZSH_NAME.
  source "$GARFIELD_INSTALL/Garfield++-env.sh"
fi
set -u
HERE="$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")"
REPO="$(cd "$HERE/../../.." && pwd)"
PROG="$HERE/full_wiretruth1800_campaign.py"
BIN="$REPO/test/mwpc-stage-c/build-wiretruth1800-analysis/wiretruth1800_distributions"
MERGER="$REPO/test/mwpc-stage-c/build-wiretruth1800-analysis/wiretruth1800_merge_full"
CAMPAIGN="${MWPC_1800_FULL:-}"
QUEUE="${MWPC_1800_PBS_QUEUE:-N}"
PARALLEL="${MWPC_1800_PARALLEL:-2}"
SHARD_CHUNKS="${MWPC_1800_SHARD_CHUNKS:-6}"
fail(){ echo "FULL 1800V: $*" >&2;exit 2; }
check_vars(){
  [[ -n "$CAMPAIGN" ]] || fail 'Export MWPC_1800_FULL to a NEW, unique test_runs folder'
  [[ "$CAMPAIGN" == "$REPO/test_runs/"* ]] ||
    fail "This full analysis must stay in the isolated analysis worktree: $REPO/test_runs/"
  [[ "$QUEUE" =~ ^[A-Za-z][A-Za-z0-9_-]*$ ]] || fail "Bad PBS queue"
  [[ "$PARALLEL" =~ ^[1-9][0-9]*$ ]] || fail "MWPC_1800_PARALLEL must be positive integer"
  [[ "$SHARD_CHUNKS" =~ ^[1-9][0-9]*$ ]] || fail "MWPC_1800_SHARD_CHUNKS must be positive integer"
  (( PARALLEL <= 4 )) || fail "Refusing concurrency >4 on shared WIPP storage; use two"
  (( SHARD_CHUNKS <= 20 )) || fail "Refusing >20 heavy original chunks per PBS shard"
}
build(){
  command -v root-config >/dev/null || fail "ROOT environment not loaded"
  mkdir -p "$(dirname "$BIN")"
  echo "ROOT_VERSION=$(root-config --version) CXX_VERSION=$(c++ -dumpfullversion)"
  c++ -std=c++17 -O2 -Wall -Wextra \
    $(root-config --cflags) "$HERE/wiretruth1800_distributions.cc" \
    $(root-config --libs) -o "$BIN"
  c++ -std=c++17 -O2 -Wall -Wextra \
    $(root-config --cflags) "$HERE/wiretruth1800_merge_full.cc" \
    $(root-config --libs) -o "$MERGER"
  python3 -m py_compile "$PROG"
  [[ -x "$BIN" && -x "$MERGER" ]] || fail "ROOT/C++ target missing"
  echo "PASS both full archive analyzers compiled; zero new Garfield muons"
}
ensure_binaries(){
  [[ -x "$BIN" && -x "$MERGER" ]] ||
    fail "First run: bash $HERE/full1800.sh build"
}
nshards(){
  python3 - "$CAMPAIGN/locked/campaign.json" <<'PY'
import json,sys
j=json.load(open(sys.argv[1]))
assert j["schema"]=="NA60_1800V_FULL_WIRE_ASYM_SHARDED_V1"
assert j["n_shards"]>=1 and j["n_chunks"]>=1
print(j["n_shards"])
PY
}
prepare(){
  check_vars;ensure_binaries
  [[ ! -e "$CAMPAIGN" ]] || fail "Campaign exists. Never alter frozen snapshot: $CAMPAIGN"
  python3 "$PROG" prepare --campaign "$CAMPAIGN" --per-shard "$SHARD_CHUNKS"
  echo "PASS full 1800V metadata freeze: $CAMPAIGN"
  echo "Next: bash $HERE/full1800.sh submit"
}
submit(){
  check_vars;ensure_binaries
  [[ -s "$CAMPAIGN/locked/campaign.json" ]] || fail "First prepare frozen source campaign"
  local n
  n="$(nshards)"
  mkdir -p "$CAMPAIGN/pbs-logs"
  echo "Submitting all $n frozen 1800V ROOT-only shards; max_parallel=$PARALLEL"
  echo "queue=$QUEUE 1CPU 4GB 8hours io=0.1 per-array-member"
  local job
  job=$(/opt/pbs/bin/qsub -q "$QUEUE" -N w18all \
    -J "0-$((n-1))%$PARALLEL" \
    -l select=1:ncpus=1:mem=4gb -l walltime=08:00:00 -l io=0.1 \
    -v "MWPC_1800_FULL=$CAMPAIGN,MWPC_1800_WORKTREE=$REPO" \
    -o "$CAMPAIGN/pbs-logs" -e "$CAMPAIGN/pbs-logs" \
    "$HERE/pbs_wiretruth1800_full_shard.sh") ||
      fail "PBS rejected full analysis array; no ROOT/Garfield jobs started"
  [[ -n "$job" ]] || fail "Empty PBS job ID"
  printf '%s\t%s\t%s\t%s\n' "$(date -u +%FT%TZ)" "$job" "$n" "$PARALLEL" \
    >> "$CAMPAIGN/pbs-logs/submission_ledger.tsv"
  echo "PBS_SUBMITTED_FULL_ARRAY=$job"
  echo "CAMPAIGN=$CAMPAIGN"
  echo "Monitor: bash $HERE/full1800.sh status"
}
status(){
  check_vars
  python3 "$PROG" status --campaign "$CAMPAIGN"
  echo "PBS USER JOBS:"
  /opt/pbs/bin/qstat -t -u "$USER" || true
}
merge_submit(){
  check_vars;ensure_binaries
  [[ -s "$CAMPAIGN/locked/campaign.json" ]] || fail "Frozen original full campaign is missing"
  [[ ! -e "$CAMPAIGN/merged" ]] || fail "Merge output exists; no overwriting"
  python3 - "$CAMPAIGN" <<'PY'
import json,sys
from pathlib import Path
root=Path(sys.argv[1])
p=json.load(open(root/"locked/campaign.json"))
missing=[j["id"] for j in p["shards"]
         if not (root/"shards"/f"shard_{j['id']:04d}"/"sealed.json").is_file()]
if missing:
    raise SystemExit(f"REFUSE MERGE: {len(missing)} PBS shards are missing: {missing}")
print(f"MERGE_GATE_PASS all {p['n_shards']} sealed shards present "
      f"for {p['n_muons']} original Geant4 muons; full SHA reread done by PBS merger")
PY
  local job
  job=$(/opt/pbs/bin/qsub -q "$QUEUE" -N w18merge \
    -l select=1:ncpus=1:mem=4gb -l walltime=08:00:00 -l io=0.1 \
    -v "MWPC_1800_FULL=$CAMPAIGN,MWPC_1800_WORKTREE=$REPO" \
    -o "$CAMPAIGN/pbs-logs" -e "$CAMPAIGN/pbs-logs" \
    "$HERE/pbs_wiretruth1800_full_merge.sh") ||
      fail "PBS rejected ROOT merge job"
  echo "PBS_SUBMITTED_FINAL_ROOT_MERGE=$job"
}
case "${1:-}" in
  build) build;;
  prepare) prepare;;
  submit|resume) submit;;
  status) status;;
  merge-submit) merge_submit;;
  *) fail 'Usage: bash full1800.sh {build|prepare|submit|status|resume|merge-submit}';;
esac
