#!/usr/bin/env bash
# Large-sample K2A campaign at the frozen 2100 V cathode geometry. Plan: K2A_CAMPAIGN_PLAN.md.
# ONLY Phase 1 (read-only preflight) is enabled. Phases 2-4 refuse to run until the preflight has been inspected and
# the next step is explicitly approved (they need new tooling and frozen acceptance criteria first).
#
#   bash run_k2a_campaign.sh selftest WORK                 SYNTHETIC tests of the preflight tools (no ROOT, no archive)
#   bash run_k2a_campaign.sh inventory                     1a  slot / seal / SHA / inventory / role check   -> $PF/inventory
#   bash run_k2a_campaign.sh moments  PF                   1b  batch birth-record moments (frozen reader, keep-going,
#                                                              chunks of $CHUNK; sealed set in a separate directory)
#   bash run_k2a_campaign.sh domain   PF                   1c  kernel-domain coverage, causes, extension cost -> $PF/domain
#   bash run_k2a_campaign.sh preflight                     1a + 1b + 1c in one new $PF
#   bash run_k2a_campaign.sh recover PF                    READ-ONLY, input metadata: D anatomy, recentering recovery,
#                                                              radial check, K256 fraction, provenance -> $PF/recenter_<ts>
#   bash run_k2a_campaign.sh supplementary PF RECOVER [N]  deterministic non-sealed phi/omega development selection
#   bash run_k2a_campaign.sh recenter-refs [MUON ...]      E1 + E2a on the reference muons (recenter/ experiment; no
#                                                              new transport; frozen replay re-run on relabelled exports)
#   bash run_k2a_campaign.sh recenter-far PF INDEX ...     E2b on non-sealed far-D muons (frozen export + T1 + grouped S0)
#   bash run_k2a_campaign.sh phase2|phase3|phase4          refused (see plan)
#
# READ-ONLY on: the archive, the selection, the validated S-ladder run (kernel, Green table, references), the repo.
# Writes only below $OUTROOT/preflight_<timestamp>/ (umask 077). No S0 export/replay, no avalanche, no Geant4, no PBS.
# Environment (defaults = the paths used by run_s_ladder.sh):
#   MWPC_ARCHIVE_2100_K256  archive root;  MWPC_SIGNAL40_OUTPUT  frozen selection directory
#   MWPC_K2A_REFERENCE_RUN  validated run (default .../s_ladder_2100/run_20261124_124740);  MWPC_K2A_TAG (r14t72edge5)
#   MWPC_BIRTHS_BIN         p3_birth_moments built from THIS commit (needs --keep-going and .event.json)
#   MWPC_K2A_OUT            output root;  MWPC_K2A_CHUNK (100 events per converter process);  MWPC_K2A_SHA_JOBS (2)
set -Eeo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../../../.." && pwd)"
set -u
umask 077
ARCHIVE="${MWPC_ARCHIVE_2100_K256:-$HOME/na6p/cathode-overnight-2100/test_runs/mwpc_stage_c/overnight_2100_weighted_k256_v1}"
SELECTION="${MWPC_SIGNAL40_OUTPUT:-$HOME/NA6PRoot-signal40-2100/test_runs/mwpc_stage_c/phase_space_2100_signal40_v1}"
RUN="${MWPC_K2A_REFERENCE_RUN:-$REPO/test_runs/mwpc_stage_c/s_ladder_2100/run_20261124_124740}"
TAG="${MWPC_K2A_TAG:-r14t72edge5}"
OUTROOT="${MWPC_K2A_OUT:-$REPO/test_runs/mwpc_stage_c/k2a_campaign_2100}"
CHUNK="${MWPC_K2A_CHUNK:-100}"; SHAJOBS="${MWPC_K2A_SHA_JOBS:-2}"; NSLOTS="${MWPC_K2A_SLOTS:-1000}"
KD="$RUN/kernel_$TAG"
fail(){ echo "K2A CAMPAIGN FAILED: $*" >&2; exit 2; }
need_dir(){ [[ -d "$1" ]] || fail "missing $2: $1"; }

kernel_id(){
  python3 - "$KD" <<'PY'
import hashlib, sys
from pathlib import Path
d = Path(sys.argv[1]); h = hashlib.sha256()
for f in ("K.npy", "KPHI0.npy", "kernel.json"):
    s = hashlib.sha256()
    with open(d / f, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 24), b""):
            s.update(b)
    h.update(s.hexdigest().encode())
print(h.hexdigest())
PY
}

do_inventory(){
  PF="$1"; need_dir "$ARCHIVE" archive; need_dir "$SELECTION" selection
  python3 "$HERE/k2a_inventory.py" --archive "$ARCHIVE" --selection "$SELECTION" --reference-run "$RUN" \
      --sha-jobs "$SHAJOBS" --slots "$NSLOTS" --out "$PF/inventory" | tee "$PF/inventory.log" || fail "inventory check (see $PF/inventory/inventory_check.json)"
}

do_moments(){
  PF="$1"; [[ -x "${MWPC_BIRTHS_BIN:-}" ]] || fail "set MWPC_BIRTHS_BIN to p3_birth_moments built from this commit"
  grep -aq -- "--keep-going" "$MWPC_BIRTHS_BIN" && grep -aq "n_primaries_zero_ions" "$MWPC_BIRTHS_BIN" ||
    fail "$MWPC_BIRTHS_BIN predates this commit (no --keep-going / event.json): rebuild it"
  [[ -s "$PF/inventory/events_nonsealed.txt" ]] || fail "run inventory first"
  for set in nonsealed sealed; do
    L="$PF/inventory/events_$set.txt"; O="$PF/moments_$set"; [[ ! -e "$O" ]] || fail "$O exists"; mkdir -p "$O/lists"
    split -d -a 3 -l "$CHUNK" "$L" "$O/lists/chunk_"
    for c in "$O"/lists/chunk_*; do
      n=${c##*/}; t0=$(date +%s)
      set +e
      if command -v /usr/bin/time >/dev/null; then
        /usr/bin/time -v -o "$O/$n.time.txt" "$MWPC_BIRTHS_BIN" --batch "$c" "$O/$n" --keep-going > "$O/$n.log" 2>&1; rc=$?
      else
        "$MWPC_BIRTHS_BIN" --batch "$c" "$O/$n" --keep-going > "$O/$n.log" 2>&1; rc=$?
      fi
      set -e
      echo "$n rc=$rc events=$(wc -l < "$c") wall_s=$(( $(date +%s) - t0 ))" | tee -a "$O/chunks.txt"
      # rc 3 = some events failed (listed as *.error.txt); any other non-zero = the process died: those events show up
      # as NOT_CONVERTED in the domain check, never silently dropped
    done
  done
  grep -h "FAILED" "$PF"/moments_nonsealed/chunk_*.log | head -20 || true
}

do_domain(){
  PF="$1"; need_dir "$KD" "validated kernel"; [[ -d "$PF/moments_nonsealed" ]] || fail "run moments first"
  KID=$(kernel_id); echo "kernel_id=$KID" > "$PF/kernel_id.txt"
  VAL=$( (ls -d "$RUN"/engine_"${TAG}"_* 2>/dev/null || true) | sort | tail -1)
  if [[ -n "$VAL" && -s "$VAL/kernel_id.txt" ]]; then
    grep -qx "kernel_id=$KID" "$VAL/kernel_id.txt" || fail "kernel $KD differs from the one validated in $VAL"
    echo "kernel id matches the validated engine run $VAL" | tee -a "$PF/kernel_id.txt"
  else
    echo "WARNING no validated engine_${TAG}_* run found to compare the kernel id" | tee -a "$PF/kernel_id.txt"
  fi
  python3 "$HERE/k2a_domain_check.py" --kernel "$KD" --kernel-id "$KID" --inventory "$PF/inventory" --selection "$SELECTION" \
      --moments "$PF/moments_nonsealed" --sealed-moments "$PF/moments_sealed" --reference-run "$RUN" --tag "$TAG" \
      --out "$PF/domain" > "$PF/domain.log" || fail "domain check"
  ( cd "$PF" && find inventory domain -type f ! -path "*/SEALED_DO_NOT_OPEN/*" -print0 | sort -z | xargs -0 sha256sum > PREFLIGHT_SHA256SUMS )
  git -C "$REPO" rev-parse HEAD > "$PF/SOURCE_GIT_SHA.txt"
  echo "return for review: $PF/inventory/inventory_check.json $PF/inventory/excluded_slots_detail.md"
  echo "                   $PF/domain/coverage.md $PF/domain/coverage.json $PF/domain/extension_estimate.json"
  echo "                   $PF/moments_*/chunks.txt  (do NOT open $PF/domain/SEALED_DO_NOT_OPEN)"
}

case "${1:-}" in
  selftest) W="${2:?WORK dir (new)}"; python3 "$HERE/tests/test_preflight_synthetic.py" "$W/tools"
    bash "$HERE/tests/test_runner_synthetic.sh" "$W/tools" "$W/runner";;
  inventory) PF="$OUTROOT/preflight_$(date +%Y%m%d_%H%M%S)"; mkdir -p "$PF"; do_inventory "$PF"; echo "PF=$PF";;
  moments) do_moments "${2:?PF dir}";;
  domain) do_domain "${2:?PF dir}";;
  preflight)
    PF="$OUTROOT/preflight_$(date +%Y%m%d_%H%M%S)"; mkdir -p "$PF"; echo "PF=$PF"
    do_inventory "$PF"; do_moments "$PF"; do_domain "$PF";;
  recover)
    PF="${2:?PF dir}"; need_dir "$KD" "validated kernel"; O="$PF/recenter_$(date +%Y%m%d_%H%M%S)"
    python3 "$HERE/k2a_recenter_recovery.py" --pf "$PF" --kernel "$KD" --out "$O" > "$O.log" || fail "recovery ($O.log)"
    echo "return: $O/recovery.md $O/recovery.json (per event: $O/per_event_recenter.csv)";;
  supplementary)
    PF="${2:?PF dir}"; RECOVER="${3:?recover output dir}"; N="${4:-8}"; O="$PF/supplementary_$(date +%Y%m%d_%H%M%S)"
    python3 "$HERE/k2a_supplementary_selection.py" --pf "$PF" --recenter "$RECOVER" --out "$O" --per-channel "$N" > "$O.log" || fail "selection"
    echo "return: $O/supplementary_dev.md $O/supplementary_dev.csv $O/SHA256SUMS";;
  recenter-refs)
    shift; MWPC_K2A_OUT="$OUTROOT" bash "$HERE/recenter/k2a_recenter_experiment.sh" refs "$RUN" "$TAG" "$@";;
  recenter-far)
    PF="${2:?PF dir}"; shift 2
    MWPC_K2A_PF="$PF" MWPC_K2A_OUT="$OUTROOT" bash "$HERE/recenter/k2a_recenter_experiment.sh" far "$RUN" "$TAG" "$@";;
  phase2|phase3|phase4)
    fail "$1 is not enabled: inspect the preflight, then approve the next step (K2A_CAMPAIGN_PLAN.md sections 5-7)";;
  *) sed -n 2,26p "$0"; exit 1;;
esac
