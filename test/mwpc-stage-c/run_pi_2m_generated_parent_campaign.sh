#!/usr/bin/env bash
# Stage-C 2M GENERATED dimuon parents per source, full 298-chamber geometry.
# benchmark = 2x50k/source in SAME frozen suite as full = 40x50k/source.
# Outputs remain preserved, validated and resumable per native ROOT chunk.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd -P)"
REPO="$(cd "$HERE/../.." && pwd -P)"
cd "$REPO"
MODE="plan"
if [[ $# -gt 0 ]]; then MODE="$1"; shift; fi
if [[ $# -gt 0 ]]; then
  echo "Usage: bash $0 plan|env-check|status|benchmark|full" >&2
  exit 2
fi
case "$MODE" in plan|env-check|status|benchmark|full) ;; *)
  echo "Usage: bash $0 plan|status|benchmark|full" >&2
  exit 2 ;;
esac

TARGET=2000000
CHUNK=50000  # Existing native Geant4 runner maximum (NOT 100k).
WORKERS="${PI_2M_WORKERS:-4}"
HOURS="${PI_2M_HOURS:-2.5}"
MIN_GIB="${PI_2M_MIN_FREE_GIB:-20}"
SUITE="$REPO/test_runs/mwpc_stage_c/pi_geant4_2m_parent_campaign_v1"

python3 - "$TARGET" "$CHUNK" "$WORKERS" "$HOURS" "$MIN_GIB" "$SUITE" <<'PY'
import math,sys
n,chunk,workers,hours,minimum,suite=sys.argv[1:]
n,chunk,workers=int(n),int(chunk),int(workers)
hours,minimum=float(hours),float(minimum)
if n!=2000000 or chunk!=50000 or not 1<=workers<=32 or (
    not math.isfinite(hours) or hours<0 or not math.isfinite(minimum) or minimum<0):
    sys.exit("ERROR: unsafe 2M campaign parameters")
print(f"Target: {n:,} GENERATED parents/channel x 3 channels; "
      f"{n*3:,} total, {n//chunk} chunks/source x {chunk:,} parents.")
print("Benchmark: first TWO 50k chunks/source = 100k/source, SAME frozen "
      "suite as full 2M/source; no re-generation of completed chunks.")
print(f"workers={workers}; per-invocation wall budget={hours} h, "
      "checked only BETWEEN complete chunks.")
print(f"Frozen campaign suite: {suite}")
print("NO Geant4 quadrant-only simplification. Garfield LUT fold is separate.")
PY

if [[ "$MODE" == "plan" ]]; then exit 0; fi
# Small, read-only environment check. Does NOT create a frozen campaign
# manifest and does NOT start Geant4 or even the 50k benchmark.
if [[ "$MODE" == "env-check" ]]; then
  NA6P_STAGEC_REQUIRE_ROOT_HADD=1
  source "$HERE/activate_pi_stagec_env.sh"
  for program in root hadd na6psim na6psim_parallel; do
    command -v "$program" || {
      echo "ERROR: missing executable $program" >&2
      exit 127
    }
  done
  echo "PASS: native simulator + ROOT/hadd dynamic-link environment ready."
  echo "No simulation and no ROOT campaign files created by env-check."
  exit 0
fi
if [[ "$MODE" == "status" ]]; then
  python3 - "$SUITE" <<'PY'
from pathlib import Path
import json,sys
suite=Path(sys.argv[1])
manifest=suite/"campaign.json"
if not manifest.exists():
    print("Campaign NOT STARTED. No frozen manifest.");sys.exit(0)
frozen=json.loads(manifest.read_text())["frozen_inputs"]
if (frozen["events_per_source"],frozen["chunk_size"])!=(2000000,50000):
    sys.exit("ERROR: frozen suite is not the expected 2M / 50k campaign")
total=0
for ch,seed0 in (("Jpsi",202610181),("Phi",202710181),
                 ("Omega",202810181)):
    done=0;elapsed=0.;clocked=0;bytes_total=0
    for i in range(40):
        path=suite/ch.lower()/f"chunk{i+1:04d}"/"completed.json"
        if not path.is_file():continue
        v=json.loads(path.read_text())
        if v.get("events")!=50000 or v.get("seed")!=seed0+i*251:
            sys.exit(f"ERROR: inconsistent completed checkpoint {path}")
        done+=1
        chunk_dir=path.parent
        bytes_total+=sum(p.stat().st_size for p in chunk_dir.rglob("*")
                         if p.is_file())
        timer=v.get("wall_seconds",{}).get("chunk_total")
        if timer is not None:
            elapsed+=float(timer);clocked+=1
    total+=done*50000
    print(f"{ch}: {done}/40 completed chunks; {done*50000:,}/2,000,000 "
          f"parents; {clocked} measured chunks, {elapsed/3600:.2f} h observed; "
          f"{bytes_total/(1024**3):.2f} GiB retained original files")
    if done:
        print(f"  CONDITIONAL disk-only linear extrapolation at SAME output-per-"
              f"event: {bytes_total/done*40/(1024**3):.1f} GiB/channel. "
              "Actual output footprint may differ.")
print(f"Completed generated parents: {total:,}/6,000,000")
PY
  exit 0
fi

python3 - "$SUITE" "$MIN_GIB" <<'PY'
from pathlib import Path
import shutil,sys
parent=Path(sys.argv[1]).parent
parent.mkdir(parents=True,exist_ok=True)
available=shutil.disk_usage(parent).free/(1024**3)
required=float(sys.argv[2])
print(f"Available: {available:.2f} GiB, required start guard {required:.2f} GiB")
if available<required:sys.exit("ERROR: insufficient disk, no Geant4 started")
print("NOTE: this is a startup guard only; inspect actual 50k ROOT size "
      "before proceeding to the full 2M/channel campaign.")
PY

python3 -m py_compile "$HERE/run_pi_geant4_campaign.py" "$HERE/test_pi_2m_campaign.py"
python3 "$HERE/test_pi_2m_campaign.py"
if [[ "$MODE" == "full" ]]; then
  python3 - "$SUITE" <<'PY'
from pathlib import Path
import json,sys
suite=Path(sys.argv[1])
for channel,base in (("Jpsi",202610181),
                     ("Phi",202710181),("Omega",202810181)):
    for i in (1,2):
        path=suite/channel.lower()/f"chunk{i:04d}"/"completed.json"
        if not path.is_file():
            sys.exit(f"FULL PAUSED: benchmark not complete, missing {path}. "
                     "Repeat benchmark to finish first 100k parents "
                     "for EACH source.")
        info=json.loads(path.read_text())
        if info.get("events")!=50000 or info.get("seed")!=base+(i-1)*251:
            sys.exit(f"FULL REFUSED: wrong frozen benchmark checkpoint: {path}")
print("PASS: all first 100k generated parents/source have validated "
      "checkpoints. Check disk and observed speed before FULL.")
PY
fi
NA6P_STAGEC_REQUIRE_ROOT_HADD=1
source "$HERE/activate_pi_stagec_env.sh"
MAX=0
if [[ "$MODE" == "benchmark" ]]; then
  MAX=2
  echo "BENCHMARK: at most first two 50k-event chunks for EACH source."
else
  echo "FULL: resume from the exact existing 2M-per-source frozen suite."
fi
python3 "$HERE/run_pi_geant4_campaign.py" \
  --mode full --suite "$SUITE" \
  --events-per-source "$TARGET" --chunk-size "$CHUNK" \
  --workers "$WORKERS" --hours "$HOURS" \
  --min-free-gib "$MIN_GIB" --max-chunks-per-source "$MAX"
if [[ -f "$SUITE/campaign.json" ]]; then
  python3 "$HERE/pi_stagec_performance_report.py" \
    --suite "$SUITE" --out "$SUITE/pi_2m_measured_performance.json"
fi
echo "Progress: bash $HERE/run_pi_2m_generated_parent_campaign.sh status"
