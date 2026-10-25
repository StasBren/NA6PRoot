#!/usr/bin/env bash
# Dedicated HIGH-STAT 500 independent Garfield microscopic trajectories per
# u/v box, per config. Uses the SAME two p3mm x FW 0/75um geometries
# as the original frozen 20-trial study, but a DIFFERENT immutable output
# root and RNG seed. Cannot silently change 6M production's frozen LUTs.
#
# USAGE
#   bash test/mwpc-stage-c/run_pi_pitch_fieldwire_500.sh check
#   bash test/mwpc-stage-c/run_pi_pitch_fieldwire_500.sh start
#   bash test/mwpc-stage-c/run_pi_pitch_fieldwire_500.sh status
#   bash test/mwpc-stage-c/run_pi_pitch_fieldwire_500.sh verify
#
# Launches the existing validated Garfield runner at low CPU priority and
# saves full, per-u-column immutable checkpoints. Stop/resume via the SAME
# 'start' command: already completed columns are checked and skipped.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd -P)"
REPO="$(cd "$HERE/../.." && pwd -P)"
BASE="$REPO/test_runs/mwpc_stage_c/pi_pitch_fieldwire_500_p3_v1"
OUTPUT="$BASE/pilot"
PIDFILE="$BASE/runner.pid"
LOG="$BASE/runner.log"
SAMPLES=500
SEED=260500
FW_UM=75
MODE="$1"
if [[ "$#" -ne 1 || ! "$MODE" =~ ^(check|start|status|verify)$ ]]; then
  echo "Usage: bash $0 check|start|status|verify" >&2
  exit 2
fi
mkdir -p "$BASE"

pid_running() {
  [[ -f "$PIDFILE" ]] || return 1
  local pid
  pid="$(cat "$PIDFILE")"
  [[ "$pid" =~ ^[0-9]+$ ]] || return 1
  kill -0 "$pid" 2>/dev/null || return 1
  # Avoid a stale PID file accidentally matching some unrelated process.
  ps -p "$pid" -o args= 2>/dev/null | grep -q 'run_pi_pitch_fieldwire_pilot.sh pilot 75 50'
}

if [[ "$MODE" == start ]]; then
  if pid_running; then
    echo "Already running with PID=$(cat "$PIDFILE")"
    echo "Follow the log: tail -f '$LOG'"
    exit 0
  fi
  if [[ -f "$OUTPUT/two_pitch_configurations.json" ]]; then
    echo "Final two-geometry summary exists. Run verify for integrity."
    exit 0
  fi
  if [[ ! -d "$REPO/test_runs/mwpc_stage_c/pi_pitch_fieldwire_pilot_v1/pilot" ]]; then
    echo "WARNING: original 20-trial baseline directory not found."
  fi
  echo "Launching 2 new Garfield LUTs: 500 independent microscopic drifts per cell,"
  echo "900,000 total single-electron Garfield trajectories, B=0, p=3mm, FW=0/75um."
  echo "Original 20-trial physics LUTs and active 6M Heed/Garfield fold unchanged."
  echo "New frozen maps: $OUTPUT"
  echo "Important: this is a SECOND CPU workload, run at nice +15."
  nohup env \
    NA60_PI_PITCH_LUT_OUT="$BASE" \
    NA60_PI_PITCH_LUT_SEED="$SEED" \
    NA60_PI_PITCH_LUT_ONLY_PITCH="3.0" \
    nice -n 15 bash "$HERE/run_pi_pitch_fieldwire_pilot.sh" \
      pilot "$FW_UM" "$SAMPLES" > "$LOG" 2>&1 < /dev/null &
  child=$!
  echo "$child" > "$PIDFILE"
  echo "Started Garfield 50-trial launcher PID=$child"
  echo "Status: bash '$HERE/run_pi_pitch_fieldwire_500.sh' status"
  echo "Log: tail -f '$LOG'"
  exit 0
fi

python3 - "$HERE" "$OUTPUT" "$MODE" "$SAMPLES" "$SEED" <<'PY'
import csv
import hashlib
import json
import sys
from pathlib import Path

here, root, mode, samples, seed = (
    Path(sys.argv[1]),Path(sys.argv[2]),sys.argv[3],
    int(sys.argv[4]),int(sys.argv[5]))
sys.path.insert(0,str(here))
from pi_pitch_fieldwire_pilot import (
    tag,DU,NV,validate_column,load_column_timing
)

expected_columns=2*round(3.0/DU)
expected_trials=NV*samples*expected_columns
complete_cols=0
total_elapsed=0.
verified_all=True
summary_check=[]
print("FROZEN p=3 mm / 500-trial Garfield LUT root:",root)
print("Samples/box:",samples,"expected columns:",expected_columns,
      "expected full Garfield electron trials:",expected_trials)
for pitch in (3.0,):
    for fw in (0.,75.):
        name=tag(pitch,fw)
        folder=root/name
        n_u=round(pitch/DU)
        manifest_path=folder/"frozen_config.json"
        files=list((folder/"columns").glob("iu_*.csv"))
        complete_cols+=len(files)
        expected_per_config=n_u*NV*samples
        config_ready=False
        if manifest_path.exists():
            m=json.loads(manifest_path.read_text(encoding="utf-8"))
            checks={
              "study":"PI_GARFIELD_PITCH_FIELDWIRE_V1",
              "pitch_mm":pitch,
              "field_wire_diam_um":fw,
              "samples_per_box":samples,
              "seed":seed,
              "mode":"pilot",
              "magnetic_field_t":0.,
              "gas":"Ar/CO2 70:30, 293.15 K, 760 Torr",
              "gap_minus_mm":3.,"gap_plus_mm":3.,
              "anode_voltage_v":1800.,"anode_diameter_um":30.,
              "du_mm":.1,"dv_mm":.2,"n_u":n_u,"n_v":30,
            }
            if any(m.get(k)!=v for k,v in checks.items()):
                raise ValueError("Inconsistent 500-trial Garfield frozen manifest: "+str(folder))
            config_ready=True
        elif files:
            raise ValueError("Columns present without immutable physics manifest: "+str(folder))
        print(f"  {name:22s}: {len(files):2d}/{n_u} columns "
              f"({len(files)*NV*samples:,}/{expected_per_config:,} trial rows); "
              f"manifest={'OK' if config_ready else 'not started'}")
        if mode=="verify":
            if not config_ready or len(files)!=n_u:
                raise ValueError("Incomplete 500-trial p3 Garfield data: "+name)
            for i in range(n_u):
                path=folder/"columns"/f"iu_{i:03d}.csv"
                validate_column(path,i,pitch,samples)
                wall=load_column_timing(path)
                if wall is None:
                    raise ValueError("Missing SHA-checked runtime sidecar: "+str(path))
                total_elapsed+=wall
            final=folder/"garfield_electrons.csv"
            if not final.is_file():
                raise ValueError("Missing assembled Garfield microtrajectories: "+str(final))
            count=0
            with final.open(newline="",encoding="utf-8") as f:
                for _ in csv.DictReader(f):
                    count+=1
            if count!=expected_per_config:
                raise ValueError(f"Corrupted {name}: assembled rows {count}/{expected_per_config}")
            summ=json.loads((folder/"summary.json").read_text(encoding="utf-8"))
            if (summ["n_electrons"]!=expected_per_config or
                summ["n_anode_collected"]+
                summ["n_field_wire_collected"]+
                summ["n_not_on_wires"]!=expected_per_config):
                raise ValueError("Summary conservation error: "+name)
            summary_check.append(name)
print(f"TOTAL: {complete_cols}/{expected_columns} completed checkpoint columns; "
      f"{complete_cols*NV*samples:,}/{expected_trials:,} saved raw microtrajectories "
      "(not counting in-progress .tmp)")
if mode=="verify":
    meta=root/"two_pitch_configurations.json"
    if not meta.exists():
        raise ValueError("Missing final two-geometry physics summary")
    final=json.loads(meta.read_text(encoding="utf-8"))
    if len(final.get("results",[]))!=2 or len(summary_check)!=2 or \
       {(r["pitch_mm"],r["field_wire_diameter_um"]) for r in final["results"]} != {(3.0,0.0),(3.0,75.0)}:
        raise ValueError("Incomplete two-geometry final summary")
    print("PASS: TWO p3 LUTs, 500 trajectories/box in each, 900,000 ""
          "raw microscopic Garfield trajectories; independent SHA-checked "
          "column checkpoints and electron conservation.")
    print(f"SUM measured independent column wall time: {total_elapsed/3600:.3f} h "
          "(this is NOT concurrent elapsed runtime)")
elif mode=="check":
    print("PASS: frozen 500-trial physics identifiers valid where checkpoints exist. "
          "No Garfield execution started.")
else:
    print("STATUS: Garfield physics follows the per-u-column progress above.")
PY

if [[ "$MODE" == status && -f "$PIDFILE" ]]; then
  if pid_running; then
    echo "RUNNING: launcher PID=$(cat "$PIDFILE")"
  else
    echo "NOT RUNNING: old launcher PID=$(cat "$PIDFILE")"
  fi
fi
if [[ "$MODE" == status && -f "$LOG" ]]; then
  echo "--- Last 12 log lines ---"
  tail -n 12 "$LOG"
fi
