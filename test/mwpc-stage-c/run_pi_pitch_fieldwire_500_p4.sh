#!/usr/bin/env bash
# Independent 4.0 mm / +/- 3 mm cathode gap / B=0 Garfield LUT500.
# 40(u) x 30(v) x 500 microscopic trials x two configurations = 1.2M.
# Does NOT touch the existing p3 LUT500, six-geometry LUT20, or 6M folds.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd -P)"
REPO="$(cd "$HERE/../.." && pwd -P)"
BASE="$REPO/test_runs/mwpc_stage_c/pi_pitch_fieldwire_500_p4_v1"
OUTPUT="$BASE/pilot"
LOG="$BASE/runner.log"
PIDFILE="$BASE/runner.pid"
MODE="${1:-}"
if [[ $# -ne 1 || ! "$MODE" =~ ^(check|start|status|verify)$ ]]; then
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
  ps -p "$pid" -o args= 2>/dev/null | grep -q 'run_pi_pitch_fieldwire_pilot.sh pilot 75 500'
}
if [[ "$MODE" == start ]]; then
  if [[ -f "$OUTPUT/two_pitch_configurations.json" ]]; then
    echo "Already assembled: verify the two p4 LUT500 configurations instead of overwriting."
    echo "bash '$HERE/run_pi_pitch_fieldwire_500_p4.sh' verify"
    exit 0
  fi
  if pid_running; then
    echo "Garfield LUT500 p4 already running. PID=$(cat "$PIDFILE")"
    echo "tail -f '$LOG'"
    exit 0
  fi
  echo "Starting NEW p4 LUT500: noFW + grounded 75um FW."
  echo "80 resumable u-columns; 1,200,000 independent Garfield microscopic drifts."
  echo "Isolated output $BASE; B=0; Ar/CO2 70:30; HV=1800 V; gap=3+3 mm."
  echo "Running at nice +15. Avoid concurrent full 6M production if CPU is constrained."
  nohup env NA60_PI_PITCH_LUT_OUT="$BASE" \
    NA60_PI_PITCH_LUT_SEED=260504 \
    NA60_PI_PITCH_LUT_ONLY_PITCH=4 \
    nice -n 15 bash "$HERE/run_pi_pitch_fieldwire_pilot.sh" pilot 75 500 \
    >"$LOG" 2>&1 </dev/null &
  echo "$!" > "$PIDFILE"
  echo "Launcher PID=$(cat "$PIDFILE"). See $LOG"
  exit 0
fi
# Validate checkpoints using the canonical existing Garfield utilities.
python3 - "$HERE" "$OUTPUT" "$MODE" <<'PY'
from pathlib import Path
import csv
import json
import sys
here, root, mode = Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3]
sys.path.insert(0, str(here))
from pi_pitch_fieldwire_pilot import validate_column, load_column_timing, tag
from pi_pitch_fieldwire_pilot import DU, NV
samples, seed, pitch = 500, 260504, 4.0
expected_cols = round(pitch / DU)
per_config = expected_cols * NV * samples
complete = 0
for fw in (0., 75.):
    name = tag(pitch, fw)
    folder = root / name
    meta = folder / "frozen_config.json"
    present = list((folder / "columns").glob("iu_*.csv"))
    complete += len(present)
    if meta.is_file():
        cfg = json.loads(meta.read_text(encoding="utf-8"))
        required = {
          "study":"PI_GARFIELD_PITCH_FIELDWIRE_V1",
          "mode":"pilot", "pitch_mm":4., "field_wire_diam_um":fw,
          "samples_per_box":samples,"seed":seed,
          "gap_minus_mm":3.,"gap_plus_mm":3.,"anode_voltage_v":1800.,
          "anode_diameter_um":30.,"magnetic_field_t":0.,
          "n_u":40,"n_v":30,"du_mm":.1,"dv_mm":.2,
          "gas":"Ar/CO2 70:30, 293.15 K, 760 Torr",
        }
        for key, value in required.items():
            if cfg.get(key) != value:
                raise ValueError(f"Wrong immutable p4 LUT500 parameter {key} in {meta}")
    elif present:
        raise ValueError(f"Garfield columns present without frozen manifest: {folder}")
    print(f"{name}: {len(present)}/{expected_cols} complete columns; "
          f"{len(present)*NV*samples:,}/{per_config:,} raw trials")
    if mode == "verify":
        if not meta.is_file() or len(present) != expected_cols:
            raise ValueError(f"Not all columns completed for {name}")
        for iu in range(expected_cols):
            file = folder/"columns"/f"iu_{iu:03d}.csv"
            validate_column(file, iu, pitch, samples)
            if load_column_timing(file) is None:
                raise ValueError(f"Missing SHA-pinned timing for {file}")
        assembled = folder/"garfield_electrons.csv"
        if not assembled.is_file():
            raise ValueError(f"No assembled LUT: {assembled}")
        with assembled.open(newline="",encoding="utf-8") as f:
            total = sum(1 for _ in csv.DictReader(f))
        if total != per_config:
            raise ValueError(f"Wrong assembled LUT rows: {total}/{per_config}")
        summary = json.loads((folder/"summary.json").read_text())
        if (summary["n_electrons"] != per_config or
            summary["n_anode_collected"] +
            summary["n_field_wire_collected"] +
            summary["n_not_on_wires"] != per_config):
            raise ValueError("Endpoint conservation failure")
print(f"Total saved columns: {complete}/80; "
      f"trials: {complete*NV*samples:,}/1,200,000")
if mode == "verify":
    meta = root/"two_pitch_configurations.json"
    if not meta.is_file():
        raise ValueError("Missing final two-configuration summary")
    result = json.loads(meta.read_text(encoding="utf-8"))
    configs = {(r["pitch_mm"],r["field_wire_diameter_um"],
                r["n_electrons"]) for r in result["results"]}
    assert configs == {(4.,0.,per_config),(4.,75.,per_config)}
    print("PASS: both p4 LUT500 physics maps, 1.2M microscopic trials.")
elif mode == "check":
    print("PASS: p4 LUT500 manifest/checkpoint preflight; no Garfield run.")
PY
if [[ "$MODE" == status ]]; then
  if pid_running; then
    echo "RUNNING: launcher PID=$(cat "$PIDFILE")"
  else
    echo "NOT RUNNING (or completed)"
  fi
  if [[ -f "$LOG" ]]; then
    echo "--- Last 12 log lines ---"
    tail -n 12 "$LOG"
  fi
fi
