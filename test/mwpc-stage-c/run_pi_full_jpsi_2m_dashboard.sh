#!/usr/bin/env bash
# Read-only COMPLETE J/psi 2M-parent six-station dashboard.
# Does NOT start/stop/modify the original full 6M fold.
# Usage: bash test/mwpc-stage-c/run_pi_full_jpsi_2m_dashboard.sh check|start|status|verify
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd -P)"
REPO="$(cd "$HERE/../.." && pwd -P)"
BASE="$REPO/test_runs/mwpc_stage_c/pi_full_jpsi_2m_dashboard_v1"
SCRIPT="$HERE/pi_full_jpsi_2m_dashboard.py"
PIDFILE="$BASE/dashboard.pid"
LOG="$BASE/dashboard_build.log"
MODE="$1"
if [[ "$#" -ne 1 || ! "$MODE" =~ ^(check|start|status|verify)$ ]]; then
  echo "Usage: bash $0 check|start|status|verify" >&2
  exit 2
fi
mkdir -p "$BASE"
active() {
  [[ -f "$PIDFILE" ]] || return 1
  local pid
  pid="$(cat "$PIDFILE")"
  [[ "$pid" =~ ^[0-9]+$ ]] || return 1
  kill -0 "$pid" 2>/dev/null || return 1
  ps -p "$pid" -o args= | grep -q 'pi_full_jpsi_2m_dashboard.py build'
}
if [[ "$MODE" == status ]]; then
  completed=0
  for st in 0 1 2 3 4 5; do
    if [[ -f "$BASE/station_summaries/MS$st.json" ]];then
      echo "AVAILABLE: full 2M J/psi station MS$st"
      completed=$((completed+1))
    else
      echo "PENDING: MS$st"
    fi
  done
  echo "Station dashboards: $completed/6 (40 completed 50k-parent chunks each)"
  for st in 0 1 2 3 4 5; do
    if [[ -f "$BASE/pages/jpsi40-bore/MS${st}_phi_atlas.html" ]];then
      echo "PHI ATLAS READY: MS$st (all physical DET modules)"
    fi
  done
  if active; then echo "RUNNING builder PID=$(cat "$PIDFILE")"
  else echo "No active dashboard builder detected";fi
  if [[ -f "$BASE/complete_dashboard_manifest.json" ]];then
    echo "FINAL MANIFEST: created (run verify for SHA validation)"
  fi
  if [[ -f "$LOG" ]];then
    echo "---- Last 12 builder log lines ----"
    tail -n 12 "$LOG"
  fi
  exit 0
fi

python3 -m py_compile "$SCRIPT" "$HERE/build_pi_100k_english_dashboard.py" \
  "$HERE/pi_full_jpsi_2m_phi_atlas.py" \
  "$HERE/test_pi_full_jpsi_2m_phi_atlas.py"
python3 "$HERE/test_pi_full_jpsi_2m_phi_atlas.py"
python3 "$HERE/test_pi_100k_english_dashboard.py"
if [[ "$MODE" == check ]];then
  python3 "$SCRIPT" check
  exit
fi
if [[ "$MODE" == verify ]];then
  python3 "$SCRIPT" verify
  exit
fi
if active;then
  echo "Builder already running PID=$(cat "$PIDFILE")"
  exit 0
fi
python3 "$SCRIPT" check
echo "BUILD: read-only 240 SHA-verified Jpsi/6-station fold checkpoints."
echo "Output: $BASE/pages/index.html (includes all-chamber phi atlas)"
echo "Original 6M Geant4 + Heed/Garfield continues uninterrupted."
nohup nice -n 15 python3 "$SCRIPT" build \
  > "$LOG" 2>&1 < /dev/null &
pid=$!
echo "$pid" > "$PIDFILE"
echo "STARTED dashboard PID=$pid"
echo "Monitor: bash '$HERE/run_pi_full_jpsi_2m_dashboard.sh' status"
echo "Live log: tail -f '$LOG'"
