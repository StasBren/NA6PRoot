#!/usr/bin/env bash
# 6M generated-parent MWPC production. Resumable native Geant4,
# verified Geant4 ROOT->CSV, common Heed + 6 frozen Garfield-LUT folding.
# Usage: bash test/mwpc-stage-c/run_pi_weekend_2m.sh check|start|status
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd -P)"
REPO="$(cd "$HERE/../.." && pwd -P)"
SUITE="$REPO/test_runs/mwpc_stage_c/pi_geant4_2m_parent_campaign_v1"
EXPORT="$REPO/test_runs/mwpc_stage_c/pi_six_geometry_2m_native_exports_v1"
FOLD="$REPO/test_runs/mwpc_stage_c/pi_six_geometry_2m_timing_v1/full"
OPER="$REPO/test_runs/mwpc_stage_c/pi_weekend_2m_operations_v1"
PIDFILE="$OPER/campaign.pid"
LOG="$OPER/campaign.log"
MODE="${1:-check}"
[[ "$#" -eq 1 && "$MODE" =~ ^(check|start|status|work)$ ]] || {
  echo "Usage: bash $0 check|start|status" >&2; exit 2;
}
cd "$REPO"
preflight(){
python3 - "$SUITE" <<'PY'
import json,sys,shutil
from pathlib import Path
suite=Path(sys.argv[1])
manifest=suite/"campaign.json"
if not manifest.is_file():
    raise SystemExit("ERROR: existing frozen Geant4 2m/source manifest missing: "+str(manifest))
inputs=json.loads(manifest.read_text())["frozen_inputs"]
for k,v in (("mode","full"),("events_per_source",2_000_000),
            ("chunk_size",50_000),("channels",["Jpsi","Phi","Omega"])):
    if inputs.get(k)!=v:
        raise SystemExit(f"ERROR: frozen {k}: {inputs.get(k)!r}, expected {v!r}")
workers=inputs.get("workers")
if not isinstance(workers,int) or not 1<=workers<=32:
    raise SystemExit("ERROR: workers missing from frozen Geant4 inputs")
complete=[]
bytes_used=0
for ch in ("jpsi","phi","omega"):
    for i in range(1,41):
        path=suite/ch/f"chunk{i:04d}"
        cp=path/"completed.json"
        if cp.exists():
            m=json.loads(cp.read_text())
            if m.get("events")!=50_000 or m.get("chunk")!=i:
                raise SystemExit("ERROR: invalid completed Geant4 marker "+str(cp))
            complete.append(cp)
            bytes_used+=sum(p.stat().st_size for p in path.rglob("*") if p.is_file())
free=shutil.disk_usage(suite).free/(1024**3)
print(f"FROZEN: 3 sources x 2,000,000 parents; 40 x 50k chunks/source; "
      f"{workers} workers; {len(complete)}/120 completed Geant4 chunks")
if complete:
    remain=(120-len(complete)) * bytes_used / len(complete) / (1024**3)
    print(f"Existing Geant4 chunks: {bytes_used/(1024**3):.1f} GiB; "
          f"estimated remaining at measured average: ~{remain:.1f} GiB.")
    if free < remain+35:
        print("WARNING: estimated remaining Geant4 storage exceeds free space "
              "minus 35 GiB guard; producer will pause at safe checkpoint.")
print(f"Linux WSL filesystem free: {free:.1f} GiB; minimum per-chunk reserve 35 GiB.")
print("Also check physical WINDOWS host-drive free space and disable sleep.")
if free<40:
    raise SystemExit("ERROR: need at least 40 GiB free to start safely.")
PY
}
status(){
if [[ -f "$PIDFILE" ]] && kill -0 "$(cat "$PIDFILE")" 2>/dev/null; then
  echo "RUNNING launcher PID=$(cat "$PIDFILE")"
else
  echo "No active launcher PID."
fi
python3 - "$SUITE" "$EXPORT" "$FOLD" <<'PY'
import sys
from pathlib import Path
suite,export,fold=map(Path,sys.argv[1:])
for label,root,glob,limit in (
 ("Geant4",suite,"*/chunk*/completed.json",120),
 ("ROOT native CSV",export,"*/chunk*/crossings_manifest.json",120),
 ("Six-LUT station folds",fold,"*/MS*/chunk*.manifest.json",720)):
 print(f"{label}: {sum(1 for _ in root.glob(glob))}/{limit}")
print("Geant4 campaign_complete.json:",
      (suite/"campaign_complete.json").is_file())
PY
if [[ -f "$LOG" ]];then echo "--- Last 25 production log lines ---";tail -n 25 "$LOG";fi
}
if [[ "$MODE" == status ]];then status;exit 0;fi
preflight
if [[ "$MODE" == check ]];then
  echo "PREFLIGHT PASS (no physics executed).";exit 0
fi
mkdir -p "$OPER"
if [[ "$MODE" == start ]];then
  if [[ -f "$PIDFILE" ]] && kill -0 "$(cat "$PIDFILE")" 2>/dev/null; then
    echo "ERROR: existing active runner PID=$(cat "$PIDFILE")";exit 4
  fi
  nohup bash "$HERE/run_pi_weekend_2m.sh" work > "$LOG" 2>&1 < /dev/null &
  child=$!
  echo "$child" > "$PIDFILE"
  echo "DETACHED: production launcher PID=$child"
  echo "Monitor: bash '$HERE/run_pi_weekend_2m.sh' status"
  echo "Live log: tail -f '$LOG'"
  echo "Windows/WSL must stay awake; every phase saves checkpoints."
  exit 0
fi
echo "START Geant4 6M production $(date -Is)"
NA6P_STAGEC_REQUIRE_ROOT_HADD=1
source "$HERE/activate_pi_stagec_env.sh"
WORKERS="$(python3 - "$SUITE/campaign.json" <<'PY'
import sys,json
print(json.load(open(sys.argv[1]))["frozen_inputs"]["workers"])
PY
)"
echo "PHASE 1/3: Native Geant4; preserve original 298-detector geometry."
python3 "$HERE/run_pi_geant4_campaign.py" \
  --mode full --events-per-source 2000000 --chunk-size 50000 \
  --workers "$WORKERS" --hours 0 --min-free-gib 35 --suite "$SUITE"
if [[ ! -f "$SUITE/campaign_complete.json" ]]; then
  echo "SAFE PAUSE: Geant4 checkpoints preserved, perhaps low disk. "
  echo "No new exports/folds started. Extend disk, then rerun start."
  exit 5
fi
echo "PHASE 2/3: ROOT exports; source+chunk+event identity retained."
python3 "$HERE/pi_2m_export_completed_chunks.py" \
  --root "$(command -v root)" --max-chunks-per-source 40
echo "PHASE 3/3: Matched six frozen Garfield LUTs and common Heed births."
python3 "$HERE/pi_2m_native_six_geometry_fold.py" fold \
  --channel all --station -1 --chunks 40
echo "PRODUCTION PASS: 120 Geant4 chunks, 120 native exports, "
echo "720 station/chunk six-Garfield folds. Strict 40-chunk final merge "
echo "and parent-bootstrap to be handled separately; 100k reports preserved."
