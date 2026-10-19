#!/usr/bin/env bash
# Physically NEW symmetric MWPC field: anode at v=0 between cathodes
# v=-3mm and v=+3mm. Default 50 Garfield drifts/start box (30k total).
#
# This is a NEW microscopic GARFIELD run, NOT a rerun of the 600k Geant4
# campaign. Results remain completely separate from legacy 2+4mm maps.
#
# Usage:
#   bash test/mwpc-stage-c/run_pi_centered_transport_box_map.sh [50] [180019]
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd -P)"
REPO="$(cd "$HERE/../.." && pwd -P)"
cd "$REPO"
SAMPLES="${1:-50}"
BASE_SEED="${2:-180019}"
if ! [[ "$SAMPLES" =~ ^[1-9][0-9]*$ && "$BASE_SEED" =~ ^[0-9]+$ ]]; then
  echo "Usage: $0 [positive-electrons-per-box] [nonnegative-seed]" >&2
  exit 2
fi
ROOT_OUT="$REPO/test_runs/mwpc_stage_c/c_transport_box_map_gm3p0_gp3p0"
mkdir -p "$ROOT_OUT"

# Existing run resumes only with exactly the same immutable parameters,
# including sample count and seed, to prevent mixing incompatible columns.
python3 - "$ROOT_OUT/scan_progress.log" "$SAMPLES" "$BASE_SEED" <<'PY'
from pathlib import Path
import re
import sys
log=Path(sys.argv[1])
if not log.exists():
    raise SystemExit(0)
content=log.read_text(errors="replace")
def previous(key):
    return re.findall(re.escape(key)+r"\s*:\s*([^\n]+)",content)
checks=[
    ("electrons / box",str(int(sys.argv[2]))),
    ("gap (-,+) [mm]","3.0, 3.0"),
    ("B [T]","0"),
    ("du,dv [mm]","0.2, 0.2"),
]
for key,wanted in checks:
    actual=previous(key)
    if actual and any(v.strip()!=wanted for v in actual):
        raise SystemExit(
           f"Existing centered response was produced under different "
           f"{key}: {actual}; use a new output directory, do not mix maps"
        )
# The underlying Garfield runner currently does not print its base seed.
# Stamp a separate immutable seed token in the wrapper output directory.
seal=log.parent/"pi_centered_run_config.json"
if not seal.exists() and any((log.parent/"columns").glob("iu_*.csv")):
    raise SystemExit("Cannot resume prior 3+3 Garfield columns with unverified base seed: "
                     "move aside this directory for provenance, do not mix")
if seal.exists():
    import json
    val=json.loads(seal.read_text())
    if (val.get("electrons_per_box")!=int(sys.argv[2]) or
        val.get("base_seed")!=int(sys.argv[3])):
        raise SystemExit("Incompatible centered map sample count or seed")
PY
python3 - "$ROOT_OUT/pi_centered_run_config.json" "$SAMPLES" "$BASE_SEED" <<'PY'
from pathlib import Path
import json
import sys
p=Path(sys.argv[1])
cfg={
 "geometry":"centered_3+3_v0_anode",
 "electrons_per_box":int(sys.argv[2]),
 "base_seed":int(sys.argv[3]),
 "B_tesla":0,
 "gap_minus_mm":3.0,
 "gap_plus_mm":3.0,
 "du_mm":0.2,
 "dv_mm":0.2,
 "wire_pitch_mm":4.0,
 "anode_voltage_V":1800,
}
if not p.exists():
    p.write_text(json.dumps(cfg,indent=2)+"\n")
else:
    old=json.loads(p.read_text())
    if old!=cfg:
        raise SystemExit("Centered-map generation settings changed: refusing resume")
PY

export C_TRANSPORT_GAP_MINUS_MM=3.0
export C_TRANSPORT_GAP_PLUS_MM=3.0
export C_TRANSPORT_B_TESLA=0
export C_TRANSPORT_DU_MM=0.2
export C_TRANSPORT_DV_MM=0.2
python3 -m py_compile \
 "$HERE/pi_validate_centered_box_map.py" \
 "$HERE/test_pi_validate_centered_box_map.py"
python3 "$HERE/test_pi_validate_centered_box_map.py"

echo "RUNNING new microscopic Garfield drift for centered 3+3 mm field."
echo "Frozen Geant4 transport and legacy 2+4 mm files will not change."
bash "$HERE/run_c_transport_box_map.sh" "$SAMPLES" "$BASE_SEED"

python3 "$HERE/pi_validate_centered_box_map.py" \
  --input "$ROOT_OUT/c_transport_box_electrons.csv" \
  --per-box "$SAMPLES" \
  --progress "$ROOT_OUT/scan_progress.log" \
  --manifest "$ROOT_OUT/pi_centered_transport_manifest.json"
echo "PASS: centered 3+3 mm Garfield electron response field validated."
