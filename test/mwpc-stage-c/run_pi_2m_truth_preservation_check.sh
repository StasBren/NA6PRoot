#!/usr/bin/env bash
# Mandatory READ-ONLY audit BEFORE the 6M generated-parent campaign.
# Never changes frozen Geant4 generator/simulator or original 100k.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd -P)"
REPO="$(cd "$HERE/../.." && pwd -P)"
SUITE="$REPO/test_runs/mwpc_stage_c/pi_geant4_2m_parent_campaign_v1"
EXPORT="$REPO/test_runs/mwpc_stage_c/pi_six_geometry_2m_native_exports_v1"
OUT="$REPO/test_runs/mwpc_stage_c/pi_weekend_2m_operations_v1"
mkdir -p "$OUT"
cd "$REPO"
NA6P_STAGEC_REQUIRE_ROOT_HADD=1
source "$HERE/activate_pi_stagec_env.sh"
python3 - "$SUITE" "$EXPORT" <<'PY'
import json,sys
from pathlib import Path
suite,export=map(Path,sys.argv[1:])
campaign=json.loads((suite/"campaign.json").read_text())["frozen_inputs"]
assert campaign["events_per_source"]==2000000
assert campaign["chunk_size"]==50000
assert campaign["channels"]==["Jpsi","Phi","Omega"]
for ch in ("Jpsi","Phi","Omega"):
 for i in (1,2):
  dir=suite/ch.lower()/f"chunk{i:04d}"
  completion=dir/"completed.json"
  if not completion.is_file():
   raise SystemExit(f"ERROR: missing 100k benchmark ROOT checkpoint: {completion}")
  m=json.loads(completion.read_text())
  required=["MCKine.root","HitsMuonSpec.root","geometry.root",
            "MuonReferencePlanes.root","na6pLayout.ini"]
  if m.get("events")!=50000 or m.get("chunk")!=i or \
      not m.get("original_worker_subruns_preserved") or \
      not m.get("verified_gas_step_tree") or \
      not m.get("verified_virtual_muon_plane_tree"):
    raise SystemExit(f"ERROR: invalid frozen native 50k parent checkpoint: {completion}")
  for name in required:
   p=dir/name
   if not p.is_file() or p.stat().st_size<=0 or name not in m["output_sha256"]:
    raise SystemExit(f"ERROR: cannot recover full physical ROOT/geometry: {p}")
  if not (dir/"raw_workers").is_dir():
   raise SystemExit(f"ERROR: original-worker track ID provenance NOT retained: {dir}")
  ex=export/ch.lower()/f"chunk{i:04d}"/"crossings_manifest.json"
  if not ex.exists():
   raise SystemExit(f"ERROR: missing validated native 100k CSV export: {ex}")
  em=json.loads(ex.read_text())
  if em["immutable_input"]["input_root_sha256"] != {
        name:m["output_sha256"][name] for name in
        ("HitsMuonSpec.root","geometry.root","MCKine.root")}:
   raise SystemExit(f"ERROR: ROOT provenance differs from native CSV: {ex}")
print("INVENTORY PASS: six 50k Geant4 checkpoints, event/kine/hit/geometry/"
      "worker truth and native CSV ROOT hashes agree.",flush=True)
PY
echo "ARCHIVE: copy actual frozen installed field maps, cuts, source and binaries"
python3 "$HERE/freeze_pi_2m_transport_inputs.py"
for source in jpsi phi omega; do
  for chunk in 1 2; do
    dir="$SUITE/$source/$(printf 'chunk%04d' "$chunk")"
    output="$OUT/truth_audit_$source.chunk$chunk.log"
    invocation="$HERE/audit_pi_muon_truth_bank.C+(\"$dir\",50000)"
    echo "AUDIT original ROOT: $source chunk$chunk"
    root -l -b -q "$invocation" > "$output" 2>&1 || {
      tail -n 60 "$output" >&2
      echo "STOP: native truth audit failed; do not start 6M Geant4" >&2
      exit 4
    }
    if ! grep -q '^MUON_TRUTH_BANK_AUDIT_PASS$' "$output"; then
      tail -n 60 "$output" >&2
      echo "STOP: ROOT did not confirm track/hit provenance" >&2
      exit 5
    fi
    grep -E '  checked|MUON_TRUTH_BANK_AUDIT_PASS' "$output"
  done
done
echo "TRUTH BANK PASS: original Geant4 ROOT preserves event entries, "
echo "per-muon mother/track kinematics, native gas entry/exit, "
echo "positive-Edep gas steps, interpolated reference planes, original "
echo "worker-local truth. No simulation performed by this audit."
