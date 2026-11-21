#!/usr/bin/env bash
# Read-only frozen muon0301 comparison. Three NEW deterministic RKF banks.
# No new Heed, microscopic avalanche, FEM, U/V, or electrons.
set -eo pipefail

ROOT_INSTALL="$HOME/na6p/install/root"
GARFIELD_INSTALL="$HOME/na6p/install/garfieldpp"
source "$ROOT_INSTALL/bin/thisroot.sh"
export LD_LIBRARY_PATH="$GARFIELD_INSTALL/lib:$ROOT_INSTALL/lib:$LD_LIBRARY_PATH"

BIN="$PWD/test/mwpc-avalanche/build-polar/mwpc_deterministic_polar_plane_bank"
PY="test/mwpc-stage-c/cathode-ion-reference/replay_deterministic_polar_bank.py"
MOBILITY="$HOME/na6p/src/garfieldpp/Data/IonMobility_Ar+_Ar.txt"
PILOT="$PWD/test_runs/mwpc_stage_c/polar_grid_2100_v1"
SOURCE="$HOME/na6p/cathode-overnight-2100/test_runs/mwpc_stage_c/overnight_2100_weighted_k256_v1/muon_0301/ion_births.root"
BASE="$HOME/na6p/cathode-independent-stereo-2100/test_runs/mwpc_stage_c/independent_muon_dt10_v1"
REF="$BASE/muon_0301_direct_dt10_waveforms.csv"
SUMMARY="$BASE/muon_0301_direct_dt10_summary.csv"
IDW="$BASE/muon_0301_fast_dt10_waveforms.csv"

for file in "$BIN" "$PY" "$MOBILITY" "$PILOT/grid6x12_dt10.root" "$SOURCE" \
            "$REF" "$SUMMARY" "$IDW"; do
  if [[ ! -s "$file" ]]; then
    echo "STOP: missing $file" >&2
    exit 1
  fi
done
if [[ ! -x "$BIN" ]]; then
  echo "STOP: bank generator not executable" >&2
  exit 1
fi

(cd "$(dirname "$SOURCE")" && test -s complete.txt && sha256sum -c SHA256SUMS) || {
  echo "STOP: original microscopic source integrity check failed" >&2
  exit 1
}

# Preflight output paths: never overwrite or silently reuse a partial run.
for spec in "6x24:6:24" "12x12:12:12" "12x24:12:24"; do
  IFS=: read -r shape nr nt <<< "$spec"
  for p in "$PILOT/grid"$shape"_dt10.root" "$PILOT/grid"$shape"_generation.log"; do
    if [[ -e "$p" ]]; then
      echo "STOP: existing bank/log $p" >&2; exit 1
    fi
  done
done
for shape in 6x12 6x24 12x12 12x24; do
  if [[ -e "$PILOT/replay_muon0301_"$shape"_windowed_v2" ]]; then
    echo "STOP: existing replay output for $shape" >&2; exit 1
  fi
done

python -m py_compile "$PY"

# Resolve angular vs radial grid resolution independently.
for spec in "6x24:6:24" "12x12:12:12" "12x24:12:24"; do
  IFS=: read -r shape nr nt <<< "$spec"
  BANK="$PILOT/grid"$shape"_dt10.root"
  LOG="$PILOT/grid"$shape"_generation.log"
  echo "=== C++ RKF grid $shape ==="
  "$BIN" --mobility "$MOBILITY" --out-root "$BANK" \
    --n-radius "$nr" --n-angle "$nt" \
    --clearance-min-um 0.00001 --radius-max-um 80 \
    --theta-min-deg 0.001 --dt-ns 10 --window-us 100 --max-step-mm 0.05 \
    > "$LOG" 2>&1 || {
      echo "STOP: failed $shape, last 20 lines of RKF log:" >&2
      tail -n 20 "$LOG" >&2
      exit 1
    }
  tail -n 3 "$LOG"
done

# Re-run original 6x12 prediction in a NEW directory to obtain the same
# per-window diagnostic as the new grids. No original result is overwritten.
for shape in 6x12 6x24 12x12 12x24; do
  BANK="$PILOT/grid"$shape"_dt10.root"
  OUT="$PILOT/replay_muon0301_"$shape"_windowed_v2"
  echo "=== ORIGINAL MUON 0301: grid $shape vs unchanged direct RKF ==="
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
    python "$PY" --grid-root "$BANK" --ion-file "$SOURCE" \
      --reference-waveforms "$REF" --reference-summary "$SUMMARY" \
      --baseline-waveforms "$IDW" --out-dir "$OUT" || {
        echo "STOP: replay failed for $shape; retained outputs" >&2
        exit 1
      }
done

echo "=== Grid convergence summary ==="
python - "$PILOT" <<'PY'
import csv
from pathlib import Path
import sys

folder=Path(sys.argv[1])
results=[]
for name in ("6x12", "6x24", "12x12", "12x24"):
    d=folder/("replay_muon0301_"+name+"_windowed_v2")
    with (d/"event_metrics.csv").open() as f:
        bymodel={r["model"]:r for r in csv.DictReader(f)}
    p=bymodel["polar_grid_bilinear_uncompressed"]
    b=bymodel["EXISTING_IDW4"]
    with (d/"window_metrics.csv").open() as f:
        bywindow={r["window"]:r for r in csv.DictReader(f)
                  if r["model"]=="polar_grid_bilinear_uncompressed"}
    results.append(dict(
        grid=name,
        joint_L2=float(p["combined_L2"]),
        IDW4_L2=float(b["combined_L2"]),
        upper_100us_error_fC=float(p["plus_window_error_fC"]),
        lower_100us_error_fC=float(p["minus_window_error_fC"]),
        first100ns_L2=float(bywindow["first_100ns"]["combined_L2_window"]),
        from100ns_to1us_L2=float(bywindow["100ns_to_1us"]["combined_L2_window"]),
        from1us_to10us_L2=float(bywindow["1us_to_10us"]["combined_L2_window"]),
        from10us_to100us_L2=float(bywindow["10us_to_100us"]["combined_L2_window"])
    ))
out=folder/"polar_convergence_summary_v2.csv"
with out.open("x",newline="") as f:
    writer=csv.DictWriter(f,fieldnames=list(results[0]))
    writer.writeheader();writer.writerows(results)
for r in results:
    print(r)
print("PASS original muon 0301 polar convergence:", out)
PY

echo "No run of locked independent muons 0408 or 0446."
