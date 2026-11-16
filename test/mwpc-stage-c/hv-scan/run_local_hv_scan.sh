#!/usr/bin/env bash
# Reproducible LOCAL laptop WSL voltage scan, ROOT-only physics analysis.
# No Python. Uses same 1000 frozen Geant4 crossings and same Heed seeds at every HV.
# First call: bash run_local_hv_scan.sh 5; extend to 20, 100, 1000 afterwards.
set -Eeuo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"
cd "$REPO"

N="${1:-5}"
CHUNK="${HV_CHUNK_SIZE:-5}"
CAP="${HV_AVALANCHE_CAP:-1000000}"
VOLTAGES=(1900 1950 2000 2050)
SELECTION="$REPO/test/mwpc-stage-c/p3_wiretruth_1000_selection.csv"
OUT="${HV_OUTPUT_DIR:-$REPO/test_runs/mwpc_stage_c/p3_hv_scan_v1}"
BUILD="$REPO/test/mwpc-avalanche/build-hvscan"
EXE="$BUILD/mwpc_p3_noFW_penning057_hvscan_v1"

if [[ ! "$N" =~ ^[0-9]+$ ]] || (( N < 1 || N > 1000 )); then
  echo "Choose 1..1000 muons, e.g. 5, 20, 100 or 1000" >&2; exit 2
fi
if [[ ! "$CHUNK" =~ ^[0-9]+$ ]] || (( CHUNK<1 || CHUNK>20 )); then
  echo "HV_CHUNK_SIZE must be in 1..20" >&2; exit 2
fi
if (( N % CHUNK != 0 )); then
  echo "For resume consistency use N divisible by HV_CHUNK_SIZE=$CHUNK" >&2;exit 2
fi
if [[ ! "$CAP" =~ ^[0-9]+$ ]] || (( CAP<200000 || CAP>5000000 )); then
  echo "HV_AVALANCHE_CAP must be in [200000,5000000]" >&2;exit 2
fi

if [[ -f "$HOME/na6p/env.sh" ]]; then source "$HOME/na6p/env.sh";fi
GARFIELD_SETUP="$HOME/na6p/install/garfieldpp/share/Garfield/setupGarfield.sh"
if [[ -f "$GARFIELD_SETUP" ]];then source "$GARFIELD_SETUP";fi
command -v root >/dev/null || { echo "ROOT is required for the final analysis." >&2;exit 2; }
command -v cmake >/dev/null || { echo "cmake is required." >&2;exit 2; }
[[ -s "$SELECTION" ]] || { echo "Frozen 1000 selection is missing: $SELECTION" >&2;exit 2; }

echo "BUILDING separate high-voltage target ..."
cmake -S "$REPO/test/mwpc-avalanche" -B "$BUILD" \
  -DGarfield_DIR="$HOME/na6p/install/garfieldpp/lib/cmake/Garfield"
cmake --build "$BUILD" --target mwpc_p3_noFW_penning057_hvscan_v1 -j2

mkdir -p "$OUT"
SOURCE_SHA="$(sha256sum "$SELECTION" | cut -d' ' -f1)"
BINARY_SHA="$(sha256sum "$EXE" | cut -d' ' -f1)"
GIT_REV="$(git rev-parse HEAD)"
LOCK_CONTENT="input_sha256=$SOURCE_SHA
binary_sha256=$BINARY_SHA
git_revision=$GIT_REV
chunk_size=$CHUNK
avalanche_cap=$CAP
voltages=1900,1950,2000,2050
geometry=P3_noFW;ArCO2_70_30;Penning_r0.57;B0"
if [[ -e "$OUT/config.lock" ]];then
  [[ "$(cat "$OUT/config.lock")" == "$LOCK_CONTENT" ]] || {
    echo "Configuration differs from sealed campaign. Change output directory or preserve original settings." >&2;exit 3;
  }
else
  printf '%s\n' "$LOCK_CONTENT" > "$OUT/config.lock"
fi
echo "HV scan: muons/point=$N, chunk=$CHUNK, avalanche cap=$CAP"
echo "Output: $OUT"

run_chunk() {
  local V="$1" idx="$2" start="$3" finish="$4"
  local tag dir tmp nrows
  printf -v tag 'chunk_%04d' "$idx"
  dir="$OUT/HV_${V}/$tag"
  nrows=$((finish-start))
  if [[ -e "$dir/complete.txt" ]];then
    (cd "$dir" && sha256sum -c SHA256SUMS >/dev/null) || {
      echo "Previously completed chunk has checksum mismatch: $dir" >&2;exit 3;
    }
    echo "SKIP verified: HV=$V $tag"; return 0
  fi
  if [[ -d "$dir" ]];then
    local quarantine="${dir}.incomplete.$(date +%Y%m%d_%H%M%S)"
    mv "$dir" "$quarantine"
    echo "QUARANTINED incomplete chunk: $quarantine"
  fi
  mkdir -p "$dir"
  tmp="$dir/input.csv"
  {
    head -n 1 "$SELECTION"
    sed -n "$((start+2)),$((finish+1))p" "$SELECTION"
  } > "$tmp"
  [[ "$(wc -l < "$tmp")" -eq "$((nrows+1))" ]] || { echo "Input split error" >&2;exit 3; }
  echo "RUN HV=$V $tag $nrows muons (source rows $start..$((finish-1)))"
  if ! "$EXE" --wiretruth-hv-p3-noFW "$tmp" "$dir/microscopic" "$V" "$CAP" > "$dir/physics.log" 2>&1; then
    echo "Garfield++ ERROR: $dir/physics.log" >&2;tail -n 30 "$dir/physics.log" >&2;exit 4;
  fi
  local cross="$dir/microscopic_crossing.csv"
  local seeds="$dir/microscopic_seed.csv"
  [[ -s "$cross" && -s "$seeds" && -s "$dir/microscopic_wire_muon.csv" ]] || {
    echo "Missing Garfield outputs: $dir" >&2;exit 4;
  }
  [[ "$(wc -l < "$cross")" -eq "$((nrows+1))" ]] || {
    echo "Wrong crossing count: $cross" >&2;exit 4;
  }
  # Reject any size-capped avalanche; never quote biased gain fluctuation.
  if ! awk -F, '
    NR==1 { for(i=1;i<=NF;i++) if($i=="capped") cap=i; if(cap==0) exit 2 }
    NR>1 { if($(cap)+0!=0) bad++ }
    END { if(NR<2||bad>0) exit 1 }
  ' "$seeds"; then
    echo "CAPPED AVALANCHE detected at HV=$V. Inspect $seeds and physics.log." >&2
    echo "Do not use this voltage in RMS/mean plots; increase cap cautiously or stop." >&2
    exit 5
  fi
  (cd "$dir" && sha256sum input.csv physics.log microscopic_*.csv > SHA256SUMS)
  printf 'HV=%s muons=%s avalanche_cap=%s\n' "$V" "$nrows" "$CAP" > "$dir/complete.txt"
  echo "PASS HV=$V $tag"
}

for V in "${VOLTAGES[@]}";do
  mkdir -p "$OUT/HV_$V"
  chunk_id=0
  for ((start=0;start<N;start+=CHUNK));do
    run_chunk "$V" "$chunk_id" "$start" "$((start+CHUNK))"
    chunk_id=$((chunk_id+1))
  done
done

# Only reached if ALL voltage chunks were successfully sealed.
for V in "${VOLTAGES[@]}";do
  vdir="$OUT/HV_$V"
  for kind in crossing seed;do
    merged="$vdir/merged_${kind}.csv"
    first=1
    : > "$merged"
    for ((idx=0;idx<N/CHUNK;idx++));do
      printf -v name 'chunk_%04d' "$idx"
      file="$vdir/$name/microscopic_${kind}.csv"
      [[ -f "$vdir/$name/complete.txt" ]] || { echo "Unsealed chunk: $name" >&2;exit 4; }
      if ((first));then head -n1 "$file" >> "$merged";first=0;fi
      tail -n +2 "$file" >> "$merged"
    done
  done
done

echo "RUN ROOT/C++ bootstrap analysis and ATLAS-style plots"
root -l -b -q "$HERE/AnalyzeHV_ROOT.C(\"$OUT\")"
echo "DONE. Results: $OUT/analysis/"
