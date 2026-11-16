#!/usr/bin/env bash
# Isolated and resumable P3 noFW 30/50um anode diameter and voltage scan
set -Eeuo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"
cd "$REPO"

MODE="${1:-pilot}"
if [[ "$MODE" != pilot && "$MODE" != full ]];then
  echo "Usage: bash $0 [pilot|full]" >&2; exit 2
fi
CAP="${HV_DIAM_CAP:-1000000}"
OUT="${HV_DIAM_OUT:-$REPO/test_runs/mwpc_stage_c/hvdiam_study_20261116_v1}"
INPUT="$REPO/test/mwpc-stage-c/p3_wiretruth_1000_selection.csv"
BUILD="$REPO/test/mwpc-avalanche/build-hvdiam"
EXE="$BUILD/mwpc_p3_noFW_penning057_hvdiam_v1"

if [[ ! "$CAP" =~ ^[0-9]+$ ]] || (( CAP < 200000 || CAP > 5000000 ));then
  echo "HV_DIAM_CAP must be between 200000 and 5000000" >&2; exit 2
fi
[[ -s "$INPUT" ]] || { echo "Missing frozen crossing input $INPUT" >&2; exit 2; }
[[ -f "$HOME/na6p/env.sh" ]] && source "$HOME/na6p/env.sh"
SETUP="$HOME/na6p/install/garfieldpp/share/Garfield/setupGarfield.sh"
[[ -f "$SETUP" ]] && source "$SETUP"
command -v root >/dev/null || { echo "ROOT not found" >&2; exit 2; }

cmake -S "$REPO/test/mwpc-avalanche" -B "$BUILD" \
  -DGarfield_DIR="$HOME/na6p/install/garfieldpp/lib/cmake/Garfield"
cmake --build "$BUILD" --target mwpc_p3_noFW_penning057_hvdiam_v1 -j2

mkdir -p "$OUT"
INPUT_SHA="$(sha256sum "$INPUT" | cut -d ' ' -f1)"
BIN_SHA="$(sha256sum "$EXE" | cut -d ' ' -f1)"
LOCK="input=$INPUT_SHA
binary=$BIN_SHA
cap=$CAP
chunk_size=5
gas=ArCO2_70_30_Penning_r0.57_B0
geometry=P3_noFW_cathodes_pm3mm"
if [[ -e "$OUT/config.lock" ]];then
  [[ "$(cat "$OUT/config.lock")" == "$LOCK" ]] || {
    echo "Configuration changed; choose NEW HV_DIAM_OUT" >&2; exit 3;
  }
else
  printf '%s\n' "$LOCK" > "$OUT/config.lock"
fi

run_chunk() {
  local v="$1" diameter="$2" start="$3"
  local end=$((start+5))
  local group="$OUT/V_${v}_D${diameter}" idx=$((start/5))
  local name dir
  printf -v name 'chunk_%04d' "$idx"
  dir="$group/$name"
  mkdir -p "$group"
  if [[ -e "$dir/complete.txt" ]];then
    (cd "$dir" && sha256sum -c SHA256SUMS >/dev/null) || {
      echo "Corrupt completed chunk $dir" >&2;exit 3;
    }
    echo "SKIP verified $dir";return 0
  fi
  if [[ -d "$dir" ]];then
    mv "$dir" "$dir.incomplete_$(date +%Y%m%dT%H%M%S)"
  fi
  mkdir -p "$dir"
  {
    head -n1 "$INPUT"
    sed -n "$((start+2)),$((end+1))p" "$INPUT"
  } > "$dir/input.csv"
  [[ "$(wc -l < "$dir/input.csv")" -eq 6 ]] || { echo "Incomplete input" >&2;exit 3; }
  echo "RUN v=$v diameter=$diameter um muons $start..$((end-1))"
  if ! "$EXE" --wiretruth-hv-diam-p3-noFW "$dir/input.csv" \
     "$dir/micro" "$v" "$CAP" "$diameter" > "$dir/run.log" 2>&1;then
    tail -n 40 "$dir/run.log" >&2;exit 4
  fi
  [[ -s "$dir/micro_crossing.csv" && -s "$dir/micro_seed.csv" &&
     -s "$dir/micro_wire_muon.csv" ]] || {
     echo "Missing Garfield outputs $dir" >&2;exit 4;
  }
  [[ "$(wc -l < "$dir/micro_crossing.csv")" -eq 6 ]] || {
     echo "Wrong muon count $dir" >&2;exit 4;
  }
  if ! awk -F, '
    NR==1 {for(i=1;i<=NF;i++){if($i=="capped") c=i;if($i=="successful") s=i;}
           if(!c||!s) bad=1;next}
    NR>1 {if(!c||!s||($c+0)!=0||($s+0)!=1) bad=1}
    END{if(bad||NR<2)exit 1}
  ' "$dir/micro_seed.csv";then
    echo "CAPPED or failed avalanche found. STOP: $dir" >&2;exit 5
  fi
  (cd "$dir" && sha256sum input.csv run.log micro_*.csv > SHA256SUMS)
  printf 'HV=%s diameter_um=%s muons=5 cap=%s\n' "$v" "$diameter" "$CAP" > "$dir/complete.txt"
  echo "PASS $dir"
}

# Same 5 primary muons in every group (paired comparison).
run_chunk 2000 30 0
run_chunk 2000 50 0
run_chunk 2100 30 0
run_chunk 2050 30 0
# Full mode extends each 30um group to 20, but keeps 50um at 5.
if [[ "$MODE" == full ]];then
  for v in 2000 2050 2100;do
    for i in 5 10 15;do run_chunk "$v" 30 "$i";done
  done
fi

for groupid in 2000:30 2050:30 2100:30 2000:50;do
  v="${groupid%%:*}"
  diameter="${groupid##*:}"
  group="$OUT/V_${v}_D${diameter}"
  n=5
  [[ "$MODE" == full && "$diameter" == 30 ]] && n=20
  for kind in crossing seed;do
    merged="$group/merged_$kind.csv"
    head -n1 "$group/chunk_0000/micro_$kind.csv" > "$merged"
    for ((j=0;j<n/5;j++));do
      printf -v chunk 'chunk_%04d' "$j"
      [[ -f "$group/$chunk/complete.txt" ]] || {
        echo "Chunk unsealed: $group/$chunk" >&2;exit 4;
      }
      tail -n +2 "$group/$chunk/micro_$kind.csv" >> "$merged"
    done
  done
done

root -l -b -q "$HERE/Analyze_HVdiam_ROOT.C(\"$OUT\")"
echo "DONE: $OUT/analysis"
