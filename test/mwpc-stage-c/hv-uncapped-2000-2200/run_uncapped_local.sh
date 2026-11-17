#!/usr/bin/env bash
# True Garfield++ AvalancheMicroscopic no-size-cutoff study.
# One real Geant4 muon PER PROCESS to contain WSL memory failure.
# Failure/timeout/oom never counts as a finite-gain observation.
set -Eeuo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"
cd "$REPO"

N="${1:-1}"
VOLTAGES=(2200 2150 2100 2050 2000) # expensive end first
TIMEOUT_S="${HV_UNCAPPED_TIMEOUT_S:-1200}"
VMEM_GB="${HV_UNCAPPED_VMEM_GB:-4}"
OUT="${HV_UNCAPPED_OUT:-$REPO/test_runs/mwpc_stage_c/p3_uncapped_2000_2200_v1}"
SELECTION="$REPO/test/mwpc-stage-c/p3_wiretruth_1000_selection.csv"
BUILD="$REPO/test/mwpc-avalanche/build-uncapped"
EXE="$BUILD/mwpc_p3_noFW_penning057_uncapped_v1"

if [[ ! "$N" =~ ^[0-9]+$ ]] || (( N<1 || N>20 ));then
  echo "Usage: bash run_uncapped_local.sh {1|5|20} (N=1..20)" >&2;exit 2
fi
if [[ ! "$TIMEOUT_S" =~ ^[0-9]+$ ]] || (( TIMEOUT_S<60 || TIMEOUT_S>86400 ));then
  echo "HV_UNCAPPED_TIMEOUT_S must be integer 60..86400" >&2;exit 2
fi
if [[ ! "$VMEM_GB" =~ ^[0-9]+$ ]] || (( VMEM_GB<2 || VMEM_GB>64 ));then
  echo "HV_UNCAPPED_VMEM_GB must be integer 2..64" >&2;exit 2
fi
command -v timeout >/dev/null || { echo "GNU timeout required";exit 2; }
command -v root >/dev/null || { echo "ROOT required for analysis";exit 2; }
[[ -s "$SELECTION" ]] || { echo "Frozen 1000 selection missing: $SELECTION" >&2;exit 2; }

[[ -f "$HOME/na6p/env.sh" ]] && source "$HOME/na6p/env.sh"
GARFIELD_SETUP="$HOME/na6p/install/garfieldpp/share/Garfield/setupGarfield.sh"
[[ -f "$GARFIELD_SETUP" ]] && source "$GARFIELD_SETUP"

cmake -S "$REPO/test/mwpc-avalanche" -B "$BUILD" \
  -DGarfield_DIR="$HOME/na6p/install/garfieldpp/lib/cmake/Garfield"
cmake --build "$BUILD" --target mwpc_p3_noFW_penning057_uncapped_v1 -j2

mkdir -p "$OUT"
SOURCE_SHA="$(sha256sum "$SELECTION" | cut -d ' ' -f 1)"
BINARY_SHA="$(sha256sum "$EXE" | cut -d ' ' -f 1)"
GIT_REV="$(git rev-parse HEAD)"
LOCK="input_sha256=$SOURCE_SHA
binary_sha256=$BINARY_SHA
git_commit=$GIT_REV
voltages=2000,2050,2100,2150,2200
avalanche_size_limit=DISABLED
per_muon_timeout_seconds=$TIMEOUT_S
process_virtual_memory_GB=$VMEM_GB
physics=P3_noFW;ArCO2_70_30;Penning_r0.57;B0;diameter_30um"
if [[ -f "$OUT/config.lock" ]];then
  if [[ "$(cat "$OUT/config.lock")" != "$LOCK" ]];then
    echo "INPUT/BINARY/GUARDS DIFFER from sealed campaign." >&2
    echo "Use a NEW HV_UNCAPPED_OUT, do not reuse existing results." >&2;exit 3
  fi
else
  printf '%s\n' "$LOCK" > "$OUT/config.lock"
fi
echo "UNLIMITED Garfield avalanche size. N=$N muons per voltage."
echo "External guard only: $VMEM_GB GB virtual address space, $TIMEOUT_S seconds/muon."
echo "Data: $OUT"

run_muon() {
  local V="$1" i="$2" tag dir
  printf -v tag 'muon_%04d' "$i"
  dir="$OUT/HV_${V}/$tag"
  if [[ -f "$dir/complete.txt" ]];then
    (cd "$dir" && sha256sum -c SHA256SUMS >/dev/null) || {
      echo "CHECKSUM FAILED at $dir" >&2;exit 3;
    }
    echo "SKIP verified V=$V $tag";return 0
  fi
  if [[ -d "$dir" ]];then
    local quarantine="$dir.incomplete.$(date +%Y%m%dT%H%M%S)"
    mv "$dir" "$quarantine"
    echo "Moved incomplete result to $quarantine"
  fi
  mkdir -p "$dir"
  {
    head -n 1 "$SELECTION"
    sed -n "$((i+2))p" "$SELECTION"
  } > "$dir/input.csv"
  [[ "$(wc -l < "$dir/input.csv")" -eq 2 ]] || {
    echo "Failed to select genuine Geant4 muon index $i" >&2;exit 4
  }
  echo "RUN V=$V index=$i; NO Garfield avalanche size limit."
  local status=0
  if (
    ulimit -v "$(( VMEM_GB * 1024 * 1024 ))"
    if [[ -x /usr/bin/time ]];then
      /usr/bin/time -v -o "$dir/resources.log" \
        timeout --signal=TERM --kill-after=20s "$TIMEOUT_S" \
        "$EXE" --wiretruth-uncapped-p3-noFW "$dir/input.csv" "$dir/microscopic" "$V"
    else
      timeout --signal=TERM --kill-after=20s "$TIMEOUT_S" \
        "$EXE" --wiretruth-uncapped-p3-noFW "$dir/input.csv" "$dir/microscopic" "$V"
    fi
  ) >"$dir/physics.log" 2>&1; then
    :
  else
    status="$?"
    printf 'UNFINISHED voltage=%s index=%s exit_status=%s\n' "$V" "$i" "$status" \
      > "$dir/UNFINISHED.txt"
    echo "ABORTED V=$V index=$i (exit=$status); physical size not measured." >&2
    echo "See $dir/physics.log and resources.log; no completed result is accepted." >&2
    tail -n 24 "$dir/physics.log" >&2
    exit 5
  fi
  [[ -s "$dir/microscopic_crossing.csv" &&
     -s "$dir/microscopic_seed.csv" &&
     -s "$dir/microscopic_wire_muon.csv" &&
     -s "$dir/microscopic_timing.csv" ]] || {
       echo "Incomplete output at $dir" >&2;exit 5;
     }
  [[ "$(wc -l < "$dir/microscopic_crossing.csv")" -eq 2 ]] || {
    echo "Crossing count mismatch $dir" >&2;exit 5;
  }
  if ! awk -F, '
    NR==1 {for(i=1;i<=NF;i++) {
      if($i=="capped")cap=i;
      if($i=="successful")ok=i;
      if($i=="ion_births")ion=i;
      if($i=="electron_endpoints")nend=i;
      if($i=="collected_e")ne=i;
    }if(!cap||!ok||!ion||!nend||!ne)bad=1;next}
    NR>1 {if(bad||($cap+0)!=0||($ok+0)!=1||($ion+0)<0||
              ($nend+0)<0||($ne+0)<0)bad=1}
    END{if(bad||NR<2)exit 1}
  ' "$dir/microscopic_seed.csv";then
    echo "Seed verification failed $dir" >&2;exit 5;
  fi
  (
    cd "$dir"
    sha256sum input.csv physics.log microscopic_*.csv > SHA256SUMS
    [[ ! -e resources.log ]] || sha256sum resources.log >> SHA256SUMS
  )
  echo "V=$V index=$i COMPLETE uncapped" > "$dir/complete.txt"
  echo "PASS V=$V muon=$i"
}

for V in "${VOLTAGES[@]}";do
  for (( i=0;i<N;i++ ));do
    run_muon "$V" "$i"
  done
done

for V in "${VOLTAGES[@]}";do
  outV="$OUT/HV_$V"
  for kind in crossing seed;do
    merge="$outV/merged_$kind.csv"
    head -n1 "$outV/muon_0000/microscopic_$kind.csv" > "$merge"
    for ((i=0;i<N;i++));do
      printf -v tag 'muon_%04d' "$i"
      [[ -f "$outV/$tag/complete.txt" ]] || {
        echo "Missing verified full muon $outV/$tag" >&2;exit 5;
      }
      tail -n +2 "$outV/$tag/microscopic_$kind.csv" >> "$merge"
    done
  done
  [[ "$(wc -l < "$outV/merged_crossing.csv")" -eq "$((N+1))" ]] || {
    echo "Merged crossing mismatch HV=$V" >&2;exit 5;
  }
done

echo "ROOT/C++ statistics from ALL five COMPLETE voltage ensembles"
root -l -b -q "$HERE/AnalyzeUncappedHV_ROOT.C(\"$OUT\")"
echo "COMPLETE: $OUT/analysis/"
