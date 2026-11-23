#!/usr/bin/env bash
# T-ladder timing ablation (T0 / T1 / per-wire time / per-muon time) on frozen 2100 V K256 archives.
# READ-ONLY on archives. Uses the frozen exporter and the frozen replay unchanged. No new Garfield++ avalanche run.
# Development muons only (development20.csv); sealed slots are dropped before any file of theirs is touched.
#
#   bash run_t_ladder.sh check                 # tool presence + synthetic-free gates on the first muon (identity bytes)
#   bash run_t_ladder.sh run [N]               # first N development muons (default all 20)
#   bash run_t_ladder.sh compare DIR           # shaped comparison of a finished run directory
#
# Required environment (no defaults are guessed for physics inputs):
#   MWPC_EXPORT_BIN     path to mwpc_stage_c_p3_export_full_spectral_paths (built on the wire branch)
#   MWPC_ION_MOBILITY   the SAME mobility file used for the frozen 2100 V replays
#   MWPC_SPECTRAL_LUT   spectral weighting LUT .npz built by p3_spectral_weighting_lut.py build
#   MWPC_NUMW_DIR       directory holding p3_replay_full_spectral_muon.py and its imports (numerical-weighting/)
# Optional: MWPC_ARCHIVE_2100_K256, MWPC_SIGNAL40_OUTPUT, MWPC_TLADDER_OUT, MWPC_HALF_STRIPS (5)
set -Eeo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../../../.." && pwd)"
set -u
umask 077
ARCHIVE="${MWPC_ARCHIVE_2100_K256:-$HOME/na6p/cathode-overnight-2100/test_runs/mwpc_stage_c/overnight_2100_weighted_k256_v1}"
SELECTION="${MWPC_SIGNAL40_OUTPUT:-$HOME/NA6PRoot-signal40-2100/test_runs/mwpc_stage_c/phase_space_2100_signal40_v1}"
DEV="$SELECTION/development20.csv"; SEALED="$SELECTION/sealed_validation20.csv"
OUTROOT="${MWPC_TLADDER_OUT:-$REPO/test_runs/mwpc_stage_c/t_ladder_2100}"
HS="${MWPC_HALF_STRIPS:-5}"
fail(){ echo "T-LADDER FAILED: $*" >&2; exit 2; }
need(){ [[ -n "${!1:-}" ]] || fail "set $1"; [[ -e "${!1}" ]] || fail "$1=${!1} does not exist"; }
idx(){ awk -F, 'NR==1{for(i=1;i<=NF;i++)if($i=="index")c=i; if(!c)exit 3; next} NF{printf "%04d\n",$c}' "$1"; }
dev_list(){ [[ -s "$DEV" && -s "$SEALED" ]] || fail "selection lists missing"
  comm -23 <(idx "$DEV" | sort) <(idx "$SEALED" | sort); }
verified(){ local d="$ARCHIVE/muon_$1"
  [[ -s "$d/ion_births.root" && -s "$d/SHA256SUMS" ]] && (cd "$d" && sha256sum -c SHA256SUMS >/dev/null 2>&1); }

one_muon(){ # $1 index, $2 run dir
  local m="$1" R="$2/muon_$1" t0 t1
  verified "$m" || { echo "skip $m: archive not verified"; return 0; }
  mkdir -p "$R"
  t0=$(date +%s.%N)
  "$MWPC_EXPORT_BIN" "$ARCHIVE/muon_$m/ion_births.root" "$MWPC_ION_MOBILITY" "$R/frozen_export" > "$R/export.log" 2>&1
  t1=$(date +%s.%N); echo "export_s=$(echo "$t1-$t0"|bc)" >> "$R/timing.txt"
  declare -A MODE=([T0]=identity [T1]=primary [TW]=wire [TM]=muon)
  for v in T0 T1 TW TM; do
    python3 "$HERE/t_time_collapse_export.py" "$R/frozen_export" "$R/$v/export" --mode "${MODE[$v]}" > "$R/$v.collapse.log"
    [[ $v != T0 ]] || grep -q '^PASS identity' "$R/$v.collapse.log" || fail "identity gate $m"
    t0=$(date +%s.%N)
    (cd "$MWPC_NUMW_DIR" && python3 p3_replay_full_spectral_muon.py "$MWPC_SPECTRAL_LUT" "$R/$v/export" "$R/$v/replay" --half-strips "$HS") > "$R/$v.replay.log" 2>&1
    t1=$(date +%s.%N); echo "replay_${v}_s=$(echo "$t1-$t0"|bc)" >> "$R/timing.txt"
    python3 "$HERE/t_prompt_electron_term.py" "$MWPC_SPECTRAL_LUT" "$R/$v/export" "$R/$v/prompt" --replay-dir "$MWPC_NUMW_DIR" --half-strips "$HS" > "$R/$v.prompt.log" 2>&1
  done
  rm -rf "$R/frozen_export/vertices.f64" "$R"/T?/export/vertices.f64   # bulky intermediates, regenerable; SHA kept in ablation.json
  echo "done $m"
}

case "${1:-}" in
  check)
    for v in MWPC_EXPORT_BIN MWPC_ION_MOBILITY MWPC_SPECTRAL_LUT MWPC_NUMW_DIR; do need $v; done
    python3 -c "import numpy, scipy" || fail "numpy/scipy"
    m=$(dev_list | head -1); [[ -n "$m" ]] || fail "no development muon"
    echo "PASS tools present; first development muon $m (sealed excluded: $(idx "$SEALED" | wc -l))";;
  run)
    for v in MWPC_EXPORT_BIN MWPC_ION_MOBILITY MWPC_SPECTRAL_LUT MWPC_NUMW_DIR; do need $v; done
    N="${2:-20}"
    RUN="$OUTROOT/run_$(date +%Y%m%d_%H%M%S)"; [[ ! -e "$RUN" ]] || fail exists; mkdir -p "$RUN"
    { echo "git_sha=$(git -C "$REPO" rev-parse HEAD)"; echo "git_dirty=$(git -C "$REPO" status --porcelain | wc -l)"
      echo "export_bin_sha256=$(sha256sum "$MWPC_EXPORT_BIN"|cut -d' ' -f1)"
      echo "mobility_sha256=$(sha256sum "$MWPC_ION_MOBILITY"|cut -d' ' -f1)"
      echo "lut_sha256=$(sha256sum "$MWPC_SPECTRAL_LUT"|cut -d' ' -f1)"
      echo "replay_sha256=$(sha256sum "$MWPC_NUMW_DIR/p3_replay_full_spectral_muon.py"|cut -d' ' -f1)"
      echo "archive_root=$ARCHIVE"; echo "sealed_list_sha256=$(sha256sum "$SEALED"|cut -d' ' -f1)"; } > "$RUN/provenance.txt"
    dev_list | head -n "$N" > "$RUN/muons.txt"
    while read -r m; do one_muon "$m" "$RUN"; done < "$RUN/muons.txt"
    echo "RUN DIR $RUN";;
  compare)
    D="${2:?run dir}"
    python3 "$HERE/t_shaped_compare.py" "$D" --variants T0,T1,TW,TM --tp 100,200 --out "$D/compare_$(date +%H%M%S)"
    python3 "$HERE/t_shaped_compare.py" "$D" --variants T0,T1,TW,TM --tp 100,200 --no-prompt --out "$D/compare_ion_only_$(date +%H%M%S)";;
  *) sed -n 2,16p "$0"; exit 1;;
esac
