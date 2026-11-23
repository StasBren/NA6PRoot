#!/usr/bin/env bash
# T-ladder timing ablation (T0 / T1 / per-wire time / per-muon time) on frozen 2100 V K256 archives.
# READ-ONLY on archives. Uses the frozen exporter and the frozen replay unchanged. No new Garfield++ avalanche run.
# Development muons only (development20.csv); sealed slots are dropped before any file of theirs is touched.
#
#   bash run_t_ladder.sh check                 # tool presence, selection lists, first development muon
#   bash run_t_ladder.sh one [INDEX]           # ONE development muon (default: first in development20.csv) -- do this first
#   bash run_t_ladder.sh run N                 # first N development muons (only after one muon has been reviewed)
#   bash run_t_ladder.sh compare DIR           # ion-only comparison FIRST, then ion + instantaneous electron step
#
# Variants: MWPC_TLADDER_VARIANTS (default T0,T1,TW; TM = one time per muon is an optional extreme control).
# Shaping 100/200 ns is PROVISIONAL (NA60+ VMM3a peaking time not yet verified); CR-RC^3 is a stand-in.
#
# Required environment (no defaults are guessed for physics inputs):
#   MWPC_EXPORT_BIN     path to mwpc_stage_c_p3_export_full_spectral_paths (built on the wire branch)
#   MWPC_ION_MOBILITY   the SAME mobility file used for the frozen 2100 V replays
#   MWPC_SPECTRAL_LUT   spectral weighting LUT .npz built by p3_spectral_weighting_lut.py build
#   MWPC_NUMW_DIR       directory holding p3_replay_full_spectral_muon.py and its imports (numerical-weighting/)
# Optional: MWPC_ARCHIVE_2100_K256, MWPC_SIGNAL40_OUTPUT, MWPC_TLADDER_OUT, MWPC_HALF_STRIPS (5),
#   MWPC_FROZEN_REPLAY_DIR  an existing frozen replay output of the SAME muon: T0 waveforms.csv must then be identical
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
VARIANTS="${MWPC_TLADDER_VARIANTS:-T0,T1,TW}"
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
  for v in ${VARIANTS//,/ }; do
    [[ -n "${MODE[$v]:-}" ]] || fail "unknown variant $v"
    python3 "$HERE/t_time_collapse_export.py" "$R/frozen_export" "$R/$v/export" --mode "${MODE[$v]}" > "$R/$v.collapse.log"
    [[ $v != T0 ]] || grep -q '^PASS identity' "$R/$v.collapse.log" || fail "identity gate $m"
    t0=$(date +%s.%N)
    (cd "$MWPC_NUMW_DIR" && python3 p3_replay_full_spectral_muon.py "$MWPC_SPECTRAL_LUT" "$R/$v/export" "$R/$v/replay" --half-strips "$HS") > "$R/$v.replay.log" 2>&1
    t1=$(date +%s.%N); echo "replay_${v}_s=$(echo "$t1-$t0"|bc)" >> "$R/timing.txt"
    if [[ $v == T0 && -n "${MWPC_FROZEN_REPLAY_DIR:-}" ]]; then
      if cmp -s "$R/T0/replay/waveforms.csv" "$MWPC_FROZEN_REPLAY_DIR/waveforms.csv"; then echo "PASS T0 waveforms identical to frozen replay" | tee -a "$R/gates.txt"
      else echo "FAIL T0 waveforms differ from $MWPC_FROZEN_REPLAY_DIR (different LUT, mobility or replay?)" | tee -a "$R/gates.txt"; fi
    fi
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
  one|run)
    for v in MWPC_EXPORT_BIN MWPC_ION_MOBILITY MWPC_SPECTRAL_LUT MWPC_NUMW_DIR; do need $v; done
    command -v bc >/dev/null || fail "bc missing (sudo apt install bc)"
    RUN="$OUTROOT/$1_$(date +%Y%m%d_%H%M%S)"; [[ ! -e "$RUN" ]] || fail exists; mkdir -p "$RUN"
    if [[ $1 == one ]]; then
      M="${2:-$(dev_list | head -1)}"; M=$(printf %04d "$((10#$M))")
      dev_list | grep -qx "$M" || fail "$M is not a non-sealed development muon"
      echo "$M" > "$RUN/muons.txt"
    else
      N="${2:?N}"; dev_list | head -n "$N" > "$RUN/muons.txt"
    fi
    { echo "git_sha=$(git -C "$REPO" rev-parse HEAD)"; echo "git_dirty=$(git -C "$REPO" status --porcelain | wc -l)"
      echo "variants=$VARIANTS"; echo "half_strips=$HS"
      echo "export_bin_sha256=$(sha256sum "$MWPC_EXPORT_BIN"|cut -d' ' -f1)"
      echo "mobility=$MWPC_ION_MOBILITY"; echo "mobility_sha256=$(sha256sum "$MWPC_ION_MOBILITY"|cut -d' ' -f1)"
      echo "lut=$MWPC_SPECTRAL_LUT"; echo "lut_sha256=$(sha256sum "$MWPC_SPECTRAL_LUT"|cut -d' ' -f1)"
      echo "replay_sha256=$(sha256sum "$MWPC_NUMW_DIR/p3_replay_full_spectral_muon.py"|cut -d' ' -f1)"
      echo "wire_branch_sha=$(git -C "$MWPC_NUMW_DIR" rev-parse HEAD 2>/dev/null || echo unknown)"
      echo "archive_root=$ARCHIVE"; echo "sealed_list_sha256=$(sha256sum "$SEALED"|cut -d' ' -f1)"; } > "$RUN/provenance.txt"
    while read -r m; do one_muon "$m" "$RUN"; done < "$RUN/muons.txt"
    echo "RUN DIR $RUN"
    echo "next: bash $0 compare $RUN";;
  compare)
    D="${2:?run dir}"
    V=$(sed -n 's/^variants=//p' "$D/provenance.txt"); V="${V:-$VARIANTS}"; S=$(date +%H%M%S)
    # 1) ion term only (exact for the frozen RKF paths + spectral LUT), 2) plus the instantaneous electron-step approximation
    python3 "$HERE/t_shaped_compare.py" "$D" --variants "$V" --tp 100,200 --no-prompt --out "$D/compare_1_ion_only_$S"
    python3 "$HERE/t_shaped_compare.py" "$D" --variants "$V" --tp 100,200 --out "$D/compare_2_ion_plus_electron_step_$S"
    echo "return: $D/compare_1_ion_only_$S/t_ladder_summary.txt and $D/compare_2_ion_plus_electron_step_$S/t_ladder_summary.txt";;
  *) sed -n 2,22p "$0"; exit 1;;
esac
