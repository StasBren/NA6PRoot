#!/usr/bin/env bash
# S-ladder (spatial-description ablation) on real 2100 V development muons, T1 timing baseline.
# READ-ONLY on archives; frozen exporter, frozen replay and frozen LUT unchanged; no new avalanche.
#
#   bash run_s_ladder.sh check
#   bash run_s_ladder.sh prepare [INDEX]   one non-sealed development muon: frozen export -> T1 export -> S0 replay
#   bash run_s_ladder.sh green RUN         Green grid paths (frozen transport code) + single-ion table, sized from RUN
#   bash run_s_ladder.sh synth RUN         S-variants from the Green table for every prepared muon
#   bash run_s_ladder.sh compare RUN       shaped comparison vs S0: 1) ion only  2) + instantaneous electron step
#
# Environment: as run_t_ladder.sh (MWPC_EXPORT_BIN, MWPC_ION_MOBILITY, MWPC_SPECTRAL_LUT, MWPC_NUMW_DIR) plus
#   MWPC_GREEN_BIN   built p3_green_grid_paths (spatial-ablation/CMakeLists.txt, FROZEN_DIR = wire branch)
# Optional: MWPC_S_NR (radial nodes, 14), MWPC_S_NTHETA (72), MWPC_S_VARIANTS, MWPC_SLADDER_OUT, MWPC_HALF_STRIPS (5)
# Shaping 100/200 ns PROVISIONAL; CR-RC^3 stand-in.
set -Eeo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TA="$(cd "$HERE/../timing-ablation" && pwd)"
REPO="$(cd "$HERE/../../../.." && pwd)"
set -u
umask 077
ARCHIVE="${MWPC_ARCHIVE_2100_K256:-$HOME/na6p/cathode-overnight-2100/test_runs/mwpc_stage_c/overnight_2100_weighted_k256_v1}"
SELECTION="${MWPC_SIGNAL40_OUTPUT:-$HOME/NA6PRoot-signal40-2100/test_runs/mwpc_stage_c/phase_space_2100_signal40_v1}"
DEV="$SELECTION/development20.csv"; SEALED="$SELECTION/sealed_validation20.csv"
OUTROOT="${MWPC_SLADDER_OUT:-$REPO/test_runs/mwpc_stage_c/s_ladder_2100}"
HS="${MWPC_HALF_STRIPS:-5}"; NR="${MWPC_S_NR:-14}"; NTH="${MWPC_S_NTHETA:-72}"
VARS="${MWPC_S_VARIANTS:-G,S1,S2,S3,SR,S2R,S3R,SW,SWW,K1,K2}"
fail(){ echo "S-LADDER FAILED: $*" >&2; exit 2; }
need(){ [[ -n "${!1:-}" ]] || fail "set $1"; [[ -e "${!1}" ]] || fail "$1=${!1} does not exist"; }
idx(){ awk -F, 'NR==1{for(i=1;i<=NF;i++)if($i=="index")c=i; if(!c)exit 3; next} NF{printf "%04d\n",$c}' "$1"; }
dev_list(){ [[ -s "$DEV" && -s "$SEALED" ]] || fail "selection lists missing"; comm -23 <(idx "$DEV" | sort) <(idx "$SEALED" | sort); }
verified(){ local d="$ARCHIVE/muon_$1"; [[ -s "$d/ion_births.root" && -s "$d/SHA256SUMS" ]] && (cd "$d" && sha256sum -c SHA256SUMS >/dev/null 2>&1); }
envs(){ for v in MWPC_EXPORT_BIN MWPC_ION_MOBILITY MWPC_SPECTRAL_LUT MWPC_NUMW_DIR "$@"; do need $v; done; }

case "${1:-}" in
  check)
    envs MWPC_GREEN_BIN
    bash "$TA/run_t_ladder.sh" check
    usage=$("$MWPC_GREEN_BIN" 2>&1 || true)
    [[ $usage == *"usage: p3_green_grid_paths"* ]] || fail "MWPC_GREEN_BIN is not p3_green_grid_paths"
    echo "PASS S-ladder dependencies";;
  prepare)
    envs
    M="${2:-$(dev_list | head -1)}"; M=$(printf %04d "$((10#$M))")
    dev_list | grep -qx "$M" || fail "$M is not a non-sealed development muon"
    verified "$M" || fail "archive muon_$M not verified"
    RUN="${MWPC_S_RUN:-$OUTROOT/run_$(date +%Y%m%d_%H%M%S)}"; R="$RUN/muon_$M"; [[ ! -e "$R" ]] || fail "$R exists"; mkdir -p "$R"
    { echo "git_sha=$(git -C "$REPO" rev-parse HEAD)"; echo "git_dirty=$(git -C "$REPO" status --porcelain | wc -l)"
      echo "mobility=$MWPC_ION_MOBILITY"; echo "mobility_sha256=$(sha256sum "$MWPC_ION_MOBILITY"|cut -d' ' -f1)"
      echo "lut_sha256=$(sha256sum "$MWPC_SPECTRAL_LUT"|cut -d' ' -f1)"; echo "export_bin_sha256=$(sha256sum "$MWPC_EXPORT_BIN"|cut -d' ' -f1)"
      echo "replay_sha256=$(sha256sum "$MWPC_NUMW_DIR/p3_replay_full_spectral_muon.py"|cut -d' ' -f1)"
      echo "sealed_list_sha256=$(sha256sum "$SEALED"|cut -d' ' -f1)"; } > "$RUN/provenance_prepare_$M.txt"
    "$MWPC_EXPORT_BIN" "$ARCHIVE/muon_$M/ion_births.root" "$MWPC_ION_MOBILITY" "$R/frozen_export" > "$R/export.log" 2>&1 || fail "exporter ($R/export.log)"
    grep -q '^PASS archive census' "$R/export.log" || fail "exporter did not PASS"
    python3 "$TA/t_time_collapse_export.py" "$R/frozen_export" "$R/T1_export" --mode primary > "$R/T1.collapse.log" || fail "T1 collapse"
    (cd "$MWPC_NUMW_DIR" && python3 p3_replay_full_spectral_muon.py "$MWPC_SPECTRAL_LUT" "$R/T1_export" "$R/S0/replay" --half-strips "$HS") > "$R/S0.replay.log" 2>&1 || fail "S0 replay"
    python3 "$TA/t_prompt_electron_term.py" "$MWPC_SPECTRAL_LUT" "$R/T1_export" "$R/S0/prompt" --replay-dir "$MWPC_NUMW_DIR" --half-strips "$HS" > "$R/S0.prompt.log" 2>&1 || fail "S0 electron step"
    rm -rf "$R/frozen_export/vertices.f64"      # T1_export keeps its own vertices (needed by synth)
    echo "PREPARED $R"; echo "RUN DIR $RUN  (prepare more muons with MWPC_S_RUN=$RUN bash $0 prepare INDEX)";;
  green)
    envs MWPC_GREEN_BIN; RUN="${2:?RUN dir}"; [[ ! -e "$RUN/green" ]] || fail "$RUN/green exists"
    python3 "$HERE/s_grid_spec.py" "$RUN"/muon_*/T1_export --nr "$NR" --half-strips "$HS" | tee "$RUN/grid_spec.txt"
    W=$(sed -n 's/^--wires \([^ ]*\) --r-um.*/\1/p' "$RUN/grid_spec.txt"); RU=$(sed -n 's/.*--r-um \([^ ]*\)$/\1/p' "$RUN/grid_spec.txt")
    DM=$(sed -n 's/.*--d-max-mm \([0-9.]*\)$/\1/p' "$RUN/grid_spec.txt")
    [[ -n "$W" && -n "$RU" && -n "$DM" ]] || fail "grid spec parse"
    t0=$(date +%s)
    "$MWPC_GREEN_BIN" "$MWPC_ION_MOBILITY" "$RUN/green_paths" --wires "$W" --r-um "$RU" --ntheta "$NTH" > "$RUN/green_paths.log" 2>&1 || fail "grid paths ($RUN/green_paths.log)"
    grep -q '^PASS green grid' "$RUN/green_paths.log" || fail "grid paths did not PASS"
    echo "green_paths_s=$(( $(date +%s) - t0 ))" >> "$RUN/timing.txt"
    python3 "$HERE/s_green_table.py" "$MWPC_SPECTRAL_LUT" "$RUN/green_paths" "$RUN/green" --replay-dir "$MWPC_NUMW_DIR" --d-max-mm "$DM" | tee "$RUN/green_table.log"
    echo "green_table_s=$(( $(date +%s) - t0 ))" >> "$RUN/timing.txt";;
  synth)
    RUN="${2:?RUN dir}"; [[ -s "$RUN/green/axes.json" ]] || fail "run green first"
    for R in "$RUN"/muon_*; do
      [[ -e "$R/S/G" ]] && { echo "skip $R (S exists)"; continue; }
      python3 "$HERE/s_synthesize.py" "$RUN/green" "$R/T1_export" "$R/S0/replay/waveforms.csv" "$R/S" --variants "$VARS" | tee "$R/synth.log"
      for v in ${VARS//,/ }; do mv "$R/S/$v" "$R/$v"; done
    done;;
  compare)
    RUN="${2:?RUN dir}"; S=$(date +%H%M%S)
    python3 "$TA/t_shaped_compare.py" "$RUN" --variants "S0,$VARS" --tp 100,200 --no-prompt --plots --out "$RUN/compare_1_ion_only_$S"
    python3 "$TA/t_shaped_compare.py" "$RUN" --variants "S0,$VARS" --tp 100,200 --plots --out "$RUN/compare_2_ion_plus_electron_step_$S"
    echo "return: $RUN/compare_*_$S/t_ladder_summary.txt, $RUN/grid_spec.txt, $RUN/green/axes.json, $RUN/muon_*/S/s_synthesize.json";;
  *) sed -n 2,16p "$0"; exit 1;;
esac
