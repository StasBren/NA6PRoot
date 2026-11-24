#!/usr/bin/env bash
# S-ladder (spatial-description ablation) on real 2100 V development muons, T1 timing baseline.
# READ-ONLY on archives; frozen exporter, frozen replay and frozen LUT unchanged; no new avalanche.
#
#   bash run_s_ladder.sh check
#   bash run_s_ladder.sh prepare [INDEX]        one non-sealed development muon: frozen export -> T1 export -> S0 replay
#   bash run_s_ladder.sh validate RUN           grid executable vs frozen exporter: same start points -> same paths
#   bash run_s_ladder.sh green RUN [TAG]        Green grid paths + single-ion table into RUN/green_TAG (never overwritten)
#   bash run_s_ladder.sh closure RUN [TAG]      G only (numerical closure) -> RUN/muon_*/G_TAG
#   bash run_s_ladder.sh closure-compare RUN TAG[,TAG2,...]   S0 vs G_TAG(s): shaped AND unshaped, ion only + electron step
#   bash run_s_ladder.sh synth RUN [TAG]        full S-variant set from green_TAG -> RUN/muon_*/<variant>_TAG
#   bash run_s_ladder.sh compare RUN [TAG]      S0 vs all variants of TAG: 1) ion only  2) + instantaneous electron step
#   bash run_s_ladder.sh kernel RUN [TAG]       harmonic kernels from green_TAG; K1/K2/K1A/K2A from cloud moments only,
#                                               timed (s_kernel.py); outputs <model>kern_TAG next to the S-variants
#   bash run_s_ladder.sh report RUN [TAG]       compare vs S0 AND vs G_TAG (ion only, + electron step) -> s_ladder_report.md
#   bash run_s_ladder.sh plan RUN [TAG]         choose representative non-sealed development muons + cost estimate
# TAG defaults to r${MWPC_S_NR:-14}t${MWPC_S_NTHETA:-72}; a different grid = a different TAG, so results coexist.
#
# Environment: as run_t_ladder.sh (MWPC_EXPORT_BIN, MWPC_ION_MOBILITY, MWPC_SPECTRAL_LUT, MWPC_NUMW_DIR) plus
#   MWPC_GREEN_BIN   built p3_green_grid_paths (spatial-ablation/CMakeLists.txt, FROZEN_DIR = wire branch)
# Optional: MWPC_S_NR (14), MWPC_S_NTHETA (72), MWPC_S_DSUB (10), MWPC_S_VARIANTS, MWPC_SLADDER_OUT, MWPC_HALF_STRIPS (5)
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
HS="${MWPC_HALF_STRIPS:-5}"; NR="${MWPC_S_NR:-14}"; NTH="${MWPC_S_NTHETA:-72}"; DSUB="${MWPC_S_DSUB:-10}"
DEFTAG="r${NR}t${NTH}"; [[ $DSUB == 10 ]] || DEFTAG="${DEFTAG}d${DSUB}"
VARS="${MWPC_S_VARIANTS:-G,S1,S1A,S2,S2A,S3,S3A,SR,S2R,S3R,SW,SWW,SRW,S3W,K1,K1A,K2,K2A}"
KMODELS="K1,K2,K1A,K2A"
suffixed(){ local o="" v; for v in ${1//,/ }; do o="$o,${v}_$2"; done; echo "${o#,}"; }
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
  validate)
    envs MWPC_GREEN_BIN; RUN="${2:?RUN dir}"
    for R in "$RUN"/muon_*; do
      V="$R/transport_validation_$(date +%H%M%S)"; mkdir -p "$V"
      python3 "$HERE/s_validate_transport.py" points "$R/T1_export" 200 "$V/points.csv"
      "$MWPC_GREEN_BIN" "$MWPC_ION_MOBILITY" "$V/paths" --points "$V/points.csv" > "$V/paths.log" 2>&1 || fail "grid tool ($V/paths.log)"
      python3 "$HERE/s_validate_transport.py" compare "$R/T1_export" "$V/points.csv" "$V/paths" | tee "$V/result.txt"
      grep -q '^PASS' "$V/result.txt" || fail "grid executable does not reproduce the frozen export ($V/result.txt)"
    done;;
  green)
    envs MWPC_GREEN_BIN; RUN="${2:?RUN dir}"; TAG="${3:-$DEFTAG}"; GD="$RUN/green_$TAG"
    [[ ! -e "$GD" && ! -e "$RUN/green_paths_$TAG" ]] || fail "$GD exists (choose another TAG)"
    python3 "$HERE/s_grid_spec.py" "$RUN"/muon_*/T1_export --nr "$NR" --half-strips "$HS" | tee "$RUN/grid_spec_$TAG.txt"
    W=$(sed -n 's/^--wires \([^ ]*\) --r-um.*/\1/p' "$RUN/grid_spec_$TAG.txt"); RU=$(sed -n 's/.*--r-um \([^ ]*\)$/\1/p' "$RUN/grid_spec_$TAG.txt")
    DM=$(sed -n 's/.*--d-max-mm \([0-9.]*\)$/\1/p' "$RUN/grid_spec_$TAG.txt")
    [[ -n "$W" && -n "$RU" && -n "$DM" ]] || fail "grid spec parse"
    t0=$(date +%s)
    "$MWPC_GREEN_BIN" "$MWPC_ION_MOBILITY" "$RUN/green_paths_$TAG" --wires "$W" --r-um "$RU" --ntheta "$NTH" > "$RUN/green_paths_$TAG.log" 2>&1 || fail "grid paths ($RUN/green_paths_$TAG.log)"
    grep -q '^PASS green grid' "$RUN/green_paths_$TAG.log" || fail "grid paths did not PASS"
    t1=$(date +%s); echo "green_paths_${TAG}_s=$((t1 - t0))" >> "$RUN/timing.txt"
    python3 "$HERE/s_green_table.py" "$MWPC_SPECTRAL_LUT" "$RUN/green_paths_$TAG" "$GD" --replay-dir "$MWPC_NUMW_DIR" --d-max-mm "$DM" --d-sub "$DSUB" | tee "$RUN/green_table_$TAG.log"
    echo "green_table_${TAG}_s=$(( $(date +%s) - t1 ))" >> "$RUN/timing.txt";;
  closure|synth)
    RUN="${2:?RUN dir}"; TAG="${3:-$DEFTAG}"; GD="$RUN/green_$TAG"; [[ -s "$GD/axes.json" ]] || fail "run: green $RUN $TAG"
    V="G"; [[ $1 == synth ]] && V="$VARS"
    for R in "$RUN"/muon_*; do
      MISS=""; for v in ${V//,/ }; do [[ -e "$R/${v}_$TAG" ]] || MISS="$MISS,$v"; done; MISS="${MISS#,}"
      [[ -n "$MISS" ]] || { echo "$R: all variants of $TAG present"; continue; }
      python3 "$HERE/s_synthesize.py" "$GD" "$R/T1_export" "$R/S0/replay/waveforms.csv" "$R" --variants "$MISS" --suffix "_$TAG" | tee -a "$R/synth_$TAG.log"
    done;;
  closure-compare)
    RUN="${2:?RUN dir}"; TAGS="${3:-$DEFTAG}"; S=$(date +%H%M%S); V="S0"
    for t in ${TAGS//,/ }; do V="$V,G_$t"; done
    python3 "$TA/t_shaped_compare.py" "$RUN" --variants "$V" --tp 100,200 --no-prompt --plots --out "$RUN/closure_1_ion_only_${TAGS//,/+}_$S"
    python3 "$TA/t_shaped_compare.py" "$RUN" --variants "$V" --tp 100,200 --plots --out "$RUN/closure_2_ion_plus_electron_step_${TAGS//,/+}_$S"
    echo "return: $RUN/closure_*_$S/t_ladder_summary.txt and the PNGs";;
  kernel)
    RUN="${2:?RUN dir}"; TAG="${3:-$DEFTAG}"; GD="$RUN/green_$TAG"; KD="$RUN/kernel_$TAG"
    [[ -s "$GD/axes.json" ]] || fail "run: green $RUN $TAG"
    [[ -e "$KD" ]] || python3 "$HERE/s_kernel.py" build "$GD" "$KD" | tee "$RUN/kernel_build_$TAG.log"
    for R in "$RUN"/muon_*; do
      python3 "$HERE/s_kernel.py" moments "$R/T1_export" "$R/cloud_moments.csv"
      for m in ${KMODELS//,/ }; do
        python3 "$HERE/s_kernel.py" apply "$KD" "$R/cloud_moments.csv" "$R/S0/replay/waveforms.csv" "$R" --model "$m" --suffix "_$TAG" --repeat 5 \
          | tee -a "$R/kernel_$TAG.log"
      done
      python3 -c "import json,sys; print('S0 frozen replay wall_s', json.load(open(sys.argv[1]))['wall_s'])" "$R/S0/replay/summary.json" | tee -a "$R/kernel_$TAG.log"
    done;;
  report)
    RUN="${2:?RUN dir}"; TAG="${3:-$DEFTAG}"; S=$(date +%H%M%S)
    V="$(suffixed "$VARS" "$TAG")"
    [[ -e "$(ls -d "$RUN"/muon_*/K2kern_"$TAG" 2>/dev/null | head -1)" ]] && V="$V,$(suffixed "$(echo "$KMODELS" | sed 's/\([A-Z0-9]*\)/\1kern/g')" "$TAG")"
    VG="G_$TAG,$(echo "$V" | tr ',' '\n' | grep -v "^G_$TAG\$" | paste -sd,)"
    for sig in ion_only electron_step; do
      NP=(); [[ $sig == ion_only ]] && NP=(--no-prompt)
      python3 "$TA/t_shaped_compare.py" "$RUN" --variants "S0,$V" --tp 100,200 "${NP[@]}" --out "$RUN/report_${TAG}_${S}/vs_S0_$sig"
      python3 "$TA/t_shaped_compare.py" "$RUN" --variants "$VG" --tp 100,200 "${NP[@]}" --out "$RUN/report_${TAG}_${S}/vs_G_$sig"
      python3 "$HERE/s_ladder_report.py" --vs-s0 "$RUN/report_${TAG}_${S}/vs_S0_$sig" --vs-g "$RUN/report_${TAG}_${S}/vs_G_$sig" \
        --tag "$TAG" --out "$RUN/report_${TAG}_${S}/s_ladder_report_$sig.md"
    done
    echo "return: $RUN/report_${TAG}_${S}/s_ladder_report_*.md (+ the vs_* directories)";;
  plan)
    RUN="${2:?RUN dir}"; TAG="${3:-$DEFTAG}"
    python3 "$HERE/s_plan_expansion.py" "$ARCHIVE" "$DEV" "$SEALED" "$RUN" --tag "$TAG" --nr "$NR" --ntheta "$NTH" | tee "$RUN/plan_$TAG.txt";;
  compare)
    RUN="${2:?RUN dir}"; TAG="${3:-$DEFTAG}"; S=$(date +%H%M%S); V="S0,$(suffixed "$VARS" "$TAG")"
    python3 "$TA/t_shaped_compare.py" "$RUN" --variants "$V" --tp 100,200 --no-prompt --plots --out "$RUN/compare_1_ion_only_${TAG}_$S"
    python3 "$TA/t_shaped_compare.py" "$RUN" --variants "$V" --tp 100,200 --plots --out "$RUN/compare_2_ion_plus_electron_step_${TAG}_$S"
    echo "return: $RUN/compare_*_${TAG}_$S/t_ladder_summary.txt, $RUN/grid_spec_$TAG.txt, $RUN/green_$TAG/axes.json, $RUN/muon_*/s_synthesize_*.json";;
  *) sed -n 2,24p "$0"; exit 1;;
esac
