#!/usr/bin/env bash
# Minimal reference experiment: does per-primary strip recentering + absolute-strip accumulation reproduce the FROZEN S0
# weighting response? Design: ../RECENTERING_DESIGN_REVIEW.md section 5. READ-ONLY on the reference run, archive, LUT.
#
#   bash k2a_recenter_experiment.sh refs RUN TAG [MUON ...]     E1 + E2a on reference muons (existing T1_export, S0,
#                                                              cloud_moments.csv, kernel_TAG): no new transport
#   bash k2a_recenter_experiment.sh far  RUN TAG INDEX ...      E2b on NON-SEALED far-D muons: frozen export + T1 (new
#                                                              files only), then the same comparisons; needs
#                                                              MWPC_EXPORT_BIN, MWPC_ION_MOBILITY, MWPC_K2A_PF (preflight)
# E1  symmetry:  whole export relabelled by m = +1 / -1 strips, frozen replay -> must equal the frozen S0 on the
#                overlapping absolute strips (tolerance 1e-10 rel L2)
# E2  recentering: per-primary groups -> frozen replay + electron step per group -> S0_local on absolute strips;
#                K2A_local = unchanged engine per group; comparisons:
#                  K2A_local vs S0_local      numerical accuracy of the recentered kernel (same windows)
#                  S0_local  vs S0 (frozen)   effect of the per-primary window choice vs the base-centred window
#                  K2A_local vs K2A (frozen kernel, base-centred) and vs S0, on the base strips
# Env: MWPC_SPECTRAL_LUT, MWPC_NUMW_DIR (frozen replay dir), MWPC_K2A_OUT (output root), MWPC_HALF_STRIPS (5).
set -Eeo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; TA="$HERE/../../timing-ablation"; SA="$HERE/../../spatial-ablation"
REPO="$(cd "$HERE/../../../../.." && pwd)"
set -u
umask 077
HS="${MWPC_HALF_STRIPS:-5}"; OUTROOT="${MWPC_K2A_OUT:-$REPO/test_runs/mwpc_stage_c/k2a_campaign_2100}"
ARCHIVE="${MWPC_ARCHIVE_2100_K256:-$HOME/na6p/cathode-overnight-2100/test_runs/mwpc_stage_c/overnight_2100_weighted_k256_v1}"
fail(){ echo "RECENTER EXPERIMENT FAILED: $*" >&2; exit 2; }
need(){ [[ -n "${!1:-}" && -e "${!1}" ]] || fail "set $1 (existing path)"; }
need MWPC_SPECTRAL_LUT; need MWPC_NUMW_DIR

replay(){ # EXPORT OUTDIR  (frozen replay + frozen electron step, unchanged)
  (cd "$MWPC_NUMW_DIR" && python3 p3_replay_full_spectral_muon.py "$MWPC_SPECTRAL_LUT" "$1" "$2/replay" --half-strips "$HS") > "$2.replay.log" 2>&1 || return 1
  python3 "$TA/t_prompt_electron_term.py" "$MWPC_SPECTRAL_LUT" "$1" "$2/prompt" --replay-dir "$MWPC_NUMW_DIR" --half-strips "$HS" > "$2.prompt.log" 2>&1
}

e2(){ # R(muon dir with T1_export, cloud_moments.csv) KD E(out dir) HAVE_S0(0/1)
  local R="$1" KD="$2" E="$3" HAVE="$4"
  python3 "$HERE/k2a_group_export.py" "$R/T1_export" "$E/groups" --mode primary | tee "$E/groups.log"
  mkdir -p "$E/S0groups"
  for g in "$E"/groups/g_*; do replay "$g" "$E/S0groups/${g##*/}" || fail "frozen replay of group ${g##*/} ($E/S0groups/${g##*/}.replay.log)"; done
  python3 "$HERE/k2a_local_sum.py" "$E/groups" "$E/S0groups" "$E/S0_local" | tee "$E/local_sum.log"
  python3 "$HERE/k2a_recentered_k2a.py" "$KD" "$R/cloud_moments.csv" "$E/K2A_local" --groups "$E/groups/groups.json" --half-strips "$HS" | tee "$E/k2a_local.log"
  python3 "$HERE/k2a_recenter_compare.py" "$E/S0_local" "$E/K2A_local" --out "$E/cmp_K2Alocal_vs_S0local.json" --label "K2A_local vs S0_local" | tee "$E/compare.txt"
  if [[ $HAVE == 1 ]]; then
    python3 "$HERE/k2a_recenter_compare.py" "$R/S0" "$E/S0_local" --out "$E/cmp_S0local_vs_S0.json" --label "S0_local vs S0" | tee -a "$E/compare.txt"
    python3 "$HERE/k2a_recenter_compare.py" "$R/S0" "$E/K2A_local" --out "$E/cmp_K2Alocal_vs_S0.json" --label "K2A_local vs S0" | tee -a "$E/compare.txt"
    [[ -d "$R/K2Akern_${TAG}" ]] && python3 "$HERE/k2a_recenter_compare.py" "$R/K2Akern_${TAG}" "$E/K2A_local" --out "$E/cmp_K2Alocal_vs_K2Abase.json" --label "K2A_local vs K2A_base" | tee -a "$E/compare.txt"
  fi
  return 0
}

MODE="${1:-}"; RUN="${2:?RUN}"; TAG="${3:?TAG}"; shift 3
KD="$RUN/kernel_$TAG"; [[ -s "$KD/kernel.json" ]] || fail "missing $KD"
X="$OUTROOT/recenter_experiment_$(date +%Y%m%d_%H%M%S)"; mkdir -p "$X"
git -C "$REPO" rev-parse HEAD > "$X/SOURCE_GIT_SHA.txt" 2>/dev/null || true
{ echo "lut_sha256=$(sha256sum "$MWPC_SPECTRAL_LUT" | cut -d' ' -f1)"; echo "replay_sha256=$(sha256sum "$MWPC_NUMW_DIR/p3_replay_full_spectral_muon.py" | cut -d' ' -f1)"
  echo "kernel=$KD"; echo "half_strips=$HS"; } > "$X/provenance.txt"
case "$MODE" in
  refs)
    MS=("$@"); [[ ${#MS[@]} -gt 0 ]] || MS=($(for d in "$RUN"/muon_*; do [[ -s "$d/S0/replay/waveforms.csv" && -s "$d/cloud_moments.csv" ]] && echo "${d##*/muon_}"; done))
    for M in "${MS[@]}"; do
      R="$RUN/muon_$M"; E="$X/muon_$M"; mkdir -p "$E"
      [[ -s "$R/T1_export/vertices.f64" && -s "$R/S0/replay/waveforms.csv" && -s "$R/cloud_moments.csv" ]] || fail "$R is not a prepared reference"
      for m in 1 -1; do        # E1
        python3 "$HERE/k2a_group_export.py" "$R/T1_export" "$E/shift_$m" --mode shift --m-U "$m" --m-V "$m" > /dev/null
        if replay "$E/shift_$m/g_${m}_${m}" "$E/S0shift_$m"; then
          python3 "$HERE/k2a_recenter_compare.py" "$R/S0" "$E/S0shift_$m" --out "$E/E1_shift_$m.json" --tol-l2 1e-10 --label "E1 shift $m" | tee -a "$E/E1.txt" || true
        else
          echo "E1 shift $m: frozen replay refused (see $E/S0shift_$m.replay.log): $(tail -1 "$E/S0shift_$m.replay.log")" | tee -a "$E/E1.txt"
        fi
      done
      e2 "$R" "$KD" "$E" 1
    done;;
  far)
    need MWPC_EXPORT_BIN; need MWPC_ION_MOBILITY; need MWPC_K2A_PF
    [[ $# -ge 1 ]] || fail "list non-sealed indices (from the supplementary selection)"
    for I in "$@"; do
      M=$(printf %04d "$((10#$I))")
      role=$(awk -F, -v s="$M" 'NR==1{for(i=1;i<=NF;i++)h[$i]=i; next} $h["slot"]==s{print $h["role"]}' "$MWPC_K2A_PF/inventory/slots.csv")
      [[ -n "$role" && "$role" != sealed_validation ]] || fail "$M is sealed or unknown (role '$role')"
      E="$X/muon_$M"; R="$E/prepared"; mkdir -p "$R"
      "$MWPC_EXPORT_BIN" "$ARCHIVE/muon_$M/ion_births.root" "$MWPC_ION_MOBILITY" "$R/frozen_export" > "$R/export.log" 2>&1 || fail "exporter $M"
      grep -q '^PASS archive census' "$R/export.log" || fail "exporter did not PASS for $M"
      python3 "$TA/t_time_collapse_export.py" "$R/frozen_export" "$R/T1_export" --mode primary > "$R/T1.log" || fail "T1 collapse $M"
      rm -f "$R/frozen_export/vertices.f64"
      python3 "$SA/s_kernel.py" moments "$R/T1_export" "$R/cloud_moments.csv" > "$R/moments.log"
      # the frozen base-centred S0 is attempted and its own verdict recorded (it refuses if ions land outside the window)
      if replay "$R/T1_export" "$R/S0"; then HAVE=1; echo "frozen base-centred S0: ran" > "$E/base_S0.txt"
      else HAVE=0; echo "frozen base-centred S0 refused: $(grep -m1 -o -E 'selected channels omit the landing strip.*|query outside declared strip phase range' "$R/S0.replay.log" || tail -1 "$R/S0.replay.log")" | tee "$E/base_S0.txt"; fi
      e2 "$R" "$KD" "$E" "$HAVE"
    done;;
  *) sed -n 2,19p "$0"; exit 1;;
esac
echo "return: $X/muon_*/{E1.txt,compare.txt,local_sum.log,base_S0.txt} and $X/muon_*/cmp_*.json"
