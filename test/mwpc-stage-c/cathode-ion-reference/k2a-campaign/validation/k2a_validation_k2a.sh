#!/usr/bin/env bash
# K2A stage of run_k2a_validation.sh for ONE muon (also called directly by the synthetic tests):
#
#   bash k2a_validation_k2a.sh V E KERNEL_DIR KERNEL_ID ION_BIRTHS_ROOT
#
# V  validation muon directory (steps.tsv, timing.tsv; T0/S0_local if the S0_T0 reference exists)
# E  far-mode experiment muon directory (prepared/cloud_moments.csv, groups/groups.json, S0_local, K2A_local)
# env K2AV_EV (k2a_strip_eval), MWPC_BIRTHS_BIN, MWPC_HALF_STRIPS (5), K2AV_REPEAT (20), K2AV_BIRTHS_REPEAT (5)
#
# Production K2A is evaluated WITHOUT any information from S0:
#   K2A      p3_birth_moments (archive -> compact moments) -> k2a_strip_eval with ITS OWN production grouping /
#            recentering rule m = floor(Dbar/P + 1/2) per primary (no --shifts).  PRIMARY physics prediction.
#   K2Aimp   the same moments with the S0 grouping imposed (--shifts groups.json).  DIAGNOSTIC only: isolates the
#            kernel / interpolation error at identical windows; never used as evidence for the prediction.
#   shift_assignment   per-primary (m_U, m_V) of the production rule vs the S0 grouping (k2a_shift_assignment.py):
#            AGREE -> K2A and K2Aimp must be bit-identical (grouping_identity); DISAGREE -> listed with charge, and
#            the effect on waveforms / observables is measured (cmp_K2Aimp_vs_K2A, shaped_refK2Aimp_ionelec).
# Equivalence tests (unchanged in meaning): K2A_cc_export (export-path moments, S0 grouping imposed, as the Python
#   K2A_local) == K2A_local (1e-12); K2A_cc_export == K2Aimp (moment path archive vs export, 1e-7); s_compare_moments.
# Comparisons (k2a_recenter_compare k2a_vs_s0): S0T1_vs_K2A, S0T0_vs_K2A (primary); S0T1_vs_K2Aimp, S0T0_vs_K2Aimp
#   (diagnostic); S0T0_vs_S0T1 (time collapse); K2Aimp_vs_K2A (grouping effect). Shaped (t_shaped_compare) on a
#   common absolute-strip label set (k2a_label_union.py pads strips a variant did not compute with zeros, recorded).
# Time grid: before any comparison, k2a_time_grid.py requires every compared waveform (S0_T1, S0_T0, K2A, K2Aimp,
#   K2A_cc_export, Python K2A_local) to share dt, first bin centre and the bins of the K2A window (window_ns of
#   k2a_strip_eval): any mismatch -> step time_grid FAIL, nothing compared; the shaped inputs are written on exactly
#   those bins, so every variant handed to t_shaped_compare has an identical time column.
set -Eeo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; CAMP="$(cd "$HERE/.." && pwd)"; CIR="$(cd "$CAMP/.." && pwd)"
RC="$CAMP/recenter"; TA="$CIR/timing-ablation"; SA="$CIR/spatial-ablation"; INL="$CIR/k2a-inline"
set -u
[[ $# -eq 5 ]] || { echo "usage: k2a_validation_k2a.sh V E KERNEL_DIR KERNEL_ID ION_BIRTHS_ROOT" >&2; exit 2; }
V="$1"; E="$2"; KD="$3"; KID="$4"; ARCHF="$5"
EV="${K2AV_EV:?set K2AV_EV}"; HS="${MWPC_HALF_STRIPS:-5}"; REP="${K2AV_REPEAT:-20}"; BREP="${K2AV_BIRTHS_REPEAT:-5}"
TF="$V/timing.tsv"
source "$HERE/k2a_validation_lib.sh"
R="$E/prepared"; G="$E/groups/groups.json"; M=$(basename "$V"); M=${M#muon_}
HAVE_T1=0; [[ -s "$E/S0_local/replay/waveforms.csv" && -s "$E/S0_local/prompt/prompt_waveforms.csv" ]] && HAVE_T1=1
HAVE_T0=0; [[ -s "$V/T0/S0_local/replay/waveforms.csv" ]] && grep -q $'^T0_S0_local_census\tOK' "$V/steps.tsv" && HAVE_T0=1

# 1. archive -> compact moments (the frozen archive reader)
mkdir -p "$V/prod"; rc=0
tm births_moments "$MWPC_BIRTHS_BIN" "$ARCHF" "$V/prod/birth_moments.csv" --repeat "$BREP" > "$V/prod/births.log" 2>&1 || rc=$?
if [[ $rc -ne 0 ]] || ! grep -q '^PASS birth moments' "$V/prod/births.log"; then
  vstep births_moments FAIL "p3_birth_moments exit $rc: $(tail -1 "$V/prod/births.log" | head -c 250)" "$V/prod/births.log"; exit 0
fi
vstep births_moments OK "$(grep -m1 '^PASS birth moments' "$V/prod/births.log")" "$V/prod/birth_moments.csv.json"
MOM="$V/prod/birth_moments.csv"
if [[ -s "$R/cloud_moments.csv" ]]; then
  rc=0; python3 "$SA/s_compare_moments.py" "$R/cloud_moments.csv" "$MOM" --report "$V/prod/compare_moments.md" > "$V/prod/compare_moments.log" 2>&1 || rc=$?
  vstep moments_archive_vs_export "$([[ $rc -eq 0 ]] && echo PASS || echo FAIL)" "$(tail -1 "$V/prod/compare_moments.log" | head -c 250)" "$V/prod/compare_moments.md"
else
  vstep moments_archive_vs_export INCONCLUSIVE "no export-path moments (far experiment incomplete)"
fi

# 2. PRIMARY: production grouping rule, nothing from S0
rc=0; tm k2a_eval "$EV" --kernel "$KD" --kernel-id "$KID" --moments "$MOM" --half-strips "$HS" --out "$V/K2A_prod" --repeat "$REP" \
        > "$V/K2A_prod.log" 2>&1 || rc=$?
evstep K2A_prod $rc "$V/K2A_prod.log"

# 3. DIAGNOSTIC: S0 grouping imposed; assignment comparison; equivalence tests
if [[ -s "$G" ]]; then
  rc=0; "$EV" --kernel "$KD" --kernel-id "$KID" --moments "$MOM" --shifts "$G" --half-strips "$HS" --out "$V/K2A_imposed" > "$V/K2A_imposed.log" 2>&1 || rc=$?
  evstep K2A_imposed $rc "$V/K2A_imposed.log"
  rc=0; python3 "$HERE/k2a_shift_assignment.py" "$G" "$V/K2A_prod" "$V/shift_assignment.json" > "$V/shift_assignment.log" 2>&1 || rc=$?
  vstep shift_assignment "$(case $rc in 0) echo AGREE;; 2) echo DISAGREE;; 3) echo INCONCLUSIVE;; *) echo FAIL;; esac)" "$(tail -1 "$V/shift_assignment.log" | head -c 300)" "$V/shift_assignment.json"
  if [[ $rc -eq 0 && -s "$V/K2A_prod/replay/waveforms.csv" && -s "$V/K2A_imposed/replay/waveforms.csv" ]]; then
    if cmp -s "$V/K2A_prod/replay/waveforms.csv" "$V/K2A_imposed/replay/waveforms.csv" && cmp -s "$V/K2A_prod/prompt/prompt_waveforms.csv" "$V/K2A_imposed/prompt/prompt_waveforms.csv"; then
      vstep grouping_identity PASS "assignments agree: production-rule and S0-imposed responses are bit-identical"
    else
      vstep grouping_identity FAIL "assignments agree but the production-rule and S0-imposed responses differ"
    fi
  fi
  rc=0; "$EV" --kernel "$KD" --kernel-id "$KID" --moments "$R/cloud_moments.csv" --shifts "$G" --half-strips "$HS" --out "$V/K2A_cc_export" > "$V/K2A_cc_export.log" 2>&1 || rc=$?
  evstep K2A_cc_export $rc "$V/K2A_cc_export.log"
  [[ -s "$E/K2A_local/replay/waveforms.csv" && -s "$V/K2A_cc_export/replay/waveforms.csv" ]] && stripcmp identity_cc_vs_python "$E/K2A_local" "$V/K2A_cc_export" 1e-12
  [[ -s "$V/K2A_cc_export/replay/waveforms.csv" && -s "$V/K2A_imposed/replay/waveforms.csv" ]] && stripcmp moment_path_archive_vs_export "$V/K2A_cc_export" "$V/K2A_imposed" 1e-7
else
  vstep K2A_imposed INCONCLUSIVE "no S0 grouping (groups.json): far experiment incomplete; diagnostic variant and equivalence tests not run"
fi

# 4. one time grid for every compared waveform (explicit; nothing is compared on a mismatch)
HK=0; [[ -s "$V/K2A_prod/replay/waveforms.csv" ]] && HK=1
HI=0; [[ -s "$V/K2A_imposed/replay/waveforms.csv" ]] && HI=1
CMPW=$(python3 -c "
import json, sys
for p in sys.argv[1:]:
    try:
        print(json.load(open(p))['window_ns']); break
    except (OSError, ValueError, KeyError):
        pass
else:
    print(7000)" "$V/K2A_prod/result.json" "$V/K2A_imposed/result.json")
GV=(); [[ $HAVE_T1 == 1 ]] && GV+=("S0T1=$E/S0_local"); [[ $HAVE_T0 == 1 ]] && GV+=("S0T0=$V/T0/S0_local")
[[ $HK == 1 ]] && GV+=("K2A=$V/K2A_prod"); [[ $HI == 1 ]] && GV+=("K2Aimp=$V/K2A_imposed")
[[ -s "$V/K2A_cc_export/replay/waveforms.csv" ]] && GV+=("K2A_cc_export=$V/K2A_cc_export")
[[ -s "$E/K2A_local/replay/waveforms.csv" ]] && GV+=("K2A_local_python=$E/K2A_local")
if [[ ${#GV[@]} -ge 2 ]]; then
  rc=0; python3 "$HERE/k2a_time_grid.py" "$V/time_grid.json" --window-ns "$CMPW" "${GV[@]}" > "$V/time_grid.log" 2>&1 || rc=$?
  if [[ $rc -ne 0 ]]; then
    vstep time_grid FAIL "$(tail -1 "$V/time_grid.log" | head -c 300)" "$V/time_grid.json"
    exit 0                                                     # no waveform comparison on different time grids
  fi
  vstep time_grid PASS "$(tail -1 "$V/time_grid.log" | head -c 300)" "$V/time_grid.json"
fi

# 5. comparisons
[[ $HAVE_T1 == 1 && $HK == 1 ]] && physcmp S0T1_vs_K2A "$E/S0_local" "$V/K2A_prod"
[[ $HAVE_T0 == 1 && $HK == 1 ]] && physcmp S0T0_vs_K2A "$V/T0/S0_local" "$V/K2A_prod"
[[ $HAVE_T1 == 1 && $HI == 1 ]] && physcmp S0T1_vs_K2Aimp "$E/S0_local" "$V/K2A_imposed"
[[ $HAVE_T0 == 1 && $HI == 1 ]] && physcmp S0T0_vs_K2Aimp "$V/T0/S0_local" "$V/K2A_imposed"
[[ $HAVE_T0 == 1 && $HAVE_T1 == 1 ]] && physcmp S0T0_vs_S0T1 "$V/T0/S0_local" "$E/S0_local"
[[ $HI == 1 && $HK == 1 ]] && physcmp K2Aimp_vs_K2A "$V/K2A_imposed" "$V/K2A_prod"

# 6. shaped observables on one time grid and one absolute-strip label set
VARS=(); [[ $HAVE_T1 == 1 ]] && VARS+=("S0T1=$E/S0_local"); [[ $HAVE_T0 == 1 ]] && VARS+=("S0T0=$V/T0/S0_local")
[[ $HK == 1 ]] && VARS+=("K2A=$V/K2A_prod"); [[ $HI == 1 ]] && VARS+=("K2Aimp=$V/K2A_imposed")
LU=1
if [[ ${#VARS[@]} -ge 2 ]]; then
  rc=0; python3 "$HERE/k2a_label_union.py" "$V/shaped_layout/muon_$M" --window-ns "$CMPW" "${VARS[@]}" > "$V/label_union.log" 2>&1 || rc=$?
  LU=$rc
  if [[ $rc -eq 0 ]]; then vstep label_union OK "$(tail -1 "$V/label_union.log" | head -c 300)" "$V/shaped_layout/muon_$M/label_union.json"
  else vstep label_union FAIL "$(tail -1 "$V/label_union.log" | head -c 300)" "$V/label_union.log"; fi
fi
if [[ $LU -eq 0 && $HAVE_T1 == 1 && $HK == 1 ]]; then
  X1="S0T1,K2A"; [[ $HI == 1 ]] && X1+=",K2Aimp"; X0="$X1"; [[ $HAVE_T0 == 1 ]] && X1+=",S0T0"
  shaped refS0T1_ionelec "$X1"
  shaped refS0T1_ion "$X1" --no-prompt
  if [[ $HAVE_T0 == 1 ]]; then shaped refS0T0_ionelec "S0T0,${X0#S0T1,},S0T1"; fi
  [[ $HI == 1 ]] && shaped refK2Aimp_ionelec "K2Aimp,K2A"
fi
exit 0
