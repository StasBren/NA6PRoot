#!/usr/bin/env bash
# Resolve (read-only) the inputs of run_k2a_validation.sh on the existing WSL setup and write them to a NEW env file
# for review before sourcing. Nothing is guessed silently: every path is checked, ambiguities are listed in the file.
#
#   source ~/t_ladder_deps.env
#   bash k2a_validation_env.sh ENV_FILE
#
# Resolves (overridable by exporting the variable beforehand):
#   K2AV_RUN     the reference run .../test_runs/mwpc_stage_c/s_ladder_2100/run_20261124_124740 with kernel_r14t72edge5
#   K2AV_TAG     r14t72edge5
#   MWPC_K2A_PF  newest .../k2a_campaign_2100/preflight_* that has inventory/slots.csv, domain/per_event.csv,
#                moments_nonsealed/ and a recenter_*/per_event_recenter.csv (the recovery output)
#   K2AV_REC     newest PF/recenter_*/ containing per_event_recenter.csv
#   MWPC_BIRTHS_BIN        ~/na6p/build-spatial/p3_birth_moments (the spatial-ablation build), else the only executable
#                          p3_birth_moments under ~/na6p/build*/; its usage text must offer --batch and --repeat
#   MWPC_K2A_INLINE_BUILD  ~/na6p/build-k2a-inline (must contain k2a_strip_eval)
#   MWPC_ARCHIVE_2100_K256 the frozen archive default of k2a_recenter_experiment.sh (must contain muon_0000)
# Searched roots: ~/na6p/test_runs and ~/na6p/*/test_runs. Exit 0 file written; 2 something missing (named).
set -Eeuo pipefail
OUTF="${1:?usage: k2a_validation_env.sh ENV_FILE}"
[[ ! -e "$OUTF" ]] || { echo "refusing: $OUTF exists" >&2; exit 2; }
miss(){ echo "NOT RESOLVED: $*" >&2; exit 2; }
for v in MWPC_SPECTRAL_LUT MWPC_NUMW_DIR MWPC_EXPORT_BIN MWPC_ION_MOBILITY; do
  [[ -n "${!v:-}" && -e "${!v}" ]] || miss "$v (source ~/t_ladder_deps.env first)"
done
shopt -s nullglob
NOTES=()
TAG="${K2AV_TAG:-r14t72edge5}"
if [[ -z "${K2AV_RUN:-}" ]]; then
  C=(); for d in "$HOME"/na6p/test_runs/mwpc_stage_c/s_ladder_2100/run_20261124_124740 "$HOME"/na6p/*/test_runs/mwpc_stage_c/s_ladder_2100/run_20261124_124740; do
    [[ -s "$d/kernel_$TAG/kernel.json" && -s "$d/kernel_$TAG/K.npy" ]] && C+=("$(cd "$d" && pwd -P)"); done
  [[ ${#C[@]} -eq 0 ]] || mapfile -t C < <(printf '%s\n' "${C[@]}" | sort -u)
  [[ ${#C[@]} -ge 1 ]] || miss "reference run run_20261124_124740 with kernel_$TAG under ~/na6p/{,*/}test_runs (export K2AV_RUN)"
  [[ ${#C[@]} -eq 1 ]] || miss "several reference runs: ${C[*]} (export K2AV_RUN to choose)"
  K2AV_RUN="${C[0]}"
fi
[[ -s "$K2AV_RUN/kernel_$TAG/kernel.json" ]] || miss "kernel_$TAG in $K2AV_RUN"
if [[ -z "${MWPC_K2A_PF:-}" ]]; then
  C=(); for d in "$HOME"/na6p/test_runs/mwpc_stage_c/k2a_campaign_2100/preflight_* "$HOME"/na6p/*/test_runs/mwpc_stage_c/k2a_campaign_2100/preflight_*; do
    r=("$d"/recenter_*/per_event_recenter.csv)
    [[ -s "$d/inventory/slots.csv" && -s "$d/domain/per_event.csv" && -d "$d/moments_nonsealed" && ${#r[@]} -ge 1 ]] && C+=("$(cd "$d" && pwd -P)")
  done
  [[ ${#C[@]} -eq 0 ]] || mapfile -t C < <(printf '%s\n' "${C[@]}" | awk -F/ '{print $NF"\t"$0}' | sort | cut -f2 | uniq)
  [[ ${#C[@]} -ge 1 ]] || miss "no complete preflight (inventory, domain, moments_nonsealed, recenter_*/per_event_recenter.csv) under ~/na6p/{,*/}test_runs/mwpc_stage_c/k2a_campaign_2100 (export MWPC_K2A_PF)"
  MWPC_K2A_PF="${C[-1]}"
  [[ ${#C[@]} -eq 1 ]] || NOTES+=("other complete preflights (not used): ${C[*]:0:${#C[@]}-1}")
fi
if [[ -z "${K2AV_REC:-}" ]]; then
  R=(); for d in "$MWPC_K2A_PF"/recenter_*/; do [[ -s "$d/per_event_recenter.csv" ]] && R+=("${d%/}"); done
  [[ ${#R[@]} -ge 1 ]] || miss "no recenter_*/per_event_recenter.csv in $MWPC_K2A_PF (export K2AV_REC)"
  K2AV_REC="${R[-1]}"
  [[ ${#R[@]} -eq 1 ]] || NOTES+=("other recovery outputs in the preflight (not used): ${R[*]:0:${#R[@]}-1}")
fi
if [[ -z "${MWPC_BIRTHS_BIN:-}" ]]; then
  if [[ -x "$HOME/na6p/build-spatial/p3_birth_moments" ]]; then MWPC_BIRTHS_BIN="$HOME/na6p/build-spatial/p3_birth_moments"
  else
    C=("$HOME"/na6p/build*/p3_birth_moments); X=(); for f in "${C[@]}"; do [[ -x "$f" ]] && X+=("$f"); done
    [[ ${#X[@]} -eq 1 ]] || miss "p3_birth_moments: ${#X[@]} executables under ~/na6p/build*/ (${X[*]:-none}); export MWPC_BIRTHS_BIN"
    MWPC_BIRTHS_BIN="${X[0]}"
  fi
fi
u=$("$MWPC_BIRTHS_BIN" 2>&1 || true)
[[ $u == *"--batch"* && $u == *"--repeat"* ]] || miss "$MWPC_BIRTHS_BIN does not start or is not p3_birth_moments with --batch/--repeat (ROOT/Garfield runtime?): ${u:0:200}"
MWPC_K2A_INLINE_BUILD="${MWPC_K2A_INLINE_BUILD:-$HOME/na6p/build-k2a-inline}"
[[ -x "$MWPC_K2A_INLINE_BUILD/k2a_strip_eval" ]] || miss "k2a_strip_eval in $MWPC_K2A_INLINE_BUILD"
ARCH="${MWPC_ARCHIVE_2100_K256:-$HOME/na6p/cathode-overnight-2100/test_runs/mwpc_stage_c/overnight_2100_weighted_k256_v1}"
[[ -s "$ARCH/muon_0000/ion_births.root" ]] || miss "archive $ARCH (export MWPC_ARCHIVE_2100_K256)"
umask 077
{ echo "# written by k2a_validation_env.sh $(date -u +%FT%TZ); REVIEW before sourcing. Read-only inputs of the K2A validation."
  for n in "${NOTES[@]}"; do echo "# NOTE: $n"; done
  echo "export K2AV_RUN='$K2AV_RUN'"; echo "export K2AV_TAG='$TAG'"
  echo "export MWPC_K2A_PF='$MWPC_K2A_PF'"; echo "export K2AV_REC='$K2AV_REC'"
  echo "export MWPC_BIRTHS_BIN='$MWPC_BIRTHS_BIN'   # sha256 $(sha256sum "$MWPC_BIRTHS_BIN" | cut -d' ' -f1)"
  echo "export MWPC_K2A_INLINE_BUILD='$MWPC_K2A_INLINE_BUILD'   # k2a_strip_eval sha256 $(sha256sum "$MWPC_K2A_INLINE_BUILD/k2a_strip_eval" | cut -d' ' -f1)"
  echo "export MWPC_ARCHIVE_2100_K256='$ARCH'"; } > "$OUTF"
cat "$OUTF"
