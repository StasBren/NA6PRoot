#!/usr/bin/env bash
#PBS -N cSig21
#PBS -m n
#PBS -l select=1:ncpus=1:mem=32gb
#PBS -l walltime=12:00:00
#PBS -j oe
# Wrapper: map 0/1 pilot array indices to actual frozen 2100V Geant4 indices.
# Physics worker is the already-reviewed separate uncapped ROOT-ion target.
set -Eeuo pipefail
: "${MWPC_SIGNAL40_SELECTION:?Set sealed selected40.csv}"
: "${MWPC_SIGNAL40_SELECTION_SHA:?Set frozen SHA256 of selection}"
: "${MWPC_DETAIL_REPO:?Set isolated physics clone}"
: "${PBS_ARRAY_INDEX:?PBS array job required}"
[[ "$PBS_ARRAY_INDEX" =~ ^[01]$ ]] || {
  echo "Only first two development-muon pilot indices 0 and 1 allowed" >&2
  exit 3
}
[[ "$(sha256sum "$MWPC_SIGNAL40_SELECTION" | cut -d' ' -f1)" ==
    "$MWPC_SIGNAL40_SELECTION_SHA" ]] || {
  echo "Selected40 manifest SHA changed after submission" >&2;exit 3;
}
# A regular CSV row is never modified: only map the existing explicit index.
# Exclude headers, validation, and the previously inspected 0301 event.
real_index="$(awk -F, '$1=="development" {print $2}' \
  "$MWPC_SIGNAL40_SELECTION" | sed -n "$((PBS_ARRAY_INDEX+1))p")"
[[ "$real_index" =~ ^[0-9]+$ ]] &&
  (( real_index >=0 && real_index <1000 && real_index != 301 )) || {
    echo "Invalid/previously-seen selected real Geant4 index" >&2;exit 3;
  }
echo "SELECTED_2100V_SIGNAL40_PILOT task=$PBS_ARRAY_INDEX original_muon=$real_index"
export PBS_ARRAY_INDEX="$real_index"
exec bash "$MWPC_DETAIL_REPO/test/mwpc-stage-c/cathode-ion-reference/pbs_cathode_ion_reference.sh"
