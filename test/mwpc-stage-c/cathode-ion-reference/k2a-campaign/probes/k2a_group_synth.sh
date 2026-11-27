#!/usr/bin/env bash
# Per-group S-ladder synthesis for ONE validation muon (e.g. 0637), from its EXISTING group exports and the frozen Green
# table: no avalanche, no ion transport, no S0 rerun. Each variant is a superposition of exact single-ion Green
# responses with one spatial simplification (s_synthesize.py, unchanged); the group's existing S0 (frozen replay at the
# chosen late step) is the reference. Feed the layout to k2a_closure_anatomy.py --layout.
#
#   bash k2a_group_synth.sh MUON_VALIDATION_DIR RUN TAG NEW_OUT [VARIANTS]
#
# MUON_VALIDATION_DIR  OUT/muons/muon_M of run_k2a_validation.sh (far_dir.txt -> E/groups/g_*, E/S0groups/g_*)
# RUN, TAG             the reference run and tag holding green_TAG (the Green table the kernel was built from)
# VARIANTS             default G,SR,S3,S3A,SW,K2A
# The group exports must still hold vertices.f64 (not pruned). Cost: Python, per ion; scales with the group's rows.
set -Eeo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; CIR="$(cd "$HERE/../.." && pwd)"; SA="$CIR/spatial-ablation"
set -u
V="${1:?usage: k2a_group_synth.sh MUON_VALIDATION_DIR RUN TAG NEW_OUT [VARIANTS]}"; RUN="${2:?RUN}"; TAG="${3:?TAG}"; OUT="${4:?NEW_OUT}"
VARS="${5:-G,SR,S3,S3A,SW,K2A}"
[[ ! -e "$OUT" ]] || { echo "refusing: $OUT exists" >&2; exit 2; }
GD="$RUN/green_$TAG"; [[ -s "$GD/axes.json" && -s "$GD/G.npy" ]] || { echo "no Green table $GD" >&2; exit 2; }
E=$(cat "$V/far_dir.txt")
umask 077; mkdir -p "$OUT"
for g in "$E"/groups/g_*; do
  n=${g##*/}; S0="$E/S0groups/$n"
  [[ -s "$g/vertices.f64" ]] || { echo "$g/vertices.f64 missing (pruned?)" >&2; exit 2; }
  [[ -s "$S0/replay/waveforms.csv" ]] || { echo "no S0 for $n" >&2; exit 2; }
  mkdir -p "$OUT/$n"; ln -s "$S0" "$OUT/$n/S0"
  ( time python3 "$SA/s_synthesize.py" "$GD" "$g" "$S0/replay/waveforms.csv" "$OUT/$n" --variants "$VARS" --suffix "_$TAG" ) \
      > "$OUT/$n/synth.log" 2>&1 || { echo "s_synthesize failed for $n ($OUT/$n/synth.log)" >&2; exit 1; }
  echo "  $n: $(grep -c . "$g/ions.csv") ion rows synthesized ($(tail -3 "$OUT/$n/synth.log" | grep real || true))"
done
echo "layout ready: python3 $HERE/k2a_closure_anatomy.py NEW_DIAG --layout $OUT --tag $TAG"
