#!/usr/bin/env bash
set -euo pipefail

EXE="./build/mwpc_phase_a_avalanche_visualization"
RENDERER="render_avalanche_gif.py"

if [[ ! -x "$EXE" ]]; then
  echo "Missing $EXE. Build test/mwpc-avalanche first."
  exit 1
fi

COMMON=(
  --pitch-mm 4
  --gap-minus-mm 2
  --gap-plus-mm 4
  --wire-diam-um 30
  --hv 1800
  --b 0
  --energy-ev 0.1
  --half-wires 4
  --avalanche-limit 5000
  --max-lines 1500
  --max-points-per-line 120
  --field-half-width-mm 3.0
  --field-u-steps 181
  --field-v-steps 181
)

echo
echo "=== Centered seed above the anode wire ==="
"$EXE" "${COMMON[@]}"   --u-mm 0.0   --v-mm 1.0   --w-mm 0.0   --output-prefix avalanche_centered

python3 "$RENDERER"   --prefix avalanche_centered   --output avalanche_centered_contours.gif   --background contours   --half-width-mm 1.25   --frames 80   --fps 18

echo
echo "=== Off-axis seed approaching the same wire ==="
"$EXE" "${COMMON[@]}"   --u-mm 0.8   --v-mm 1.0   --w-mm 0.0   --output-prefix avalanche_offaxis

python3 "$RENDERER"   --prefix avalanche_offaxis   --output avalanche_offaxis_contours.gif   --background contours   --half-width-mm 1.5   --frames 80   --fps 18

echo
echo "Created:"
echo "  avalanche_centered_contours.gif"
echo "  avalanche_offaxis_contours.gif"
echo
echo "Optional field-arrow rendering:"
echo "  python3 render_avalanche_gif.py --prefix avalanche_centered --output avalanche_centered_field.gif --background both --half-width-mm 1.25"
