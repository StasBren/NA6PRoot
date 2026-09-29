#!/usr/bin/env bash
set -euo pipefail

# Generate three presentation GIFs from the same detector configuration:
#   1) normal incidence close to one wire;
#   2) normal incidence at the midpoint between two wires;
#   3) a 20-degree muon tilted in the local u-v plane.
#
# The electric-field streamlines are seeded only in a narrow band around the
# actual muon / Heed-cluster region. The same renderer therefore adapts to all
# three geometries without filling the whole chamber with irrelevant lines.

EXE="./build/mwpc_heed_visualization"
RENDERER="render_heed_gif.py"

if [[ ! -x "$EXE" ]]; then
  echo "Missing $EXE. Build the mwpc-avalanche targets first."
  exit 1
fi

COMMON=(
  --pitch-mm 4
  --gap-minus-mm 2
  --gap-plus-mm 4
  --wire-diam-um 30
  --hv 1800
  --b 0
  --momentum-gev 10
  --sample-mode balanced
  --max-electrons 60
  --max-per-cluster 8
  --rng-seed 12345
  --max-path-points 250
  --field-u-steps 721
  --field-v-steps 361
)

RENDER_COMMON=(
  --background stream
  --v-min-mm 0.0
  --v-max-mm 4.1
  --stream-density 0.70
  --stream-focus track
  --stream-seeds 11
  --stream-band-mm 0.22
  --muon-frames 18
  --cluster-frames 24
  --drift-frames 70
  --hold-frames 18
  --fps 18
)

echo
echo "=== Near-wire event: u0 = 0.4 mm ==="
"$EXE" "${COMMON[@]}" \
  --u0-mm 0.4 \
  --output-prefix heed_nearwire

python3 "$RENDERER" \
  --prefix heed_nearwire \
  --output heed_nearwire.gif \
  --u-min-mm -0.8 \
  --u-max-mm 1.6 \
  "${RENDER_COMMON[@]}"

echo
echo "=== Midpoint event: u0 = 2.0 mm ==="
"$EXE" "${COMMON[@]}" \
  --u0-mm 2.0 \
  --output-prefix heed_midpoint

python3 "$RENDERER" \
  --prefix heed_midpoint \
  --output heed_midpoint.gif \
  --u-min-mm -0.5 \
  --u-max-mm 4.5 \
  "${RENDER_COMMON[@]}"

echo
echo "=== Angled event: theta = 20 deg toward +u, u0 = 0.4 mm ==="
"$EXE" "${COMMON[@]}" \
  --u0-mm 0.4 \
  --theta-deg 20 \
  --phi-deg 0 \
  --output-prefix heed_angle20

python3 "$RENDERER" \
  --prefix heed_angle20 \
  --output heed_angle20.gif \
  --u-min-mm -0.5 \
  --u-max-mm 4.5 \
  "${RENDER_COMMON[@]}"

echo
echo "Created:"
echo "  heed_nearwire.gif"
echo "  heed_midpoint.gif"
echo "  heed_angle20.gif"
