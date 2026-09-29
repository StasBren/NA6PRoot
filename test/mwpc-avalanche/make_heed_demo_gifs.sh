#!/usr/bin/env bash
set -euo pipefail

# Generate two presentation GIFs from the same detector configuration:
#   1) a muon close to one wire;
#   2) a muon exactly at the midpoint between two wires.
#
# These are illustrative single events, not a controlled event-by-event
# comparison of identical microscopic ionisation histories.

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
  --max-electrons 40
  --max-per-cluster 3
  --rng-seed 12345
  --max-path-points 250
)

echo
echo "=== Near-wire event: u0 = 0.4 mm, visible wires 0 and 4 mm ==="
"$EXE" "${COMMON[@]}"   --u0-mm 0.4   --output-prefix heed_nearwire

python3 "$RENDERER"   --prefix heed_nearwire   --output heed_nearwire.gif

echo
echo "=== Midpoint event: u0 = 2.0 mm, visible wires 0 and 4 mm ==="
"$EXE" "${COMMON[@]}"   --u0-mm 2.0   --output-prefix heed_midpoint

python3 "$RENDERER"   --prefix heed_midpoint   --output heed_midpoint.gif

echo
echo "Created:"
echo "  heed_nearwire.gif"
echo "  heed_midpoint.gif"
