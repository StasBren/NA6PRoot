#!/usr/bin/env bash
# Compare Heed primary ionisation for symmetric 5 mm and 6 mm gas gaps.
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
NTRACKS=20000
if [ "$#" -ge 1 ]; then NTRACKS="$1"; fi
if ! [[ "$NTRACKS" =~ ^[1-9][0-9]*$ ]]; then
  echo "Usage: bash run_phaseA_heed_gap_scan.sh [positive_tracks] [output_directory]" >&2
  exit 2
fi
if [ "$#" -gt 2 ]; then
  echo "Usage: bash run_phaseA_heed_gap_scan.sh [positive_tracks] [output_directory]" >&2
  exit 2
fi

if [ "$#" -eq 2 ]; then
  OUTDIR="$2"
else
  OUTDIR="$HERE/runs/heed_gap_scan_$(date +%Y%m%d_%H%M%S)"
fi
mkdir -p "$OUTDIR"
OUTDIR="$(cd "$OUTDIR" && pwd)"

EXE="$HERE/build/mwpc_phase_a_heed_stats"
if [ ! -x "$EXE" ]; then
  echo "Missing $EXE. Source Garfield, then run:" >&2
  echo "  cmake -S \"$HERE\" -B \"$HERE/build\" -DGarfield_DIR=\"\$HOME/na6p/install/garfieldpp/lib/cmake/Garfield\"" >&2
  echo "  cmake --build \"$HERE/build\" -j4 --target mwpc_phase_a_heed_stats" >&2
  exit 1
fi

# Heed-only scan. Fixed gas, particle, angle, wire pitch and HV.
# The chamber is symmetric about the wire plane in both cases.
# Note phaseA_heed_statistics.cxx starts 0.1 mm inside the +v cathode;
# effective nominal ionisation path lengths are about 4.9 and 5.9 mm.
for item in "5:2.5" "6:3.0"; do
  IFS=: read -r GAP HALF <<< "$item"
  PREFIX="$(printf '%s/heed_%smm' "$OUTDIR" "$GAP")"
  echo
  echo "=== Heed gas thickness $GAP mm ($HALF+$HALF mm) ==="
  "$EXE" \
    --tracks "$NTRACKS" \
    --momentum-gev 10 \
    --theta-deg 0 \
    --phi-deg 0 \
    --u0-mm 1 \
    --gap-minus-mm "$HALF" \
    --gap-plus-mm "$HALF" \
    --pitch-mm 4 \
    --wire-diam-um 30 \
    --hv 1800 \
    --output-prefix "$PREFIX" \
    > "$PREFIX.log" 2>&1
  tail -n 16 "$PREFIX.log"
done

python3 "$HERE/analyze_phaseA_heed_gap_scan.py" \
  --five "$OUTDIR/heed_5mm_tracks.csv" \
  --six "$OUTDIR/heed_6mm_tracks.csv" \
  --output-dir "$OUTDIR"

echo
echo "Results: $OUTDIR"
echo "See summary.csv, low_tail.csv, comparison.txt and (if matplotlib is installed) PNG plots."
