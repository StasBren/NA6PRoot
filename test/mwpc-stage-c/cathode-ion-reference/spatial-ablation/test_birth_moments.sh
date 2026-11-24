#!/usr/bin/env bash
# SYNTHETIC test of the birth-record -> moments path (no ROOT): BirthMoments.hh core vs the established
# export -> T1 collapse -> s_kernel.py moments path, then kernel waveforms from both moment files.
#   bash test_birth_moments.sh WORKDIR [LUT.npz NUMW_DIR]
# With LUT and NUMW_DIR it also builds a small SYNTHETIC-transport Green table + kernels covering the births and checks
# that the K2A kernel waveforms from the two moment files agree to full precision (ion term and electron step).
set -Eeuo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; TA="$HERE/../timing-ablation"
W="${1:-$(mktemp -d)}"; mkdir -p "$W"; [[ -z "$(ls -A "$W")" ]] || { echo "workdir not empty" >&2; exit 2; }
g++ -O2 -std=c++17 -Wall -Wextra -o "$W/core" "$HERE/tests/birth_moments_core_test.cc"
python3 "$HERE/tests/make_synthetic_births.py" "$W/syn" 7
python3 "$TA/t_time_collapse_export.py" "$W/syn/frozen_export" "$W/syn/T1_export" --mode primary > /dev/null
python3 "$HERE/s_kernel.py" moments "$W/syn/T1_export" "$W/mom_export.csv"
"$W/core" "$W/syn/births.csv" "$W/syn/frame.csv" "$W/mom_births.csv"
python3 "$HERE/s_compare_moments.py" "$W/mom_export.csv" "$W/mom_births.csv" --tol 1e-9 --report "$W/compare_moments.md" | tail -16
grep -q '^PASS cloud moments agree' "$W/compare_moments.md"
if [[ $# -ge 3 ]]; then
  LUT="$2"; NUMW="$3"
  python3 "$HERE/s_grid_spec.py" "$W/syn/T1_export" --nr 8 > "$W/spec.txt"
  RU=$(sed -n 's/.*--r-um \([^ ]*\)$/\1/p' "$W/spec.txt"); WI=$(sed -n 's/^--wires \([^ ]*\) --r-um.*/\1/p' "$W/spec.txt")
  DM=$(sed -n 's/.*--d-max-mm \([0-9.]*\)$/\1/p' "$W/spec.txt")
  python3 "$HERE/synthetic_paths.py" grid "$W/grid" --wires="$WI" --r-um="$RU" --ntheta 16 > /dev/null
  python3 "$HERE/s_green_table.py" "$LUT" "$W/grid" "$W/green" --replay-dir "$NUMW" --d-max-mm "$DM" --ntau 120 | tail -1
  python3 "$HERE/s_kernel.py" build "$W/green" "$W/kern"
  # labels for apply: any replay-format file with the strip frame of this export (only the header is used)
  python3 - "$W" <<'PY'
import sys, csv
W = sys.argv[1]; r = next(csv.DictReader(open(f"{W}/syn/T1_export/ions.csv")))
bu, bv = int(r["upper_index"]), int(r["lower_index"])
lab = [f"U_{bu + o}" for o in range(-5, 6)] + [f"V_{bv + o}" for o in range(-5, 6)]
open(f"{W}/labels.csv", "w").write("time_center_ns," + ",".join(lab) + "\n" + ",".join(["5"] + ["0"] * len(lab)) + "\n")
PY
  mkdir -p "$W/kx" "$W/kb"
  python3 "$HERE/s_kernel.py" apply "$W/kern" "$W/mom_export.csv" "$W/labels.csv" "$W/kx" --model K2A --repeat 1
  python3 "$HERE/s_kernel.py" apply "$W/kern" "$W/mom_births.csv" "$W/labels.csv" "$W/kb" --model K2A --repeat 1
  python3 "$HERE/s_compare_waveforms.py" "$W/kx/K2Akern" "$W/kb/K2Akern" --tol 1e-9 --report "$W/compare_waveforms.md" | tail -8
fi
echo "PASS birth-moments synthetic test; work dir $W"
