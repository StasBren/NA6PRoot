#!/usr/bin/env bash
# Real TrackHeed timing with the frozen 2100 V reference executable, Heed only (--validate-heed: no avalanche, no drift).
# READ-ONLY. Uses the frozen 1000-crossing selection; rows whose archive index is in sealed_validation20.csv are dropped
# before anything runs (assumption, stated in the output: archive index muon_NNNN = 0-based data row NNNN of the
# selection, as written by the overnight runner).
#
#   bash run_heed_timing.sh [N_SMALL N_LARGE REPEATS]      default 50 400 3
#
# Per-crossing cost = slope of wall time between N_SMALL and N_LARGE crossings, so process start-up and the one-time
# Heed initialisation (cross-section tables on the first NewTrack) fall into the intercept. Each crossing builds a new
# TrackHeed + sensor (as in production HeedPrimaries), with delta-electron transport enabled.
set -Eeo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../../../.." && pwd)"
if [[ -f "$HOME/na6p/env.sh" ]];then source "$HOME/na6p/env.sh";fi
if [[ -f "$HOME/na6p/install/garfieldpp/share/Garfield/setupGarfield.sh" ]];then
  source "$HOME/na6p/install/garfieldpp/share/Garfield/setupGarfield.sh"; fi
set -u
umask 077
BIN="${MWPC_HEED_BIN:-$REPO/test/mwpc-avalanche/build-cathode-2100-wsl/mwpc_p3_noFW_penning057_cathode_ion_reference_v1}"
SEL="${MWPC_HEED_SELECTION:-$REPO/test/mwpc-stage-c/p3_wiretruth_1000_selection.csv}"
SELECTION="${MWPC_SIGNAL40_OUTPUT:-$HOME/NA6PRoot-signal40-2100/test_runs/mwpc_stage_c/phase_space_2100_signal40_v1}"
SEALED="$SELECTION/sealed_validation20.csv"
OUTROOT="${MWPC_HEED_TIMING_OUT:-$REPO/test_runs/mwpc_stage_c/heed_timing_2100}"
N1="${1:-50}"; N2="${2:-400}"; REP="${3:-3}"
fail(){ echo "HEED TIMING FAILED: $*" >&2; exit 2; }
[[ -x "$BIN" ]] || fail "binary $BIN missing (build target mwpc_p3_noFW_penning057_cathode_ion_reference_v1; set MWPC_HEED_BIN)"
[[ -s "$SEL" && -s "$SEALED" ]] || fail "need $SEL and $SEALED"
(( N1 >= 10 && N2 > N1 && N2 <= 900 )) || fail "need 10 <= N_SMALL < N_LARGE <= 900"
OUT="$OUTROOT/run_$(date +%Y%m%d_%H%M%S)"; mkdir -p "$OUT"
awk -F, 'NR==1{for(i=1;i<=NF;i++)if($i=="index")c=i; if(!c)exit 3; next} NF{print $c+0}' "$SEALED" | sort -n > "$OUT/sealed_rows.txt"
[[ $(wc -l < "$OUT/sealed_rows.txt") -eq 20 ]] || fail "expected 20 sealed indices"
# header + non-sealed data rows (0-based data row r = file line r+2)
awk 'NR==FNR{s[$1]=1; next} FNR==1{print; next} !((FNR-2) in s)' "$OUT/sealed_rows.txt" "$SEL" > "$OUT/nonsealed.csv"
head -n $((N1+1)) "$OUT/nonsealed.csv" > "$OUT/heed_${N1}.csv"; head -n $((N2+1)) "$OUT/nonsealed.csv" > "$OUT/heed_${N2}.csv"
{ echo "git_sha=$(git -C "$REPO" rev-parse HEAD)"; echo "binary=$BIN"; echo "binary_sha256=$(sha256sum "$BIN"|cut -d' ' -f1)"
  echo "selection_sha256=$(sha256sum "$SEL"|cut -d' ' -f1)"; echo "sealed_list_sha256=$(sha256sum "$SEALED"|cut -d' ' -f1)"
  echo "assumption=archive index NNNN == 0-based data row NNNN of the selection"
  echo "cpu=$(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2- | xargs)"; echo "nproc=$(nproc)"
  echo "garfield=${GARFIELD_INSTALL:-${GARFIELD_HOME:-unknown}}"; } > "$OUT/provenance.txt"
echo "n,rep,wall_s,user_s,sys_s,clusters,electrons" > "$OUT/timing.csv"
for n in "$N1" "$N2"; do for r in $(seq 1 "$REP"); do
  TIMEFORMAT='%R %U %S'   # bash builtin time: wall, user, sys (no /usr/bin/time needed)
  { time "$BIN" --validate-heed "$OUT/heed_${n}.csv" > "$OUT/log_${n}_${r}.txt" 2>&1; } 2> "$OUT/t_${n}_${r}.txt" \
    || fail "binary failed, see $OUT/log_${n}_${r}.txt"
  grep -q '^PASS: true TrackHeed' "$OUT/log_${n}_${r}.txt" || fail "no PASS line in $OUT/log_${n}_${r}.txt"
  read -r wall usr sys < "$OUT/t_${n}_${r}.txt"
  cl=$(sed -n 's/.*generated \([0-9]*\) clusters and \([0-9]*\) low-energy.*/\1/p' "$OUT/log_${n}_${r}.txt")
  el=$(sed -n 's/.*generated \([0-9]*\) clusters and \([0-9]*\) low-energy.*/\2/p' "$OUT/log_${n}_${r}.txt")
  echo "$n,$r,$wall,$usr,$sys,$cl,$el" >> "$OUT/timing.csv"; echo "n=$n rep=$r wall=${wall}s user=${usr}s"
done; done
python3 - "$OUT/timing.csv" "$N1" "$N2" <<'EOF' | tee "$OUT/heed_timing_summary.txt"
import csv, statistics as st, sys
rows = list(csv.DictReader(open(sys.argv[1]))); n1, n2 = int(sys.argv[2]), int(sys.argv[3])
def med(n, k): return st.median(float(r[k]) for r in rows if int(r["n"]) == n)
def get(n, k): return [float(r[k]) for r in rows if int(r["n"]) == n]
for k in ("wall_s", "user_s"):
    a, b = med(n1, k), med(n2, k)
    s = (b - a) / (n2 - n1)
    lo = (min(get(n2, k)) - max(get(n1, k))) / (n2 - n1); hi = (max(get(n2, k)) - min(get(n1, k))) / (n2 - n1)
    print(f"{k}: median {a:.3f} s ({n1} crossings), {b:.3f} s ({n2}); per-crossing slope {1e3*s:.3f} ms [range {1e3*lo:.3f}, {1e3*hi:.3f}]; intercept {a - s*n1:.3f} s")
e1, e2 = med(n1, "electrons"), med(n2, "electrons"); c1, c2 = med(n1, "clusters"), med(n2, "clusters")
s = (med(n2, "wall_s") - med(n1, "wall_s")) / (n2 - n1)
pe = (e2 - e1) / (n2 - n1)
print(f"Heed primaries per crossing (marginal) {pe:.2f}, clusters {(c2 - c1)/(n2 - n1):.2f}; wall per primary electron {1e6*s/pe:.2f} us")
print("measured: real TrackHeed, frozen 2100 V executable, single process; includes per-crossing TrackHeed/sensor construction")
EOF
echo "OUTPUT $OUT (return heed_timing_summary.txt, timing.csv, provenance.txt)"
