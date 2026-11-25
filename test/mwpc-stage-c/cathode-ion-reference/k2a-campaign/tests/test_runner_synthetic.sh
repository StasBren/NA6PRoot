#!/usr/bin/env bash
# SYNTHETIC end-to-end test of run_k2a_campaign.sh preflight (inventory -> chunked keep-going conversion -> domain),
# with a fake converter (shell script) on the fake archive of test_preflight_synthetic.py. No ROOT, no real data.
#   bash tests/test_runner_synthetic.sh WORK_OF_test_preflight_synthetic NEW_WORK
set -Eeuo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; S="${1:?}"; T="${2:?}"
[[ ! -e "$T" ]] || { echo "$T exists"; exit 2; }; mkdir -p "$T"
cp -r "$S/archive" "$S/selection" "$T/"
(cd "$T/archive/muon_0025" && sha256sum ion_births.root microscopic_wire_muon.csv > SHA256SUMS)   # repair the corrupted slot
mkdir -p "$T/run/kernel_tt"; for m in 0000 0001; do mkdir -p "$T/run/muon_$m/S0/replay"; : > "$T/run/muon_$m/S0/replay/waveforms.csv"; done
cp "$S/kernel/kernel.json" "$T/run/kernel_tt/"
python3 -c "import numpy as np; np.save('$T/run/kernel_tt/K.npy', np.zeros(3, np.float32)); np.save('$T/run/kernel_tt/KPHI0.npy', np.zeros(3))"
cat > "$T/fakebin" <<'FB'
#!/usr/bin/env bash
# fake p3_birth_moments: --keep-going n_primaries_zero_ions
[[ "$1" == --batch ]] || exit 9
mkdir "$3"; n=0
while read -r id path; do
  if [[ "$id" == muon_0013 ]]; then echo "synthetic failure" > "$3/$id.error.txt"; continue; fi
  printf 'primary,wire,t_ns,t_rms_ns,N_fC,r_eff_um,D_U_mm,D_V_mm,C1,S1,C2,S2,A,min_clearance_um,base_U,base_V,t_raw_rms_ns,rows\n0,0,30,3,10,30,0.1,0.2,0,0,0,0,0,5,1,2,3,100\n' > "$3/$id.moments.csv"
  s=${id#muon_}; s=$((10#$s)); ch=(Jpsi phi omega)
  echo "{\"source_channel\": \"${ch[$((s%3))]}\", \"station\": $((s%6)), \"det_id\": $((2*(s%6))), \"pdg\": 13, \"n_primaries\": 10, \"ion_births_exact\": 1000, \"rows_read\": 100}" > "$3/$id.moments.csv.event.json"
  n=$((n+1))
done < "$2"
[[ -n "$(ls "$3"/*.error.txt 2>/dev/null)" ]] && { echo "BATCH FAILURES"; exit 3; }
echo "PASS birth moments batch events=$n"
FB
chmod +x "$T/fakebin"
MWPC_K2A_SLOTS=30 MWPC_ARCHIVE_2100_K256="$T/archive" MWPC_SIGNAL40_OUTPUT="$T/selection" MWPC_K2A_REFERENCE_RUN="$T/run" \
MWPC_K2A_TAG=tt MWPC_BIRTHS_BIN="$T/fakebin" MWPC_K2A_OUT="$T/out" MWPC_K2A_CHUNK=7 bash "$HERE/run_k2a_campaign.sh" preflight > "$T/preflight.log"
P=$(ls -d "$T"/out/preflight_*)
grep -q "| all | 24 | 23 " "$P/domain/coverage.md"
grep -q "muon_0013 | pool | phi | 1 | READ_FAILED | READ_FAILED: synthetic failure" "$P/domain/coverage.md"
grep -q "rc=3" "$P/moments_nonsealed/chunks.txt"
! grep -q SEALED_DO_NOT_OPEN "$P/PREFLIGHT_SHA256SUMS"
ENVS=(MWPC_K2A_SLOTS=30 MWPC_ARCHIVE_2100_K256="$T/archive" MWPC_SIGNAL40_OUTPUT="$T/selection" MWPC_K2A_REFERENCE_RUN="$T/run" MWPC_K2A_TAG=tt MWPC_K2A_OUT="$T/out")
env "${ENVS[@]}" bash "$HERE/run_k2a_campaign.sh" recover "$P" > "$T/recover.log"
RD=$(ls -d "$P"/recenter_*/ | head -1); grep -q "## 2. Coverage recovery" "$RD/recovery.md"
env "${ENVS[@]}" bash "$HERE/run_k2a_campaign.sh" supplementary "$P" "$RD" 3 > "$T/supp.log"
SD=$(ls -d "$P"/supplementary_*/ | head -1); (cd "$SD" && sha256sum -c --quiet SHA256SUMS)
! grep -q -E "^(0003|0004|0011)," "$SD/supplementary_dev.csv"
echo "PASS runner recover + supplementary (synthetic): recovery report, selection excludes the sealed slots"
echo "PASS runner preflight (synthetic): chunked keep-going conversion, failure visible, sealed blinded, checksums written"
