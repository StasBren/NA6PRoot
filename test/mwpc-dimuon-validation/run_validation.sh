#!/usr/bin/env bash
set -euo pipefail

NEVENTS="${1:-20000}"
WORKERS="${2:-4}"

ROOT_DIR="$(git rev-parse --show-toplevel)"
GEN="${ROOT_DIR}/test/genDimuonBgEvent.C"
HOOK="${ROOT_DIR}/test/mwpc-dimuon-validation/mwpcDimuonHooks.C"
KIN_ANALYSIS="${ROOT_DIR}/test/mwpc-dimuon-validation/analyzeMWPCDimuon.C"
HIT_ANALYSIS="${ROOT_DIR}/test/mwpc-dimuon-validation/plotMWPCHitDensity.C"
ANGLE_ANALYSIS="${ROOT_DIR}/test/mwpc-dimuon-validation/analyzeMWPCIncidentAngles.C"
OUT_BASE="${ROOT_DIR}/test_runs/mwpc_dimuon"

mkdir -p "${OUT_BASE}"

if ! command -v na6psim >/dev/null 2>&1; then
  echo "na6psim is not in PATH. Source the NA6PRoot environment first." >&2
  exit 1
fi
if ! command -v na6psim_parallel >/dev/null 2>&1; then
  echo "na6psim_parallel is not in PATH. Source the current NA6PRoot install first." >&2
  exit 1
fi

# Compile the generator and hook once before parallel workers start. This avoids
# several ACLiC processes trying to build the same macro at the same time.
# NA6PRoot canonicalizes keyval.output_dir during configuration loading, so the
# directory must already exist for this direct na6psim smoke run.
SMOKE_DIR="${OUT_BASE}/_compile_smoke"
rm -rf "${SMOKE_DIR}"
mkdir -p "${SMOKE_DIR}"
na6psim \
  -n 1 \
  -r 20260914 \
  -g "${GEN}+(1,\"Jpsi\",0.,5.,0.,6.,false)" \
  -u "${HOOK}+" \
  --doDigitization false \
  --configKeyValues "beam.energyPerNucleon=40;keyval.output_dir=${SMOKE_DIR}"
rm -rf "${SMOKE_DIR}"

run_channel() {
  local channel="$1"
  local seed="$2"
  local out="${OUT_BASE}/${channel}"

  echo
  echo "============================================================"
  echo "Running ${channel}: ${NEVENTS} events, ${WORKERS} workers"
  echo "Output: ${out}"
  echo "============================================================"

  rm -rf "${out}"
  na6psim_parallel \
    --workers "${WORKERS}" \
    -n "${NEVENTS}" \
    -r "${seed}" \
    -g "${GEN}+(1,\"${channel}\",0.,5.,0.,6.,false)" \
    -u "${HOOK}+" \
    --doDigitization false \
    --configKeyValues "beam.energyPerNucleon=40;keyval.output_dir=${out}"

  # Parent/daughter momentum distributions and decay-closure checks.
  root -l -b -q "${KIN_ANALYSIS}+(\"${out}\",\"${channel}\")"

  # Direct Geant4 MWPC sensitive-gas hit density. Every chamber crossing is
  # counted separately at its actual position; 1 cm bins resolve the 3 cm seams.
  root -l -b -q "${HIT_ANALYSIS}+(\"${out}\",\"${channel}\",1.)"

  # Stage C0: incidence-angle phase space at MS0-MS3. The analysis uses the
  # entrance momentum stored in each NA6PMuonSpecHit and reports both the
  # current analytic 320x240 cm MNP33 aperture and the conservative 245x240 cm
  # useful-aperture fiducial.
  root -l -b -q "${ANGLE_ANALYSIS}+(\"${out}\",\"${channel}\")"
}

run_channel Jpsi 20260915
run_channel Omega 20260916
run_channel Phi 20260917

echo
echo "All channels finished."
echo "Kinematics: ${OUT_BASE}/{Jpsi,Omega,Phi}/plots/"
echo "Hit density: ${OUT_BASE}/{Jpsi,Omega,Phi}/plots_direct_hits/"
echo "Incident angles: ${OUT_BASE}/{Jpsi,Omega,Phi}/plots_angles/"
