#!/usr/bin/env bash
set -euo pipefail

# End-to-end PI-style fast MWPC response workflow.
#
# Uses:
#   1) realistic MNP33 crossing kinematics from c_pi_track_sample.csv;
#   2) Poisson cluster counts with density from the 2+4 mm Heed library;
#   3) empirical Heed cluster electron multiplicities;
#   4) empirical Garfield single-electron transport responses from the complete
#      2+4 mm transport-box map.
#
# The workflow stops before avalanche-gain fluctuations, cathode induction and
# electronics.
#
# Usage:
#   bash test/mwpc-stage-c/run_c_pi_fast_response.sh
#
# Optional environment:
#   C_PI_MAX_TRACKS=0        0 means all tracks in c_pi_track_sample.csv
#   C_PI_CLUSTER_SEED=180101
#   C_PI_TRANSPORT_SEED=181001
#   C_PI_REGENERATE_CLUSTERS=0  set 1 to force regeneration

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
STAGE_DIR="${ROOT_DIR}/test/mwpc-stage-c"

TRACK_DIR="${ROOT_DIR}/test_runs/mwpc_stage_c/c_pi_tracks"
HEED_DIR="${ROOT_DIR}/test_runs/mwpc_stage_c/c_heed_library_2p4"
TRANSPORT_DIR="${ROOT_DIR}/test_runs/mwpc_stage_c/c_transport_box_map_gm2p0_gp4p0"
OUT_DIR="${ROOT_DIR}/test_runs/mwpc_stage_c/c_pi_fast_response_2p4"

TRACKS="${TRACK_DIR}/c_pi_track_sample.csv"
HEED_TRACKS="${HEED_DIR}/heed_2p4_tracks.csv"
HEED_CLUSTERS="${HEED_DIR}/heed_2p4_clusters.csv"
TRANSPORT="${TRANSPORT_DIR}/c_transport_box_electrons.csv"

CLUSTERS="${OUT_DIR}/c_pi_clusters_heed2p4.csv"
PREFIX="${OUT_DIR}/c_pi_fast"
SUMMARY="${OUT_DIR}/c_pi_fast_summary.csv"
UBINS="${OUT_DIR}/c_pi_fast_u_bins.csv"

MAX_TRACKS="${C_PI_MAX_TRACKS:-0}"
CLUSTER_SEED="${C_PI_CLUSTER_SEED:-180101}"
TRANSPORT_SEED="${C_PI_TRANSPORT_SEED:-181001}"
REGENERATE_CLUSTERS="${C_PI_REGENERATE_CLUSTERS:-0}"

mkdir -p "${OUT_DIR}"

for f in "${TRACKS}" "${HEED_TRACKS}" "${HEED_CLUSTERS}" "${TRANSPORT}"; do
  if [[ ! -f "${f}" ]]; then
    echo "Missing required input: ${f}" >&2
    exit 2
  fi
done

# Guard against the historical one-column transport-map bug. The complete
# 2+4 mm map is 20 x 30 = 600 boxes.
python3 - "${TRANSPORT}" <<'PY'
import csv, sys
path = sys.argv[1]
boxes = set()
n = 0
with open(path, newline="") as f:
    for r in csv.DictReader(f):
        boxes.add((int(r["iu"]), int(r["iv"])))
        n += 1
if len(boxes) != 600:
    raise SystemExit(
        f"Transport map incomplete: found {len(boxes)} boxes, expected 600."
    )
print(f"Transport map check: {len(boxes)} boxes, {n} electron trials")
PY

echo
echo "=== 1/3 Prepare PI clusters from realistic tracks + Heed library ==="
if [[ -f "${CLUSTERS}" && "${REGENERATE_CLUSTERS}" != "1" ]]; then
  echo "Reusing existing cluster sample: ${CLUSTERS}"
else
  python3 "${STAGE_DIR}/generate_c_pi_poisson_clusters.py" \
    --tracks "${TRACKS}" \
    --output "${CLUSTERS}" \
    --gap-minus-mm 2 \
    --gap-plus-mm 4 \
    --heed-tracks-csv "${HEED_TRACKS}" \
    --heed-clusters-csv "${HEED_CLUSTERS}" \
    --heed-reference-length-mm 5.9 \
    --cluster-electron-model heed-empirical \
    --max-tracks "${MAX_TRACKS}" \
    --seed "${CLUSTER_SEED}"
fi

echo
echo "=== 2/3 Fold clusters through Garfield box responses ==="
python3 "${STAGE_DIR}/fold_c_pi_clusters_with_transport.py" \
  --clusters "${CLUSTERS}" \
  --transport-electrons "${TRANSPORT}" \
  --output-prefix "${PREFIX}" \
  --wire-pitch-mm 4 \
  --seed "${TRANSPORT_SEED}"

echo
echo "=== 3/3 Analyse wire coordinate and arrival asymmetry ==="
python3 "${STAGE_DIR}/analyze_c_pi_fast_response.py" \
  --events "${PREFIX}_event_summary.csv" \
  --output-summary "${SUMMARY}" \
  --output-u-bins "${UBINS}" \
  --u-bin-mm 0.2

echo
echo "PI fast-response workflow COMPLETE"
echo "  clusters : ${CLUSTERS}"
echo "  events   : ${PREFIX}_event_summary.csv"
echo "  wires    : ${PREFIX}_wire_summary.csv"
echo "  summary  : ${SUMMARY}"
echo "  u bins   : ${UBINS}"
