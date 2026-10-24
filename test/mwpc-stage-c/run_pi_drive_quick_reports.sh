#!/usr/bin/env bash
# Rapid curated Stage-C Drive publication: a handful of small status,
# physics diagnostics and numeric summaries only. No full ROOT, 120 CSVs,
# gzip shards, multi-MB merged station JSON, huge HTML or plots.
#
# Stops on persistent Drive upload errors; does NOT delete files.
# Can be used while local Heed/Garfield six-geometry folds are still running.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd -P)"
REPO="$(cd "$HERE/../.." && pwd -P)"
BASE="$REPO/test_runs/mwpc_stage_c"
REMOTE="${NA60_DRIVE_REMOTE:-na60drive:}"
[[ "$REMOTE" == *: ]] || {
  echo "ERROR: remote must be a root such as na60drive:" >&2;exit 2;
}
command -v rclone >/dev/null || {
  echo "ERROR: rclone not installed" >&2;exit 3;
}
MODE="${1:-copy}"
[[ "$#" -eq 1 && "$MODE" =~ ^(check|copy)$ ]] || {
  echo "Usage: bash $0 check|copy" >&2;exit 2;
}
echo "QUICK research publication only: $REMOTE"
echo "Authoritative physics and original checkpoint files REMAIN on WSL."
python3 "$HERE/pi_2m_lightweight_progress_snapshot.py"
# Explicit allowlist intentionally avoids large HTML and sparse grouped
# gzip physics checkpoints. Full high-volume results stay local until
# the separate cloud-backup campaign succeeds.
FILES=(
  "pi_weekend_2m_operations_v1/lightweight_progress_snapshot.json"
  "pi_weekend_2m_operations_v1/native_gas_path_audit.json"
  "pi_weekend_2m_operations_v1/runtime_comparison.json"
  "pi_weekend_2m_operations_v1/runtime_comparison.txt"
  "pi_six_geometry_2m_timing_v1/full/anode_timing_analysis_en/summary.txt"
  "pi_six_geometry_2m_timing_v1/full/anode_timing_analysis_en/station_summary.csv"
  "pi_six_geometry_2m_timing_v1/full/anode_timing_analysis_en/station_anode_timing.json"
  "pi_six_geometry_2m_timing_v1/full/pages/english_dashboard_manifest.json"
  "pi_pi_meeting_20261025_v1/selected_data.json"
  "pi_pi_meeting_20261025_v1/README.txt"
  "pi_pi_meeting_20261025_v1/01_time_six_station.png"
  "pi_pi_meeting_20261025_v1/01_time_six_chamber.png"
  "pi_pi_meeting_20261025_v1/02_time_overlay_station.png"
  "pi_pi_meeting_20261025_v1/03_numbered_anodes_mean_rms_chamber.png"
  "pi_pi_meeting_20261025_v1/04_electron_approach_vs_muon_theta_u_heatmaps.png"
  "pi_pi_meeting_20261025_v1/04_electron_approach_vs_muon_theta_w_heatmaps.png"
  "pi_pi_meeting_20261025_v1/05_electron_approach_angle_slices_theta_u.png"
  "pi_pi_meeting_20261025_v1/05_electron_approach_angle_slices_theta_w.png"
  "pi_pi_meeting_20261025_v1/06_chamber_vs_station.png"
  "pi_pi_meeting_20261025_v1/07_time_per_physical_anode_heatmaps.png"
)
rclone lsf "$REMOTE" --max-depth 1 >/dev/null || {
  echo "ERROR: Drive remote not reachable" >&2;exit 4;
}
count=0
for rel in "${FILES[@]}"; do
  file="$BASE/$rel"
  if [[ ! -f "$file" ]]; then
    echo "SKIP not yet created: $rel"
    continue
  fi
  bytes="$(stat -c %s "$file")"
  if (( bytes > 8388608 )); then
    echo "SKIP unexpectedly large report (>8 MiB): $rel ($bytes B)"
    continue
  fi
  echo "REPORT: $rel ($bytes B)"
  if [[ "$MODE" == "copy" ]]; then
    rclone copyto "$file" "${REMOTE}test_runs/mwpc_stage_c/$rel" \
      --checksum --transfers 1 --drive-chunk-size 8M --tpslimit 2 \
      --timeout 60s --contimeout 20s --retries 1 \
      --low-level-retries 1 -v
    count=$((count+1))
  fi
done
echo "PASS: quick $MODE processed $count reports, no local or Drive deletion."
