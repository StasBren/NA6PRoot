#!/usr/bin/env bash
# Read existing complete full-campaign native QA outputs; make a shareable,
# conservation-checked text report. Does not invoke ROOT or Geant4.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd -P)"
REPO="$(cd "$HERE/../.." && pwd -P)"
cd "$REPO"
python3 -m py_compile "$HERE/summarize_pi_native_full_qa.py"
python3 "$HERE/summarize_pi_native_full_qa.py"
echo
echo "Send this small text file to ChatGPT for PI physics interpretation:"
echo "  $REPO/test_runs/mwpc_stage_c/pi_native_full_quality/pi_native_QA_readout.txt"

# Deliver ALL finished PI QA products in ONE small, portable file for review.
# Include all histogram bins so we can inspect Edep/angles rather than
# trusting only the pre-summarized averages. Original CSVs remain untouched.
OUT="$REPO/test_runs/mwpc_stage_c/pi_native_full_quality"
python3 - "$OUT" <<'PY'
from pathlib import Path
from zipfile import ZipFile, ZIP_DEFLATED
import sys
d=Path(sys.argv[1])
members=("pi_native_QA_readout.txt",
         "pi_native_all_QA_report.json",
         "pi_native_all_station_quality.csv",
         "pi_native_all_chamber_quality.csv",
         "pi_native_all_histograms.csv")
archive=d/"pi_native_PI_review.zip"
with ZipFile(archive,"w",ZIP_DEFLATED,compresslevel=6) as z:
    for member in members:
        source=d/member
        if not source.is_file():
            raise FileNotFoundError(source)
        z.write(source,arcname=member)
print("PASS: PI review package:",archive)
print("Please attach THIS ZIP file in ChatGPT")
PY
