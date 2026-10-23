#!/usr/bin/env bash
# Six-field native Garfield/Heed electron fold over checkpointed 50k-parent
# Geant4 per-source CSV. NO Geant4 or Garfield simulation performed here.
#
# smoke: 200 real Jpsi chunk0001 MS0 muon hits; separately checkpointed.
# fold:  all first two 50k chunks x 3 sources x 6 stations; resumable.
# merge: only after ALL 36 source/station/chunk shards pass and keep exact
#        moments, produce NEW 100k parents/source PI V3 + per-anode pages.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd -P)"
REPO="$(cd "$HERE/../.." && pwd -P)"
cd "$REPO"
MODE="smoke"
if [[ $# -gt 0 ]]; then MODE="$1"; shift; fi
if [[ $# -ne 0 ]]; then
  echo "Usage: bash $0 smoke|fold|merge" >&2
  exit 2
fi
python3 -m py_compile \
 "$HERE/pi_six_geometry_native_timing.py" \
 "$HERE/pi_2m_native_six_geometry_fold.py" \
 "$HERE/test_pi_2m_native_six_geometry_fold.py"
python3 "$HERE/test_pi_2m_native_six_geometry_fold.py"
case "$MODE" in
 smoke)
   echo "REAL small test: first 200 muon gas hits in Jpsi MS0, "
   echo "from NEW 50k native ROOT chunk0001. No Garfield/G4 rerun."
   python3 "$HERE/pi_2m_native_six_geometry_fold.py" smoke \
     --channel Jpsi --station 0 --chunks 1 \
     --max-muon-hits-per-station 200
   echo "PASS: real 200-hit new-chunk native six-geometry smoke."
   ;;
 fold)
   echo "REFOLD: six frozen Garfield response maps, shared Heed births."
   echo "2 50k-parent chunks x 3 sources x 6 stations; resumable."
   echo "Full 100k/source histograms only after ALL shards finish."
   python3 "$HERE/pi_2m_native_six_geometry_fold.py" fold \
     --channel all --station -1 --chunks 2
   ;;
 merge)
   echo "STRICT 100k/channel merge from all 36 native response shards."
   echo "No per-electron pseudo-independent statistical error bars."
   python3 "$HERE/pi_2m_native_six_geometry_fold.py" merge --chunks 2
   ;;
 *)
   echo "Usage: bash $0 smoke|fold|merge" >&2
   exit 2
   ;;
esac
