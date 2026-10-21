#!/usr/bin/env bash
# PI native Geant4 gas-hit campaign for the CURRENT installed 298-chamber
# material geometry, plus full original MCKine and Edep-step data.
#
# A smoke run is required before going long (also precompiles ROOT macros).
#
#   bash test/mwpc-stage-c/run_pi_geant4_campaign.sh smoke
#   bash test/mwpc-stage-c/run_pi_geant4_campaign.sh parallel-smoke
#   bash test/mwpc-stage-c/run_pi_geant4_campaign.sh full
#
# Optional full arguments: events_per_source chunk_size workers hours
#   bash test/mwpc-stage-c/run_pi_geant4_campaign.sh full 30000 3000 4 2.5
#
# The full command is resumable: completed chunks are never overwritten.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
cd "$REPO"
# Always use the installed native executable AND libsimLib built from THIS
# checkout. The user's default ~/na6p/install/NA6PRoot may be legacy.
# The activation affects this subshell only, not the user's login profile.
NA6P_STAGEC_REQUIRE_ROOT_HADD=1
source "$HERE/activate_pi_stagec_env.sh"
python3 -m py_compile "$HERE/run_pi_geant4_campaign.py" || {
 echo "Geant4 campaign Python syntax check failed" >&2
 exit 2
}
MODE="smoke"
if [[ $# -gt 0 ]]; then MODE="$1"; shift; fi
case "$MODE" in
 smoke)
   if [[ $# -gt 0 ]]; then
     echo "Smoke takes no extra arguments" >&2
     exit 2
   fi
   python3 "$HERE/run_pi_geant4_campaign.py" --mode smoke
   ;;
 parallel-smoke)
   if [[ $# -gt 0 ]]; then
     echo "Parallel smoke takes no extra arguments" >&2
     exit 2
   fi
   python3 "$HERE/run_pi_geant4_campaign.py" --mode parallel-smoke
   ;;
 full)
   EVENTS="30000"
   CHUNK="3000"
   WORKERS="4"
   HOURS="2.5"
   if [[ $# -gt 0 ]]; then EVENTS="$1"; shift; fi
   if [[ $# -gt 0 ]]; then CHUNK="$1"; shift; fi
   if [[ $# -gt 0 ]]; then WORKERS="$1"; shift; fi
   if [[ $# -gt 0 ]]; then HOURS="$1"; shift; fi
   if [[ $# -gt 0 ]]; then
     echo "Too many full arguments: [events_per_source] [chunk_size] [workers] [hours]" >&2
     exit 2
   fi
   python3 "$HERE/run_pi_geant4_campaign.py" \
     --mode full --events-per-source "$EVENTS" \
     --chunk-size "$CHUNK" --workers "$WORKERS" --hours "$HOURS"
   ;;
 *)
   echo "Usage: bash $0 smoke | parallel-smoke | full [events_per_source] [chunk_size] [workers] [hours]" >&2
   exit 2
   ;;
esac
