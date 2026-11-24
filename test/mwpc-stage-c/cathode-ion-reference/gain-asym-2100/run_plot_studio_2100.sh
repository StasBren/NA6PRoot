#!/usr/bin/env bash
# NA60+ 2100 V Plot Studio: existing V0 CSV outputs only; no physics simulation.
set -Eeo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../../../.." && pwd)"
VENV="${NA60_PLOTSTUDIO_VENV:-$HOME/na6p/venvs/na60-plotstudio-2100}"
DATA="${NA60_PLOTSTUDIO_DATA_ROOT:-$REPO/test_runs/mwpc_stage_c/2100_gain_asym}"
MODE="${1:-check}"
fail(){ echo "NA60+ Plot Studio: $*" >&2; exit 2; }
case "$MODE" in
 install)
   command -v python3 >/dev/null || fail "python3 missing"
   mkdir -p "$(dirname "$VENV")"
   if [[ ! -f "$VENV/bin/activate" ]];then
     python3 -m venv "$VENV" || fail "venv unavailable. Install your distro's python3-venv package."
   fi
   "$VENV/bin/python" -m pip install --upgrade pip
   "$VENV/bin/python" -m pip install -r "$HERE/requirements_plot_studio_2100.txt"
   echo "PASS Plot Studio Python dependencies installed in $VENV";;
 check)
   [[ -x "$VENV/bin/python" ]] || fail "run '$0 install' first"
   "$VENV/bin/python" -m py_compile "$HERE/plot_studio_2100.py" "$HERE/per_anode_census_2100.py" "$HERE/test_per_anode_2100.py"
   "$VENV/bin/python" -m unittest discover -s "$HERE" -p "test_per_anode_2100.py" -q
   "$VENV/bin/python" -c 'import streamlit,plotly,pandas,numpy; print("PASS Streamlit",streamlit.__version__,"Plotly",plotly.__version__)'
   [[ -d "$DATA" ]] || fail "no analysis results found in $DATA; run gain-asym pilot first"
   count=0
   shopt -s nullglob
   for d in "$DATA"/pilot_*/census "$DATA"/population_*/census;do
     [[ -s "$d/v0_muons.csv" && -s "$d/v0_primaries.csv" && -s "$d/v0_fragments.csv" && -s "$d/v0_wires.csv" ]] || continue
     echo "FOUND: $d"; ((count+=1))
   done
   ((count>0)) || fail "no complete census outputs (pilot or population) found"
   echo "PASS $count saved analysis folders";;
 run)
   [[ -x "$VENV/bin/streamlit" ]] || fail "run install first"
   [[ -d "$DATA" ]] || fail "no data directory $DATA"
   export NA60_PLOTSTUDIO_DATA_ROOT="$DATA"
   PORT="${2:-8501}"
   [[ "$PORT" =~ ^[0-9]+$ ]] || fail "port must be numeric"
   echo "Starting read-only NA60+ Plot Studio at http://localhost:$PORT"
   echo "Input data: $NA60_PLOTSTUDIO_DATA_ROOT"
   exec "$VENV/bin/streamlit" run "$HERE/plot_studio_2100.py" \
     --server.address 127.0.0.1 --server.port "$PORT" \
     --browser.gatherUsageStats false;;
 *)
   fail "Usage: bash $0 install|check|run [port]";;
esac
