# p3_noFW only: 1000 true Geant4 J/psi muons with full TrackHeed + Garfield

## Scientific scope and exact benchmark

This is a **single geometry** native physics benchmark, NOT a LUT replay or
synthetic Poisson cluster approximation. The 1000 tracks are selected from
the previously frozen **50,000 generated J/psi parent event** Geant4 CSV
`Jpsi/chunk0003/geant_muon_crossings.csv` through
`crossing_provenance("Jpsi",3)` with SHA verification. One eligible
gas crossing is chosen per **distinct generated parent event**; exactly
1000 parent event IDs are selected via deterministic SHA ranking, independent
of Garfield outcomes. This is **one muon crossing per distinct parent event,
not both muons/all chamber crossings from each J/psi decay**. Native
Geant4 momentum/chord/angles/station/detector are retained verbatim. The
source is eligible-crossing/parent-sampled, **not an unweighted all-crossings
detector flux**, and it is not intended for final detector efficiency
measurement. Default source is chunk 3, avoiding existing old/new120 and
the separate development 2400 source in chunk 2. Source availability on
the user's WSL must be verified; if missing, choose another completed,
SHA-verified original Jpsi chunk >=3 *before preparing*.

**Physics**: original real TrackHeed::NewTrack/GetCluster/GetElectron;
one true 0.1eV thermal start `AvalancheMicroscopic::AvalancheElectron`
per conduction electron, exactly one p3_noFW geometry
(pitch 3mm, 30um anodes, 0 field wire, cathode planes +/-3mm, Ar:CO2
70:30 at 760 Torr, 293.15K, anode +1800V, B=0). Garfield
microscopic avalanche cap is the original 200000; caps and unsuccessful
outcomes are *counted and reported*, never silently treated as a correct
upper-tail sample. A 128-point unbiased spatial/time ion-birth reservoir
per avalanche is retained, but the **full ion creation census** and
sectors are kept without subsampling. **Not all ion birth coordinates**
and **not the full electron current** are saved; this is a true full
microscopic computation for measuring cost, NOT a cathode signal simulation.

This is an **isolated fork** of frozen
`phaseB_native_muon_3p3_spacetime_v1.cxx`, altered ONLY to run
`GEOMS.front() = p3_noFW` in its main `--muon-cache-physics-p3-noFW`
mode and to add native per-muon/per-avalanche timing outputs. Existing
four-geometry C++ and SHA-sealed campaigns are unchanged. A separate
new CMake target is built; old native binary is not recompiled by the
command below. A compiler smoke and the native first-20 check
are **mandatory**. The new code has NOT been built/executed on WSL
by the assistant.

## First gate: 20 real muons only, then inspect time/errors

~~~bash
cd ~/na6p/src/NA6PRoot-mwpc
git switch feature/mwpc-stage-c
git pull --ff-only
source "$HOME/na6p/env.sh"
source "$HOME/na6p/install/garfieldpp/share/Garfield/setupGarfield.sh"
set -euo pipefail

cmake -S test/mwpc-avalanche -B test/mwpc-avalanche/build \
  -DGarfield_DIR="$HOME/na6p/install/garfieldpp/lib/cmake/Garfield"
cmake --build test/mwpc-avalanche/build \
  --target mwpc_phase_b_native_p3_noFW_jpsi_bench_v1 -j2

PY=test_runs/mwpc_stage_c/.venv_pi_slides/bin/python
SCRIPT=test/mwpc-stage-c/pi_p3_noFW_jpsi1000_native_benchmark_v1.py
"$PY" -m py_compile "$SCRIPT"
"$PY" "$SCRIPT" selftest
"$PY" -u "$SCRIPT" prepare --source-chunk 3
"$PY" -u "$SCRIPT" preflight

"$PY" -u "$SCRIPT" run --max-chunks 1 --hours 2 \
  | tee "$HOME/na6p/p3_noFW_jpsi1000_first20.log"
"$PY" "$SCRIPT" status
~~~

**Do not proceed** if the C++ compile, first source selection, preflight,
first 20 native physics events or native SHA audit fails. In particular,
no `prepare` retry after success: all 1000 source IDs are frozen.
The runner creates one non-overwriting chunk folder per 20 independent
parents, all outputs + SHA `complete.json`; it refuses reruns of an
unsealed chunk. A `.writer.lock` advisory `flock` excludes two
simultaneous instances of this runner.

## 100 muons then the rest — only after first 20 PASS

~~~bash
# After 1 sealed 20-parent chunk: four more chunks -> 100 muons.
"$PY" -u "$SCRIPT" run --max-chunks 4 --hours 6 \
  | tee "$HOME/na6p/p3_noFW_jpsi1000_next80.log"
"$PY" "$SCRIPT" status

# Only if 100 muons give tolerable caps/failures and cost:
# 45 chunks remaining -> 1000 muons. Soft time cutoff BETWEEN chunks.
"$PY" -u "$SCRIPT" run --max-chunks 45 --hours 18 \
  | tee "$HOME/na6p/p3_noFW_jpsi1000_remaining900.log"
"$PY" "$SCRIPT" status
"$PY" "$SCRIPT" audit
~~~

The `status` and `audit` modes do **not** call Garfield. `status`
checks SHA of the exact source rows and each completed chunk seal,
but skips big output file hashing for speed. `audit` **rechecks all
individual saved file hashes**; run it at completion.
A live/unsealed final chunk is reported but NOT counted; inspect its
`physics.log` and running process before any resume. The runner
automatically skips previously complete, SHA-sealed chunks on resume.

## Native timings and farm decision

Each chunk contains:

- `muon_reference_muon_timing.csv`: one row per actual Geant4 muon;
  Heed time, Garfield+output time, total per-muon time.
- `muon_reference_timing.csv`: every single actual Garfield
  `AvalancheMicroscopic` call, `wall_seconds`, electron/ion counts,
  `success`, `capped`.
- `complete.json`: the actual wall-clock time **of the full C++
  subprocess**, including fresh gas initialization and native CSV I/O;
  validated total numbers of true Heed electrons and Garfield avalanches,
  SHA digest of every physics file.

`status` reports measured seconds/muon and a **provisional sequential
1000-muon estimate** with the actual 20-parent per-process startup cost
from the observed chunks. Use 100 muons (five different source-mixed
chunks) for a more stable estimate than 20. Native timings exclude Python
SHA audits/selection, CMake builds and human scheduling, so compare full
shell `time -p` if including all throughput overhead.

For a farm: multiply observed CPU seconds per *true TrackHeed electron*
by source-specific expected electron count, include measured gas-start
overhead per batch, and benchmark shared-filesystem I/O. One parent
event / chunk must not be counted as 2 independent muons. Parallel speedup
is not assumed linear; separately measure per-core efficiency, memory,
quota, software deployment, and data transfer.

Current script explicitly reports `n_capped_seeds` and
`n_unsuccessful_seeds` for diagnosis. Scientific results from
any right-censored/capped trials cannot be interpreted as an uncapped
full Garfield distribution. This is a runtime feasibility campaign.

## Immutable storage

`test_runs/mwpc_stage_c/pi_p3_noFW_jpsi1000_native_v1/`.

No modifications to old D500, A20/A100, full-cloud bank, 9000-cloud
run, old/new120 or 2400-parent native four-geometry references. No
Geant4 rerun, new muon generator, weighting fields, strip readout, or
new LUT training in this campaign.
