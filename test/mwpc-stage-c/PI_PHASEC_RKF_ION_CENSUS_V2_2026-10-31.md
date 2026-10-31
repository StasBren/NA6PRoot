# Phase-C V2 — RKF estimated avalanche ion population on 100 exact D500 starts

**Status 2026-10-31:** V2 standalone C++ and Python runner committed; checked user-installed Garfield header signature supports `void DriftLineRKF::GetAvalancheSize(double& ne, double& ni) const`. Offline Python synthetic statistical self-test passed. **Not yet compiled/physically run against user's installed Garfield++; do not claim physical accuracy.**

## Why this test exists

Frozen Phase-C V1 reported `RKF::GetGain()` from integrals of Townsend coefficient, which is **not identical** to the original `AvalancheMicroscopic` output `ion_births`. V2 records both `rkf_electrons_estimate` and `rkf_ions_estimate` through the **confirmed local header** `GetAvalancheSize(double&, double&)`, without changing V1 source, gas, 100 exact D500 starts, V1 output, 84k original A20 trials, 16k A100 trials or 120 real muon full microscopic archive. This is a different, **macroscopic** model's estimated particle population, not one random microscopic avalanche nor ion-cloud spatial coordinates.

### New files

- `test/mwpc-avalanche/phaseC_rkf_ion_census_pilot_v2.cxx`: V1 geometry and transport unchanged; adds two macroscopic RKF output fields per start.
- `test/mwpc-avalanche/CMakeLists.txt`: registers **separate** target `mwpc_phase_c_rkf_ion_census_pilot_v2`; original V1 target remains untouched.
- `test/mwpc-stage-c/pi_run_rkf_ion_census_v2.py`: SHA-guarded run, comparison, selftest and status. Preserves V1 original code SHA recorded in the immutable plan and verifies gasfile/output/manifests. V2 refuses overwrite and writes to a separate directory.
- This document and updated handoff.

## Exact WSL commands

```bash
(
  set -euo pipefail
  cd ~/na6p/src/NA6PRoot-mwpc
  git pull --ff-only
  python3 test/mwpc-stage-c/pi_run_rkf_ion_census_v2.py selftest

  # Reconfigure old CMake cache to register the new target.
  cmake -S test/mwpc-avalanche -B test/mwpc-avalanche/build
  cmake --build test/mwpc-avalanche/build \
    --target mwpc_phase_c_rkf_ion_census_pilot_v2 -j2

  python3 -u test/mwpc-stage-c/pi_run_rkf_ion_census_v2.py buildcheck
)
```

If and only if this passes, execute:

```bash
(
  set -euo pipefail
  cd ~/na6p/src/NA6PRoot-mwpc
  python3 -u test/mwpc-stage-c/pi_run_rkf_ion_census_v2.py run
  python3 -u test/mwpc-stage-c/pi_run_rkf_ion_census_v2.py compare
  python3 test/mwpc-stage-c/pi_run_rkf_ion_census_v2.py status
)
```

**Do not rerun** the V1 commands `prepare`, `gas` or `run`. V2 reuses the already successful, SHA-pinned V1 100 exact starts, its B=0 Ar70/CO2 30 Magboltz gas table, and existing original microscopic A100 trials; it does not make more avalanches or new gas tables.

If compilation fails, copy full native compiler error; only V2 files may be adjusted, never frozen V1 C++. If V2 fails to reproduce the V1 geometry/drift `GetGain()` and `GetLoss()` for the same coordinates, the code stops before presenting a comparison. Invalid RKF estimates (negative, nonfinite) are rejected; if all estimated ion populations are 0 while gain is nonzero, **investigate** how `DriftLineRKF::GetAvalancheSize` is populated with `EnableSignalCalculation(false)` on the installed Garfield version, not assume zero physical ions.

## Output schema and interpretation

Separate root: `test_runs/mwpc_stage_c/pi_RKF_ion_census_B0_pilot_v2`.

- `rkf_estimated_ions_100_exact_starts.csv`: both `rkf_electrons_estimate` and `rkf_ions_estimate`, original V1 `rkf_mean_gain`, `rkf_loss`, geometry, start coords, transport endpoint, time, physically in-metal indicator and old illustrative proxy.
- `matched_100_micro_vs_RKF_ions.csv`: exact `(geometry,box_id,sample)` join with source `trials_extra_[geometry].csv`, including `ion_births`, `garfield_ni`, `collected_e`, estimated RKF electrons/ions, per-event residual.
- `summary.json`: per-geometry and aggregate means, zero fractions, sample-level Pearson, MAE and run time; comparisons are exploratory and **not** a calibration or speed-up measurement.
- `run_receipt.json` and `comparison_manifest.json`: SHA of sources, original gas, original V1 runs, new C++ binary and outputs.

**Physics caution:** User's installed header confirms a count API, not that the macroscopic RKF's count has the exact same collision-level definition/variance as `AvalancheMicroscopic`. Only 25 microscopic realizations per geometry exist in this test; a difference is not necessarily systematic error at a fixed position. V1 results suggest at least 100 successful RKF transport traces, but no starts inside metal were encountered, so metal handling remains untested. No cathode strip signals, ion positions, magnetic field, full 120-muon reconstruction or unseen Geant4 source is represented here.

Before interpreting zero estimated RKF `ni`, check local `~/na6p/src/garfieldpp/Source/DriftLineRKF.cc` for whether count accumulation requires signal calculation. Source and installed headers determine the exact semantic; do NOT silently enable signal/ion-tail and still claim the same V1 numerical method.

## External context
- Garfield++ User Guide: https://garfieldpp.web.cern.ch/doxygen/UserGuide.pdf
- Garfield++ discussion of deterministic macroscopic RKF vs microscopic avalanche fluctuations: https://root-forum.cern.ch/t/help-on-garfield-avalanchemicroscopic-and-avalanchemc/53678
