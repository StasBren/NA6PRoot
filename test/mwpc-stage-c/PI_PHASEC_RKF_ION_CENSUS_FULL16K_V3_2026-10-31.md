# NA60+ MWPC Phase-C V3 — complete 16,000-start matched RKF ion census

**Status:** V3 standalone C++ source and independent Python analysis code committed and Python offline selftest passed. **No V3 actual Garfield++ compilation or execution has happened on the user's WSL yet.** V1, V2, 84k original Garfield A20, 16k original Garfield A100, 120 original full muon references and prior gasfile are immutable.

## Scientific reason

The user confirmed the 100-start Phase-C V2 RKF `GetAvalancheSize(double& ne,double& ni)` benchmark completed. For 25 starts per geometry, `mean_RKF_ions / mean_microscopic_Garfield_ion_births` was 1.464 / 1.761 / 0.946 / 1.359 for p3_noFW / p3_FW75 / p4_noFW / p4_FW75, respectively; 17/100 microscopic events had zero ion births but 0/100 macroscopic RKF expected ion populations were zero. Intra-geometry per-start correlations were weak, which is not proof of physics failure because one stochastic microscopic realization does **not** equal an expected response at that coordinate. Small n=25 and high-variance selected cells make calibration premature.

**V3** instead runs only fast RKF on **all 16,000 saved A100 original D500 start positions**, each with one *already existing* microscopic Garfield avalanche, grouped into the **same 200 selected spatial cells with 80 new starts per cell**. This makes cell-mean calibration and cell-wise variance much more meaningful. Still **source-/cell-selection conditioned** on original A20 variance and 120 initial Heed positions; NOT external independent-muon validation.

## New repository assets

- `test/mwpc-avalanche/phaseC_rkf_ion_census_full16k_v3.cxx` — exactly the successful V2 RKF geometry/model, with only input-size cap raised to 16,000. Records direct RKF estimated electron/ion counts, gain and drift outcome; no new microscopic Garfield++ physics.
- CMake target `mwpc_phase_c_rkf_ion_census_full16k_v3` — separate from V1 and V2; do not modify previous sources.
- `test/mwpc-stage-c/pi_run_rkf_ion_census_full16k_v3.py` — immutable archive verification, native parser check, run, read-only paired A100 comparison and SHA-sealed outputs. Performs explicit coordinate matching and *reproduces all 100 previous V2 ion estimates* from the same D500 start positions.
- The original V2 C++ contains an **illustrative** randomly sampled Pólya proxy column; V3 inherits this extra output but **does not use it in any quantitative analysis**. No claim that it is ion births or a calibrated stochastic model.

### Exact WSL steps

First build and validate. The already-prepared V1/V2 input and gas are reused. Do NOT rerun `prepare`, `gas`, either V1/V2 `run`, or any microscopic physics:

```bash
(
set -euo pipefail
cd ~/na6p/src/NA6PRoot-mwpc
git pull --ff-only
python3 -u test/mwpc-stage-c/pi_run_rkf_ion_census_full16k_v3.py selftest
cmake -S test/mwpc-avalanche -B test/mwpc-avalanche/build
cmake --build test/mwpc-avalanche/build --target mwpc_phase_c_rkf_ion_census_full16k_v3 -j2
python3 -u test/mwpc-stage-c/pi_run_rkf_ion_census_full16k_v3.py buildcheck
)
```

Expected: `PASS C++ ORIGINAL D500 exact-start parser (16000 rows)`. Then:

```bash
(
set -euo pipefail
cd ~/na6p/src/NA6PRoot-mwpc
python3 -u test/mwpc-stage-c/pi_run_rkf_ion_census_full16k_v3.py run
python3 -u test/mwpc-stage-c/pi_run_rkf_ion_census_full16k_v3.py compare
python3 -u test/mwpc-stage-c/pi_run_rkf_ion_census_full16k_v3.py status
)
```

**Original input:** `test_runs/mwpc_stage_c/pi_A100_targeted_50cells_per_geometry_v1/selected_exact_old_D500_starts.csv` (16k exact coordinates). **Original source:** four `trials_extra_<geom>.csv` (each 4k true microscopic trials). These are SHA verified against the published A100 physics manifest `2a8520c61dc38dc2b32aa31319c3665e36d9ed86f9119e8f3737f32849d7da82`. V1 existing Magboltz gas SHA pinned `f6f6a8a561d52c788c7fe8168256d3add78bdca131172fcf9526159707e1e693`. V1 result and V2 result are separately pinned by known SHA. No new gas table is generated.

**Exclusive new output root:** `test_runs/mwpc_stage_c/pi_RKF_ion_census_B0_A100all16k_v3/`. Contains:
- `rkf_estimated_ions_16k_exact_starts.csv`: RKF direct macro ne/ni for all exact positions; inherited *unused* illustrative Pólya proxy column.
- `run_receipt.json`: identity/SHA of original gas, physics, C++ source, compiled binary and new RKF output.
- `matched_16k_micro_vs_rkf_ions.csv`: exact start-wise join on `(geometry,box_id,sample)` and verified coordinates.
- `cell_80_micro_vs_80_rkf.csv`: 200 independent selected-cell summaries, each with 80 original micro trials.
- `summary.json`: per geometry 4000-start averages, 50-cell Pearson correlation, **cluster bootstrap across cells** for RKF/micro mean ratio, and two-way held-out cell calibration diagnostic (25 train / 25 evaluation cells per fold).
- `comparison_manifest.json`: SHA-sealed comparison outputs. Each stage refuses overwrite.

**Important numerical guard:** comparison **stops** rather than fabricates a physical zero if RKF returns a nonmetal transport failure or a physical-metal start that prevents like-for-like cell-mean comparison. Classify such starts and fix the statistical estimand separately if they occur.

## Physical interpretation / limitations

- `GetAvalancheSize(ne, ni)` from macroscopic RKF is an **expected population estimate** from Townsend/attachment gas table, not an event-level random microscopic `AvalancheMicroscopic` collision history.
- 80 *different original D500 starts* per cell, not 80 identical exact-coordinate microscopic replicas: the cell mean integrates both position distribution inside that bin and collision randomness. Correlation of the two cell means tests broad variations across cells, not exact pointwise conditional expectation. This is enough to assess suitability as a fast average-gain approximation but not a fully calibrated stochastic replacement.
- Spatial cell fold calibration uses the same 200 old-A20-selected cells and 16k source simulation; it is an internal test, not wholly unseen detector geometry/Geant4 muons. Avoid deployment conclusions based on this alone.
- Original D500 drift-only success/collected flags are not overwritten or equated to microscopic avalanche zero; 17% zeros on original 100 are real microscopic outcomes.
- Magboltz `ncoll=5` 40-point field grid may require stability checks before interpreting persistent physical bias.
- B=0, no ion drift, no `A_{ion}` sectors and no induced U/V cathode signal in RKF V3. Even if cell-mean ion counts agree, this does not constitute a strip-signal simulator.

**Native V3 results are pending.** Do not assert improvement, 16k correlation, crossvalidated gain correction or speedups until WSL output is seen.

## Native 16k completion and sealed reporting bug (2026-10-31)

User WSL completed compilation, `buildcheck`, `run`, `compare`, `status`: 16,000 RKF physics estimates generated and SHA-sealed, 16,000 exact A100 microscopic starts matched, 100 shared V2 values verified, 200 cells × 80 starts verified. User's reported `status` confirms all output artifacts exist. Original `summary.json` is **internally inconsistent**: it says `n_transport_ok=0`, `n_nonmetal_fail=16000`, but all 200 cell-level diagnostics require 80 successful RKF outcomes per cell and return nonempty means, correlations and crossfit values.

Root cause confirmed in the original `summary_for`: it compared `rkf_transport_ok` and `rkf_metal` (integers after constructing `trial`) with **string literals** `'1'`, so it misreported every physical trial as failed. `run` and `compare` actual physical pair construction, cell means, ratios, cell bootstrap and crossfit used integer conversion and completed successfully. **Do not delete or modify the original V3 `summary.json` or its SHA-sealed source; do not rerun `run` or `compare`.** The main Python script's `summary_for` was fixed to normalize both integer/string flags, with regression `selftest` (commit `ee7b6bf`). Original summary remains a frozen incorrect historical report.

**Separate full 16k read-only audit:** `test/mwpc-stage-c/pi_audit_rkf_v3_16k_status.py`. It verifies original V3 output receipt and comparison SHA, all 16k raw RKF vs paired CSV identities and actual independent microscopic Garfield original trials, all 200 stored cell means, original cell correlations and fold calibration, and recomputes corrected transport and ion summary. Writes exclusively to a **new** sibling directory `test_runs/mwpc_stage_c/pi_RKF_ion_census_B0_A100all16k_v3_status_audit_v1` and refuses overwrite. Local offline synthetic `selftest` PASS, plus a full-size 16k synthetic integration with the real 16k original micro trial files and fully fabricated RKF responses PASS. **User's native real-data audit has NOT run yet.**

To run **only offline audit** in WSL:

```bash
(
set -euo pipefail
cd ~/na6p/src/NA6PRoot-mwpc
git pull --ff-only
python3 test/mwpc-stage-c/pi_run_rkf_ion_census_full16k_v3.py selftest
python3 test/mwpc-stage-c/pi_audit_rkf_v3_16k_status.py selftest
python3 -u test/mwpc-stage-c/pi_audit_rkf_v3_16k_status.py audit
)
```

If audit prints PASS, use `test_runs/mwpc_stage_c/pi_RKF_ion_census_B0_A100all16k_v3_status_audit_v1/corrected_summary.json`, **not** original v3 `summary.json`, when reporting overall transport success and RKF mean ion count. It also creates `audit_note.md` and `audit_manifest.json` with SHA of new reports.

**Scientific preliminary numbers from user-reported original cell means (independent of the string/int bug):** estimated RKF/microscopic mean ion ratios: p3_noFW ~1.273, p3_FW75 ~1.333, p4_noFW ~1.334, p4_FW75 ~1.321. Cell-mean Pearson correlations only 0.0498/0.0830/0.1133/0.1801. A simple per-geometry training-cell ratio rescale decreased average held-out cell MAE from 157.7→59.1, 1906.0→427.5, 490.1→160.4, 1738.3→420.6, respectively. These are **internal selected-cell measurements, not independent Geant4 muon validation**; still not a calibrated stochastic avalanche response. Zero microscopic ion births: 2708/16000. Do not interpret any numerical `n_transport_ok=0` or nonexistent `mean_rkf_ions` from the old report. Reconfirm with real native audit before publication.
