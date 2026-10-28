# NA60+ / MWPC — first TRUE-MUON avalanche B1 (±3 mm cathodes)

**28 October 2026 · code in feature/mwpc-stage-c · C++ compile and Garfield results NOT YET verified on WSL**

## Objective

One original validated **Geant4 gas crossing from the J/ψ muon sample** is the physical observation.
Read its real local entry/exit, signed arrival angles, **already exported** local
muon momentum, full parent/hit/det/station IDs. Use **TrackHeed** once on that
real muon chord to create an actual collection of primary electron clusters.
Replay that exact collection in four *separate microscopic Garfield avalanches*:

| key | pitch | field wire | gaps | role |
| --- | --- | --- | --- | --- |
| **p3_noFW** | **3 mm** | none | **3+3 mm** | **baseline** |
| p3_FW75 | 3 mm | grounded 75 µm | 3+3 mm | control |
| p4_noFW | 4 mm | none | 3+3 mm | control |
| p4_FW75 | 4 mm | grounded 75 µm | 3+3 mm | control |

Anode diameter 30 µm, +1800 V, Ar/CO₂ 70:30, T=293.15K,
P=760 Torr, B=0T. **Not gain-matched; this is fixed-HV physics.**

For every Heed electron seed, the microscopic Garfield ionization callback
captures births from **all generations** in that seed's avalanche. We sum
raw ion birth counts across **every** seed of one **original muon crossing**
before forming A. This is not average per-primary-electron A, and not an
induced cathode signal:

- `A_ion_hemi=(N_up-N_down)/(N_up+N_down)` over the entire 360° anode surface.
- `A_ion_TB=(N_top45-N_bottom45)/(N_top45+N_bottom45)` over
  top [45,135) and bottom [-135,-45) degree sectors. Others are excluded.
- Both use the **nearest physical anode** (a midpoint field wire is not
  an anode) in a bounded, translated local 9-wire analytic-field model.
- Denominator zero = blank / undefined; failed/zero avalanches are retained.
- Also save amplified electron endpoint A and effective full-hit gain
  `N_collected_amplified_e / N_true_Heed_primary_e`. This is not A_Q(T);
  actual cathode induced signal is a subsequent B3 task.

## Source integrity

Original Geant4 files are unmodified:
`test_runs/mwpc_stage_c/pi_six_geometry_2m_native_exports_v1/jpsi/chunk0001/geant_muon_crossings.csv`.

`pi_2m_native_six_geometry_fold.py::crossing_provenance` SHA-checks the
original CSV and original ROOT metadata. The macro
`export_c_geant_muon_crossings.C` already wrote:
`pu_in_GeV,pv_in_GeV,pw_in_GeV` plus both signed muon angles.
**No new ROOT export or fixed 10-GeV surrogate momentum is needed.**
The input selector uses a tiny deterministic per-station/per-sign sample
within 16 mm u-edge protection, and rejects momentum–gas-chord direction
misalignment. This deliberately **is not a sample reproducing the physical
muon incidence-angle PDF**. The old `full34` Geant4+empirical Heed
marked-Poisson study remains a different frozen data product.

## Step 1. Build just the NEW target on WSL

```bash
cd ~/na6p/src/NA6PRoot-mwpc
git switch feature/mwpc-stage-c
git pull --ff-only
bash test/mwpc-stage-c/run_pi_native_muon_b1.sh build
```

This runs offline regression tests and compiles ONLY
`mwpc_phase_b_native_muon_3p3` with the existing
`~/na6p/install/garfieldpp`. There is no new Geant4 run.

## Step 2. Select exactly TWO verified genuine muon gas crossings

```bash
bash test/mwpc-stage-c/run_pi_native_muon_b1.sh prepare2
bash test/mwpc-stage-c/run_pi_native_muon_b1.sh check2
```

Selection and SHA provenance:
`test_runs/mwpc_stage_c/pi_muon_avalanche_b1_v1/smoke2_chunk0001/`.
It chooses two physical crossing rows using a reproducible hash rank in
signed angle/station strata; run output is not co-mingled with gain pilot100.

## Step 3. Real full chain smoke: 2 muons × 4 geometries

```bash
bash test/mwpc-stage-c/run_pi_native_muon_b1.sh smoke2
bash test/mwpc-stage-c/run_pi_native_muon_b1.sh status
```

The non-overwriting output directory is:
`test_runs/mwpc_stage_c/pi_muon_avalanche_b1_v1/run_smoke2_chunk0001/`.

Important files:
- `native_muon_b1_crossing.csv`: **one original physical muon crossing × one geometry**; four rows per muon, signed angles, Heed multiplicity, full-hit A_hemi/A_TB, gain, zeros/failures/caps.
- `native_muon_b1_seed.csv`: each Heed electron and its own avalanche bookkeeping, angle counts and gain.
- `native_muon_b1_ion_birth_sample.csv`: a **bounded, unbiased reservoir of up to 100 real ion-birth positions per seed**, using an **independent RNG** that does not perturb Garfield. Each sampled birth carries `seed_ion_births_total` and `reservoir_weight` to allow weighted spatial plots across unequal avalanche sizes. Physics asymmetries use **every ionization birth**, not the sample.
- `native_muon_crossings_input.csv`, `physics.log`, `run_config.json`, `summary.json`, `manifest.json`: inputs, human-readable physics trace, true source hashes and immutable SHA-pinned outputs.

An interrupted/failed run remains on disk for diagnosis and is NEVER
automatically retried into the same output folder. Report compilation errors,
failed seed IDs, limits or identity mismatches before broadening sample.
If no births are recorded from all four geometries the audit fails closed.

## Step 4. Upload the two-crossing result for direct ChatGPT reading

Existing WSL remote is **`na60drive:`**, already authorized and tested.
No reconfiguration or 252 GiB re-upload needed.

```bash
# Default is DRY RUN, no cloud write.
python3 test/mwpc-stage-c/pi_publish_run_to_drive.py \
  --source test_runs/mwpc_stage_c/pi_muon_avalanche_b1_v1/run_smoke2_chunk0001 \
  --run-id muon_B1_true_Heed_smoke2_20261028

# Execute only after inspecting the dry run:
python3 test/mwpc-stage-c/pi_publish_run_to_drive.py \
  --source test_runs/mwpc_stage_c/pi_muon_avalanche_b1_v1/run_smoke2_chunk0001 \
  --run-id muon_B1_true_Heed_smoke2_20261028 \
  --upload
```

The file-limited publisher reuses
`na60drive:avalanche_3p3_runs`, copies only selected small CSV/log/manifests
and refuses overwriting a different remote file.

## Step 5. Eight physical crossings after TWO succeed

```bash
bash test/mwpc-stage-c/run_pi_native_muon_b1.sh prepare8
bash test/mwpc-stage-c/run_pi_native_muon_b1.sh check8
bash test/mwpc-stage-c/run_pi_native_muon_b1.sh smoke8
```

Output `run_smoke8_chunk0001/` is separate and retains 2-crossing samples;
**the first two of the eight repeat the smoke2 source crossing IDs**, so do not
treat them as an independent replicate. Eight crossings are not enough for a
physics RMS(A) versus entrance-angle result; they test the full processing
chain. After both modes are verified, add representative native-muon samples
with correct angle-bin exposure/sampling weights, parent-correlated errors,
finite detector U/W boundaries and true-cathode B3 signal.

**Status:** the code and scripts have been committed to GitHub but this
document does not claim successful Garfield build, physics output, or a
new B1 scientific measurement.


## 28 Oct 2026: initial smoke2 failure — CSV newline format (FIXED IN SOURCE)

The first native B1 build and source selection PASSED. The original
`run_smoke2_chunk0001/physics.log` recorded:

```
B1 FAILURE: Wrong real Geant4 muon input column schema
```

**Root cause:** Python `csv.DictWriter` produces `CRLF` line endings
(`\\r\\n`) even on WSL/Linux. The original C++ `std::getline` strips
the `\\n` but retains `\\r`, so an exact string header comparison fails.
The underlying 2 selected real Geant4 crossings and immutable selection
CSV are valid, and no Garfield avalanche physics had begun.

**Fix:** C++ `ReadCrossings` now explicitly strips terminal CR from
both header and each row; the wrong-header exception now prints actual
and expected headings. C++ `--validate-input` checks the original file
**without initializing Garfield gas or performing any physical simulation**.
The Python `check2` mode invokes this C++ test before every B1 smoke.

**Preserve the original failed run unchanged**. The corrected physics
attempt writes to `run_smoke2_crlf_fix1_chunk0001/`:

```bash
cd ~/na6p/src/NA6PRoot-mwpc
git pull --ff-only
bash test/mwpc-stage-c/run_pi_native_muon_b1.sh build
bash test/mwpc-stage-c/run_pi_native_muon_b1.sh check2
bash test/mwpc-stage-c/run_pi_native_muon_b1.sh smoke2retry1
```

If successful, upload this NEW directory via the existing `na60drive:`
project-root remote:

```bash
python3 test/mwpc-stage-c/pi_publish_run_to_drive.py \
  --source test_runs/mwpc_stage_c/pi_muon_avalanche_b1_v1/run_smoke2_crlf_fix1_chunk0001 \
  --run-id muon_B1_true_Heed_smoke2_crlf_fix1_20261028 \
  --upload
```

Do not execute the old `smoke2` mode against its existing directory,
and do not delete, mutate or re-upload the incomplete old run.


## 28 Oct 2026: second smoke2 failure — non-ionisable Heed source (FIXED IN SOURCE)

The second attempt `run_smoke2_crlf_fix1_chunk0001/physics.log` recorded:

```text
ComponentAnalyticField::CellCheck:
    All potentials in the cell are the same.
    There is no point in going on.
B1 FAILURE: No ionisable medium at initial position in
Garfield::TrackHeed::NewTrack(...)
```

**Root cause**: The deliberately geometry-independent Heed primary
ionisation source was represented by a `ComponentAnalyticField` with
two identical-potential planes (0 V / 0 V), no anode wires. Garfield's
analytic cell check correctly rejected this degenerate electrostatic
model, leaving `TrackHeed` without an ionisable gas medium. This
happened **before any ionisation/avalanche physics**, and is NOT a
negative result about the actual ±3 mm MWPC geometry or field wires.

**Fix**: the Heed source now uses `ComponentConstant` with an explicitly
bounded gas volume (`SetArea` + `SetMedium`), a weak 1 V/cm
bookkeeping field, and sensor bounds. This *source* field does not
replace or alter any of the FOUR proper 1800 V wire fields used by
`AvalancheMicroscopic`. The initial muon position/direction still
comes from its real frozen Geant4 chord/momentum, and the SAME Heed
electrons are reused across all four geometries.

`check2` now runs a real C++ `--validate-input` mode that not only
parses the frozen CSV but verifies `Sensor::GetMedium` and
`IsIonisable` at **each actual muon start**, without running
TrackHeed ionisation or any avalanches. Required new output:

```text
PASS: parsed 2 original Geant4 muon crossings (CSV CRLF/LF supported)
PASS: Heed gas source medium is ionisable at all 2 muon starts (±3 mm)
```

**Third attempt: NEW immutable output**, preserving BOTH prior failed
runs verbatim:

```bash
cd ~/na6p/src/NA6PRoot-mwpc
git pull --ff-only
bash test/mwpc-stage-c/run_pi_native_muon_b1.sh build
bash test/mwpc-stage-c/run_pi_native_muon_b1.sh check2
bash test/mwpc-stage-c/run_pi_native_muon_b1.sh smoke2retry2
```

Only after `smoke2retry2` ends with a fully validated manifest and
a physical output should it be published:

```bash
python3 test/mwpc-stage-c/pi_publish_run_to_drive.py \
  --source test_runs/mwpc_stage_c/pi_muon_avalanche_b1_v1/run_smoke2_gas_medium_fix2_chunk0001 \
  --run-id muon_B1_true_Heed_smoke2_gas_medium_fix2_20261028 \
  --upload
```

This is an implementation validation; neither smoke with 2 physical
crossings nor the fixed-point pilot100 constitutes a measured
muon-angle-dependent detector asymmetry.


## 28 Oct 2026: C++ compilation failed — installed Sensor::GetMedium API (SOURCE PATCHED)

The local user's installed Garfield++ header explicitly declares:

```cpp
Garfield::Medium* Garfield::Sensor::GetMedium(double x, double y, double z);
```

The third attempt did not compile because the previous code incorrectly
used the *different* 4-argument API pattern `GetMedium(x,y,z,medium)`.
The subsequent `check2` correctly refused to use a stale executable,
and **no third B1 simulation was started**.

The constructor now calls exactly the installed 3-argument,
pointer-returning method and checks `medium && medium->IsIonisable()`.
Code review additionally found and fixed:
- `TrackHeed::NewTrack` is checked for returned false.
- The exact Garfield++ User Guide (section 5.1.1) states that when
  `EnableDeltaElectronTransport()` is enabled, `GetElectron`
  provides low-energy **conduction-electron positions**, but the
  returned kinetic energy and directions are not meaningful. The
  previous code incorrectly fed these numbers to `AvalancheElectron`.
  The corrected code seeds thermalized electrons at **0.1 eV** with
  zero initial direction (isotropic sampling), identically for all
  four detector configurations. This is a stated modelling approximation.
- A dedicated `--validate-heed` mode now runs **REAL TrackHeed**
  on the unchanged original two Geant4 muons, prints cluster and
  low-energy electron counts, and fails if all ionization is absent.
  **No avalanches are executed by that preflight.**
- `check2` first validates genuine C++ CSV parsing and the gas-medium
  pointer, then runs the Heed-only check. A successful Python synthetic
  test is no longer mistaken for successful Heed physics.
- `build` writes a **source SHA256 / executable SHA256** stamp only
  after CMake completes; `check2` verifies it. This handles noisy WSL
  file modification timestamps / clock-skew warnings without accepting
  an old executable by accident.

**Validation status:** the header signature is supported by the actual
compiler diagnostics and the Garfield++ manual. The changed files have
been reviewed through connected GitHub. **The current ChatGPT container
does not contain the user's installed Garfield++ libraries or headers**;
Google Drive is storage and does not execute compiled C++.
No successful compile or TrackHeed smoke has yet been observed for this
latest revision; this is an explicit unverified boundary, not a claim
of completed physical testing.

Safest immediate sequence: first compile, then Heed-only preflight.
Do not start full avalanches until BOTH pass:

```bash
cd ~/na6p/src/NA6PRoot-mwpc
git pull --ff-only
bash test/mwpc-stage-c/run_pi_native_muon_b1.sh build
bash test/mwpc-stage-c/run_pi_native_muon_b1.sh check2
```

If successful, separately run the already-created immutable
`smoke2retry2` mode, which remains unused after the previous
compiler failure. Both previous failed run folders remain unchanged.
