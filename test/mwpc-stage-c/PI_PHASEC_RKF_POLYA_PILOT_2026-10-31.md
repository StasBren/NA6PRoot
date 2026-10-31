# MWPC Phase-C: separate macroscopic RKF + Pólya *pilot* (2026-10-31)

**Status:** Source committed and tested syntactically against isolated Garfield API stubs; no live native Garfield++ compilation, gas-table generation, or physical RKF execution has yet occurred on the user's WSL. Neither a verified gain-speedup nor a calibrated replacement of the 100,000 Garfield microscopic avalanches has been demonstrated. No old physics data is overwritten.

## Physical scope

The pilot copies the **exact analytic geometry** of `phaseB_native_muon_3p3_avalanche.cxx`: Ar/CO₂ 70/30 at 293.15 K, 760 Torr; anodes 30 μm diameter at +1800 V, 3/4 mm pitch; optional field wires 75 μm diameter at 0 V; cathodes ±3 mm at 0 V; B=0; 9-anode finite patch, original mm→cm coordinate conversion and original D500 exact start coordinates. The original source C++ is NOT changed.

Unlike microscopic `AvalancheMicroscopic`, `DriftLineRKF` needs an **Ar/CO₂ Magboltz swarm transport gas table** (drift velocity, Townsend, attachment). The pilot generates an independent 40-log-grid (E=100..300000 V/cm) gasfile with `GenerateGasTable(ncoll)`; it never replaces the original microscopic gas cross-section model or assumes this is the identical physics process.

`DriftLineRKF::GetGain()` / `GetLoss()` provide a **macroscopic gain/loss estimate**. The separately sampled positive Gamma/Pólya random number with optional zero gate is labelled `polya_proxy_ion_births`; **it is an uncalibrated proxy, NOT Garfield microscopic `ion_births` nor a measured cathode signal**. `SetGainFluctuationsPolya` controls gain fluctuations in RKF *signal calculation* and must not be silently equated with the microscopic ion birth census. The proxy's default theta=0.4 and p0=0.17 are illustrative, motivated by old A20 distributions but NOT fitted to RKF per-start predictions. Further full 120-muon validation is required.

**Metal starts:** Exact anode and field-wire diameter geometry is used for physical zero start classification; other failed drift starts are recorded as **nonmetal transport failures**, never silently turned into physical zero avalanches. The program has explicit checks for nonfinite/excessive gains and invalid numerical input. It never creates new Geant4 or TrackHeed primaries.

## Actual first WSL attempt and CMake cache fix (user-confirmed 2026-10-31)

The user pulled through `ef8e8a3` and attempted `cmake --build test/mwpc-avalanche/build --target mwpc_phase_c_rkf_polya_pilot -j2`. Make returned `No rule to make target`. This indicates the existing build tree did not know the new CMake target; **no C++ compilation began**. The user separately ran `prepare`, which **PASSED**, creating the immutable 100 exact original D500 start plan. The later `buildcheck`, `gas`, `run` and `compare` modes failed because no RKF executable or output exists. This was a dependent build failure, not evidence for physical inaccuracies. Do **not** re-run `prepare`, wipe build caches, delete the new plan or modify old Garfield sources. The fix is `cmake -S test/mwpc-avalanche -B test/mwpc-avalanche/build` followed by building the target. Stop on any compiler error and inspect the actual user's installed Garfield API; the earlier mock-header syntax test is not native compilation. The user has not yet confirmed a successful reconfiguration/build.

## Safe first native WSL sequence

```bash
cd ~/na6p/src/NA6PRoot-mwpc
git pull --ff-only
python3 -m py_compile test/mwpc-stage-c/pi_run_rkf_polya_pilot.py

# The previous build directory predates the newly added CMake target.
# Reconfigure it explicitly. This reuses its existing cached Garfield settings;
# it does NOT delete existing binaries or change original Garfield physics.
cmake -S test/mwpc-avalanche -B test/mwpc-avalanche/build
cmake --build test/mwpc-avalanche/build --target mwpc_phase_c_rkf_polya_pilot -j2

# Verify 40 SHA-sealed A100 physical chunk manifests and prepare 100 fixed starts.
python3 -u test/mwpc-stage-c/pi_run_rkf_polya_pilot.py prepare
python3 -u test/mwpc-stage-c/pi_run_rkf_polya_pilot.py buildcheck

# Distinct Magboltz macroscopic table; do not replace any original A20/A100 gas files.
python3 -u test/mwpc-stage-c/pi_run_rkf_polya_pilot.py gas --ncoll 5

# 25 physical identical start coordinates per field, under separate output root.
python3 -u test/mwpc-stage-c/pi_run_rkf_polya_pilot.py run --theta 0.4 --pzero 0.17
python3 -u test/mwpc-stage-c/pi_run_rkf_polya_pilot.py compare
python3 test/mwpc-stage-c/pi_run_rkf_polya_pilot.py status
```

No `prepare` or `run` command is idempotent. Outputs are immutable and refuse overwrite. On native build failure, STOP and preserve the error; do not alter original archived microscopic source or `Garfield::` API semantics to force a superficially successful response. Fix the isolated pilot only after inspecting the actual installed Garfield headers.

The pilot output root is
`test_runs/mwpc_stage_c/pi_RKF_Polya_separate_B0_pilot_v1`.

## Validation gates

- `prepare` calls `validate_extra()` to check SHA and checkpoint integrity of all **40 completed A100 chunks**, original 84k A20 and original 120-muon full reference. It chooses 25 distinct target cells per geometry and one *unused original D500 start* per cell; selection uses deterministic randomization, **not the corresponding ion birth outcomes**.
- `buildcheck` invokes actual compiled C++ `--check` for all 100 pilot exact starts and validates source/unit/CSV schemas without Garfield drift.
- `gas` refuses overwrites and reports Magboltz gasfile SHA. Note `--ncoll 5` is a pilot gas table; numerical stability/accuracy must be studied versus additional gas-table grid/collision samples before making detector claims.
- `run` refuses overwrites, saves exact starts, RKF physical transport result, failures, sampled proxy, individual times and SHA-linked receipt. Other microscopic sources remain bytewise intact.
- `compare` joins exactly `(geometry,box_id,sample)` to 16k existing fully microscopic trials, checks their SHA. It reports preliminary means, failures and time; **do not treat a 100-seed comparison as a calibrated 120-muon predictive-interval result**.

If gasfile creation fails or RKF returns NaN/overflow, the underlying Garfield transport-table/settings or installed version need investigation (not a reason to silently fill missing responses with 0).

## Public reference

- Garfield++ documentation: https://garfieldpp.web.cern.ch/doxygen/UserGuide.pdf — drift-tube example shows `SetFieldGrid`, `GenerateGasTable`, `WriteGasFile`, `LoadGasFile`; Chapter 6 explains RKF gain and optional Pólya.
- Garfield++ maintainers explain macroscopic versus microscopic transport and need for a gasfile: https://root-forum.cern.ch/t/li-in-a-cylindrical-drifttube/62635
- Garfield++ Q&A: https://root-forum.cern.ch/t/help-on-garfield-avalanchemicroscopic-and-avalanchemc/53678

## Next physics decision

First require native WSL **compiled run + no unexpected nonmetal transport failures + sensible RKF gain and runtime**. Then estimate any calibration of RKF mean gain, Pólya theta and zero probability on a *separate calibration set*, and only then test **unseen independent** microscopic starts / full Geant4–Heed muons with matched definitions of ion counts, observable uncertainty, and speed. The old 120-muon positions were used in choosing the A100 50-cell targets, so their reanalysis is **source-conditioned**, not a blinded independent external test.
