# Isolated p=3 mm Garfield drift-method benchmark

Objective: find the CHEAPEST response model that reproduces the
observable electron-anode collection, phi_hit, drift-time distributions
needed for the NA60+ MWPC LUT. This is independent of the previous
20-/500-trial frozen Garfield maps and all native Geant4/Heed folds.

## Current reference / physics

Same electric geometry as original c_transport_box_map.cxx:
- Ar/CO2 70:30, T=293.15 K, P=760 Torr, no magnetic field.
- Anodes at 3 mm pitch, 30 um diameter, +1800 V.
- Cathode planes v=-3 and +3 mm, 0 V.
- With or without grounded 75 um-diameter field wires at half pitch.
- Finite analytic field lattice of 9 anode wires, and matching 8 FW
  in the FW configuration, with periodic-cell starts near the center.
- Exact old anode hit criterion: a drift endpoint is closer than 1.5
  anode radii to a finite anode, with Garfield status != -7.
- Single initial electron; NO avalanche multiplication, strip
  induction, discriminator, or electronics.

Three algorithms are invoked independently on exactly identical
birth coordinates in the same field and gas:
1. AvalancheMicroscopic::DriftElectron: individual electron-molecule
   collisions. This is our original LUT physics.
2. AvalancheMC::DriftElectron: macroscopic drift velocity with
   stochastic diffusion, from a precomputed Magboltz gas table.
3. DriftLineRKF::DriftElectron with avalanche disabled: mean
   deterministic drift path, diffusion OFF, therefore useful as a
   qualitative field-line baseline but not an eventwise microscopic
   response model.

The MC and RKF methods need a precomputed Magboltz swarm transport
table; it is NOT produced by the old microscopic Garfield maps.
The gas table must cover electric fields reaching ~100s kV/cm near
the 15 um-radius anodes, not just drift-gap E. The small benchmark
starts with 16 logarithmic field values between 100 and 500,000 V/cm,
with two times 10^7 Magboltz collisions per point by default.
THIS grid and the stepping size need separate convergence tests.
Do not equate a coarse-grid pilot with a precision gas transport model.
Gas-table preparation can be expensive and is measured separately.

## WSL build

Only a new target is built; do not reconfigure a production
Geant4 tree or modify c_transport_box_map.cxx.

    cd ~/na6p/src/NA6PRoot-mwpc
    git pull --ff-only
    cmake -S test/mwpc-stage-c -B test/mwpc-stage-c/build
    cmake --build test/mwpc-stage-c/build --target pi_garfield_drift_methods_benchmark -j2

If GarfieldConfig.cmake is not found, use the same CMAKE_PREFIX_PATH
as the existing Stage-C transport build. Do not install a different
Garfield version to work around missing environment variables.

## Garfield runtime shared libraries (WSL)

CMake can compile and link this executable while the runtime ELF loader
still cannot locate libTree.so from the independently installed ROOT.
The working Stage-C pitch pilot explicitly resolves ROOT and Garfield libs.

After building, in the **same WSL terminal** as the subsequent prepare/run:

    export ROOTSYS="$HOME/na6p/install/root"
    export LD_LIBRARY_PATH="$ROOTSYS/lib:$HOME/na6p/install/garfieldpp/lib:$HOME/na6p/install/garfieldpp/lib64:${LD_LIBRARY_PATH:-}"

    BIN=test/mwpc-stage-c/build/pi_garfield_drift_methods_benchmark
    ldd "$BIN" | grep 'not found' || echo "PASS: all libraries found"

Do not run prepare unless the loader reports no missing shared libraries.
If the ROOT install location differs, resolve libTree.so from the same
ROOT prefix that built Garfield; do not mix multiple ROOT installations.
These exports apply only to the current shell and need repeating in a
new terminal. See run_pi_pitch_fieldwire_pilot.sh for the validated
dynamic ROOT library discovery procedure.

## Prepare ONE independent gas table

    BASE=test_runs/mwpc_stage_c/pi_garfield_drift_methods_v1
    mkdir -p "$BASE"
    BIN=test/mwpc-stage-c/build/pi_garfield_drift_methods_benchmark

    "$BIN" --mode prepare --gas-file "$BASE/ar70_co230_B0_transport.gas" --ncoll 2

This writes a .gas transport COEFFICIENT INPUT, NOT a response LUT.
It refuses to overwrite an existing gas file. Reuse it for all
later method/step comparisons. This avoids pretending the gas-table
precomputation was free. You may opt for --ncoll 1 as a purely
diagnostic (LOWER-QUALITY) first compilation/runtime test. Such a
table is not an accurate benchmark for final physics conclusions.

## Same-start pilot: 5 um MC steps

    "$BIN" --mode run --gas-file "$BASE/ar70_co230_B0_transport.gas" \
      --trials-per-box 4 --mc-step-um 5 \
      --output "$BASE/matched_5um.csv"

    PY=test_runs/mwpc_stage_c/.venv_pi_slides/bin/python
    "$PY" test/mwpc-stage-c/pi_garfield_drift_methods_report.py \
      "$BASE/matched_5um.csv" --out "$BASE/report_5um"

The sparse test uses 25 spatial boxes (5 u locations x 5 v
locations) with 4 independent starts per box: 100 physical starting
points, each passed to each of THREE algorithms, for EACH of TWO
geometries (600 drifts). This is not a substitute for the 900x500
full resolution map. The two gas geometries are run separately.

The CSV is preserved rather than replaced. Every result includes
anode/FW/other collection classification, impact phi, drift time
and only timed Garfield drift invocation. Program startup, gas-table
loading and gas-table GENERATION are excluded from transport-speed
ratios. They matter for total workflow runtime and must also be
tracked, especially for small LUTs.

The Python report validates identical start coordinates and computes:
- Anode collection fractions for both gas/electrode geometries.
- Mean and RMS of anode-collected electron drift times.
- Coarse 12-sector phi distributions (30 degrees per sector), with
  PNG overlay and total variation relative to microscopic reference.
- Mean and median CPU-wall use per electron, speedup from measured
  per-electron transport time only.

## Diagnose 5 um collection discrepancy BEFORE judging AvalancheMC

Our first 100-start/geometry numbers gave very different anode fractions:
microscopic 0.91 / 0.83, MC 0.59 / 0.57, RKF 1.00 / 1.00.
This is far too large to explain by the ~55-60x measured MC speed.
We must distinguish physically lost electrons from Garfield endpoint and
step-size classification, before replacing any LUT.

Existing matched_5um.csv ALREADY contains the exit status, endpoint
coordinates and original proximity-based anode/FW/other flags. The
read-only diagnostic below prints all endpoint status counts and
distances to actual anodes, field wires and cathodes, plus a same-start
comparison against the microscopic results. It does NOT relabel hits.

    cd ~/na6p/src/NA6PRoot-mwpc
    git pull --ff-only
    BASE=test_runs/mwpc_stage_c/pi_garfield_drift_methods_v1
    python3 -m py_compile test/mwpc-stage-c/pi_garfield_drift_methods_diagnose.py
    python3 test/mwpc-stage-c/pi_garfield_drift_methods_diagnose.py \
      "$BASE/matched_5um.csv" --out "$BASE/diagnose_5um.json"

If the MC endpoints cluster just beyond the 1.5-anode-radius criterion
(22.5 um from the wire center), that suggests an endpoint/step geometry
artifact. Many status -3 or missing endpoints instead suggest a numerical
or transport-table problem. Endpoints by the cathodes or field wires
indicate genuine different trajectories/collection basins. No diagnosis
can be concluded from status code alone.

## Why the 5 um MC efficiency was low: status -7 = ATTACHMENT

Diagnostic results from 100 matched starts/config:

| Geometry | Microscopic -7 (attachment) | AvalancheMC -7 | RKF -7 |
|---|---:|---:|---:|
| No FW | 9/100 | 41/100 | 0/100 |
| FW 75 um | 17/100 | 43/100 | 0/100 |

Garfield's -7 endpoint means **electron attachment** (capture in the gas),
NOT an out-of-bounds drift and NOT a failure of our 1.5-wire-radius
geometrical hit criterion. See the Garfield++ User Guide termination codes
and developer explanation at:
https://root-forum.cern.ch/t/help-on-garfield-avalanchemicroscopic-and-avalanchemc/53678

The RKF solver integrates the deterministic mean drift path and does NOT
sample single-electron attachment; its 100% collection is therefore
NOT a direct prediction for actual electron survival.

### Diagnostic physics intervention: disable MC attachment ONLY

DO NOT alter the frozen microscopic LUT or original matched_5um.csv.
Compile the benchmark new target after pulling the code:

    git pull --ff-only
    cmake --build test/mwpc-stage-c/build \
      --target pi_garfield_drift_methods_benchmark -j2

Then run the SAME 100 matched start sites/config, same 5 um stepping and
gas table, with only the macroscopic MC attachment process switched OFF:

    BASE=test_runs/mwpc_stage_c/pi_garfield_drift_methods_v1
    BIN=test/mwpc-stage-c/build/pi_garfield_drift_methods_benchmark
    GAS="$BASE/ar70_co230_B0_transport.gas"
    "$BIN" --mode run --gas-file "$GAS" \
      --mc-attachment off --trials-per-box 4 --mc-step-um 5 \
      --output "$BASE/matched_5um_no_attachment.csv"

    PY=test_runs/mwpc_stage_c/.venv_pi_slides/bin/python
    "$PY" test/mwpc-stage-c/pi_garfield_drift_methods_report.py \
      "$BASE/matched_5um_no_attachment.csv" \
      --out "$BASE/report_5um_no_attachment"
    python3 test/mwpc-stage-c/pi_garfield_drift_methods_diagnose.py \
      "$BASE/matched_5um_no_attachment.csv" \
      --out "$BASE/diagnose_5um_no_attachment.json"

The new CSV records `mc_attachment=off`. Reporting remains compatible
with the previous on-by-default CSV, which did not record this flag.
The benchmark never silently overwrites gas tables or result CSV files.
This test isolates the role of tabulated attachment. **MC with attachment
OFF is an unphysical diagnostic control for Ar/CO2; do not use its
collection efficiency to claim model agreement.**

If the MC collected-anode fraction grows to near the RKF baseline, the
dominant disagreement is the attachment modelling, not the geometry
classification. We then need to examine the Magboltz attachment
coefficient eta(E) and its coarse high-field interpolation separately.
If it remains low, examine diffusion/step numerics and stopping positions.

## Step-size numerical convergence test

Only after first 5 um run and diagnostic complete, run:

    "$BIN" --mode run --gas-file "$BASE/ar70_co230_B0_transport.gas" \
      --trials-per-box 4 --mc-step-um 2 \
      --output "$BASE/matched_2um.csv"

    "$PY" test/mwpc-stage-c/pi_garfield_drift_methods_report.py \
      "$BASE/matched_2um.csv" --out "$BASE/report_2um"

The SAME random start sites are used at 2 and 5 um; stochastic
realizations and per-electron outcomes need not match. Compare
ensemble distributions with statistical uncertainty, not individual
trajectories. Also later test 1 um and the high-field Magboltz grid.

CAUTION: this is an exploratory pilot (100 starts/config), so
estimates of 30-degree angular PDFs are imprecise; a visually
reasonable plot alone does not validate MC or RKF. Increase to
20-50 starts per spatial box after confirming a stable gas table,
API compatibility and step-size convergence. Preserve those output
CSVs under separate filenames. The RKF method is deterministic
per start and is NOT a substitute for Garfield fluctuations.

No large-scale Garfield LUT generation, Geant4 or native Heed folds
are started by this benchmark.
