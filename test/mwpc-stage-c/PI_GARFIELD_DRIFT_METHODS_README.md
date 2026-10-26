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
