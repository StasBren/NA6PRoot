# Fast Garfield LUT: reuse frozen microscopic LUT500 start positions

This is the executable Stage-0/Stage-2 branch of
PI_FAST_LUT_VALIDATION_ROADMAP.md. It does **not** replace the
production LUT. It never reruns microscopic Garfield/Geant4/Heed.

Both frozen 450,000-trial p3 microscopic CSVs are verified row-by-row,
including each full 30x30x500 grid and the immutable-source SHA256 from
the existing 100k-parent analysis. From each of 900 cells, four distinct
microscopic birth coordinates are chosen for **fit** and four distinct
ones for **validation**; total 14,400 starts across two geometries.

Fast MC uses these exact same physical starts and new independent
random seeds. The existing microscopic end points, status, angular
arrival, and time values remain the reference.

## 1. Update, test, compile (no physics run)

    cd ~/na6p/src/NA6PRoot-mwpc
    git pull --ff-only
    PY=test_runs/mwpc_stage_c/.venv_pi_slides/bin/python
    "$PY" -m py_compile test/mwpc-stage-c/pi_fast_lut_select_reference.py \
      test/mwpc-stage-c/pi_fast_lut_existing_starts_report.py \
      test/mwpc-stage-c/test_pi_fast_lut_existing_starts.py
    "$PY" test/mwpc-stage-c/test_pi_fast_lut_existing_starts.py
    cmake -S test/mwpc-stage-c -B test/mwpc-stage-c/build
    cmake --build test/mwpc-stage-c/build \
      --target pi_fast_lut_existing_starts_mc -j2

    export ROOTSYS="$HOME/na6p/install/root"
    export LD_LIBRARY_PATH="$ROOTSYS/lib:$HOME/na6p/install/garfieldpp/lib:$HOME/na6p/install/garfieldpp/lib64:$LD_LIBRARY_PATH"
    BIN=test/mwpc-stage-c/build/pi_fast_lut_existing_starts_mc
    ldd "$BIN" | grep 'not found' || echo "PASS: loader resolved"

## 2. Validate frozen micro LUTs, select reference births

    BASE=test_runs/mwpc_stage_c/pi_fast_lut_reference_v1
    "$PY" test/mwpc-stage-c/pi_fast_lut_select_reference.py \
      --fit-per-cell 4 --validation-per-cell 4 \
      --output "$BASE/selected_micro_reference.csv"

Expected counts: 450,000 complete source microscopic rows for each
geometry; 7,200 selected birth sites per geometry; 14,400 total.
Output: selected_micro_reference.csv and a SHA-verified
selected_micro_reference.manifest.json.

The selector intentionally fails on any frozen input-SHA mismatch,
incomplete 900x500 spatial grid, duplicate trial, or changed experiment
configuration. Do NOT bypass checks to force progress. Existing
result files are never silently overwritten.

## 3. One MC run, step=5um, attachment OFF

    GAS=test_runs/mwpc_stage_c/pi_garfield_drift_methods_v1/ar70_co230_B0_transport.gas
    "$BIN" --starts "$BASE/selected_micro_reference.csv" \
      --gas-file "$GAS" --mc-attachment off --mc-step-um 5 \
      --output "$BASE/mc_no_attach_5um.csv"

    "$PY" test/mwpc-stage-c/pi_fast_lut_existing_starts_report.py \
      --reference "$BASE/selected_micro_reference.csv" \
      --mc "$BASE/mc_no_attach_5um.csv" \
      --out "$BASE/mc_transport_step5_report.json"

This is a diagnostic of the FAST un-killed drift, NOT a realistic
survival/attachment model. The report compares the microscopic and MC
anode collection, 36 x 10deg phi histograms, arrival times,
mean MC drift CPU runtime, region-by-region biases, and independent
fit/validation summaries. Paired collection uncertainty is clustered
by LUT spatial cell.

**Important:** phi/time distributions are conditional on each
model's *own* collected electrons: they have different survivor sets.
Differences cannot be assigned to diffusion alone while the original
micro attachment is active. Also, four trials per region-cell are not
enough to demonstrate 2% PDF precision. This is step-convergence
triage, not final model certification.

## 4. Optional: confirm step convergence after 5um succeeds

    "$BIN" --starts "$BASE/selected_micro_reference.csv" \
      --gas-file "$GAS" --mc-attachment off --mc-step-um 2 \
      --output "$BASE/mc_no_attach_2um.csv"

    "$BIN" --starts "$BASE/selected_micro_reference.csv" \
      --gas-file "$GAS" --mc-attachment off --mc-step-um 1 \
      --output "$BASE/mc_no_attach_1um.csv"

    "$PY" test/mwpc-stage-c/pi_fast_lut_existing_starts_report.py \
      --reference "$BASE/selected_micro_reference.csv" \
      --mc "$BASE/mc_no_attach_5um.csv" \
           "$BASE/mc_no_attach_2um.csv" \
           "$BASE/mc_no_attach_1um.csv" \
      --out "$BASE/mc_transport_steps_report.json"

LUT source file names and hashes are defined in
pi_fast_lut_select_reference.py. If your old frozen files differ,
investigate why rather than changing the hashes.

## Later physics decision

The second phase will integrate the Magboltz attachment hazard
eta(E) along stored attachment-free MC paths, then fit a minimal k
factor using ONLY the fit set, and freeze it for validation.

If attachment-corrected MC fails phi, time or eventwise RMS(A)
requirements on both geometries, it will NOT be promoted to production,
regardless of raw speedup.
