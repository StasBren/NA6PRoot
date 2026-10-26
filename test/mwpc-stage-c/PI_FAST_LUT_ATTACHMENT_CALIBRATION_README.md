# Fast LUT attachment-hazard calibration: FIRST PASS, not production

## Rationale
The validated 5/2/1 um MC comparison on 14,400 matched source sites
showed minimal shifts of mean drift time (<0.1 ns), with 36-bin phi TV
differences of order 0.04-0.06, similar to finite-microscopic reference
noise. This does NOT establish exact equivalence, because the MC with
attachment disabled has a different surviving electron population.

The read-only Magboltz eta(E) probe showed high-field attachment
peaking near 7.19 /cm in the gas table and ~6.93 /cm in the sparse
field sample. Chamber sampled fields were ~0.4-22 kV/cm. A sparse
sampling is NOT an attachment probability or path field maximum;
very close to anode there are higher E values.

Next goal: test whether a SINGLE k in

    eta_eff(E) = k * eta_Magboltz(E)

can bring fast stochastic MC closer to the frozen Microscopic LUT.

## Crucial hazard subtlety

AvalancheMC stochastic stepping uses drift advance plus diffusion
displacement, so geometric chord length between points is *not*
necessarily equal to the nominal physical drift step. We record:
- hazard_eta_vdt = sum[ eta(E_mid) * |v_d(E_mid)| * dt ]
- hazard_eta_chord = sum[ eta(E_mid) * |delta_x_MC| ]

Both are approximate alternatives. The *vdt* hazard is the primary
calibration candidate; chord is a sensitivity diagnostic. Neither
has been demonstrated equal to the internal Garfield++ ComputeAlphaEta
and ComputeGainLoss integration. Before choosing the production model
we MUST run a native MC-attachment-on control at k=1 on the
same input population and compare predicted losses to actual losses.

We also record path steps, E sampled range, path lengths and the
additional CPU microseconds used to extract hazards. No microscopic
recalculation or gas table generation occurs.

## 1. Pull, synthetic tests, new CMake target

    cd ~/na6p/src/NA6PRoot-mwpc
    git pull --ff-only
    PY=test_runs/mwpc_stage_c/.venv_pi_slides/bin/python
    "$PY" -m py_compile \
      test/mwpc-stage-c/pi_fast_lut_attachment_fit.py \
      test/mwpc-stage-c/test_pi_fast_lut_attachment_fit.py
    "$PY" test/mwpc-stage-c/test_pi_fast_lut_attachment_fit.py

    cmake -S test/mwpc-stage-c -B test/mwpc-stage-c/build
    cmake --build test/mwpc-stage-c/build \
      --target pi_fast_lut_attachment_hazards -j2

    export ROOTSYS="$HOME/na6p/install/root"
    export LD_LIBRARY_PATH="$ROOTSYS/lib:$HOME/na6p/install/garfieldpp/lib:$HOME/na6p/install/garfieldpp/lib64:$LD_LIBRARY_PATH"

## 2. Un-killed stochastic MC drift with hazard extraction

Uses the *same* 5 um MC algorithm, exact same 14,400 selected birth
coordinates AND exact same RNG seeds as mc_no_attach_5um.csv. The
new output is separate and refuses overwrite. The analysis later
checks every endpoint and drift time against the previous 5um CSV
to detect accidental RNG or physics changes.

    BASE=test_runs/mwpc_stage_c/pi_fast_lut_reference_v1
    GAS=test_runs/mwpc_stage_c/pi_garfield_drift_methods_v1/ar70_co230_B0_transport.gas
    BIN=test/mwpc-stage-c/build/pi_fast_lut_attachment_hazards
    "$BIN" \
      --starts "$BASE/selected_micro_reference.csv" \
      --gas-file "$GAS" \
      --mc-attachment off \
      --mc-step-um 5 \
      --output "$BASE/mc_no_attach_5um_hazards.csv"

If a trajectory reaches an invalid field segment or its hazard cannot
be computed, STOP and inspect. Do not silently substitute H=0.

## 3. One-parameter fit, independently assess both geometries

    "$PY" test/mwpc-stage-c/pi_fast_lut_attachment_fit.py \
      --reference "$BASE/selected_micro_reference.csv" \
      --hazards "$BASE/mc_no_attach_5um_hazards.csv" \
      --baseline-mc "$BASE/mc_no_attach_5um.csv" \
      --out "$BASE/attachment_k_first_pass.json"

The tool computes 3 solutions for each hazard convention:
- k_joint: fit on both geometry FIT sets;
- k_noFW: fit using only noFW FIT starts;
- k_FW75: fit using only FW75 FIT starts.

The k values are then frozen and compared on both independent
VALIDATION sets. This reveals whether the factor transfers across
geometry. Global efficiency is fitted by construction, so the
important out-of-fit metrics are:
- validation global collection bias;
- 6x6 coarse-grid spatial efficiency residuals;
- 36-bin weighted phi angular TV;
- survival-weighted mean drift time;
- drift + hazard CPU/electron (not just drift CPU).

The fit is currently a deterministic expected survival
S(k,H)=exp(-kH). This is the correct way to optimize without MC
reruns, but a real production LUT must draw an independent
Bernoulli(S) sample per electron (separate pinned RNG stream)
to retain event fluctuations.

### Fundamental limitations

This is the FIRST exploratory calibration. A k that matches global
efficiency is NOT validated without spatial, angular, timing and
event-level checks. The reported normal-approximation confidence
intervals are not yet the final cell-clustered intervals.
The fast kernel is not to be silently promoted to production.

Before promotion we still need: k=1 comparison with native Garfield
attachment-on, sensitivity to Magboltz E-table resolution (currently
only 16 logarithmic field grid points), fit/validation reproducibility,
and the existing native Geant4/Heed event-level fold.
