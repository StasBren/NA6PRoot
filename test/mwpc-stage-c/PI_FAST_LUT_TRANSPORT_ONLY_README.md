# NA60+ MWPC: current active workstream is TRANSPORT ONLY

**Decision date: 2026-10-26.**
Do NOT work on electron attachment calibration in this workstream.
Do NOT alter prior microscopic LUT500, Geant4/Heed outputs, Magboltz
tables or the saved attachment calculations.

Objective: fastest model for electron drift-time and impact-angle
distributions around the anode at p=3 mm, B=0, Ar/CO2 70:30,
anode 1800 V, +/-3 mm cathode planes, with/without grounded
75 um field wire. Garfield::AvalancheMC, EnableAttachment(false)
and a preliminary 5 um step is the fast candidate.

Already produced:
- Frozen microscopic LUT500, both geometries, SHA256 verified.
- selected_micro_reference.csv, 7200 matched birth sites / geometry,
  split into 3600 fit and 3600 validation.
- mc_no_attach_5um.csv, 7199/7200 noFW anode endpoints,
  7200/7200 FW75. Native MC-only mean drift CPU 0.732 / 1.104 ms/e.
- mc_no_attach_2um.csv, mc_no_attach_1um.csv and noise-floor tests.
- Step 5 -> 1 um shifts mean MC drift by <0.1 ns; angular TV
  differences have the scale of reference finite statistics.
This does NOT yet establish agreement of absolute drift-time
distributions with the microscopic model.

## Do now: read-only detailed transport audit, no new Garfield runs

    cd ~/na6p/src/NA6PRoot-mwpc
    git pull --ff-only

    PY=test_runs/mwpc_stage_c/.venv_pi_slides/bin/python
    "$PY" -m py_compile \
      test/mwpc-stage-c/pi_fast_lut_transport_only_audit.py \
      test/mwpc-stage-c/test_pi_fast_lut_transport_only_audit.py

    "$PY" test/mwpc-stage-c/test_pi_fast_lut_transport_only_audit.py

    BASE=test_runs/mwpc_stage_c/pi_fast_lut_reference_v1
    "$PY" test/mwpc-stage-c/pi_fast_lut_transport_only_audit.py \
      --reference "$BASE/selected_micro_reference.csv" \
      --mc "$BASE/mc_no_attach_5um.csv" \
      --out "$BASE/transport_only_5um_detailed_v1.json"

Generated, and never silently overwritten:
- transport_only_5um_detailed_v1.json: fit AND validation separately,
  36-bin hit-angle PDFs, time mean/RMS and 10/50/90/99 percentiles,
  two-sample time CDF/KS, anode wire probabilities and 6x6 coarse
  birth-position maps. Every input includes its SHA256.
- transport_only_5um_detailed_v1.png: validation common-survivor
  angle PDFs, drift-time CDFs and 6x6 mean-time bias heatmaps.

Mean drift bias for common birth sites uses a 900-cell cluster
bootstrap (500 repeats) for an approximate 95% interval. This is a
diagnostic of response contrast, not proof of full physical equivalence.

## Essential selection caveat

The reference AvalancheMicroscopic had attachment active, so only
surviving electrons have a measured anode phi/time.

We therefore provide two DIFFERENT observables:

1. Own-survivor comparison: micro collected electrons vs MC collected
   electrons, with different survival populations.
2. Common-birth survivor comparison: choose the SAME initial birth
   coordinates for which BOTH micro and MC reach anode, and compare
   their phi/time distributions, not pair-by-pair stochastic
   collision histories.

Common-birth still selects microscopic tracks by attachment outcome.
Neither yields a truly attachment-free microscopic reference.

We do NOT compare the raw anode efficiencies 0.82/0.79 micro and ~1.0
MC-off as an error in fast electron drift: the physical scopes differ.

## Our next physics decision

After reviewing the time bias, 10/50/90/99 percentiles, KS distances,
wire-index PDF and spatial map:
- If transport metrics are stable and fit the project precision needs,
  use the fast-MC-off kernel only for drift and intrinsic angular
  transport observables (explicitly mark attachment=off in metadata).
- If substantial time or spatial bias exists, examine the Magboltz
  drift velocity/diffusion tables and local electric fields first,
  not an ad-hoc attachment factor.
- Before broad generalization or detector-response claims, obtain
  an independent small microscopic transport-only control, not just
  surviving microscopic electrons.

Deferred for later gain/charge modeling:
electron attachment, effective Townsend coefficient alpha - eta,
avalanche multiplication fluctuations, induced signal, electronics
threshold, per-event detected charge and RMS(A).

Do not silently substitute the attachment-free transport LUT into
the existing detector-efficiency/Geant4-Heed event-folding pipeline.
