# Next stage: implement and validate physical cathode-strip response

**Status:** planned. This document explicitly distinguishes *existing prototypes* from validated strip-signal physics. No background simulation is assumed.

## M0 — baseline freeze, voltage decision and input audit

- Start from known 1800 V P3 noFW geometry as a reproducible reference, but **do not select the new cathode-response HV yet**. Current laptop/farm scan is at 1900, 1950, 2000, 2050 V.
- Wait for actual ROOT/C++ voltage-scan results: gain distributions including zeros, RMS(g)/mean(g), RMS(g|g>0)/mean(g|g>0), RMS(G_mu)/mean(G_mu), signed A distributions, capping and runtime. Decide a usable voltage only once results and physical-model limitations (space-charge absence) are understood.
- Freeze the Geant4 crossing / Heed seed and Garfield++ code/config/versions with immutable sha and manifest. Keep an independently held-out validation sample. Compare all six station/detector occupancies and geometry before mixing samples.
- Output a short, frozen run/physics summary, not just an arbitrary PDF plot.
- **No change to ongoing high-voltage laptop run or archived 1800 V reference.**

## M1 — prove the electric field and strip weighting fields

**Existing code to inspect (not assume production-ready):**
- \`test/mwpc-avalanche/phaseB_strip_weighting.cxx\` : ideal strip weighting prototype. Earlier configuration used a strip family *perpendicular* to anode wires (segmented in w); it is **NOT** automatically the near-parallel stereo strip geometry specified in group guidance.
- \`test/mwpc-avalanche/phaseC_p3_noFW_two_cathode_ion_reference_v1.cxx\`: actual P3 noFW drift field, ion RKF and two WHOLE cathode-plane signals, **not individual stereo strips**.
- \`test/mwpc-avalanche/phaseB_real_avalanche_signal.cxx\`: exploratory microscopic ion-birth capture, single-event ion drift and strip signals; inspect simplifications before adapting.

Technical steps:
1. Derive upper and lower cathode electrode geometries from current \`NA6PMWPCChamber\` and team input. Parametrize strip width/pitch/gap, both orientations, finite boundaries, board details and station variations. Provisional width 1.7–2 mm and tan(alpha)=0.1 to wire; do not assume upper/lower angles equal or opposite without proof.
2. Inspect whether Garfield's analytic strip weighting models support the necessary arbitrary stereo angle. If not, use a verified geometric transform or 3D field-map approach, **not a fake strip label or unsupported rotation**. Explicitly document why the approach is valid.
3. Distinguish real E-field (\`anodeHV\` / actual cathode biases) from weighting field (\`Phi_w=1\` on one readout electrode, \`Phi_w=0\` elsewhere). For analytic Garfield cells, auxiliary non-degenerate plane voltages may be required for \`CellCheck\`; these **do not represent physical detector voltages**.
4. Validate \`Phi_w\` at the selected strip/adjacent strips/cathodes and \`E_w=-grad(Phi_w)\`. Test sign, symmetry and readout-label handling. Repeat at off-center u/w, both cathodes, near the strip gap and board edge. Check partition of weighting potential where applicable.
5. Make a ROOT plot of the 2D/3D weighting-field response and its strip dependence. **No Python plots or computations.**

**Deliverable:** \`weighting_field_checks.root\`, ROOT-generated PNG/PDF plots, CSV of electrode geometry and a C++ unit test with explicit expected symmetries and numerical convergence.

## M2 — one-ion transport / current sanity test

1. Load valid positive-ion mobility for the Ar/CO2 gas or explicitly document and bound the effect of a temporary approximation. The existing two-cathode Phase C code labels its Ar+ in Ar mobility approximation *provisional*, so do not claim mixture-validated ion drift times.
2. Use \`Garfield::DriftLineRKF\` in the real electric field with \`EnableSignalCalculation(true)\`; sample the induced current on a **physical strip** using its weighting field; keep electron and ion signals distinct.
3. Check the Shockley–Ramo convention \`i_c(t) = -q v(t) dot E_weight,c(r(t))\` with consistent polarity. Integrated charge must match \`-q(Phi_w(end)-Phi_w(start))\` according to chosen current sign convention.
4. Check dt (e.g. 1,2,5,10 ns if realistic), numerical RKF step, endpoint status, cathode arrival, finite recording window and integrated current. If the ion drifts beyond the record window, label and quantify missing charge; never silently truncate.
5. Test both readout sides and several ion start azimuths / w positions. Verify expected sign/geometry/symmetry and charge-sharing behavior.

**Deliverable:** reproducible Garfield++ C++ single-ion benchmark plus ROOT-plotted \`i_c(t)\`, integrated charge and uncertainty/convergence tests.

## M3 — many ions from a microscopic avalanche

- First create a **small direct reference** containing individual ion-birth positions \`(u,v,w,t)\` from FULL Garfield++ ionization callbacks at the chosen HV; save sampled ROOT/CSV with provenance. The 1000-muon wire truth provides only an <=128-point reservoir per seed, and the HV scan exports none. Neither is a complete high-fidelity ion cloud.
- Use direct per-ion RKF and weighting field for a feasible small reference. Track every ion start / time and both cathode polarities.
- Prototype unbiased **weighted macroions** using sampled ion births (weight = exact births / retained sample) and/or precomputed signal transfer kernels over \`(u,v,w,t)\` and wire geometry. Compare against direct per-ion reference for K=32,64,128,256,512,1024 (where feasible); report waveform metric convergence and errors. Never assume 128 is enough.
- Verify that weighting and adding macroion currents is valid only within the linear fixed-field/no-space-charge approximation.
- Preserve correlations among ion birth angle, longitudinal position, time and gain; different anodes must not be merged by global asymmetry \`A_mu\`.

**Deliverable:** microscopic vs macroion/kernel overlay plots and quantified errors in integrated strip charge, peak current, early charge, rise time, and hit position.

## M4 — complete single-muon cathode signal

- Assign actual collecting-wire IDs based on Garfield electron endpoints. Maintain separate ion-birth association to nearest wire (those counts need not match).
- Sum all primary-avalanche contributions with original \`t_Heed + t_ionBirth + t_ionDrift\` or correctly defined absolute times. Use correct per-wire u/w geometry and upper/lower stereo segmentation.
- Include prompt electron contribution or explicitly quantify conditions in which it can be neglected.
- Evaluate for one-anode and multi-anode (>5% second-wire collected electron charge) cases. The 1000-wire-truth sample found 192/1000 multi-anode, and 143/192 with opposite local A signs; **this is a validation stress test, not a user input to the simulator**.
- Compare with a direct Garfield full-event signal on a manageable sample, especially at large |theta_u| and near wire boundaries.

**Deliverable:** per-strip current \`i_c(t)\`, induced charge, event-wise wire and strip geometry metadata, ROOT plots of waveforms, and comparisons per angle/station.

## M5 — front-end electronics and reconstruction interface

- Implement configurable linear electronics shaping \`V_c(t)=(i_c*h_c)(t)+noise\`, pedestal, threshold, digitization/ADC/TDC and realistic dead-time **only after** specification is provided. Never invent electronics parameters.
- Compare strip multiplicity, peak amplitude, integrated charge, time-over-threshold, hit resolution and coordinate residuals vs true muon crossing, split by station, angle and number of anodes.
- Interface with \`NA6PRoot\` /\`NA6PMuonSpec\`/\`NA6PMWPCChamber\` preserving local-to-global and det_id metadata; keep optimized geometry configurable because of ongoing mechanical design decisions.

**Deliverable:** C++/ROOT integration tests and documented mock-signal vs reconstruction checks.

## M6 — fast inference, optionally Fast5

- **First** produce validated Garfield reference strip currents. Then test Fast4-generated ion clouds as sources of the same signal chain.
- Alternative Fast5-L idea: donor-muon **complete multi-strip waveform** or current-template library, conditioned on impact phase \`u mod 3 mm\`, signed theta_u/theta_w, gas path, station, momentum (only if shown useful), and perhaps retained Heed multiplicity. Preserve whole-event inter-wire correlations; avoid reconstructing A_k from a single global A_mu.
- Donor sets for fit vs held-out test must be disjoint, with special test sample at large angles and edge positions. Use whole-muon bootstrap.
- Previous Fast5-T prototype was ~3.4x faster than Fast4 **response-only** in one benchmark, but had inaccurate local A correlations. Do **not** claim 10x is achieved, nor that response-only speed is detector end-to-end speed.
- Publish runtime/accuracy trade-off on the same laptop/CPU and code configuration after including signal generation and electronics.

## Acceptance checklist (must be demonstrated, not asserted)

- One-ion current obeys weighting-potential integral within predeclared numerical tolerance.
- Correct multi-electrode geometry, sign, time units and time integration.
- Macroion/sampling convergence relative to true direct-ion reference; report random and systematic uncertainty.
- Realistic per-wire Q/A/ion-birth and multi-anode correlations preserved; no global-to-local asymmetry fabrication.
- Per-station, signed incident-angle, charge-sharing (>5%), initial u phase and production-channel comparisons.
- ROOT-only analysis, ATLAS-like styling, clearly marked \`NA60+ Simulation\`, axes/units/uncertainties; no Gaussian assumption for gain.
- Updated GitHub source, reproducible runs, manifest, README and commands, held-out validation, known limitations.
