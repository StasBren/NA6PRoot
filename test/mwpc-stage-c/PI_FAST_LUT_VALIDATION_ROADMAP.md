# NA60+ MWPC — fast Garfield LUT validation roadmap

Status: **research plan, not a validated replacement**. Current branch
feature/mwpc-stage-c. Frozen 20- and 500-trial microscopic LUTs and the
Geant4/Heed native folds remain immutable.

## Decision and use case

Goal: find the lowest-cost B=0, p=3 mm, Ar/CO2 70:30 transport model
that reproduces the microscopic single-electron collection, azimuth
phi_hit, arrival-time distributions AND resulting event-level up/down
asymmetry. Aim for a substantial end-to-end LUT generation speedup.

**Shortlist:** AvalancheMC (stochastic diffusion) with a carefully
calibrated attachment mechanism. DriftLineRKF stays a deterministic
field-line control, not an event-fluctuation surrogate. The old
AvalancheMicroscopic is the reference. Old pilot timings:
micro 54-69 ms/electron; MC 0.88-1.26 ms/electron (attachment on);
RKF 4.39-7.27 ms/electron. A second run varied by multiple factors,
so benchmark larger homogeneous batches before making speed claims.

Current 100-start observations:
- Micro anode survival = 0.91 (noFW), 0.83 (FW75).
- MC native table attachment on = 0.59, 0.57; lost electron status -7.
- MC attachment off = 1.00, 1.00; same geometry and start coordinates.
The loss difference is caused by attachment in this small sample.
Neither bare attachment-off MC nor a factor calibrated on 100 starts
is a valid physics replacement.

## Stage 0 — immutable input baseline and scientific acceptance gates

Inputs: existing 2 x (900 spatial cells x 500 independent microscopic
electrons) CSVs at:
test_runs/mwpc_stage_c/pi_pitch_fieldwire_500_p3_v1/pilot/
  p3p0_anode_only/garfield_electrons.csv
  p3p0_field75p0um/garfield_electrons.csv

Validate the frozen_config.json identifiers and SHA before analysis.
These files include TRUE u_start_mm, v_start_mm, u_end_mm, v_end_mm,
wire_index, phi_hit_deg, drift_time_ns, and attachment/collection
statuses. Reuse these data; do not regenerate microscopic electrons.

Create a deterministic, stratified, fixed held-out split of original
Garfield start sites covering all 900 bins, anode-adjacent, field-wire-
adjacent, both gas halves, and central-gap regions. Use separate
microscopic replicates/cells for fitting vs final validation; also
hold out an entire geometry in leave-one-geometry-out stress tests
when judging transferability. Same birth sites for comparisons do NOT
imply matched Monte Carlo collision histories.

Predeclare physics metrics and uncertainty estimators:
- Per-configuration and per-(u,v)-region anode collection and
  attachment probability, field-wire/cathode/missed endpoints.
- The 36-bin 10-degree normalized phi_hit PDF with error bands;
  harmonics and upper/lower fractions if useful.
- Conditional drift-time CDF and selected quantiles/mean/RMS;
  survivor selection must be treated consistently.
- Geant4/Heed-folded phi distributions and EVENT-level A mean/RMS,
  not only single-electron averages.
- Timing: runtime for a fixed throughput, excluding and including
  separate Magboltz table preparation, serialization, initialization.
Evaluate 95% uncertainty intervals using independent microscopic
trajectories per cell and (for events) clustered/parent bootstrap.
Do NOT count reused LUT lookups as independent physics.

**Provisional practical gate targets** (can be amended once *before*
the first calibration, but not retrospectively to make a model pass):
  global anode probability within 1 percentage point of microscopic,
  each high-occupancy spatial region within its statistically
  resolvable error band, 36-bin phi TV <= 0.02 (plus an explicit
  error-band comparison), mean drift-time bias <= 1 ns or <= 3%,
  and event-level RMS(A) error <= 0.015 (absolute).
All criteria must hold on BOTH noFW and FW75 held-out sets.
To call an effect established, account for finite LUT statistics and
paired/clustered sampling. Precision numbers are project targets,
not measured performance.

## Stage 1 — explain the current attachment discrepancy

Already added read-only programs:
- pi_garfield_attachment_field_probe.cxx
- pi_garfield_attachment_field_report.py

Run once on existing .gas table to inspect eta(E), alpha(E),
diffusion/drift and local electric field in both geometries.
Check the electric-field interval near 15 um anode radii (including
out-of-range interpolation). Fix any demonstrably broken gas-table
grid/configuration before attempting an empirical correction.
The table was prepared on 16 logarithmic E points, 100-500000 V/cm;
test finer high-field grids only if the observed field distribution
and eta(E) require it. Benchmark gas-table preparation cost separately.

Output: attachment_eta_vs_field.png and JSON, plus a succinct
diagnosis: transport-table issue / numerical step issue /
physical discrepancy not explained yet.

## Stage 2 — isolate un-killed stochastic transport

Run MC WITH attachment disabled on a stratified 2-10k-site subset
taken from the existing microscopic 500-trial CSV (do not rerun
micro). Record all endpoint coordinates, time and wire ID.
Test MC fixed distance step = 5, 2, 1 um; check convergence and
CPU dependence, especially near anodes/FW. Optional adaptive small
steps near wires if numerically justified.

Since attachment-off and microscopic-on have different survivor
ensembles, compare endpoint physics both unconditionally where
possible and in a separate attachment-adjusted validation. The
attachment control alone cannot validate the full 36-bin PDF.

GO/NOGO: if MC un-killed transport has seriously wrong destination
wire, phi-shape, or timing after step convergence, an attachment
scale factor CANNOT fix the model. Investigate diffusion, transport
table (velocity/DL/DT), field geometry; retain microscopic default.

## Stage 3 — minimal attachment correction; no brute-force scan

Let eta_tab(E) be the macroscopic attachment coefficient in 1/cm.
Start with a single global dimensionless multiplier k>=0:
    eta_eff(E) = k * eta_tab(E).

For an MC path with attachment disabled, store
    H_path = integral eta_tab(E(s)) ds
using the actual stochastic MC drift trajectory (NOT an arbitrary
constant length), and its terminal electrode/end coordinates.
Conditional survival approximation:
    S_k(path) = exp(-k * H_path).

Crucially: compute MANY candidate k using the SAME stored no-attach
MC trajectories + H_path, WITHOUT rerunning Garfield for each k.
For a discrete LUT, sample an independent Bernoulli(S_k) survival
for each electron with a reproducible distinct RNG stream.
Do not merely multiply event electron counts by S_k:
event-level fluctuations and phi-selection would be lost.

Fit k only on training cells/trajectories, e.g. grouped binomial
likelihood for microscopic reference attachment counts as a
function of start position and relevant observables. Use held-out
data to judge both geometries and spatial regions. A single k fit
to noFW MUST predict FW75 without refitting to claim
transferability; also report the reverse experiment.

If a global k fails, try a physically interpretable 2-regime
k(E) (low vs high electric field) and penalize complexity.
Do not freely fit per-cell factors then claim a portable gas model.
As a LAST RESORT, a geometry-specific smoothed survival map
P_survive(u,v) is an honest response surrogate for THAT geometry,
not a universal modification of the attachment coefficient.

Implementation choice: a custom post-drift Bernoulli survival
layer or an explicitly versioned scaled attachment medium; never
modify original .gas file and do not edit Garfield++ internals.
Check against an independently implemented attachment-on MC
control at k=1, because discrete numerical attachment algorithms
may differ from the continuous exponential survival approximation.

## Stage 4 — blind quantitative validation

After calibration, freeze k (or the minimal correction formula).
Evaluate on untouched held-out microscope sites AND on both
complete spatial response geometries (3 mm noFW / FW75).
Compare collection maps, field-wire attachment/misses,
36-bin phi and timing, and angular hemisphere asymmetry.

Important: successful reproduction of the 2 global collection
probabilities is insufficient. Rejection if conditional phi,
time, or high-occupancy start-position PDFs fail declared gates,
even when overall efficiency agrees.

Repeat a limited real J/psi MS0 native Geant4/Heed fold with the
fast model on the frozen input events and compare with the
previous microscopic-folded results. Evaluate eventwise A and
its RMS using correct sampled Poisson-Heed events, not analytic
expected Heed occupancy.

Output: one validation report with pass/fail for EACH predetermined
metric, confidence intervals and all input SHA/physics settings.

## Stage 5 — fastest correct LUT production strategy

Only after physics validation:
- Choose fastest MC step size that meets the validated tolerances.
- Exploit both exact B=0 reflectivity checks u->-u, v->-v.
  For the central symmetric field, a 30x30 spatial grid could need
  15x15 truly simulated cells (225) instead of 900. Explicitly
  check transformed wire indices, phi angles and edge semantics.
  Mirrored copies are CORRELATED, not additional trials.
  Check any future B!=0 configuration separately; do not assume
  both symmetries survive magnetic field.
- Explore adapting trials per cell to its Heed occupancy and
  target error (under-sampling an important rare region is forbidden).
- Benchmark actual production wall time for a fixed number of
  uncertainty-qualified independent microscopic/MC source trials,
  including initialization, table load, file IO and manifest checks.
- Report end-to-end speedup rather than multiplying nominal factors
  uncritically.

For the ALREADY simulated fixed geometries, a direct lookup in the
existing 500-trial microscopic LUT is the fastest operational choice.
A calibrated fast generator is mainly useful for NEW wire pitches,
field-wire dimensions, HV or gas settings; changes beyond validated
physics scope require renewed checks.

## Stage 6 — promotion and safe ongoing use

The approved fast-LUT generator gets a new explicit model identifier,
independent output namespace, CLI selection flag and automated tests.
Pin Garfield source/ABI, gas-table SHA, field configuration,
attachment model k/version, integration step, RNG streams and
calibration/validation sample hashes.

Run a 1-column smoke test -> 900-cell small control -> full LUT
only after the frozen acceptance report PASS.
Retain original microscopic LUTs as immutable reference.
No production Geant4/Heed job silently changes its transport kernel.

If no candidate passes ALL physics gates, publish a clear
negative result and keep microscopic LUT500. Speed never justifies
silently breaking the physical observables.

## Recommended immediate next action

Complete the EXISTING attachment eta(E) probe (no new trajectories).
After reviewing it, implement a tiny stratified-start extractor
from frozen microscopic CSV and run MC without attachment at
5/2/1 um on those existing reference birth positions.
Only THEN implement/calibrate k on the saved no-attachment paths.
