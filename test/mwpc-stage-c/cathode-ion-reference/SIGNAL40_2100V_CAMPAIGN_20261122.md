# 2100 V MWPC / DiCE: phase-space-complete cathode signal development campaign

**Status, 2026-11-22:** Isolated C++/ROOT survey and two-muon PBS pilot code
are committed. No cluster jobs were submitted from this branch and the new
selector has not been compiled/run with the real 997 ROOT archives in this chat.
Do not confuse a code commit with a completed physics analysis.

## Why this is not another 1800 V run

- Frozen numerical geometry is P3 **no field wires**, nine analytic anodes,
  **3 mm pitch, 30 μm diameter, v=±3 mm**, Ar/CO₂ 70:30,
  Penning Ar r=0.57, B=0, **2100 V**, no Garfield avalanche-size limit.
- The *as-built* Prototype-3 documented elsewhere has a **4 mm pitch,
  asymmetric 2+4 mm gas gap and 1800 V**. The numerical geometry is not a
  validated hardware representation, even if its source calls it \`p3_noFW\`.
- Previous 1800 V 3-channel 30k Geant4 campaign was an earlier physics
  reference. No output, manifest, binary or active PBS job is modified here.
- The more recent 2100 V run already has **997 SHA-verified weighted K256
  full Garfield/TrackHeed muon ROOT archives out of 1000** as reported in
  2026-11-21 work. Three slots are **1 failed and 2 unsealed**, not physics.
  Actual archive paths and state must be verified on the machine before use.
- Source archive (WSL/laptop) expected by this new selector:
  \`$HOME/na6p/cathode-overnight-2100/test_runs/mwpc_stage_c/overnight_2100_weighted_k256_v1/\`.
  It uses actual Geant4 and Heed seeds, exact Garfield ion callback census,
  true collecting anode counts, K256 weighted positions/times and per-wire
  physical identifiers. This is already a 2100 V simulation.
- Therefore a blind re-run of those 997 at the same voltage wastes CPU and
  may threaten reproducibility. We FIRST select existing events and use the
  997 for statistics. New Garfield jobs must answer a missing-data question.

## Scientific event selection (C++/ROOT only)

The \`mwpc_survey_and_select_2100_signal\` target opens every candidate
\`muon_NNNN/ion_births.root\` and \`microscopic_wire_muon.csv\` and verifies
that exact ROOT sum of per-primary births/collection equals physical-wire
muon sums. It produces:

1. \`archive_inventory.csv\`: original slot/channel/station/detector,
   signed theta_u and theta_w, impact phase in 3 mm lattice,
   momentum, path, Heed primaries, ion births, actual collected electrons,
   positive-wire count, second-wire collected fraction, mean gain.
2. \`channel_station_coverage.csv\`: all 18 channel×station cells:
   accepted archive counts, impact phase, angle sign, multiwire (>5% second
   wire), high-gain decile. Note some cells may have sparse occupancy.
3. \`survey.root\`: ROOT histograms and a frozen \`selected40\` metadata tree.
4. \`selected40.csv\`: **20 new J/ψ development** and **20 sealed independent
   channel validation events** (10 φ + 10 ω). All are already simulated at
   2100 V as weighted K256. Greedy deterministic purpose-specific coverage
   rewards each station, det_id, anode-centered/transition/near-boundary
   phase, both theta signs, large angles, >5% two-wire sharing, high-gain
   decile, extremes of |w| and long gas path.
5. \`excluded_slots.txt\`: any unsealed/missing candidates.
6. \`previously_inspected_not_sealed.txt\`: 0301, 0408, 0446, 0435, 0511,
   0101 are *not* eligible for the sealed validation. Event 0301 was heavily
   examined for U/V development and is **debugging only**. These exclusions
   avoid post-selection information leakage.

**Critical statistical distinction:** the selected 40 events are a
purposively enriched DEVELOPMENT/validation diagnostic and do **not**
have known inclusion probabilities as an unbiased population sample. Compute
station/channel occupancy, gain distributions, RMS/mean, charge-sharing
fractions and uncertainties from all verified 997 archives, not from 40.
Validation events are held unopened for threshold/tuning/model choices.
No claim is made that 20/20 covers every rare joint condition; the
18-cell coverage table highlights gaps.

Geometry mapping: station-0 anode phase is +1.5 mm, stations 1–5 phase
zero. \`impact_phase_mm = u_in_mm - (round((u_in_mm-phase)/3)*3+phase)\`
in [-1.5, 1.5]. Distinguish patch-relative Garfield coordinates from
physical chamber/global-wire coordinates. \`w\` magnitude is only a proxy
for longitudinal extremes; true PCB/board-edge clearance requires final
mechanical geometry, not an invented bound.

## Direct use on the same WSL machine that holds the 997 archives

\`\`\`bash
cd "$HOME"
git clone --depth 1 \
  --branch feature/p3-2100v-phase-space-signal-campaign-20261122 \
  https://github.com/StasBren/NA6PRoot.git \
  "$HOME/NA6PRoot-signal40-2100"

cd "$HOME/NA6PRoot-signal40-2100"
RUNNER=test/mwpc-stage-c/cathode-ion-reference/run_signal40_selection.sh

# Uses local established ~/na6p Garfield++/ROOT env. ROOT/C++ only.
bash -n "$RUNNER"
bash "$RUNNER" prepare
bash "$RUNNER" survey

OUT="$PWD/test_runs/mwpc_stage_c/phase_space_2100_signal40_v1"
head -n 8 "$OUT/channel_station_coverage.csv"
head -n 8 "$OUT/selected40.csv"
cat "$OUT/selection_scope.txt"
cat "$OUT/FROZEN_SOURCE_SHA256.txt"
\`\`\`

The runner checks actual SHA256 checksums for all 40 selected ROOT archive
folders; the 997-pool discovery step checks previously reported seal-file
presence and exact electron/ion ROOT-vs-CSV counts, **not full SHA rehash
of all 997**. For a full population cryptographic re-audit, rerun the
existing campaign's \`status\` verifier separately.

## Physical cathode response scope and required retainable information

**Retain from all 997:** \`muons\`, \`primaries\`, \`ions\`, real
\`microscopic_wire_muon.csv\`, timestamps and physical/global anode IDs,
original Geant4 entrance/exit and momentum/angles, Heed cluster/primary,
Garfield RNG seeds, exact callback birth counts and sample inclusion
weights; keep each source crossing intact. K256 is sampled per Heed primary;
there is NO reconstruction of individual ions from global hemisphere A.
Preserve zero-gain and zero-ion primaries.

**For the selected 40 U/V studies:** response should use both U and V
electrode banks, a fixed **10 ns** initial sampling grid, full ion-drift
continuation to roughly **100 μs** for total integrated charge as needed,
a provisional **22-channel window (~11 strips per plane)** where a
verified field covers them; use **0–1000 ns** as the prompt diagnostic
window. Do NOT mistake 1000 ns for the complete slow-ion collection.
Keep per-strip waveforms, drift status, births, time origins and
physical hit position. Electron prompt signal and Heed *initial*
positive-ion contribution are distinct; the latter is not included
in avalanche-ion callbacks. No electronics response or final thresholds
until VMM settings and shaping are fixed.

**Do not conflate existing models:**
- Two infinite whole cathode planes are NOT separate physical U/V strips.
- The wire-free analytic strip model is NOT the chamber's real weighting
  field. A wire-inclusive 3D FEM U/V map exists as a candidate and has
  convergence issues in thin-wire field gradients. No claim that all
  997 waveforms have physically validated strip accuracy.
- Provisional positive ion mobility Ar+ in Ar is NOT a mixture-validated
  Ar/CO₂ 70:30 transport model.
- With fixed-field uncapped Garfield++ avalanches, space charge/discharge
  is not modeled. Numerical uncapped gain is not proof of stable operation.

An earlier phase-resolved 10ns whole-plane kernel matched the same
weighted-ion inputs at ~1.355% waveform L2 on one *development* muon;
the three-muon whole-plane check was ~4.81%. These measurements are
**not** physical stereo-strip uncertainties.

## Optional new 2100 V Garfield++ high-detail pilot (NO bulk jobs)

Purpose: on two previously unseen *development* J/ψ events selected above,
repeat **exact same Geant4 crossing and original Heed RNG seed**, but
record K512 weighted ion births per primary. Only if this small independent
pilot is accepted should complete (all-ion) archives be attempted on
exceptional events, followed by full-vs-K256/K512/K1024 waveform
comparisons. No cap, no synthetic replacement; failures are unsealed.

1. Transfer only the small \`selected40.csv\` from WSL to WIPP after SHA
   verification (do not transfer 997 large ROOT files just for scheduling).
2. On WIPP, create a **separate** clone of this branch, set a NEW output
   directory and run the inherited \`manage_cathode_ion_reference.sh prepare\`.
   It builds the identical uncapped Garfield++ physics source plus ROOT
   archives and freezes Geant4 and binary hashes.
3. Inspect \`qstat -Qf\` for permitted memory/walltime. The prior
   32GB/12h request was rejected; avoid blindly resubmitting. Set
   explicit \`MWPC_DETAIL_QUEUE\`, \`MWPC_DETAIL_MEM\`,
   \`MWPC_DETAIL_WALL\` to supported limits.
4. \`submit_signal40_2100_pilot.sh\` submits **only an array 0-1%1**,
   maps those array positions to the two frozen original Geant4 development
   indices and calls the inherited checked PBS worker. No 1000-run and
   no touching historical outputs.

For example, after preparing and copying the selection file:

\`\`\`bash
ROOT=/storage/cgrp/stanis/mwpc_penning057/phase_space_2100_v1
export MWPC_DETAIL_OUTPUT_DIR="$ROOT"
export MWPC_SIGNAL40_SELECTION="$HOME/selected40.csv"
# Set these only to values confirmed accepted by qstat -Qf:
export MWPC_DETAIL_QUEUE='<queue accepted for memory and time>'
export MWPC_DETAIL_MEM='<e.g. 16gb if allowed and scientifically tolerable>'
export MWPC_DETAIL_WALL='<e.g. 08:00:00 if allowed>'
bash test/mwpc-stage-c/cathode-ion-reference/submit_signal40_2100_pilot.sh
\`\`\`

No new calculation is performed by this command sequence until PBS has
accepted the explicit pilot, which we have not attempted here.
The PBS worker writes \`complete.txt\` only after validating all ROOT
and per-wire exact electron/ion sums, preserved true ion births,
original event identity and SHA256. Memory/time failures are FAILED
unsealed data, not truncated physical gains.

## When would a larger 2100 V run make scientific sense?

The old 1800 V run reached around 130–138 twenty-muon blocks per channel
(approximately 8.1k total). Repeating those thousands blindly at 2100 V
is expensive and unnecessary for the first response test because
997 *already simulated at 2100 V* supply broad summary statistics and
a physics-accurate weighted-ion starting population.

A larger 2100 V Geant4 crossing campaign should be an explicitly
conditional expansion for **missing channel×station×angle×phase×multiwire
cells or reconstruction/ghost/efficiency tails**, with a documented
~10^4-crossing truth goal if the reconstruction validation actually
requires it. A larger fast library must not be built by 100-fold
recycling each of 997 donors; use held-out-conditioned generative
validation. This expanded step is NOT authorized/submitted by the tiny
2-muon PBS pilot scripts.

## Provenance / checks

- Commit and input-selection hash frozen before any new physical job.
- Original K256 archives remain read-only.
- Physics analysis and plots: C++/ROOT, honest **NA60+ Simulation**.
- Per-waveform signal validation: electron+ion, 10ns bins, U/V
  multiplicity and sharing, peak/time, local strip residuals, 0–1000ns
  early charge and full late-tail integrated charge and completion status.
- Population metrics from all sealed archives; stage-specific
  test subset cannot be used to optimize strip pitch, weighting geometry,
  ion mobility, interpolation kernels, electronics, and then claimed as
  independent validation.
