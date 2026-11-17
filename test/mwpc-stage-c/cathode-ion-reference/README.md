# 2100 V uncapped cathode ion-birth reference — separate production target

**Status (2026-11-17): code committed; CMake/real Garfield++ compilation,
PBS physics runs and numerical convergence are NOT YET verified.**
Do not claim a validated induced strip current from these reference archives.

## Why this is a separate mode

Original immutable main campaign:
\`feature/p3-hv2100-uncapped-pbs1000-20261117\`,
\`test/mwpc-stage-c/hv2100-uncapped/\`, output
\`/storage/cgrp/stanis/mwpc_penning057/hv2100_uncapped_1000_v1/\`.
Its compiled 2100V binary and source/input SHA are not touched.

This branch adds:
- \`mwpc_p3_noFW_penning057_cathode_ion_reference_v1\`:
  isolated copy of the actual uncapped Garfield++ C++ source, using the
  same **real Geant4+TrackHeed** selection and seed derivation;
  \`DisableAvalancheSizeLimit()\` remains active. No new physics cap.
- \`CathodeIonArchive.hh\`: buffered ROOT TTree ion birth archive,
  independently seeded uniform-reservoir and complete callback streaming.
- \`verify_cathode_ion_root\`: ROOT/C++ physics QA, true per-Heed-primary
  callback census, physical-wire electron endpoint and ion sum comparisons.
- \`audit_geant4_1000_root\`: original 1000 input's station/channel/angle/
  phase/momentum/path distributions and honestly labeled ROOT plots.
- \`select_detailed_1000_root\`: **252** exactly stratified random
  Geant4-index choices with within-stratum sampling probabilities and
  **24** separate non-representative geometric stress cases.
- \`compare_sampling_to_full_root\`: per-primary reservoir K=256/512/1024
  versus exactly the same original muon with all callback births, reports
  finite-population sampling error and local anode birth associations.
- \`pbs_cathode_ion_reference.sh\` / \`manage_cathode_ion_reference.sh\`:
  separate PBS pilot and SHA-sealed results; no broad production submit.

## Physical modes

**Mode A:** the existing (unchanged) 1000-muon 2100V uncapped reference:
exact anode endpoints/gain/ion counts and sector moments. Ion-birth
reservoir is **zero**, so Mode A cannot generate individual cathode-strip
currents by itself. No simulated Mode A results should be claimed before
the existing PBS resource issue is resolved.

**Mode B:** \`weighted 256\`, \`weighted 512\`, \`weighted 1024\`.
Each Heed primary with N callback ion births retains exactly
\`n=min(N,K)\` ion records chosen by standard reservoir sampling without
replacement, with \`inclusion_probability=n/N\` and
\`inverse_inclusion_weight=N/n\` (when n>0). Original N is exact. A zero
ion birth primary has n=0, and is **not dropped** from the primaries tree.
The sampling \`std::mt19937_64\` is independent of the microscopic Garfield
random-number stream; it starts from a deterministic derived seed,
without requesting additional random numbers from Garfield.

**Mode C:** \`full 0\`. Every positive-ion creation callback is streamed
directly to compressed ROOT TTree \`ions\` with unit weight. No reservoir
limit/cap and no silent switching to sampling. Large files and 2100V
avalanches may exhaust PBS limits: these are FAILED, **unsealed** muons.
The pilot is only two original indices with at most one running at once.

Each ROOT archive represents **exactly one original Geant4 crossing** and
contains:
- \`muons\`: crossing identity, channel/station/det_id, native Geant4
  in/out coordinates and p-vector, momentum, path, angles, original
  Heed seed, physical anode phase/base and patch-to-chamber shift
- \`clusters\`: accepted TrackHeed cluster positions and times
- \`primaries\`: primary seed/cluster/position/time; actual collected
  electron endpoint count; exact callback births, number archived,
  selection weight and simulation successful flag
- \`ions\`: per-ion callback ordinal, Garfield seed, primary/cluster,
  local and physical chamber \`u,v,w\` [cm], creation time [ns], actual
  nearest physical anode, weight and inclusion probability
- \`electron_endpoints\`: all electron starting/ending points, times and
  endpoint statuses **only for first Heed primary in Mode C**
- \`electron_control\`: the first primary's \`0.05 ns\`-binned
  **whole-cathode-plane** prompt electron current control for ±v,
  not an individual physical stereo strip.

The \`microscopic_wire_seed.csv\` and \`microscopic_wire_muon.csv\`
from the independent detailed target retain all actual collected-electron
per-anode and exact ion-birth censuses; ROOT/C++ checker cross-validates.
\`microscopic_ion_birth_sample.csv\` is legacy-compatible and **empty** in
this independent reference; use the actual ROOT \`ions\` tree for ion
coordinates/times.

## Derivation of patch to chamber coordinate

Garfield analytic field uses \`x=u_patch\`, \`y=v_absolute\`, \`z=w_relative\`
(all cm), with nine local anodes \`k=-4..4\` at \`x=k*0.3cm\`.

Let \`u_in_mm\` be original chamber-local Geant4 entry coordinate.
The station has \`N_anodes = round(width/pitch)\` with
\`width=307.2mm\` for station 0 and \`512mm\` for stations 1..5,
\`pitch=3mm\`. Station 0 has even \`N_anodes\`: physical phase 1.5mm.
Stations 1..5 have odd \`N_anodes\`: phase 0mm.

Compute \`base=round((u_in_mm-phase_mm)/3mm)\`,
\`u_phase_cm=0.1*(u_in_mm-(phase_mm+3*base))\` and
\`u_wire0_cm=0.1*(phase_mm+3*base)\`.
For each birth:
\`u_chamber_cm=u_wire0_cm+u_patch_cm\`,
\`v_chamber_cm=v_absolute_cm\`,
\`w_chamber_cm=0.1*w_in_mm+w_relative_cm\`.
Its physical nearest-anode global index is \`base+round(u_patch_cm/0.3)\`.
Global wire center is \`phase_mm+3*(base+local)\` mm.
Checks reject inconsistent maps before a sealed physics result.

The callback does NOT label positive-ion gas species. The later
ion-transport layer must supply and document a mixture-dependent
mobility/ion-assignment assumption. Existing whole-plane reference is
not a physical stereo-strip weighting-field solution.

## Connect to WIPP and prepare only the independent reference

The previous **32 GB / 12 h** PBS pilot of the original 2100V campaign
was rejected with \`Job violates queue and/or server resource limits\`.
Do not misinterpret it as Garfield physics failure. Inspect actual queue
limits and choose an allowed queue/resources before submitting ANY new job.
Never run full uncapped avalanche computations on the login node.

\`\`\`bash
ssh -o ServerAliveInterval=60 stanis@wipp-an1.weizmann.ac.il

# Safe separate clone: don't update the pinned HV2100 or previous 1800V repos.
git clone --depth 1 \
 --branch feature/p3-2100v-cathode-ionbirth-reference-20261117 \
 https://github.com/StasBren/NA6PRoot.git \
 "$HOME/NA6PRoot-cathode-ion2100"
cd "$HOME/NA6PRoot-cathode-ion2100"

MANAGER=test/mwpc-stage-c/cathode-ion-reference/manage_cathode_ion_reference.sh
WORKER=test/mwpc-stage-c/cathode-ion-reference/pbs_cathode_ion_reference.sh
bash -n "$MANAGER" && bash -n "$WORKER"

# Can compile and audit the Geant4 input on the WIPP login; no physics runs.
bash "$MANAGER" prepare

# Inspect PBS queue and server rules before requesting 32 GB / 12 h.
qstat -Qf | grep -E '^Queue:|resources_max|resources_min|enabled|started'
qstat -Bf | grep -E 'resources_max|resources_default|max_.*run'
\`\`\`

If a queue genuinely accepts the requested 32 GB RAM and 12 hours,
configure it explicitly, e.g.:
\`\`\`bash
export MWPC_DETAIL_QUEUE='<verified-compatible-queue>'
export MWPC_DETAIL_MEM='32gb'
export MWPC_DETAIL_WALL='12:00:00'
bash "$MANAGER" pilot-weighted     # 2 actual Geant4 muons, K=256
\`\`\`
Pilot tasks create
\`/storage/cgrp/stanis/mwpc_penning057/cathode_2100_uncapped_v1/weighted_k256/muon_NNNN\`.
Only if they complete and resource use is acceptable:
\`\`\`bash
bash "$MANAGER" pilot-full        # two same muons, ALL callback ion births
bash "$MANAGER" weighted512       # same two muons, reservoir 512
bash "$MANAGER" weighted1024      # same two muons, reservoir 1024
bash "$MANAGER" status
\`\`\`
Each PBS task has one original muon; job 0 and job 1 are **not** guaranteed
representative across channel/station. Actual 252 stratified muon indices
are in:
\`.../geant4_input_audit/weighted_252_original_indices.csv\`.
No automatic 200+ job submit exists until resource pilot is audited.

Once the first full and sampled archives are actually SHA-sealed:
\`\`\`bash
ROOT=/storage/cgrp/stanis/mwpc_penning057/cathode_2100_uncapped_v1
CHECK=./test/mwpc-avalanche/build-cathode-2100/compare_sampling_to_full_root
for K in 256 512 1024;do
 "$CHECK" "$ROOT/full_k0/muon_0000/ion_births.root" \
  "$ROOT/weighted_k$K/muon_0000/ion_births.root"
done
\`\`\`
A completed per-muon worker writes:
\`input.csv\`, \`ion_births.root\`, \`microscopic_*.csv\`,
\`physics.log\`, \`root_audit.log\`, \`resource_time.txt\`,
\`SHA256SUMS\`, \`metadata.txt\`, \`complete.txt\` ONLY on success.
Failures leave an unsealed directory with logs/FAILED marker and cannot
be silently overwritten. The manager's \`status\` verifies SHA of every
completed task.

## Physical validation limits and next decisions

- ROOT checker requires exact ion birth callback and per-wire
  electron/ion sums for EVERY primary, physical anode phase consistency,
  K/N inclusion weights and first-primary direct electron control.
- \`compare_sampling_to_full_root\` calculates conditional true versus
  weighted spatial/time means, per-nearest-wire ion populations and
  explicit finite-population standard errors for matching full/sampled
  events. Those comparisons must be run and reviewed: **not yet passed**.
- The 252-event geometry-stratified selection is representative within
  explicitly recorded station/channel/signed-theta-u/anode-boundary strata.
  Separate 24 geometric stress cases are *not population-representative*.
  Heed multiplicity, gain-tail and actual multi-anode charge-sharing
  strata cannot be constructed from Geant4 input alone. They must be
  measured at 2100V **after a real accepted pilot**; phi/omega remain
  held-out for any Fast4/Fast5 donor model calibrated using J/psi.
- The original G4 selection does not freeze physical *w* board edges:
  \`edge_w_distance_mm=-1\` means UNKNOWN, not a true edge distance.
- The archive does NOT store a positive-ion mixture/species, run ion RKF,
  invent weighting-field geometry or calculate real stereo strip current.
  Those are separate next-stage Garfield++ tasks once mechanics/mobility
  and strip weighting fields are validated.
- Current model omits space charge, streamer/discharge, photon feedback.
  **Uncapped fixed-field 2100V gain is not proof of stable detector HV.**

In all follow-on plots use **ROOT/C++**, ATLAS-like style, honest
\`NA60+ Simulation\` labels, no Python physics/statistical analysis.
