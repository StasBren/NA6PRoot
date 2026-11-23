# Claude review handoff — 1800V native per-wire charge/birth distributions
2026-11-23, NA60+/DiCE MWPC, `StasBren/NA6PRoot`

## Working agreement

ChatGPT implements independent analysis tools on GitHub. The user runs
them at WIPP and returns real logs/ROOT summaries/PNG figures. ChatGPT
interprets results and decides whether a focused Claude code review or
scientific cross-check is needed. Do **not** invent Garfield/ROOT results,
alter running PBS arrays, modify 1800V production source/output, or
pretend the physical U/V strip solver is already calibrated.

## Source of truth

New analysis branch:
`feature/p3-1800v-wiretruth-8k-analysis-20261123`.

Use ONLY existing *fully completed, receipt-backed* per-wire simulation
from `/storage/cgrp/stanis/mwpc_penning057_wiretruth30k_v1/results/`
(jpsi, phi, omega).
Source physical production commit `599e29522edb6fcec9dfcca2904308cb3186edce`
and binary sha256
`2cf92cb27d3a180cd26d3ab5417f57293473864f6d685ee43bacf931712aa8f1`.
Frozen original 20-muon Geant4 input for chunk N is
`/storage/cgrp/stanis/mwpc_penning057/original2m_CH_10000_v1/chunks/chunk_NNNN.csv`.
The stopped/held Jpsi/Phi and growing Omega make actual count a
**snapshot**, approximately 8,000, not a synthetic 8,000 exactly.

Read `README.md` plus
`wiretruth1800_distributions.cc`,
`run_wiretruth1800_analysis.sh`, `pbs_wiretruth1800_analysis.sh`.
The Bash/python component is metadata/hash integrity ONLY; numerical
physics/statistics are C++/ROOT.

## Physics definitions — avoid common pitfalls

- Each physical muon crossing is one statistical event. 20 muons per
  chunk, not 20 primary electrons. Nonzero per-wire rows outnumber muons.
- `Q_{muon,k}` is *actual collected Garfield electron endpoints* on
  physical anode `k`. `N_{ion,muon,k}` is exact ionisation callback
  positive-ion births assigned to the **nearest birth anode**; they
  are not ion drift endpoints, and need not equal Q.
- Per-primary `microscopic_wire_seed.csv` and per-muon
  `microscopic_wire_muon.csv` represent different aggregation levels.
  The new ROOT file saves `muons`, `wires`, `seeds`,
  and `seed_wires` TTrees for both.
- `Q2/Qtotal` is the second-largest collected-wire share. This is not
  a second-strip signal; it is not guaranteed to be an adjacent wire.
  `Q2/Q1` has a different denominator and must not be confused.
- Full zero-gain muons remain in multiplicity and sharing-rate
  denominators. Ion births with no collected endpoints remain genuine.
- Positive-wire-only gain log histograms are **conditional Q>0** and must
  be labelled accordingly; they do not estimate absolute detection
  efficiency.
- Station-0 nearest global wire phase is **1.5 mm** for pitch 3 mm;
  stations 1..5 use phase 0 mm. Analyze real Geant4 signed
  theta_u and theta_w and original impact relative to nearest anode.
- The snapshot has *exact* per-wire full-population ion births and
  endpoint totals; sampling weights of 1800V ion-birth reservoir are
  NOT needed for these quantities.
- All plots: **NA60+ Simulation**, 1800V P3 noFW numerical model,
  Ar/CO2 70:30, Penning r=.57, symmetric ±3mm cathodes,
  30µm wires at 3mm pitch. NOT the verified as-built chamber.

## Updated workflow: 60-muon smoke test before ~8k full analysis

The analysis branch additionally supports
`bash test/mwpc-stage-c/wiretruth1800/run_wiretruth1800_analysis.sh pilot`.
It submits only a 1-core/4GB/1h OpenPBS `w1800pil` job: one frozen,
SHA256-verified 20-muon chunk from each of Jpsi, Phi and Omega
(60 real source Geant4 muons total). The independent C++ code validates
per-primary, per-wire and per-muon charge, birth-sector totals and
reported asymmetries before exporting ROOT plots.

Wait for `Exit_status=0`, inspect `asymmetry_denominators.csv`,
`asym_gain_vs_A_hemi_wire_all.png` and
`asym_A_hemi_vs_local_wire_all.png`. If all correct, the new output
prefix with `bash .../run_wiretruth1800_analysis.sh submit` selects all
completed chunks (roughly 8k original muons); do not alter or resume
the production arrays just for statistics.

WIPP November 23 evidence: GCC 13.1.0 / ROOT 6.43.01 compiled the
analyzer successfully (warnings only). `qsub: Access to queue is denied`
was caused by an inaccessible hard-coded `normE` queue, now replaced by
`N`, with `MWPC_1800_PBS_QUEUE` override. New PBS science result still
requires actual successful run; do not infer any asymmetry values from build.

The code has been committed but never executed at WIPP in this chat.
A syntactically valid source cannot establish that it compiles against
local ROOT or that the physics conservation gates pass.
Do not report real gain/asymmetry values until user shares actual output.

## Current deliverables

1. `channel_station_summary.csv`, independent muon denominator with
   Q2_share>5% binomial SE and exact charge/ion totals;
2. `wiretruth1800_8k.root`, original crossing + physical wire rows;
3. `hLogQWire.png`, `hLogIonWire.png`, `hNWire.png`,
   `hSecond.png`, `hNear.png`, `hA.png`;
4. `share_probability_vs_phase.png`,
   `share_probability_vs_theta.png`;
5. `all_per_wire_fraction_collected_electrons.png` and
   `all_per_wire_fraction_ion_births.png`;
6. per-channel per-wire fraction heatmaps and SHA-sealed input manifest.

## Validation needed before physics interpretation

The source C++ has **not yet been compiled on WIPP**. First user runs:
```bash
bash test/mwpc-stage-c/wiretruth1800/run_wiretruth1800_analysis.sh build
```
then submits isolated `N` 4GB/4h PBS analysis:
```bash
bash test/mwpc-stage-c/wiretruth1800/run_wiretruth1800_analysis.sh submit
```
Review actual job `Exit_status=0` and SHA256 output check.
Do not report a measured charge-sharing result from code alone.

The first meaningful question for the actual 8k data:
**Does multi-anode charge sharing increase near half-pitch boundaries
and at large signed incident theta_u, and is that dependence
consistent across the 6 stations and Jpsi/Phi/Omega?**
Check for nearest-wire Q anomalies, non-adjacent Q2, high-ion zero-Q
cases and high-gain tails.

If/when code review is needed, focus on:
- handling of exact frozen crossing IDs, station-aware physical wire
  centre, 9-wire boundaries, and per-seed to per-muon aggregation;
- accidental use of ion birth nearest-wire labels as actual
  collecting anode endpoints;
- ROOT histogram numerator and denominator consistency, conditioning
  on positive collected charge, snapshot truncation/selection biases;
- whether the 3mm nine-wire local patch underestimates high-angle
  or far-boundary tails and makes Q2 overly dependent on field limits;
- quantitative uncertainty of channel/station binomial probabilities.
Do not swap in a toy strip weighting field or reinterpret the ion
hemisphere asymmetry as measured cathode U/V signal.

## New analysis requested: gain versus true ion-birth asymmetry per anode

The same analysis branch now includes
`wiretruth1800_asymmetry.hh`, integrated into the ROOT/C++ executable.
Do not run Garfield again or invent plots before completed WIPP output.

**Exact definitions:** For every native `microscopic_wire_muon.csv` physical
anode row, `A_ion_hemi = (ion_up-ion_down)/(ion_up+ion_down)`,
and `A_ion_TB=(ion_top-ion_bottom)/(ion_top+ion_bottom)`.
Hemisphere `up/down` means birth azimuth above/below the anode x-axis;
`top/bottom` are the restricted 45°–135° and -135°– -45° sectors.
Exact boundary/side counts are excluded from the respective denominators.
When denominator is zero, A is undefined (NaN); never set to zero.
Anode ions are assigned by **ion birth location**, electrons by
**actual Garfield collection endpoint**. These wire populations differ.

**Gain denominators:**
- `G_{muon}=Q_{muon}/N_{Heed,muon}` (effective gain);
- `G_{muon,wire}=Q_{muon,wire}/N_{Heed,muon}` (effective contribution
  of one collecting anode);
- `G_{primary,wire}=Q_{primary,wire}` (collected endpoints produced from
  one initial Heed electron);
- `Q=0` is a physically valid response and must stay in the ROOT trees
  and `asymmetry_denominators.csv`, though it has undefined log10 gain.

New PNGs under output include
`asym_gain_vs_A_hemi_wire_all.png`,
`asym_gain_vs_A_hemi_muon_all.png`,
`asym_gain_vs_A_hemi_primary_wire_all.png`,
`asym_gain_vs_A_TB_wire_all.png`,
`asym_A_hemi_vs_A_TB_all.png`,
`asym_A_hemi_vs_local_wire_all.png`,
`asym_A_hemi_vs_signed_impact_all.png`,
`asym_A_hemi_vs_theta_u_all.png`,
`asym_A_hemi_vs_theta_w_all.png`,
`asym_muon_sharing_vs_A_hemi_all.png`,
`asym_charge_fraction_vs_A_hemi_wire_all.png`,
`asym_ion_births_vs_A_hemi_wire_all.png`,
and `asym_zero_collected_ion_wires_A_hemi_all.png`.
Every channel has equivalent figures, and 18 station/channel subsets
have station-resolved gain-A and A-local-wire plots plus all ROOT histograms.

The baseline reads only files independently locked to the original
1800 V Geant4/Heed/Garfield SHA-verified campaign.
`asymmetry_denominators.csv` specifies the actual per-wire/per-muon/
per-primary-wire population of each graph.

**Interpretation caution:** h2 densities have one entry per ion-bearing
physical wire (or per primary wire), not per independent muon, except
specifically labeled muon plots. A large-gain primary can contribute
several wires. TProfile ROOT bin errors are not event-cluster bootstrap
uncertainties. The comparison to 2100V needs common crossing/geometry
and real Garfield, not a plot from weighted hypothetical strip ions.

Before model decisions, request WIPP job `Exit_status=0`,
`asymmetry_denominators.csv`, and 4–6 key PNGs.
Scientific hypothesis: whether per-wire collected gain and Q2 sharing
retain a conditional dependence on ion-birth asymmetry once anode
offset, station, incident theta_u/theta_w, and channel are accounted for.

## Future detector strategy

First obtain a coherent current-geometry baseline (even though numerical
P3 reference differs from actual as-built hardware), then parameterize
strip width/pitch, stereo angle, wire pitch, cathode gap and HV in a
configuration-lock scheme. Changing electrostatic geometry or field
typically requires new weighting/drift-field solves, and can also require
new avalanche physics; do not reuse one geometry's response blindly
for another's strip pitch/width.

When the user returns data, ChatGPT compares these real 1800V
distributions to validated 2100V data only with aligned input/geometry
and decides which missing phase-space cells or uncertainty justify
new microscopic jobs.

## Scheduler fix (confirmed from original production)

On 2026-11-23, WIPP's `N` queue was shown to be a routing queue with
`resources_min.io=0.1` and destinations longE/shortE/normE/gpuE.
Original successful 1800V worker `pbs_wiretruth30k_chunk_v1.sh`
explicitly requested `#PBS -l io=0.1` in addition to `#PBS -q N`.
The first analysis `normE` submission was denied access; the second
`N` submission omitted IO and was rejected by queue resource limits.
The new launcher and analysis worker now request `io=0.1` explicitly.
No analysis result has yet been reported; wait for a real PBS job ID and
completed SHA-verified pilot before interpreting gain-versus-asymmetry plots.


## Verified real 60-muon outcome; next performance gate

See [REAL_PILOT60_RESULTS_20261123.md](REAL_PILOT60_RESULTS_20261123.md)
for source-anchored counts, actual visual findings and appropriate
limits on gain-versus-ion-hemisphere-A statements.

The user's genuine WIPP 60-muon job `5588743.pbs` exited 0,
with 35m36s walltime but only 7s CPU. Source manifest and exact
wire/primary electron and ion conservation passed. At muon level 15/60
had exactly two collecting anodes; 10/60 had Q2/Qtotal>5%;
75 physically ion-bearing anode rows, 2661 per-Heed-primary/anode
ion-bearing rows and only one ion-bearing primary-wire with zero Q.
Total effective gain Q/N_Heed=3225.04 at 1800 V. No physical cathode
current, reconstruction efficiency, or calibrated sign of U/V waveform
follows from this result. The supplied result tar omitted the TTree ROOT
file and 18 other named products; all 112 included files passed SHA.

Because CPU/wall ratio suggests I/O wait, the GitHub analysis launcher
now supports `pilot180` (3 original 20-muon chunks per channel,
180 actual muons, resource 1 core/4 GB/4h/io=.1) and optional node
TMPDIR staging of ROOT and PNG output. Its stdout stage timestamp
markers isolate bottlenecks while preserving the full metadata/sha
verifications.

Do not advise a single full ~8k 4h PBS job before the 180 resource
gate. If physical correlation is requested now, insist on original
per-event ROOT tree; visual scatter cannot give slope, Pearson r,
bootstrap confidence interval, or channel/station adjustment.


## 2026-11-23 decision: run ALL completed original 1800V ~8k muons

The user specifically requested **all ~8000 now**, without an intermediate
180-muon performance pilot. ChatGPT implemented a new full-sample
**read-only ROOT-only** PBS archive campaign. It selects every completed
SHA-receipt-backed original 20-muon WIPP chunk at preparation time
(actual number recorded, not assumed exactly 8k), then freezes its
channel + original crossing index set.

Critical GitHub files:
- FULL_8K_CAMPAIGN_20261123.md (complete WIPP launch/merge instructions)
- full1800.sh (build/prepare/submit/status/resume/merge-submit)
- full_wiretruth1800_campaign.py (receipt+source SHA provenance, PBS-local
  copies, sealed independent shard outputs, strict full merge gate)
- pbs_wiretruth1800_full_shard.sh and pbs_wiretruth1800_full_merge.sh
- wiretruth1800_distributions.cc (new --no-plots flag for shards; **no
  physics changes** to gain/A definitions)
- wiretruth1800_merge_full.cc (real C++/ROOT TFileMerger, exact
  crossing_id uniqueness, summed integer counters, independently
  recomputed full-population Q2 fraction and uncertainty, final
  channel and station plot rendering)

Source physics is historical 1800V P3 noFW 3mm pitch symmetric +/-3mm
analytic nine-wire model, not the as-built chamber and NOT a 2100V
uncapped simulation. All collected charges are genuine Garfield electron
wire endpoints and ion A values are true callback positive-ion *birth*
sector counts. Keep the distinction absolute.

WIPP pilot job 5588743.pbs (60 genuine muons, 35m36s wall, 7s CPU,
exit 0) demonstrated I/O domination, motivating default 6 original
20-muon chunks per PBS array member, at most 2 simultaneous members,
N routing queue with io=0.1, 1 CPU, 4GB, 8h per member. New stage
copies all five input files to node TMPDIR while simultaneously verifying
their frozen SHA, runs identical ROOT physics without redundant per-shard
PNGs, seals each separate shard on shared storage. The final PBS job
first validates ALL shards, merges TTrees and TH1/TProfile distributions,
checks original crossing IDs unique and exact primary/wire population,
and generates final PNGs once.

No full simulation/analysis job is claimed launched or successful at
the moment of this handoff. User runs the commands and returns PBS IDs
and actual outputs. Before reporting a physical correlation from these
anodes, check independent-muon clustering, sample conditioning by
station/source/impact/angle and original ion count.

The full 1800V campaign MUST NOT modify held/original Garfield production
arrays, old verified source chunks, 2100V ion/strip solver, or source
geometry. If a PBS task errors, correct the source on a *new frozen
analysis campaign* rather than silently mixing different ROOT executable
SHA versions within one selection.
