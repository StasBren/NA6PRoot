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
then submits isolated `normE` 4GB/4h PBS analysis:
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
