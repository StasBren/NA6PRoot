# Actual 1800 V per-anode gain/asymmetry pilot — 2026-11-23

**Status:** Measured ROOT/C++ 60-muon WIPP pilot. No additional Garfield++
simulation. Exact numbers below are from the user's **uploaded** PBS analysis
`channel_station_summary.csv`, `asymmetry_denominators.csv`,
`analysis_readme.txt`, `provenance.txt` and physically generated
`asym_*.png` in
`wiretruth1800_asym_pilot_results.tar.gz`. These source files are not
committed to GitHub; only the findings and provenance are recorded here.
Do not mistake the text for an executable or data file.

## Physical and computational provenance

- Physics campaign: 1800 V / P3 noFW numerical nine-wire, 3 mm anode pitch,
  30 μm wire diameter, ±3 mm symmetric cathodes, Ar/CO₂ 70:30, Penning r=.57.
- Original code SHA `599e29522edb6fcec9dfcca2904308cb3186edce`.
- Original frozen binary SHA256:
  `2cf92cb27d3a180cd26d3ab5417f57293473864f6d685ee43bacf931712aa8f1`.
- Pilot analysis source commit:
  `41d8dec957713d5f0a08edfdd31a8cc6fb48ca01`.
- ROOT 6.43.01, GCC 13.1.0, PBS `N` queue, `io=0.1`.
- Three SHA-verified original 20-muon chunks (one each Jpsi, Phi, Omega).
- PBS job `5588743.pbs`: exit status 0, wall 00:35:36,
  CPU 00:00:07, peak reported memory 547392 kb.
- Runner's terminal `analysis_readme.txt` records exact
  per-physical-wire source vs Heed primary vs Garfield collected
  electron and positive-ion birth conservation as PASS.
- The user's transfer tar is a **selected bundle**:
  112 SHA256-listed files are present and checksum-correct;
  19 entries in the original output SHA256SUMS were intentionally not
  included (the source ROOT TTree, input manifest and general plots).
  Thus `sha256sum -c SHA256SUMS` cannot pass **within this partial tar**,
  even though its 112 supplied files have no mismatches.
  Do not claim any quantitative ROOT-tree correlation was calculated here.

## Measured physical totals

| | All | Jpsi | Phi | Omega |
|---|---:|---:|---:|---:|
| Independent Geant4 muons | 60 | 20 | 20 | 20 |
| Original Heed primaries | 3148 | 1169 | 977 | 1002 |
| Actual collected e endpoints | 10,152,431 | 3,752,167 | 3,167,121 | 3,233,143 |
| Exact born positive ions | 10,149,698 | 3,751,145 | 3,166,300 | 3,232,253 |
| Total Q / total Heed primaries | 3225.04 | 3209.72 | 3241.68 | 3226.69 |
| Single collecting anode | 45 | 15 | 14 | 16 |
| Two collecting anodes | 15 | 5 | 6 | 4 |
| Three-or-more collecting anodes | 0 | 0 | 0 | 0 |
| Q2/Qtotal > 5% | 10 | 4 | 2 | 4 |
| Completely zero-charge muons | 0 | 0 | 0 | 0 |

Muon-level probability of collecting on exactly 2 anodes is 15/60
=25%; fraction with Q2>5% is 10/60 =16.67%. Quoted independent-muon
binomial standard error for the latter is 0.0481 (4.81 percentage points),
without systematics or channel/angle/impact conditioning.
This **does not** establish comparable true physical strip multiplicity.

`Q_{collected} - N_{born} = 2733`; individual collected e endpoints
and ion births remain distinct quantities. No signed U/V strip signal
can be inferred from this approximate near equality.

## Asymmetry populations and interpretation

- Full-muon rows with defined A_hemi: 60; gain+defined A_hemi: 60.
- Ion-bearing physical muon×anode rows: 75; both A_hemi and A_TB defined
  for all 75; all 75 have positive collected Q.
- Individual Heed-primary×anode ion-bearing rows: 2661;
  A_hemi defined for all; 2660 have positive collected Q; **1** has
  born ions but zero collected electrons. It is kept in ROOT and
  denominator QA but excluded from log10 gain plots.
- Thus the muon×wire gain vs A scatter has just 75 contributing
  physical-wire observations, correlated within the 60 muons.
  Per-primary×wire scatter has 2660 positive-gain observations,
  also correlated by muon.
- Visual inspection of the **actual** uploaded plot
  `asym_gain_vs_A_hemi_wire_all.png` shows the majority of
  positive-wire gain contributions near log10(G_wire) ≈3.5, with a
  handful of much smaller secondary-anode contributions.
  **No compelling monotonic gain–A trend** can be established with
  only 75 populated physical-wire entries. No quantitative slope,
  p-value or correlation coefficient has been calculated, because
  the pilot ROOT tree was not in the supplied tar.
- The per-Heed-primary×anode plot shows much broader gain variation
  and many A values near ±1; these are real native exact angular
  sector counts, not K256 sampling noise. However, finite ion birth
  counts, station, geometry and per-muon effects must be conditioned
  before concluding a physical relationship.
- A_hemi and A_TB are strongly related by angular-sector definitions;
  a correlation between them is not by itself a new physics result.
- The A versus local anode plot is concentrated around zero offset
  (most muons only collect on one wire) and near ±1 for some second
  wires. Near-zero occupancy away from the first neighboring wire.
- Geometry-angle/impact panels are sparse at 60 muons; station-level
  bins frequently have only 1–7 independent muons. Do NOT claim
  a sign-dependent asymmetry effect yet.
- ROOT TProfile errors are not clustered-by-muon uncertainties.

## Critical computational implication

Despite only 60 muons, PBS used **35m36s wall time and 7s CPU**.
This indicates dominant storage/I/O or waiting time, NOT
a 35-minute numerical calculation. A naive 8,000-muon single
`4h` PBS job would be unsafe. Running the full original sample
unbatched could require tens of hours at the pilot's apparent
walltime/chunk rate; only a stage-timed pilot can distinguish
source hashing, C++ loading, PNG writes and final SHA seals.

The branch's current runner was updated **after** this real pilot
to timestamp SHA source preparation / ROOT physics / output persistence,
write ROOT/PNG products to PBS node-local `TMPDIR` before immutable
copyback, and offer an intermediate `pilot180` mode (3 verified
original 20-muon chunks per source channel = 180 muons).
The new runner changes are not yet exercised on WIPP.

## Decision and next gate

1. **No more Garfield++ simulations for this analysis.**
2. Run the new 180-muon intermediate ROOT-only PBS pilot
   to measure whether the node-local output staging improves elapsed
   time and whether channel/station asymmetry coverage increases.
   The run must be SHA-sealed and exit 0.
3. Inspect the per-stage timestamps and `asymmetry_denominators.csv`,
   compare the population vs the verified 60-muon pilot and
   validate no change in physical definitions.
4. **Do not launch all ~8k in a single 4-hour PBS job** until
   storage performance and batch size are measured. Plan chunk-index
   partitioned jobs + exact ROOT merge with event identity and
   duplicate protection if full-sample walltime is too long.
5. For a quantitative A–gain correlation fit with independent-muon
   cluster resampling, the user must transfer `wiretruth1800_8k.root`
   from the pilot or future full result; the current tar omitted it.
   Condition correlations on ion-birth statistics, incident u/angles,
   station, local collecting wire, and channel. Avoid data-driven
   cherry-picking from 60 muons.

## Claude handoff

The correct question for Claude is whether conditional correlation of
true per-wire A_hemi/A_TB with G_primary,wire and G_muon,wire (electron
endpoint counts) remains after holding primary ion count/geometry
constant. Do NOT substitute wire ion births for collected electron
charge, average over all wire rows as if statistically independent, or
mistake a 1800V original Garfield cap of 200000 (zero capped seeds
in the verified data) for the **uncapped 2100V** campaign.
