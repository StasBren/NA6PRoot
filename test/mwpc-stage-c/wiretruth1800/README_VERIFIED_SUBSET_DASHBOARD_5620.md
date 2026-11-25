# 1800 V MWPC wiretruth: verified 5620-muon subset / dashboard input

This branch creates a **new, independent ROOT merge** from exactly 48 previously
SHA-sealed C++/ROOT shards of the archived 1800 V campaign. It does not
regenerate Garfield++ avalanches or change original physics archives.

Original frozen source: `~/NA6PRoot-wiretruth1800-analysis/test_runs/wiretruth1800_FULL_20261123T145716Z`.
Excluded missing shard IDs: 0, 5, 7, 16, 17, 27, 32.
Expected retained events: **5620 = 2060 Jpsi + 2520 Phi + 1040 Omega**.
The exclusion is not random; avoid presenting it as the complete 6460-muon dataset.

## Job stages
1. Create a *separate* worktree from `feature/p3-1800v-5620-subset-dashboard-20261125`.
2. Export `MWPC_1800_SUBSET="$PWD/test_runs/wiretruth1800_SUBSET_5620_<unique_timestamp>"`.
3. Run `bash test/mwpc-stage-c/wiretruth1800/subset_1800v.sh prepare`.
4. Review prepared subset count and provenance.
5. Submit the one-time **ROOT-only** PBS job with `bash .../subset_1800v.sh submit`; monitor via `bash .../subset_1800v.sh status` and the separate PBS log.
6. Inspect `<subset>/merged/` only after PASS and SHA checksum verification.

The worker stages each required sealed ROOT/CSV input to PBS local scratch while
verifying SHA-256; invokes the **unchanged original frozen C++ merger**, checks
distinct crossing IDs and original-primary counts, plots, writes new provenance
and SHA, then copies output to the new subset area. The original 1800 V
`merged` directory is never written.

## Intended dashboard sections (C++/ROOT, ATLAS-style)
- **Population gain**: total ion births and collected electron endpoint count
  per muon; total collected electrons per Heed primary using ratio of sums;
  distinguish zero-gain population. From
  `channel_station_summary.csv`, compute ion gain as
  `total_ion_births/total_heed_primaries` and electronic gain as
  `total_collected_e/total_heed_primaries`. The `muons`, `seeds`,
  `wires`, `seed_wires` TTrees support additional distributions.
- **Charge sharing**: active collecting-wire multiplicity, `Q2/Qtotal`,
  second-wire fraction >5%, signed impact position, muon incidence angle,
  nearest-anode assignment and wire charge fractions.
- **Up/down ion asymmetry only**:
  `A_UD=(N_up-N_down)/(N_up+N_down)`, defined only for nonzero
  `N_up+N_down`; display its distribution and dependence on electronic
  gain, local wire, signed impact, signed angles and charge sharing.
  The old merger contains legacy top/bottom plots, but **do not use them in
  the dashboard**.
- **Longitudinal ion development**:
  `sigma_w_birth_mm` and `mean_w_birth_global_mm` are stored for
  primary-wire and aggregate-wire ROOT rows; their dependence on gain can
  be plotted in a separate ROOT postprocessing macro. Distinguish
  per-avalanche spread from multiple spatially separated primary electrons.
- **Angular ion development**: angular *asymmetry* is recoverable from
  exact up/down full-census counts, but circular angular *spread* cannot
  be reconstructed exactly from only those counts. Use original
  per-primary ion-birth reservoir samples with sampling weights or extend
  the production output to retain circular moments. Never label a
  hemispheric asymmetry as an angular standard deviation.
- **Later response dashboard**: ion-source packets and per-wire positions
  may feed weighting-field / Shockley--Ramo cathode U/V signal models,
  but stored ROOT charge distributions are **not** induced strip currents.

## Scope disclaimer
Existing C++ merger produces population and up/down ROOT plots plus some
legacy top/bottom PNG files; dedicated gain-vs-longitudinal-spread and
circular-angular-spread graphs are **not yet implemented**. The first step
is to obtain one SHA-verified merged ROOT data product, then add these new
derived visualizations without rerunning Garfield++.
