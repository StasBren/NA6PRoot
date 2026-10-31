# NA60+ MWPC — protected spatial LUT with shared finite-training bootstrap

**31 October 2026 — completed, source-conditioned physics diagnostics.**

## Goal and inputs

Evaluate protected spatial-neighbor pooling of original A20 and targeted A100 empirical Garfield++ avalanche response libraries, including the *finite number of independent avalanche outcomes used to train each cell*.

Physical provenance (all hash-verified): 84,000 independent old20 full microscopic Garfield++ outcomes, 16,000 new80 targeted outcomes (50 selected cells in each of 4 geometries), one full-microscopic reference outcome per real parent×geometry, and the frozen 5,881 TrackHeed electrons of 120 distinct Geant4 parent muons. Total 480 correlated physical reference cases, four geometries per parent. **No new Garfield++, Geant4, TrackHeed, or Magboltz simulation.**

## Method

Within **each global bootstrap repetition**, resample old20 outcomes with replacement **once per detector cell**. Reuse this common resampled training cell for *every* muon and every primary electron that visits the cell. Independently bootstrap additional new80 outcomes in A100-selected cells; respect their 20% old20/80% new80 mixture, sharing the old20 component with A20 to reduce comparison noise (common random numbers). Conditional on the globally bootstrapped training libraries, sample one whole 5-component correlated Garfield ion-census vector per primary. Evaluate A20, A100, guarded neighbor-pooled A20 and guarded neighbor-pooled A100 on exactly the same primary coordinates. Interval is central 5th–95th percentiles of this marginal shared-library/bootstrap predictive distribution. Four detector geometry outcomes are treated as one Geant4-parent cluster when bootstrapping confidence intervals for changes in coverage.

Guard/kernel fixed before evaluation: nearest 5×5 original A20 donor cells, Gaussian distance weighting exp[-(du²+dv²)/2], 25% local / 75% neighbors. Pooling prohibited across v=0 ion hemisphere; wire-adjacent cells below 0.35 mm (center-to-center), and near-axis opposite-u-sign donor cells are excluded. No periodic wrapping. No ion spatial clouds, ion drift, induced U/V strip signal, avalanche time or B≠0.

## Result — full 5,000 global bootstrap repetitions (random seed 20261031)

| Model | Nominal 90% interval covered (480) | Empirical coverage | Mean absolute predictive-mean error (ions) | Mean 90% width (ions) |
|---|---:|---:|---:|---:|
| Original A20 | 406/480 | 84.58% | 22,391 | 87,085 |
| Targeted A100 | 399/480 | 83.13% | 20,858 | 79,477 |
| **Guarded pooled A20** | **417/480** | **86.88%** | **18,750** | 78,888 |
| Guarded pooled A100 | 414/480 | 86.25% | 18,843 | 77,613 |

Parent-clustered bootstrap 95% confidence intervals on *paired coverage differences*: guarded A20 minus A20 **+2.29 percentage points**, CI **[0,+4.58] pp**; guarded A100 minus A100 **+3.13 pp**, CI **[+0.83,+5.42] pp**. A20 coverage improvement is *borderline* (interval touches 0); must not be overstated. Both pooled approaches still fall below the nominal 90%.

Ion hemisphere asymmetry empirical 90%-interval coverage: A20 90.21%, A100 89.79%, guarded A20 90.63%, guarded A100 90.21%. Top/bottom asymmetry: 88.33%, 88.33%, 89.38%, 88.96%, respectively. These are **ion-sector counts, not induced detector signals**.

Per geometry (covered / 120):

| Configuration | A20 | A100 | Guarded A20 | Guarded A100 |
|---|---:|---:|---:|---:|
| p3_noFW | 98 | 95 | 101 | 100 |
| p3_FW75 | 104 | 102 | 106 | 104 |
| p4_noFW | 102 | 99 | 105 | 105 |
| p4_FW75 | 102 | 103 | 105 | 105 |

Independent pseudo-random sensitivity check: another 3,000 global bootstrap repetitions, seed 7, yielded covered counts A20 409, A100 401, guarded A20 421 and guarded A100 418; same qualitative result but Monte Carlo quantile noise affects exact integers and coverage-delta CIs.

**Cross-check:** prior independent A20 vs A100 *non-pooled hierarchical bootstrap* produced 407/480 and 399/480 respectively (published previous work). New 406/480 and 399/480 are consistent given different seeds/repetition counts.

## Interpretation, limits and next experiment

The spatially pooled model reduces predictive-mean error of ion births by ~16% relative to A20 after including finite-library uncertainty, while modestly improving nominal coverage. With global training bootstrap, the old A20 coverage rises from 76.7% in conditional-only replay to 84.6%; therefore *the incremental pooling contribution to interval coverage is only +2.3 pp*, not the +9 pp observed in the conditional-only comparison. Additional A100 statistics **do not clearly improve pooled performance**, which weakens the rationale for immediately filling all 4,200 cells to K=100.

Cautions: (1) these 120 Heed source coordinates contributed to A100 cell selection; (2) spatial pooling hyperparameters were chosen after diagnostics on high-variance/exposure-selected A100 cells; (3) global bootstrap is an approximation to library-estimation uncertainty, **not** an exact Bayesian posterior or independent calibration; (4) only one random microscopic truth per parent×geometry; (5) guarded 0.35mm criterion is defined at cell-center resolution and is not a field-line/metal-surface proof; (6) sector counts do not establish the distribution of actual U/V cathode strip induced charge.

**Decision:** hold new expensive microscopic Garfield campaign while preparing genuine *new 100–120 Geant4 parent muons* with new frozen TrackHeed source and full microscopic Garfield truth, analyzed only after model/hyperparameters are frozen. Do not use the new parents to choose A100 cells or retune pooling before the blinded comparison. Compare empirical 90%-PI coverage and predictive-mean error by geometry, track distance from anode/field-wire; then determine whether a 48k targeted or 320k uniform full Garfield avalanche expansion is justified. Costs previously measured at ~0.35–0.46 s / original microscopic avalanche, not a guarantee for new runs.

## Code and reproducibility

`test/mwpc-stage-c/pi_spatial_lut_muon120_hier_v1.py` (separate from original old physics).

In user's original WSL repository after pulling latest branch:

```bash
cd ~/na6p/src/NA6PRoot-mwpc
git pull --ff-only
python3 test/mwpc-stage-c/pi_spatial_lut_muon120_hier_v1.py selftest
python3 -u test/mwpc-stage-c/pi_spatial_lut_muon120_hier_v1.py run --native --reps 5000
```

Native run refuses overwrites, uses frozen source verifier, and writes only under `test_runs/mwpc_stage_c/pi_A100_targeted_50cells_per_geometry_v1/diagnose_guarded_spatial_pool_global_shared_training_bootstrap_v1/` (unless given another new `--output` path).

For SHA-checked published input archive instead of WSL original chunks:

```bash
python3 pi_spatial_lut_muon120_hier_v1.py run --published-dir /mnt/data --output /mnt/data/NEW_HIER_RUN --reps 5000
```

Both use original old20/new80 micro ion births and five ion-sector counters, not synthetic physical input. The synthetic selftest is distinct and is not evidence of physical validity. The code was exercised in archive mode with full 5,000 physical-source-conditioned global bootstrap replicates and with different seed 7 × 3,000 replicates.