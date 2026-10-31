# NA60+ MWPC: existing fast avalanche models and independent tests (2026-10-31)

## Research question

Can a cheap model replace event-by-event Garfield++ `AvalancheMicroscopic` while preserving the physical distribution of ion births and eventually the U/V cathode-strip signal?

**Scope:** Four existing zero-magnetic-field wire-chamber geometries: `p3_noFW`, `p3_FW75`, `p4_noFW`, `p4_FW75`. Gas Ar/CO2 70/30. All new tests below are read-only analyses of **original frozen A20 (84,000) and A100 added trials (16,000)**. No additional Garfield/TrackHeed/Geant4 physics was executed.

## Credible, available precedents (checked online)

1. **Garfield++ `DriftLineRKF` + `SetGainFluctuationsPolya(theta, mean)`**: standard supplied Garfield++ functionality. The mean multiplication follows Townsend and attachment coefficients along drift lines; Pólya samples gas-gain fluctuations. It can simulate induced signals and an ion tail if appropriate mobility/weighting-field information is supplied. **This is the strongest reusable on-machine alternative to a hand-built LUT**, but is not microscopic event-by-event physics and must be validated against full Garfield++ and, ultimately, detector measurements. Garfield++ User Guide: https://garfieldpp.web.cern.ch/doxygen/UserGuide.pdf ; detailed author clarification: https://root-forum.cern.ch/t/help-on-garfield-avalanchemicroscopic-and-avalanchemc/53678 . A concrete user-side code example: https://root-forum.cern.ch/t/signal-calculation-of-low-energy-particles-and-fission-fragments-in-a-micromegas-geometry/45235
2. **GTS / a fast and parametric digitization for triple-GEM detectors**: published simulation factorizes primary ionization, diffusion, gas gain (fit Pólya), induction, electronics, and validates against measurements. The development pattern directly matches our goals, though GEM detector physics differs and a reusable production code package for **our MWPC** was **not verified**. R. Farinelli et al., https://arxiv.org/abs/1904.06142 ; slides https://indico.cern.ch/event/773049/contributions/3474785/attachments/1931837/3211586/20191024_GTS_v3.pdf . The paper explicitly reports experimental tuning of gain/diffusion, so parameterisation alone is not guaranteed to reproduce measured MWPC response.
3. **Garfield++ accelerated microscopic calculation**: CUDA acceleration for `AvalancheMicroscopic`, described in JINST 20 (2025) P12019, https://doi.org/10.1088/1748-0221/20/12/p12019 . Promising for future high-statistical campaigns, contingent on installed version, compilation, hardware, and validation; a GPU cannot by itself fix inappropriate spatial response binning.
4. **CERN DRD1 Garfield overview:** analytic electric field for 2D wire chambers, HEED, Magboltz, charge induction via weighting fields, electronics convolution are built-in capabilities; no need to reimplement all that in isolation. https://drd1.web.cern.ch/garfield/ .
5. **Other methods:** `AvalancheMC` is macroscopic stochastic drift/avalanche, and `AvalancheGrid` is tailored to parallel-plate RPCs. These are not drop-in replacements for our wire-chamber geometry without physics review. `AvalancheGrid` DRD1 presentation: https://indico.cern.ch/event/1489983/contributions/6363875/attachments/3021531/5331074/Slides.pdf .

## Two independent hold-out diagnostics performed on actual 16k new trials

The **original A20 20 avalanches per cell were used to rank the 200 targets**. To avoid reusing them to test smooth spatial dependence, spatial cross-validation uses *only* the independently generated **new80** from each target cell. The Pólya experiment uses **original20** for fitting and **new80** for external trial-level testing. Exact source CSV SHA256 were checked against the frozen published source manifests.

### Test A: does smooth subcell position improve gain prediction?

Within each selected 0.1×0.2 mm cell, predict either `log(1 + ion_births)` or the probability of zero ion births using: constant cell mean, linear in exact start `(u,v)`, or quadratic in `(u,v)`; ridge penalties suppress noise. Fourfold cross-validation on only the new80 per cell, 16,000 total held-out response observations across folds. Difference is reduction in **out-of-sample squared loss** versus constant-cell model (positive = improvement).

| Geometry | Linear gain-log MSE reduction | Quadratic gain-log MSE reduction | Linear zero-response Brier improvement |
|---|---:|---:|---:|
| p3_noFW | **−2.40%** | −2.24% | −2.50% |
| p3_FW75 | −1.63% | −1.88% | −1.82% |
| p4_noFW | −1.40% | −1.35% | −1.38% |
| p4_FW75 | −2.43% | −2.15% | −2.38% |
| Overall 200 cells | **−1.95%** | **−1.90%** | **−2.01%** |

A 4,000-draw bootstrap over **200 independent selected cells** yields linear log-gain loss-reduction 95% percentile interval **[−2.45%, −1.42%]**. Therefore within these selected cells neither naïve smooth model beats constant-on-cell on these specific held-out loss metrics. This does **NOT** disprove finer, nonlinear physical structure near wires, gradients **between** different cells, source position distribution shift between D500 and TrackHeed, or conditional ion spatial/timing distribution dependence. The 16k new independent avalanches also lack direct simulated U/V signals.

### Test B: does a physically motivated Pólya gain law work on unseen avalanches?

For each selected cell fit a **zero-mass plus continuous Gamma/Pólya conditional on positive gain** using the 20 original Garfield avalanches. The zero mass has a Jeffreys-half-count smoothing prior, and the conditional positive-Gamma is fit by maximum likelihood. Contrast 5th–95th percentile intervals with raw empirical A20 5th–95th percentiles using **the independent new80 in each selected cell**, 16,000 hold-out observations.

| Geometry | 90%-nominal empirical A20 coverage | Pólya/zero mixture coverage |
|---|---:|---:|
| p3_noFW | 91.43% | 94.90% |
| p3_FW75 | 90.28% | 94.00% |
| p4_noFW | 91.50% | 94.05% |
| p4_FW75 | 90.65% | 95.33% |
| **Overall** | **90.96%** | **94.57%** |

Bootstrap over 200 selected cells for the difference Pólya minus empirical: **[+2.73,+4.54] percentage points** (95% percentile). The Gamma **shape** median is around 1.34–1.50 by geometry (equivalent Pólya θ approximately 0.34–0.50 in the common shape=θ+1 convention). This is a descriptive fit, not a universal gas property.

**Interpretation caution:** Roughly 15–18% of original/heal-out single-electron trials yield zero ions. Thus the 5th-percentile bound in a zero-atom mixture is often exactly zero, while the 95th bound cuts off ~5% on the high side. **Even a perfectly specified two-sided nominal 90% interval may have ~95% actual coverage with an atom at zero.** Hence 94.57% is not by itself evidence Pólya is a better-calibrated 90% predictor than empirical A20, or a fit to the entire gain distribution. Evaluate separate upper-tail calibration and proper scoring rules before choosing a distribution. The empirical A20 has ~6.3% upper-tail misses and Pólya about 4.8% on the new80 holdout.

### Key scientific inference

The independent new80 data show that the original A20 sampling gives approximately nominal single-electron predictive coverage (90.96% on 200 *selected* cells), while original full-muon validation gives ~77–79% gain coverage for a nominal 90% interval. **The main mismatch emerges at the muon response aggregation/conditionalization layer**, not unambiguously from an incorrect single-electron mean or a readily fixed smooth within-bin gradient. Shared finite training-cell uncertainty, correlated multiple electrons per cell, large ion-gain tails, exact metal-start exclusion, and *different distributions of starts inside each spatial cell for D500 and TrackHeed* remain specific hypotheses. The tests do not show each is independently causal.

## Recommended next executable detector-level test

1. **Do not discard existing LUT or add more avalanche simulation yet.** Preserve all chunk-level outputs, including ion-birth spatial/time reservoirs.
2. **Build a separate `DriftLineRKF + Polya` B=0 prototype** using the exact same `ComponentAnalyticField` used by existing native MWPC Garfield C++ driver (wire positions, voltages, planes, gas), and the frozen 120 parent-muon Heed starts. Do not modify full physical reference or simulate a new Geant4 source yet.
3. **Require explicit zero-response probability** from electron loss to attachment, anode/field-wire losses, inactive starts, plus conditional Pólya gain *only if the electron reaches an amplification region*. Do not use a Pólya distribution alone to force every electron to avalanche. Be precise whether comparing **ion births** versus Garfield's final electron count; they are not exactly the same observable.
4. Validate at per-primary and per-muon scales: p_zero, means, variance, 5th/50th/95th, tails, coverage, source-angle and station strata, geometry, and runtime. Do not tune on the 120 full Garfield reference and claim unseen validation. Once functioning, evaluate on newly independent parent muons and/or experimental data.
5. Later compute ions' trajectories, weighting fields of both U/V strip planes, Ramo current and electronics response. The present input A20/A100 trials contain ion-birth counts and coarse sectors, **not the full induced U/V signal**.

## Reproducibility

Two read-only Python scripts were created, executed on original SHA-verified data, and committed in `StasBren/NA6PRoot`, `feature/mwpc-stage-c`:

- `test/mwpc-stage-c/pi_A100_subcell_spatial_diagnostic.py`
- `test/mwpc-stage-c/pi_A100_polya_holdout_20v80.py`

Both allow `--old-dir`, `--new-dir`, `--out-dir`; output to new directories only and refuse to overwrite existing outputs. The local downloaded source archive has separate A20/A100 physics rows in one `/mnt/data` folder, whereas WSL's original directories are separate. As usual, the same 120 muon *source coordinates* were used to select targets, so neither test is fully representative of an independent unseen muon flux.