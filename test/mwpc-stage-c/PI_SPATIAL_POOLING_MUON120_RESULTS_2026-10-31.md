# NA60+ MWPC — wire-guarded neighbor-pooled avalanche LUT evaluated on real Geant4 muons

**Date:** 2026-10-31. **Status:** Archived-data analysis completed and SHA-sealed; source inputs were not modified. Native WSL replay not yet independently confirmed. No new Geant4, TrackHeed, Garfield++ or Magboltz physics was performed.

## Experimental question

Does spatial neighbor-pooling of physically simulated Garfield++ avalanche outcomes improve the **summed response of an entire real muon**, relative to the old per-cell A20 and targeted A100 libraries? More importantly, does that improvement survive a conservative guard around electrodes, anode/field-wire boundaries and the `v=0` ion-hemisphere boundary?

### Inputs and provenance

- 120 genuine independent Geant4 parent muons (four correlated detector geometries per parent, **480** muon×geometry comparisons).
- Exactly **5,881** saved original TrackHeed primary electrons (reused for each geometry: **23,524** primary×geometry instances).
- 4,200 spatial cells across all geometries, each with **20** full microscopic Garfield++ avalanche outcomes: 84,000 original trials.
- 200 specially targeted cells (50/geometry) have **80 additional independent Garfield++** avalanches, making 100 each: another 16,000 trials. A100 target selection used the **same** 120 muons' Heed positions and old A20 variances, not the full Garfield validation responses.
- Frozen original one-realization-per-muon full microscopic Garfield++ ion and hemisphere-sector counts. All source CSV files checked against their published SHA-256 digests. Old20/noFW and FW75 have different analytic fields; never mix geometries.

### Candidate lookup distributions (six models)

- **A20**: independent uniform draws from the local original twenty full Garfield samples, with real 5-component ion census `(births, up, down, top, bottom)` sampled jointly as one physical avalanche.
- **A100**: draws from all hundred samples in the 50 A100-selected cells per geometry, otherwise the original 20.
- **Guarded pooled A20/A100**: independently for each Heed electron, draw from its local A20/A100 cell with probability **25%** or from the **old-A20** distributions of nearby 5×5 cells with probability **75%**; neighbor weights `exp(-(du²+dv²)/2)` in spatial bin-index units. The guard never mixes positive and negative local `v` across 0, avoids cross-wire `u` sign near the central anode and forbids smoothing if the bin center lies within **0.35 mm** of the anode or any field wire. It does not wrap periodic cell boundaries. This is a simple conservative geometrical heuristic (not an exact field-line separatrix), whose guard distance was not physically calibrated. Only full Garfield outcomes from the **same geometry** may enter.
- **Unguarded pooled A20/A100**: same smoothing probabilities and radius but without electrode/`v=0` exclusions; retained only as an unsafe diagnostic comparison, not a proposed implementation.

The 5-dimensional birth-sector tuple is sampled **jointly** to retain within-avalanche angular correlations, rather than sampling total ion count and asymmetry separately. A20 and A100 are coupled using **common random numbers**: unchanged cells generate **bitwise-identical** A20/A100 realizations, including for 158 muon×geometry cases that have no primary electron in a selected A100 cell. Every model uses independent avalanche draws for different primary electrons, conditional on the frozen empirical response pools.

The radius and mixing fraction were previously explored using **the same targeted A100 selected cells** for independent single-avalanche holdout. They must be treated as *development-selected*, not independent confirmation on the 120 muons.

### Results: 2,000 Monte Carlo response replicas per muon×geometry

| Model | Covered among 480 nominal 90% intervals | Coverage | Mean absolute error of predicted muon ion yield |
|---|---:|---:|---:|
| A20 original | 368/480 | **76.67%** | 22,463 |
| A100 targeted | 378/480 | **78.75%** | 20,869 |
| **Guarded A20 pooling** | **411/480** | **85.63%** | **18,726** |
| Guarded A100 pooling | 410/480 | 85.42% | 18,797 |
| Unguarded A20 pooling | 409/480 | 85.21% | 18,891 |
| Unguarded A100 pooling | 409/480 | 85.21% | 18,785 |

Coverage means the single full microscopic Garfield++ reference for that physical muon falls in the model's central 5th-to-95th percentile interval. **Nominal 90% is not achieved**, and a single reference realization is not the conditional mean of a muon. The mean absolute error is computed on the sampled predictive **mean**, not on matching stochastic avalanche samples event-by-event. Because geometries have very different gain scales, the pooled MAE must not be interpreted as an equal-error metric for each geometry.

#### Coverage by geometry (all have 120 parents)

| Geometry | A20 | A100 | Guarded A20 | Guarded A100 |
|---|---:|---:|---:|---:|
| p3_noFW | 86/120 = 71.67% | 91/120 = 75.83% | 98/120 = 81.67% | 99/120 = 82.50% |
| p3_FW75 | 99/120 = 82.50% | 97/120 = 80.83% | 105/120 = 87.50% | 104/120 = 86.67% |
| p4_noFW | 90/120 = 75.00% | 92/120 = 76.67% | 103/120 = 85.83% | 102/120 = 85.00% |
| p4_FW75 | 93/120 = 77.50% | 98/120 = 81.67% | 105/120 = 87.50% | 105/120 = 87.50% |

**Paired 120-parent-cluster bootstrap**, preserving the correlation between all four geometries of the same parent:

- Guarded A20 vs original A20: **+8.96 percentage points**, 95% interval **[+6.25, +12.08] pp**.
- Guarded A100 vs targeted A100: **+6.67 pp**, 95% interval **[+4.17, +9.38] pp**.

The tested heuristic exposed 22,825 of 23,524 primary×geometry instances to guarded neighborhood pooling; 699 were protected as local-only. The small difference between guarded and unguarded coverage on this sample **does not** establish safety around electrode surfaces or for induced signals. Total-ion sector 90%-interval coverage also improves but is not an induced electrical signal: overall hemisphere-asymmetry coverage A20 **85.0%**, guarded A20 **89.4%**; top/bottom-asymmetry coverage A20 **84.4%**, guarded A20 **88.3%**.

### Statistical and physical limitations

1. **Not external validation.** These same 120 muon primary coordinates informed the original target-cell ranking. Neighbor radius/weights were explored using the variance-targeted A100 sample. Their apparent performance may be optimistic, even though the actual full Garfield muon response was not used to select those cells.
2. **Conditional intervals only.** The 2,000 response replicas resample frozen finite-size response libraries. They do **not** include uncertainty in the per-cell library itself (the separate hierarchical bootstrap previously raised A20 coverage from ~77% toward ~85%). We therefore cannot assert that pooling is statistically calibrated simply because 85.6% is closer to 90%.
3. **Spatial guard is heuristic.** The actual source coordinate and donor avalanche coordinates can differ inside the 0.1×0.2 mm² bins; geometry-wire distance is computed only at bin centers. Physical field lines, collection efficiency, B field, and near-wire atypical tails deserve more targeted validation.
4. **Avalanche ion counts ≠ strip response.** No ion drift, Shockley–Ramo weighting field, cathode U/V strips or electronics are simulated. Preserving 5-dimensional angular census is useful but insufficient to construct real induced waveforms. For future full Garfield runs preserve ion birth positions and times / reservoir rather than counts alone.
5. **A100 + pooling adds little on these muons.** This may mean finite-cell noise or spatial mis-specification currently dominates, but should not be generalized to all unselected cells.

### Research decision

**Do not launch a days-long uniform A100/A200 simulation yet.** First do (a) a finite-library **training bootstrap** for the guarded pooled response using shared cell uncertainty, (b) inspect behavior close to anodes/field wires, and (c) reserve **new fully independent parent muons** (e.g. 100–120), generated from Geant4→Heed→full Garfield with independently frozen seeds and no re-tuning, to validate intervals and muon response.

If still useful, run an immutable targeted next 48,000 full microscopic avalanches to bring roughly 200 top-ranked cells per geometry to at least 100 samples; this is estimated to cost **~5–6 hours** at previously measured 0.35–0.46 seconds per avalanche. For all 4,200 cells to reach A100 requires 320,000 new avalanches (~31–41 wall hours if sequential), and to A200 requires 740,000 (~3–4 days). These are pilot-runtime estimates, not fixed production promises.

### Reproduce using the original WSL native archives

The new script is committed to:

`test/mwpc-stage-c/pi_spatial_lut_muon120_crn_v1.py`

```bash
cd ~/na6p/src/NA6PRoot-mwpc
git pull --ff-only
python3 test/mwpc-stage-c/pi_spatial_lut_muon120_crn_v1.py selftest
python3 -u test/mwpc-stage-c/pi_spatial_lut_muon120_crn_v1.py run --native --reps 2000
```

The script uses the project's existing native SHA-verified A100/120-muon reader and writes only under a new sibling output root `compare_frozen_120_real_muons_guarded_spatial_pooling_v1`. It refuses to overwrite existing outputs. The **archived-data evaluation has run and passed locally**, the distinct native WSL path has **not** yet been confirmed by the user.