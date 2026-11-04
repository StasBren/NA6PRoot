# NA60+ MWPC — blinded new120 fixed-guard LUT validation (2026-11-04)

**Actual WSL native physics and blind comparison completed.** Original 120 parent tracks were excluded by source (source/chunk/event) before selecting 120 different original Geant4 J/psi parent events from the already-generated, SHA-validated chunk-2 source. These are newly evaluated parent events and new microscopic Garfield++ histories, not newly generated Geant4. Selection was balanced across 6 stations × two signed incidence-angle strata, not flux-unbiased.

- Reference: 120 new parents, 2,580 TrackHeed clusters, 6,193 accepted electrons, 24,772 full AvalancheMicroscopic calls over four B=0 geometries; 12/12 original SHA-sealed 10-parent blocks. Physically in-metal zero-ion anode starts in one completed block were strictly verified without rerunning physics. `REFERENCE_SEALED True`.
- Training: original immutable A20 84k and A100 targeted 16k single-electron true microscopic avalanche pools. Frozen *before seeing new120 outcomes*: protected 5×5 local spatial pooling; 25% target-cell empirical pool / 75% Gaussian-weighted neighboring original A20 pools; no crossing v=0, wire-proximity veto and near-wire sign safeguards.
- Analysis: 3,000 shared-global-cell-training bootstrap iterations, 480 linked (120 physical parent × four geometry) observations; all model comparisons use exactly identical new parent truth. This bootstrap approximates finite LUT training uncertainty; it is not a fully calibrated posterior.

| Metric | A20 | A100 | guarded A20 | guarded A100 |
|---|---:|---:|---:|---:|
| Nominal 90%-PI covers | 420/480 | 416/480 | **432/480** | **434/480** |
| Coverage | 87.500% | 86.667% | **90.000%** | **90.417%** |
| MAE of predicted mean ion births | 21,746.464 | 21,612.590 | **18,480.082** | **18,467.451** |
| Mean nominal 90%-PI width | 87,276.607 | 86,259.070 | **78,803.038** | **78,668.103** |
| Hemisphere asymmetry interval covers | 86.250% | 86.042% | **89.167%** | 88.958% |
| Top/bottom asymmetry interval covers | 88.958% | 88.750% | **89.167%** | **89.167%** |

Prespecified 120-parent-cluster bootstrap 95% percentile CI for changes in empirical ion birth PI coverage:
- Guarded A20 minus plain A20 = **+2.50 percentage points**, **[+0.2083, +4.7917]** pp.
- Guarded A100 minus A100 = **+3.75 pp**, **[+1.6667, +6.0417]** pp.

Guarded A20 lowers predictive mean MAE by **15.02%** and mean PI width by **9.71%** relative to A20, *while increasing coverage*, so this is not trivial interval inflation. Guarded A100 offers only 2 additional covered cases / 480 (+0.417 pp), 12.63 fewer absolute mean-error ions and 134.94 fewer nominal PI width ions than guarded A20; **no compelling incremental benefit** demonstrated for A100 targeted sampling atop pooling.

**Interpretation:** Blinded, source-disjoint parent validation supports protected spatial A20 pooling as the leading fast mean/ion-sector avalanche emulator for this detector model. Observed 90% coverage does **not** prove exact calibration on arbitrary flux, B, gas or wire geometries: all 4 configurations for each parent are correlated, new parents use the same Geant4 source/generator and a balanced signed-angle/station selection, and data reflect only one true microscopic Garfield draw per Heed primary. Two paired comparisons and many output metrics are exploratory; do not inflate significance. Global cell training bootstrap is approximate. No induced U/V strip signal, ion drift, weighting field or full ion birth spatio-temporal distribution has been validated.

**Recommendation:** Keep **guarded A20** as primary reference and guarded A100 as robustness comparison; do **not** spend days expanding every LUT cell at present. Next run read-only diagnostics on immutable 480 per-parent rows: per geometry, parent-clustered paired MAE and interval width, source exposure to A100 target cells, near-wire and tail/zero checks. Then implement/validate a physically consistent ion birth cloud, ion drift and Shockley–Ramo induced U/V readout model. Unmodified C++ physics, A20/A100 archives, original old120 reference and new120 reference remain SHA sealed.

Authoritative new WSL files: `test_runs/mwpc_stage_c/pi_external_new120_guarded_lut_blind_v1/new_full_garfield_120/manifest.json` and `.../blind_evaluation_fixed_guard_v1/{per_external_muon_4_models.csv,summary.json,manifest.json}`. All output numbers above were pasted by user from the actual WSL run; assistant has **not** independently inspected the 480-row CSV and does not claim any per-geometry subgroup findings yet.
