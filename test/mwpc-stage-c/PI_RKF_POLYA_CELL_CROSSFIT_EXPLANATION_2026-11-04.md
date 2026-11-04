# NA60+ MWPC — RKF + zero-inflated Pólya against physical A20 (4 Nov 2026)

This is an **exploratory count-only test**, not a replacement for the current
SHA-frozen `guarded_A20`. All source files stay immutable.

## What is being tested?

For an exact electron start with geometry g and initial position (u,v),
`DriftLineRKF` computes a deterministic/macroscopic estimate `m_RKF(u,v,g)`
using drift-line integration with Townsend and attachment coefficients.

The candidate is a **new, explicit statistical model**, NOT a direct native
claim that `DriftLineRKF::SetGainFluctuationsPolya` correctly counts microscopic
ion births:

- probability `p0_g`: zero microscopic ion births;
- otherwise `Y | (Y>0, u,v,g) ~ Gamma(k_g, mean = c_g*m_RKF/(1-p0_g))`.

The model thus has `E[Y|u,v,g] = c_g*m_RKF(u,v,g)`, when `m_RKF>0`.
It is a minimal *hurdle Gamma/Pólya surrogate for ion_births*, NOT a
simulated full avalanche or a distribution of collected electrons. The
Pólya shape convention is `theta = k-1` for positive gain; parameters
are **fitted**, not guessed. We do not claim the Gamma law holds for every
position or for the rare tails.

Controls:
1. `A20`: archived 20 original microscopic Garfield ion counts / local cell.
2. `geom_Polya`: same Gamma+zero fit per geometry, **without RKF**.
   If RKF does not outperform this control, local RKF shape information
   adds little predictive value.

## Existing physics — no new simulation

- Original A20: **84,000** actual microscopic avalanches / 4,200 cells.
- Selected A100 extra: **16,000** microscopic avalanches in 200 cells,
  80 new starts/cell. These 16,000 positions are exactly matched to the
  existing completed `DriftLineRKF` V3 run.
- The old A100 target-cell selection depended on A20 and earlier 120
  muons. Therefore these 200 cells are not a random sample of the full
  detector and cannot establish generalization to new Geant4 parents.
- Original 120 + independent newer 120 muon reference and their physical
  data are **not used for fitting** and not rerun.
- A reader verifies the physical archives' SHA seals and every original
  sample identity. No source C++ changes.

## Validation (no data leakage)

For each of four geometries, deterministically partition the 50 target
cells into 25 training and 25 validation cells. Fit `(p0,c,k)` using
the 80 microscopic new-A100 outcomes in each training cell. Test on the
80 complete new-A100 outcomes in each withheld cell. Swap halves and
combine scores. **No cell's 80 new micro responses enter its own RKF
fit**, but old A20 responses from that same cell may be used as the
A20 baseline, consistent with the original20→new80 independence.

Report per geometry, with 95% **cell-clustered** bootstrap:
- `cell_mean_absolute_error`: |mean 80 truth - mean 80 predictions|,
  then average over validation cells (not a single-draw error);
- `one_draw_MAE`: exploratory, has irreducible avalanche randomness;
- `zero_Brier`: mean squared error of zero-probability forecasts;
- `PI90_coverage`: fraction inside nominal [q05,q95];
- `upper95_exceedance` and `upper99_exceedance`: upper-tail calibration;
- `q95_pinball` and `q99_pinball`: proper quantile scores (smaller is
  better); `PI90` alone is misleading when ~17% true avalanches are zero.
- paired delta for RKF+Pólya minus both controls, bootstrap by CELL.

**No rule declares RKF good from one global mean ratio or from nominal
coverage alone.** If it does not decisively improve proper scores on
unseen cells, there is no reason to replace guarded_A20.

## Run on WSL

```bash
cd ~/na6p/src/NA6PRoot-mwpc
git switch feature/mwpc-stage-c
git pull --ff-only

PY=test_runs/mwpc_stage_c/.venv_pi_slides/bin/python
"$PY" -m py_compile test/mwpc-stage-c/pi_rkf_polya_cell_crossfit_v1.py
"$PY" test/mwpc-stage-c/pi_rkf_polya_cell_crossfit_v1.py selftest

# Existing V3 forensic audit; only if not run previously.
# This writes one separate new sibling summary, no simulation.
# If already run, SKIP and do not overwrite.
# "$PY" test/mwpc-stage-c/pi_audit_rkf_v3_16k_status.py audit

"$PY" -u test/mwpc-stage-c/pi_rkf_polya_cell_crossfit_v1.py audit --bootstrap 3000
```

The new program never invokes Garfield++. It creates **only**:

`test_runs/mwpc_stage_c/pi_RKF_Polya_cell_crossfit_v1/`

with `summary.json`, `fold_parameters.json`,
`per_cell_results.json` and `manifest.json`. Existing output
directories are never overwritten.

If the specialized venv path does not exist or does not contain NumPy
and SciPy, use an existing Python interpreter with both. Do not install
a new Garfield++ version just to run this analysis.

## Next experimental decision

Regardless of whether count-only RKF+Pólya is competitive, it **cannot**
produce the real conditional ion birth space/time cloud. That is why
the joint-response microscopic bank remains the current preferred
development path for a future cathode-ready avalanche simulator.
Before allocating multi-day full Garfield computation, first validate
the 128-point reservoir against an unsampled microscopic birth-cloud
reference, and test that ion count/angle/space/time are sampled *jointly*
per avalanche. Do not unfreeze or modify ongoing 2400-parent source and
its one complete SHA checkpoint. Keep `guarded_A20` frozen.
