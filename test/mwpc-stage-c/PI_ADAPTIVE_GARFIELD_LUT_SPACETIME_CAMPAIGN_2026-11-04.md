# NA60+ MWPC — staged microscopic Garfield++ campaign and replacement LUT

**Date:** 2026-11-04. This is a *new, separate* research campaign. NO native Garfield jobs have been executed by ChatGPT. The script passes Python syntax checks and a **synthetic selection/chunking** selftest in the analysis environment; C++ compilation, source availability, Geant4/Heed preflight, and physics run on the user's WSL **remain to be verified**.

## 1. GitHub files added

- `test/mwpc-avalanche/phaseB_native_muon_3p3_spacetime_v1.cxx`: exact copy of previous `phaseB_native_muon_3p3_avalanche.cxx` with only the independent-random-reservoir cap raised from **5 to 128 birth positions/times per actual primary-electron avalanche** in `--muon-cache-physics`; preserves full ion sector counting and the original C++ target. Coordinates `(u,v,w,t)` are weighted by `N_ion_births / n_saved` from the reservoir. The physical ionization process and RNG engine are not intended to change; validate whole-census identity with a matched source test if necessary.
- `test/mwpc-avalanche/CMakeLists.txt`: adds only new `mwpc_phase_b_native_muon_3p3_spacetime_v1` build target. Original target remains.
- `test/mwpc-stage-c/pi_adaptive_muon_spacetime_campaign_v1.py`: fixed source selection, 10-parent physics chunks, strict anode/FW in-metal zero-ion validation via original trusted bytecode with isolated augmentation, file hashes, staged paired parent-cluster bootstrap report.

## 2. WSL commands (first 10 physical muons only)

```bash
cd ~/na6p/src/NA6PRoot-mwpc
git switch feature/mwpc-stage-c
git pull --ff-only
python3 test/mwpc-stage-c/pi_adaptive_muon_spacetime_campaign_v1.py selftest

# Reuse already-installed native Garfield++/Heed environment.
source "$HOME/na6p/env.sh" 2>/dev/null || true
source "$HOME/na6p/install/garfieldpp/share/Garfield/setupGarfield.sh" 2>/dev/null || true
cmake -S test/mwpc-avalanche -B test/mwpc-avalanche/build \
  -DGarfield_DIR="$HOME/na6p/install/garfieldpp/lib/cmake/Garfield"
cmake --build test/mwpc-avalanche/build \
  --target mwpc_phase_b_native_muon_3p3_spacetime_v1 -j2

# Freeze source and the reserved validation parent IDs *before* physics.
python3 -u test/mwpc-stage-c/pi_adaptive_muon_spacetime_campaign_v1.py prepare --n 2400

# 10 real Heed-only tracks but no full microscopic avalanches.
python3 -u test/mwpc-stage-c/pi_adaptive_muon_spacetime_campaign_v1.py preflight

# First 10 independent muon parents *only*, all 4 fields; inspect results.
python3 -u test/mwpc-stage-c/pi_adaptive_muon_spacetime_campaign_v1.py run --max-chunks 1 --hours 8
python3 test/mwpc-stage-c/pi_adaptive_muon_spacetime_campaign_v1.py status
```

**Do not initiate the full multi-day campaign before inspecting the first chunk**. The new target is not yet WSL-compiled/tested. If a block fails validation or native physics, the original raw output and logs are left in a unique, **unsealed** directory, and `run` refuses to redo or skip it. Inspect and repair using a new explicit protocol, not `rm -rf`, and never silently accept non-metal Garfield failures. Do not edit committed Python/C++ while a prepared campaign with code/SHA locks is active. Do not run `prepare` twice. No Google Drive upload is performed.

If the first chunk PASSes, resume **in separate bounded blocks**:

```bash
# For example, finish 120 development muons total by running 11 more chunks.
python3 -u test/mwpc-stage-c/pi_adaptive_muon_spacetime_campaign_v1.py run --max-chunks 11 --hours 8
python3 -u test/mwpc-stage-c/pi_adaptive_muon_spacetime_campaign_v1.py report --reps 2000

# Later: extend to 600, 1200, 1920 development muons without changing input selection.
# 'max-chunks' = number of additional NEW 10-muon blocks for that invocation.
# After development is frozen, finish the final 480 reserved parents if justified.
```

With `--n 2400`, parents 1–1920 comprise the development cohort; 1921–2400 are reserved as an independent event-level validation set. Ordinary `report` restricts itself to first 1920 even if later blocks exist. The final reserved report is locked until all 240 chunks have completed and explicitly requested using `report --final-holdout`. Once the final 480 are inspected, do **not** tune the model on them and then call them independent validation.

All original old120 + new120 **240 parent event IDs are excluded** and SHA-verified. Exactly one original J/ψ parent event is used per new selected muon. All four geometries share each muon's original Heed primary list. Source: already-generated, frozen Geant4 J/ψ `chunk0002`. Geometry B=0, 1800V, p=3/4mm, with/without FW75, cathode ±3mm, Ar/CO2 70/30, T=293.15K, P=760 Torr. The chosen sample deliberately has equal station×signed-angle populations and **does not have the native source-angle distribution**; the report presents both equal-stratum estimates and a source-eligible-stratum-count-weighted **proxy**, not a truly unstratified physical flux sample. The latter needs a separate source-level selection with acceptance/weights.

## 3. Statistical stopping rule (physics descriptive statistics, not LUT validity)

The script computes independently within each geometry:
- Equal-parent-weight mean collected effective gain (`sum collected avalanche electrons / Heed primary electrons` per muon), plus separately the pooled electron-weighted gain.
- RMS of **whole-muon ion-birth** `A_hemi` and `A_TB` (undefined values left missing), not the primary-electron endpoint asymmetry.
- Paired `RMS(A_TB)_FW75 - RMS(A_TB)_noFW` at both pitches.
- Prespecified parent-bootstrap 95% CI, resampling source parent IDs **within station×angle stratum**; four geometries always remain paired. Also stratum-reweighted source-proxy moments, effective sample size and exploratory gain q50/q90/q95/q99.

Initial **illustrative** decision tolerances (change ONLY before freezing a new campaign): relative 95%-CI halfwidth **≤10%** of mean effective gain in each geometry; 95%-CI halfwidth **≤0.03** on both asymmetry RMS; 95%-CI halfwidth **≤0.03** on both paired FW-minus-noFW `A_TB` RMS differences; require both balanced and source-proxy estimates. These are *not* thresholds delivered by Sasha or experiment resolution. Once RMS/gain converge, require **additional** criteria for zero-ion fraction, large-tail quantiles, near-wire sectors, theta/station slices, absolute ion cloud shape and induced U/V waveform. The script intentionally does not claim they have converged.

For a mean with standard deviation sigma and target 95% halfwidth delta, rough independent-parent count is `n ≈ (1.96 sigma/delta)^2`. For a rare tail of probability p, N independent parents yields only N*p expected tail events; ~100 events above a 1% threshold require ~10000 parents *before* further angular/station stratification. Parent clustering, source reweighting, non-Gaussian gain and true Garfield stochasticity can all make simple estimates optimistic. Bootstrap convergence versus `n=120, 600, 1200, 1920` matters more than choosing 2000 arbitrarily.

## 4. What a true ion-cloud / cathode LUT must be

A20/GuardedA20 exists now for a **joint 5-number ion sector census** `(N_birth, N_up, N_down, N_top, N_bottom)` conditional on 0.1×0.2mm spatial cells, not for charge readout. The physics response from one Heed primary at true `(u,v,w,t)` should be represented by one **jointly sampled** microscopic outcome:

`R = (N_e_collected, N_ion_births, N_up, N_down, N_top, N_bottom, {weighted birth (delta_u,delta_v,delta_w,delta_t) points}, electron endpoints, failure/cap flags)`.

The initial Heed primary's own ion and prompt electron transport also contribute to the final induced signal. Shift a canonical periodic cell response to the physically struck anode, with longitudinal z/w and time referenced to that primary's original start; treat FW gas/metal and boundary cells separately. Preserve the stochastic dependence of gain and spatial/time pattern; do **not** independently sample a scalar Polya gain and angular bins or multiply an average shape by random gain without checking that factorisation against Garfield.

For each exact Geant4 muon crossing: generate genuine Heed clusters/electrons once; for each Heed electron select from the conditional joint-response library for the same field and same local seed coordinates; combine full ion populations across the primary electrons, then drift ions in the correct field and apply each cathode strip's weighting potential. The current 128-point equal-probability reservoir is a **candidate compressed representation**, not proven accurate for rare local signal features; test against a small exact unsampled/full-birth reference or an equivalent high-resolution binned ion-cloud census. Keep 0-ion branch and rare large avalanches. Trained LUT trials and validation muon full Garfield **must not share** exact stochastic outcomes; independent parent selection is necessary but does not automatically validate each under-sampled spatial cell.

Uniform 100/4200 cells across four fields requires **320k new** actual microscopic avalanches beyond the already existing 100k, typically 31–41 hours based on actual earlier runtime and excluding enlarged reservoir I/O. But rather than blindly simulating every cell to the same depth, allocate more samples to conditionally important/unstable cells and sparse FW/near-wire sectors according to *ion cloud / readout error*, not solely total ion-birth variance. Require separate holdout at the **cell/sample** level and **full parent** level. For example, if gain coefficient of variation per cell is near 1, obtaining a 10% 95%-CI halfwidth on mean gain of one cell by direct sampling alone takes ~385 independent avalanche trials, not 20. Cloud-level response may be more difficult and remains unvalidated.

**Acceptance of replacement LUT:** held-out full Garfield should agree in zero-ion rate, mean and tails of collected gain, full ion birth sectors and angular distributions, and measured error/coverage conditional on angle, seed phase, station, wire proximity and selected detector geometry; **ultimately** use identical Heed primary starts, then independently propagate exact Garfield cloud and LUT reconstructed cloud through the **same physical ion drift, U/V strip weighting-field Shockley–Ramo, time-integration/electronics**. Compare induced strip charge, charge sharing, centroid bias, time pickoff, hit efficiency, threshold crossing and tails. If these don't agree within physically motivated detector tolerances, the count/sector LUT is not yet an adequate strip-signal replacement.

## 5. Current limitations

- Actual C++ WSL compilation, 10-event preflight and Garfield jobs have **NOT** been executed by the assistant; they require user's WSL.
- The separate 128-point reservoir per avalanche uses independent sampling RNG but must be checked for retained cloud fidelity; not enough for arbitrary waveform claims.
- Seed `(u,v,w,t)` and field geometry models are fixed to the existing small periodic 9-anode analytic patch; physical chamber edges, B field and gas-condition uncertainties are NOT validated.
- Stage results are source-/geometry-specific; no prototype gain/ion mobility/Penning/field-wire voltage calibration. Sample size does not remove model systematic errors.
- A statistically converged full Garfield parent ensemble answers **what the microscopic model predicts**, not necessarily what the actual chamber measures. Reliable replacement requires a *signal-level* physical reference and ideally measurements.

## 6. Verified nonphysics tests and PI figures

The synthetic 120-parent test is available as `test/mwpc-stage-c/test_pi_adaptive_muon_spacetime_campaign_v1.py`. It checks independent parent selection, 12 strata, source-proxy weights, paired four-geometry bootstrap, gain/ion-sector metrics and immutable report files, **without running Garfield**. Run after `git pull`:

```bash
python3 test/mwpc-stage-c/test_pi_adaptive_muon_spacetime_campaign_v1.py
```

After 12 or more SHA-valid completed 10-parent native Garfield blocks, run `report` and then create three PI-facing real-muon figures, strictly from the **development cohort** (reserved validation never read):

```bash
python3 -u test/mwpc-stage-c/pi_adaptive_muon_spacetime_campaign_v1.py report --reps 2000
python3 -u test/mwpc-stage-c/pi_adaptive_muon_sasha_figures_v1.py build
```

Outputs: `test_runs/mwpc_stage_c/pi_adaptive_muon_spacetime_v1/plots_development_NNNN/{01_mean_effective_gain_CI.png,02_true_muon_ion_asymmetry_distributions.png,03_weighted_ion_birth_azimuth.png,manifest.json}`. The third plot is normalized **per ionization birth via inverse-inclusion reservoir weights**, not per muon and not an induced cathode-charge observable. The first two concern physical muon crossing outcomes. These plots are research illustrations and should be cross-checked against the independently validated ion-sector count summaries.

The report also contains *exploratory*, parent-bootstrap estimates for `ion_gain` (new avalanche ionization births per Heed primary), `zero_frac` (zero **collected avalanche electron** fraction, not necessarily zero ion births) and `cathode_sector_frac` (fraction of created avalanche ions born in top/bottom ±45° sectors), each for balanced-stratum and approximate source-stratum weighting. Signal readiness **is not assessed** by these proxies.
