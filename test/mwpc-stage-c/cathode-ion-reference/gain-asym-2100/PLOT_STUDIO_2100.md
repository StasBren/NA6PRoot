# NA60+ Plot Studio — 2100 V, up to 997 archived Garfield muons

Read-only interactive plotting dashboard for original 2100 V Garfield++ microscopic statistics, restoring the editing workflow of our original Plot Studio. It reads the **C++/ROOT V0 census CSVs** and never reruns Garfield, ion drift, or cathode response. All plot output is generated interactively. Scientific data remain immutable.

## WSL commands

From your existing checkout `~/na6p/wt-gain-asym-2100`:

```bash
DASH=test/mwpc-stage-c/cathode-ion-reference/gain-asym-2100/run_plot_studio_2100.sh
bash "$DASH" install
bash "$DASH" check
bash "$DASH" run
```

Open **http://localhost:8501** in Windows, Ctrl+C to stop. If that port is occupied: `bash "$DASH" run 8502`.

Python packages are isolated in `~/na6p/venvs/na60-plotstudio-2100` so your existing CERN ROOT/Garfield installs remain untouched.

The dashboard automatically discovers and offers all complete `pilot_*` and `population_*` outputs under `test_runs/mwpc_stage_c/2100_gain_asym/`. If you only have the five-muon pilot, the UI prominently labels it a software test, NOT a population result. To prepare full statistics, in a **separate terminal** run the already committed C++/ROOT analysis:
```bash
RUNNER=test/mwpc-stage-c/cathode-ion-reference/gain-asym-2100/run_gain_asym_2100.sh
bash "$RUNNER" population
```
After that completes, refresh the dashboard and select the `population_*` dataset.

## Available sections

1. Gas gain: per-primary ion births, per-muon ion and collected-electron gains, total counts and 2D gain correlations.
2. Physical anode charge sharing: collecting-wire multiplicity, second-wire fraction f2, P(f2>0.05) vs angle and f2 vs ionization gain. Muon-level collected-wire census closure is checked before using each event.
3. Avalanche angular structure: K256-estimated circular R and sigma_phi versus gain, and their 2D population distributions. These are sampled angle *moments*, not exact full-ion angular density.
4. Along-wire structure: exact microscopic within-fragment longitudinal RMS sigma_w versus gain, angle and asymmetry, not full-muon longitudinal spread.
5. Up/down asymmetry: signed A and mean |A| versus exact ion gain, impact position, signed angle, plus 2D densities at primary and muon levels.

For every plot, edit title, horizontal/vertical names, axis scales and limits, histogram binning/normalization, color, line/marker width, fonts, legend and muon-block bootstrap bars. The explanation panel details what is plotted, exactly how computed, and limitations. Export standalone interactive HTML, displayed-data CSV, and presentation-quality PNG or vector PDF (Kaleido).

## Verified conventions

```
phi = atan2(v, u - nearest_anode_u)
up: 0 < phi < pi; down: -pi < phi < 0
A_UD = (N_up - N_down)/(N_up + N_down)
G_ion(muon) = ion births_exact / Heed primaries
G_collected(muon) = collected electrons_exact / Heed primaries
mean G across all muons = sum of counts / sum of Heed primaries
```

Births on the dividing axis are excluded only from A's denominator, not ion yield. A is undefined if N_up+N_down=0; zero-ion primaries remain in gain histograms. Exact A uses original Garfield per-wire counter census in `microscopic_wire_seed.csv`; K256 sampling is used only where identified as such. Each avalanche or fragment gets one observation in binned trends; errors use 160 **independent muon-block resamples** to preserve correlations between primaries from one muon.

**All 997** includes the previous sealed validation subset: valid for descriptive histograms, NOT an independent holdout for model calibration. Idealized 2100 V model: 3-mm anode pitch, 30-µm diameter, cathodes ±3 mm, no field wires, Ar/CO2 70/30, Penning transfer 0.57, B=0. Not the as-built prototype geometry.

## Missing inputs are identified, not guessed

The compact V0 exports do not contain exact all-ion azimuth histograms, full primary-relative delta-w distribution, complete drift trajectories for all 997, U/V weighting potentials, induced strip signals or collected cathode charge. The dashboard does not synthesize these. A separate ion-birth ROOT export or positive-ion drift simulation would be needed.

**Troubleshooting:** If `check` finds no outputs, run the five-muon C++ pilot first; if static PNG/PDF fails, HTML and CSV remain downloadable. Errors from pip can be resolved within the isolated venv without reinstalling CERN ROOT.
