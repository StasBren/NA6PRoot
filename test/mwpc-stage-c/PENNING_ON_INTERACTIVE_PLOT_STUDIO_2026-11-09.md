# Penning-ON MWPC interactive Plot Studio

**Purpose:** fast interactive scientific graph edits without repeated Python
code changes, re-running Geant4/TrackHeed/Garfield or re-reading original
ion-birth coordinate CSV files. Works locally in a browser via Streamlit
and Plotly. Source is exclusively the SHA-sealed Penning-ON 1,000 genuine
J/psi muons (p3_noFW, Ar/CO2 70:30, Ar Penning transfer 0.57).

The dashboard Python script is:
test/mwpc-stage-c/pi_p3_noFW_penning057_interactive_dashboard_v1.py

## Safe install using GitHub fetch (does NOT merge dirty older plotting files)

```bash
cd ~/na6p/src/NA6PRoot-mwpc
git fetch origin feature/mwpc-stage-c
SCRIPT=test/mwpc-stage-c/pi_p3_noFW_penning057_interactive_dashboard_v1.py
test ! -e "$SCRIPT" || { echo 'Already present: preserve/inspect existing script'; exit 1; }
git show "FETCH_HEAD:$SCRIPT" > "$SCRIPT"
PY=test_runs/mwpc_stage_c/.venv_pi_slides/bin/python
"$PY" -m py_compile "$SCRIPT"
"$PY" "$SCRIPT" selftest
"$PY" -m pip install streamlit plotly
```

**First prepare (once only):**

```bash
"$PY" -u "$SCRIPT" prepare | tee ~/na6p/penning_on_dashboard_prepare.log
```

The 'prepare' mode rechecks all 50 original Penning-ON SHA seals and
physics-log Penning markers, source identity, all exact per-seed
electrons / ion-birth counts / upper-lower exact asymmetry counters,
and for sampled ion births verifies that all reservoir weights
restore the full ion census. It produces a NEW compact compressed cache:

test_runs/mwpc_stage_c/pi_p3_noFW_jpsi1000_penning057_v1/
interactive_dashboard_v1/dashboard_data.npz
interactive_dashboard_v1/provenance.json

No original physics file, previous graph or backup is modified.
Never delete original 50 native checkpoints to save disk space:
the cache is derived and finite-resolution, not lossless physics.

**Launch browser app:**

```bash
"$PY" -m streamlit run "$SCRIPT" --server.address localhost --server.port 8501
```

Then visit **http://localhost:8501** from the Windows browser.
Stop the app with Ctrl+C; cached Penning data remain on disk.

## Plots (all Penning ON only)

1. Avalanche collected-electron gain, zero included: exactly-zero
   separate [0,1) bin, variable positive bin widths 10–500 electrons.
2. Avalanche collected-electron gain conditional on gain>0.
3. Ion birth azimuth, **equal weight per ion-producing avalanche**.
4. Ion birth azimuth, **ion weighted** via exact birth count/128 points.
5. Exact per-avalanche upper/lower ion-birth asymmetry A.
6. Exact per-real-muon hemisphere asymmetry A.
7. Physical ion-birth delay (relative to each initiating Heed electron's
   stored creation time), weighted by number of created ions.
8. Per-avalanche sampled birth-time span (>= lower bound; max 128 points).
9. Mean |A| versus avalanche collected gain (exact asymmetry counts).
10. Ion birth position heatmap (nearest-anode u offset and v, weighted).
11. Mean collected gain versus initial Heed electron starting u,v.
12. Ion birth delay versus azimuth 2D weighted joint distribution.

## Editable settings, without regenerating any physical data

- Title, legend/colorbar text, x- and y-axis titles.
- Linear/logarithmic axes where mathematically valid.
- X/Y limits, curve colour, legend visibility and legend position.
- 1D stepped curve or bars, curve width, font sizes.
- Gain bin width, angular bin width, asymmetry bins, time grouping.
- Per-bin probability, probability density or raw (weighted) counts.
- Heatmap linear or log colour scaling.
- Reset one chart's controls to defaults without modifying data.
- Plotly interactive hover/zoom, camera PNG/SVG/JPEG export.
- Download self-contained HTML plot or plotted histogram CSV.

Different plot types retain independent controls inside the running
Streamlit session. 'Reset this plot' restores its defaults.

## Scientific caveats

- The angular and spatial birth plots are based on at most 128
  uniformly drawn positions per avalanche and real ion-count weights.
  They are NOT exact per-ion trajectories or later drifting-ion clouds.
- The asymmetry is exact (callback ion birth counts in upper/lower sectors).
- The starting electron is a genuine TrackHeed primary. No real cathode
  induced-current waveform has yet been simulated (Shockley–Ramo, ion
  mobility/drift, electronics convolution remain future steps).
- 'Time' means physical Garfield ion-birth time minus the initiating
  Heed electron start time, in nanoseconds, NEVER wall-clock CPU timing.
- Time bins are intentionally nonuniform; use probability density when
  comparing bin heights. A zero-delay bin exists separately and thus
  log-X cannot show it without a scientifically meaningful exclusion.
- Zero-electron gain events still exist with Penning enabled.
- All observations refer to 3mm symmetric p3_noFW geometry at 1800V.
  This is not the 4mm asymmetric (2+4mm) 2026 prototype.

## If preparation stops

Do not restart Garfield. Read the exact Python exception. The exporter
intentionally rejects missing SHA seals, capped avalanches, unpaired
starting electrons, inconsistent ion-birth reservoir sizes/weights,
and impossible negative creation-time offsets. Any such failure needs
investigation rather than hiding bad data in a graph.


## New: physics interpretation above every plot

The dashboard now displays an expanded "Understanding this plot" panel
**immediately above** the selected figure (all 12 modes). Each panel
gives:
- What the axes and the source distribution actually represent.
- How Garfield, TrackHeed, exact sector counters or the 128-per-avalanche
  random ion-birth reservoir generated those numbers.
- A derived observation from the loaded, verified Penning-ON cache
  (e.g. mean/RMS/zero fraction, signed asymmetry/mean |A|, angular
  axial-to-transverse contrast or approximate ion-weighted time
  percentiles), with limitations stated explicitly.

These paragraphs use the physical cached data; they do not claim to
have reconstructed the induced cathode waveform or individual drift
trajectories. The interpretation is outside the downloadable figure,
so slide-quality PNG/SVG exports remain uncluttered. Uncheck
"Show physics interpretation above plot" in the sidebar to hide it.

**No need to re-run 'prepare'.** This addition changes ONLY the
local Streamlit presentation. Re-fetch the one dashboard Python file,
run selftest and restart Streamlit. To force a white browser theme:

```bash
"$PY" -m streamlit run "$SCRIPT" --server.address localhost --server.port 8501 --theme.base light
```

