# Full J/ψ source dashboard: two million generated parents

The ongoing six-million-parent source campaign already completed 40 frozen J/ψ 50k-parent native Geant4/Heed/Garfield folds **for all six MS stations** (240 station/chunk shards). This independent dashboard publishes those **finished** numerical results without starting more transport, changing the original 20-trial Garfield LUTs, modifying the currently-running Phi/Omega shards, or altering original 100k comparison outputs.

## Read-only dashboard builder

```bash
cd ~/na6p/src/NA6PRoot-mwpc
git pull --ff-only
bash test/mwpc-stage-c/run_pi_full_jpsi_2m_dashboard.sh check
bash test/mwpc-stage-c/run_pi_full_jpsi_2m_dashboard.sh start

# Monitor at any time:
bash test/mwpc-stage-c/run_pi_full_jpsi_2m_dashboard.sh status

# After all six stations complete:
bash test/mwpc-stage-c/run_pi_full_jpsi_2m_dashboard.sh verify
explorer.exe test_runs/mwpc_stage_c/pi_full_jpsi_2m_dashboard_v1/pages
```

Open `pages/index.html` locally (ideally Firefox/Chrome/Edge; pages embed their physical data and work offline). Each of the six station links contains an interactive XY physical chamber map, numbered anodes, full 10 ns birth-to-anode electron drift PDF, signed muon entrance angle slices θu/θw, Garfield last-step approach angle relative to the inward wire normal and wire-surface impact-position azimuth φ (not a direction azimuth), and six pitch/field-wire geometries. There is also a detailed spectrum dashboard for every station, with logarithmic Y and time/angle selections.

## New XY dashboard controls: multiple DETs, six histogram panels and rebin

The existing `MS0_layout_v3.html` ... `MS5_layout_v3.html` and their
`*_en.html` siblings now support **additive physical-chamber selection**.
Click a chamber on the XY map to select only it. Hold Ctrl, Cmd or Shift
and click to add/remove chambers; alternatively open the checkbox list
under the map. Select all, clear all and layer A/B/C/D filters are available.
An empty chamber selection produces no data rather than silently reverting
to the entire station. Single-DET physical-anode selection remains enabled;
for multiple DETs the anode selector is disabled to avoid mixing distinct
physical wires with the same index.

Each selected chamber contributes its saved chamber-level counters and
histograms **exactly once**. Primary births, gas crossing segments,
collected electrons, 31-bin time histograms, 18-bin electron incidence
histograms, 36-bin φ_hit histograms and time moments are summed.
The aggregate drift mean and centered RMS are recalculated from native
`sum_drift_ns` and `sum_drift_ns2`, not averaged by detector. The native
10 ns timing histogram retains the >=300 ns overflow in its denominator.

The time plot modes are: six separate histograms (default), six overlaid
lines, six stepped histograms, and bars for one selected configuration.
The electron arrival angle plot also has six separate histograms (default),
for both last-step incidence to normal and wire-surface φ_hit.
All six mini-panels use the same Y range; log-Y is supported for time.
Grouping can only use **native or coarser bins**: 10/20/50 ns for drift
time, and 10°/20°/30° for electron angles. No fake 5 ns or 5° interpolation.

To upgrade a previously built/verified 2M dashboard, repeat the regular
`run_pi_full_jpsi_2m_dashboard.sh start`. If the station summary and
merged JSON pass SHA256 verification, the runner re-renders **only the
outdated layout HTML** and recomputes its publication SHA, without
rerunning the original 240 Geant4/Heed/Garfield fold checkpoints or
modifying `merged/jpsi/MS*.json`.

## Wire-surface electron impact-position φ_hit: every station and DET

The main dashboard `pages/index.html` now links to a six-station polar atlas.
Each station page `pages/jpsi40-bore/MS0_phi_atlas.html` through `MS5_phi_atlas.html`
includes *every physical chamber* and an `Entire station` option. Chambers
without finite-φ collected electrons remain visibly empty. Six polar histograms
compare 3, 3.5, and 4 mm pitches with and without the grounded 75 μm field wire.
The selected Geant4 muon entrance angle can be all, signed θu, or signed θw.
The common radial scale is normalized by the finite-φ electron count in
each selected configuration and angle bin. The page also displays the
field-wire vs no-field-wire histogram total-variation distance and the
separate absolute anode collection fractions.

The machine-readable `phi_atlas_stats.csv` contains results for every
(station, DET, pitch) combination, including station-wide sums. New atlas
pages and CSV are included in the complete SHA256 manifest and Drive export.
Raw physical-anode groups are excluded to avoid double counting; the
per-chamber and per-angle phi histograms must conserve station-wide counts.
Run the normal dashboard `start` command again to add this atlas to a
previously built dashboard: already SHA-verified station JSON and legacy
HTML are reused without rerunning Geant4, Heed, Garfield, or frozen folds.

φ_hit is the electron ENDPOINT POSITION on the anode circumference,
not the electron velocity direction or incidence angle. The 2M value counts
generated J/ψ parents, not independent electron tracks. No parent-event
bootstrap error bars are claimed, and the original B=0 / 20-trial-per-cell
Garfield limitations remain.
Both **English** `*_en.html` and Russian `.html` variants are generated; the dashboard landing page and numerical CSV are in English. Existing source/station/DET shapes and physical anode numbers are used directly; no interpolation or modeled curves are substituted.

## How the science is verified

1. All 240 original chunk/station SHA256 manifests are required, with exactly 50,000 **generated parent events per Geant4 chunk**, correct J/ψ source and physical station IDs, complete group counts and saved gzip hashes.
2. On build, all source gzip files are SHA256-verified and their original per-chunk full Heed/Garfield six-LUT fingerprints checked, including the frozen 20-independent-microscopic-trials-per-box Garfield provenance.
3. All raw histogram counts and *unrounded electron moments* are summed; no re-sampling. Per-anode → chamber → station time histograms and both signed muon-angle bin sets are conserved across the entire source.
4. For each station, one merged JSON, four independent Russian/English HTML pages and SHA-pinned compact statistics are retained in the new, isolated dashboard folder. A successful rerun verifies finished outputs and skips them.
5. `verify` checks all six station JSON/HTML SHA256 digests, the 40 input manifest SHA256 values for each station, the portal, φ_hit atlas pages, and the summary CSV files.

## Data scope

**2,000,000 generated J/ψ parent events** is the denominator of source event generation, *not* the number of recorded muons, station gas crossings, independent Garfield trajectories or anode electrons. Separate counters show:

- Real native muon gas-crossing segments per station (the same muon can cross multiple stations).
- Heed primary electron births per station.
- For each of six simulated pitch/field-wire geometries: anode-collected electrons, collection fraction, known drift time, mean and centered RMS.
- Original **20 truly independent microscopic Garfield drift trials per spatial cell** reused statistically for the Heed electron population. The **new 50-trial LUT run is separate** and does not retroactively alter this dashboard.

Limitations: no field on (B=0 analytic reference); no avalanche, signal induction, electronics threshold, or parent-event bootstrapped uncertainties; the 75 μm grounded field-wire diameter is only a study assumption.

## Publishing in Google Drive

After `verify` passes:

```bash
bash test/mwpc-stage-c/run_pi_drive_full_jpsi_2m_dashboard.sh check
bash test/mwpc-stage-c/run_pi_drive_full_jpsi_2m_dashboard.sh copy
bash test/mwpc-stage-c/run_pi_drive_full_jpsi_2m_dashboard.sh verify
```

Publishes only standalone HTML pages (including φ_hit atlases), the portal, numerical summaries, phi_atlas_stats.csv and a complete manifest to `na60drive:test_runs/mwpc_stage_c/pi_full_jpsi_2m_dashboard_v1`. Large native ROOT files, original Heed fold .gz checkpoints and merged source JSON remain on local WSL. HTML previews inside Google Drive might not execute JavaScript: download/open the HTML files locally for full interactivity.
