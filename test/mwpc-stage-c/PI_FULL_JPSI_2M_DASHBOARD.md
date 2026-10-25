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

Open `pages/index.html` locally (ideally Firefox/Chrome/Edge; pages embed their physical data and work offline). Each of the six station links contains an interactive XY physical chamber map, numbered anodes, full 10 ns birth-to-anode electron drift PDF, signed muon entrance angle slices θu/θw, Garfield last-step normal approach angle and φ, and six pitch/field-wire geometries. There is also a detailed spectrum dashboard for every station, with logarithmic Y and time/angle selections.

Both **English** `*_en.html` and Russian `.html` variants are generated; the dashboard landing page and numerical CSV are in English. Existing source/station/DET shapes and physical anode numbers are used directly; no interpolation or modeled curves are substituted.

## How the science is verified

1. All 240 original chunk/station SHA256 manifests are required, with exactly 50,000 **generated parent events per Geant4 chunk**, correct J/ψ source and physical station IDs, complete group counts and saved gzip hashes.
2. On build, all source gzip files are SHA256-verified and their original per-chunk full Heed/Garfield six-LUT fingerprints checked, including the frozen 20-independent-microscopic-trials-per-box Garfield provenance.
3. All raw histogram counts and *unrounded electron moments* are summed; no re-sampling. Per-anode → chamber → station time histograms and both signed muon-angle bin sets are conserved across the entire source.
4. For each station, one merged JSON, four independent Russian/English HTML pages and SHA-pinned compact statistics are retained in the new, isolated dashboard folder. A successful rerun verifies finished outputs and skips them.
5. `verify` checks all six station JSON/HTML SHA256 digests, the 40 input manifest SHA256 values for each station, the portal and the summary CSV.

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

Publishes only the standalone HTML pages, the portal, station numerical summaries and a complete manifest to `na60drive:test_runs/mwpc_stage_c/pi_full_jpsi_2m_dashboard_v1`. Large native ROOT files, original Heed fold .gz checkpoints and merged source JSON remain on local WSL. HTML previews inside Google Drive might not execute JavaScript: download/open the HTML files locally for full interactivity.
