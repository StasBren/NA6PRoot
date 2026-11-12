# Additional NA60+ Penning057 30k summaries

Source: `pi_penning057_extract_additional_summaries_30k_v1.py`.
Run on local WSL with existing read-only SSHFS mount. No Garfield re-simulation.

```bash
cd ~/na6p/src/NA6PRoot-mwpc
git pull --ff-only
mountpoint ~/na6p/cluster_mwpc
PY=python3
SCRIPT=test/mwpc-stage-c/pi_penning057_extract_additional_summaries_30k_v1.py
OUT=$PWD/test_runs/mwpc_stage_c/extra_30k_summaries_v1
mkdir -p "$OUT"
$PY -m py_compile "$SCRIPT"

# Pilot: first 20 real Jpsi muons:
$PY -u "$SCRIPT" extract --root ~/na6p/cluster_mwpc --out "$OUT" --channel Jpsi --max-chunks 1

# Full archive: 500 20-muon chunks per channel; resumes already completed checkpoints.
set -o pipefail
$PY -u "$SCRIPT" extract --root ~/na6p/cluster_mwpc --out "$OUT" --channel All --report-every 1 2>&1 | tee -a "$OUT/extract.log"
$PY "$SCRIPT" status --out "$OUT"
$PY -u "$SCRIPT" combine --out "$OUT"
ls -lh "$OUT"/summary_csv_gz/*
```

Saved tables (five gzip CSVs and provenance JSON):

- `avalanche_longitudinal.csv.gz`: per original channel/crossing/Heed-primary exact start `w0_mm`, exact avalanche ion/gain counts, and sample-derived longitudinal mean, population std, 5th/95th percentile of ion BIRTH w. A maximum of 128 uniform birth samples exists per avalanche; estimates are NOT exact for larger avalanches.
- `crossings_3d.csv.gz`: real 3D Geant4 entrance and exit, momentum, path length and deposited energy, with matched TrackHeed cluster/electron counts.
- `heed_clusters.csv.gz`: all TrackHeed clusters including ZERO accepted primaries, original/accepted multiplicities, 3D centre positions and times.
- `avalanche_wall_times.csv.gz`: original microscopic avalanche steady-clock elapsed wall seconds, linked to exact primary and channel. Not process CPU seconds.
- `muon_wall_times.csv.gz`: original Heed, Garfield+output and complete-muon elapsed wall seconds. Not process CPU seconds.

Each verified 20-muon chunk produces one small local ZIP checkpoint. Interrupted extraction resumes without restarting processed chunks. All original eight SHA-256 seals per chunk are verified while streaming the native CSV once. Source archive remains read-only.

Not available from existing 30k native campaign: matched electron arrival time, longitudinal position, azimuth and energy on a pre-multiplication surface. Those require a new small calibration with matched trajectories. The archived assigned starting energy of 0.1 eV must not be described as a measured arrival energy. Defensible speedup also requires a matched fast-model benchmark.
