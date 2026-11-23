# Full archived 1800 V MWPC per-anode gain/asymmetry analysis

## Scope

This is the full available frozen original WIPP 1800 V run, approximately 8,000
Geant4 muons, with the exact event count established and locked at preparation.
Every complete original 20-muon chunk with a valid verified_report.json is
included. New production chunks completed after preparation are not silently
added. No Garfield++ simulation is launched.

The historical numerical geometry: P3 noFW 9-wire analytic model, 3mm
wire pitch, 30um diameter, symmetric cathode ±3mm, 1800V, Ar/CO2 70:30,
Penning 0.57. This is not the final as-built MWPC chamber configuration.

## WIPP operations

Run on wipp-an1 in the separate, previously cloned worktree.

```bash
cd "$HOME/NA6PRoot-wiretruth1800-analysis"
git fetch origin feature/p3-1800v-wiretruth-8k-analysis-20261123
git merge --ff-only FETCH_HEAD

FULL=test/mwpc-stage-c/wiretruth1800/full1800.sh
bash -n "$FULL"

# Necessary because C++ now accepts --no-plots and the full ROOT merger is new
bash "$FULL" build

# Create this once and retain it throughout analysis/merge/retrieval
export MWPC_1800_FULL="$PWD/test_runs/wiretruth1800_FULL_$(date -u +%Y%m%dT%H%M%SZ)"
export MWPC_1800_PBS_QUEUE=N
export MWPC_1800_PARALLEL=2
export MWPC_1800_SHARD_CHUNKS=6

# Freeze every original source chunk that was complete at this instant
bash "$FULL" prepare

# Start analysis of EVERY selected original muon, not an intermediate pilot
bash "$FULL" submit
```

Preparation prints exact n_muons, chunks and shards. For roughly 400
completed original chunks, expect ~67 independently restartable array tasks
because shard boundaries are kept within channels. Each full shard contains
up to six original 20-muon chunks = 120 independent Geant4 muons. The PBS
array uses N routing queue, 1 CPU/4GB/8h/io=0.1 for each task and throttles
to two active tasks. Exact duration is unknown; do not assume the pilot's
35-minute walltime scales linearly.

## Internal integrity/physics workflow

Every shard checks frozen code/Garfield binary and source selection locks,
rechecks its own immutable source manifest, reads original completed
verified_report.json and compares its SHA-256 to the prepared manifest,
copies the original Geant4 CSV and four Garfield CSV files to PBS node-local
TMPDIR while calculating and verifying every recorded SHA-256, then runs
the C++/ROOT physical-anode analyzer on those local immutable copies.

Each worker creates ROOT TTrees (muons, wires, seeds, seed_wires), ALL
original gain/ion/asymmetry histograms, independent-muon and per-wire QA,
and CSV summaries. Flag --no-plots suppresses *per-shard PNG generation*,
which previously caused many writes through slow shared storage. No physics
information is omitted; final PNGs are generated once from merged ROOT.

Every shard uses a unique output folder, writes SHA256SUMS and provenance,
checks its copyback, then seals results by atomic rename. Original data
is never modified. Repeated successfully sealed shards are verified
and skipped.

## Monitoring, retry and final merge

```bash
bash "$FULL" status
/opt/pbs/bin/qstat -t -u "$USER"
```

Status prints the number of completed original muons and missing shard IDs.
If all PBS array members have finished and some shards failed, examine
the corresponding PBS logs and any stale attempts/claim_NNNN owner first.
Once there is no running original shard still owning that claim, rerun:

```bash
bash "$FULL" resume
```

This requeues the same frozen range, but already completed SHA-sealed
shards are skipped. Do not change scientific code, commit HEAD, source
selection or compiled binaries after the prepare lock. If code changes
are required, start a fresh campaign directory.

Once status shows ALL original muons processed:

```bash
bash "$FULL" merge-submit
```

The separate PBS merge job checks every source result SHA, copies all
ROOT+summary files to node-local scratch, combines ROOT TTrees and
TH1/TH2/TProfile histograms, validates no duplicate original
crossing_id and counts, sums original per-shard integers, and
recomputes aggregate Q2>5% probabilities, binomial uncertainties,
and mean Q2 fraction over positive-charge muons. No averaging shard
percentage estimates.

The merged output is $MWPC_1800_FULL/merged with:
- wiretruth1800_full.root: all original crossing, wire, primary, primary-wire
  TTrees and all channel/station histograms
- channel_station_summary.csv: exact additive physical counters
- asymmetry_denominators.csv: ion birth A defined, zero Q vs positive Q
- merged_root_qa.txt: unique crossing IDs, channel/primary/wire census
- full_gain_vs_A_hemi_wire_all.png and all channel variants
- full_gain_vs_A_hemi_primary_wire_all.png
- full_gain_vs_A_TB_wire_all.png
- full_A_hemi_vs_local_wire_all.png
- full_A_hemi_vs_signed_impact_all.png
- full_A_hemi_vs_theta_u_all.png
- full_A_hemi_vs_theta_w_all.png
- full_muon_sharing_vs_A_hemi_all.png
- full_collected_e_per_wire_all.png
- 18 channel×station gain-A and asymmetry-wire plots
- SHA256SUMS, frozen_campaign.json and sealed.json.

After the merge PBS job exits 0:

```bash
MERGED="$MWPC_1800_FULL/merged"
cat "$MERGED/merged_root_qa.txt"
cat "$MERGED/channel_station_summary.csv"
cat "$MERGED/asymmetry_denominators.csv"
(cd "$MERGED" && sha256sum -c SHA256SUMS)

cd "$MERGED"
tar -czf "$HOME/wiretruth1800_FULL_plots_and_stats.tar.gz" \
  *.csv merged_root_qa.txt frozen_campaign.json full_*.png sealed.json
```

On local WSL (not WIPP), transfer to Windows Downloads:

```bash
scp stanis@wipp-an1.weizmann.ac.il:~/wiretruth1800_FULL_plots_and_stats.tar.gz \
  /mnt/c/Users/stanm/Downloads/
```

Upload the bundle to ChatGPT. The combined ROOT file may be transferred
separately for rigorous per-muon clustered correlation fits.

## Scientific interpretation limits

The plot axes use G_muon,wire = actual Garfield Q_collected,wire /
N_Heed,muon; an individual Heed primary-wire gain is its exact Q
from one initial Heed electron. Ion birth asymmetries are:
A_hemi = (N_up-N_down)/(N_up+N_down) and
A_TB = (N_top-N_bottom)/(N_top+N_bottom).
When the denominator is zero, asymmetry is undefined, never zero.
A wire can have positive ions but zero electron collection; its log
gain is undefined, and it is retained in ROOT and QA tables.
Per-wire observations are correlated within the same original muon,
so regression uncertainty must cluster by parent Geant4 crossing.
Native ion-birth-nearest-anode assignment differs from electron
endpoint anode and from actual U/V cathode strip response.

No new 2100V analysis or geometry scan is implied by the 1800V results.
For future variations (strip pitch/width, geometry, HV), use a
validated configuration-locked weighting/drift-field and avalanche
model rather than reinterpreting this source as actual hardware.
