# 28 Oct 2026 — executable next step: new ±3 mm avalanche gain / safe Google Drive publication

## Fixed scope

Four analytic-field geometries, cathodes always at **v = −3 mm / +3 mm**:

| ID | Anode pitch | Grounded 75 μm field wire | Role |
| --- | ---: | --- | --- |
| **p3_noFW** | **3 mm** | No | **Primary model** |
| p3_FW75 | 3 mm | Yes | FW effect at 3 mm |
| p4_noFW | 4 mm | No | Pitch effect |
| p4_FW75 | 4 mm | Yes | Pitch × FW effect |

Other fixed inputs: 30 μm anodes at +1800 V, Ar/CO₂ (70:30), 293.15 K, 760 Torr, B=0 T. New Phase-A gain smoke uses the same explicit **single-electron initial seed (u,v)=(0, +1 mm)** and same per-event **random seed integers** across configurations. This is **not** a representative muon/Heed ionization ensemble or calibrated experimental gain.

Old `test/mwpc-avalanche` Phase A/B outputs made for **2+4 mm** gaps are methodological baselines only, not 3+3 physics numbers. They are NEVER overwritten.

## 1. Changes made to existing Garfield engine

Existing `test/mwpc-avalanche/phaseA_gain_fluctuations.cxx` now supports:
- `--field-wire-diam-um 0|75` with grounded midpoint wires `u=(i+0.5)*pitch`, v=0 and explicit geometry-overlap guard;
- `--base-seed N` for a reproducible Garfield RNG sequence per seed event; omitted = legacy random engine behavior;
- `--avalanche-limit N` with old 200,000 default preserved;
- original no-FW 4 mm **2+4 mm** defaults still intact for historical tests.

For each incoming electron seed the existing output columns are preserved:
`event,avalanche_electrons,avalanche_ions,collected_on_wire,attached,other_endpoints,mean_collection_time_ns`.
The appropriate effective single-seed collected gain is `G_eff=collected_on_wire`,
**not** `avalanche_electrons+avalanche_ions`. The code flags possible avalanche-size-limit censoring. Its wire endpoint check counts electrons collected on anodes, not grounded field wires.

## 2. Build and run on the local WSL machine

```bash
cd ~/na6p/src/NA6PRoot-mwpc
git switch feature/mwpc-stage-c
git pull --ff-only

# Rebuild the ONE existing Garfield executable with the new field wires.
bash test/mwpc-stage-c/run_pi_avalanche_3p3_gain.sh build

# Local synthetic tests, source/binary checks, NO physics.
bash test/mwpc-stage-c/run_pi_avalanche_3p3_gain.sh check

# FIRST PHYSICS: 5 single-electron avalanches per geometry, only 20 total.
bash test/mwpc-stage-c/run_pi_avalanche_3p3_gain.sh smoke5

bash test/mwpc-stage-c/run_pi_avalanche_3p3_gain.sh status

# Only after reviewing 3+3 geometry, possible size caps, and 5-seed results:
bash test/mwpc-stage-c/run_pi_avalanche_3p3_gain.sh pilot100
```

Separate outputs are created under:

`test_runs/mwpc_stage_c/avalanche_3p3_gain_v1/smoke5/`
and
`test_runs/mwpc_stage_c/avalanche_3p3_gain_v1/pilot100/`.

Each output stores four independent CSVs/logs, `run_config.json`,
`summary.json`, `manifest.json` and input source/binary SHA256,
run code commit and collection/zero/cap statistics.
An already existing output directory is NEVER overwritten; inspect
an interrupted run instead of automatically restarting it.

For scientific interpretation, if any cap-like events occur,
do not treat the mean as an uncensored gas gain. The location (0,+1 mm)
is deliberately fixed for this initial controlled comparison; it
does not imply the full range of Heed electron starting positions.

## 3. One-time Google Drive upload setup on WSL

The ChatGPT-connected Google Drive account **does not transfer OAuth
credentials to the local WSL shell**. WSL must have its own `rclone`
installation and explicitly authorized Google Drive remote.

```bash
sudo apt-get update
sudo apt-get install rclone
rclone config
rclone listremotes
```

During `rclone config`, create a Google Drive remote, e.g. `gdrive`,
using the user-approved OAuth flow. Choose appropriate Drive access
(not appfolder-only / drive.file-only for pre-existing folders).
Never commit OAuth tokens or the rclone configuration to GitHub.

The destination folder **already exists and was verified through the
connected Google Drive**:
[NA60plus_MWPC_Research / avalanche_3p3_runs](https://drive.google.com/drive/folders/1Xbh45V86BhVu-WLYSM5qtosqPEG363M7).

The publisher explicitly sets the remote root to this folder by ID,
so it does not depend on navigating the Shared Drive root correctly.

### Dry run (no Google Drive write)

```bash
python3 test/mwpc-stage-c/pi_publish_run_to_drive.py \
  --source test_runs/mwpc_stage_c/avalanche_3p3_gain_v1/smoke5 \
  --run-id avalanche_3p3_smoke5_20261028 \
  --remote gdrive
```

### Explicit upload after verifying the dry-run plan

```bash
python3 test/mwpc-stage-c/pi_publish_run_to_drive.py \
  --source test_runs/mwpc_stage_c/avalanche_3p3_gain_v1/smoke5 \
  --run-id avalanche_3p3_smoke5_20261028 \
  --remote gdrive --upload
```

`rclone copy --immutable` never deletes remote files or modifies differing
existing files; a post-upload `rclone check --one-way` confirms readback.
It copies an **immutable local publication bundle** of only small
summaries, manifests, plots and bounded CSV/logs; raw `crossings.csv.gz`,
large ZIP/HTML, etc are intentionally excluded. It stages under
`test_runs/mwpc_stage_c/drive_publications_v1/`; the same unchanged
bundle can be reused from dry-run to final upload, but changed sources
and colliding run IDs are refused.

The newly created Drive folder can then be read by ChatGPT directly,
without attaching a ZIP to the conversation.

## 4. Next scientific steps after the gain smoke

1. Validate the four gains with matching seed populations; study
   cap rate, zero-gain fraction, mean and RMS. Do not draw conclusions
   from five seeds; run a larger sample after the first audit.
2. Adapt the existing **Stage-B** ionisation-birth angle scan to four
   symmetric 3+3 configurations; retain **both** definitions:
   `A_hemi` (all upper/lower half-circumferences) and
   `A_TB` (only ±45° cathode-facing sectors).
3. Record the full **zero-response population**, including failed
   avalanches (old Stage-B scan silently filters some failed seeds).
   Ion-birth A and avalanche endpoint A are separate observables.
4. Compute actual time-dependent cathode charge `A_Q(T)` with
   ion drift and weighting fields for both ±3 mm cathodes;
   it need not follow `A_TB` monotonically.
5. Only then fold microscopic gain and induced charge with selected
   real Geant4+Heed muon gas crossings, maintaining proper parent/crossing
   identifiers and uncertainties.

**No new Garfield physics has been run via GitHub or through ChatGPT.**
Actual compilation, numerical output and distribution verification will
occur on the user's WSL machine after the above commands.
