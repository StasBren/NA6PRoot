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

## 3. Google Drive publication — EXISTING rclone configuration

**Previously configured and tested 24 October 2026.**
The user already has a working WSL rclone OAuth remote `na60drive:`
whose **root folder is NA60plus_MWPC_Research**, ID
`1Svbp3E230YmVk53D0M1AHwTGo6cBxy38`.
Its custom OAuth client was validated, and quick-results reports were
successfully uploaded. Do **NOT** rerun `rclone config`, create a new
`gdrive` remote, or upload the entire ~252 GiB dataset again.

Existing repository helpers:
- `run_pi_drive_archive.sh` (large archive; DO NOT use for this study).
- `run_pi_drive_results_only.sh` (selected result folders).
- `run_pi_drive_quick_reports.sh` (small reports, previously tested).

New destination folder was created under that exact project root:
[avalanche_3p3_runs](https://drive.google.com/drive/folders/1Xbh45V86BhVu-WLYSM5qtosqPEG363M7).

Check connection, **no upload**:

```bash
rclone lsf na60drive: --max-depth 1
rclone lsf na60drive:avalanche_3p3_runs/ --max-depth 1
```

### Dry run, existing `na60drive` remote

```bash
python3 test/mwpc-stage-c/pi_publish_run_to_drive.py \
  --source test_runs/mwpc_stage_c/avalanche_3p3_gain_v1/smoke5 \
  --run-id avalanche_3p3_smoke5_20261028
```

### Explicit upload after reviewing dry-run

```bash
python3 test/mwpc-stage-c/pi_publish_run_to_drive.py \
  --source test_runs/mwpc_stage_c/avalanche_3p3_gain_v1/smoke5 \
  --run-id avalanche_3p3_smoke5_20261028 \
  --upload
```

The default is now `--remote na60drive`; an explicit `--remote` override
is only needed if the user's local remote name actually changes.
The script writes to
`na60drive:avalanche_3p3_runs/avalanche_3p3_smoke5_20261028/`
relative to the **already configured project root**, without overriding
the remote's `root_folder_id`.

`rclone copy --immutable` never deletes remote files or modifies differing
existing files; a post-upload `rclone check --one-way` confirms readback.
Transfers are intentionally throttled (`--transfers 1 --tpslimit 2`).
The immutable local publication bundle contains only small summaries,
manifests, plots and bounded CSV/logs; raw `crossings.csv.gz`, large ZIP
and HTML are excluded by design. The same unchanged bundle can be reused
from dry-run to actual upload.

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
