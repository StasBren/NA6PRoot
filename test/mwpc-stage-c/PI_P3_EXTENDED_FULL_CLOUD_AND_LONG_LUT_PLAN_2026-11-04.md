# NA60+ MWPC p3 — status, extended full-ion cloud study, then long conditional-LUT campaign

**Protocol date:** 2026-11-04. This is a new independent archive, NOT an update
to the first sealed 48-full-cloud physics pilot. Only `p3_noFW` and
`p3_FW75`; decisions, errors, confidence intervals and calibration
are reported **separately** for both geometries.

## Exactly where we are

1. **Source**: actual Geant4 muon -> TrackHeed primary electron positions,
   and full Garfield++ `AvalancheMicroscopic` are working.
2. **Count+angle LUT**: frozen A20 18,000 independent microscopic avalanches
   per 3-mm geometry, and targeted A100 4,000 extra per 3-mm geometry.
   Fixed guarded-neighbor A20 passed on older independent muons. The
   gamma/Pólya count-only benchmark is close or slightly better on mean
   scores, but `gamma-guarded` CRPS differences on the 120 new parent
   muons have cell/parent-bootstrap CIs including zero. No significant
   model replacement decision. `DriftLineRKF + Pólya` adds essentially
   no proven local count benefit for these selected domains.
3. **Full cloud (48)**: two 24-start per-geometry source sets,
   99,388 *exact ionization birth callbacks* from 48 independent
   full microscopic avalanches, verified against exact C++ moments,
   angular counts and complete SHA-sealed point records.
   Each ion-birth point is (u,v,w,t) of a positive ion being CREATED,
   **not** a drifting ion position or induced cathode charge.
   Previous reservoir sizes 128 and 512 show statistically meaningful
   improvement at 512 on 12-angle x 6-time TV, mean radius and q90 birth
   time, but only 18/19 independent nonzero avalanches exist.
4. **Next, THIS program**: amplify exact-cloud COMPRESSOR validation to
   `480 p3_noFW + 480 p3_FW75` independent full microscopic avalanches,
   with source region/sign stratification, multiple independent source
   positions in very near-wire D500 cells. Compare 32/128/256/512/**1024**
   versus each same exact microscopic full cloud. Do not call any
   induced-cathode model.
5. **Later actual long-run**: only after separately selecting a sufficiently
   accurate cloud representation for the two p3 geometries, train
   conditional spatial-time micro-avalanche generators on a much
   larger SHA-frozen matched grid. Proposed first learning curve is
   20/50/100 independent avalanches per 0.1 x 0.2 mm D500 cell,
   emphasizing near anode, near field wire, narrow transport
   boundaries and nonzero/tail events. That would be roughly
   36k/90k/180k total microscopic trials for BOTH p3 geometries at
   20/50/100 samples per 900 cells/geometry, before adaptive additions.
   These are SCALE OPTIONS, not yet approved production-job counts.
   Reserve spatially separate cells and real independent Geant4/TrackHeed
   parents as holdouts. Validate *joint* distribution of number of ions,
   zero fraction, angular pattern, radius, w coordinate, birth-time
   and their correlations. Then, and ONLY then, discuss ion drift and
   Shockley–Ramo cathode signals.

## Why this extended job is NOT the long production LUT

- This run tests how many `(u,v,w,t)` birth points preserve the shape
  of an ALREADY KNOWN avalanche. It **does not train** an emulator
  that generates new avalanche clouds from arbitrary Heed electron starts.
- Two geometries use separate source selection/analysis. The D500
  source is deliberately stratified, NOT representative of the
  muon incidence PDF. `mean Nion` in the source is not the actual
  muon mean.
- Full-cloud *compression* comparisons are paired to the same exact
  Garfield event; `64` or `32` random reservoirs are not new
  independent avalanches.
- Near-wire cells are few. Many independent initial electrons from
  the same spatial cell can be used, but generalization uncertainty
  must be bootstrapped by **spatial cell**, not by ion count.
- A 512-point statistical reservoir does not mean only 512 ions
  were physically born. Exact `N_ion` is stored separately, and
  each selected point represents `N_ion / min(N_ion,512)` births.

## Independent source and physics controls

A dedicated script `pi_p3_extended_full_cloud_campaign_v1.py` reads
**original frozen D500 initial electron records**, only passes 100–139
of the reproducible 500-start permutation, never the A20/A100 sample
passes 0–99. It also SHA-checks and excludes the 48 earlier full-pilot
start IDs. The **same exact** `phaseB_native_muon_joint_cloud_lut_v1.cxx`
compiled executable and unaltered Garfield physics is used. No
changes to C++ or any old SHA manifest.

The original 900 cells/geometry are classified by actual electron
position: `near_anode`, `transition`, `bulk`, and, for
`p3_FW75` only, `near_field`, separately for `v>0` or `v<0`.
480 initial-electron starts/geometry means 80 per noFW stratum
(6 strata) and 60 per FW75 stratum (8 strata). Across each stratum
we first spread selections over as many distinct spatial cells as
available, then sample independent new starts within cells.

The source is sealed before physics. Runs are in separate
`test_runs/mwpc_stage_c/pi_p3_extended_full_cloud_v1/`. Each
8-electron chunk retains the FULL native Garfield ion birth CSV,
exact birth moments, all trial outcomes, Garfield timing and a
SHA manifest. A run refusing a partial SHA checkpoint cannot silently
overwrite it. All chunks use one untouched C++ executable hash from
preflight. No network jobs, no old-archive edits.

## Commands (from project WSL)

### A. Pull and test software, NO physics yet

```bash
cd ~/na6p/src/NA6PRoot-mwpc
git switch feature/mwpc-stage-c
git pull --ff-only
PY=test_runs/mwpc_stage_c/.venv_pi_slides/bin/python

"$PY" -m py_compile \
  test/mwpc-stage-c/pi_p3_extended_full_cloud_campaign_v1.py \
  test/mwpc-stage-c/pi_p3_extended_cloud_readout_v1.py
"$PY" test/mwpc-stage-c/pi_p3_extended_full_cloud_campaign_v1.py selftest
"$PY" test/mwpc-stage-c/pi_p3_extended_cloud_readout_v1.py selftest
"$PY" test/mwpc-stage-c/test_pi_p3_extended_full_cloud_campaign_v1.py
```

The source regression builds two 900-cell x 40-pass synthetic datasets
in temporary storage, checks both 480 and 960 selections and
preparation/manifest SHA. It does **not** call Garfield++.

### B. Freeze source and validate existing native parser, then ONE 8-start block

```bash
"$PY" -u test/mwpc-stage-c/pi_p3_extended_full_cloud_campaign_v1.py prepare --per-geometry 480
"$PY" -u test/mwpc-stage-c/pi_p3_extended_full_cloud_campaign_v1.py preflight
"$PY" -u test/mwpc-stage-c/pi_p3_extended_full_cloud_campaign_v1.py run --max-chunks 1 --hours 1
"$PY" test/mwpc-stage-c/pi_p3_extended_full_cloud_campaign_v1.py status
```

The previous C++ binary is used as-is. **DO NOT rebuild** after the
preflight, as its hash is frozen in source and checkpoint metadata.
The first block should be `VERIFIED CHUNKS 1 / 120`. If not, stop:
inspect raw `physics.log` and never overwrite any partial/failed
chunk or interpret it as zero-gain.

### C. RUN remaining 119 chunks (after confirmed first SHA checkpoint)

```bash
nohup nice -n 15 "$PY" -u \
  test/mwpc-stage-c/pi_p3_extended_full_cloud_campaign_v1.py \
  run --max-chunks 119 --hours 18 \
  > "$HOME/na6p/p3_extended_full_cloud_v1.log" 2>&1 < /dev/null &
echo "Background WSL PID: $!"
```

Read-only status:

```bash
tail -n 35 "$HOME/na6p/p3_extended_full_cloud_v1.log"
"$PY" test/mwpc-stage-c/pi_p3_extended_full_cloud_campaign_v1.py status
```

The native runner checks at each complete 8-trial boundary for
`--hours 18`; if the budget expires at a boundary, repeat the same
`run` command without changing its source selection or binary.
Run **only one** writer process at any time. Stopping WSL/suspending
the machine *inside* a chunk can leave an intentionally unsealed
chunk, in which case the program stops rather than silently reusing
or deleting it. There is no automated cleanup/recovery.

This 960-avalanche COMPLETE CLOUD job is a calibration-size run, likely
much smaller than the future long training run. Native runtime and
disk growth should be inferred from first completed blocks, not
promised in advance.

### D. After 120/120 chunks SHA PASS: 5-cap analysis, then box-bootstrap readout

```bash
"$PY" -u test/mwpc-stage-c/pi_p3_extended_full_cloud_campaign_v1.py \
  analyze --reps 32
"$PY" -u test/mwpc-stage-c/pi_p3_extended_cloud_readout_v1.py \
  report --bootstrap 3000 \
  > "$HOME/na6p/p3_extended_cloud_report_v1.log"
tail -n 100 "$HOME/na6p/p3_extended_cloud_report_v1.log"
```

Analysis output stays separately under
`pi_p3_extended_full_cloud_v1/analysis_complete_v2/`.
The readout SHA-checks every native physics chunk and all analysis
outputs. It includes each geometry and near-wire/bulk stratum,
five reservoir caps, exact/uncompressed cases, bootstrap by
independent **spatial cells**, ion-count-weighted sensitivity, and
paired 128→512 and 512→1024 comparisons.

## Decision gate before large production run

For each geometry independently review (i) coarse joint angle/time TV,
(ii) ion birth radial-moment and upper time-quantile errors, (iii)
the undefined/rare one-sided ion asymmetry cases, (iv) compressed-only
and larger-gain cases, and (v) CIs clustered by spatial cell. If 512
still fails our accepted precision, choose 1024 or a better conditional
representation BEFORE creating tens of thousands of clouds. Existing
C++ supports full birth mode and 512 birth mode; a production 1024 mode
would be added as its own isolated target only after gate. Keep
(N_{ion}), source coordinates and birth times together **per
microscopic trial** and preserve covariance between count and cloud
shape. Count-only Gamma/Pólya may be a separate fast proposal but
cannot generate full ion clouds alone.

**No conclusion that 512 or 1024 is cathode-signal-ready is justified
by these compression-only tests.**

## 2026-11-04: LIVE status race on chunk_0013 (13 sealed / 120)

Observed after background `nohup` (PID originally **142387**):
`chunk_0000` ... `chunk_0012` each printed `PASS SHA FULL BIRTHS`;
the writer next printed `FULL MICROSCOPIC ... chunk 13`.
A simultaneous `status` process raised

```text
ValueError: Unsealed/incomplete chunk must be inspected:
.../full_cloud_chunks/chunk_0013
```

This **does not show Garfield failed**. The original `status` called
`core.completed()`, a deliberately strict all-chunks SHA verifier
also used by final analysis. The asynchronous writer creates
`chunk_0013` and its `physics.log` first, then later creates
`complete.json`. The status query raced in that normal interval.

**Patch:** the NEW extended campaign's `status()` now uses a
read-only checkpoint directory scan:
- audits every sequential **sealed** chunk using the original
  complete.json, source data and all physical file SHA;
- reports the single newest **unsealed** chunk as
  `ACTIVE or INTERRUPTED; not counted`, without treating the
  presence of a directory as proof of successful simulation;
- still fails closed on gaps, two unfinished chunks, a sealed
  later chunk after an unsealed predecessor, or a changed sealed SHA;
- **does not alter** `core.completed()`, `core.run()`,
  C++ target, preflight binary hashes, original 48-cloud pilot or
  any physical output. Synthetic tests now include the live-reader race.

Run (safe while the original one writer is still running):

```bash
cd ~/na6p/src/NA6PRoot-mwpc
git switch feature/mwpc-stage-c
git pull --ff-only
PY=test_runs/mwpc_stage_c/.venv_pi_slides/bin/python
"$PY" test/mwpc-stage-c/pi_p3_extended_full_cloud_campaign_v1.py selftest
ps -p 142387 -o pid,ppid,stat,etime,cmd || true
tail -n 25 "$HOME/na6p/p3_extended_full_cloud_v1.log"
"$PY" test/mwpc-stage-c/pi_p3_extended_full_cloud_campaign_v1.py status
```

**Do not launch another `run` process while PID 142387 or its
Garfield child is still active.** If the PID has exited and an
unsealed chunk remains, do not retry blindly: inspect
`chunk_0013/physics.log`, checkpoint files and process exit
reason. Keep them for exact forensic recovery. Only after all
120 chunks are `VERIFIED` should `analyze --reps 32` run.

## 2026-11-04: All 960 true Garfield avalanches PASS; selftest fixture correction

WSL ended successfully with **120/120 SHA-sealed blocks**, **960/960
independent D500 initial electron avalanches**, **480 per each of the two
p3 geometries**, and **2,981,533 exactly recorded ionization-birth
points**. The background job ended normally with
`PASS ALL COMPLETE 960 independent exact-start p3 avalanches`.
The updated read-only `status` independently verified all 120 native
chunk SHA receipts. No missing or unfinished physical chunks exist.

One Python `selftest` failed after the status-reader code update:
the synthetic test created `chunk_0003` unsealed, then created
`chunk_0002` also unsealed, but incorrectly EXPECTED the scanner
to accept the two simultaneous unsealed physical chunks. The scanner
correctly rejected this as an invalid campaign state. In the
`selftest`, this is now an EXPECTED failure case; when
`chunk_0002` is sealed, the remaining `chunk_0003` is correctly
reported as the single possible in-flight chunk. A final all-sealed
case was added as well.

**Only Python selftest was repaired.** No changes to the executable,
physical source, frozen input, physics records, SHA receipts, strict
chunk scanner or first 48 full-cloud pilot.

Next WSL commands:

```bash
cd ~/na6p/src/NA6PRoot-mwpc
git switch feature/mwpc-stage-c
git pull --ff-only
PY=test_runs/mwpc_stage_c/.venv_pi_slides/bin/python
"$PY" test/mwpc-stage-c/pi_p3_extended_full_cloud_campaign_v1.py selftest
"$PY" test/mwpc-stage-c/pi_p3_extended_full_cloud_campaign_v1.py status
# Analysis reads only the already-completed exact birth points:
"$PY" -u test/mwpc-stage-c/pi_p3_extended_full_cloud_campaign_v1.py analyze --reps 32
"$PY" -u test/mwpc-stage-c/pi_p3_extended_cloud_readout_v1.py report --bootstrap 3000 \
  | tee "$HOME/na6p/p3_extended_cloud_report_v1.log"
```

If `analysis_complete_v2/` already exists from an earlier attempt,
DO NOT delete it. Inspect the directory and its log first; the analyzer
refuses overwrite. The analysis and report must be completed before
fixing any final 512/1024 compression format or starting the much
larger conditional full-LUT learning-curve experiment.
