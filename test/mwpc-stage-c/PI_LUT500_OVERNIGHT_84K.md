# NA60+/MWPC Stage B1 — Overnight 84,000 independent avalanches on frozen LUT500

**Run design, 29 October 2026.** Nightly unattended WSL calculation.
**Scientific scope:** ONLY stages 1–2 (microscopic avalanche timing and conditional electron-response LUT). No electronics, induced strips, or claimed physical muon-angle dependence.

## Physics

Use the ORIGINAL SHA-pinned, full-coverage 500-electron-per-cell Garfield **drift** LUTs, one for each of:

| geometry | old cells | new independent AvalancheElectron starts |
| --- | ---: | ---: |
| p3_noFW | 900 | 18,000 |
| p3_FW75 | 900 | 18,000 |
| p4_noFW | 1,200 | 24,000 |
| p4_FW75 | 1,200 | 24,000 |
| **total** | **4,200** | **84,000** |

For each (geom,box), select exactly **20 different original sample indices** among old 500 by a SHA256-derived per-cell offset and a step of **137 modulo 500**, so no duplicates. Every new avalanche starts from the old CSV's **exact original electron start** at 0.1 eV, isotropic direction, and the appropriate genuine analytic 1800V geometry. NO old endpoint is reused as an avalanche start; no old collection flag is applied as a survival gate. The 500-sample source files remain read-only. Pilot is not driven by Heed or true Geant4 muon angular distributions.

There are **20 complete passes** through all 4,200 spatial cells, random-shuffled deterministically with separate seeds per pass, so that an early time-budget stop leaves full past passes and an approximately geographically balanced partial pass. 210 blocks of 400 microscopic trials each. Each block is saved to an immutable `complete.json` SHA-protected directory and validates the exact original starting coordinates, ion sectors, one-to-one trial identity, and Garfield time. An interrupted/incomplete block is NOT silently discarded or overwritten. `--resume` explicitly continues completed chunks, but requires the same C++ source/executable and frozen LUT bytes. The run does not launch uncontrolled repeated loops.

**Full ion birth count and angle census** are retained for EVERY trial, alongside gain and endpoint diagnostics. Only the POSITION reservoir is reduced from 100 to 5 random ion-birth locations per avalanche for disk size. Reservoir RNG is separate from Garfield's microscopic RNG; individual spatial samples carry weight N_birth/n_sample. This preserves unbiased pooling of sampled positions within each avalanche, but does NOT preserve every ion position. A true comparison of microphysics (full100 vs compact5 reservoir) on the SAME eight original starts is required as a preflight, checking exact equality of 39 trial CSV fields before the night run proceeds.

## Runtime budget

Measured real WSL `full432` mean 0.34664s/avalanche, 432 successful samples, no caps. 84,000 × 0.34664 seconds ≈ 8.09 **CPU hours**, excluding process startup, selection, disk and publication. Budget **8.75 wall-clock hours** (31500s), with a 5-minute finalization reserve; stop after the current 400-trial chunk on deadline, then aggregate whatever full chunks finished and upload. Actual times vary with gain, geometry, and outliers, so exact 84k completion is NOT guaranteed. The last valid checkpoint and all outputs are retained; capped high-tail avalanches are explicitly marked as right-censored and warned about in summary, not silently used as unbiased tail samples. If interrupted, an incomplete chunk is kept for forensics, not automatically overwritten.

Output root (never overwrite historical results):

`test_runs/mwpc_stage_c/pi_lut500_avalanche_overnight_v1/run_20percell_firstpass_v1/`

- `chunks/chunk_0000` … `chunk_0209`: original exact start subset CSV, source Garfield stdout/stderr, full single-seed trial CSV, compact birth spatial reservoir, timing CSV, SHA-stamped `complete.json`
- `trials_{geom}.csv` and `timing_{geom}.csv`: ALL per-trial gain and clock metrics, not summary-only
- `birth_reservoir_{geom}_part{0,1}.csv`: 5 spatial ion-birth samples per nonzero avalanche (or fewer if ion count <5), ion-count weighting; split for Drive size cap
- `cell_statistics.csv`: ALL 4,200 spatial cells with number of fresh trials, zero fraction, gain mean/standard deviation/median, timing mean and ion-weighted angular asymmetries
- `summary.json`, `manifest.json`, `run_config.json`: status, coverage, source hashes, complete/partial coverage, right-censored trials, Gaussian/zero caveats

## Commands

**Manual preflight / observe PASS, then launch overnight:**

```bash
cd ~/na6p/src/NA6PRoot-mwpc
git pull --ff-only
bash test/mwpc-stage-c/run_pi_native_muon_b1.sh build
python3 test/mwpc-stage-c/pi_run_lut500_overnight.py smoke
nice -n 15 python3 -u test/mwpc-stage-c/pi_run_lut500_overnight.py overnight --hours 8.75 --upload
```

The `build` runs offline Python regression checks (and forcibly recompiles ONLY the dedicated C++ B1 target to avoid WSL clock skew). The `smoke` verifies the genuine compact5 mode against full100 with identical physical trial CSVs for 2 original LUT500 starts per each of the four configurations. No overnight compute begins unless smoke verdict is PASS. If smoke is re-run in an unchanged, already-existing immutable smoke directory, it refuses to overwrite; proceed directly to overnight only after inspecting the saved `smoke_verified.json`.

**Single unattended command (background, preflight blocks overnight on errors):**

```bash
cd ~/na6p/src/NA6PRoot-mwpc
git pull --ff-only
nohup bash -c '
  cd "$HOME/na6p/src/NA6PRoot-mwpc" &&
  bash test/mwpc-stage-c/run_pi_native_muon_b1.sh build &&
  python3 -u test/mwpc-stage-c/pi_run_lut500_overnight.py smoke &&
  nice -n 15 python3 -u test/mwpc-stage-c/pi_run_lut500_overnight.py overnight --hours 8.75 --upload
' > "$HOME/na6p/lut500_overnight_20261029.log" 2>&1 < /dev/null &
echo "Night PID $!"
```

**Check morning progress:**

```bash
tail -n 35 "$HOME/na6p/lut500_overnight_20261029.log"
```

**Power warning:** WSL jobs are suspended if WINDOWS itself sleeps; plugging the PC into power and configuring Windows to NOT sleep for the night is necessary. Screen may lock/off; the machine must remain awake.

**Manual restart after unexpected termination:** inspect incomplete block and contact the assistant instead of deleting or overwriting files. If all existing chunks are complete, invoke `overnight --resume --hours 8.75` to continue (and append `--upload` if wanted). An incomplete saved block requires explicit forensic recovery before resumable script can move past it.

## Google Drive

`--upload` invokes EXISTING `pi_publish_run_to_drive.py` with the pre-configured and upload-tested `na60drive:` remote; creates a uniquely named immutable run folder below `avalanche_3p3_runs`. It hashes local data before transfer, uploads with `rclone --immutable --checksum`, and independently checks the transferred files. The local run can have hundreds of raw checkpoint folders; the Drive publication contains the full per-geometry per-trial CSV, sampled ion-birth positions, cell-level summaries, config and SHA manifest, but not nested raw checkpoint logs.

Google Drive folder: https://drive.google.com/drive/folders/1Xbh45V86BhVu-WLYSM5qtosqPEG363M7

This C++ overnight mode has NOT been compiled or physically run on the user's WSL at document creation. Do not claim it was.


## 29 Oct 2026: smoke failed in Python-only CSV identity check; RECOVER, DO NOT DELETE

User confirmed successful `build` and full 84,000 frozen source-position extraction, but compact5 smoke raised:

`ValueError: LUT identities and endpoint flags no longer paired`

Root cause is isolated in `pi_run_avalanche_lut_sidecar.py::validate_outputs`: overnight Python generated original `random_seed` as an `int`, while the same integer from C++ output is necessarily a `str` loaded with `csv.DictReader`. Previous check compared `"123" != 123` and falsely rejected the row. Fixed to `int(...) != int(...)` and added a regression test for both equal-different-types and unequal numerical seeds.

The C++ process itself had already returned success; the error happened only AFTER it completed, in Python output validation. The new explicit `recover-smoke` command **does not repeat or modify any original eight Garfield trial CSVs**. It requires the exact original selected source keys, all raw CSV output identities and full ion birth count conservation, and C++ PASS log, then creates a SHA-stamped `complete.json`. Next it proceeds with the originally intended full100-vs-compact5 SAME-input Garfield physical equivalence comparison; only upon success writes `smoke_verified.json`. The original incomplete smoke directory is deliberately preserved, without `rm`, `mv`, `sync`, `--force`, or manual CSV edits.

**Continue from the observed error**:

```bash
cd ~/na6p/src/NA6PRoot-mwpc
git pull --ff-only
python3 test/mwpc-stage-c/test_pi_avalanche_lut_sidecar.py
python3 test/mwpc-stage-c/test_pi_lut500_overnight.py
python3 test/mwpc-stage-c/pi_run_lut500_overnight.py recover-smoke
```

**Require** `PASS recovered 8 ORIGINAL ... without rerun` and `PASS compact5 smoke AND exact 100-vs-5 ... physics identity`. If anything fails, STOP and preserve diagnostic files; no night run.

Only then:

```bash
nohup nice -n 15 python3 -u test/mwpc-stage-c/pi_run_lut500_overnight.py overnight --hours 8.75 --upload > "$HOME/na6p/lut500_overnight_20261029.log" 2>&1 < /dev/null &
echo "Night PID: $!"
```

This Python-only recovery patch DOES NOT change C++ source/binary; recompilation is NOT needed after the previous successful B1 build. New overnight physics on 84k starts remains unrun/unverified until user logs and Drive publication arrive.
