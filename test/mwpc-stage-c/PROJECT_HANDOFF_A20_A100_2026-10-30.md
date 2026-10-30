# Project continuity memo, 2026-10-30 — NA60+ MWPC

Read this memo after editing/forking an earlier ChatGPT message. Code and raw physics files may survive a conversation fork even if later chat messages do not. This memo records status, not new physical execution.

## Sources and definitions

- WSL: ~/na6p/src/NA6PRoot-mwpc
- GitHub repository: StasBren/NA6PRoot
- Branch: feature/mwpc-stage-c
- DriftLUT20 (D20) = historical **20 electron DRIFT trials** per spatial cell, NOT 20 microscopic avalanche calculations.
- DriftLUT500 (D500) = historical **500 electron DRIFT trials** per cell, NOT 500 microscopic avalanche calculations.
- AvalancheLUT20-from-Drift500 (A20@D500) = 20 genuinely independent microscopic full Garfield++ avalanches per 0.1x0.2mm cell, started at selected exact D500 electron START locations. 4,200 cells x20 = **84,000 original full avalanche calculations**.
- Targeted A100@D500: 100 independent microscopic avalanches only in selected important cells; original 84k avalanche files and all old D20/D500 drift tables remain immutable.
- Physical quantity for now = ion births and their hemispheric/top-bottom sector census, NOT collected cathode U/V strip signal, ion drift or B-field effects. Four geometries: p3_noFW, p3_FW75, p4_noFW, p4_FW75; B=0.

## Original complete 120-muon independent validation

120 distinct original Geant4 J/psi parent events and physical muon tracks from frozen source chunk2, disjoint from previous 20-muon pilot chunk1. Exact TrackHeed produced 2,731 physical clusters, 5,881 accepted conduction electron starts; these SAME Heed electrons were used for all four full Garfield geometries. Total 23,524 true primary×geometry full Garfield microscopic attempts, zero capped. Eleven p4_FW75 full-reference primary electrons were physically born inside the 75um field-wire metal and therefore correctly produced zero ion births, with Garfield successful=0, not a computational avalanche cap. Original raw unsuccessful flags were NOT changed.

Original immutable Google Drive folders:
- Full real Geant4/Heed/Garfield: https://drive.google.com/drive/folders/1f8DsQwDbD6Kyz42mVZFbNyqZnz00M7sb
- 120-muon A20@D500 comparison: https://drive.google.com/drive/folders/176fe_fEUghfm3giCFD4Hp5eBXoSaUtYn
- Original 84k avalanche data: https://drive.google.com/drive/folders/1R-GxGNkwH8LWqqYlDKHyPhc1jm56ubiB

For geometries p3_noFW, p3_FW75, p4_noFW, p4_FW75 in that order, ratio average A20/full gain = 0.97863, 1.00114, 1.01689, 1.00631. Nominal per-muon 90% predictive coverage = 70.83%, 82.50%, 74.17%, 76.67%, i.e. mean gain good but uncertainty distribution systematically undercovers. A_hemi / A_TB correlations roughly 0.84–0.92, not validated as exact per-avalanche asymmetries. For one Heed electron per original cluster, 90% predicted cluster interval covers ~89.6–91%; for >=11 cluster electrons, ~67.4–74.7%.

Timing benchmark: original 11/12 full 10-muon chunks timed, 7,036.83 seconds for 110 actual Geant4 muons in 4 geometries (~63.97 s/muon) including TrackHeed, gas initialization and CSV output. Original recovered eighth chunk has NO recorded wall time; do not invent it. Fast warm A20 conditional 120 muon x four geometry pass ~0.01482s and original 84k cold load/SHA ~0.532s. Not a like-for-like end-to-end benchmark because Heed source is already cached in LUT replay.

## Why A20 undercovers: real Heed electron clustering

81.07% of original true Heed primaries share their 0.1x0.2mm A20 cell with another electron in the SAME muon; 44.55% lie in cells holding >=10 electrons within one muon. Of 224 actual Heed clusters with >=4 accepted electrons, 222 have at least half their electrons in one cell. Finite original 20-avalanche empirical cell mean error is COMMON to all those same-cell electrons, with model variance sum_cells n_cell(muon)^2 * s_cell(A20)^2 / 20, not just proportional to n_cell. This explains much but not all of the undercoverage: there may be within-cell spatial gradients, rare big avalanches and exact wire-metal start effects.

Target ranking uses only original source Heed start coordinates and the old A20 per-cell sample SD (NOT original full-Garfield muon gain outcomes), by sum over 120 muons of n(muon,cell)^2 * s_cell^2 /20. Top 50 cells per geometry explain ~72.8–76.6% of estimated finite-A20 mean-uncertainty contribution ON THIS SAMPLE.

## NEW 16,000 additional avalanche plan — CODE ONLY, NO NEW PHYSICS YET

50 top cells for each of four geometries (200 total). Add 80 independent full Garfield microscopic avalanches for each (16,000 NEW true avalanches) at previously UNUSED original D500 starts. Selected 200 cells become A100@D500 with original 20 plus 80; all other 4,000 cells remain A20@D500. Forty immutable SHA-checkpointed 400-trial native C++ batches, resumable without overwriting previous physics. Old D20/D500 drift data, old A20 84k avalanche data and full 120-muon reference must never be rewritten.

GitHub paths:
- test/mwpc-stage-c/pi_A100_targeted_50cell_campaign.py (prepare, preflight, run, finish, publish, status)
- test/mwpc-stage-c/test_pi_A100_targeted_campaign.py (offline synthetic check)
- test/mwpc-stage-c/pi_compare_A100_targeted_vs_A20.py (compare, publish; replays SAME original 120 Heed sources without new Heed/full muon physics)
- test/mwpc-stage-c/PI_A100_TARGETED_A20_UNDERCOVERAGE_V1.md (detailed methods)
- test/mwpc-stage-c/pi_run_lut500_overnight.py (ORIGINAL FROZEN selector, DO NOT EDIT)

## Already identified/fixed bug in the new A100 scripts

First user's WSL attempt on old commit 91e733a stopped during OFFLINE test with ValueError "Bad avalanche pass" because original pi_run_lut500_overnight.selected_j only accepts pass IDs 0..19. This guard is intentional for original A20 and must NOT be modified. New targeted helper additional_original_D500_sample(g,cell,p) handles passes 20..99 with the same SHA-stable original offset and coprime stride137:

    original_D500_index = (offset_for(geom,box) + 137*pass_id) % 500

Original 20 indices and new 80 indices are disjoint for all 4,200 original cells. The user HAS NOW RUN the updated WSL synthetic offline test and it PASSED (all 4,200 cells). The same user's next step prepare FAILED due to a SECOND simple Python interface bug: the targeted prepare() wrongly did `_,meta_original,src=payload_rows()` even though original holdout `payload_rows()` returns exactly TWO values (rows, manifest), whereas `verify_full()` returns THREE (rows, stats, SHA). This was fixed in GitHub (commit ce9774b9) to `src,meta_original=payload_rows()` and `ref_rows,full_stats,_=verify_full()`, and cross-check all 120 physical muon IDs. A new opt-in WSL real-source-only test `test_pi_A100_targeted_campaign.py --real-source-check` now asserts these actual interfaces and verifies original input SHA/schema without launching Garfield physics. The previous preparation failed BEFORE creating the new A100 ROOT directory; no clean-up is necessary. User has NOT YET confirmed the revised prepare/preflight on WSL, and no new A100 full avalanches have been run.

## Follow-up static audit: exhaustive C++ preflight (2026-10-30)

- GitHub commit 46eb88d2d9b490ad0a1596618c48abbded47d97a changes `pi_A100_targeted_50cell_campaign.py` to validate **ALL 40 x 400 = 16,000 prepared original D500 rows** with the actual native `--validate-lut-input` C++ parser, rather than only the first and last block. C++ hard limit is 432 rows/call, so all blocks use 400.
- The `run` mode now calls this full `preflight()` before ANY new Garfield computation, including on resumes. The original D500, 84k A20 and 120-muon physics remain unchanged.
- Source-readback verified matching C++ parser flag, 18-column header and success marker. An isolated mocked dispatch test checked exactly 40 batches / 16k row indices; **this is NOT a successful user-WSL actual `--real-source-check`, `prepare`, or native C++ preflight**.
- WSL local access is not available from the assistant's tool runtime. The user must execute the exact next commands below on their WSL and report actual output; until then A100 physics remains NOT STARTED. Do not mark any WSL gate PASS on the strength of the isolated static/mock check.

## Published Drive diagnostics: completed A100/A20 (2026-10-30)

**CONFIRMED publications through authenticated Google Drive connector:**
- A100 aggregated physical trial folder: https://drive.google.com/drive/folders/1cqY2mNpdOBhu4rO6DyFz529nYf8bSdSo
- A20 vs A100 comparison folder: https://drive.google.com/drive/folders/1ppAtbBOvAWrn6rK3dn9Ww0yf0wNQh1C7
- Original A20 84k folder: https://drive.google.com/drive/folders/1R-GxGNkwH8LWqqYlDKHyPhc1jm56ubiB
- Downloaded all four new A100 aggregated trials, all four original A20 trials, per-muon comparison, target list and cell rank; **each downloaded byte file SHA256 matches its original published manifest digest**. A100 `summary.json` SHA also matched its published manifest. All 16k new avalanches and 40/40 chunks complete, zero capped/nonmetal censored and zero recorded metal starts on new 16k, new physics wall-clock sum 7287.985900988278 s (~2.024 h). Do not confuse this with 11 previously observed physically in-metal start electrons in the original 120-muon p4_FW75 reference.

**Read-only new scientific diagnostics on 480 frozen muon×geometry comparison rows:**
- Transitions missed→covered / covered→missed: p3_noFW **9/5**, p3_FW75 **4/5**, p4_noFW **7/7**, p4_FW75 **8/2**, total **28/19**. Pooled covered 365/480 → 374/480 = **+1.875 percentage points** (480 not independent parents).
- 20,000-resample *paired cluster bootstrap by 120 physical Geant4 parent muons* (resample all 4 geometries together), seed 20261030, gives 95% percentile CI **[-0.833,+4.792] percentage points**, contains zero. This is source-conditioned and does not account for ranking selection or replay Monte Carlo uncertainty.
- Old full Garfield outside A20 interval, lower/upper: p3_noFW 13/22, p3_FW75 8/13, p4_noFW 14/17, p4_FW75 16/12. New outside A100 interval lower/upper: p3_noFW 10/21, p3_FW75 10/12, p4_noFW 12/19, p4_FW75 10/12. No single uniform lower-versus-upper miss pattern.
- Median old→new 90%-PI widths (ion births): p3_noFW 11160.85→11467.95, p3_FW75 117828.60→118131.875, p4_noFW 30281.475→29675.775, p4_FW75 113159.30→110949.05. Average absolute per-muon mean-gain prediction error falls on all four geometry samples, but intervals remain poorly calibrated.
- **Natural no-target control:** 38/32/43/45 of 120 muons for p3_noFW/p3_FW75/p4_noFW/p4_FW75, respectively, touch *no* A100-selected cells; old and new underlying response pools for these muons are identical. Even among these unchanged observations, independently drawn 200-replica A20/A100 predictions flip 2 misses to covered and 2 covered to misses in aggregate. Common-random-number (CRN) replay is needed to isolate sampling noise.
- Source-selected Heed electron fraction across the 480 geometry comparisons has median about **6.7%**, 25th percentile **0**, 75th percentile **32.2%**; 158/480 observations touch no selected A100 cells. The top50 estimated finite-A20 cell-mean variance share (73–77%) must not be equated with fraction of electrons (26.5–29.2%) or total predictive error.

**Independent old20 versus new80 full Garfield ion birth samples in the 200 target cells:**
- Original first 20 sample mean exceeds additional independent 80 sample mean in **140/200 cells** (by geometry p3_noFW 37/50, p3_FW75 34/50, p4_noFW 37/50, p4_FW75 32/50).
- Frozen 120-Heed-electron-count-weighted per-selected-cell mean ion births (old20→additional80): p3_noFW 622.31→574.64 (**−7.66%**), p3_FW75 6031.27→5694.24 (**−5.59%**), p4_noFW 1653.18→1483.37 (**−10.27%**), p4_FW75 5732.17→5358.90 (**−6.51%**). The drop is compatible with noisy old20 estimates / positive-skew gain tail / outcome-blind high-SD cell selection, but it does NOT uniquely establish a cause.
- Selected-cell original ranking score reconstructed replacing old20 `sd^2` with independent new80 `sd^2` (fixed original source n-cell weights) is 85.1%/98.8%/98.1%/89.2% of the old scores for p3_noFW/p3_FW75/p4_noFW/p4_FW75 respectively. These are *within-target score ratios*, NOT overall error reduction or validation of 90% coverage.

**Next steps**: preserve all existing physics and associated directory, then run new *purely statistical* CRN A20/A100 replay with more than 200 replicas; hierarchical cell-pool bootstrap to propagate `n_cell^2` correlated training uncertainty; check within-bin start dependence and true heavy-tail behavior; then independent third Geant4 source sample for external generalization. No new Garfield calculations authorized for this diagnostic.

**Archival caution:** Published A100 folder currently has aggregates and targeting files, **not the 40 chunk directories including physics.log, ion-birth reservoirs, and per-start timing CSVs**. Leave original WSL SHA-sealed chunks intact; make another immutable complete archive before any cleanup, especially for eventual spatial cloud studies.

## User-confirmed A100 vs A20 comparison (2026-10-30)

User pasted the complete successful WSL `pi_compare_A100_targeted_vs_A20.py compare` output and `summary.json`. It reports `PASS A100 target-cell sensitivity analysis`, 120 physical muons, 5,881 primary TrackHeed electrons, 200 conditional LUT replicas per muon×geometry, original A20 84k unchanged, SHA-sealed A100 16k unchanged, original full 120 Garfield unchanged. Comparison `status=complete` and targeted A100 full manifest SHA256 `2a8520c61dc38dc2b32aa31319c3665e36d9ed86f9119e8f3737f32849d7da82`; original A20 SHA256 `8ba639e3f082192d716d328775e13743ca72f9035ebfe703b1f722a202fb006c`; full 120 reference SHA256 `dfffdeabfac813f339a1746480a99e19aa78ee8f376b2a8683f42e2fc393a763`.

Measured A20 -> targeted A100 90%-interval coverage (out of 120 true full Garfield muons per geometry):
- p3_noFW: 85/120 (70.833%) -> 89/120 (74.167%), Δ +4/120 (95% paired bootstrap [-3.33,+9.17] percentage points); A100/full mean-gain ratio 0.960195 vs old 0.978626; target Heed fraction 0.291787.
- p3_FW75: 99/120 (82.5%) -> 98/120 (81.667%), Δ -1/120 (95% [-5.83,+4.17] pp); A100/full mean-gain ratio 0.989470 vs old 1.001144; target Heed fraction 0.292297.
- p4_noFW: 89/120 (74.167%) -> 89/120 (74.167%), Δ 0 (95% [-6.67,+6.67] pp); A100/full mean-gain ratio 0.991593 vs old 1.016894; target Heed fraction 0.267301.
- p4_FW75: 92/120 (76.667%) -> 98/120 (81.667%), Δ +6/120 (95% [0,+10.00] pp); A100/full mean-gain ratio 0.990905 vs old 1.006311; target Heed fraction 0.265091.
- Naive pooling: 365/480 (76.04%) -> 374/480 (77.92%), +9/480 (+1.875 pp); **NOT 480 independent muons**, same 120 sources reused over four geometries, and source coordinates informed target selection. Do NOT claim significant overall improvement without a 120-parent clustered paired test; per-geometry bootstrap intervals all include 0.
- For A_hemi coverage A20->A100: p3_noFW 84.17->87.50%; p3_FW75 85.83->85.83%; p4_noFW 84.17->86.67%; p4_FW75 81.67->81.67%.
- For A_TB coverage A20->A100: p3_noFW 84.17->81.67%; p3_FW75 81.67->85.00%; p4_noFW 83.33->85.83%; p4_FW75 81.67->82.50%. No uniformly consistent calibration improvement.

**Scientific interpretation**: increasing from 20 to 100 independent microscopic avalanches in 50 important cells/geometry only weakly improved conditional gain prediction; all geometries still under nominal 90%, and the changes are not individually compelling. The top-50 73–77% score is a *noisy, A20-estimated finite-training-mean variance share*, NOT a share of total predictive error. Only 26.5–29.2% of actual Heed electrons lie in priority cells. Rankings are source-coordinate-conditioned (not fully unseen validation). A20 and A100 original comparison scripts use different Monte Carlo replay random seeds and only 200 replicas, so compare changes include some stochastic replay noise. The A100 mean gain falls in all four geometries; this may be a finite-training / selected-sample effect, NOT proven physics. Eleven original p4_FW75 physically in-metal Heed starts remain unmodelled by the pure empirical LUT and are not to be overwritten or silently reclassified.

**Next steps**: preserve all on-WSL outputs; publish both A100 final physical data and SHA-sealed comparison with existing `pi_A100_targeted_50cell_campaign.py publish` and `pi_compare_A100_targeted_vs_A20.py publish` after user command. Next use per-muon comparison CSV to diagnose missed interval above/below, paired gained/lost 90%-coverage, interval widths, selected-source fraction and mean bias, and calculate 95% aggregate CI with bootstrap clustered by 120 *parent muons*, not 480 independent geometry rows. Consider same-random-number A20/A100 replays with more than 200 replicas to isolate LUT changes from replay Monte Carlo error. Do NOT generate another wave of Garfield physics until these diagnostics discriminate training uncertainty, within-cell spatial gradients, gain tails and field-wire metal starts.

## User-confirmed full A100 physics completion (2026-10-30)

- User pasted actual WSL output from `pi_A100_targeted_50cell_campaign.py status`: `PLAN prepared 16000`, `TARGETED GARFIELD 40/40 blocks`, `SEALED COMPLETE True`.
- Same user's run log tail reports successful SHA-sealed checkpoint messages for chunks 12..40, including `PASS TARGETED 40/40, +400 independent Garfield avalanche starts` and the final message `PASS TARGETED A100 full microscopic Garfield complete, 16,000 new original D500 starts, SHA-sealed without modifying original 84k or 120 real muons`.
- Log excerpt shows `physically embedded metal starts 0` for chunks 12..40; do not generalize this to chunks 1..11 without their complete logs or `summary.json`.
- **A100 NEW PHYSICS COMPLETE per user-provided native execution and final SHA manifest.** Assistant has not independently read WSL chunk manifests/outputs; source verification is the next comparison's job. Do NOT rerun `prepare`, `run` or `finish` and do NOT delete the new A100 output directory.
- Next safe step: `python3 -u test/mwpc-stage-c/pi_compare_A100_targeted_vs_A20.py compare`; this re-verifies all 40 SHA chunks, original 84k A20 and frozen 120 physical Geant4/TrackHeed/Garfield reference, then performs 200 stochastic LUT replays per muon and geometry. This is NOT new Garfield physics. Check `summary.json`, then publish only after successful comparison and integrity checks.
- Important inference limit: 120-muon *coordinates* were used for target ranking. Treat comparison as source-conditioned A100 sensitivity, not independent generalization. For external validation use unseen third Geant4 source sample.

## User-confirmed WSL gate outcome (2026-10-30, after pull to 9f6225a)

- User pasted actual shell output proving `py_compile` succeeded and `test_pi_A100_targeted_campaign.py --real-source-check` printed PASS ORIGINAL REAL WSL INPUTS and PASS SYNTHETIC/OFFLINE ONLY, including all 4,200 spatial cells and original A20 / D500 / 120 Geant4 + Heed SHA contracts.
- `pi_A100_targeted_50cell_campaign.py prepare` ACTUALLY PASSED, creating `test_runs/mwpc_stage_c/pi_A100_targeted_50cells_per_geometry_v1` and SHA-sealed 16,000 new source starts; target 50 scores cover 76.554%, 73.354%, 73.394%, 72.804% of modeled finite-A20 mean variance for p3_noFW, p3_FW75, p4_noFW, p4_FW75.
- C++ preflight ACTUALLY PASSED only first and last 400 rows because that WSL pull ended at `9f6225a`, BEFORE subsequent GitHub fix `46eb88d2` introduced ALL-40-block preflight and a preflight gate in `run`. WSL status: PLAN prepared 16000, TARGETED GARFIELD 0/40, SEALED COMPLETE False.
- **Do NOT run `prepare` again**: it correctly refuses overwrite of existing prepared root. Next action: WSL `git pull --ff-only` to include exhaustive preflight fix, then `python3 test/mwpc-stage-c/pi_A100_targeted_50cell_campaign.py preflight`; require ALL 40 blocks / 16,000 rows PASS before launching new Garfield physics. No new A100 avalanches reported.

## Exact immediate NEXT WSL command (NO new Garfield avalanches)

    cd ~/na6p/src/NA6PRoot-mwpc
    git pull --ff-only && \
    python3 -m py_compile test/mwpc-stage-c/pi_A100_targeted_50cell_campaign.py test/mwpc-stage-c/pi_compare_A100_targeted_vs_A20.py test/mwpc-stage-c/test_pi_A100_targeted_campaign.py && \
    python3 test/mwpc-stage-c/test_pi_A100_targeted_campaign.py --real-source-check && \
    python3 -u test/mwpc-stage-c/pi_A100_targeted_50cell_campaign.py prepare && \
    python3 test/mwpc-stage-c/pi_A100_targeted_50cell_campaign.py preflight

Proceed with NO new physical Garfield calculations unless BOTH PASS PREPARED and actual compiled C++ input-parser PASS appear. If failure: inspect full trace and fix code, do NOT delete original physics or disable SHA checks.

After all safe checks pass, run on user's WSL:

    nohup nice -n 15 python3 -u test/mwpc-stage-c/pi_A100_targeted_50cell_campaign.py run --hours 6.75 --max-chunks 40 > "$HOME/na6p/A100_targeted_50x4_20261030.log" 2>&1 < /dev/null &

Check status:

    python3 test/mwpc-stage-c/pi_A100_targeted_50cell_campaign.py status
    tail -n 25 ~/na6p/A100_targeted_50x4_20261030.log

After 40/40 SHA-sealed complete:

    python3 -u test/mwpc-stage-c/pi_compare_A100_targeted_vs_A20.py compare
    python3 test/mwpc-stage-c/pi_A100_targeted_50cell_campaign.py publish
    python3 test/mwpc-stage-c/pi_compare_A100_targeted_vs_A20.py publish

Evaluate original A20 vs selected A100/mixed A20 mean gain, 90% gain intervals, asymmetry A_hemi/A_TB and timed speed; 120 source coordinates informed target ranking, hence it is a source-conditioned recheck, not independent detector-wide generalization. For final external validation use a newly unseen third physical Geant4 source sample. Avoid inventing numerical improvement or claiming unrun WSL physics complete.
