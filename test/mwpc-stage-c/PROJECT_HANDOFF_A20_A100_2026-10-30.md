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
