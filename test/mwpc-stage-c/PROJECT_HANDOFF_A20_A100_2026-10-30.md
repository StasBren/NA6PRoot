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

Original 20 indices and new 80 indices are disjoint for all 4,200 original cells. The new GitHub test exhaustively checks this; correction committed in branch feature/mwpc-stage-c, latest fix commit 6327ce09. An isolated arithmetic check of all 4200 cells passed, but actual UPDATED WSL offline test has not yet been reported. Because user's command chained with &&, prepare and preflight never started in the failure: no recovery or file deletion is needed.

## Exact immediate NEXT WSL command (NO new Garfield avalanches)

    cd ~/na6p/src/NA6PRoot-mwpc
    git pull --ff-only && \
    python3 -m py_compile test/mwpc-stage-c/pi_A100_targeted_50cell_campaign.py test/mwpc-stage-c/pi_compare_A100_targeted_vs_A20.py test/mwpc-stage-c/test_pi_A100_targeted_campaign.py && \
    python3 test/mwpc-stage-c/test_pi_A100_targeted_campaign.py && \
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
