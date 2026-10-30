# A20@D500 undercoverage and targeted A100 augmentation

## Physical diagnosis from 120 actual Geant4 muons

120 independent physical J/psi muons, 2,731 accepted Heed clusters, 5,881
conduction electrons, 4 Garfield geometries, 23,524 full seed simulations.
A20@D500 = only 20 independently simulated full Garfield avalanches per
0.1x0.2 mm spatial cell, selected among 500 ORIGINAL DRIFT points/cell.
There are 4,200 cells and 84,000 original microscopic avalanches.

Nominal 90% conditional gain coverage on 120 new muons:
- p3_noFW 70.83%, p3_FW75 82.50%,
- p4_noFW 74.17%, p4_FW75 76.67%.

Cluster coverage declines sharply as accepted Heed electron multiplicity
increases. Single-electron cluster coverage approximately 89.6-91.0%;
clusters with >=11 primaries only 67.4-74.7%.

Source coordinate mapping to the true exact C++ station wire phase confirms:
- 81.07% of Heed primaries in a real muon occupy a 0.1x0.2mm cell with
  at least one other electron of that SAME muon.
- 44.55% occupy a cell with at least ten electrons of that muon.
- Of 224 physical Heed clusters with >=4 accepted electrons, 222 have
  at least half the electrons in the same cell.

Repeated use of the SAME uncertain 20-trial empirical cell mean produces
a shared nuisance error whose variance scales quadratically with
the number of real Heed electrons n in that cell:

    Var_trained_cell_mean_per_muon ≈ sum_cell n(muon,cell)^2 s(cell)^2 / 20

By comparison the stochastic independent-response variance scales as n.
This conditional finite-training error explains much, but not necessarily
all, of the observed undercoverage; within-bin position dependence,
rare high-gain avalanches and exact field-wire metal starts also matter.

Adding the modeled finite-cell-mean variance to per-muon LUT variance
as an explanatory (NOT independently validated) sensitivity check changes
the observed standard deviation of normalized gain residuals:
- p3_noFW 1.427 -> 1.211
- p3_FW75 1.278 -> 1.089
- p4_noFW 1.371 -> 1.160
- p4_FW75 1.252 -> 1.054.

## Outcome-blind target selection

The new script ranks cells for each geometry using ONLY immutable
source Heed coordinates and the ORIGINAL 20 original Garfield-trial
per-cell standard deviations, NEVER the 120 reference Garfield outcomes:

    score(geometry,cell) = sum over real muons n(muon,cell)^2 * s_A20(cell)^2 / 20

Highest 50 cells per geometry account for 76.55%, 73.35%, 73.39%,
72.80% of the respective estimated finite-training cell-mean variance
in p3_noFW, p3_FW75, p4_noFW, p4_FW75 on this source sample.
This is a source-exposure-ranked test, not a full detector angular PDF.

Selected targets:
- 50 highest-scoring original cells per geometry (200 total).
- 80 NEW independent full Garfield avalanche runs per cell (16,000 total).
- Original SHA-offset + stride137 original 500 DRIFT starting point
  indices from passes 20..99. Original 20 full-avalanche passes 0..19
  are excluded, so there is no duplicated original D500 trial index.
- These 200 cells become A100@D500; all other 4,000 original
  cells stay A20@D500. The original 84k physics files remain untouched.
- 40 separate 400-source Garfield chunks with physical completion,
  exact original-input hashes, resumability and strict no-overwriting
  of incomplete results. Any capped result stops physics validation.
  Only strictly proven nontransporting electron born within an original
  75um-diameter field-wire metal can be included as genuine zero response.
- Previously measured full Garfield + Heed time ~0.33 s per
  microscopic avalanche gives a VERY ROUGH ~1.5h runtime estimate
  for 16k new starts. Measure the actual wall time instead of
  treating this extrapolation as a timing guarantee.

New source:
- test/mwpc-stage-c/pi_A100_targeted_50cell_campaign.py
- test/mwpc-stage-c/pi_compare_A100_targeted_vs_A20.py
- test/mwpc-stage-c/test_pi_A100_targeted_campaign.py

## WSL commands: first prepare, then physics, then compare

All commands from repository root.

    cd ~/na6p/src/NA6PRoot-mwpc
    git pull --ff-only
    python3 -m py_compile test/mwpc-stage-c/pi_A100_targeted_50cell_campaign.py test/mwpc-stage-c/pi_compare_A100_targeted_vs_A20.py test/mwpc-stage-c/test_pi_A100_targeted_campaign.py
    python3 test/mwpc-stage-c/test_pi_A100_targeted_campaign.py --real-source-check
    python3 -u test/mwpc-stage-c/pi_A100_targeted_50cell_campaign.py prepare
    python3 test/mwpc-stage-c/pi_A100_targeted_50cell_campaign.py preflight
    python3 test/mwpc-stage-c/pi_A100_targeted_50cell_campaign.py status

The source test verifies BOTH the synthetic 4,200-cell disjoint
sampling logic and the ORIGINAL already existing 120 Geant4/TrackHeed/
Garfield source read contracts, 84k empirical A20 input SHA/schema,
and D500 source headers WITHOUT running new microscopic physics.
The earlier Python unpacking error in prepare() was fixed: original
payload_rows() returns exactly (rows, manifest), NOT three values.

This PREPARATION verifies original SHA and scans original frozen 500
drift-trial records/cell; it does not launch any Garfield physics.
Only proceed if it prints PASS PREPARED and 16,000 frozen new starts.
Do not run physics if there is an exception.

    nohup nice -n 15 python3 -u \
      test/mwpc-stage-c/pi_A100_targeted_50cell_campaign.py run \
      --hours 6.75 --max-chunks 40 \
      > "$HOME/na6p/A100_targeted_50x4_20261030.log" 2>&1 < /dev/null &
    echo "New microscopic Garfield PID: $!"

Monitor without altering physics:

    tail -n 25 ~/na6p/A100_targeted_50x4_20261030.log
    python3 test/mwpc-stage-c/pi_A100_targeted_50cell_campaign.py status

Full run creates 40 immutable SHA checkpoints of 400 original D500 starts
each. Repeat the same RUN command after a normal between-chunk time stop.
DO NOT delete or overwrite any incomplete chunk after a fault; preserve
the original log and CSVs for investigation. The script creates a final
SHA manifest only when all 16,000 new full-Garfield starts are verified.

On SEALED COMPLETE:

    python3 -u test/mwpc-stage-c/pi_compare_A100_targeted_vs_A20.py compare
    python3 test/mwpc-stage-c/pi_A100_targeted_50cell_campaign.py publish
    python3 test/mwpc-stage-c/pi_compare_A100_targeted_vs_A20.py publish

Compare the original 120 *same* independently generated Geant4 muon
trajectories and saved, unchanged Heed electron starts, 200 new random
conditional responses per muon+geometry, without rerunning true
full Garfield on muons. Compare both mean gain and conditional 90%
gain + angular birth asymmetry coverage. Report paired bootstrap
over 120 original independent parents. The 120 source Heed COORDINATES
guided target selection, so this is not a fully blind new test;
for definitive external generalization use an untouched new source
sample. Exact field-wire-metal veto is NOT silently retrofitted
into either original A20 or new A100 empirical sampler.

## Status

Scripts committed in GitHub feature/mwpc-stage-c. The user has NOT
yet executed new 16k original Garfield avalanches on their WSL.
No claim of improved calibration/coverage or confirmed C++ execution
is made before user runs source validation and physical campaign.
