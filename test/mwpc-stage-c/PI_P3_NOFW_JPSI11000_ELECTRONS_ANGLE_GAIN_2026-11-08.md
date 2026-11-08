# p3_noFW: 11,000 J/psi muons — microscopic electrons, ion-birth angles and gain

Read-only analysis: 1000 original + 10,000 additional distinct J/psi parents;
588,560 real microscopic Garfield avalanche calls. All 550 chunk receipts
and consumed physics seed/birth files are checked against SHA hashes,
including the completed high-Heed recovery chunk 0426. Nothing rewritten.

Three figures:

1. 01_electrons_per_avalanche.png — distribution of exact native
   GetAvalancheSize(ne) and anode-collected endpoint counts per primary
   TrackHeed electron. Positive-count histogram uses log axes; CDF includes
   zero-avalanche cases. These two electron counters are NOT interchangeable.
2. 02_ion_birth_angles_around_anode.png — phi=atan2(v,u-nearest_3mm_anode)
   ion-birth probability in 5-degree bins, ion-weighted and separately
   avalanche-weighted, plus circular-mean avalanche direction. This comes
   from up to 128 randomly saved ion births per avalanche, each weighted by
   true full birth count / saved count. Per-bin azimuth is a statistical
   estimator, NOT an exact 333M-ion position census; R<0.2 centroids
   are excluded and reported.
3. 03_mean_collected_gas_gain.png — effective mean gas gain:
   total collected electron endpoints / total genuine Heed initial
   electrons = mean collected_e per one starting electron; also
   conditional mean gain versus exact start coordinates u,v (30x30 bins).
   Number of electron-ion pair creations per starting electron is given
   separately. The ratio collected electrons / created pairs is a
   charge-balance diagnostic, NOT gas gain.

Native operational definition: collected_e counts electron endpoints
within 1.5x the physical anode radius of nearest wire. For a 30um
anode diameter this is a 22.5um radial cutoff.

Run on WSL after syncing branch:

~~~bash
cd ~/na6p/src/NA6PRoot-mwpc
git switch feature/mwpc-stage-c
git pull --ff-only
PY=test_runs/mwpc_stage_c/.venv_pi_slides/bin/python
SCRIPT=test/mwpc-stage-c/pi_p3_noFW_jpsi11000_gain_angle_plots_v1.py
"$PY" -m py_compile "$SCRIPT"
"$PY" "$SCRIPT" selftest
set -o pipefail
"$PY" -u "$SCRIPT" plot | tee ~/na6p/p3_noFW_jpsi11000_gain_angle_plots.log
~~~

Original large birth-position sample files take time to read and SHA-check.
No new Garfield physics is executed.

New non-overwriting output:
test_runs/mwpc_stage_c/pi_p3_noFW_jpsi11000_gain_angle_v1/
  01_electrons_per_avalanche.png
  02_ion_birth_angles_around_anode.png
  03_mean_collected_gas_gain.png
  angular_distribution.csv
  mean_gain_map.csv
  summary.json

Note that these concern births and collection, not positive-ion drift,
Shockley-Ramo induced whole-cathode signals or U/V strip reconstruction.
Only hemisphere asymmetry is considered in the project's other analysis.

## Linear-scale revision (2026-11-08)

The source file `pi_p3_noFW_jpsi11000_gain_angle_plots_v1.py` was updated
in place. The **physics inputs are unchanged**. Electron multiplicity
is plotted with **linear X and Y**, no percentile truncation, including
zero-electron avalanches. Legend of both histogram and CDF reports the
full-sample mean and ROOT-style RMS (`np.std(ddof=0)`) separately for
`garfield_ne` and `collected_e`.

Quick redraw of ONLY requested electron-number figure, avoiding the
60,780,188-line ion-birth reservoir CSV parsing:

~~~bash
cd ~/na6p/src/NA6PRoot-mwpc
git pull --ff-only
PY=test_runs/mwpc_stage_c/.venv_pi_slides/bin/python
SCRIPT=test/mwpc-stage-c/pi_p3_noFW_jpsi11000_gain_angle_plots_v1.py
"$PY" -m py_compile "$SCRIPT"
"$PY" "$SCRIPT" selftest
"$PY" -u "$SCRIPT" electrons-linear
~~~

New immutable output:
`test_runs/mwpc_stage_c/pi_p3_noFW_jpsi11000_electrons_linear_v2/`
(figure and `statistics.json`). SHA source file checks are retained;
no new Garfield runs and no edits to old images.

For ALL three figures using the updated styling (including a linear
gain-map colorbar) run `"$PY" -u "$SCRIPT" plot`. It creates a second,
distinct output directory
`test_runs/mwpc_stage_c/pi_p3_noFW_jpsi11000_gain_angle_linear_v2/`.
This full mode must read all 60.8M source birth-reservoir rows and is
much slower, but still runs no Garfield simulations.
