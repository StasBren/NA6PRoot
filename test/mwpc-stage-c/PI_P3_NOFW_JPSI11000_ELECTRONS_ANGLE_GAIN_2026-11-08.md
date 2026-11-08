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

## Zero-gain conditional mean check

For all original 588560 full-Garfield calls, the recorded summary shows:
- Collected electron mean over ALL primary seeds = 566.6484436591002
- Zero collected-electron outcomes = 91708
- Active collected-electron outcomes = 496852
- Conditional mean collected_e GIVEN collected_e>0 =
  566.6484436591002 * 588560 / 496852 = **671.2393388775732**
- Garfield ne mean over ALL seeds = 566.3984385619138
- Zero Garfield ne outcomes = 91699
- Conditional mean ne GIVEN ne>0 = **670.9310350379684**

These are mathematically exact consequences of the already-available
sample means and event counts, up to printed floating point precision.
The unconditional mean is the physically correct effective average
gain per genuine TrackHeed starting electron when modeling actual
charge response; conditional mean only describes successful avalanches.

To verify numerically from every original SHA-sealed seed CSV and
regenerate the linear plot with both means in the legend (without
rehashing 60M ion births or running Garfield):

~~~bash
cd ~/na6p/src/NA6PRoot-mwpc
git pull --ff-only
PY=test_runs/mwpc_stage_c/.venv_pi_slides/bin/python
SCRIPT=test/mwpc-stage-c/pi_p3_noFW_jpsi11000_gain_angle_plots_v1.py
"$PY" -m py_compile "$SCRIPT"
"$PY" "$SCRIPT" selftest
"$PY" -u "$SCRIPT" electrons-linear
~~~

New, nonoverwriting output
`test_runs/mwpc_stage_c/pi_p3_noFW_jpsi11000_electrons_linear_v3_nonzero_mean/`
has `statistics.json` containing both all-inclusive and conditional
means/RMSs, plus their zero fractions.

## Histogram fix v5 — no CDF, separate zero bin, linear and log Y

The previous user-local patch added a return schema without the
`RMS_given_positive` field, causing a `KeyError` **after** writing
its histogram files. The GitHub source now has the correctly matching
producer/consumer schema (including positive-only RMS).

Mode `electrons-hist` reads the ORIGINAL 550 SHA-sealed native 20-muon
source chunks, including replay 0426, but checks only source selection,
crossing and seed CSVs: **does not scan/read the expensive ion-birth
reservoir** and does not execute Garfield.

Two PNGs, **NO CDF**:
- `01_electrons_per_avalanche_with_zero_bin.png`: explicitly separate
  single-electron-count zero bin `[0,1)`, and nonzero bins `[1,...]`.
  Probability per bin with normalization over ALL 588560 avalanches.
- `01_electrons_per_avalanche_with_zero_bin_without_zero_bin.png`:
  removes bin `[0,1)`; normalizes over positive-count avalanches only.

Each PNG contains **left linear-Y** and **right logarithmic-Y** panels,
both with linear X over full observed counts. Legend contains
mean(all), mean(>0), RMS(all) and number of zero outcomes, separately
for Garfield ne and anode-collected electrons. JSON adds
`RMS_given_positive`, positive counts and zero fractions.

Output is new, immutable:
`test_runs/mwpc_stage_c/pi_p3_noFW_jpsi11000_electrons_hist_v5/`.
The full `plot` mode creates another independent
`pi_p3_noFW_jpsi11000_gain_angle_hist_v5/` with the same two
electron histograms plus ion-birth azimuth and conditional gain maps.

Because an earlier assistant instructed making local direct edits to
the same tracked script, ordinary `git pull --ff-only` **may refuse**
to update the dirty tracked file. First back up the existing local
script outside the Git checkout, then fetch the GitHub branch and
extract the exact reviewed file from `FETCH_HEAD` (see chat commands).
No original Garfield data or existing previous plot outputs are touched.
