# Penning ON r=0.57: SAME 1000 real J/psi muons — controlled campaign

Original reference: pi_p3_noFW_jpsi1000_native_v1 (Penning OFF).
New output: pi_p3_noFW_jpsi1000_penning057_v1 (Penning ON).
The new runner copies the original 1000 true Geant4 inputs BYTE FOR BYTE;
it does not replace or modify any original baseline file.

Only physics change: the previously validated independent Penning pilot
executable enables gas.EnablePenningTransfer(0.57, 0., "ar").
Geometry is the original p3_noFW: 3mm pitch, 30um wires, cathodes
at +/-3mm, +1800V, Ar/CO2 70:30, B=0. It is NOT the actual
2026 prototype with 4mm pitch and asymmetric 2+4mm cathode gaps.

### 0. Fetch ONLY the new script, without merging dirty plotting files

```bash
cd ~/na6p/src/NA6PRoot-mwpc
git fetch origin feature/mwpc-stage-c
SCRIPT=test/mwpc-stage-c/pi_p3_noFW_jpsi1000_penning057_campaign_v1.py
test ! -e "$SCRIPT" || { echo 'Inspect already-existing script; refusing overwrite'; exit 1; }
git show "FETCH_HEAD:$SCRIPT" > "$SCRIPT"
PY=test_runs/mwpc_stage_c/.venv_pi_slides/bin/python
"$PY" -m py_compile "$SCRIPT"
"$PY" "$SCRIPT" selftest
```

The dependencies (original frozen 1000 runner, successful Penning
3-muon pilot C++ source+binary, and Python pilot verifier) must already
exist locally. If not, do not attempt a partial physics campaign.

### 1. Prepare and preflight, then the FIRST 20 (not 1000!)

```bash
"$PY" "$SCRIPT" prepare
"$PY" "$SCRIPT" preflight
set -o pipefail
"$PY" -u "$SCRIPT" run --max-chunks 1 --hours 24 | tee ~/na6p/p3_noFW_penning1000_first20.log
"$PY" "$SCRIPT" gate
"$PY" "$SCRIPT" status
```

The first 20 muons are fully recomputed with true Garfield Penning ON.
Every chunk validates full microscopic counts, uncapped avalanche size,
original Heed electron identities, cluster indices, start position,
energy, time, original Garfield seeds, and 128-point/avalanche reservoir.
The gate requires gain ON/OFF ratio between 2 and 30 plus zero
size-capped events. This catches inadvertent wrong source/binary but
is NOT an experimental calibration of the chamber gain.

The runner refuses to start more than 20 new muons before the first
successful SHA-sealed chunk. If the first chunk fails, keep all files;
never delete or rerun unsealed physics without root-cause analysis.

### 2. Only after passing the first-chunk gate: remaining 980 muons

```bash
"$PY" -u "$SCRIPT" run --max-chunks 49 --hours 24 | tee ~/na6p/p3_noFW_penning1000_remaining.log
"$PY" "$SCRIPT" status
"$PY" "$SCRIPT" audit
```

49 existing/remaining chunks are skip-audited or computed and then
SHA sealed. User's original 15-hour Penning OFF physics is never changed.
Input geometry, path, Heed starts and random seeds are compared exactly.
Hours argument is a safe stop ceiling, not a duration prediction.

### 3. Plot gain and angular ion-birth distributions after completion

```bash
"$PY" -u "$SCRIPT" plots | tee ~/na6p/p3_noFW_penning1000_plots.log
```

Plots are in
test_runs/mwpc_stage_c/pi_p3_noFW_jpsi1000_penning057_v1/plots_collected_gain_and_ion_birth_angles_v1/

1. 01_penning_on_avalanche_gain_with_zero.png: collected gain all
   true Heed electrons; first zero-only bin [0,1), positive bins 50
   collected electrons wide; separate linear-Y and log-Y panels;
   legend gives mean(all), RMS(all).
2. 02_penning_on_avalanche_gain_positive_only.png: zero bin absent,
   normalized over collected gain >0; legend gives mean and RMS
   conditional on successful collection.
3. 03_penning_on_ion_birth_angular_equal_avalanche.png: only ORANGE
   curve, 5-degree bins around nearest 3mm anode; every ion-producing
   avalanche has total statistical weight 1, legend present.
4. 04_penning_on_ion_birth_angular_ion_weighted.png: only ORANGE
   curve with each ion birth given equal statistical weight, via
   128-point per avalanche unbiased spatial reservoir.
5. angular_distribution.csv and summary.json (paired Penning
   ON and OFF gain, exact total ion births, zero-gain fractions, SHA).

Angular histograms are estimates from 128 uniformly sampled ion
birth positions per avalanche, not a census of every ion position.
Full ion birth counts and their source SHA receipts are audited.
No ion drift, induced cathode charge, or U/V strip signals are included.

ATTENTION: 1000 independent Geant4 J/psi parents; individual
Garfield avalanches within the same muon share initial conditions.
