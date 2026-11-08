# NA60+ MWPC p3_noFW: diagnosing low microscopic Garfield gain

The completed immutable Penning-OFF 11k real Geant4/TrackHeed muon
campaign has:
- 588,560 independent microscopic one-electron Garfield avalanche calls
  (correlated under only 11,000 independent parent muons).
- mean anode-collected electrons per original Heed seed: 566.648.
- mean given nonzero collected e: 671.239.
- about 15.58% zero-collected avalanches.

These values were already checked from exact all-11k native CSVs and do not
constitute an estimate of the expected Penning-ON physical Ar/CO2 gain.

## Confirmed omission in the sealed production source

The source 'phaseB_native_muon_p3_noFW_jpsi_runtime_v1.cxx' has gas:
  gas.SetComposition("ar",70.,"co2",30.);
  gas.SetTemperature(293.15);gas.SetPressure(760.);
  gas.SetMaxElectronEnergy(200.);gas.Initialise(false);
and does NOT enable Penning excitation-to-ionisation transfer.

Garfield++ Ar:CO2 microscopic documented setting:
  gas.EnablePenningTransfer(0.57, 0., "ar");
Compare UserGuide, Ar/CO2 example, and 2026 Garfield ROOT Forum
"Penning effect for Drift Tube for (Ar/CO2)".
This physical omission can suppress gain by large factors. Actual ratio in
this MWPC at 1800V is UNKNOWN until the independent native test is run.

## Geometry caveat — important for comparison with prototype

Current 'p3_noFW' archive:
- anode pitch 3 mm, diameter 30 microns
- cathodes v=+/-3 mm, +1800 V anodes, Ar/CO2 70/30.

Local 'New_chamber_design(1).pdf' p.2 Table 1:
- Prototype 2025: wire diameter 30 microns, p=4 mm, 3+3 mm
  cathode gaps, HV=2.1 kV.
- Prototype 2026: wire diameter 30 microns, p=4 mm, 2+4 mm
  cathode gaps, HV=1.8 kV.
So the reference at 1800 V with 3-mm pitch and symmetric 3+3 mm
cathode gaps is not directly one of these prototype configurations.
DO NOT copy nominal expected gains across different voltage/geometry
without measured gain-vs-HV reference.

## Safe Penning-only matched pilot

An independent C++ source and compiled target were added:
test/mwpc-avalanche/phaseB_native_muon_p3_noFW_jpsi_penning057_pilot_v1.cxx
test/mwpc-avalanche/build/mwpc_p3_noFW_jpsi_penning057_pilot_v1

This source is **byte-for-byte identical** to the frozen p3_noFW original
except for three provenance comments and ONE gas modification above.
The Python pilot performs exact source diff audit at runtime. No original
source/compiled binary/11k SHA data can be changed by the script.

It selects three existing **fully simulated and sealed** real Geant4
muon parents with 25..60 original Heed primaries from archived old-1000
chunk_0000. The old Penning-OFF seed CSV is read, NOT rerun.

Build/test:

~~~bash
cd ~/na6p/src/NA6PRoot-mwpc
git pull --ff-only
source ~/na6p/env.sh
source ~/na6p/install/garfieldpp/share/Garfield/setupGarfield.sh
cmake -S test/mwpc-avalanche -B test/mwpc-avalanche/build \
  -DGarfield_DIR="$HOME/na6p/install/garfieldpp/lib/cmake/Garfield"
cmake --build test/mwpc-avalanche/build \
  --target mwpc_p3_noFW_jpsi_penning057_pilot_v1 -j2

PY=test_runs/mwpc_stage_c/.venv_pi_slides/bin/python
SCRIPT=test/mwpc-stage-c/pi_p3_noFW_penning_gain_pilot_v1.py
"$PY" -m py_compile "$SCRIPT"
"$PY" "$SCRIPT" preflight --muons 3
set -o pipefail
"$PY" -u "$SCRIPT" run --muons 3 | tee ~/na6p/p3_noFW_penning057_gain_pilot.log
~~~

Paired exact Check in Python before producing comparison:
- Same true Geant4 muon sources, original 20-parent SHA baseline.
- Same original Heed primary electron identity, cluster, position and time.
- Same Garfield seed and electron starts for each individual pair.
- Uncapped exact collected_e and ion_births for original/penning-on.
- Computes all-seed mean, mean conditional positive and zero fraction.
- Difference ONLY the Penning gas setting.
- Any change in source/primaries, capped avalanche, or incomplete output
  prevents interpretation and leaves raw output for physical diagnosis.
- Important: pilot only 3 J/psi muon parents. Estimate effect size and
  computational cost; not yet a precision claim or experimental calibration.

Report (new directory refuses overwrite):
test_runs/mwpc_stage_c/pi_p3_noFW_jpsi_penning057_gain_pilot_v1/comparison.json

If run errors, DO NOT remove or overwrite partial output. Review:
test_runs/mwpc_stage_c/pi_p3_noFW_jpsi_penning057_gain_pilot_v1/physics.log

Do NOT repair the 11k Penning-off campaign by a trivial gain multiplier or
by recoding 333 million original ion births. Enabling Penning can change
per-avalanche gain distributions, ion birth positions, the fraction of
zero-electron outcomes, and therefore ion-birth asymmetry. A new validated
physics sample is required for calibrated Penning-on predictions.

Next if Penning improves gain but does not match prototype: hold gas,
fix geometry 3mm vs 4mm pitch and 3+3 vs 2+4 cathodes in separately
named C++ variant; scan actual HV using small matched seed samples;
compare with a genuine measured gain-vs-HV reference.
