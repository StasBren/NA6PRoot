# Six-geometry 50-trial Garfield LUT (separate study)

## Scientific definition

Six **new** independent single-electron Garfield response grids:
`anode pitch = 3.0, 3.5, 4.0 mm`, each with and without a **75 μm** grounded midpoint field wire.

Physics held constant relative to the frozen 20-trial pilot: Ar/CO₂ 70:30 at 293.15 K and 760 Torr; +1800 V anodes, 30 μm anode diameter; cathode planes at -3 and +3 mm, 0 T; `du=0.1 mm`, `dv=0.2 mm`, 30 bins in v. Field-wire diameter is a study assumption, not a validated NA60+ engineering value.

For every (pitch, field-wire) configuration, **50 Garfield microscopic drift trials per position box**. Grid sizes: 900, 1050, 1200 boxes for 3, 3.5, 4 mm pitch respectively, giving **315,000 distinct Garfield drift trials** across all six geometries; each output is analysed conditionally on anode collection. These are not 315,000 Geant4 muons or 315,000 independent experimental detection events.

## Safety and isolation

- Retain and never overwrite the original 20-trial maps at `test_runs/mwpc_stage_c/pi_pitch_fieldwire_pilot_v1/pilot/`.
- New maps live solely at `test_runs/mwpc_stage_c/pi_pitch_fieldwire_50_v1/pilot/`.
- A dedicated RNG base seed `260050` and immutable `frozen_config.json` per configuration record the new sampling conditions and source/binary SHA256.
- The original 6M parent Geant4 + Heed folding pipeline **continues on the original 20-trial maps**. This run does not alter its frozen physics or existing results.
- The Garfield jobs run with `nice +15` priority so the existing Geant4/Heed CPU work has preference; two workloads nevertheless share machine resources.
- The six-configuration run is checkpointed one u-column at a time; 210 complete columns total, 1500 Garfield electrons each. On restart, previously saved columns are checked before being skipped. Do not delete or edit checkpoint files while the job is running.

## Run

```bash
cd ~/na6p/src/NA6PRoot-mwpc
git pull --ff-only

bash test/mwpc-stage-c/run_pi_pitch_fieldwire_50.sh check
bash test/mwpc-stage-c/run_pi_pitch_fieldwire_50.sh start

# Any time:
bash test/mwpc-stage-c/run_pi_pitch_fieldwire_50.sh status
tail -f test_runs/mwpc_stage_c/pi_pitch_fieldwire_50_v1/runner.log

# After successful completion:
bash test/mwpc-stage-c/run_pi_pitch_fieldwire_50.sh verify
bash test/mwpc-stage-c/run_pi_drive_quick_reports.sh copy
```

The Drive command above uploads only summaries; the large raw Garfield tables remain local until deliberately archived.

## Interpretation

The new LUT is an independent **response lookup ensemble**. Merely generating it does not update the preliminary plots or the frozen 6M-fold samples. To create new plots for the PI, the **same saved Geant4/Heed inputs** must be re-folded with the 50-trial maps into an isolated analysis namespace, then compared with 20-trial maps. Do not directly change the 6M original fold pipeline or represent reuse of single-electron LUT draws as independently simulated microscopic trajectories.

Comparison to 20 trials should include cell-level uncertainty/variance, station mean drift, collection fraction, and **full** conditional electron incidence and azimuth distributions vs signed muon entrance angles. Source sample size and independence of Garfield trials must remain separately stated.
