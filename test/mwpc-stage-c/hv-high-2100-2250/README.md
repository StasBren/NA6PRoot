# P3 noFW high-voltage exploratory scan: 2100, 2150, 2200, 2250 V

**Purpose:** test the voltage dependence of *collected-electron gas gain* and relative gain fluctuations in full Garfield++ `AvalancheMicroscopic`. This is **not a stable-operation recommendation**. High-field microscopic avalanches omit space-charge, streamers, discharge and related collective effects.

**Physics setup:** 30 um anode diameter; 3 mm pitch; cathodes +/-3 mm; no field wires; Ar/CO2 70/30; 293.15 K, 760 Torr; Penning Ar r=0.57, B=0. Other configurations frozen. The new high-HV source differs from the earlier HV source by its allowed upper-voltage bound only and a safety comment. It is a separate CMake target, output tree and Git branch.

## WSL commands

```bash
git clone --depth 1 --branch feature/p3-hv-high-2100-2250-20261116 \
  https://github.com/StasBren/NA6PRoot.git ~/na6p/p3-hv-high

# First run exactly 1 muon at each voltage, starting with the most demanding 2250 V.
bash ~/na6p/p3-hv-high/test/mwpc-stage-c/hv-high-2100-2250/run_local_high_hv.sh 1

# Only after confirming runtime and NO capped avalanches:
bash ~/na6p/p3-hv-high/test/mwpc-stage-c/hv-high-2100-2250/run_local_high_hv.sh 5

# Expand to 20 muons at each voltage (80 full muon simulations total)
bash ~/na6p/p3-hv-high/test/mwpc-stage-c/hv-high-2100-2250/run_local_high_hv.sh 20
```

The same frozen selection of 1000 Geant4 crossings used in the 1k wire-truth run is reused in the same order for every voltage, including the same TrackHeed and Garfield seed logic. The 1/5/20 modes reuse completed one-muon chunks with SHA256 audit, not recomputation. Do NOT run this at the same time as the old 1900-2050 V scan on a resource-limited laptop.

**Before running:** ROOT and Garfield++ installations from prior P3 noFW study must be available in WSL. The runner builds the independent CMake target `mwpc_p3_noFW_penning057_hvhigh_v1`. It does not touch old build directories or scan outputs.

**Safeguards:** by default the run sets `EnableAvalancheSizeLimit(1000000)` and hard-stops when any seed is flagged capped or failed. Per-muon runtime is capped at 1200 seconds by `timeout` (the run will fail, leave an incomplete chunk quarantined on retry, and will not produce fake physics plots). To increase run timeout use a **new output directory** with e.g. `HV_HIGH_TIMEOUT_SECS=3600 HV_HIGH_OUTPUT_DIR=... bash ...`, because the per-run configuration is sealed. Memory consumption and CPU cost can rise rapidly at 2250 V; extending the cap is not necessarily physically meaningful.

Results under `~/na6p/p3-hv-high/test_runs/mwpc_stage_c/p3_hv_high_2100_2250_v1/`:
- one-muon chunks with per-wire collecting-electron and ion-birth census, log files, checksums;
- `analysis/gain_voltage_summary.csv`: mean collected gain, zero fractions, all-positive RMS/mean and per-muon effective gain; uncertainties from resampling *muons* (not treating all avalanche electrons as independent);
- `analysis/mean_gain_vs_HV.pdf/png`;
- `analysis/gain_RMS_over_mean_vs_HV.pdf/png`;
- `analysis/positive_gain_distributions.pdf/png`;
- `analysis/signed_asymmetry_distributions.pdf/png`.

Plots from 1 or 5 muons are for technical pilots ONLY and must not be interpreted as accurate signed-asymmetry or gain-fluctuation distributions.

## Does a diffusion step suppress our electron gain?

No macroscopic diffusion step is configured: this experiment uses `Garfield::AvalancheMicroscopic`, which follows microscopic electron collision histories from `MediumMagboltz`. The `AvalancheMC` time/distance/collision-step parameters are a distinct approximation, not used here. If the gain disagrees with experiment, check actual fields, wire diameter and pitch, collection vs ionization definitions, attachment, gas pressure, effective Penning model, avalanche capping, and missing collective effects before inferring inadequate stepping.

## Status

GitHub files and source changes were verified. The assistant could not compile Garfield++ or ROOT on the user's WSL installation; code is **not claimed to have run** until the user supplies logs.
