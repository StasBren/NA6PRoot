# NA60+ MWPC HV scan — native Garfield++, full per-anode truth, ROOT-only analysis

## Experimental design

Five anode voltages: **1700, 1750, 1800, 1850, 1900 V**. The fixed chamber is
P3/no field wires, 3-mm anode pitch, 30-micron anode diameter, cathodes at
v=+/-3 mm, Ar/CO2 70/30, Penning transfer r=0.57, B=0.

Each voltage reuses exactly the **same 1,000 genuine Geant4 muon crossings**
from the prior wire-truth study, split 334 J/psi + 333 phi + 333 omega. The
frozen pseudo-random selection covers MS0-MS5 and is not enriched in incident
angle. An identical Heed RNG seed is used at all voltages. The runner compares
Heed primary and cluster CSV bytes across voltages and stops on mismatch.
It is a *paired-input* scan, NOT five independent muon samples.

The new HV-aware microscopic engine is an isolated CMake target:

`mwpc_p3_noFW_penning057_hvscan_v1`

and is called as

`mwpc_p3_noFW_penning057_hvscan_v1 --wiretruth-hv-p3-noFW CHUNK.csv PREFIX VOLTAGE_V`.

The validated reference 1800-V source and outputs are untouched. At every
voltage Garfield++ computes full microscopic multiplication and exact
collecting-wire electron endpoints plus an ion-birth census and moments by
anode. Each chunk is accepted only if ROOT validates all seed/muon/per-wire
charge/ion sums, checks n_Heed and detects no capped avalanches. Failures
are deliberately NOT silently censored.

## WSL / local launch

First clone or fetch the branch:

```bash
# New clean worktree, does not affect existing NA6PRoot worktrees.
git clone --depth 1 --branch feature/p3-hv-scan-20261116 \
    https://github.com/StasBren/NA6PRoot.git ~/na6p/p3-hv-scan
cd ~/na6p/p3-hv-scan

# Build and run 20 muons per voltage. Inspect no cap/failures before full run.
bash test/mwpc-stage-c/hv_scan_v1/run_hv_scan_v1.sh pilot

# Run remaining 49 chunks / voltage (resume-safe):
bash test/mwpc-stage-c/hv_scan_v1/run_hv_scan_v1.sh run

# ROOT-only computations and ATLAS-style figures:
bash test/mwpc-stage-c/hv_scan_v1/run_hv_scan_v1.sh analyse
```

The scan creates distinct directories under
`test_runs/mwpc_stage_c/hv_scan_wiretruth_1000_v1/voltages/V####/chunk_####/`.
There are 50 chunks of 20 muons for each voltage. Re-running the script skips
chunks with valid `ROOT_AUDIT_PASS.txt` seal; it never overwrites unsealed
partial output. `status` shows sealed chunk counts.

For a shorter first test, use
`NA60_HV_VOLTAGES="1750 1800 1850" bash .../run_hv_scan_v1.sh pilot`.
The full five-voltage scan is still the standard plan.

## PBS cluster option

Build and run the pilot in a cluster worktree first with `NA60_HV_OUT`
pointing at a *new* shared storage directory, not the previous 30k run.
Only after a clean pilot, submit production chunks 1..49 per voltage:

```bash
export NA60_HV_OUT=/storage/cgrp/stanis/p3_hv_scan_wiretruth_1000_v1
# Pilot and CMake build must already have completed in the same worktree.
qsub -J '0-244%12' -v "NA60_HV_OUT=$NA60_HV_OUT,NA60_HV_REPO=$PWD" \
  test/mwpc-stage-c/hv_scan_v1/pbs_hv_scan_chunk.sh
```

Confirm PBS array syntax, concurrency cap, resource requests, and fair-use
policy with your site's queue. The wrapper does NOT build 245 times or alter
completed runs. `NA60_HV_OUT` is mandatory on sites without shared home paths.

## ROOT analysis outputs

`reports/voltage_summary_ROOT.csv` with muon-bootstrap errors and
`reports/voltage_station_summary_ROOT.csv` plus ATLAS-style PDF/PNG and
`voltage_distributions.root`:

- Effective muon gain G_mu distributions and mean vs voltage.
- Single-avalanche gain fluctuation RMS(g)/mean(g) **including zeros**.
- Separate RMS(g|g>0)/mean(g|g>0), and zero-gain fraction.
- Muon effective-gain fluctuation RMS(G_mu)/mean(G_mu).
- Signed A_mu distributions, RMS(A_mu), mean |A_mu| and mean signed A_mu.
- Multi-anode events (second anode >5% of collected charge) and local |A_k|.
- Station-separated mean gain, gain fluctuations, signed asymmetry RMS, and
  multi-anode fractions at each voltage.

**Statistical unit for uncertainty:** the complete muon is bootstrap-resampled
with all its Heed avalanches together, not billions of correlated ion births.
400 bootstrap resamples; reported error is the standard deviation of each
bootstrap estimator (approximately one standard error).

## Is 1,000 muons / voltage enough?

Yes as an exploratory *overall* voltage scan. At 1800 V the observed
RMS(G_mu)~560 and mean~3265 imply an illustrative relative statistical
standard error of the mean of ~0.54% for N=1000 (before correlations/selection
weights). The multi-anode fraction near 19% has about 1.2 percentage-point
binomial SE (95% interval roughly +/-2.4 points). For individual stations,
large-angle tails, or subtle RMS/mean changes, 1000 total will not suffice.
If needed, expand to targeted kinematics and/or ~1000 muons per station,
not immediately to 10k per channel at every HV.

At higher V, the avalanche-size cap (200,000) can censor the gain tail. The
run rejects capped results rather than treating them as physical gain samples;
a high-voltage point that caps requires a new validated simulation protocol.
Voltage changes electron drift physics: **the fixed LUT500 built at 1800 V
must not simply be reused at 1700 or 1900 V without validation or regeneration**.

## What is still missing?

The ROOT/C++ integrity and analysis macros are supplied and script syntax
was checked. The C++ target has not been compiled or executed in this ChatGPT
environment (Garfield++ and ROOT are absent). Run the 20-muon pilot and share
its ROOT audit and compilation log to confirm compatibility before the full
campaign. Station counts for the original 30k (not this selected 1k) remain
handled by the separate 30k ROOT station audit in `na60-dice-mwpc-simulation`.
