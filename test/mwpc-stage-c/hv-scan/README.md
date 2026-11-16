# NA60+ P3 noFW voltage scan (ROOT-only analysis)

This **separate** branch runs microscopic Garfield++ avalanches at **1900, 1950,
2000 and 2050 V**, using the **same 1,000 frozen Geant4 crossings and the same
Heed seeds at each voltage**. Inputs are the 334 J/psi, 333 phi and 333 omega
crossings used in the prior 1000-muon per-anode reference (all six stations).

## Physics

- Reference: P3 noFW; anode pitch 3 mm, diameter 30 um; cathodes at +/-3 mm.
- Ar/CO2 = 70/30, Penning Ar r=0.57, B=0.
- Only anode high voltage is changed for this scan.
- Every Heed primary gets a full Garfield++ AvalancheMicroscopic call.
- Real collecting-wire IDs, exact per-wire ion birth census/moments, seed gain,
  global and local asymmetries are preserved.
- The independent 128-ion-position reservoir is **disabled in HV mode** to
  avoid enormous I/O; it is not needed for the voltage-dependent observables.
- Avalanche cap is **1,000,000** by default in HV mode (old reference uses
  200,000). The runner aborts if even one cap is reached. NEVER interpret a
  capped run as a true gain distribution.
- Space charge, streamers and discharge are NOT simulated. Large-gain voltage
  results are microscopic-model predictions, not operating voltage advice.

## Laptop (WSL)

Clone in a separate directory so previous runs are untouched:

```bash
git clone --depth 1 \
  --branch feature/p3-hv-scan-20261116 \
  https://github.com/StasBren/NA6PRoot.git \
  ~/na6p/p3-hv-scan

bash ~/na6p/p3-hv-scan/test/mwpc-stage-c/hv-scan/run_local_hv_scan.sh 5
```

The runner loads `~/na6p/env.sh` and
`~/na6p/install/garfieldpp/share/Garfield/setupGarfield.sh` if available,
then builds a **separate** `build-hvscan` CMake directory. It requires ROOT,
CMake and the Garfield++ installation already used for the prior wire-truth run.

This first command runs **5 muons at each of four voltages**. If all pass and
runtime/memory/capped checks look safe, grow to 20, 100 and finally 1000 muons
per voltage using the same command with a new integer:

```bash
bash ~/na6p/p3-hv-scan/test/mwpc-stage-c/hv-scan/run_local_hv_scan.sh 20
bash ~/na6p/p3-hv-scan/test/mwpc-stage-c/hv-scan/run_local_hv_scan.sh 100
# Only after measuring realistic total runtime:
bash ~/na6p/p3-hv-scan/test/mwpc-stage-c/hv-scan/run_local_hv_scan.sh 1000
```

Previously sealed chunks are checksum-verified and reused. Incomplete chunks
are quarantined. Do not change `HV_CHUNK_SIZE`, `HV_AVALANCHE_CAP` or checkout
a different source commit mid-campaign. The default chunk size is **5 muons**.

To change the cap or start a distinct run, select a NEW output folder, e.g.:

```bash
HV_AVALANCHE_CAP=2000000 \
HV_OUTPUT_DIR="$HOME/na6p/p3-hv-scan/test_runs/mwpc_stage_c/p3_hv_scan_cap2m" \
bash ~/na6p/p3-hv-scan/test/mwpc-stage-c/hv-scan/run_local_hv_scan.sh 5
```

Larger caps consume more CPU and memory; they do NOT make large-gain physics
more realistic. A frequently capped 2050 V point may be impractical on a
laptop. Consider the farm if the 5-muon pilot is slow.

## Outputs

`test_runs/mwpc_stage_c/p3_hv_scan_v1/HV_<voltage>/chunk_XXXX/`

Each chunk contains `input.csv`, `physics.log`,
`microscopic_crossing.csv`, `microscopic_seed.csv`,
`microscopic_wire_muon.csv`, `microscopic_wire_seed.csv`, timing CSVs and
`complete.txt` + `SHA256SUMS`.

ROOT/C++ analysis outputs:
`test_runs/mwpc_stage_c/p3_hv_scan_v1/analysis/`

- `gain_voltage_summary.csv`
- `mean_gain_vs_HV.pdf/png`
- `gain_RMS_over_mean_vs_HV.pdf/png`
- `positive_gain_distributions.pdf/png`
- `signed_asymmetry_distributions.pdf/png`

The ROOT bootstrap resamples **whole muons**, preserving their internally
correlated avalanches, so the reported error bars are not based on the false
assumption that ~53,000 avalanches are independent across muons.

The plots show **RMS/mean**, not Gaussian fit parameters. It is evaluated for
all primary seeds, positive-response seeds only, and effective muon gain.
Signed asymmetry is plotted as a distribution; RMS(A)/mean(A) is not meaningful
because the signed mean is near zero.

Important: the 5-muon pilot provides only a smoke-test; do not use its error
bars to claim small voltage effects. The 1000-muon study is the first
statistically meaningful comparison.

## Status and monitoring

```bash
tail -f ~/na6p/p3-hv-scan/test_runs/mwpc_stage_c/p3_hv_scan_v1/HV_2050/chunk_0000/physics.log
```

If 2050 V is extremely slow or a capped avalanche occurs, preserve logs and
stop. Do not turn off cap checks merely to complete the scan.
