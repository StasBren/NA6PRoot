# High-voltage and anode-diameter Garfield++ test

**Separate from the currently running 1900–2050 V scan.**

## Exact requested study

Same Geant4 crossings (from the 1000-muon frozen selection) and seeded TrackHeed at every point:

| Anode voltage | Diameter | Muons | Purpose |
|---:|---:|---:|---|
| 2000 V | 30 um | 20 | reference |
| 2050 V | 30 um | 20 | gain / RMS/mean voltage dependence |
| 2100 V | 30 um | 20 | upper-voltage study, check caps and physical validity |
| 2000 V | 50 um | 5 | paired *first five* muons vs 30 um |

**Pilot** first runs five muons at all four configurations. **Full** resumes these and extends 30 um points to 20 muons (50 um stays five). The 30 vs 50 um comparison uses the identical first five muons at 2000 V, *not* 20 vs 5.

New source: \`test/mwpc-avalanche/phaseB_native_muon_p3_noFW_penning057_hvdiam_v1.cxx\`.
New executable: \`mwpc_p3_noFW_penning057_hvdiam_v1\`.

- True Garfield anode wire diameter is set to 0.003 cm for 30 um or 0.005 cm for 50 um via a runtime argument.
- The electron-endpoint *collection radius* changes consistently with the wire diameter.
- Wire pitch 3 mm, cathode separation 6 mm, no field wires, Ar/CO2 70:30, T=293.15 K, pressure 760 Torr, Penning r=0.57, B=0 are unchanged.
- Each primary Heed electron gets a complete Garfield::AvalancheMicroscopic simulation; original per-wire electron/ion census retained.
- The avalanche cap defaults to 1,000,000; if any cap or failure occurs, stop and **do not use the truncated gain distribution**.
- High-gain results are without space charge or discharge physics, so 2100 V predictions are not operating voltage recommendations.

## Run on your WSL laptop

From your existing terminal:

\`\`\`bash
git clone --depth 1 --branch feature/p3-hv-wire50-check-20261116 \
  https://github.com/StasBren/NA6PRoot.git \
  ~/na6p/p3-hv-wire50

bash ~/na6p/p3-hv-wire50/test/mwpc-stage-c/hv-diameter-study/run_hvdiam_local.sh pilot
\`\`\`

Once pilot 2100 V performance, cap rates and logs look reasonable:

\`\`\`bash
bash ~/na6p/p3-hv-wire50/test/mwpc-stage-c/hv-diameter-study/run_hvdiam_local.sh full
\`\`\`

Alternatively use \`git pull --ff-only\` inside \`~/na6p/p3-hv-wire50\` when later fixes land.

Results remain under
\`~/na6p/p3-hv-wire50/test_runs/mwpc_stage_c/hvdiam_study_20261116_v1/\`,
**not** the original voltage-scan directory. Each five-muon chunk has a checksum file and seal; completed chunks are verified before being reused. A corrupt or capped chunk stops the campaign. A configuration change requires a **new output directory** via \`HV_DIAM_OUT=/separate/path\`.

The script builds a separate CMake target (\`build-hvdiam\`) and then executes ROOT/C++ macro
\`Analyze_HVdiam_ROOT.C\`.

ROOT outputs in \`analysis/\`:
- \`gain_summary.csv\`: per-primary mean gain including zeros, positive-only mean, zero fraction, RMS/mean on both definitions, and muon-effective gain; whole-muon bootstrap uncertainties.
- \`mean_gain_vs_voltage_30um.pdf/png\`
- \`gain_fluctuations_vs_voltage_30um.pdf/png\`
- \`positive_gain_distributions_30um.pdf/png\`
- \`diameter_2000V_5muon_comparison.pdf/png\` (paired same five inputs).

## Interpretation

20 muons produce approximately 1,000 individual Heed avalanches at each 30 um voltage. It is a *preliminary* voltage scan, not a precision result, especially if the microscopic gain distribution has long tails. Five muons at 50 um are only a geometry-sensitivity pilot; do not extrapolate precise relative changes or stable operation.

Use ROOT only (ATLAS-like plots marked NA60+ Simulation), without Gaussian fitting; for signed asymmetry with mean near zero use RMS(A), mean|A| and the full distribution instead of RMS(A)/mean(A).

Compare with experimental mean gain only after matching wire diameter, chamber geometry, operating gas and temperature/pressure. The Penning factor has NOT been retuned to fit experiment.
