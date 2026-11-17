# NA60+ MWPC: truly uncapped microscopic avalanche-size scan (P3 noFW)

Requested 2026-11-17. Run **20 genuine Geant4 muons at EACH of 2000, 2050,
2100, 2150, 2200 V**, with identical frozen 1000-crossing selection prefix,
identical Heed seeds at each voltage, Penning Ar r=0.57, anode diameter 30 um,
anode pitch 3 mm, cathode planes +/-3 mm, Ar/CO2 70/30, 293.15 K and 760 Torr.

**This is an UNLIMITED-SIZE Garfield++ AvalancheMicroscopic simulation:**
- New, isolated executable:
  \`test/mwpc-avalanche/phaseB_native_muon_p3_noFW_penning057_uncapped_v1.cxx\`
- The dedicated executable explicitly calls
  \`avalanche.DisableAvalancheSizeLimit();\`, and has **no**
  \`avalanche.EnableAvalancheSizeLimit(...)\` call.
- That is different from the existing 2100-V and 1900–2050-V simulations,
  which used avalanche-size limits (typically 1,000,000).
- New CSV field \`electron_endpoints\` records the exact number of electron
  trajectory endpoints per individual Heed primary.
- Other per-primary output columns include \`ion_births\` (exact ionization
  callback count), \`garfield_ni\`, \`garfield_ne\`,
  \`collected_e\` (electron endpoints actually hitting a physical anode),
  \`attached\`, \`other\`, \`successful\`, \`capped\` (always 0 here).
- **Avalanche size** is not necessarily identical to collected-electron gain.
  The latter excludes attachment and endpoints that do not reach the anode.
- The independent small ion birth *reservoir* is disabled to reduce I/O,
  but the full callback ion-birth count and per-wire census remain.

## Important: no finite computer can guarantee that an unlimited avalanche completes

At large HV, unbounded microscopic multiplication can exhaust RAM/CPU or
become physically invalid through collective space charge, discharge and
streamer physics NOT modeled here. This calculation is therefore exploratory,
not a validated description of a stable operating MWPC.

To protect the WSL laptop:
- Each **ONE muon is a separate process**, with a **1200-s wall-time guard**
  and **4-GB virtual address-space guard** by default.
- The process-level guards **never artificially cap the avalanche**. A timeout,
  bad allocation, nonzero Garfield error, or a missing/partial output is
  marked **UNFINISHED** and is NOT included as a measured finite avalanche.
- An incomplete result is quarantined and no ROOT histogram is produced
  unless *all five* voltage ensembles completed successfully.
- Reusing a completed muon requires a SHA256 validation; if changing process
  limits or the executable, use a **new output directory**.
- To allow a larger process virtual-memory budget on a machine with enough
  physical RAM: \`HV_UNCAPPED_VMEM_GB=8\`. For longer per-muon execution:
  \`HV_UNCAPPED_TIMEOUT_S=3600\`. These require a fresh \`HV_UNCAPPED_OUT\`.

A runtime/virtual-memory guard is **not a physics size limit**; it can only
abort and **censor** an unfinished observation. It is incorrect to quote a
finite gain, a gain tail or a mean for aborted/censored events.

## WSL laptop commands

\`\`\`bash
git clone --depth 1 --branch feature/p3-uncapped-hv-2000-2200-20261117 \
  https://github.com/StasBren/NA6PRoot.git ~/na6p/p3-uncapped-hv

# One muon per voltage first; test memory/runtime with 2200 V first.
bash ~/na6p/p3-uncapped-hv/test/mwpc-stage-c/hv-uncapped-2000-2200/run_uncapped_local.sh 1

# After checking runtime, output and memory, extend to 5.
bash ~/na6p/p3-uncapped-hv/test/mwpc-stage-c/hv-uncapped-2000-2200/run_uncapped_local.sh 5

# Final: 20 per voltage, 100 full muons in total.
bash ~/na6p/p3-uncapped-hv/test/mwpc-stage-c/hv-uncapped-2000-2200/run_uncapped_local.sh 20
\`\`\`

Do not run this simultaneously with another heavy Garfield++ campaign on a
memory-constrained laptop.

For a fresh run with more memory (ONLY if your laptop can actually afford it):

\`\`\`bash
HV_UNCAPPED_VMEM_GB=8 HV_UNCAPPED_TIMEOUT_S=3600 \
HV_UNCAPPED_OUT="$HOME/na6p/uncapped_hv_8GB" \
bash ~/na6p/p3-uncapped-hv/test/mwpc-stage-c/hv-uncapped-2000-2200/run_uncapped_local.sh 1
\`\`\`

## Results: ROOT-only physics and ATLAS-style plots

\`~/na6p/p3-uncapped-hv/test_runs/mwpc_stage_c/p3_uncapped_2000_2200_v1/analysis/\`

- \`uncapped_avalanche_summary.csv\`: mean collected gain including zeros,
  RMS/mean all and positive only, mean ion births and RMS/mean, mean Garfield
  ion and electron counts, maximum collected gain, maximum ion births,
  maximum electron endpoints and tail fractions >1e6 and >1e7 ion births.
  All uncertainty estimates bootstrap COMPLETE muons.
- \`uncapped_mean_gain_vs_voltage.pdf/png\`
- \`uncapped_mean_ion_births_vs_voltage.pdf/png\`
- \`uncapped_gain_fluctuations_vs_voltage.pdf/png\`
- \`uncapped_collected_gain_distributions.pdf/png\`
- \`uncapped_ion_birth_distributions.pdf/png\`

The 1- and 5-muon stages are technical probes, not high-precision fluctuation
measurements. Exact gas gain may differ from experimental measurements even with
the cap disabled, because of physics-model and detector-design differences.
The **working-point decision of 2100 V remains unchanged** unless the user
explicitly revisits it. This scan is a numerical cross-check of the avalanche
cutoff, not an automatic recommendation to operate at 2200 V.

**Status**: sources committed on GitHub. Build/execution on the user's laptop
must be confirmed by actual logs; no simulated numerical results are claimed.
