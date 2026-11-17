# Copy-paste prompt for a separate ChatGPT project chat

We are continuing the NA60+ / DiCE MWPC **cathode-signal response simulation**. Our earlier conversations are long; **do not assume you remember them accurately**.

**First action (mandatory):** access GitHub repository
\`StasBren/NA6PRoot\`, branch
\`feature/mwpc-cathode-response-handoff-20261116\`, folder
\`docs/mwpc-cathode-response/\`.
Read, not merely summarize from file names:
\`README.md\`, \`REQUIREMENTS_AND_STATUS.md\`,
\`DATA_CONTRACT.md\`, \`SIGNAL_IMPLEMENTATION_PLAN.md\`, and \`AGENTS.md\`.
Also inspect the *actual* C++ sources linked there:
\`test/mwpc-avalanche/phaseB_strip_weighting.cxx\`,
\`phaseC_p3_noFW_two_cathode_ion_reference_v1.cxx\`,
\`phaseB_real_avalanche_signal.cxx\`,
\`phaseB_native_muon_p3_noFW_penning057_wiretruth_v1.cxx\`
and the HV-scan source and runner.

Scientific goal: from real Geant4/TrackHeed muons plus verified microscopic Garfield ion birth data, compute **physical induced current on individual upper/lower cathode strips**, then electronics-shaped response and reconstruction observables. Maintain real collecting-wire IDs and local ion-cloud correlations; global muon asymmetry \`A_mu\` must not be used as a proxy for actual local per-wire ion clouds or cathode charge.

Constraints:
1. **Physics analysis and every plot must use C++/ROOT and ROOT ATLAS-like visual style**. No Python calculations/plots. Present plots honestly as \`NA60+ Simulation\`.
2. Do not modify historical frozen 1800 V reference, 1000 wire-truth data or ongoing HV scan. Create isolated branch, target and output folders. Work incrementally through GitHub; give me simple WSL bash commands. Never claim to have executed Garfield++ on my laptop or HPC unless logs actually prove it.
3. Current baseline: P3 noFW, anode pitch 3 mm, 30 um wire diameter, cathodes at +/-3 mm, Ar/CO2 70:30, Penning r=0.57, B=0. **2100 V is the user-selected working point** (reported gain ~80,000, not yet independently audited); use branch `feature/p3-2100V-only-20261117` for the dedicated working-point cross-check. Other voltage scans are historical. Keep HV configurable, check avalanche cap and missing space-charge effects.
4. From Sasha's confirmed guidance: width ~1.7–2 mm and \`tan(alpha)=0.1\` to anode wire are **provisional strip design values**, not final specifications; detector-response code belongs with our group and must remain adjustable. Distinguish the old perpendicular-to-wire strip weighting prototype from the intended stereo geometry.
5. User-relayed Sasha tasks: gain vs HV; RMS/mean (not Gaussian fit); full distributions; channel-dependent azimuth differences; station-by-station checks; slide-9 asymmetry vs voltage; ROOT-only ATLAS-like plots; 30k muon crossing selection audit.
6. The old 30k Garfield ensemble (10k J/psi + 10k phi + 10k omega) lacked collecting-wire IDs. The completed independent 1000 wire truth includes them; 192/1000 events have a second anode with >5% charge, and 143/192 have opposite ion-birth A on the two leading wires. The *new* 30k per-wire campaign status must be verified.
7. The 1000 wire-truth archive's \`ion_birth_sample.csv\` stores <=128 sampled births per seed; HV scan stores no sampled birth points. Per-wire \`Q_k,A_k,mean_w,sigma_w,mean_t,sigma_t\` are NOT sufficient to produce arbitrary current waveforms. Plan a small full-ion-birth direct reference and validate weighted macroions / kernels before making fast signal claims.
8. Existing Phase C two-cathode-ion code calculates current on WHOLE cathode planes, NOT physical stereo strips, and its ion mobility approximation is explicitly provisional. Distinguish physical electric field from strip weighting field, and verify Ramo integrated charge, signs, times, angular/longitudinal distributions, convergence and station/edge cases.
9. Keep Fast4 as detailed surrogate; Fast5 is postponed until per-anode and signal validation. Consider a conditional *whole-muon induced-current template library* after that; do not claim another order-of-magnitude speedup without measured same-hardware benchmarks.

**Immediate task:** read and critically audit existing Phase B/C signal code and implement the FIRST validated milestone: one positive ion inducing current on a *properly defined cathode strip* at the chosen/reference voltage, with a ROOT-only waveform plot, numerical convergence and Ramo charge-integral check. If 2100 V gain/cap results are not yet validated, keep the source HV configurable and use 1800 V only as a reproducible numerical control. Do not call the 2100 V gain or detector stability experimentally confirmed. Explain what is known, what is assumed, and what must be measured or calibrated (e.g. Ar/CO2 ion mobility, stereo geometry/electronics). Then provide code via GitHub and a minimal WSL command to build/run. Continue stepwise; don't jump straight to full-muon waveform claims.

Before generating code, identify concrete risks/bugs in the current sources and present a short prioritized plan. After changes, cite Git commits/paths and actual completed local checks. Do not ask me to re-explain context already provided in the handoff documents.
