# NA60+ MWPC — cathode-signal response: cross-chat handoff

**Prepared 2026-11-16.** This directory is the durable, version-controlled source of truth for starting the cathode-signal stage in another conversation. It is a *plan*, not a claim that strip signals are already implemented or validated.

## Read these files in order

1. [REQUIREMENTS_AND_STATUS.md](REQUIREMENTS_AND_STATUS.md): which instructions come from group guidance/user, which facts are confirmed, and unresolved decisions.
2. [DATA_CONTRACT.md](DATA_CONTRACT.md): what Garfield++ really records; what is missing for waveform generation.
3. [SIGNAL_IMPLEMENTATION_PLAN.md](SIGNAL_IMPLEMENTATION_PLAN.md): prioritized physics/engineering implementation, tests and acceptance gates.
4. [NEXT_CHAT_PROMPT.md](NEXT_CHAT_PROMPT.md): complete prompt for a new ChatGPT project chat.
5. [AGENTS.md](AGENTS.md): standing instructions for AI assistants doing work in this directory.

**Do not treat other chats as automatically synchronized.** A new assistant must fetch and read these *actual files* before coding or offering a research decision.

## Relevant verified code

- 1k per-wire Garfield++ source:
  \`test/mwpc-avalanche/phaseB_native_muon_p3_noFW_penning057_wiretruth_v1.cxx\`
- Voltage scan, in progress on laptop at 1900, 1950, 2000, 2050 V:
  \`test/mwpc-stage-c/hv-scan/\` and \`test/mwpc-avalanche/phaseB_native_muon_p3_noFW_penning057_hvscan_v1.cxx\`
- Existing strip-weighting field prototype:
  \`test/mwpc-avalanche/phaseB_strip_weighting.cxx\`
- Existing full-plane two-cathode ion-only signal prototype:
  \`test/mwpc-avalanche/phaseC_p3_noFW_two_cathode_ion_reference_v1.cxx\`
- Exploratory paired electron/ion and realistic avalanche signal:
  \`test/mwpc-avalanche/phaseB_real_avalanche_signal.cxx\`, \`phaseB_electron_ion_cloud_signal.cxx\`
- Separate ROOT-only station/gain audit branch in repository
  \`StasBren/na60-dice-mwpc-simulation\`: \`feature/root-only-mwpc-analysis-20261115\`.
- Finished preliminary 1k reference has true wire IDs, and archived 30k reference did not. New 30k wire-truth *may be running separately*; verify run manifests rather than asserting completion.

## Next scientific milestone

Select an operating voltage *after* completing/assessing the voltage scan and cap diagnostics. Demonstrate **one honest, physically consistent cathode-strip current** from a verified ion-birth input and weighting field, checking charge conservation and numerical convergence. Then extend to both cathode sides, real geometry, electron contribution, electronics and muon-level validation. Do not jump straight from six per-wire moments to a claim of realistic waveforms.

## Reference geometry (baseline, NOT final hardware)

P3 noFW, 3 mm anode pitch, 30 um wire diameter, cathodes at v= +/-3 mm, Ar/CO2 70:30, Penning Ar r=0.57, B=0, nominal historical 1800 V. Actual 1900–2050 V operating point undecided. Strip width around 1.7–2 mm and tan(alpha)=0.1 to wire are *provisional design guidance*. Station-specific mechanics and stereo orientation must be configurable.
