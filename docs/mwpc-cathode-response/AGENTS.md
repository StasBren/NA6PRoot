# Instructions for AI coding assistants working on MWPC cathode signal

This directory is the project handoff and scientific source of truth. Always read README, REQUIREMENTS_AND_STATUS, DATA_CONTRACT and SIGNAL_IMPLEMENTATION_PLAN before modifications.

- **Never fabricate simulation or validation results.** Distinguish code written, static checks, local compile, and actual Garfield++ physics runs.
- Use **C++/ROOT for every physics/statistical calculation and ROOT for every graph**. ATLAS-like fonts/layout, honest NA60+ Simulation label; no Python analysis or plotting.
- User has WSL Garfield++ and ROOT; AI tool execution may not. Work through GitHub with isolated branches and offer one clear bash command. Respect user's running high-voltage scan.
- Preserve immutable 1800 V reference and frozen 1000/30k datasets. No changes to their binaries, result folders or identities. Add new targets/config.
- Always separate per-Heed-primary, per-physical-collecting-wire, nearest-wire ion birth, per-muon, and per-cathode-strip observables.
- Physical anode electron collection is assigned by Garfield electron endpoint. Nearest-wire ion birth association is **not** a collecting endpoint. Signed A is not induced cathode charge.
- Use conditional ion-birth distribution and weighting fields for signal. Do not invent per-anode moments from A_mu or use six moments as a surrogate for full waveform.
- Differentiate drift-only LUT500 (30x30 cells x 500 trials, historical reference voltage) from full AvalancheMicroscopic. LUT500 must not be reused at a changed HV without validation or recomputation.
- The user selected **2100 V as the working development point** on 2026-11-17 (user-reported approximate gain 80,000). Do not treat it as physically validated until uncapped runs and missing space-charge effects are evaluated. The active single-voltage verification branch is `feature/p3-2100V-only-20261117`. Do not automatically extend to 2150/2200/2250 V.
- Strip width approximately 1.7–2 mm and tan(alpha)=0.1 to wire are design guidance, subject to optimization. Old Phase B strips run perpendicular to wire; don't silently adopt them as stereo readout.
- Existing Phase C two-cathode-ion signal is whole-cathode and has provisional ion mobility; don't call it individual-strip truth.
- Unit tests: Ramo current-integral vs weighting potential, field geometry, polarity, endpoint and time-window, weighted-ion convergence, charge/birth conservation, station and multi-wire cases.
- Document assumptions, inputs, versions, citations, output checksums, held-out samples, numerical and statistical errors and measured wall times.
- Prefer detailed Fast4 validation first; Fast5-L whole-muon template strategy is a later hypothesis, not a proven deployment.
