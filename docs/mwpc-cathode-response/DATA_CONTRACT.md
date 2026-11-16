# Data contract for the next cathode-signal chat

## Coordinate and charge conventions

- Chamber-local u = transverse to anode wires; v = normal to cathode planes; w = along the wires.
- Real wire pitch 3 mm, physical wire labels are **global/station-aware**. Local nine-wire Garfield patch indices are only intermediate.
- \`collected_e[k]\`: number of electron endpoints that actually reached physical anode k. \`Q_anode,k = -e * collected_e[k]\` as signed physical electron charge (or positive gain magnitude e*N for plots). This is NOT a measured cathode-strip charge.
- \`ion_births[k]\`: number of ionization births associated with the nearest physical anode to the birth position. This **does not assert** that the associated avalanche electrons were collected by that same anode.
- \`A_k=(N_up,k-N_down,k)/(N_up,k+N_down,k)\` uses ion-birth positions; 0/0 is undefined, not 0.
- Time \`t_birth_ns\` and \`t_collected_ns\` have different physical meanings; preserve primary Heed time and absolute origins.

## Existing CSVs in the 1000-muon wire-truth run

- \`*_heed_primaries.csv\`: crossing_id, primary_index, cluster_index, u_rel_cm, v_abs_cm, w_rel_cm, t_ns, and raw Heed direction fields.
- \`*_wire_seed.csv\` / \`*_wire_muon.csv\`: crossing_id, station, theta_u_deg, theta_w_deg, primary_index/cluster_index, \`wire_index_global\`, \`wire_u_global_mm\`, \`collected_e\`, \`ion_births\`, \`ion_up\`, \`ion_down\`, \`ion_boundary\`, \`ion_top\`, \`ion_bottom\`, \`A_ion_hemi\`, \`A_ion_TB\`, \`mean_w_birth_global_mm\`, \`sigma_w_birth_mm\`, \`mean_t_birth_ns\`, \`sigma_t_birth_ns\`, \`mean_t_collected_ns\`, \`sigma_t_collected_ns\`. Do not assume listed means reconstruct complete probability distributions.
- \`*_ion_birth_sample.csv\`: sampled positive-ion creation points/times per seed, not the full billions of ion births; includes \`reservoir_weight\` approximately \`N_ion / N_sample\` for unbiased intra-seed pooling. Reservoir size <=128 in regular wire-truth reference. The **HV-scan fork sets the reservoir size to zero** in HV mode to save resources; it cannot be used as a waveform source on its own.
- \`*_seed.csv\` / \`*_crossing.csv\`: original Garfield++ gain/ion census/seed outcomes and original Geant4 crossing identities.

## Essential input for high-fidelity cathode currents

At minimum a verified weighted list \`[(u_j,v_j,w_j,t_j,weight_j,ion_species)]\` or a validated equivalent transfer-kernel representation. Also:
- Geometry and physical electric field at chosen HV; cathode-specific readout strip definitions, weighting field per electrode; chamber boundary and stereo conventions.
- Positive-ion transport model, calibrated ion mobility/ion mixture and correct time evolution; optional electron prompt current.
- Per-muon linkage, each electrode ID, time-window and sampling, clear physical units, environment/config checksums.
- Independent detailed Garfield++ event subset or high-stat direct ion current reference for waveform validation.

**Six per-wire summary moments \`{Q_k,A_k,mean_w,sigma_w,mean_t,sigma_t}\` are not sufficient to reconstruct arbitrary cathode waveforms.** A sum of signed A values can hide opposing clouds. Weighted reservoirs may lose fine details: test convergence explicitly.

## Minimum validation/crosschecks

1. Sum_k collected_e[k] equals collected electron endpoint count for each primary and muon.
2. Sum_k ion_births[k] equals callback ion births exactly (and zero/undefined cases are preserved).
3. Reconstruct muon A from exact ion up/down counts and compare to stored global A.
4. For ion sampling, compare macroion count/moments/angles, conditional on start/wire/gain; quantify reservoir variance vs size.
5. For strip waveform, check signed integrated induced charge via weighting-potential difference per ion and sum over strips vs expected boundaries; check dt/step-size/window convergence.
6. For complete event, separate true physical cathode strip response, artificial ideal infinite strips, and whole-plane cathode currents.
7. Physics-channel validation must not train on held-out phi/omega data; per-station and multi-anode (>5%) tests mandatory.

## Sources of truth

Use GitHub source and frozen manifest **before** relying on claims in earlier assistant messages. In particular inspect:
\`test/mwpc-avalanche/phaseB_native_muon_p3_noFW_penning057_wiretruth_v1.cxx\`,
\`test/mwpc-avalanche/phaseC_p3_noFW_two_cathode_ion_reference_v1.cxx\`,
\`test/mwpc-avalanche/phaseB_strip_weighting.cxx\`.
Check current local/farm runs with hashes, not just a copied status log.
