# Decision: reuse immutable LUT20/LUT500 for Stage B1 avalanche pilot — 29 Oct 2026

**Conclusion from code-level inspection:** YES, reuse all four existing ±3 mm p3/p4 noFW/FW75 single-electron Garfield *drift* grids as **start-point trial indexes** and provenance. NO, do not treat their terminal anode endpoints as the physical starting state of an avalanche or assume a new AvalancheElectron trajectory is an exact continuation of the original drift.

## Sources reviewed in branch feature/mwpc-stage-c

- `test/mwpc-stage-c/c_transport_box_map.cxx`, primary CSV writer and `Garfield::AvalancheMicroscopic::DriftElectron` loop
- `test/mwpc-stage-c/pi_pitch_fieldwire_pilot.py`, 20-per-cell configuration and immutable frozen config
- `test/mwpc-stage-c/pi_six_geometry_native_timing.py::load_map/response_for_births/birth_sample`
- `test/mwpc-stage-c/pi_lut_20_vs_500_convergence.py` and `pi_lut500_muon_crossing_asymmetry.py`
- `test/mwpc-avalanche/phaseB_native_muon_3p3_avalanche.cxx` and `test/mwpc-stage-c/pi_benchmark_native_muon_b1.py`

### Actual frozen raw LUT schema

`garfield_electrons.csv` has `iu,iv,box_id,u_center_mm,v_center_mm,u_start_mm,v_start_mm,w_start_mm,sample,random_seed,drift_ok,collected,status,wire_index,wire_u_mm,u_end_mm,v_end_mm,w_end_mm,drift_time_ns,phi_hit_deg,cos_phi_hit,sin_phi_hit,hit_side_u,incidence_normal_deg,last_step_u_mm,last_step_v_mm,last_step_w_mm,last_step_energy_eV,last_dir_u,last_dir_v,last_dir_w,field_wire_collected,field_wire_index,field_wire_u_mm`.

The Python high-level `load_map` **discards** most of these fields and retains only `(collected, wire_index, drift_time_ns, incidence_normal_deg, phi_hit_deg)`. For an avalanche add-on, read the ORIGINAL raw CSV, not only the five-column NumPy `pool`. The raw CSV has exact **randomly sampled within-cell original starts** and frozen per-trial identities; therefore no need to regenerate a new *start grid*.

Physics identity is compatible across original 20 and 500 series: Ar/CO2 70:30, T=293.15K, P=760 Torr, 30 µm anode +1800V, ±3mm cathodes, grounded 75 µm midpoint FW if enabled, B=0. Pitches p=3mm and p=4mm. `du=0.1mm,dv=0.2mm` with 30 v bins and 30 or 40 u bins per pitch. Raw C++ generator uses `DriftElectron`, not `AvalancheElectron`, with a starting low-energy electron at 0.1 eV and direction (0,0,0). Geometry model is a finite 9-anode analytic patch, not full chamber edges.

| pitch | cells per geometry | raw LUT20 trials per geometry | raw LUT500 trials per geometry |
| --- | --- | --- | --- |
| 3 mm | 900 | 18,000 | 450,000 |
| 4 mm | 1,200 | 24,000 | 600,000 |
| **TOTAL FOUR geometries** | **4,200** | **84,000** | **2,100,000** |

These counts come from code and frozen layouts; the full WSL data files have **not** been inspected directly from ChatGPT. The read-only `pi_audit_avalanche_lut_reuse.py` is provided for a true, local audit of all eight source CSVs and SHA identities before physics jobs.

## Physics issue: what can and cannot be attached

A drift electron reaching the ANODE is **already past** the multiplication region; `last_step_*` is the last microscopic step near termination, not the state before Townsend multiplication. Launching `AvalancheElectron` from `u_end,v_end` or last-step point would omit the avalanche gain that developed upstream of that point.

**A: VALID AND SIMPLE first prototype:** for a selected frozen LUT key `(geometry,old_lut_csv_sha256,box_id,sample)`, take the **original stored exact initial coordinate** `(u_start,v_start,w_start)` and start a NEW `AvalancheElectron` simulation in the same geometry, with an **independent namespace RNG seed**, 0.1eV, isotropic start. Record `N_ion_birth`, `N_collected_anode`, N endpoint outcomes, `A_ion_hemi`, `A_ion_TB` and ion birth spatial reservoir, as well as full per-trial runtime and cap/failure status. New result is **paired by initial position and geometry only**; it is an independent microscopic realization and need not reach the same wire or reproduce the original drift endpoint. The new avalanche's terminal gain/attachment must be taken from this new run, not forced to equal the original old `collected` flag. Preserve old terminal fields as comparison/provenance columns only.

**B: EFFICIENT FUTURE PRODUCTION:** split at a physically justified **pre-avalanche** cylinder `r_switch` where local `alpha(E)` integrated along the remaining path becomes appreciable. Export crossings of this surface (position, energy/direction/time, wire ID, survival status) in a NEW drift kernel and draw a conditional near-wire microscopic avalanche response. The original raw LUT endpoints do **not** contain the state on that cylinder, so they cannot directly provide an exact split. An approximate `P(G,A_ion | phi_hit,incidence)` conditioned on old endpoint status needs explicit validation against FULL new `AvalancheElectron` starts and should never be called an exact paired physical trajectory.

**DO NOT double-condition electron survival:** gating on the old `collected==1` and then sampling an unconditional full-start avalanche response would count attachment/loss processes twice. Also never multiply all surviving electrons by an independent average gain and claim to preserve spatial gain–angle correlations. Keep zero avalanches and FW losses.

The original `LUT20/500` already supplies **electron-drift endpoint and arrival-time distributions**. The first avalanche extension will be an *independent correlated trial table* with original start identity plus full new avalanche response; the two can be statistically compared but are not the same particle history. A later validated split can truly compose drift and near-wire multiplication.

## Two immediate stages only

**Stage 1 measured cost:** run the already committed `run_pi_native_muon_b1.sh build`, then `pi_benchmark_native_muon_b1.py --seeds 1` (8 microscopic avalanches from the two real Heed muons), and optionally `--seeds 3` (24 avalanches). These are engineering timings only, neither full-muon A nor a physics sample. This code has not yet been observed compiling/running on the user's WSL host. Confirm timestamps, caps, seed wall times, failure rate, gas initialization overhead; do not extrapolate aggressively from 8/24 samples. Only advance after benchmark and local pilot pass.

**Stage 2 read-only source audit:** `python3 test/mwpc-stage-c/pi_audit_avalanche_lut_reuse.py --trials all`, checks all four geometries in LUT20 and LUT500, raw full geometry configs, complete row/cell/trial coverage, unique sample IDs, correct original start bounds, anode vs FW collection and SHA256. With `--export`, write a NEW immutable `lut_exact_start_index.csv` in `test_runs/mwpc_stage_c/pi_avalanche_lut_reuse_v1/frozen_exact_starts_all`, selecting 3 u positions × 6 v positions × 3 frozen trials per geometry per LUT (=432 total index rows); original eight LUT directories are never modified.

Next, after these tests: create a new **avalanche sidecar runner** which consumes these selected original exact-start rows, uses the same B1 full microscopic avalanche engine and outputs another SHA-sealed directory; not implemented or claimed run yet. Compare sidecar endpoints (distribution, not individual-history exact match) with old drift, timing costs, gain distribution with zeros and angular ion-birth metrics. Consider an adaptive cell allocation weighted by real Heed birth occupancy and high-variance gain/FW regions before scaling. Stratified rows oversample extrema, so they require weights for physical muon averages. Existing older large-scale crossing summaries generated from a marked-Poisson Heed surrogate, not a direct `TrackHeed` run per muon; keep those source models distinct.

**NOT to do now:** run all 84k or 2.1 million microscopic avalanches, overwrite raw LUT20/500 files, infer physical angular trends from a few artificial single-electron seeds, run B3 ion drift or U/V strip electronics until B1 pilot fidelity is established.
