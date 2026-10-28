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


## 29 Oct 2026: local timing AND frozen-LUT audit passed

The user completed `pi_benchmark_native_muon_b1.py --seeds 3` on WSL:
- 24 direct Garfield `AvalancheElectron` calculations on 6 true TrackHeed electrons per geometry; wall 7.17s total including gas init.
- Mean seconds / avalanche: `p3_noFW=.097`, `p3_FW75=.447`, `p4_noFW=.109`, `p4_FW75=.510`.
- Max per avalanche: `.149,.824,.302,1.094` seconds respectively.
- Rough full 2-muon + 4-geometry 81-Heed-electron extrapolation: `94.47s` excluding init and Heed; **NOT validated as actual full-run timing**.
- The 24-row console output does NOT show whether any seeds hit the cap; inspect the `summary.json` / raw timing rows before interpreting gain or runtime tails.

The user also ran `pi_audit_avalanche_lut_reuse.py --trials all --export`. ALL eight complete original LUTs passed full file validation with **84,000** LUT20 and **2,100,000** LUT500 original independent Garfield drift rows. Stratified subset index of **432** starts created in a new immutable directory; 54 starts from each (LUT20/LUT500 × 4 geometry) stratum.

Old anode-collected counts (NOT ion avalanche): LUT20 `14875,14258,19383,19041` and LUT500 `368796,356644,484762,474787` respectively `p3_noFW,p3_FW75,p4_noFW,p4_FW75`. Old explicit `field_wire_collected` counts were zero in every geometry; this is an observational property of the old drift-classification fields, **not proof that field wires have no effect on electron transport**. Review original status/other-endpoint classifications separately as needed.

These 432 independent full-start avalanches might cost `2 × 54 × (.097 + .447 + .109 + .510) ≈ 126s` at the 24-sample means, **not a performance guarantee** because the input electron-start distribution is deliberately not muon-weighted, and tails/size caps may matter.

## New implementation: original exact-start microscopic ion-birth sidecar

New C++ mode `--lut-sidecar index.csv output_prefix` in the **same** `phaseB_native_muon_3p3_avalanche.cxx` and same Garfield `RunGeometry` microphysics as the verified native-muon B1 benchmark. Unlike native muons, exact-start coordinates are ALREADY in the wire-centered local unit cell, so the LUT sidecar explicitly sets entrance phase to **zero**, avoiding accidental chamber anode-parity offset. No `TrackHeed` is used in this sidecar because these are Garfield-only original grid start electrons, not physical muon clusters.

New Python runner `test/mwpc-stage-c/pi_run_avalanche_lut_sidecar.py` verifies every original LUT source SHA, selected exact starts, the current C++ executable SHA-stamped build, and a dedicated C++ input parse preflight (`--validate-lut-input`). It strictly validates: one row and microscopic timing per original start; exact (u,v,w) positions; independent new RNG; original old outcome labels retained for metadata **but NEVER used as a gate**; ion birth and sector counts; reservoir size min(Nions,100), reservoir weight and provenance. Results are stored **ONLY** under a new immutable directory. Original LUT20/500 files and all earlier physics results remain read-only.

New two-stage pilot:
```bash
cd ~/na6p/src/NA6PRoot-mwpc
git pull --ff-only
bash test/mwpc-stage-c/run_pi_native_muon_b1.sh build
python3 test/mwpc-stage-c/pi_run_avalanche_lut_sidecar.py check
python3 test/mwpc-stage-c/pi_run_avalanche_lut_sidecar.py smoke8
```
The first `smoke8` uses exactly ONE original saved start per `(LUT20/LUT500 × four geometries)`; output:
`test_runs/mwpc_stage_c/pi_avalanche_lut_reuse_v1/run_smoke8_independent_full_start_avalanche_v1/`.

Only if smoke8 C++ compile/run and output validation PASS, no size caps or failures, consider `python3 test/mwpc-stage-c/pi_run_avalanche_lut_sidecar.py full432`. This is a separate explicit action that NEVER starts automatically. The run has all 432 frozen exact-start identities; again NOT true per-muon asymmetry. This newest C++ sidecar mode has NOT been compiled or executed on WSL yet. User must run a new `build` after `git pull` because the C++ source changed.

Do not treat a new avalanche as a continuation of the original old drift. Do not gate full avalanche on old `collected` because that would count loss twice. The sidecar is a NEW independent microscopic response drawn from the same exact starting coordinate.
