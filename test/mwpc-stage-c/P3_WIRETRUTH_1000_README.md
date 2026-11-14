# Independent native Garfield++ collecting-wire truth, 1000 muons

**Purpose:** finally validate the Fast4/Fast5-T physical assignment of real collected electrons to one or several MWPC anodes, and save the per-anode ion-birth space/time moments needed for the next cathode-response step.

## Input and fixed physics

- 1000 genuine Geant4 muon gas crossings frozen from the prior 30k archive: 334 Jpsi, 333 Phi, 333 Omega, one crossing per previously selected parent, six stations. SHA256-ranked random selection with no angle enrichment; the exact 1000 crossings are bundled in `p3_wiretruth_1000_selection.csv` and its signed provenance manifest. The 30k TrackHeed electron *positions* were not retained in that CSV; this new reference uses deterministic **new TrackHeed streams**, and saves the actual electron coordinates for exact subsequent fast-model replay.
- A separate C++ executable cloned from `phaseB_native_muon_p3_noFW_jpsi_penning057_highheed_v1.cxx` (the Penning-ON full microscopic campaign source). The old source, validated library and research data are not modified.
- Geometry: p3_noFW, 3 mm pitch, 30 um anode diameter, cathodes v=+/-3 mm, 1800 V, Ar/CO2 70:30, Penning r=0.57, B=0. All genuine TrackHeed conduction electrons are propagated with full microscopic Garfield++ `AvalancheElectron`, with the same 200,000 cap and 0.1 eV thermal start convention as previous reference. Any capped event is fatal for reference validation.
- The reference is 1000 *muon crossings* selected from old Geant4 records; not a fresh Geant4 transport run and not a measurement of detector trigger/efficiency.

## New exact physical-anode observables

Original native files (`_crossing.csv`, `_seed.csv`, `_heed_primaries.csv`, `_heed_clusters.csv`, `_ion_birth_sample.csv`, `_timing.csv`, `_muon_timing.csv`) are retained. New:

- `microscopic_wire_seed.csv`: each (muon, Heed primary, physical anode) with **actual collected electron endpoints** on that wire, **exact** ion-birth census classified by their nearest physical anode, signed upper/lower and top/bottom birth asymmetry, exact ion-birth w/time means and std, and electron collection-time mean/std. Zero-gain seeds remain represented in the original `_seed.csv`.
- `microscopic_wire_muon.csv`: the same fields aggregated for each actual muon and physical anode, including charge across multiple wires.
- Ion birth classified by nearest wire is not the same thing as the originating seed's destination wire, and neither equals an induced cathode charge. We keep both separately.
- Global physical wire indices account for the **station-dependent 1.5 mm shift of station 0**; do not aggregate local `-1,0,+1` wire IDs across different muons.
- First/second wire charge fractions, angular dependence, correlation with Heed primary number and muon channel can now be evaluated directly against these full microscopic Garfield++ reference outputs.

## Execution (WSL / Ubuntu)

Work in an isolated Git worktree using branch `feature/p3-wiretruth-1000-20261114`, from `StasBren/NA6PRoot`. Do **not** switch or reset the existing `~/na6p/src/NA6PRoot-mwpc` working tree. The user has previously installed Garfield++ under `~/na6p/install/garfieldpp`.

First do a 20-muon native pilot if desired:

```bash
NA60_WIRETRUTH_PILOT_ONLY=1 bash test/mwpc-stage-c/run_p3_wiretruth_1000.sh
```

Then complete all 1000 (will skip any SHA-sealed completed chunks):

```bash
bash test/mwpc-stage-c/run_p3_wiretruth_1000.sh
```

An end-to-end run compiles only the new target, verifies C++ input/Heed, runs the first 20 real muons, verifies exact wire conservation, and then resumes the remaining batches of 20; all outputs, source inputs and binary hashes are recorded in `complete.json`. An unfinished chunk is never silently overwritten; inspect its log if interrupted.

Results:

`test_runs/mwpc_stage_c/p3_noFW_penning057_wiretruth_1000_v1/`

The runner also creates `wiretruth_compact_1000_muons.zip` in its results folder: it includes the per-anode truth and actual Heed primaries for model replay, omitting the large ion-birth-reservoir CSV (which is retained locally). Upload this compact ZIP to ChatGPT for analysis.

Reports after completion:

- `reference_wire_summary.md`
- `reference_wire_sharing_by_angle.csv`
- 50 chunk directories, each with per-anode truth and a SHA manifest.

## Important limitations and next comparison

A native Garfield/Geant4 run cannot be executed through the GitHub connector itself; this code must compile and run on the user's Garfield-equipped WSL machine. Neither an answer nor a previous cross-model agreement constitutes a Garfield collecting-anode validation before this new simulation completes.

The **same fresh Heed primaries** must be passed into Fast4 and Fast5-T for a fair matched-muon comparison. Their existing 30k aggregate-input adaptation is insufficient for exact 1000-event replay. Once the new outputs are back, construct a thin adapter for these new per-primary columns and test the per-wire charge distribution directly against the full microscopic truth, including Q2 fraction, n_anodes and angle bins. The proposed second 10k per channel run should be decided using statistical uncertainty and conditional multi-anode discrepancies in the first 1000.
