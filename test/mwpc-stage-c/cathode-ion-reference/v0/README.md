# V0: archive-only statistics for the 2100 V K256 campaign

Reads `muon_NNNN/ion_births.root` and the exact per-(primary, wire) census `microscopic_wire_seed.csv` written by the
frozen Garfield++ source (`phaseB_native_muon_p3_noFW_penning057_cathode_ion_reference_v1.cxx`). Read-only. No
avalanche, no ion transport, no weighting field.

## Status
| Item | State |
|---|---|
| Source compiles (g++ 13, C++17, no ROOT) | yes |
| ROOT reader | syntax-checked against stub headers only; not yet compiled against real ROOT |
| Synthetic estimator checks | pass (see `run_v0_archive_stats.sh selftest`) |
| Real 2100 V archives | not yet run |

## Conventions verified in the repository
* `CathodeIonArchive.hh` is byte-identical on the overnight, independent-muon and v0-pilot branches.
  Branch types: `/I` Int_t, `/L` Long64_t, `/D` Double_t, `crossing_id`/`source_channel`/`storage_mode` std::string.
  The reader checks every required branch type and fails on any mismatch.
* `ions.local_wire_index = llround(u_patch_cm / pitch_cm)`: the anode nearest to the ION BIRTH. The collecting wire is
  known only from `collected_e` in the census CSV. The reader re-derives the index from `u_patch_cm` and fails if it
  disagrees.
* Each avalanche is started at the Heed primary time (`AvalancheElectron(u, v, w, e.t, ...)`), so
  `t_birth_ns - t_primary_ns` is the delay after the primary was created, including drift.
* `primaries.u_patch_cm` already includes the patch phase, i.e. it is in the same frame as the ions.
* Census CSV `crossing_id` is the same string as ROOT `muons.crossing_id`; its hemisphere convention
  (up: 0 < phi < 180 deg, exact 0/180 excluded) and nearest-wire attribution match the ROOT ion rows.
* Weights: the reader requires `inverse_inclusion_weight = ion_births_exact / ion_births_saved` per row and
  `ion rows per primary = ion_births_saved`, as `Plane2100Archive.hh` does.

## Run (WSL)
```bash
R=test/mwpc-stage-c/cathode-ion-reference/v0/run_v0_archive_stats.sh
bash $R prepare      # cmake target mwpc_v0_archive_stats
bash $R selftest     # synthetic checks
bash $R schema       # branch types of the first development archive
bash $R pilot        # 20 J/psi development archives (development20.csv)
bash $R population   # all SHA-verified archives except the 20 sealed ones (+ plots/ if matplotlib and pandas exist)
bash $R plots DIR    # re-plot an existing V0 output directory (v0_plots.py)
# optional 1800 V columns: export MWPC_LUT500_DIR=.../pi_pitch_fieldwire_500_p3_v1/pilot/p3p0_anode_only
```
Sealed slots are dropped before any file of theirs is opened or checksummed; the binary also refuses to run without
the sealed list. Return `v0_summary.txt`, `v0_manifest.txt`, `provenance.txt`, `selection_note.txt`.

## Outputs
`v0_summary.txt` sections: [1] primary census (ion births vs collected electrons, zero-ion primaries kept),
[2] nearest-birth vs collecting wire and reservoir coverage of wire fragments, [3] nested birth-time decomposition,
[4] K256 error measured against the exact census, [5] conditional descriptor tables, [6] rank correlations,
[7] effective-number / independence tests, [8] 2100 V transport/timing map by primary start position
(optionally next to the 1800 V LUT500, never merged; see `VOLTAGE_CONSISTENCY_AUDIT.md`). CSVs: `v0_fragments.csv`, `v0_primaries.csv`, `v0_wires.csv`
(additive per-wire moments for E3 calibration), `v0_muons.csv`, `v0_transport_map.csv`.

## Plots (`v0_plots.py`, added 2026-11-23)
Eight figure files plus `v0_plots_summary.txt`: delay vs start (u0, v0); A vs signed v0 and |du| (maps and profiles);
P(|A| > 0.9); conditional correlation profiles; within-avalanche ion-birth times (low/high-gain examples from
`v0_ion_examples.csv`, written with `--ion-dump-muons N`); onset spread within a muon and within a collecting wire;
collecting anode vs anode nearest the PRIMARY START; 1800 V LUT500 vs 2100 V (labelled, comparison only).
Bands are 68% muon-bootstrap intervals; cells with fewer than 20 primaries are blank.

`v0_plots_summary.txt` also has the repaired correlation estimator `rho_w` (between-primary residual correlation
within muons, leave-muon-out conditioning on 12 x 12 signed start cells x yield tercile, flagged when outside the
admissible range). Synthetic checks: independent law, rho_w(A) = +0.0015 [-0.011, +0.015]; correlated law
(`--rho-side 0.25`, expected about 0.2 to 0.24 from the generating law), rho_w(A) = +0.212 [+0.178, +0.247].
Residual synthetic bias on the delay estimate: -0.022 [-0.040, -0.003], so differences below about 0.03 are not resolved.

New CSV columns: `v0_primaries.csv` start_nearest_wire, exact_delay_ns, exact_sigt_prim_ns, exact_onset_ns,
collected_frac_max_wire; `v0_fragments.csv` exact_onset_ns, exact_sigt_collected_ns.

Fix (2026-11-23): `argmax_collected_wire` / `argmax_birth_wire` were set to global wire 0 for a primary whose census
rows all had zero count. The 52 / 42,879 count in section [2] should be re-checked with this build.
