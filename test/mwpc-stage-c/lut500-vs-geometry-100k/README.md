# NA60+/DiCE — 100k muons per channel/station: LUT500 vs geometric midline

Prepared 18 November 2026. **Status: source package prepared and synthetic tests executed; the actual physics runs have NOT been executed in this environment.**

## Scientific question

For native frozen Geant4 muon gas crossings from the **2-million-parent/channel** campaign, compare assignment of the same source electrons and the same gain weights to physical MWPC anodes by (a) nearest-anode geometric midline split and (b) conditional drift LUT500. Study MS0, MS3, MS5 for each of J/ψ, φ, ω (nine 100,000-crossing samples), then determine per-muon Q₂/Q_total and the number of *significant collecting physical anodes* at 1%, 5%, 10% charge thresholds.

**LUT500 supports >2 anodes per muon.** It samples Δwire ∈ {-1,0,+1} for each primary electron *relative to that electron's own starting wire cell*. The global physical ID is `k_base + Δwire`. A muon crossing many wire cells may therefore involve many anodes. No arbitrary cap of two anodes occurs in this implementation. A *single* electron migrating past a second neighboring wire is not represented by this LUT.

## Inputs and provenance

- Existing WSL repository branch `feature/mwpc-stage-c` with original Geant4 source chunks and `test/mwpc-stage-c/pi_2m_native_six_geometry_fold.py`. We read the **SHA-verified** `geant_muon_crossings.csv` exports via `crossing_provenance(channel,chunk)` and never rerun Geant4.
- `NA60_LUT500.zip` (original frozen 900 cell × 500 drift trajectories). Each cell has 0.1 mm in u × 0.2 mm in v, for ±1.5 mm in u and ±3 mm in v. The LUT was generated **without Penning transfer**, so it is treated as a **drift-only conditional-on-collection** kernel.
- Existing Heed **calibration** `c_heed_library_2p4/heed_2p4_tracks.csv` and `heed_2p4_clusters.csv` via `read_heed`. The fast study samples a marked-Poisson cluster process *along the actual, native Geant4 chord*, exactly the type of approximation previously used in `pi_six_geometry_native_timing.py::birth_sample`. **This is NOT a true TrackHeed microscopic simulation for each of the 900k crossings.** Actual per-crossing TrackHeed replay requires a separate Heed-only C++ export, because the 2M Geant4 pool does not contain frozen Heed seeds.
- Optional but strongly recommended: J/ψ calibration `model_calibration_jpsi_v1.npz` (keys `cell`, `gain`), sampling **one shared empirical Penning-ON electron gain for both geometric and LUT paths**. Zero-gain seeds remain zero for both. Do not apply an extra Penning-OFF LUT survival veto after sampling a donor.

The native detector uses 3 mm pitch, ±3 mm cathode gaps, 30 μm anodes, 1800 V, Ar/CO₂ 70:30, B=0 T, no field wires as a common **P3 noFW model**. Geometry indexes each finite wire lattice with the station-correct half-pitch phase (`N_MS0=102`, `N_MS3=171`, `N_MS5=171`) in the local-chamber coordinate. It is **not a final physical claim about MS5's coarse MuID readout segmentation**.

## Standard GitHub → WSL Bash execution

The study is committed in the `test/mwpc-stage-c/lut500-vs-geometry-100k/` directory of the dedicated GitHub study branch. **No manual file copy or ZIP download is necessary.** The branch is based on `feature/mwpc-stage-c`, so its existing native Geant4 exporters and Heed calibration interfaces are available. Git switching preserves the local untracked `test_runs/` data but will stop if tracked changes would be overwritten; do not force-reset your current work.

Once the GitHub branch has been fetched and checked out in `~/na6p/src/NA6PRoot-mwpc`, run the following from the repository root:

```bash
cd ~/na6p/src/NA6PRoot-mwpc
STUDY=test/mwpc-stage-c/lut500-vs-geometry-100k
bash "$STUDY/run.sh" test
bash "$STUDY/run.sh" export
bash "$STUDY/run.sh" select
bash "$STUDY/run.sh" pilot
# After reviewing the 9 pilot summaries:
bash "$STUDY/run.sh" full
```

`export` SHA-verifies and resumes frozen Geant4-to-CSV exports for up to 40 x 50,000-parent chunks per channel. No Geant4 transport is repeated. `select` chooses 100,000 eligible crossings in each of the nine channel/station strata; it **fails rather than silently shortens a sample**. The `pilot` takes the first nested 10,000 of each saved 100k sample, so it is strictly paired to the full run.

The runner first reuses any SHA-verified files in `~/na6p/studies/lut500_vs_geometry_100k/frozen_inputs/`, otherwise automatically retrieves the exact historical `NA60_LUT500.zip` and `model_calibration_jpsi_v1.npz` from the validated Fast4/Fast5 Google Drive folder using the existing WSL `na60drive:` rclone remote. It SHA-256 verifies both inputs against the frozen archive manifest before the first physics event. **No re-generation of LUT500 and no new Geant4 transport are performed.** The cached originals are reused on later runs. If rclone is unavailable, point explicitly at the original files:

```bash
LUT_ZIP=/absolute/path/NA60_LUT500.zip \
CALIBRATION_NPZ=/absolute/path/model_calibration_jpsi_v1.npz \
bash "$STUDY/run.sh" pilot
```

Do not use a different LUT calibration silently. Requires the existing Heed-calibration CSVs under `test_runs/mwpc_stage_c/c_heed_library_2p4/`. The scripts use NumPy and optionally matplotlib for plots; install them in your existing WSL Python environment if needed.

## What is saved

`~/na6p/studies/lut500_vs_geometry_100k/` contains:

- `selection_manifest.json` with source exposure, selection policy, valid/rejected/overlap counts, and SHA-256 per source sample.
- 9 source-station `selected_*_100000.csv.gz` files, one **physical muon per station**, longest valid gas chord if several overlap-module hits exist.
- `comparison_10000_donor/` and `comparison_100000_donor/` (separate, non-overwriting output folders): 9 `events_*.csv.gz`, 9 `summary_*.json`, 9 `angle_*.png`, a consolidated `comparison_table.csv` and `comparison_table.md`, and `run_manifest.json`.
- `n_positive` and `n_gt{1,5,10}pct` counts on both geometric and LUT anodes, no multiplicity cap. Multiplicity histogram contains 0/1/2/3/4/5+; summary also reports max positive-wire multiplicity.
- Parent-clustered paired 95% uncertainty on the difference in >5%-second- and third-anode rates, using original `(chunk,event)` cluster IDs.

Important: one selected event is one **physical muon crossing at a station**, not one generated parent and not one detector strip/channel firing. Events with no amplified collected charge remain in the denominator; near-edge out-of-aperture charge is separately counted. No ion drift, weighting fields, cathode strip charge or ASIC thresholds are modeled.

## Guardrails and scientific limitations

1. The full source ROOT campaign and any unexported chunk data are on the user's WSL/Garfield machine, **not** mounted in this chat's container. The scripts have not been executed there, and the nine new physics percentages are not yet known.
2. For a fair comparison of the **wire-assignment mechanism alone**, geometric and LUT variants use the exact same electron starts and **the exact same gain draws**. Differences are caused by drift-wire assignment. Previous results comparing cluster-nearest geometry versus donor-weighted Fast4 also changed the gain model and must not be conflated with this controlled comparison.
3. The sample is native Geant4 and is representative of its incidence distribution, but the primary-ionization realizations are a **marked-Poisson Heed-calibrated approximation**, not per-track TrackHeed truth. If we want truly exact TrackHeed, instrument a Heed-only export for the selected gas chords and replay it for both methods.
4. For high-angle tracks across many physical wire cells, LUT can represent ≥3, ≥4, and more collecting anodes. For a **single electron** moving >1 cell away, LUT500's stored offsets cannot represent that motion; count and flag such configurations in an extended LUT study if needed.
5. Multiple crossings of one physical muon within the same station are reduced to the longest valid gas chord. That rule is explicit to avoid double-counting a physical muon as two independent observations; if station-level total signal from overlapping modules is the target, the definition must change.
6. A full-truth study should also validate station/chamber-specific geometry assumptions for MS5 (MuID), and check whether rates conditional on the 5% threshold are robust to the underlying gain fluctuations and zero-gain events.
## Postrun fine-angle diagnostics (18 November 2026)

After the 100k full run, use:

```bash
bash test/mwpc-stage-c/lut500-vs-geometry-100k/run.sh diagnose
```

This is a **read-only** scan of the nine existing `events_*_MS*.csv.gz` and `summary_*_MS*.json` files. It verifies exact per-event/summary agreement for >5% two-anode and three-anode rates and hashes the frozen inputs. It writes only to a NEW `comparison_100000_donor/angle_tail_diagnostics_v1/` directory and refuses to overwrite existing diagnostic output. No physics, sampling, LUT/Geant4/Heed, or gain draw is rerun.

Key products are `angle_bins_fine.csv` (5° bins from 20–50°, then wider to 180°), `high_angle_muons_ge30.csv` (individual original track IDs, inclinations, source cells, anode multiplicity, charge fractions and edge losses), `three_plus_anodes_5pct.csv` (all muons with >=3 significant anodes in either method), and `standardized_0_20.csv` (rates reweighted to the same pooled 0–20° angle-bin mix). The latter isolates differences in broad-angle composition **descriptively**; it cannot establish a causal detector effect, and uncertainty bars are not derived there. The rare >=30° bin should not be treated as uniformly distributed or as a validation of detector geometry without inspecting individual muon tracks.
