# Independent validation: A20@D500 with 120 new real Geant4 muons

## Naming — IMPORTANT

- **DriftLUT500 / D500**: 500 ORIGINAL ELECTRON DRIFT TRIALS per (u,v) spatial cell. These do **not** include 500 new microscopic full Garfield avalanche realizations per cell.
- **AvalancheLUT20 / A20**: 20 independent NEW FULL Garfield microscopic avalanche calculations per cell, seeded at selected ORIGINAL D500 drift START coordinates. There are 4,200 cells and precisely 84,000 completed avalanche trials in total.
- **A20@D500**: the empirical full-avalanche model being validated. DO NOT refer to this as "LUT500 avalanche statistics". The old file/directory names containing LUT500 are archived physics provenance and must not be changed.

## Physical validation question

Does the original A20@D500 single-electron empirical-response model predict the distribution of **muon-level total ion births, zero-ion electron outcomes, birth hemisphere/top-bottom asymmetries, and Heed-cluster total ion birth fluctuations** for new genuinely simulated GEANT4 muons, relative to fully microscopic Garfield++?

*This does NOT increase the per-cell A20 training statistics; it tests out-of-sample generalization.*

## Independent holdout

- Select **120 truly different Geant4 J/psi parent events** from archived native crossing **chunk 2** (the previous 20-muon B2 pilot used chunk 1). Only one muon track may be chosen from any parent event; neither one track crossing multiple stations nor both muons from the same dimuon are counted as independent observations.
- Balanced 10 muons per one of 12 station/sign(theta_u) strata: MS0..MS5 × {theta_u < 0, theta_u >= 0}. Deterministic SHA-ranked candidate selection. NO invented muon energy, track or angle.
- The already used nine-anode analytic Garfield patch is valid over a restricted local transverse domain; require source physical chord |delta_u| <= 8 mm and |delta_w| <= 40 mm to leave margin for local coordinates. Thus selected tracks are a **domain-limited stratified sample, not the true incidence-angle PDF**. This cannot measure genuine detector-wide efficiency.
- Original Geant4 CSV and original ROOT hashes are read-only, verified against the original frozen chunk2 manifest. Newly selected 120 source identities and source SHA locked in their own input manifest.

## True microphysics and 12 resumable checkpoints

- The EXACT same original TrackHeed clusters and conduction-electron identities are generated ONCE for each native Geant4 muon. Every accepted original electron then drives one full Garfield++ AvalancheElectron in each of four fields p3_noFW, p3_FW75, p4_noFW, p4_FW75 (anodes +1800V, ±3mm cathodes, Ar/CO2 70:30, B=0). Full ionization-birth census preserved. Compact 5-ion spatial reservoir does NOT change full ion birth/sector counts.
- C++ B1 accepts 1..24 muon crossings, hence this new source is run in **12 separately SHA-protected physical chunks × 10 muons**. Each chunk must have C++ physics completion marker, the original cluster/primary correspondence and full per-muon ion census independently validated before it receives a `complete.json`. Any interrupted/incomplete chunk remains untouched for forensic diagnosis; completed chunks are SHA-checked and SKIPPED on resume.
- A completed run aggregates all original 12 physical source results into separate root-level CSVs, verifies all source identities, 120 true Geant4 source tracks/parents, all Heed clusters, every primary Garfield full avalanche and full four-geometry muon response, then creates final SHA manifest.

## Stochastic replay and uncertainty

- Uses precisely the **same stored physical Heed electron list**, never reruns Heed or Geant4 in the fast comparator.
- Reuses 84k original A20@D500 trial CSVs with SHA checks and **exactly 20 independent avalanche realizations per 0.1mm × 0.2mm bin**.
- For 120 original muons × 4 geometries × 200 conditional A20 replications = **96,000 synthetic muon-response draws** (these are computational reuses of a FIXED A20 table, not extra independent Garfield avalanches).
- Per muon: full Garfield total ion births, A20 conditional mean/sd/5th/50th/95th percentiles, fraction of full reference outcomes inside nominal 90% intervals, number of zero-birth Heed electrons, and A_hemi / A_TB measured from total ion birth sector counts (NOT the average of per-electron asymmetries).
- Per original Heed cluster: original position/time and electron count; full Garfield cluster total versus A20 conditional distribution. Full cluster totals are verified to sum exactly to the whole-muon full response.
- Per geometry: out-of-sample mean gain ratio and **1,500 paired bootstrap resamples across DISTINCT Geant4 J/psi parent events**, full-muon gain SD versus A20 total SD (between-fixed-muon-Heed variance + within-conditional-avalanche variance), 90% coverage, asymmetry bias and bootstrap intervals.
- Gaussian/statistical significance is not assumed from 200 A20 replications. All such replicates share only 20 original avalanche outcomes per cell and do NOT reduce the A20 finite-training uncertainty. Additional 120 independent full-Garfield muons check real validation variability, not A20 training adequacy at the few-percent level.

## One-time physical field-wire start recovery of the already completed original chunk_07 (29 Oct 2026)

This is ONLY for the actual stopped original run with original 7 SHA-complete
chunks and originally completed eighth **native C++ physics** execution,
whose Python reference validator stopped on eleven failed starts in a
p4_FW75 75um field-wire metal, real muon
\`jpsi40-bore|2|4138|1|2|1|38\`, primary indices 18..28.

Cause: \`ComponentAnalyticField::AddWire\` receives a **diameter** of 75 um,
i.e. a wire **radius** of 0.0375 mm, with p4 field-wire centres at
\`u = (...,-2,+2,...) mm\`, \`v=0\`. All observed 11 Heed electron
starting positions lie *inside the conductor*. Microscopic Garfield++
cannot initiate gas-electron transport at a point inside a metal wire;
it returns \`successful=0, garfield_ne=1, garfield_ni=0,
ion_births=0, capped=0\`. This is a physically **absorbed or embedded
initial electron with zero gas avalanche**, NOT a capped or unphysical
cascade. Original raw CSV statuses must remain unchanged; no claim is
made that Garfield's return value is a successful transport.

The Python physical validator now accepts ONLY this strict geometrically
confirmed case, including a C++ electron endpoint recorded as `other=1`:
C++ classifies any non-collected and non-attached endpoint as `other`,
so an electron born inside the metal can have one such endpoint even when
no avalanche is transported. Both `other=0` and `other=1` are permitted
only under the exact same strict inside-wire and zero-ion conditions;
`other>1` remains forbidden. This fixes the earlier false rejection by
the initial recovery validator that mistakenly demanded `other=0`.
It continues to reject:
- any nonzero avalanche ions, endpoints, collection/attachment or
  angular-sector counts in a supposed in-wire no-transport seed,
- any cap; every failed/no-transport seed outside field-wire metal;
- start points at anodes or other media unless separately validated.

Run after previous process exited and BEFORE resuming, without
rebuilding C++ or editing original physics files:

    cd ~/na6p/src/NA6PRoot-mwpc
    git pull --ff-only
    python3 -m py_compile test/mwpc-stage-c/pi_run_muon_lut_validation.py test/mwpc-stage-c/pi_holdout120_A20_from_D500.py
    python3 test/mwpc-stage-c/test_pi_muon_lut_validation.py
    python3 test/mwpc-stage-c/test_pi_holdout120_A20_from_D500.py
    python3 test/mwpc-stage-c/pi_holdout120_A20_from_D500.py recover07
    python3 test/mwpc-stage-c/pi_holdout120_A20_from_D500.py status

\`recover07\` asserts ORIGINAL input/C++ executable SHA, ALL seven
existing sealed checkpoints, original unchanged chunk_07 Geant4 source,
physical C++ end-of-run marker, all original Heed/ion CSV schemas,
no caps and precisely **11 physically in-wire zero-ion failed seeds**
in \`p4_FW75\`. It NEVER calls Garfield. If and only if all conditions
pass it writes the single missing \`chunk_07/complete.json\` containing
SHA256 of every original C++ result and a conspicuous recovery note.

The original Python exited before recording original eighth-chunk wall
clock, so no elapsed-time value is inferred or fabricated. Its timing
will be marked missing, and whole-120 measured runtime ratios will
remain unavailable unless separately and legitimately benchmarked.
The remaining original source chunks 08..11 will have genuine
monotonic wall/CPU timing. Every original chunk_00..07 CSV stays untouched.

**Physical modeling caveat for follow-up:** A20@D500 is currently a
0.1mm×0.2mm spatial-bin empirical sampler. If a Heed primary starts
inside a 75um field-wire metal, exact geometry dictates zero gas avalanche;
the empirical spatial bin might include both metal and gas starts.
The original A20@D500 sampler is NOT automatically geometrically vetoed
by this change, so future holdout results should explicitly inspect
any possible LUT-vs-Garfield bias caused by embedded-metal electron
starts. A dedicated documented metal veto may be needed in an upgraded
surrogate; do not conceal this discrepancy or silently retrofit the
historical 84k physical LUT.

## Reproducible runtime benchmark (BOTH simulation methods)

The 120-muon holdout additionally benchmarks **time**, with physically
comparable primary-electron populations and explicit accounting for setup.

- **Native full microscopic Garfield++/TrackHeed**: every 10-muon C++
  subprocess is timed with `time.perf_counter()` (monotonic wall) and
  `resource.RUSAGE_CHILDREN` (user/system CPU). Each immutable
  `chunks/chunk_XX/complete.json` includes measured elapsed seconds,
  CPU seconds, original Heed primaries, and full microscopic
  primary-electron×geometry avalanche count. Includes shared Heed
  generation, gas initialization, 4 physical Garfield fields and raw
  disk output; excludes Python provenance validation, selection and Drive.
  The reference `summary.json` and
  `full_garfield_chunk_timing.csv` combine all genuinely timed
  checkpoints.
- **Fast A20@D500 LUT**: separately measure cold loading + SHA checking
  the archived 84k original avalanche rows, and the warm conditional
  lookup + RNG + ion-count/angular-sector/Heed-cluster aggregation over
  **96,000 fixed-Heed muon×geometry replicas** (120×4×200).
  `benchmark_A20D500_replay_by_geometry.csv` reports total warm wall
  time, number of stored-electron draws, time per electron and time per
  muon×geometry replica. The detailed benchmark is also written into
  the comparison `summary.json` and hashed `manifest.json`.
- Compare ONE complete physical muon/4-geometry stochastic response:
  `speedup = sum_native_full_wall_120muons_4fields /
  (sum_warm_LUT_wall_120muons_4fields_200rep / 200)`.
  This speedup deliberately describes **native full Garfield+shared
  Heed+gas init+CSV I/O versus LUT response with Heed already cached**.
  Therefore it is NOT an end-to-end speedup for identical Geant4/Heed
  pipelines. To measure such a speedup later, Heed costs must be included
  in BOTH paths; for now the benchmark reports the scopes separately.
  Do not confuse wall time, user CPU, and total processing including
  preflight/validation/Drive.
- If original physical checkpoint(s) were generated BEFORE timing
  instrumentation, their `physics_timing` is missing. They remain
  immutable; the final benchmark explicitly reports
  `n_chunks_without_timing` and sets any complete-run speedup
  to unavailable rather than fabricating wall times.
- No new Garfield C++ rebuild is required by this Python benchmark;
  the original physics model and all old measured physics files stay
  SHA protected.

These metrics allow later profiling before trying safe optimizations
(preload LUT once, vectorize lookup and aggregation, compile sampling C++,
reduce CSV overhead, batch common Heed source and gas initialization).
Optimizations must preserve common Heed primary identities, zero-birth
distribution and angular counts, and must be revalidated against the
unaltered full-Garfield reference.

## Exact WSL execution

First update and run only safe/offline checks:

```bash
cd ~/na6p/src/NA6PRoot-mwpc
git pull --ff-only
python3 -m py_compile test/mwpc-stage-c/pi_holdout120_A20_from_D500.py test/mwpc-stage-c/test_pi_holdout120_A20_from_D500.py
python3 test/mwpc-stage-c/test_pi_holdout120_A20_from_D500.py
python3 test/mwpc-stage-c/pi_holdout120_A20_from_D500.py status
```

Next select **original** independent J/psi chunk2 Geant4 muons and validate their actual compiled C++ Heed source. This uses the EXISTING SHA-pinned successful B1 binary, and no rebuild is needed unless its source/binary SHA fails (the new holdout changes only Python scripts):

```bash
python3 test/mwpc-stage-c/pi_holdout120_A20_from_D500.py prepare
python3 -u test/mwpc-stage-c/pi_holdout120_A20_from_D500.py preflight
```

**Do not run the physical overnight if offline or preflight checks fail.**

After both pass, the whole long sequence can be run unattended. It will stop safely BETWEEN complete 10-muon chunks at a soft time limit of 7h45m; no partial reference or comparison is published. Repeating the EXACT command later checks existing chunk source hashes, skips completed ones, resumes physical Garfield on untouched new chunks and publishes only when all 120 source muons have finished. When all 12 pass, automatically compare A20@D500 with the full Garfield response and publish both original full-reference and comparison bundles to the existing Drive archive using `rclone --immutable --checksum` and independent readback.

```bash
nohup nice -n 15 python3 -u \
  test/mwpc-stage-c/pi_holdout120_A20_from_D500.py pipeline \
  --hours 7.75 --replicates 200 --upload \
  > "$HOME/na6p/muon_holdout120_A20D500_20261029.log" 2>&1 < /dev/null &
echo "Holdout PID: $!"
```

To check live state without altering physical data:

```bash
tail -n 20 ~/na6p/muon_holdout120_A20D500_20261029.log
python3 test/mwpc-stage-c/pi_holdout120_A20_from_D500.py status
```

If no full reference and new comparison output exists after safe partial time cutoff, run the SAME pipeline command again: it will reuse completed SHA-verifiable physical checkpoints. Do not manually delete existing incomplete chunks; inspect their physics.log and preserve original evidence.

If automatic rclone publication fails after physical calculations are complete, physical results and the A20 comparison are still retained on WSL. After fixing the rclone issue, simply re-run the SAME pipeline command; it will verify existing SHA manifests and only retry immutable publication, not rerun Garfield.

### Publication destinations

Parent: [NA60+ MWPC Google Drive avalanche archive](https://drive.google.com/drive/folders/1Xbh45V86BhVu-WLYSM5qtosqPEG363M7)

- `real_muon120_holdout_chunk2_full_Garfield_v1`: source identity, 120-muon full Garfield input, original TrackHeed clusters/electrons, all independent full microscopic ion birth and angular sector responses, source SHA manifests, 12 per-chunk SHA provenance records (large raw chunks remain on WSL).
- `real_muon120_holdout_A20_from_D500_compare_v1`: `muon_comparison.csv`, `heed_cluster_comparison.csv`, `muon_A20D500_replicas.csv`, statistical summary and original SHA manifest.

**CURRENT STATUS: CODE STAGED IN GITHUB; NO new physical 120-muon Garfield execution or holdout result has been observed.** Never claim predictive accuracy until real source validation, compiled Heed preflight, 12 Garfield SHA checkpoints, comparison, and rclone readback complete.
