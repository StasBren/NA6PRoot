# K2A fast cathode response: validation plan at the frozen 2100 V geometry

## Objective (clarified 2026-11-25)

K2A is meant to become the **fast cathode-response component inside the complete NA60+ simulation chain**:
avalanche generator → compact cloud moments → compiled K2A engine → electronics → digitization.
- The 997-muon archive is a **validation and benchmarking dataset**, not the scientific deliverable.
- A full offline 997-event response export is **not** a goal.
- 997 detailed S0 references will **not** be built.
- During K2A validation the 2100 V field, cathode geometry and provisional electronics stay unchanged. The geometry
  study stays paused.

**Current step:** the existing read-only preflight only (section 8). It establishes kernel coverage and characterizes
the available events: domain causes, wire/D/time coverage, dev-vs-population features, and the K256 reservoir
fraction (stored rows / exact births). The next implementation step is decided from its results.

**Priorities after the preflight**, in order, each to be planned in detail once the preflight has been reviewed:

1. **Representative development validation.**
   - Extend K2A numerical (vs G) and physical (vs S0) validation on a development subset chosen from the preflight
     coverage.
   - Include non-sealed φ/ω events where the J/ψ-only development set does not span the population.
   - The sealed 10 φ + 10 ω stay untouched until the criteria are frozen.
2. **In-chain interface.**
   - Design the interface between the existing avalanche generator and the compact-moment accumulator (the
     `BirthMoments.hh` per-(primary, wire) sums).
   - Moments are accumulated in memory as births occur, so the simulation writes no `ion_births.root`, no trajectory
     export and no CSV intermediates.
3. **Full-birth vs K256 moments.**
   - Accumulate the moments from the complete ion-birth population inside the generator, before reservoir sampling.
     The per-cloud sums are additive, so this needs no stored ions.
   - Quantify full-birth vs K256-weighted moments on a few muons, on moments first and then on K2A waveforms and
     observables.
   - The preflight's reservoir-fraction column shows where sampling acts.
4. **Electronics and digitization.** Connect the compiled engine to the electronics and digitization stages behind
   the existing interchangeable `Electronics` interface. The provisional shaper stays until the VMM3a configuration is
   fixed.
5. **Full-chain validation.**

Sections 2–7 below were written for the earlier "997-event processing campaign" framing.
- They remain valid as a toolbox: reuse map, coverage limits, cost scaling, observables, draft criteria, sealed-set
  protocol.
- **Section 7 (offline 997 production) is superseded**: the 997 are processed only as a benchmark and characterization
  pass, not as a deliverable.
- **Section 5's development set is to be revised** under priority 1.

**Status.**
- Only the preflight is implemented and enabled; it has been tested on SYNTHETIC inputs only.
- The runner refuses `phase2/3/4`.
- Unchanged: geometry, electronics stand-in, avalanche archives, sealed set, Geant4/Heed, the frozen references, the
  frozen baseline (f117e2a).

## 0. Inputs (existing, read-only)

| input | path (WSL defaults used by the runners) | role |
|---|---|---|
| 1000-slot K256 archive | `$HOME/na6p/cathode-overnight-2100/test_runs/mwpc_stage_c/overnight_2100_weighted_k256_v1` | 997 sealed `muon_NNNN/ion_births.root`; 3 slots are excluded (recorded as 1 failed, 2 unsealed) |
| audited inventory and frozen selection | `$HOME/NA6PRoot-signal40-2100/test_runs/mwpc_stage_c/phase_space_2100_signal40_v1` | `archive_inventory.csv` (997 rows: channel, station, det_id, pdg, Heed primaries, exact ion births, …), `excluded_slots.txt`, `development20.csv` (20 J/ψ), `sealed_validation20.csv` (10 φ + 10 ω), `previously_inspected_not_sealed.txt`, `FROZEN_SOURCE_SHA256.txt` |
| validated reference run | `test_runs/mwpc_stage_c/s_ladder_2100/run_20261124_124740` | the five S0/G/K2A references, `green_paths_r14t72edge5`, `green_r14t72edge5`, `kernel_r14t72edge5`, `engine_r14t72edge5_*` (kernel id), `timing.txt` |

## 1. Reuse

| need | existing tool (unchanged unless marked) |
|---|---|
| archive reading, census, weights, primary association, frame check | frozen `Plane2100Archive.hh` through `p3_birth_moments` (*changed*: `--keep-going` batch, `.event.json` metadata; moments CSV unchanged) |
| cloud moments directly from ion births | `p3_birth_moments --batch` + `BirthMoments.hh` (validated 3.5e-15 vs export path) |
| K2A evaluation | `engine/cathode_engine.py` + `kernel_eval.cc` (compiled; validated == loop/vector; 0.70–10.5 ms for 6–115 clouds), `engine_cli.py` |
| S0 reference for new development muons | `run_s_ladder.sh prepare` (frozen export → T1 → frozen replay + electron step), one muon at a time |
| common Green grid / table / kernels | `s_grid_spec.py`, `p3_green_grid_paths`, `s_green_table.py`, `s_kernel.py` via `run_s_ladder.sh green / kernel` |
| G closure and K-variants on new muons | `run_s_ladder.sh closure / synth / kernel` |
| shaped and unshaped observables vs S0 and vs G | `t_shaped_compare.py` (rel L2 shaped and unshaped, window charge, U/V ratio, lead strip, peak/time, threshold times, multiplicity, flips, second/lead, centroid u/w), `s_ladder_report.py` (per event vs S0 and vs G) |
| exact birth-vs-export equivalence | `s_compare_moments.py`, `s_compare_waveforms.py` |
| freeze manifests | `s_freeze_validation.py` |

## 2. New tools

| tool | phase | status |
|---|---|---|
| `k2a_inventory.py` | 1a | **in this commit**. All 1000 slots: seal files, full `sha256sum -c`, inventory/excluded partition, frozen-selection checksums, roles, batch lists. The sealed set is identified by its index column only. |
| `p3_birth_moments --keep-going` + `<id>.moments.csv.event.json` | 1b | **in this commit**. One failing archive no longer stops the batch; its verbatim error is kept. Muon metadata (crossing id, channel, station, det, pdg, seed, angles, p, path, K, primaries, zero-ion primaries, exact/saved rows) is written for the cohort summary. |
| `k2a_domain_check.py` | 1c | **in this commit**. Per cloud: the exact engine support tests (wire, r_eff, D_U/D_V for all ±5 strips, t vs window, finiteness). Per event: status and cause. Coverage by event, ion weight and cloud (overall, per channel, station, role). Cross-check vs `archive_inventory.csv`. Patch-edge-wire warning. Representativeness of dev20 vs the population. Extension requirement and cost from **non-sealed events only**. Sealed results blinded. |
| `run_k2a_campaign.sh` | 1 | **in this commit**: `selftest`, `inventory`, `moments` (chunked, 100 events per process; a crashed chunk shows up as NOT_CONVERTED), `domain`, `preflight`; `phase2/3/4` refused |
| `k2a_grid_extension.py` | 2a | to write after preflight review, and only if the preflight shows an extension is needed. It derives a **strict superset** of the r14t72edge5 grid: same radial nodes extended at the same log step, same D spacing, added wires. Old-node kernel values are therefore reproduced, so old kernel vs extended kernel on the five references can be checked to ~1e-6. It uses the union of the preflight requirement and the exact `s_grid_spec.py` output of the 20 development exports. |
| `k2a_phase2_compare.py` | 2b | per event × plane × Tp, from the `t_shaped_compare` jsonl (vs S0 and vs G); max errors; dependence on cloud count, collected charge (avalanche size), wires ≥ 5 %, station, angle, \|D\|; criteria evaluation |
| `k2a_freeze.py` | 2c | freeze record: criteria json, kernel_id, transport/weighting ids, engine config (model, half-strips, window, dt, shaper), git sha, development results sha. Phases 3/4 refuse without it. |
| `k2a_phase3_sealed.sh` | 3 | opens the sealed set only through the freeze record: S0 + K2A for the 20, reported per channel, no exclusion |
| `k2a_production.py` | 4 | engine (compiled) over the 997 moments of a fresh, hash-checked conversion. Per-event NPZ (U/V ion + electron-step currents) and per-event observables with the provisional shaper; `cohort_summary.csv/.json` |

## 3. Expected coverage limitations of the current kernel (`kernel_r14t72edge5`)

The kernel was built from the **union of five development muons**. Its support is therefore what those five
populated (`s_grid_spec.py`):

1. **Wires.** Only the patch wires k = rint(u/3) on which the five muons had ion births. Events with clouds on other
   wires are rejected by the engine. Likely causes are large |θ_u|, multi-wire sharing, and impact phase near a cell
   boundary. Expected to be the dominant out-of-domain cause.
2. **D range.** The supported range is \|D\|max(5 muons) + 5·1.7 + 0.5 mm, reduced by the ±5-strip window. Clouds far
   along w from the entry strip (large |θ_w|, long gas path) fall outside.
3. **r_eff range.** The radial grid spans the five muons' ion radii, and the kernel uses the cloud r_eff (log mean).
   Clouds are rarely outside. Separately, an innermost ion inside rn[0] matters only for G/S0 synthesis (ion level) and
   is reported as `G_ION_R_INNER`, not as a K2A cause.
4. **Time window (7000 ns).** A cloud born after the window is **dropped silently by the engine**, so the preflight
   lists it as out of domain. Expected to be rare (late δ-electron clouds).
5. **Field-model edge (not a kernel limit).** Clouds on |k| = 4 sit next to the edge of the frozen nine-wire analytic
   patch. They are reported as a warning.
6. **Channel coverage of the validation.** The development set is 20 J/ψ only, so Phase 2 cannot measure channel
   dependence. φ/ω enter only in the sealed set (Phase 3). The preflight prints input-feature distributions (clouds,
   charge, wire multiplicity, \|D\|, t_max) of dev20 against the non-selected J/ψ, φ and ω populations. That shows
   whether φ/ω occupy regions the development set never tested.
   - Option (needs your approval, not in the plan by default): add non-selected φ/ω events from the pool as
     supplementary development muons. They would be chosen by the same deterministic coverage rule, never from the
     sealed set.

The numbers are produced by the preflight. Nothing is clamped or extrapolated, and every out-of-domain event appears
with its cause and charge fraction.

## 4. Cost estimates

All Phase 2–4 entries are **expected** costs, scaled from measured runs. The preflight measures Phase 1 itself.

| step | compute | disk | basis |
|---|---|---|---|
| 1a inventory + full SHA of 997 archives | I/O-bound: total archive bytes / disk read rate (reported as `total_root_bytes`, `wall_s`); 2 parallel `sha256sum` | < 1 MiB | read-only |
| 1b batch moments, 997 events | read+validate is the five-muon converter cost × 997; moments µs–ms. Expected minutes to tens of minutes; recorded per chunk (`chunks.txt`, `/usr/bin/time -v`, `batch_timing.json`) | ≈ 10–50 kB/event → < 50 MiB | five-muon `batch_timing.json` |
| 1c domain check | seconds (CSV only; kernel hashed once) | < 10 MiB | — |
| Green extension (if needed) | paths ∝ new (wire, r) start points; table and kernel ∝ wires × r × D nodes. `extension_estimate.md` gives MiB and seconds scaled from the measured r14t72edge5 build. Engine resident RAM = K.npy size. | same scaling | measured build × factor |
| Phase 2: S0 for 15 new development muons | per muon: frozen export (RKF, not timed in the five-muon prepare; measure on the first) + replay 16–275 s; sequential, one muon at a time on the laptop | ≈ size of the five-muon `T1_export` + `S0` per muon | five-muon S0 `summary.json` |
| Phase 2: G, kernels, K-variants, comparisons | ≈ five-muon closure/kernel times × 4 | as five-muon run | — |
| Phase 3: S0 for 20 sealed | as Phase 2, per muon | as Phase 2 | — |
| Phase 4: 997 fast responses | conversion as 1b + compiled engine 0.7–10.5 ms/event + shaping/observables. Expected minutes; single process | 2 planes × 11 strips × 700 bins × 2 terms × float32 ≈ 0.12 MB/event → ≈ 120 MB + summaries | measured engine latency |
| 997 full S0 exports | **not planned**. The per-muon S0 cost (export + 16–275 s replay, plus T1_export/S0 disk per muon) × 997 buys nothing the 40 S0 references and the domain check do not already give. It would be justified only if Phase 3 shows a failure mode correlated with an input feature that the 40 references do not span. | — | — |

## 5. Phase 2: development validation (20 J/ψ, 5 existing references + 15 new)

1. **Kernel decision.** If the preflight shows non-sealed events out of domain, build the strict-superset common grid
   once (2a), covering all non-sealed requirements and the 20 development exports. Then rebuild the Green table and
   kernels with a new tag. The kernel that passes Phase 2 is the one frozen for Phases 3–4.
   - Regression check: extended kernel vs `kernel_r14t72edge5` on the five references, K2A waveforms ≤ 1e-6 relative.
     If the nodes are a strict superset, the only differences should be floating-point.
2. **S0 expansion.** For the 15 new development muons, run `run_s_ladder.sh prepare` one muon at a time (frozen
   export, T1, frozen replay, electron step), then `closure` (G vs S0), `kernel`, `births` (birth moments == export
   moments) and `engine` (compiled == export path).
3. **Comparisons.** Both K2A vs G (kernel compression error) and K2A vs S0 (total error), ion-only and with the
   electron step, at Tp 100/200 ns (provisional CR-RC³, absolute thresholds from the reference), with
   `t_shaped_compare` + `s_ladder_report` + `k2a_phase2_compare.py`.

**Observables** (per event, per plane, per Tp):
- unshaped rel L2 (current and cumulative charge) and shaped rel L2;
- window charge per plane and U/V charge ratio;
- lead-strip identity, lead peak and peak time, second/lead ratio and per-strip peak changes (≥ 5 % strips);
- threshold crossing times at f = 2, 5, 10 %, multiplicity and above/below flips;
- centroid w and stereo u (U = w − Au, V = w + Au), default and fixed-membership centroids.

Per-event **maximum** errors are reported, with their dependence on cloud count, collected charge, wires ≥ 5 %, station,
θ_u/θ_w and \|D\| (rank correlations plus a table sorted by worst error). Channel dependence is not measurable here
(J/ψ only).

**Proposed acceptance criteria.** These are drafts, not the geometry-study 2× rule. The values are set from detector
needs and finalized and frozen after Phase 2, before any sealed file is opened. Each holds per event and plane for
**K2A vs S0**; K2A vs G must also be reported.

| observable | proposed limit | rationale |
|---|---|---|
| lead-strip identity | identical in every plane, unless the reference top-two peaks are within 1 % (counted as ambiguous, listed) | a lead change moves the cluster by a strip |
| multiplicity / hit flips at f = 5 %, 10 % | 0 flips, except strips whose reference peak is within ±2 % of threshold (marginal, listed) | discrete readout |
| threshold crossing time (f = 10 %) | \|Δt\| ≤ 1 ns | ≪ 10 ns bins and shaper rise |
| window charge per plane | \|ΔQ\|/Q ≤ 0.5 % | gain-calibration scale |
| U/V charge ratio | \|Δ\| ≤ 1 % | five-muon worst ≈ 0.7 % (informative only, not 2×) |
| second/lead (sharing) | \|Δ\| ≤ 0.005 absolute | centroid lever arm |
| w (along wires) | \|Δw\| ≤ 10 µm | ≤ 10 % of an expected ~100 µm strip-centroid resolution |
| u (stereo) | \|Δu\| ≤ 50 µm | u lever arm 1/(2A) = 5; ≤ 10 % of an expected ~0.5 mm stereo resolution |
| shaped rel L2 per plane | ≤ 1 % (diagnostic ceiling) | waveform-level sanity; the observables above decide |

**Stopping conditions in Phase 2:**
- Any read, census, metadata or equivalence failure.
- An S0 replay failure.
- A kernel-id, transport_id or weighting_id mismatch.
- Any development event out of domain of the candidate kernel.
- Any criterion failure.

The response to a stop is diagnosis on development data only. If the model changes, Phase 2 is redone. The sealed set
is never used.

## 6. Phase 3: sealed validation (10 φ + 10 ω)

The runner refuses to start without the freeze record (criteria, kernel_id, ids, engine/shaper config, git sha).
Steps:
- S0 for the 20 sealed (as in Phase 2) and K2A with the frozen kernel.
- Results reported **separately per channel**, plus all events.
- Every failure stays in the table, including any event out of domain and any criterion miss.
- No retuning. If the model is revised afterwards, the sealed set is spent: it is reported as such and no longer
  called independent.

## 7. Phase 4: 997 fast-response processing (SUPERSEDED: benchmark/characterization only, see Objective)

- Fresh conversion of all 997 with the frozen binary; moments hashed and compared with the preflight moments, which
  must be identical.
- Engine with the frozen kernel, compiled.
- Per event:
  - an NPZ with U/V ion and electron-step currents;
  - observables with the provisional shaper at 100 and 200 ns (charges, U/V ratio, lead strip, multiplicity,
    threshold times, centroid w/u).

`cohort_summary.csv` + `.json` hold, per event:
- status (OK / OUT_OF_DOMAIN with cause / READ_FAILED / NO_IONS), channel, station, det, role;
- clouds, wire multiplicity (any and ≥ 5 %), collected charge;
- the observables above;
- runtime (read, moments, engine, shaping) and peak RSS;
- config hashes: kernel_id, transport_id, weighting_id, engine config, git sha, archive sha.

Outliers are flagged by robust z (median/MAD) per channel × station on charge, U/V ratio, multiplicity and runtime.
They are listed, not removed.

**Stop if:**
- any configuration hash mismatches;
- more than 1 % of events fail technically (investigate before continuing);
- any archive SHA changed.

## 8. Phase 1: exact WSL commands (read-only)

Review first. The commands only read the archive, the selection and the validated run. They write only below
`test_runs/mwpc_stage_c/k2a_campaign_2100/preflight_<timestamp>/`.

```bash
cd ~/NA6PRoot                      # your checkout of StasBren/NA6PRoot
git fetch origin
git checkout dev/na60-2100v-validation-20261123
git pull --ff-only origin dev/na60-2100v-validation-20261123
git log --oneline -1               # the commit named in the reply

# 0. environment as for run_s_ladder.sh; MWPC_BIRTHS_BIN = your p3_birth_moments build (rebuild: new options)
export MWPC_BIRTHS_BIN=<path to your spatial-ablation build>/p3_birth_moments
cmake --build "$(dirname "$MWPC_BIRTHS_BIN")" --target p3_birth_moments
C=test/mwpc-stage-c/cathode-ion-reference/k2a-campaign
RUN=test_runs/mwpc_stage_c/s_ladder_2100/run_20261124_124740
ls -d "$RUN"/kernel_r14t72edge5 "$RUN"/engine_r14t72edge5_*        # the validated kernel and engine runs

# 1. synthetic self-test (no archive)
bash $C/run_k2a_campaign.sh selftest /tmp/k2a_selftest_$(date +%Y%m%d_%H%M%S)

# 2. regression: the rebuilt converter reproduces a frozen reference moments file byte for byte (one muon, read-only)
REF=$(ls "$RUN"/muon_*/births_*/birth_moments.csv | sort | tail -1); M=$(basename "$(dirname "$(dirname "$REF")")")
NEW=/tmp/k2a_regress_${M}_$(date +%H%M%S).csv
"$MWPC_BIRTHS_BIN" "$HOME/na6p/cathode-overnight-2100/test_runs/mwpc_stage_c/overnight_2100_weighted_k256_v1/$M/ion_births.root" "$NEW"
cmp "$NEW" "$REF" && echo "PASS converter regression $M (moments byte-identical to $REF)"

# 3. the read-only preflight (inventory + full SHA, batch moments, kernel domain)
MWPC_K2A_REFERENCE_RUN="$PWD/$RUN" MWPC_K2A_TAG=r14t72edge5 bash $C/run_k2a_campaign.sh preflight 2>&1 | tee /tmp/k2a_preflight.log
```

Return for inspection (do **not** open `domain/SEALED_DO_NOT_OPEN/`):
- `inventory/inventory_check.json`
- `inventory/excluded_slots_detail.md`
- `domain/coverage.md`
- `domain/extension_estimate.json`
- `moments_nonsealed/chunks.txt`
- the tail of `/tmp/k2a_preflight.log`

The preflight can be split into its stages, `inventory`, then `moments PF`, then `domain PF`. Set `MWPC_K2A_SHA_JOBS=1`
if the laptop is busy.
