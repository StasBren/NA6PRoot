# K2A architecture review before any kernel extension: strip frame, recentering, near-wire grid, development set

**Status (2026-11-25).** This is a design review and tooling only.
- No Green table, kernel, frozen file, archived result or sealed file is touched.
- The tools ran on SYNTHETIC/toy data only. The real numbers come from running them on your preflight directory
  (section 8); I have only your summary of that preflight.

## 0. Objective (unchanged, explicit)

The production interface is **in memory, inside the simulation chain**:

```
Geant4 crossing -> Heed primaries -> Garfield avalanche (ion-birth callbacks)
   -> per-(primary, wire) moment accumulator, fed per birth, BEFORE K256 sampling
   -> cloud moments (N, t, r_eff, D_U, D_V, C1..S2, A) + muon strip frame
   -> compiled K2A engine (per-primary strip recentering, absolute U/V strip indices)
   -> electronics (shaper / VMM model) -> digitization (thresholds, ADC/TDC per strip) -> hits
```

- `ion_births.root`, trajectory exports and CSV moments exist **only** to validate this component.
- The 997 archives are a validation and benchmark set.
- The geometry study stays paused. The 2100 V field, cathode geometry and provisional electronics are unchanged.

## 1. What D_U, D_V and base_U, base_V mean, and why D reaches −14 … +21 mm

From the frozen exporter (`p3_export_full_spectral_paths.cc`), reproduced exactly by `BirthMoments.hh`. Lengths are
in cm inside, mm out; A = 0.1, P = 1.7 mm.

```
base_U = llround((w_in - A u_in) / P)          base_V = llround((w_in + A u_in) / P)        (muon ENTRY point)
D_U    = 10 (w + w_in - A patch_shift - P base_U)    D_V = 10 (w + w_in + A patch_shift - P base_V)
strip coordinate of an ion on patch position u:  eta_U = D_U - A u_patch,  eta_V = D_V + A u_patch
```

- **base_U, base_V** are the strips under the muon **entry point**. The 11-strip output (±5) is centred there once
  per muon.
- **D** is an ion's (or cloud's) longitudinal position relative to that entry strip.
  - Substituting base gives D_U ≈ 10 w_rel + A (u_in − patch_shift) ± P/2. The entry-offset term is at most a
    few tenths of a mm, because the patch is centred on the entry wire.
  - So |D| ≫ P means the cloud itself lies many strips along w from the entry point. The likely causes are an
    inclined track crossing the 6 mm gap (Δw ≈ 6 |tan θ_w|), long gas paths and δ-electrons.
  - `k2a_recenter_recovery.py` tests this on your events: the Spearman correlation of the D span with 6 |tan θ_w|,
    the top-12 table with θ_w, path and m, and the maximum entry-offset term.
- **Consequence 1: the output window is the wrong physical object for such events.** A cloud with |D| ≈ 15–21 mm sits
  9–12 strips from the entry strip. Most of its induced signal is on strips **outside** the base-centred 11-strip
  output, for S0 and K2A alike. The preflight's D failures are therefore not only a kernel-size problem.
- **Consequence 2: the frozen S0 cannot compute such events with a base-centred window either.** It refuses with
  "query outside declared strip phase range" (LUT phase range ±12 mm) or with "selected channels omit the landing
  strip; widen or recenter". This was demonstrated on a toy far event (`tests/test_recenter_synthetic.sh`). None of
  the five references had this situation, so they could not show it.
- A direct D extension of the Green table would compute the response of far-away strips to such clouds. That costs
  ~3.6 GiB, and it still misses the cloud's own strips.

## 2. Per-primary recentering with absolute-strip accumulation

**Rule.** For each Heed primary p:
1. m_U(p) = ⌊D̄_U(p)/P + ½⌋ and m_V(p) likewise, with D̄ the charge-weighted mean over the primary's clouds.
2. Evaluate its clouds with D′ = D − mP on strips base + m + o, for o = −H…H.
3. Add the results onto **absolute** strip labels U_j, V_j.

Whole primaries are used so that each group has an integer exact census, which the frozen replay needs for its
reference.

**Conditions for this to be exact, and where each one comes from:**

| requirement | holds in the frozen model because | checked by |
|---|---|---|
| ion paths independent of a w translation | the 2100 V field is 2D (w-invariant); the replay asserts no longitudinal motion | model property (stated) |
| weighting potential of strip j depends on (u, v, η − jP) only | the LUT is built and queried in the phase coordinate η, with identical, contiguous, infinitely many strips | **E1**: whole export relabelled by ±1 strip, frozen replay, overlap vs frozen S0 (toy: 5.9e-16) |
| kernel response depends on D − oP only | the engine indexes `jd − d_sub·o` on a D axis of spacing P/d_sub, so a shift by m strips is exactly m·d_sub nodes | K2A_local vs base-centred K2A on an in-domain muon (toy: 2.3e-15) |
| linear superposition | Shockley–Ramo is linear in charge; clouds are independent | absolute-label sum (`k2a_local_sum.py`) |
| charge conservation | the frozen replay checks per group that every ion lands on a window strip and the census matches; groups partition the primaries | `local_sum.json`: Σ collected = exact × e (toy: 1.3e-15) |
| truncation | every primary still sees ±H strips, now centred on its own signal; the tail beyond is the whole-plane sum-rule tail (7.7e-4 for 13 strips, baseline LUT check) | E2: S0_local vs base-centred S0 on in-domain muons |

**Finite-geometry limitations.** Recentering inherits these from the frozen model; it does not create or remove
them, but it would make far-w events computable where the model is less trustworthy.
- infinite strip length;
- no PCB edges, frames, spacers or connectors;
- no interstrip gaps;
- the 9-wire analytic field patch in u (clouds on |k| = 4 are flagged);
- the LUT phase range ±12 mm, which local windows keep satisfied: |η| ≲ 0.85 + 5·1.7 + 1.35 mm ≈ 10.7 mm;
- the w-invariance itself.

Near real strip ends or board edges, neither S0 nor K2A is valid. That needs the final mechanical geometry, and the
preflight's |w| statistics only proxy for it.

**Output format consequence.** Per plane, an event produces a sparse set of absolute strips: the union of the local
windows. Digitization needs exactly this (channel IDs), and it removes the arbitrary entry-centred window. Engine
cost scales with clouds × (2H+1), unchanged.

**Effect on the Green table.** With recentering, the D axis no longer grows with the event population. It only has to
cover |D′| ≤ P/2 plus the primary's own spread plus ±H strips. That is the current r14t72edge5 D range or less. The
extension then reduces to the missing wires and radial nodes: size ∝ wires × r-nodes × D-nodes. The re-estimate goes
into `recovery.md` and `extension_estimate.json` after review.

## 3. Coverage recovered by recentering alone (input metadata only)

`k2a_recenter_recovery.py` re-runs the engine's own support tests per cloud after recentering, on the 977 non-sealed
events (sealed ids are refused). For each out-of-domain event it gives the **remaining** causes:
- WIRE: needs Green paths for that wire;
- R_EFF_LOW / R_EFF_HIGH: needs radial nodes;
- T_WINDOW: needs a longer window or tau axis;
- D_U/V_RECENTERED: a primary spread wider than the window.

It also reports:
- recovered event count and charge, by channel;
- union strip widths;
- charge whose local strip lies > 5 strips from the entry strip.

I cannot give the numbers without your files. Expectation, to be confirmed: almost all D-only failures are recovered,
leaving the WIRE and radial cases you listed.

## 4. Near-wire radial treatment

- The kernel grid is log-spaced, from r_lo just above 15 µm (half the innermost observed clearance) to r_hi.
- One log step inward from r_lo lands at r_lo/step, below 15 µm, **inside the anode** (toy grid: 11.66 µm). The
  extension estimate in 47ac8ef had this defect.
- Fix, in this commit: inner nodes geometric in the **clearance** c = r − 15 µm,
  c_k = c_min (c_0/c_min)^(k/n),
  - c_0 is the existing innermost clearance;
  - c_min is half the smallest required clearance, as `s_grid_spec.py` already does for the outer rule;
  - n is the number of steps needed at the existing ratio.
- Every node stays in gas. Engine interpolation (in log r over any monotone node set) is unaffected.
- Physically the limit is regular: an ion born at clearance c → 0 starts on the anode surface, where every strip's
  weighting potential is 0 (the anode is a conductor), and the RKF path from the surface is well defined.
- New grid paths are needed only for the added (wire, r) start points.
- Outward extension keeps the existing log step; nothing bounds it except the cell (r ≲ 1.5 mm).
- **R_EFF_LOW clouds need inspection before any extension.** r_eff below a node 0.0x µm above the wire means
  essentially all of the cloud's ions are at the surface. Check their row counts and clearances in
  `cloud_violations.csv` and `per_event_recenter.csv`.
- **R_EFF_HIGH clouds may break the single-r_eff K2A approximation.** They are wide or far clouds, so they need a
  G-vs-S0 check, not just more nodes.

## 5. Minimal comparison experiment against the frozen S0 weighting response

Tools: `recenter/` (group export, local S0 sum, recentered K2A, comparison) and the runner modes `recenter-refs` and
`recenter-far`.

| step | events | what | pass condition (proposed) | cost |
|---|---|---|---|---|
| E1 symmetry | the five references (or a subset) | whole T1 export relabelled by m = ±1, frozen replay, overlap vs frozen S0 | rel L2 ≤ 1e-10 on every overlapping strip | 2 frozen replays per muon (16–275 s each) |
| E2a recentering, in-domain | the five references | per-primary groups → frozen replay + electron step per group → S0_local; K2A_local with the unchanged engine and kernel | census exact (≤ 1e-9); S0_local vs S0 on the base strips: the effect of the window choice, reported; K2A_local vs K2A base: ≤ 1e-12 if all m = 0, otherwise reported; K2A_local vs S0_local: within the five-muon K2A vs S0 range (0.15–0.20 % L2) | ≈ one frozen replay per muon |
| E2b recentering, far D | 2–4 non-sealed far-D picks of the supplementary selection | frozen export + T1 (new files), then E2. The frozen base-centred S0 is attempted and its refusal recorded. | census exact; K2A_local vs S0_local not worse than on the in-domain references; anything else is a finding, not tuning | per muon: export (to be measured) + one replay |

No kernel is rebuilt, no frozen file changes, and nothing is sealed. Toy results (coarse toy grid, not physics):
- E1: 5.9e-16;
- S0_local equals S0 exactly for a single group;
- K2A_local equals K2A base to 2.3e-15;
- far toy event: base-centred S0 and engine both refuse, while grouped S0_local conserves the census
  (rel 1.3e-15) and K2A_local vs S0_local error is the same as on the unshifted toy (6.5 %; coarse toy kernel).

## 6. Supplementary development selection (non-sealed φ/ω)

`k2a_supplementary_selection.py` is deterministic (fixed hash ties) and uses input metadata only.
- Candidates are roles `pool` (optionally `previously_inspected`) with channel φ or ω. The sealed set and the 20 J/ψ
  development muons are removed by role before any feature is read.
- Per channel (default 8), it picks in order: typical (median), p95 clouds, max clouds, multiwire, far-D (recentering
  needed), one per cause still out of domain after recentering, and patch-edge if any.
- It then fills by max-min distance in standardized features: log clouds, log charge, wires ≥ 5 %, dominant share,
  |θ_u|, |θ_w|, D span, union strips, K256 fraction, plus station one-hot.
- Outputs are CSV, md and json with input hashes and SHA256SUMS.
- The far-D picks feed E2b. The full set later gets S0/G references (grouped S0 where far) for priority 1.

## 7. Provenance of your coverage.md and the K256 stored fraction

No coverage.md file reached this conversation, so I could only use your summary.
- The K256 reservoir line was added to the domain check in **78e6f1f**. A coverage.md without it was written by the
  47ac8ef domain check.
- `run_k2a_campaign.sh` records HEAD at the domain stage in `PF/SOURCE_GIT_SHA.txt`. Check it with
  `cat $PF/SOURCE_GIT_SHA.txt` and `grep -c "K256 reservoir" $PF/domain/coverage.md`.
- `recover` derives the stored-fraction statistics from the existing `*.event.json`. It needs no SHA pass and no
  conversion. It also prints the provenance (section 5 of recovery.md).

## 8. WSL commands (read-only on everything that exists)

```bash
cd ~/NA6PRoot && git fetch origin && git checkout dev/na60-2100v-validation-20261123 \
  && git pull --ff-only origin dev/na60-2100v-validation-20261123 && git log --oneline -1
C=test/mwpc-stage-c/cathode-ion-reference/k2a-campaign
export MWPC_K2A_REFERENCE_RUN="$PWD/test_runs/mwpc_stage_c/s_ladder_2100/run_20261124_124740" MWPC_K2A_TAG=r14t72edge5
PF=$(ls -d test_runs/mwpc_stage_c/k2a_campaign_2100/preflight_* | sort | tail -1); echo "$PF"
cat "$PF/SOURCE_GIT_SHA.txt"; grep -c "K256 reservoir" "$PF/domain/coverage.md"
# A. recovery estimate + radial check + K256 fraction (seconds, metadata only)
bash $C/run_k2a_campaign.sh recover "$PF"
RD=$(ls -d "$PF"/recenter_*/ | sort | tail -1)
# B. supplementary phi/omega development selection (metadata only)
bash $C/run_k2a_campaign.sh supplementary "$PF" "$RD" 8
```
Return `$RD/recovery.md`, `$RD/recovery.json` and `$PF/supplementary_*/supplementary_dev.md`.

The experiment (E1/E2a, then E2b) starts only after you review those. It runs frozen replays, about one to three per
reference muon. Its commands are `run_k2a_campaign.sh recenter-refs` and `recenter-far "$PF" <far-D indices>`, with
`MWPC_SPECTRAL_LUT` and `MWPC_NUMW_DIR` (and `MWPC_EXPORT_BIN`, `MWPC_ION_MOBILITY` for far).

## 9. Hardened comparison runner (2026-11-25, after the metadata review)

Metadata review results (yours):
- recentering recovers 149 of 197 out-of-domain non-sealed events, so coverage goes from 780/977 to 929/977;
- no D violation remains;
- 46 R_EFF_HIGH and 2 WIRE events are left;
- the K256 reservoir stores about 0.6 % of the exact births (median);
- the φ/ω supplementary set is approved.

The runner `recenter/k2a_recenter_experiment.sh` was hardened before any real run.

**Status per step, never implied by "the command finished".** Each step writes a line to `muon_*/steps.tsv`.
- E1 (each of the two relabellings) is **PASS / FAIL / INCONCLUSIVE**:
  - a frozen-replay refusal ("selected channels omit the landing strip" or "query outside declared strip phase
    range") is INCONCLUSIVE, with the reason recorded;
  - any other replay error is FAIL.
- E1 acceptance (`thresholds.json` → `e1`), per plane and term:
  - unshaped, cumulative and shaped (100/200 ns) rel L2 ≤ 1e-10;
  - window charge and U/V ratio ≤ 1e-10;
  - max per-strip charge difference / lead ≤ 1e-10;
  - same lead strip;
  - **exactly** 2H+1−|m| common strips; fewer is INCONCLUSIVE, more is FAIL.
- E2 comparisons are **OK / WARN / FAIL / INCONCLUSIVE**:
  - documented warning thresholds (`thresholds.json` → `e2`) produce WARN;
  - numerical breakage produces FAIL: nonfinite values, inconsistent time grids, an identity that should hold
    failing, or a census error;
  - nothing to compare produces INCONCLUSIVE: missing files, fewer than 3 common strips, zero reference charge.

**Cases that cannot pass silently.**
- Zero common strips, NaN/inf and missing files are explicit non-OK statuses (tested).
- `worst_L2_rel` is `null` unless every value is finite.
- The exit code is 0 only for PASS/OK.

**E2 comparisons, each with full diagnostics** (rel L2 unshaped/cumulative/shaped, ΔQ/Q, Δ(U/V), lead strip, max
per-strip ΔQ/lead, charge on strips present on only one side = missing strip tails):

| comparison | kind |
|---|---|
| K2A_local vs S0_local | k2a_vs_s0 |
| S0_local vs frozen S0 | identity when every primary has m = 0, otherwise the window-choice diagnostic |
| K2A_local vs frozen S0 | k2a_vs_s0 |
| K2A base vs frozen S0 | k2a_vs_s0, the validated baseline for context |
| K2A_local vs K2A base | identity or window-choice |

The S0_local census and the K2A_local evaluation are steps of their own.

**Aggregate.** `k2a_recenter_aggregate.py` writes `results.json` and `results.md`. The overall result is PASS only if,
for every muon:
- E1 PASSes;
- E2 has no FAIL or INCONCLUSIVE (WARN is listed for inspection);
- every reference file read is byte-identical afterwards (sha256 before/after of T1_export, S0, cloud_moments.csv and
  K2Akern_TAG).

All outputs go to a new `recenter_experiment_<timestamp>/`; the runner refuses an output path inside the reference run.

**Environment.** Variables come from the validated setup (`source ~/t_ladder_deps.env`) and are **verified, not
assumed**:
- the LUT sha256 must equal both the Green table's `lut_sha256` and the `lut_sha256` in
  `provenance_prepare_<M>.txt` of each reference muon;
- the replay code sha256 must equal its `replay_sha256`;
- the replay module must import from `MWPC_NUMW_DIR`;
- the kernel id must equal the validated `engine_TAG_*/kernel_id.txt`;
- far mode also requires the exporter and mobility sha256 to equal the references' provenance, and the exporter must
  start and print its usage (ROOT/Garfield runtime present).

ROOT is not needed for `refs`: replay, electron step and engine are Python/numpy.

Synthetic test `tests/test_recenter_synthetic.sh` (toy run):
- a wrong LUT aborts;
- E1 PASSes at 6e-16;
- the single-group identities PASS;
- the four negative compare cases (disjoint strips, NaN, a 1e-6 perturbation under E1, missing files) all end
  non-OK;
- the far toy event behaves as before;
- references are unchanged.

**Order of the real runs:**
1. `envcheck` and `refs` on development references 0000 and 0858, then inspect.
2. Only then `far` on ω 0560 and φ 0997 (both far-D events recovered by recentering).
3. φ 0133 is deferred: it still has R_EFF_HIGH and cannot be evaluated with the unchanged kernel.

No Green table is enlarged and the 2100 V field is unchanged.

### 9.1 Final corrections before the reference runs (review of 9af160a)

**Completeness in the aggregate.** `k2a_recenter_aggregate.py` now requires the content of the mode. Missing pieces
give INCONCLUSIVE, never PASS.
- **refs:**
  - exactly `E1_shift_1` and `E1_shift_-1`, both PASS;
  - all eight E2 steps: groups, S0_local census, K2A_local, K2Alocal_vs_S0local, S0local_vs_S0, K2Alocal_vs_S0,
    K2Abase_vs_S0, K2Alocal_vs_K2Abase;
  - reference checksum lists present, non-empty and identical.
- **far:**
  - no E1 steps;
  - `base_centred_S0` present (a refusal, INCONCLUSIVE, is acceptable);
  - the four E2 core steps, plus S0local_vs_S0 and K2Alocal_vs_S0 when the base-centred S0 ran.
- **Every comparison step** must point to its `cmp_<step>.json`. A missing report is INCONCLUSIVE. A report status
  that disagrees with the step status is FAIL.
- **Every muon** in `env_check.json` needs an output directory, and outputs not listed there are flagged.
- **Duplicate steps** give INCONCLUSIVE.
- **Overall:** FAIL > INCONCLUSIVE > PASS. Exit code 0 only for PASS (1 FAIL, 3 INCONCLUSIVE).
- **Synthetic negatives** (`tests/test_recenter_synthetic.sh`, 13 variants on copies of a real toy experiment). Each
  gives the expected FAIL or INCONCLUSIVE, and the unmodified copy PASSes:
  - no E1, or only one E1 relabelling;
  - one E2 comparison missing, or no E2 at all;
  - empty `steps.tsv`;
  - a deleted report;
  - a status mismatch;
  - a missing or changed reference checksum list;
  - a muon missing;
  - environment not PASS;
  - wrong mode.

**Bounded-memory group export.** `k2a_group_export.py` previously stacked every selected path in memory; in shift
mode it held the whole export several times.
- It now reads the source with explicit `pread`s (no mmap) and writes each ion path's block as it goes.
- The per-block arithmetic is unchanged (copy as <f8, column 3 −= P·m_U, column 4 −= P·m_V), and the (t, u, v)
  checksum is accumulated over the same bytes.
- `tests/test_group_export_equivalence.py` compares against the 9af160a implementation on synthetic exports with
  strip columns spread from −14 to +21 mm and integer per-primary census. Every output file is **byte-identical** in
  both modes (vertices, ions.csv, complete.txt, groups.json; 70–76 files).
- Peak RSS:

  | synthetic export | mode | new | old (9af160a) |
  |---|---|---|---|
  | 57 MiB | primary | 103 MiB | 164 MiB |
  | 57 MiB | shift | 105 MiB | 283 MiB |
  | 229 MiB | primary | 123 MiB | 559 MiB |
  | 229 MiB | shift | 131 MiB | 1038 MiB |

  The remaining growth is the ions.csv row list, not vertices.
- Disk: about 3 × T1_export per reference muon (two relabelled copies plus the groups). The runner checks free
  space before starting a muon.
- The frozen replay still reads each export the way it did for the original S0 runs; it is unchanged.

### 9.2 E2b (far-D) analysis, after the 0858 reference result

Your reference results, for context:
- 0858: E1 PASS (≈ 2e-16); two groups, g_0_0 and g_-1_0; census 3.6e-14;
- S0_local vs S0 and K2A_local vs K2A base: both 1.94e-5;
- K2A_local vs S0_local: unshaped ≤ 0.531 %, shaped 0.19–0.25 %. The 0.5 % unshaped warnings are also present in the
  original K2A, so they do not come from recentering;
- integrity PASS.

**Far mode keeps the four effects separate.** Each muon gets `muon_*/far_analysis.{json,md}` and a separate
`truncation.json`:

| # | effect | how it is measured | reference |
|---|---|---|---|
| 1 | K2A approximation error | K2A_local vs S0_local: same local windows, same frozen weighting, same paths. Total, **and per recentering group** (`K2A_groups/<g>` vs `S0groups/<g>`), so a far group's error cannot hide behind a near one | the frozen S0 of the same sub-export |
| 2 | recovered signal outside the original entry-centred ±5 strips | S0_local and K2A_local per plane: signed charge inside/outside base ±5; Σ\|q_strip\| outside as a fraction; collected-ion (landing) charge outside, from the frozen replay's channels.csv; lead strip and its distance from the entry strip | the frozen base-centred S0 cannot provide it (it refuses), so this is signal the old 11-strip definition lost |
| 3 | finite local-window truncation | exact final-state sum rule at each ion's birth with the frozen LUT and the frozen whole-plane potential (`fit(0, 24, 0)`, the replay's own check): charge beyond each primary's ±5 = Σ q (Σ_win φ − Φ). The ±6 ring is reported where the LUT phase range allows, with out-of-range ions counted, never extrapolated. Also the frozen replay's subset strip-sum gap per group | integrated (final) charge only; transient signal beyond the window is bounded only by the ±6 ring and is stated as such |
| 4 | unsupported geometry / weighting assumptions | shared by S0 and K2A: \|u\| vs the LUT cross-section 13.5 mm, local-window LUT phase margin vs 12 mm, charge on \|k\| = 4 patch wires, longitudinal extent over which infinite strips and a w-invariant field are assumed (strip ends, PCB edges, frames and gaps are not modelled). K2A-specific: per-cloud rms of D and of ln r around the cloud means the kernel uses | flags (WARN) per `thresholds.json` → `e2b`; the e1/e2 values are unchanged |

**The base-centred S0 refusal is expected.** It is recorded as `EXPECTED_REFUSAL`, which is not a failure. That status
is valid only for `base_centred_S0` in far mode; anywhere else it is FAIL. A replay error of any other kind is FAIL.
Group replays, the census, local K2A, truncation and analysis failures are never suppressed.

**Far-mode completeness** (aggregate):
- `base_centred_S0` (OK or EXPECTED_REFUSAL);
- E2 groups, census, K2A_local, K2Alocal_vs_S0local, E2_truncation, E2_far_analysis, with reports whose status
  matches the step;
- the two base comparisons, when the base S0 ran;
- archive sha256 before == after.

**Environment for far mode:**
- LUT and replay as before;
- the exporter and mobility sha256 must equal those recorded in **every** `provenance_prepare_*.txt` of the reference
  run;
- the exporter must start and print its usage;
- the archive slot must pass its own `SHA256SUMS`;
- the slot must not be sealed (checked from the preflight `slots.csv`);
- a minimum free-disk check.

A wrong mobility file and a sealed slot are refused in the synthetic test.

**Afterwards (not started).**
1. Compact-moment accumulation inside the avalanche generator, from all births before K256 sampling.
2. Moving the per-primary recentering and absolute-strip accumulation into the compiled evaluator (`kernel_eval.cc`):
   output as a sparse (first strip, n strips, samples) per plane.
