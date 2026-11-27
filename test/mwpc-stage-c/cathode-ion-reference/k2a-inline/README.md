# K2A inline full-birth moments (production integration, milestone 1)

Goal: form the K2A cloud moments of the **complete** positive-ion birth population inside the Garfield++ avalanche
path, before K256 reservoir sampling, without `ion_births.root`, trajectory exports or CSV in production.

Unchanged: the frozen generator source and binary, Garfield physics, Heed, the random-number streams, the reservoir,
the archive writer, the 2100 V field, the cathode geometry, the weighting references, the Green table, the kernels,
the sealed sample and the provisional electronics.

## Integration point

Frozen generator: `test/mwpc-avalanche/phaseB_native_muon_p3_noFW_penning057_cathode_ion_reference_v1.cxx`
(sha256 `1436a715…f874f`, untouched).

- `OnIonisation` → `BirthAccumulator::Add(u, v, w, t)` sees every birth before the reservoir. It already calls the
  optional `fullBirthSink(u, v, w, t, ordinal)`.
- Per Heed primary `i`, `RunGeometry` creates the `BirthAccumulator`. In full mode it sets the sink to
  `Archive::Ion`, runs `AvalancheElectron`, then reports `births.total`.

The new copy `phaseB_native_muon_p3_noFW_penning057_k2a_inline_v1.cxx` adds 8 marked `K2A-INLINE` blocks:

| block | where | what |
|---|---|---|
| header | top | description |
| include | after the archive include | `MomentAccumulator.hh`, hook state |
| sink | `RunGeometry`, before `gBirth=&births` | `BeginPrimary(i)`; the existing sink is **wrapped**, not replaced: archive `Ion()` first, then `AddCallbackBirth(u, v, w, t, ordinal)` |
| end_primary | after the frozen outside-patch abort | `EndPrimary(births.total)` (the generator's own census) |
| options | `main`, before the frozen parser | `--k2a-moments-out PREFIX`, `--max-heed-primaries N` (stripped from argv) |
| frame | `main`, per crossing | `MuonFrame{.1*a[0]-phase, a[0], a[2]}`: the same double expressions as `Archive::Muon` |
| subset | the `RunGeometry` call | onlyIndices = 0..N-1 (the frozen call is kept verbatim as `//@frozen`); `Finish` + write |
| exit | end of `main` | exit code 4 if the moments carry a blocking status |

The allowed-diff checker `check_inline_generator_diff.py` removes the blocks, restores the `//@frozen` lines, and
requires the result to be byte-identical to the frozen source with the recorded sha256. It also requires exactly
these 8 block names and no K2A code outside them.

**`--max-heed-primaries N` (TEST ONLY).**
- Heed is generated in full and unchanged.
- Only primaries 0..N-1 are avalanched, through the frozen `onlyIndices` path. Their original indices `i`, and hence
  their Garfield seeds `1 + (heed_seed + (i+1)·104729) mod 2e9`, are kept.

**Outside-patch semantics are unchanged.**
- In full-archive mode, `Archive::Ion` aborts at once.
- Otherwise, `birthsOutsidePatch` aborts after the avalanche.
- No K2A output is written in either case.

## Interface (`MomentAccumulator.hh`, header-only, no ROOT/Garfield)

```cpp
K2AInline::MomentAccumulator acc(K2AInline::MuonFrame{patch_shift_cm, u_in_mm, w_in_mm});
acc.BeginPrimary(original_heed_index);                    // strictly increasing
acc.AddCallbackBirth(u_cm, v_cm, w_cm, t_ns, ordinal);    // or AddBirth(..., weight)
acc.EndPrimary(callback_births);                          // generator census
K2AInline::Event ev = acc.Finish();                       // clouds (primary, wire), primary records, status
K2AInline::WriteOutputs(prefix, ev);                      // moments file only if ev.UsableForCathodeSignal()
```

**Definitions.** The definitions are exactly those of `spatial-ablation/BirthMoments.hh`, expression by expression,
in the frozen strip frame: tan 0.1, 0.17 cm pitch, base strips from the muon entry.
- A cloud is (primary, wire = `rint(u_mm/3)`), with ties rounded to **even**. The generator's census rule
  `llround(u/0.3)` rounds half away from zero. The two differ only exactly on a half-cell boundary; T0 tests all 8
  boundaries.

**Timing.**
- `t_ns` is the T1 time, the primary's weight-mean birth time. It is bitwise equal to the batch code's per-primary
  mean.
- `t_rms_ns` is the collapsed spread, 0 by definition.
- The raw spread is kept and never collapsed: `t_raw_mean_ns` and `t_raw_rms_ns` per cloud, and `t_raw_rms_ns` per
  primary.

**Electron step.** No extra data is needed. The validated K2A electron term is the instantaneous step at the T1 time on
the same moments.

**Statuses.** A blocking status suppresses the engine-facing `PREFIX.moments.csv`; only diagnostics are written.
- `INSIDE_ANODE`: r ≤ 15 µm. Blocking.
- `NONFINITE`: a non-finite coordinate or time, or weight ≤ 0. Blocking.
- `OUTSIDE_PATCH`: \|k\| > 4. Blocking.
- `CENSUS_MISMATCH`. Blocking.
- `SEQUENCE_ERROR`. Blocking.
- `ZERO_ION_PRIMARY`: informational only.

Engine-side domain statuses (wire, r_eff, D, time window) belong to the evaluator (milestone 2).

**Memory** is one per-wire accumulator for the current primary plus one record per cloud and per primary. Nothing per
birth is kept. `accumulator_state_bytes_peak` is reported.

**Outputs.**
- `PREFIX.inline.json`: status, census, frame, memory.
- `PREFIX.clouds.csv` and `PREFIX.primaries.csv`: all records, with statuses.
- `PREFIX.moments.csv` and `.frame.json`, only if the event is usable. These use the `BirthMoments::Write` schema, so
  existing K2A tools read them unchanged.

## Tests

| test | where | what |
|---|---|---|
| T0 | `tests/moment_accumulator_t0.cc` (no deps) | incremental vs `BirthMoments::Compute` on 430 k synthetic births (unit and K256 weights, 100 primaries incl. zero-ion ones), half-cell ties, all statuses, writer |
| synthetic E2E | `tests/test_t1_synthetic.sh` | frozen + inline sources compiled against **stub** Garfield/ROOT (`tests/stubs`, toy avalanche, NOT physics); allowed-diff negatives; frozen vs inline outputs byte-identical; full T1 runner; negatives |
| **T1** | `run_t1_inline.sh SLOT OUT [N]` | real 2100 V avalanche of one non-sealed muon, first N Heed primaries; see below |
| T2 | inside the T1 runner | K256 vs full-birth moments, **diagnostic only** |

**T1 runs** (all from the same build):

| run | mode | moments |
|---|---|---|
| R1 | full archive | in-memory moments, in **one** run |
| R2 | R1 repeated | — |
| R3 | moments only (production path) | yes |
| R0 | R3 without moments | no (memory baseline) |
| R4 | weighted K256 archive | yes |

**Two independent full-population references for R1:**
- (a) `p3_birth_moments` on R1's full archive, through the frozen reader and `BirthMoments.hh`;
- (b) numpy per-ion values with `math.fsum` per cloud, on `k2a_archive_dump` output. The dump reads raw TTrees and
  shares no code with (a).

`k2a_inline_compare.py` checks census, seeds, identity, frame, statuses, moments, timing and charge. Every tolerance is
derived in its header (u = 2^-53; recursive-summation bound with a factor-2 margin plus 4 ulp per libm value).
Bitwise agreement with (a) is reported separately per field.

**Required for PASS:**
- R1 = R2 for moments and archive content;
- R3 = R1 and R4 = R1 for moments;
- R4's reservoir rows are a bitwise subset of R1's births;
- R0 and R3 have identical per-seed physics CSVs.

**Reported separately, never part of the verdict:**
- T2;
- the cross-check against the frozen production K256 archive of the same muon (seeds, census, reservoir rows;
  MATCH/MISMATCH depends on the build environment).

## Observation from T0

With `-O3 -march=native`, the spatial moments differ from the batch code by rounding: FMA contraction is applied
differently to the two inlined copies. With the default x86-64 flags, or `-ffp-contract=off`, they are bitwise equal.
The T1 PASS criterion is the tolerance bound; bitwise agreement is reported, not required.

# Milestone 2: compiled absolute-strip K2A evaluator

Goal: from the in-memory full-birth moments (`K2AInline::Event`), produce sparse U/V strip waveforms on absolute
strip labels in compiled code. Nothing in the production path uses K256, ROOT or Python.

Unchanged: the kernel (`s_kernel.py` directory, identity checked), `spatial-ablation/engine/kernel_eval.cc` (compiled
in, not modified), the 2100 V field, the weighting LUT, and the validated recentering convention.

## API (`K2AStripEvaluator.hh`, `K2AKernelIO.hh`; header-only apart from `kernel_eval.cc`)

```cpp
auto k = K2AStrip::KernelTable::Open(KERNEL_DIR, KERNEL_ID);       // once per process; K.npy mapped read-only
K2AStrip::Evaluator eval(*k, {"K2A", 5, 7000.0, 10.0});            // model, half strips, window ns, dt ns
K2AStrip::Result r = eval.Evaluate(event);                         // event: K2AInline::Event (milestone 1)
// r.usable, r.status; r.U / r.V: strips (ascending absolute), ion and electron (strip-major, fC/ns, 700 bins)
// provenance: r.clouds (primary, wire, m_U, m_V, local D, status), r.shifts (Dbar, m, tie margin), r.groups
```

**Kernel identity.**
- `kernel_id` is the sha256 of the three file hashes (K.npy, KPHI0.npy, kernel.json), the same value as
  `cathode_engine.KernelTable`. A mismatch refuses.
- `transport_id` and `weighting_id` are read from the axes.

**Per-cloud arithmetic.** This is `CathodeResponseEngine._prepare` operation for operation:
- the coefficients, including the hemisphere term `(A - a0)/(s_perp·sm)`;
- log-r interpolation using `searchsorted − 1`, clipped;
- the D index and fraction.

`kernel_eval.cc` is then called once per group, exactly as the Python `compiled` implementation calls it.

**Recentering (E2a/E2b convention of `k2a_recentered_k2a.py`).**
- Each primary gets m = floor(Dbar/P + 1/2), where Dbar is the N-weighted mean D of its clouds. The sum uses numpy's
  pairwise order, so integer shifts reproduce the Python rule exactly.
- Groups are formed by (m_U, m_V) and evaluated in ascending order, with D' = D − mP and base + m.
- Each group's result is added on its absolute labels.
- An explicit shift table, such as an E2b `groups.json`, can be imposed. Any difference from the production rule is
  counted.

**Statuses.**
- Blocking:
  - `INPUT_BLOCKED`: the moment event itself is blocked;
  - `WIRE_NOT_IN_KERNEL`;
  - `R_EFF_BELOW_GRID` and `R_EFF_ABOVE_GRID`;
  - `D_U_OUTSIDE_TABLE` and `D_V_OUTSIDE_TABLE`: the local D' over all ±5 strips;
  - `NONFINITE_MOMENT`;
  - `NONPOSITIVE_CHARGE`;
  - `T_NEGATIVE`;
  - `SHIFT_UNDEFINED`;
  - `ENGINE_INDEX_ERROR`.

  If any is set, there is no waveform at all. The per-cloud statuses and the blocked charge are still returned.
- Informational: `T_AFTER_WINDOW` and `SHIFT_NEAR_TIE`.

Nothing is clamped or dropped.

## Tools

| file | purpose |
|---|---|
| `k2a_strip_eval.cc` | CLI: evaluates a moments CSV, writes replay format + clouds/shifts/groups CSV + `result.json` (kernel verify/load time, cold/warm evaluate, peak RSS, mapped vs `--preload`); `--shifts groups.json`; `--displace-event` (test only) |
| `k2a_strip_compare.py` | waveform-only comparator: labels, time grid, L2, per-strip L2, per-strip and total window charge, U/V ratio; `--label-shift` |
| `run_m2_adapter.sh` | WSL validation on the real kernel with **explicit required coverage**: E2b 0560 and 0997 (offline reproduction, displacement symmetry), the accepted T1 run of 0000, and ≥ 1 further T1-validated multiwire/displaced run; all three benchmark modes checked; anything missing is INCONCLUSIVE |
| `k2a_select_regression.py` | picks a non-sealed (slot, N) for a real multiwire / displaced-primary full-birth run, from the preflight K256 moments; N = last included original primary index + 1 (zero-ion primaries are absent from the moments), census checked against `event.json` |
| `tests/test_strip_eval_synthetic.sh` | toy-kernel tests A–H (SYNTHETIC, not physics) |
| `tests/strip_eval_inmemory_test.cc` | births → accumulator → Event → evaluator in memory, == CSV path bitwise |

**Agreement with the Python engine.** On the toy kernel the adapter matches the Python engine to rel L2 ≤ 7e-15 for
the loop, vector and compiled implementations. It is not bitwise: the hemisphere constants cos(nθ)·sign(sin θ) are
exactly zero analytically, but their computed value is rounding noise of about 1e-15 that depends on the BLAS
summation order. The required tolerance is therefore 1e-12, not bitwise equality.
