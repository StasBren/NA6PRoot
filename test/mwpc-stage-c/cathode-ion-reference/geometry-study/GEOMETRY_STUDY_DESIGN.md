# Controlled cathode strip-geometry study: design and feasibility (2100 V, five development muons)

**PAUSED (2026-11-25):** kept in the repository, not to be run or developed further for now. Current priority: the
large-sample K2A campaign at the existing geometry, `../k2a-campaign/K2A_CAMPAIGN_PLAN.md`.

**Status (2026-11-24).** This commit contains tooling and SYNTHETIC/LUT-level checks only. No real muon has been run under
the alternative geometry. No geometry scan has been run. The sealed 20-muon set, the archives, the frozen
exporter/replay/LUT/Green tables, the electronics stand-in and avalanche generation are all untouched. No frozen file
is edited: the alternative geometry is applied through a run-time overlay (section 4).

**First alternative.** Stereo tan α changes from 0.10 to 0.15 on both planes (U = +A, V = −A). The strip pitch stays
1.7 mm along w, the strips stay contiguous (width = pitch), and everything else is unchanged.
`geometries/stereo015.json`, geometry_id `3ed7d01b832d…`.

## 1. Audit: where each readout parameter enters

| parameter | frozen value | enters (frozen wire branch unless noted) | coupled to the 2100 V drift field / ion transport? | changeable here |
|---|---|---|---|---|
| stereo tan A | 0.1 | `p3_fourier_mfs_mode_pilot.py:13,36,40`: slab mode `exp(−ikAu)`, `q=\|k\|√(1+A²)`. `p3_fourier_strip_reconstruction_pilot.py:19–20` (bare strip, imported A). `p3_spectral_weighting_lut.py:47,83,99,155`: phase η = w − Au, table, meta, lookup. `p3_replay_full_spectral_muon.py:149`: landing strip. Exporter `p3_export_full_spectral_paths.cc`: strip frame (`.1`). Ours: `BirthMoments.hh`, `cathode_engine.moments_from_births`, `s_green_table.py:89` (range check), `t_shaped_compare.py:39`, `t_centroid_diagnostic.py:24`, `run_t_ladder.sh:86` (assert) | **no** (model assumption, section 2) | **yes** (first alternative) |
| strip pitch P (along w) | 1.7 mm | `p3_fourier_strip_reconstruction_pilot.py:13` (“pitch AND width”), `bare` d = P/2. LUT `:38` (strip-indicator Fourier coefficients 2 sin(kP/2)/(kL)), `:98,142,216`. Replay `:67` (offsets w − kP), `:150–151`. Exporter (`.17`). Ours: `BirthMoments.hh`, `moments_from_births`, `s_grid_spec.py:24`, `s_green_table.py:61` (dD = P/d_sub), `t_shaped_compare.py:38` | no (same assumption) | yes (second option, contiguous only) |
| strip width | = P | the frozen model has no separate width; ideal contiguous strips | a gap needs a boundary model for the gap surface (also a field question) | **no** |
| strip indexing | idx = llround((w_in ∓ A u_in)/P), offsets −5…+5 | exporter frame (`upper_index`/`lower_index`); replay `--half-strips 5`; engine `half_strips=5` | no | follows A, P (recomputed) |
| U/V plane assignment | U on +v cathode (+A), V on −v (−A, mirrored x → −x) | one LUT serves both planes by mirror symmetry (`query()`, `s_green_table`) | no | both planes change together; a single-plane change would need a second LUT (not built) |
| anode pitch | 3 mm | mode pilot `PITCH`, `CENTERS`; LUT `PITCH`, `UMAX`; exporter field; wire index rint(u/3) | **yes**: field and paths | no |
| anode radius | 15 µm | mode pilot `R`; LUT near table; replay `path_for` gas check | **yes** | no |
| anode–cathode gap | 3 mm | mode pilot `GAP`, `H`; LUT bulk grid; replay endpoint check | **yes** | no |
| number of wires | 9 | mode pilot `CENTERS`; LUT `j` clip 0…4 | **yes** | no |
| LUT numerics | L 48 mm, kmax 6/mm, N 24, phase ±12 mm, 241 phases, bulk 0.15 mm, near 32×48 | `p3_spectral_weighting_lut.py:22–24`, build defaults | no | **kept fixed** (only one physical parameter changes) |
| Green / kernel numerics | r14 × θ72, d_sub 10, τ 200 nodes | `run_s_ladder.sh`, `s_green_table.py` | grid paths = transport | kept; grid paths **reused** |

`weighting_geometry.py` refuses any geometry file that changes a transport-coupled row. It also refuses width ≠ pitch.

**Hazard found in the frozen code.** `WeightingLUT.__init__` checks `strip_mm` but not `stereo_tan`, so a LUT built
for another stereo angle loads silently. Every geometry-study tool now checks the table meta (stereo, strip, anode
pitch, radius, gap, model) against the geometry before use (`check_table`). The synthetic pipeline test confirms
that the alternative geometry refuses the baseline LUT.

## 2. Model assumptions (explicit)

The drift field and the Green-grid trajectories are reused unchanged. That is valid only because of the following
model assumptions. They are properties of the frozen model, not of a real detector.

1. **Cathodes are continuous equipotentials.** The 2100 V drift field is the frozen 9-wire analytic field with the
   cathode planes at ±3 mm as uniform grounded planes. Segmentation does not appear in the field model at all.
   Changing stereo angle or contiguous pitch therefore leaves the field and every RKF trajectory unchanged by
   construction. Physically this requires all strips to be at the same DC potential (front-end inputs at virtual
   ground, no strip biased differently) and interstrip gaps negligible compared with the 3 mm gap.
2. **Ideal contiguous strips, width = pitch.** There is no gap region, so no gap boundary condition is needed for
   the weighting field. With real gaps both fields change near the cathode: the drift field over about one gap
   width, and the weighting field through the boundary condition. That case is **not** representable here and is
   not attempted.
3. **2D field, infinitely long strips, B = 0.** Ions do not move along w (the replay checks this). A change of
   strip frame is therefore a constant shift of the strip-relative coordinate per plane (section 4.2). No ion is
   moved or re-drifted.
4. **Quasi-static Shockley–Ramo with an ideal (zero-impedance) readout.** There is no strip-to-strip cross-talk or
   capacitive network, and the electronics response is the unchanged PROVISIONAL CR-RC³ stand-in.
5. **Pitch is defined along w (the wire direction) at fixed u**, as in the frozen model. At fixed P the
   perpendicular strip pitch is P/√(1+A²): 1.6916 mm at A = 0.10 and 1.6812 mm at A = 0.15 (−0.6 %). The first
   alternative therefore changes the angle at fixed w-pitch, not at fixed perpendicular pitch. This is stated so
   that the result is not misread as a pure rotation.
6. **Avalanche positions and charges are the archived ones** (T1 timing, unchanged), and the same electronics are
   used. Thresholds in the physics comparison are absolute fC-equivalent values from the baseline S0 lead peak, so a
   geometry change cannot move the threshold.

Under assumptions 1–3, the weighting potential is the only physics that changes. Changing tan α alters three things:
the transverse wavenumber of every Fourier mode (√(1+A²)), the phase η = w − Au, and the strip frame (base indices
and strip-relative D). The field, the paths, the Green grid start points, the wires, the radii and the cloud
geometry (r, θ harmonics, T1 times) are identical. This is verified on synthetic data: the transport-side moments
of the C++ core are byte-identical between geometries, and only D_U/D_V move.

## 3. Why stereo first

- It is weighting-only under the stated assumptions, and every hard-coded A is reached by the overlay or by a new
  option.
- It has a direct readout observable: the stereo coordinate u = P (c_V − c_U)/(2A) has a lever arm ∝ A. At
  A = 0.15 the centroid error maps into u 1.5× less, while charge sharing per plane changes only through the
  phase/wavenumber.
- Feasibility of the numerical ranges, checked before choosing:
  - LUT phase range |η| ≤ 12 mm: replay queries reach about 0.85 + 5·1.7 + A·13.5 ≈ 11.4 mm at A = 0.15.
  - Green D range, `s_green_table` requires d_max + A·max|u| ≤ 12 mm: the margin shrinks by 0.05·max|u| (≈ 0.3 mm
    for |u| ≤ 6 mm). This is a hard check at run time. If it fails, the remedy is a smaller window (half-strips 4)
    as a separately tagged run, not a change of the LUT phase range.
- Pitch-only is the natural second step. It changes the strip frame and the D grid (dD = P/d_sub), and it changes
  the readout window in mm at a fixed half-strip count, so it is a larger change to the comparison.

## 4. Implementation (files)

All files are new under `geometry-study/` unless marked *(changed)*. Each changed file has defaults equal to the
frozen values, with the regression checked below.

| file | role |
|---|---|
| `weighting_geometry.py` | Geometry record (json → dataclass), `geometry_id` = sha256 of the canonical physics json, transport-coupled constants asserted. `apply()` imports the four frozen modules under their own names and sets `A`/`P` once per process; a second geometry in the same process is refused. `check_table()`, frozen-module sha256. |
| `geometries/baseline.json`, `geometries/stereo015.json` | baseline (no-op) and the first alternative |
| `geom_run.py` | runs the unchanged frozen LUT builder, frozen replay, electron-step tool or Green-table tool under a geometry; checks the LUT meta first; writes `geometry.json` sidecars (geometry id, LUT sha, frozen module sha) |
| `geom_restrip_export.py` | re-expresses an archived export in the new strip frame. Vertex columns 3–4 shift by a constant per plane; `upper/lower_index` are recomputed with the exporter formula; columns t, u, v are byte-identical (hash recorded and checked); the source is read-only |
| `geom_lut_check.py` | independent LUT checks: interpolation vs direct spectral solution, truncation vs refined (L 96, N 32), bounds, anode/strip boundary values, partition of unity → whole-plane potential (segmentation-independent invariant) |
| `geom_physics_compare.py` | same muons, two geometries: U/V charge, U/V ratio, lead peak, threshold time, multiplicity, second/lead, reconstructed w/u vs the charge-weighted ion-birth truth; for S0 (exact) and K2A (fast). Reports the fast-model error on the geometry effect |
| `geom_manifest.py` | identities: geometry_id, weighting_id (must differ), transport_id (must be identical), baseline LUT/exports intact (sha recomputed), kernel ids |
| `run_geometry_study.sh` | runner: selftest, overlay-identity, check, lut, prepare, green, closure, kernel, births, engine, report, hashes |
| `tests/test_restrip.py`, `tests/test_core_geometry.sh`, `tests/test_pipeline_geometry.sh` | SYNTHETIC tests |
| `spatial-ablation/BirthMoments.hh` *(changed)* | `StripGeometry{stereo_tan, pitch_cm}`, defaults = frozen; same double arithmetic |
| `spatial-ablation/p3_birth_moments.cc` *(changed)* | `--stereo-tan`, `--strip-pitch-mm` (single and batch); writes `<out>.frame.json` (patch_shift, u_in, w_in, bases, geometry) — the frame needed to restripe an export |
| `spatial-ablation/tests/birth_moments_core_test.cc` *(changed)* | optional geometry arguments |
| `spatial-ablation/engine/cathode_engine.py` *(changed)* | `moments_from_births(..., stereo_tan, strip_pitch_mm)` |
| `spatial-ablation/engine/engine_cli.py` *(changed)* | refuses moments whose `.frame.json` geometry differs from the kernel's |
| `spatial-ablation/s_grid_spec.py` *(changed)* | `--strip-pitch-mm` (D range) |
| `timing-ablation/t_shaped_compare.py` *(changed)* | `--strip-pitch-mm`, `--stereo-tan` for the stereo reconstruction; recorded in the summary json |

### 4.1 Reuse of ion transport: what is valid

| product | reused for the alternative? | why |
|---|---|---|
| archived births, avalanche, T1 timing | yes | not geometry dependent |
| T1 export paths (t, u, v) | yes, unchanged bytes | field unchanged (assumption 1); restriping changes only the strip columns |
| Green-grid RKF paths `green_paths_TAG` | **yes** (same files, `transport_id` identical) | start points (wire, r, θ) come from u, v only. `green` stops if the grid spec (wires, radii) recomputed from the restriped exports differs from the baseline grid spec |
| weighting LUT | **no**: new, own `weighting_id` | A enters the mode solution and the phase |
| Green table, PHI0, kernels | **no**: rebuilt from the reused paths with the new LUT | they contain φ |
| S0 reference | **no**: recomputed, exact frozen replay, new LUT, restriped export | independent reference |

### 4.2 Strip frame algebra

Exporter, in cm (mm out): idx_U = llround((w_in − A u_in)/P), wU = 10 (w + w_in − A·patch − P·idx_U). The V plane
uses +A. Re-striping from (A0, P0) to (A1, P1) adds the constant ΔU = 10 [(A0 − A1)·patch − (P1·idx_U′ − P0·idx_U)],
and ΔV likewise with the opposite A sign. This needs only the muon frame (patch_shift_cm, u_in_mm, w_in_mm). The new
`.frame.json` from `p3_birth_moments` provides it, read from the archive by the frozen reader. The tool cross-checks
that the frozen formula reproduces the export's indices before writing anything.

## 5. Results so far (labelled)

**LUT level: real numerical builds, no muon data.** They ran on the 2-core sandbox, with the frozen builder at
default settings.

| check | frozen-builder default LUT (baseline) | overlay, stereo 0.15 |
|---|---|---|
| sha256 | `ed458136e3e2…` (built here from the frozen code) | `a185886df7f4…` |
| baseline overlay rebuild | sha256 **identical** `ed458136e3e2…`, arrays and meta bit-identical | — |
| LUT − direct spectral solution, same settings (interpolation), max \|Δφ\| | 1.51e-5 | 1.51e-5 |
| direct (L 48, N 24) − refined (L 96, N 32) (truncation), max \|Δφ\| | 3.56e-8 | 3.54e-8 |
| φ range at 251 held-out points × 9 phases | 1.16e-4 … 0.9799 | 1.18e-4 … 0.9798 |
| φ on the anode surface, max \|φ\| | 2.6e-9 | 2.6e-9 |
| φ on the strip plane − indicator | 0 | 0 |
| Σ 13 strips − whole-plane potential (finite-sum tail) | 7.68e-4 | 8.00e-4 |
| build wall time | 202 s (earlier) / 205 s (overlay) | 207 s |

The sandbox LUT was built from the frozen code with this environment's numpy. `overlay-identity` repeats the no-op
proof against your own frozen LUT on WSL. If your LUT was built with another numpy/scipy, it reports the max
difference instead of bit identity.

**SYNTHETIC (toy transport, not NA60+ physics):**

- `test_restrip.py`:
  - identity restrip byte-identical;
  - stereo 0.15 restrip vs the exporter formula evaluated directly: 3.3e-16 mm, with base V −19 → −18 (an index
    change is exercised);
  - transport columns byte-identical;
  - numpy birth path vs restriped export: D_U/D_V 4.2e-16 mm.
- `test_core_geometry.sh`:
  - C++ core with explicit (0.1, 1.7) is byte-identical to the default;
  - core vs numpy at 0.15: 1.9e-16;
  - transport-side moments (t, N, r_eff, C1…S2, A) bit-identical between geometries; D columns move.
- `test_pipeline_geometry.sh`, one toy muon on a coarse toy grid (nr 8 × nθ 16):
  - G vs S0: 0.48 % baseline vs 0.51 % stereo 0.15;
  - K2A vs S0: 4.6 % vs 4.4 % (coarse toy grid; the real r14×θ72 grid gave 0.15–0.20 %);
  - the LUT/geometry mismatch is refused;
  - the engine refuses moments of another strip frame, and matched moments reproduce s_kernel K2A;
  - identities PASS: transport unchanged, weighting changed;
  - toy geometry effect: u − truth moves +0.728 mm (S0) vs +0.727 mm (K2A), and the U/V charge ratio moves −6.8e-5
    vs −7.7e-5.
- Regression of the changed tools at default settings, on toy data:
  - `t_shaped_compare` per-muon output and summary text are byte-identical to HEAD;
  - `s_grid_spec` output is identical;
  - `test_birth_moments.sh` and `engine/test_engine.sh` PASS.

## 6. Procedure and acceptance gates (proposed; nothing beyond gate 0 has run on real data)

| gate | step (runner mode) | pass criterion (proposed, for your approval) |
|---|---|---|
| 0 | `selftest`, `overlay-identity` | synthetic PASS; baseline overlay LUT bit-identical to your frozen LUT (or max \|Δ\| reported and explained) |
| 1 | `lut` | `geom_lut_check` PASS: interpolation ≤ 5e-5 (baseline 1.5e-5), truncation ≤ 1e-6 (baseline 3.6e-8), bounds, anode ≤ 1e-6; baseline LUT checked by the same tool |
| 2 | `prepare` (5 muons) | restrip OK (index cross-check, transport hash); exact S0 replay PASSes its own frozen checks (maximum principle, landing strip, census) |
| 3 | `green`, `hashes` | grid spec (wires, radii) identical to the baseline; `transport_id` identical; `weighting_id` differs; baseline LUT and exports intact |
| 4 | `closure` | G vs S0 per muon and plane (shaped, 100/200 ns, ion only and + electron step) ≤ 2 × the baseline G closure of the same muon, and ≤ 1e-3 |
| 5 | `kernel`, `report` | K2A vs S0 ≤ 2 × baseline worst (baseline 0.15–0.20 % L2), U/V charge ratio error ≤ 1 %; K2A vs G reported |
| 6 | `births`, `engine` | birth moments vs restriped export: PASS (1e-9); K2A/K2 waveforms PASS; compiled/vector engine == export path |
| 7 | `report` (geometry effect) | reported, no criterion: charge sharing, U/V ratio, waveforms, threshold crossing, multiplicity, w/u vs truth, for S0 and K2A, and the fast model's error on each effect |

A scan (more angles or pitches) is proposed only after gates 1–6 pass for this geometry.

## 7. Expected computational cost (five development muons)

Measured here are marked “measured (sandbox)”. Everything else is an expectation from the baseline records, not a
measurement.

| step | cost | basis |
|---|---|---|
| LUT build | 205–207 s, 21 MiB | measured (sandbox, 2 cores) |
| LUT checks (alt + frozen) | 36 s each | measured (sandbox) |
| overlay-identity (optional) | ≈ 205 s, 21 MiB | measured (sandbox) |
| frame (p3_birth_moments) per muon | ≪ 1 s | baseline births timing |
| restrip per muon | one copy of `T1_export/vertices.f64` (I/O only); disk = size of the baseline T1_export | expected |
| S0 exact replay per muon | ≈ the baseline S0 wall time, 16–275 s per muon | expected: same queries, same LUT size |
| electron-step term | seconds per muon | expected |
| Green grid paths | **0** (reused) | — |
| Green table | ≈ baseline `green_table_TAG_s` (BASE/timing.txt); size ≈ baseline, ∝ number of D nodes | expected; toy grid: 38 → 45 MiB because the D range followed the restriped columns |
| G closure, kernels, K-apply, births, engine | ≈ baseline (minutes in total) | expected |
| new disk total | ≈ 21 MiB LUT + 5 × (T1_export + S0) + one Green table + one kernel set | expected |

Nothing needs PBS. No new avalanche, Geant4 or Heed step is involved.

## 8. Not done / limitations

- No real-muon result exists for the alternative geometry; gates 1–7 are for you to run.
- Strip gaps, a single-plane change and transport-coupled parameters are excluded by construction.
- The physics comparison uses the PROVISIONAL shaper and the ideal readout (assumption 4).
- Baseline closure values for gate 4 are read from your existing baseline reports; the criteria are proposals.
