# Stage C: 3D grounded-anode-wire weighting potential for two-sided stereo strips

**Development stage:** a new numerical, physically wire-aware *thin-wire approximation* is committed, but must be compiled and checked on user WSL. It does **not yet** simulate microscopic avalanches, drifting ions or electrons, strip currents, electronics or any VMM response.

## Design source and explicit provisional assumptions

Source: NA60+/DiCE `New_chamber_design(1).pdf`, especially §§4.2 and 5–6. The design note defines MuTr strip width **along the anode-wire direction**, considers 1.5–2.0 mm (nominal 1.7 mm), and defines strip directions perpendicular to the wire with ± stereo angle satisfying tan(α)=0.1. A simplified **two-sided** strip readout is preferred over the complicated single-cathode interlaced pattern. Symmetric cathode gap **3+3 mm** is still being considered, with **2.5+2.5 mm** raised as a potential alternative. Do not substitute Prototype-3's as-built 4mm pitch, 2+4mm asymmetric gaps, 1.8kV without new explicit geometry assumptions.

The present *numerical reference* separately uses the frozen Garfield P3 noFW 9-anode ideal field: 3 mm wire pitch, 30 μm wire diameter, ±3 mm cathodes, 2100 V. These are not necessarily the final production chamber parameters. For the readout boundary calculation the absolute operating voltage does **not** appear in the weighting-field solve; anode metal is grounded for weighting potential. The 1.7mm **strip pitch** is merely assumed equal to 1.7mm width, i.e. contiguous zero-gap copper, until an interstrip copper gap is specified. No dielectric or real PCB stack is modeled.

Coordinates are u across wires, v perpendicular to plates (anodes v=0), w along wires. Upper U strip centerline: w − 0.1u = n p; lower V: w + 0.1u = n p. The strip potential is three-dimensional, and the opposing inclinations can be independently checked under reflection (u,v)→(−u,−v).

## Mathematical method

Boundary condition for selected strip k: φ_k=1 on that strip, φ=0 on all other conducting electrode areas, opposite cathode, and **all anode wires**. We decompose

   φ_k(u,v,w) = φ_k^(two_plate_strip)(u,v,w) + Σ_(wire j, segment m) λ_jm G_D[(u,v,w),(u_j,0,w′_m)]

where φ^(two_plate_strip) is the exact ideal infinite two-parallel-plate strip harmonic solution from the earlier validated `StereoStripBoundaryApprox.hh`. The correction uses line-charge segments along each anode and the **two-plane Dirichlet Green function G_D**, built by charge-image pairs whose exact infinite series vanishes on the upper and lower plates.

Discretize each 30 μm diameter wire into constant-line-density intervals along w; choose unknown λ on each wire/segment and collocate at four angular quadrants of the real 0.015mm-radius wire cylinder. Solve a dense shared LU system that drives the **mean** potential at each sampled wire section to zero. This models real anode shielding in 3D, but is **not exact cylinder surface BEM**. Check the *unconstrained* 12-angle wire-surface residual rather than declaring boundary perfection.

Default numeric settings: 9 wires, ±10mm w correction region, 1.0mm wire-charge segments, 12 Green image pairs to each side, 180 line-density unknowns, 1.7mm strip width & pitch, gap ±3mm. The line charge source support is finite along w even though actual wires are physically long; a convergence study must extend w support and halve the charge-segment size. Green image truncation is another approximation: evaluate against changed image count if needed.

Importantly this field solution is not equivalent to running `Garfield::ComponentAnalyticField` with a tilted strip. Instead use Garfield's existing *whole-cathode* weighting potential for the exact same 9-wire ideal model as an independent diagnostic reference: sum the many adjacent stereo strip solutions and compare to the screened whole plate potential, then compare that to the Garfield analytic full-plane weighting values at several (u,v,w) points.

## Numerical tests and what PASS means

- Selected cathode strip at unit potential; unselected parts and opposite plane zero.
- Actual 30μm anode cylinders nominally ground; separately report max residual over 12 azimuthal surface samples and several central w positions (not only averaged collocation).
- U versus V opposite stereo reflection under (u,v)↦(−u,−v) with unchanged w.
- Sum of 23 contiguous upper strips should reproduce the same numerical solver's **wire-screened whole-upper-cathode** potential in the central detector patch.
- Diagnostic comparison of numerical whole-cathode weighting potential to `Garfield::ComponentAnalyticField::WeightingPotential` for the same 9-wire reference. This is the most important physics fidelity comparison: report the worst ABSOLUTE potential discrepancy; it is **not** automatically accepted just because a boundary test passes.
- With `--convergence`, refine 1.0mm line segments to 0.5mm, extend charge-source correction by 4mm each side, and independently increase Green image count by 8 on each side. Report maximum differences at fixed interior points for both a selected strip and whole plate. Check the numerical 3D Laplacian away from electrodes as a separate harmonicity test.

**A PASS currently means only numerical boundary, symmetry, partition and broad refinement checks passed within the coded pilot thresholds. It does NOT certify physical strip current, precise charge sharing, or complete field convergence.** Required later: finite PCB/inactive edges, actual copper gaps, traces, dielectric layers, realistic wire segmentation, systematic tolerance convergence, and credible 3D field reference with ions and direct weighting-potential charge closure.

## Safe WSL instructions — NEW isolated worktree

Do not switch, clean, pull, compile or modify the running microscopic 1000-muon checkout. Create a new worktree solely for this field calculation:

~~~bash
cd ~/na6p/cathode-overnight-2100
git fetch origin feature/p3-2100v-wire-aware-stereo-weighting-20261120
git worktree add --detach ~/na6p/cathode-wire-field-2100 FETCH_HEAD
cd ~/na6p/cathode-wire-field-2100
export GARFIELD_INSTALL="$HOME/na6p/install/garfieldpp"
export LD_LIBRARY_PATH="$GARFIELD_INSTALL/lib:$HOME/na6p/install/root/lib:${LD_LIBRARY_PATH:-}"
cmake -S test/mwpc-avalanche -B test/mwpc-avalanche/build-wire-field \
  -DCMAKE_PREFIX_PATH="$GARFIELD_INSTALL" \
  -DROOT_DIR="$HOME/na6p/install/root/cmake"
cmake --build test/mwpc-avalanche/build-wire-field \
  --target mwpc_wire_aware_stereo_strip_pilot -j2
~~~

Then validate and produce sampled 3D field maps and Garfield whole-plane controls:

~~~bash
BIN="$PWD/test/mwpc-avalanche/build-wire-field/mwpc_wire_aware_stereo_strip_pilot"
OUT="$PWD/test_runs/mwpc_stage_c/wire_aware_stereo_v1"
mkdir -p "$OUT"
"$BIN" --gap-mm 3 --strip-width-mm 1.7 --strip-pitch-mm 1.7 \
  --tan-stereo 0.1 --wire-pitch-mm 3 --wire-radius-mm 0.015 \
  --wires 9 --halfspan-w-mm 10 --segment-step-mm 1 \
  --image-pairs 12 --convergence \
  --out-prefix "$OUT/wire_aware_3plus3_v1" \
  > "$OUT/run.log" 2>&1
tail -n 45 "$OUT/run.log"
cat "$OUT/wire_aware_3plus3_v1_checks.csv"
ls -lh "$OUT/wire_aware_3plus3_v1_potential.png"
~~~

The output PNG contains the upper U0 and U1 strip weighting potentials, the lower V0 strip, and the differential **shielding correction** compared with the old wire-free model. All are potentials, not current waveforms. The CSV reports boundary residuals, partition closure, 3D reflection symmetry, Garfield whole-plane deviations and segmentation/finite-w extent convergence.

### Independent control for candidate 2.5+2.5mm

After the first case is scientifically acceptable, run a separate prefix with `--gap-mm 2.5`; do not reuse a previously created prefix or pretend the 2100V microscopic reference holds in a redesigned chamber. The change is a geometry sensitivity study only.

## First actual WSL result — preliminary numerical validation

Compiled successfully on WSL (2026-11-20). With the default candidate geometry, 180 unknown longitudinal wire-line densities and `--convergence`, the completed log reports:

~~~text
wire_surface_residual U=0.003788416 V=0.003789194 whole=0.002538913
stereo_reflection_error=0.00006533331
sum_strips_minus_screened_whole=0.000005026515
max_3D_Laplace_residual_per_mm2=0.00003258118
max_abs_whole_screened_potential_minus_Garfield_analytic=0.0008872376
segment_refinement_change=0.0001592416
wire_source_extent_change=0.0005785199
image_count_change=0.0006713097
~~~

These are absolute *weighting potential* differences and are NOT guaranteed relative accuracy on the individual strip response or local weighting electric fields. The Garfield comparison sampled a small set of interior points. A ~0.0038 wire-surface deviation could be serious very near the anode wire. The program PASS reflects loose diagnostic thresholds, not physical accuracy certification.

### Required next test: near-anode field gradients

New optional `--near-wire` run compares both screened **whole-cathode** potential and the vector weighting field (`E_u`, `E_v` from finite differences) to Garfield's independent analytic 9-wire whole-cathode reference. Samples three central wires, radial distances 1.1, 1.5, 2, 5, 20 times the physical 15μm wire radius, four polar angles and two longitudinal positions (240 configurations including upper/lower electrode references). An independent *physical strip* reference does not yet exist, so this test cannot establish actual strip-current accuracy on its own. Selected strip potential is sampled as a useful extra near-wire check.

Run a NEW output prefix to avoid overwriting the already successful convergence run:

~~~bash
cd ~/na6p/cathode-independent-stereo-2100
git fetch origin feature/p3-2100v-wire-aware-stereo-weighting-20261120
git merge --ff-only FETCH_HEAD
cmake --build test/mwpc-avalanche/build-independent \
  --target mwpc_wire_aware_stereo_strip_pilot -j2
BIN="$PWD/test/mwpc-avalanche/build-independent/mwpc_wire_aware_stereo_strip_pilot"
OUT="$PWD/test_runs/mwpc_stage_c/wire_aware_stereo_v1"
"$BIN" --gap-mm 3 --strip-width-mm 1.7 --strip-pitch-mm 1.7 \
  --tan-stereo 0.1 --wire-pitch-mm 3 --wire-radius-mm 0.015 \
  --wires 9 --halfspan-w-mm 10 --segment-step-mm 1 \
  --image-pairs 12 --near-wire \
  --out-prefix "$OUT/wire_aware_nearwire_v1" \
  > "$OUT/nearwire.log" 2>&1
tail -n 25 "$OUT/nearwire.log"
head -n 10 "$OUT/wire_aware_nearwire_v1_near_wire.csv"
~~~

Before proceeding to physical induced strip currents, assess potential and *gradient* differences by radius, then decide whether the approximate thin-wire collocation suffices or requires improved angular surface treatment and finer near-wire resolution.

## Next numerical gate: near-wire *individual stereo-strip* gradient sensitivity (not Garfield whole cathode)

The newly provided 240-point reference scan CSV shows excellent **screened full-plane** agreement with Garfield: ~0.1109% combined vector-gradient L2, maximum local vector-gradient relative error ~0.548%, max absolute potential difference 3.849e-4. These are not independently validated U or V STRIP field gradients; the near-anode selected strip weighting-potential magnitude can itself be ~0.0038–0.01, similar to the prior 0.0038 wire-surface residual. It would be premature to use the existing thin-wire strip gradient to predict strip induced currents.

New `mwpc_stereo_strip_gradient_convergence` compares the actual *three-component* weighting electric field E_w=-∇φ_w of **separate U0 and V0 stereo strips** at radial shells 1.1R, 2R and 5R, across three anode wires, four polar angles and two longitudinal positions. It recomputes the U and V strip potentials from scratch under four independent perturbations of the field discretization:

- grounded cylinder angular collocation **4→8 equally spaced surface angles** (previously 4-angle mean only);
- longitudinal wire-segment length **1.0→0.5 mm**;
- finite correction-charge source half-span **10→14 mm**;
- Dirichlet image cutoff **12→20 pairs**.

For each perturbation and U/V side, report near-wire **vector gradient relative L2**, maximum local relative gradient error (can be unstable when the reference gradient is tiny), maximum absolute gradient difference in mm⁻¹, and absolute weighting-potential sensitivity. These compare **two numerical approximations** of the same physical model, not a strip-level experimental or independent BEM reference. An apparent gradient convergence is a necessary but not sufficient condition for real U/V signal validity.

Separate worktree already exists at ~/na6p/cathode-independent-stereo-2100. User can update it without touching sealed muon archives or prior successful output prefixes:

~~~bash
cd ~/na6p/cathode-independent-stereo-2100
git fetch origin feature/p3-2100v-wire-aware-stereo-weighting-20261120
git merge --ff-only FETCH_HEAD
cmake -S test/mwpc-avalanche -B test/mwpc-avalanche/build-independent \
  -DCMAKE_PREFIX_PATH="$HOME/na6p/install/garfieldpp" \
  -DROOT_DIR="$HOME/na6p/install/root/cmake"
cmake --build test/mwpc-avalanche/build-independent \
  --target mwpc_stereo_strip_gradient_convergence -j2
BIN="$PWD/test/mwpc-avalanche/build-independent/mwpc_stereo_strip_gradient_convergence"
OUT="$PWD/test_runs/mwpc_stage_c/wire_aware_stereo_v1"
mkdir -p "$OUT"
"$BIN" --gap-mm 3 --strip-width-mm 1.7 --strip-pitch-mm 1.7 \
  --tan-stereo 0.1 --wire-pitch-mm 3 --wire-radius-mm 0.015 \
  --wires 9 --halfspan-w-mm 10 --segment-step-mm 1 \
  --image-pairs 12 --out-prefix "$OUT/strip_nearwire_gradconv_v1" \
  > "$OUT/strip_nearwire_gradconv.log" 2>&1
tail -n 35 "$OUT/strip_nearwire_gradconv.log"
cat "$OUT/strip_nearwire_gradconv_v1_summary.csv"
~~~

All 4 alternative models are solved anew, and failure must be investigated rather than ignored. Do not conflate the output with Garfield's separate exact whole-plane reference, or with physical U/V strip currents. If gradient convergence is poor, improve wire-surface boundary discretization; if stable, next controlled experiment is single-ion induced-charge closure with sampled trajectories and a *physically justified* independent charge-sharing accuracy assessment.

## Measured individual U/V gradients: 1 mm wire segmentation NOT converged

User ran complete independent near-wire U/V strip gradient sensitivity on 2026-11-20. U and V strip field-vector relative L2 difference was **6.0335%** after halving wire-source element length 1→0.5mm; max local relative gradient change **11.74%**. Near-wire maximum absolute potential change **0.004109**. By contrast changing wire-source extent ±10→±14mm changes vector gradient L2 ~0.000202%, image cutoff 12→20 changes ~0.0362%, and angular quadrature mean 4→8 angles has near-zero effect (~6.63e-11 relative gradient L2). The 4-to-8 result does NOT mean anode cylinder surface potential is properly enforced because mean-zero constraints with only axial line monopoles cannot impose pointwise zero around a tilted-strip wire.

The historical program printed `PASS` to indicate calculations completed. **It was not an actual field-convergence pass.** New code now reports explicit `NOT_CONVERGED` when the gradient L2 changes by more than 2% (a screening diagnostic, not calibrated detector tolerance) and accepts `--only-variation segment_half` for a focused high-cost next test.

Next use the already-existing analysis worktree and update only it, not the running 1000-muon microscopic campaign. Dense matrix at 0.25mm segment spacing has 720 unknown wire elements, so this is a relatively costly numerical computation. Run with reduced CPU priority and only one perturbation:

~~~bash
cd ~/na6p/cathode-independent-stereo-2100
git fetch origin feature/p3-2100v-wire-aware-stereo-weighting-20261120
git merge --ff-only FETCH_HEAD
cmake --build test/mwpc-avalanche/build-independent \
  --target mwpc_stereo_strip_gradient_convergence -j2
BIN="$PWD/test/mwpc-avalanche/build-independent/mwpc_stereo_strip_gradient_convergence"
OUT="$PWD/test_runs/mwpc_stage_c/wire_aware_stereo_v1"
nice -n 10 "$BIN" \
  --gap-mm 3 --strip-width-mm 1.7 --strip-pitch-mm 1.7 \
  --tan-stereo 0.1 --wire-pitch-mm 3 --wire-radius-mm 0.015 \
  --wires 9 --halfspan-w-mm 10 \
  --segment-step-mm 0.5 --image-pairs 12 \
  --only-variation segment_half \
  --out-prefix "$OUT/strip_nearwire_seg05_vs025_v1" \
  > "$OUT/strip_nearwire_seg05_vs025.log" 2>&1
tail -n 25 "$OUT/strip_nearwire_seg05_vs025.log"
cat "$OUT/strip_nearwire_seg05_vs025_v1_summary.csv"
~~~

**Stop and reconsider the discretization** if the 0.5→0.25mm gradient differences remain large. Even if this single test looks stable, pointwise grounded wire-cylinder boundary errors and actual strip PCB/microgeometry remain unvalidated. Because the present unknowns are line-source monopoles with only azimuthal *mean* wire-boundary constraints, angular-surface multipoles or a proper 3D surface BEM/FEM may be necessary; increasing the 4/8 angular sample count alone cannot fix angular boundary errors. Do NOT infer 0.11% individual strip accuracy from the whole-plane Garfield reference, and do NOT proceed to strip current/RKF calculations until these issues are resolved.
