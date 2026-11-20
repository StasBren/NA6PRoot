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
- With `--convergence`, refine 1.0mm line segments to 0.5mm and extend the charge-source correction region by 4mm each side. Report maximum changes at fixed interior points for both a selected strip and whole plate. Need also study image truncation if significant.

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
