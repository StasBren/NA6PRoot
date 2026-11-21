# NA60+ candidate P3 2100V: next numerical/physics gates for true U/V stereo-strip response

Status at 2026-11-21: original 1000 frozen microscopic muon inputs are being recovered after 120 ROOT libGraf3d.so loader failures; 877 of the original 1000 slots are SHA/physics-verified, 1 earlier failed and 2 unsealed, 120 library-failure attempts archived for safe retry. Recovery runs in the frozen overnight checkout and MUST NOT be modified by ongoing development.

Nine-wire symmetric 3+3mm candidate FEM baseline built in a separate analysis checkout has completed: Gmsh reported 1,845,601 nodes and 1,502,492 elements; Elmer drift and U/V central-strip weighting solutions are present, and eight input/map files passed SHA256. This is completed numerical solve, NOT convergence certification.

## Gate A — physical FEM strip fields

Reuse completed baseline at `test_runs/mwpc_stage_c/p3_2100_fem_field_baseline_v1/map`; don't rebuild or overwrite. The repo already has `test/mwpc-stage-c/fem_stereo_strip_near_wire_probe.cxx` and `test/mwpc-stage-c/compare_fem_stereo_near_wire.py` for read-only 480 near-wire sample points, U/V and all three finite-difference components, symmetry and finite difference h-vs-h/2 control. Evaluate baseline near-anode shells at r/R=1.1, 1.5, 2, 5, 20 across three central wires, eight angles, w=0/1.1mm. Check valid Garfield ComponentElmer interpolation, including actual Elmer tetrahedron element order; Gmsh requests quadratic -order 2, therefore an unsupported assumption of P1 elements should not drive design decisions without checking mesh.element types. h vs h/2 agreement INSIDE an FEM element by itself does not establish spatial mesh convergence.

Once overnight muon recovery is complete and system memory permits, run an independent full Elmer field-only `wire_refine` solve with a separate fresh output path, varying only `wire_size_mm` 0.012→0.006mm. Compare identical mesh-independent sample positions using `compare_fem_stereo_near_wire.py`; report max absolute potential shift, absolute field difference mm^-1, vector-gradient L2 with weak-channel amplitude floor and max local relative difference, separately by shell and U/V. Check observed boundary values, strip handedness and u→−u, v→−v, U↔V reflection. The current 2% rule is a SCREENING criterion, not an experimentally established error budget. A second mesh comparison alone is insufficient to certify physical U/V signal accuracy; near-strip-edge mesh effects and domain patch-size convergence remain.

## Gate B — one ion and one complete muon

After above, import field maps via existing `Garfield::ComponentElmer`, retaining physical ion transport from `ComponentAnalyticField`/`DriftLineRKF` as previously benchmarked. Start with a single representative positive ion and validate pointwise induced weighting potential, time-binned current and the Ramo endpoint identity `Q_k = q[phi_w,k(r_end) − phi_w,k(r_start)]` (using project's existing sign convention). A near-perfect endpoint identity is a consistency check of current integration, NOT independent validation of the FEM potential. Use the same reference trajectories across numerical meshes to quantify strip Q sensitivity.

Next, load one existing already-verified K256 ion-birth reservoir muon WITHOUT regenerating its microscopic avalanche, apply stored weights, shift by true birth times, and form U/V waveforms at 10ns along with the prompt electron contribution from microscopic electron path archive if that data is available and accessible. Ion current must preserve physically complete drift trajectories and explicitly separate in-window 100us integral from later 100-200us endpoint Ramo charge; no late ions may be silently dropped. Initially central U0/V0 maps are usable ONLY for central channels. To construct the physical multi-strip cluster, independently generate adjacent k strip maps or validate translation of the central map against explicit k=+/-1 maps under finite Neumann sidewall boundaries; never assume exact periodicity for the finite patch. Retain weak strips when computing global metrics but use a nonzero floor for relative per-strip errors. Inspect U/V slope convention with u!=0 and w!=0.

Charge accounting: the sum of weighting potentials for an exhaustive set of independently biased conductor electrodes (including anode wires and gap/guard conductors if modeled) equals one in the ideal electrostatic boundary problem, and sum of their weighting fields zero. Do NOT assume `sum(U strips)+sum(V strips)=1`, or that a sum of all induced endpoint charges equals the transported ion charge. For each full cathode, if its modeled surface is completely partitioned into readout strips, `sum(k in cathode) phi_w,k = phi_w,whole cathode` and corresponding per-bin current sum (same geometry, complete coverage); compare with independent Garfield analytic whole-plane only allowing for FEM finite-bank sidewall/model differences. Current physics config is width=pitch with zero physical gap, an idealization, not actual PCB.

## Gate C — fast U/V response library and ~1000 muons

Tabulate **multi-channel** ion kernels from sampled positive-ion starts parameterized by transverse position around nearest anode (e.g. u-u_wire, v or rho,theta) AND longitudinal stereo-strip phase (w modulo proposed 1.7mm along-wire pitch). Keep a wire-index/edge validity class or prove spatial stationarity; finite bank and finite FEM patch break naive translational invariance. PLUS and MINUS strip phases are `w−tan(alpha)u` and `w+tan(alpha)u`, respectively; they are correlated, not arbitrary independent variables. Exact infinite-chamber strip k translation is along w by k*pitch, but in the finite FEM box this is an APPROXIMATION requiring an explicit k=±1 map control.

Start with precomputed single-ion response kernels including several neighboring strips on BOTH cathodes (number determined by measured charge tail, not an assumed 5–7). Interpolate using training templates frozen BEFORE testing on independent original microscopic muons. Reuse identical K256 sampled ion reservoirs and birth-time shifts for direct RKF vs fast model; assess combined 10ns waveform L2, first100ns charge, first100us charge, late endpoint tails, per-strip charge-sharing fractions, strip centroid u/w and coverage. Existing ~4.81% whole-cathode combined waveform L2 on three held-out muons is a baseline reference, NOT guaranteed strip-level performance or a required sub-percent numerical limit.

Once three or more genuinely held-out muons demonstrate satisfactory full channel-by-channel and centroid accuracy, apply fast library to currently recovered ~1000 microscopic muons. Report VERIFIED, FAILED, UNSEALED and ABSENT categories and never call the campaign 1000 fully simulated muons until verified. Map storage/use and throughput must be monitored; do not rerun microscopic avalanches. Electronics/noise/VMM absent from this stage.

## Safe baseline read-only FEM field probe (run when concurrent muon jobs permit memory)

```bash
cd ~/na6p/cathode-independent-stereo-2100
export LD_LIBRARY_PATH="$HOME/na6p/install/root/lib:$HOME/na6p/install/garfieldpp/lib:${LD_LIBRARY_PATH:-}"
cmake -S test/mwpc-stage-c -B test/mwpc-stage-c/build-fem-probe -DGarfield_DIR="$HOME/na6p/install/garfieldpp/lib/cmake/Garfield"
cmake --build test/mwpc-stage-c/build-fem-probe --target mwpc_stage_c_fem_stereo_near_wire_probe -j2
MAP="$PWD/test_runs/mwpc_stage_c/p3_2100_fem_field_baseline_v1/map"
TEST="$PWD/test_runs/mwpc_stage_c/p3_2100_fem_validation_v1"
mkdir -p "$TEST"
nice -n 10 test/mwpc-stage-c/build-fem-probe/mwpc_stage_c_fem_stereo_near_wire_probe --map-dir "$MAP" --out-csv "$TEST/baseline_near_wire.csv" --tag baseline > "$TEST/baseline_probe.log" 2>&1
tail -n 20 "$TEST/baseline_probe.log"
```

No command in this plan deletes any saved muon or FEM data. If a CSV already exists the probe explicitly refuses to overwrite it; choose a new output prefix.
