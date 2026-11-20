# Audit of historical NA6PRoot stereo-strip field solvers vs current 3D thin-wire pilot

Reviewed 2026-11-20, after new near-anode U/V gradient convergence test 1.0→0.5mm (6.03% vector gradient L2) and 0.5→0.25mm (2.5715% L2). This is a **source-code audit** of earlier committed solver implementations, not validation of historical test_runs artifacts.

## Findings: a full wire-surface FEM implementation already exists

Earlier Stage-B/C1 strip-sharing implementation: `test/mwpc-avalanche/MWPCStereoReadout.h` uses Garfield analytic full-wire cell for overall cathode induced charge, but approximates distribution over individual cathode strips with an analytic planar-strip model. Stage-C full cathode gap omits anode-wire conductor boundary in the **strip** weighting solve. Legacy half-gap solves put an artificial wire-plane boundary at v=0 and cannot describe opposite-half-gap carriers. Segmented charges are RESCALED to match full whole-cathode coupling. Historical C1 limitations are clearly documented in `test/mwpc-stage-c/numerical-weighting/README.md` and `test/mwpc-stage-c/README.md`, which say the analytic hybrid strip sharing is unstable for weak stereo differences. This is related to, but **not the identical numerical segmentation error** found in our new thin-wire line-density model.

Earlier Stage-C2 numerical solution is already implemented: `test/mwpc-stage-c/numerical-weighting/generate_weighting_model.py`, `test/mwpc-stage-c/numerical-weighting/weighting_config_c2b.json`, `test/mwpc-stage-c/c2_numerical_weighting_probe.cxx`, and `test/mwpc-stage-c/MWPCNumericalStereoReadout.h`. Uses real cylindrical anode holes cut from Gmsh 3D gas mesh via OpenCASCADE BooleanDifference; Elmer finite-element Dirichlet solve puts anode metal boundaries at φ=0, selected strip at φ=1, all other cathode regions at φ=0. Output results are imported by Garfield `ComponentElmer`. Explicit Gmsh mesh size controls around wire and cathode planes; numerical-domain Neumann side walls need patch-size convergence. `MWPCNumericalStereoReadout.h` uses translated single-central-strip maps and Ramo endpoint differences for readout, rejecting nonfinite or out-of-domain map queries; unlike Stage C1, no sidewise strip normalization rescales charge. This is a stronger physical *wire surface* formulation than our new axial monopole line-element method, which enforces only an angular MEAN potential on wire cylinders. It should become the basis for future production U/V weighting fields, subject to mesh and boundary tests.

**Crucial unverified status:** repository source describes C2 FEM smoke, full-event and edge-scans, but the actual local `test_runs/` map/results are not present in GitHub and have not been checked in this audit. `c2_numerical_weighting_probe.cxx` only enforces LOOSE smoke criteria (selected cathode potential >0.5, near-wire potential <0.2, reflected profile discrepancy <0.02), NOT 0.1% near-wire gradient convergence. C2 full-event reconstruction had documented unusual edge +/-u behavior and noisy strip-template estimators, but those could arise from stochastic avalanches, geometry or reconstruction and are not evidence of the new 6% longitudinal collocation discretization error. Don't claim the old FEM was validated at the needed accuracy until its actual mesh, numerical field gradients and stored outputs are checked.

## Important geometry mismatches before reusing legacy C2 code

1. Older `weighting_config_c2b.json` uses ±2.5mm gaps, 4mm anode pitch, **5** explicit anodes (`half_wires=2`), 30μm diameter, ±12mm across, ±12mm along and 1800V. Current reference uses candidate ±3mm, **3mm** pitch, **9** anodes (`half_wires=4`), 30μm diameter. Anode voltage does NOT influence the unit weighting boundary solve (it matters for physical drift map only). Legacy `validate()` in the generator requires `u_half_span_mm > (half_wires + 0.30)*wire_pitch_mm`: a 9-wire 3mm model needs u halfspan >12.9mm, not legacy 12mm. Both u/w domains need independent enlargement/refinement tests.

2. Stereo slope/channel label mismatch: historical `generate_weighting_model.py:strip_expression` defines *minus* projection `w*cos(α)-u*sin(α)` and *plus* projection `w*cos(α)+u*sin(α)`. Current `StereoStripBoundaryApprox.hh` assigns *PLUS* projection `w*cos(α)-u*sin(α)` and *MINUS* `w*cos(α)+u*sin(α)`. These are opposite naming assignments. Physical slopes are equivalent up to cathode-label swap under symmetric geometry, but when combining or comparing U/V charges the label conventions must be mapped explicitly and checked on nonzero u positions (a u=0 test will miss sign mismatch).

3. Definition of 1.7mm pitch/width differs slightly: historical C2 Dirichlet mask is `abs(w cosα +/- u sinα - stripIndex*P) <= width/2`, so P and width are perpendicular-to-strip distances. The current pilot's P=1.7mm and width=1.7mm are measured **along anode-wire direction at fixed u**. To reproduce 1.7mm along-wire spacing/width in the legacy projection, use projected P=1.7*cosα and projected width=1.7*cosα, where cosα=1/sqrt(1+0.1²) ~0.995037. This 0.5% scale change should not be silently ignored in systematic comparisons.

4. The old FEM `wire_size_mm=0.020`, `wire_refine_radius_mm=0.25`, `cathode_size_mm=0.20` controls are a STARTING MESH, not demonstrated fine enough for near-anode U/V field gradients at 1.1 wire radius (wire radius=0.015mm). Verify actual tetrahedron order, cylinder surface resolution, strip-boundary Dirichlet sharpness and refinement convergence.

5. Historical `run_c2_weighting_smoke.sh`, `run_c2_full_numerical.sh` and `run_c2_u_cell_scan.sh` contain `rm -rf` on their respective fixed run outputs. DO NOT run old scripts unmodified against valuable local test_runs. First inspect environment/tools, existing C2 map results, then build an isolated no-overwrite map generator and validation runner in a dedicated worktree.

## Recommended next research action

**STOP** brute-force refinement of the current thin-wire 720-element collocation field as the main path. Use it as a fast exploratory/independent mathematical control. Reuse and update the old Gmsh/Elmer **explicit conductor** model for the new candidate geometry, run two standalone meshes with tighter near-cylinder surface and strip-edge refinement, compare local U/V weighting vector gradients at identical 1.1R/2R/5R field points, verify wire-surface Dirichlet φ=0 and U+V channel orientation, and compare whole-cathode weighting to Garfield as an independent reference. Only after field convergence and electrode/microgeometry fidelity are established, advance to isolated ion Ramo strip charge, then full physical U/V ion current. No electronics.

Suggested read-only local preflight:

~~~bash
for cmd in gmsh ElmerGrid ElmerSolver; do
  command -v "$cmd" || printf 'missing: %s\n' "$cmd"
done
find ~/na6p -maxdepth 9 -type f -name 'weight_plus_0.result' -print 2>/dev/null | head -n 15
~~~
