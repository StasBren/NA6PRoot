# P3 2100V candidate — direct 3D Elmer U/V weighting-field gradient gate

## UPDATE: critical wire faceting actually measured 2026-11-21 — PRIORITY CHANGED

The original 1.85M-node Elmer P2 map has mesh type 510 and wire-surface triangle type 306, but corner and mid-edge node audits found median radial mid-edge offset **−0.890395 μm**, extreme **−2.221988 μm**, versus a physical near-wire test clearance of only **1.5 μm at r=1.1R**. This is geometric faceting, not merely a 12→6 μm mesh-size convergence problem. The original `Mesh.SecondOrderLinear=1` must not be treated as exact cylinder geometry. No production U/V currents using this map.

**FIRST controlled study (after overnight muon recovery)**: at the SAME 12μm nominal wire mesh resolution, generate a fresh curved-wire P2 FEM field via the updated `run_p3_2100_fem_field_only.sh curved /absolute/new/output` runner. It opts into `Mesh.SecondOrderLinear=0` and `Mesh.HighOrderOptimize=1` without modifying prior baseline data or defaults. This is a numerical comparison, NOT completed or validated WSL output. Audit the **new** `mesh.boundary` wire radii using `audit_elmer_wire_surface.py`, confirm Elmer solves and no curved-element Jacobian/interpolation failure with Garfield `ComponentElmer`. Sample original and curved maps at identical points using `mwpc_stage_c_fem_stereo_near_wire_probe` (tag `baseline` for original, `curved` for new). Compare with `compare_fem_stereo_near_wire.py --study curvature` and matching configs and fresh output CSV. A successful curved-node audit is necessary but NOT proof of field convergence.

**SECOND, later**: if curved P2 compatibility and field sampling pass, compare 12→6μm at fixed curved geometry, with a fresh `curved_refine` map. Do NOT perform the previous straight-sided `wire_refine` as the principal next physics step; it conflates element-size change with geometrical polygonization error. Do NOT run large Elmer solves concurrently with the frozen overnight 120-muon recovery. The baseline 0.11% **whole-cathode analytic** agreement and 2.57% **line-source strip** sensitivity do not establish accuracy of the FEM individual strips.

Updated 2026-11-21. This is the required numerical gate **before** any physical U/V strip current from the weighted real-muon ion archives. No electrons, positive-ion drift, avalanche simulation or electronics are run by the commands in this document.

## Scientific context and existing immutable baseline

- The archived `test_runs/mwpc_stage_c/p3_2100_fem_field_baseline_v1` is already Gmsh/Elmer complete with both selected-strip field results and SHA256 checks OK. Its source geometry is the nine-wire 2100V P3 reference: 3mm pitch, 30μm diameter, cathodes ±3mm, selected plus/minus strip width and pitch **1.7mm measured along the wire**, stereo slopes ±0.1. Mesh has 1,845,601 nodes and 1,502,492 elements.
- The Gmsh/Elmer map uses genuine grounded cylindrical anode conductor surfaces, not the approximately zero-angular-mean line-source approximation whose individual-strip gradients showed 6.03% and 2.57% refinement changes.
- The baseline numerical mesh input has `wire_size_mm=0.012`. Our first controlled refinement modifies ONLY this to 0.006mm. It does not alter drift HV, physical geometry, strip projection, channel sign convention, other mesh settings or artificial side-wall positions.
- For the finite Elmer patch, zero flux on artificial sidewalls remains a potential systematic. Two-mesh comparisons test near-wire discretization, **not** finite-patch size or PCB/material realism.
- The existing microscopic K256 1000-muon overnight campaign is a SEPARATE frozen worktree. Do not switch, merge, build or otherwise change the campaign worktree or archived ion data.

## First: build read-only near-wire map probe on the independent checkout

~~~bash
cd ~/na6p/cathode-independent-stereo-2100
git fetch origin feature/p3-2100v-wire-aware-stereo-weighting-20261120
git merge --ff-only FETCH_HEAD
export GARFIELD_INSTALL="$HOME/na6p/install/garfieldpp"
export LD_LIBRARY_PATH="$HOME/na6p/install/root/lib:$GARFIELD_INSTALL/lib:${LD_LIBRARY_PATH:-}"
cmake -S test/mwpc-stage-c -B test/mwpc-stage-c/build-fem-strip-probe \
  -DCMAKE_PREFIX_PATH="$GARFIELD_INSTALL"
cmake --build test/mwpc-stage-c/build-fem-strip-probe \
  --target mwpc_stage_c_fem_stereo_near_wire_probe -j2
BASE="$PWD/test_runs/mwpc_stage_c/p3_2100_fem_field_baseline_v1"
VALID="$PWD/test_runs/mwpc_stage_c/p3_2100_fem_validation_v1"
mkdir -p "$VALID"
(cd "$BASE" && sha256sum -c SHA256SUMS)
PROBE="$PWD/test/mwpc-stage-c/build-fem-strip-probe/mwpc_stage_c_fem_stereo_near_wire_probe"
"$PROBE" --map-dir "$BASE/map" --tag baseline \
  --out-csv "$VALID/baseline_nearwire.csv" \
  > "$VALID/baseline_probe.log" 2>&1
tail -n 20 "$VALID/baseline_probe.log"
head -n 4 "$VALID/baseline_nearwire.csv"
~~~

The C++ probe uses exactly 480 fixed locations: 3 central wires × 5 radial shells (1.1R, 1.5R, 2R, 5R, 20R) × 8 azimuthal angles × 2 longitudinal w positions × 2 cathodes. Every shell is outside the grounded 15μm wire radius; the closest shell is just 1.5μm above the surface. The code reads the genuine Gmsh/Elmer `plus_strip_0` and `minus_strip_0` results through Garfield++ `ComponentElmer` and reports φ and all **three** vector weighting-field components E_w=−∇φ in 1/mm; derivatives are evaluated by finite differences at h and h/2, strictly outside the wire.

Checks: individual strip target/opposite cathode potential sanity; finite potential/gradients; 3D U/V mirror-symmetry deviation; finite-difference step sensitivity. A reading of 0 for `WeightingPotential` can also occur when an interpolation point is outside the mesh: inspect the sampling logs and value patterns; **do not** take a completed sampling program alone as proof of validity.

## Second: after checking memory availability / overnight muon recovery, solve a separate refined field

An additional quadratic Elmer 3D mesh can consume considerable RAM/CPU; do not compete with active microscopic 120-event overnight recovery on a RAM-limited laptop. Preserve baseline; use a NEW output directory. Confirm enough free memory and disk before launching.

~~~bash
cd ~/na6p/cathode-independent-stereo-2100
free -h
df -h "$HOME"
RUNNER=test/mwpc-stage-c/numerical-weighting/run_p3_2100_fem_field_only.sh
REFINED="$PWD/test_runs/mwpc_stage_c/p3_2100_fem_field_wire006_v1"
nice -n 10 bash "$RUNNER" wire_refine "$REFINED" \
  > "$PWD/p3_2100_fem_field_wire006_v1.log" 2>&1
tail -n 20 "$PWD/p3_2100_fem_field_wire006_v1.log"
(cd "$REFINED" && sha256sum -c SHA256SUMS)
~~~

The runner **refuses to overwrite** preexisting output. A successful Elmer solve and checksum verification do not imply mesh convergence.

## Third: sample refined U/V fields and compare identical locations

~~~bash
cd ~/na6p/cathode-independent-stereo-2100
BASE="$PWD/test_runs/mwpc_stage_c/p3_2100_fem_field_baseline_v1"
REFINED="$PWD/test_runs/mwpc_stage_c/p3_2100_fem_field_wire006_v1"
VALID="$PWD/test_runs/mwpc_stage_c/p3_2100_fem_validation_v1"
PROBE="$PWD/test/mwpc-stage-c/build-fem-strip-probe/mwpc_stage_c_fem_stereo_near_wire_probe"
"$PROBE" --map-dir "$REFINED/map" --tag wire_refine \
  --out-csv "$VALID/refined_nearwire.csv" \
  > "$VALID/refined_probe.log" 2>&1
python3 test/mwpc-stage-c/compare_fem_stereo_near_wire.py \
  --baseline "$VALID/baseline_nearwire.csv" \
  --refined "$VALID/refined_nearwire.csv" \
  --baseline-config "$BASE/input_config.json" \
  --refined-config "$REFINED/input_config.json" \
  --out-summary "$VALID/mesh_gradient_comparison.csv"
cat "$VALID/mesh_gradient_comparison.csv"
~~~

The comparator rejects geometry/stereo mismatches and uncontrolled changes to mesh settings: only wire_size_mm=0.012→0.006 is allowed. It reports relative vector-gradient L2 and worst local pointwise gradient difference per U/V side and radial shell, max absolute potential/gradient changes, and the internal h vs h/2 finite-difference consistency of each mesh. A 2% screen is merely a diagnostic, NOT an accuracy claim. If near-wire meshes fail to stabilize, repeat targeted refinements/inspect Gmsh surface mesh and/or improve FEM boundary representation.

## Only after a satisfactory 3D field gate

1. **One muon:** load one already-audited K256 ion-birth archive; reuse the independent 2100V physical electric field/ion mobility and full RKF transport. On each ion path evaluate each selected strip weighting-potential difference and accumulated positive-ion induced charge in consistent 10ns bins, tracking late arrivals beyond 100μs and enforcing unit-ion Ramo endpoint closure. Do not confuse with reconstructed collected electrons or full electronics pulses.
2. **Independent checks:** verify U/V symmetry, charge signs/partition and that summed sufficiently numerous strips reproduce the wire-aware whole cathode only within demonstrated finite active strip coverage; compare field and charge under mesh refinement. Add prompt electron contribution using saved paths or properly rerun electron microscopic transport if archived data is insufficient.
3. **Fast response for ~1000 muons:** train a new strip-specific ion-response kernel on a controlled set of start locations/trajectories and check held-out muon direct-vs-fast U/V waveforms before applying fast kernels to the completed microscopic campaign. The historical 256 whole-plane templates are NOT individual-strip kernels and cannot be silently reused as such.
4. Every muon index carries explicit state: SHA-audited, failed, incomplete, late-endpoint tracking. Reporting an unfinished campaign as 1000 good U/V events is prohibited. Keep front-end electronics outside this stage.
