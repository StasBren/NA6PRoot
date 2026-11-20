# Stage C: independent-muon 10 ns check + candidate two-sided U/V stereo geometry

**Draft code: not yet compiled or executed on the user's WSL.** Physics outputs are unknown until the test runs. NO electronics response work is included.

## Source design basis and scope

The NA60+/DiCE design draft New_chamber_design(1).pdf, sections 4.2 and 6, considers a MuTr strip width measured ALONG the anode wire, specifically w=1.7mm, with lines crossing the wires at pi/2 +/- alpha, tan(alpha)=0.1. It discusses two-sided simple strips with either 3+3mm or potentially 2.5+2.5mm symmetric gaps. **The intended geometry is not finalized or as-built.**

Numerical design candidate here: width along wire = 1.7mm, strip-pitch-along-wire = 1.7mm (IDEAL contiguous strips, with zero copper-to-copper gap in this pilot), tan(alpha)=0.1, symmetric 3+3mm gap. The true copper gap, dielectric/PCB layers, guard electrodes, traces, finite module extent and optimal gap remain unresolved.

The separate frozen Garfield model uses nine ideal wires, 3mm anode pitch, 30um wire diameter, 2100V, symmetric +/-3mm cathode planes, Ar/CO2 70:30 and a provisional Ar+ in Ar mobility. This is not the built prototype.

Coordinates: u across anode wires, w along the anode wires, v normal to the cathodes, with anode plane v=0. Upper strip n centerline at w=n*p+tan(alpha)*u; lower strip n centerline at w=n*p-tan(alpha)*u, i.e. U coordinate sU=w-0.1u, V coordinate sV=w+0.1u. Strip width/pitch w and p are defined along w at fixed u. Cross-wire coordinate reconstruction u=(sV-sU)/(0.2), which magnifies differential strip-position errors by factor 5. No reconstruction uncertainty model is included yet.

## Wire-free analytic weighting field: geometry diagnostic ONLY

StereoStripBoundaryApprox.hh is a correct analytic 3D harmonic solution for a single infinite tilted strip on one conducting plane, with an opposing grounded conducting plane, **but with ALL anode wires omitted**. This is a boundary and rotation validation model, not yet a real MWPC strip weighting field. It deliberately does not compute physical induced strip currents: including anodes (grounded in the weighting solve), finite PCB strips/insulation and electrodes needs a wire-aware 3D BEM/FEM/neBEM field solution.

The geometry diagnostic checks upper/lower strip orientations, conductor boundary potential, reflection under u -> -u and v -> -v, 3D Laplace harmonicity, potential gradients and ideal-strip sum rules. For contiguous strips, summing all upper strip potentials reproduces the **wire-free** whole-plane linear potential, not the original Garfield nine-wire whole-plane response. Outputs include a ROOT histogram, a PNG four-panel potential map and CSV sample grid. Absolutely no VMM, electronics or physical strip signal is computed.

## Independent-muon fixed-template validation

Choose three frozen, already-audited K256 microscopic muons, none equal to the original muon0 256-template training source:

| Event index | Source | Station | Collected e / Heed primary | Saved K256 ion rows |
|---|---|---:|---:|---:|
| 0301 | phi40-bore | 0 | 44,351.85 | 4,096 |
| 0408 | jpsi40-bore | 3 | 29,818.18 | 7,424 |
| 0446 | omega40-bore | 2 | 36,315.79 | 5,376 |

The new runner verifies completed event seal, per-event SHA256, ROOT physical audit and uncapped Garfield evidence, then computes at 10ns over 100us: (a) weighted direct ion RKF for each independently sampled ion reservoir, and (b) the **same exact ion reservoir** with the **frozen original 256 ion-start sites** and four-nearest-template model, mirrored only at response lookup. It does not resample ions or generate new avalanches. The C++ audit saves 50/100/200/400/1000/100000ns charge comparisons, raw 10ns L2, rebinned 100ns L2 and first-1000ns ROOT/PNG overlays. The three muons are drawn from different physical channels and stations, rather than the training muon.

## WSL run — make a NEW isolated checkout; do not alter active overnight checkout

~~~bash
cd ~/na6p/cathode-overnight-2100
git fetch origin feature/p3-2100v-independent-muons-stereo-boundary-20261120
git worktree add --detach ~/na6p/cathode-independent-stereo-2100 FETCH_HEAD
cd ~/na6p/cathode-independent-stereo-2100
export GARFIELD_INSTALL="$HOME/na6p/install/garfieldpp"
export LD_LIBRARY_PATH="$GARFIELD_INSTALL/lib:$HOME/na6p/install/root/lib:${LD_LIBRARY_PATH:-}"
cmake -S test/mwpc-avalanche -B test/mwpc-avalanche/build-independent \
  -DCMAKE_PREFIX_PATH="$GARFIELD_INSTALL" \
  -DROOT_DIR="$HOME/na6p/install/root/cmake"
cmake --build test/mwpc-avalanche/build-independent --target \
  mwpc_plane_response_from_ions_root \
  mwpc_plane_muon_kernel_compare \
  mwpc_audit_independent_muon_kernel_10ns \
  mwpc_stereo_strip_boundary_pilot -j2
~~~

Read-only microscopic-input verification and three new independent runs:

~~~bash
RUNNER=test/mwpc-stage-c/cathode-ion-reference/run_independent_3muon_dt10.sh
bash "$RUNNER" preflight
bash "$RUNNER" run > "$PWD/independent3_run.log" 2>&1
tail -n 25 "$PWD/independent3_run.log"
OUT="$PWD/test_runs/mwpc_stage_c/independent_muon_dt10_v1"
cat "$OUT/independent3_events.csv"
cat "$OUT/independent3_windows.csv"
ls -lh "$OUT/independent3_early.png"
~~~

Read-only, wire-free strip geometry tests and visualization:

~~~bash
GEOM="$PWD/test/mwpc-avalanche/build-independent/mwpc_stereo_strip_boundary_pilot"
"$GEOM" --gap-mm 3 --strip-width-mm 1.7 \
  --strip-pitch-mm 1.7 --tan-stereo 0.1
mkdir -p "$PWD/test_runs/mwpc_stage_c/strip_boundary_v1"
"$GEOM" --gap-mm 3 --strip-width-mm 1.7 \
  --strip-pitch-mm 1.7 --tan-stereo 0.1 \
  --out-prefix "$PWD/test_runs/mwpc_stage_c/strip_boundary_v1/stereo_wirefree"
~~~

## Interpretation and next geometry milestone

Original benchmark K512 muon0, 256 template positions: upper first-100ns charge error from -20.24% to +0.206%, lower from -21.92% to +0.177%, when using 10ns instead of 100ns template shifting. The native 10ns waveform L2 remained 2.94%/3.45%. **These are single-muon findings.** Independent test outcomes must not be inferred in advance.

After confirming or revising spatial interpolation generalization on the three independent muons, a full stereo-strip signal requires 3D wire-aware weighting field: all non-readout electrodes and anodes grounded, individual strip at unit potential, realistic PCB geometry, boundary conditions, mesh convergence and strip/whole-plane sum rules. Only after this validation should we integrate true ion trajectories to produce physical U/V currents. No electronics in this project phase.
