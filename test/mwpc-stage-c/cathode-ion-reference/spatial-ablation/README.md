# S-ladder: minimum spatial description of the ion cloud for the U/V strip signals (2100 V)

Timing baseline: **T1** (one creation time per Heed primary). TW/TC remain faster candidates, not baselines.
Reference **S0** = the frozen replay of the T1 export: archived ion positions, their own RKF trajectories, frozen
spectral weighting LUT. Shaping 100/200 ns is PROVISIONAL (CR-RC^3 stand-in).

## Why a single-ion Green's table, and why it is physically consistent
The field is static and 2D (no ion motion along the wire). An ion's trajectory therefore depends only on its start
point (wire k, radius r, azimuth theta) and is independent of its creation time and of its longitudinal position w.
A shift in w changes only the strip-relative coordinate D passed to the weighting LUT (exact translation symmetry of
the 2D field + infinite strips in the frozen model). So every ion's ion-term signal is exactly
`q G_s(k, r, theta, D; t - t_birth)`, with G computed once per start point from that point's OWN RKF trajectory.

A different spatial description = a different set of start points, each carrying its own trajectory. We never move
a charge along a path computed for another start point. Interpolation between grid start points splits the charge
over neighbouring start points (positive weights); harmonic truncation uses signed weights on the angular grid
(a moment-matched surrogate distribution, total charge exact).

`p3_green_grid_paths.cc` compiles the frozen exporter source unchanged (its `main` renamed) and drifts single ions
from a polar grid with the same gas, 9-wire 2100 V field, cathodes, sensor area and RKF step. `s_green_table.py`
evaluates the frozen LUT on those paths: `G[s,k,r,theta,D,tau]` (ion term per unit charge) and `PHI0` (integral of the
avalanche-electron term, for the approximate electron step).

## Ladder (all T1 timing; `s_synthesize.py`)
| Code | Spatial description | Question |
|---|---|---|
| S0 | archived ions, own RKF paths (frozen replay) | reference |
| G | each ion at its own (r, theta, D) from the Green table | **closure gate**: table resolution vs S0 |
| S1 | each ion azimuthally uniform at its own r and D | is angular structure needed at all? |
| S2 | angular harmonics n <= 1 per ion (= cloud first moments resolved in r, D) | first moments sufficient? |
| S3 | harmonics n <= 2, r and D per ion | second harmonic needed? |
| SR | one radius per cloud, r_eff = exp<log r>, full angle and D | compact radial (alone) |
| S2R, S3R | r_eff + n <= 1 / n <= 2 | compact radial + harmonics |
| SW | one D per cloud (charge-weighted mean) | longitudinal: one bin per avalanche |
| SWW | one D per wire (all primaries) | longitudinal: one bin per wire |
| K1, K2 | r_eff + n <= 1 (2) + one D per cloud = (N, <cos>, <sin>, [<cos2>, <sin2>], r_eff, D) per cloud | candidate kernel input |

Cloud = (Heed primary, wire nearest ion birth). Every variant preserves the total ion weight exactly.
The comparison reuses `timing-ablation/t_shaped_compare.py` with S0 as reference: normalized L2 of signed shaped
waveforms, leading-strip peak / time, threshold times and hit flips at thresholds fixed from S0, multiplicities, U/V
charge and amplitude ratios, stereo w, u.

## Closure assessment (practical, not a fixed 0.1% gate)
G vs S0 measures only the numerical error of the Green-table description (grid in r, theta, D, tau; same ions, same
T1 times, same frozen LUT). Spatial-model error is every S-variant minus G, not minus S0. Procedure:
1. **Hard consistency checks** (enforced by `s_synthesize.py` for G, the run stops otherwise): plane window charge
   within 1% of S0, no strip with |Q| >= 1% of the largest whose sign differs, same leading strip in both planes;
   no ion outside the radial grid, all wires covered, D range inside the table; channel labels identical to S0.
2. **Run G alone** (`closure`) and report, per plane: unshaped current and cumulative-charge L2, shaped L2 at 100 and
   200 ns, leading-peak change, threshold times, hit flips, U/V ratios and w, u under every centroid definition.
3. **Decide by the question being asked.** Call eps_G(X) the closure error of observable X and Delta_k(X) the
   difference variant k - G. A ranking or rejection of variant k is unambiguous when |Delta_k(X)| >= 3 eps_G(X)
   for the observables that decide it (and the two variants being compared differ by more than 3 eps_G). Examples:
   eps_G = 0.3% easily rejects a 10% distortion; it cannot separate 0.4% from 0.7%.
   For acceptance, eps_G must also be below the observable tolerance that matters (e.g. a fraction of the
   expected u, w resolution, sub-ns threshold timing, no hit flips attributable to the grid).
4. **Refine only if a decision is ambiguous**: one targeted refinement (MWPC_S_NR, MWPC_S_NTHETA or MWPC_S_DSUB),
   under a new TAG, compared with `closure-compare RUN TAG1,TAG2`. The two closures give a convergence estimate for
   that grid direction; stop when the decision no longer changes.
Results from different grids never overwrite each other (directories and variant names carry the TAG).

## Files
* `p3_green_grid_paths.cc`, `CMakeLists.txt` — grid paths with the frozen transport (build against the wire branch).
* `s_grid_spec.py` — wires, radius grid and D range from the muon's actual ion births. Innermost node =
  r_min - min(0.005 r_min, (r_min - 15 um)/2): never above the innermost ion, never on or inside the anode; full
  precision radii; prints the clearance distribution and charge per radial interval and fails unless coverage is
  complete. (Fix 2026-11-24: the old lower bound 15.05 um excluded 1.99% of the charge of muon 0000.)
* `test_edge_coverage.sh LUT NUMW [DIR]` — SYNTHETIC regression: ions 0.001 um above the wire are covered, the
  strict no-extrapolation check still rejects the old 15.05 um grid.
* `s_green_table.py` — Green table with the frozen LUT and the replay's own path validation.
* `s_synthesize.py` — the variants above (`--variants G` for closure only, `--suffix _TAG`); consistency checks.
* `s_validate_transport.py` — grid executable vs frozen exporter on the same start points (`validate`).
* `run_s_ladder.sh` — check / prepare / validate / green / closure / closure-compare / synth / compare.
* `synthetic_paths.py` — SYNTHETIC toy transport for pipeline tests only.

## WSL (after review)
```bash
cd ~/na6p/wt-validation && git pull
SA=test/mwpc-stage-c/cathode-ion-reference/spatial-ablation
source ~/t_ladder_deps.env            # MWPC_EXPORT_BIN, MWPC_ION_MOBILITY, MWPC_SPECTRAL_LUT, MWPC_NUMW_DIR
source ~/na6p/env.sh 2>/dev/null; source ~/na6p/install/garfieldpp/share/Garfield/setupGarfield.sh 2>/dev/null
cmake -S $SA -B ~/na6p/build-spatial -DFROZEN_DIR=$HOME/na6p/wt-wire/test/mwpc-stage-c/cathode-ion-reference
cmake --build ~/na6p/build-spatial -j4         # prints the frozen exporter source sha256
export MWPC_GREEN_BIN=~/na6p/build-spatial/p3_green_grid_paths
bash $SA/run_s_ladder.sh check
bash $SA/run_s_ladder.sh prepare               # first non-sealed development muon; prints RUN DIR
bash $SA/run_s_ladder.sh validate <RUN>        # 200 archived start points: grid tool must reproduce the frozen paths
bash $SA/run_s_ladder.sh green <RUN>           # TAG r14t72
bash $SA/run_s_ladder.sh closure <RUN>
bash $SA/run_s_ladder.sh closure-compare <RUN> r14t72
# only if the closure leaves a decision ambiguous, one refinement, e.g.:
#   MWPC_S_NR=20 bash $SA/run_s_ladder.sh green <RUN>; MWPC_S_NR=20 bash $SA/run_s_ladder.sh closure <RUN>
#   bash $SA/run_s_ladder.sh closure-compare <RUN> r14t72,r20t72
# when the closure is interpretable:
bash $SA/run_s_ladder.sh synth <RUN>
bash $SA/run_s_ladder.sh compare <RUN>
```
Return `transport_validation_*/result.txt`, `closure_*/t_ladder_summary.txt` (+ PNGs), `grid_spec_*.txt`,
`green_*/axes.json`, `green_paths_*.log`, `timing.txt`, `muon_*/s_synthesize_*.json`, `provenance_prepare_*.txt`.

## Synthetic pipeline test (toy transport, real spectral LUT; NOT physics)
One toy muon (960 ions, 40 primaries), grid 3 wires x 12 radii x 72 angles, D step 0.17 mm, 200 tau nodes; table
439 MiB, 125 s on a 2-core VM; all 11 variants synthesized in 2.6 s. Tp = 100 ns, ion only, L2 vs S0 (U / V):
G closure 0.16% / 0.31%; S1 10% / 9.3%; S2 0.68% / 0.78%; S3 0.65% / 0.60%; SR 2.6% / 5.7%; SW 0.16% / 0.31%;
K1 3.0% / 5.7%. Here a 0.3% closure is ample to reject S1 or SR-type errors but would not separate S2 from S3
(0.7% vs 0.6%), which is exactly the case where the closure procedure above asks for one refinement. These numbers only show that the tools run and resolve differences; the toy radial time law is not physics.

## Dependency separation for later design studies
The S-ladder already factorizes the response the way a configurable engine needs:
* **Transport** (gas, E/B field, wire geometry, voltage): the grid paths `green_paths_*`. Valid only for the frozen
  9-wire 2100 V zero-B field; its 2D symmetry (path independent of w and t) is what allows one path per (k, r, theta).
* **Electrode weighting** (strip pitch, width, gap, stereo angle, segmentation): the LUT applied in `s_green_table.py`.
  New strip geometry with the same drift field = new LUT + rebuilt table from the SAME grid paths.
* **Electronics** (impulse response, thresholds): applied only in `t_shaped_compare.py`; G/PHI0 and all variant
  waveforms are unshaped, so a new shaper reuses everything.
* **Placement / incidence**: enters only through the ion start points (archived or generated) and the strip frame.
A B field, a 3D field or finite strips break the w-translation and/or 2D path symmetry; then paths must be indexed by
the extra coordinates. The frozen P3 tools stay hard-coded; the production engine should key caches on explicit,
versioned configuration hashes (transport, weighting, electronics).

## Limits
* Same frozen-model limits as the T-ladder: ideal 9-wire plane, provisional Ar+ mobility, ideal contiguous strips,
  10 ns replay grid, no Heed primary ions, approximate electron step.
* The grid covers the wires and radii of the prepared muons only; adding muons may need a wider grid (`green` again).
* `p3_green_grid_paths.cc` is syntax-checked here against stub Garfield/ROOT headers only; first real build is on WSL.
