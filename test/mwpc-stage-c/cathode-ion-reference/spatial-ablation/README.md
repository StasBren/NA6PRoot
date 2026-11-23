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

**Gate:** G must agree with S0 far better than the effects under test (target L2 < 1e-3 and threshold-time shifts
well below 1 ns). If not, refine the grid (MWPC_S_NR, MWPC_S_NTHETA, `--d-sub`, `--ntau`) before reading S1..K2.
S5 (shaped kernels and their runtime vs the replay) follows only after the ladder identifies the description.

## Files
* `p3_green_grid_paths.cc`, `CMakeLists.txt` — grid paths with the frozen transport (build against the wire branch).
* `s_grid_spec.py` — wires, radius grid and D range from the muon's actual ion births.
* `s_green_table.py` — Green table with the frozen LUT and the replay's own path validation.
* `s_synthesize.py` — the variants above; writes replay-format waveforms and electron-step files.
* `run_s_ladder.sh` — check / prepare / green / synth / compare; development muons only.
* `synthetic_paths.py` — SYNTHETIC toy transport for pipeline tests only.

## WSL (after review)
```bash
cd ~/na6p/wt-validation && git pull
SA=test/mwpc-stage-c/cathode-ion-reference/spatial-ablation
source ~/t_ladder_deps.env            # MWPC_EXPORT_BIN, MWPC_ION_MOBILITY, MWPC_SPECTRAL_LUT, MWPC_NUMW_DIR
source ~/na6p/env.sh 2>/dev/null; source ~/na6p/install/garfieldpp/share/Garfield/setupGarfield.sh 2>/dev/null
cmake -S $SA -B ~/na6p/build-spatial -DFROZEN_DIR=$HOME/na6p/wt-wire/test/mwpc-stage-c/cathode-ion-reference
cmake --build ~/na6p/build-spatial -j4
export MWPC_GREEN_BIN=~/na6p/build-spatial/p3_green_grid_paths
bash $SA/run_s_ladder.sh check
bash $SA/run_s_ladder.sh prepare                # first non-sealed development muon; prints RUN DIR
bash $SA/run_s_ladder.sh green   <RUN DIR>      # ~3 wires x 14 radii x 72 angles RKF paths + table
bash $SA/run_s_ladder.sh synth   <RUN DIR>
bash $SA/run_s_ladder.sh compare <RUN DIR>
```
Return `compare_*/t_ladder_summary.txt`, the PNGs, `grid_spec.txt`, `green/axes.json`, `green_paths.log`,
`timing.txt`, `muon_*/S/s_synthesize.json`, `provenance_prepare_*.txt`.

## Synthetic pipeline test (toy transport, real spectral LUT; NOT physics)
One toy muon (960 ions, 40 primaries), grid 3 wires x 12 radii x 72 angles, D step 0.17 mm, 200 tau nodes; table
439 MiB, 125 s on a 2-core VM; all 11 variants synthesized in 2.6 s. Tp = 100 ns, ion only, L2 vs S0 (U / V):
G closure 0.16% / 0.31%; S1 10% / 9.3%; S2 0.68% / 0.78%; S3 0.65% / 0.60%; SR 2.6% / 5.7%; SW 0.16% / 0.31%;
K1 3.0% / 5.7%. The closure is above the 1e-3 target on V: on real data refine the grid before interpreting
variants. These numbers only show that the tools run and resolve differences; the toy radial time law is not physics.

## Limits
* Same frozen-model limits as the T-ladder: ideal 9-wire plane, provisional Ar+ mobility, ideal contiguous strips,
  10 ns replay grid, no Heed primary ions, approximate electron step.
* The grid covers the wires and radii of the prepared muons only; adding muons may need a wider grid (`green` again).
* `p3_green_grid_paths.cc` is syntax-checked here against stub Garfield/ROOT headers only; first real build is on WSL.
