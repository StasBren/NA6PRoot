# T-ladder timing ablation (2100 V, frozen RKF paths + frozen spectral replay)

Question: how far can avalanche ion-creation times be aggregated without changing the shaped, signed U/V strip
signals at (provisional) 100 and 200 ns peaking time?

The ion paths come from the frozen exporter (`mwpc_stage_c_p3_export_full_spectral_paths`, Garfield RKF, static
field, mobility table given on its command line). In a static field with a time-independent mobility a path depends
only on its start point, so collapsing a birth time is an exact rigid time shift of that path. No new RKF, no new
avalanche, no change to the replay or the LUT.

| Variant | Export transform (group -> one ion-yield-weighted mean birth time) | Ladder rung |
|---|---|---|
| T0 | identity; byte-identical to the frozen export (checked) | microscopic reference |
| T1 | Heed primary (each avalanche keeps its ion positions; primaries keep distinct times) | T1 |
| TW | anode nearest each ION BIRTH (`round(u/3 mm)` = archive `local_wire_index`) | one time per wire |
| TC | DOMINANT COLLECTING anode of the ion's primary (exact census `microscopic_wire_seed.csv`) | one time per collecting wire |
| TM | whole crossing (optional, extreme control) | one time per muon |

TW and TC differ only for ions whose birth-nearest anode is not their primary's dominant collecting anode; V0 reported
52 / 42,879 primaries at the dominant-wire level (count from before the zero-count argmax fix; recheck). TC records the local-to-global wire offset it derived and
the weighted fraction of ions it moves.

Signal, in this order:
1. **Ion only**: the frozen replay. For the frozen RKF paths and spectral LUT this is the time-dependent
   Shockley–Ramo ion signal (chargeSign = +1).
2. **Ion + instantaneous electron step** (`t_prompt_electron_term.py`): the avalanche-electron signal replaced by its
   time integral `+q phi_k(r_birth)` in the birth-time bin. An APPROXIMATION, not an electron-trajectory calculation;
   the archives hold no avalanche-electron trajectories.

Shaping: unit-peak CR-RC^3 **stand-in** (`tau = Tp/3`) unless `--impulse` gives a real VMM3a response. 100 and 200 ns
are PROVISIONAL until the NA60+ front-end configuration and transfer function are verified.

Comparison (`t_shaped_compare.py`): normalized L2 of signed shaped waveforms per plane, leading-strip peak and peak
time, leading-edge times at thresholds f x (T0 lead peak) with f = 0.02, 0.05, 0.1, 0.2, 0.5 applied as the SAME
absolute value to every variant, strip multiplicity and hit flips per threshold, all-strip peak changes, plane charge
and amplitude, U/V ratios, stereo w and u from U/V centroids. `--plots` writes per-muon waveform and residual figures.

## Dependencies (none rebuilt silently)
| Item | Where it comes from | If missing |
|---|---|---|
| exporter | target `mwpc_stage_c_p3_export_full_spectral_paths`, `test/mwpc-stage-c/CMakeLists.txt` on `feature/p3-2100v-wire-aware-stereo-weighting-20261120` | compile that unchanged source (Garfield + ROOT); compiling does not change physics |
| mobility | `$HOME/na6p/src/garfieldpp/Data/IonMobility_Ar+_Ar.txt` (named in COMPACT_KERNEL_PILOT.md, PLANE_RESPONSE_PILOT.md, MUON_PLANE_WAVEFORM.md as the provisional Ar+ in Ar table) | do not substitute; any other table changes the physics |
| replay + LUT builder | `test/mwpc-stage-c/numerical-weighting/` on the same wire branch | -- |
| spectral LUT | an existing `.npz` with `meta.model = P3_9wire_infinite_plates_ideal_contiguous_strips` | `p3_spectral_weighting_lut.py build NEW.npz` with defaults: same model and grid (meta printed by `check`); not a new physics assumption |
| previous frozen replay outputs | `waveforms.csv` + `channels.csv` + `summary.json` | optional; with one, `MWPC_FROZEN_REPLAY_DIR` makes T0 a byte-level reproduction gate |

`locate_t_ladder_deps.py` scans `~/na6p` and `~/NA6PRoot*`, lists every candidate with SHA256 (and rows / exact_ions
for previous exports and replays so they can be paired with a muon), and writes `t_ladder_deps.env` for review.

## WSL
```bash
cd ~/na6p/wt-validation && git pull                         # branch dev/na60-2100v-validation-20261123
TA=test/mwpc-stage-c/cathode-ion-reference/timing-ablation
(cd ~ && python3 ~/na6p/wt-validation/$TA/locate_t_ladder_deps.py)   # writes ~/t_ladder_deps_report.txt, ~/t_ladder_deps.env
# wire branch (exporter, replay, LUT builder), unchanged:
git -C ~/na6p/src/NA6PRoot-mwpc fetch origin feature/p3-2100v-wire-aware-stereo-weighting-20261120
git -C ~/na6p/src/NA6PRoot-mwpc worktree add ~/na6p/wt-wire FETCH_HEAD
# only if no exporter was found (same Garfield/ROOT environment as the 2100 V builds):
source ~/na6p/env.sh 2>/dev/null; source ~/na6p/install/garfieldpp/share/Garfield/setupGarfield.sh 2>/dev/null
cmake -S ~/na6p/wt-wire/test/mwpc-stage-c -B ~/na6p/wt-wire/build-stage-c && \
  cmake --build ~/na6p/wt-wire/build-stage-c --target mwpc_stage_c_p3_export_full_spectral_paths -j4
source ~/t_ladder_deps.env                                  # after review; or export the four variables by hand:
#   MWPC_EXPORT_BIN  MWPC_ION_MOBILITY  MWPC_SPECTRAL_LUT  MWPC_NUMW_DIR=~/na6p/wt-wire/test/mwpc-stage-c/numerical-weighting
bash $TA/run_t_ladder.sh check
bash $TA/run_t_ladder.sh one                                # first non-sealed development muon, T0,T1,TW,TC
bash $TA/run_t_ladder.sh compare <RUN DIR printed by one>   # 1) ion only  2) + electron step; plots included
# only after the one-muon result is reviewed:
bash $TA/run_t_ladder.sh run 5 && bash $TA/run_t_ladder.sh compare <RUN DIR>
```
Return both `compare_*/t_ladder_summary.txt`, the PNGs, `provenance.txt`, `muon_*/gates.txt`, `muon_*/timing.txt`
and `muon_*/T?/export/ablation.json`.

## Limits
* 10 ns replay grid: threshold times finer than a few ns are interpolation, identical for all variants.
* Provisional Ar+ in Ar mobility; ideal nine-wire, infinite-plate, contiguous-strip spectral weighting (not detector
  validated); Heed primary positive ions are not in the archives.
* `make_synthetic_export.py` exists only for pipeline tests; its outputs are not physics.
