# T-ladder timing ablation (2100 V, frozen RKF paths + frozen spectral replay)

Question: can ion-birth times be collapsed (per primary, per wire, per muon) without changing the shaped, signed U/V
strip signals at 100 and 200 ns peaking time?

The ion paths come from the frozen exporter (`mwpc_stage_c_p3_export_full_spectral_paths`, Garfield RKF, static
field). In a static field with a time-independent mobility the path depends only on its start point, so collapsing a
birth time is an exact rigid time shift of that path. No new RKF, no new avalanche, no change to the replay.

| Variant | Export transform | Ladder rung |
|---|---|---|
| T0 | identity (byte-identical to the frozen export; checked) | reference |
| T1 | all ions of a Heed primary at the primary's ion-weighted mean birth time | T1 |
| TW | all ions on one nearest-birth anode at the wire's mean time | timing part of T3/T4 |
| TM | all ions of the muon at one time | extreme control |

T2, T3 and T4 also replace the spatial cloud with moments. That is a kernel test (harmonic kernel vs this replay) and
is the next experiment; it needs no new tools beyond the kernel itself.

Signal, in this order:
1. **Ion only**: the frozen replay. For the frozen RKF paths and spectral LUT this is the time-dependent
   Shockley–Ramo ion signal.
2. **Ion + instantaneous electron step** (`t_prompt_electron_term.py`): the avalanche-electron signal replaced by its
   time integral `+q phi_k(r_birth)` in the birth-time bin (frozen convention chargeSign = +1). This is an
   APPROXIMATION of the true time-dependent electron current; the archives hold no avalanche-electron trajectories. Shaping:
unit-peak CR-RC^3 **stand-in** (`tau = Tp/3`) unless `--impulse` gives a real VMM3a response; the VMM3a user guide
specifies a 3rd-order shaper (1 real + 2 complex-conjugate poles) and 25/50/100/200 ns peaking times but no pole
values. 100 and 200 ns are PROVISIONAL settings until the NA60+ front-end configuration is verified.

## Files
* `t_time_collapse_export.py SRC OUT --mode identity|primary|wire|muon`: writes `ablation.json` (source SHA256, RMS shift).
* `t_prompt_electron_term.py TABLE EXPORT OUT --replay-dir NUMW`: electron term on the replay grid and channels.
* `t_shaped_compare.py ROOT --out DIR`: per-strip peaks, peak times, threshold times, multiplicities, centroids,
  w/u; paired differences vs T0 with a muon bootstrap. `--no-prompt` = ion term only.
* `run_t_ladder.sh check|run [N]|compare DIR`: development muons only; sealed slots are never touched.
* `make_synthetic_export.py`: SYNTHETIC straight-line paths, only for pipeline tests.

## WSL
```bash
cd <worktree of this branch>
# The exporter and replay live on feature/p3-2100v-wire-aware-stereo-weighting-20261120 (not on this branch); use a checkout of that branch, unchanged:
WIRE=<checkout of the wire branch>
export MWPC_EXPORT_BIN=<build>/mwpc_stage_c_p3_export_full_spectral_paths
export MWPC_ION_MOBILITY=<same mobility file as the frozen replays>
export MWPC_NUMW_DIR=$WIRE/test/mwpc-stage-c/numerical-weighting
export MWPC_SPECTRAL_LUT=<lut>.npz   # or: (cd $MWPC_NUMW_DIR && python3 p3_spectral_weighting_lut.py build ~/lut_default.npz)  ~4 min, 21 MiB
R=test/mwpc-stage-c/cathode-ion-reference/timing-ablation/run_t_ladder.sh
bash $R check
export MWPC_FROZEN_REPLAY_DIR=<existing frozen replay output of that muon, if any>   # optional T0 identity gate
bash $R one            # FIRST: one development muon, variants T0,T1,TW
bash $R compare <RUN DIR printed by one>   # ion-only first, then + electron step
# only after reviewing one muon:
bash $R run 5           # 5 muons x 3 variants (T0,T1,TW; add TM via MWPC_TLADDER_VARIANTS). Replay wall time scales with rows: measured 19.5-21.2 s for 1,200 synthetic rows on a 2-core cloud VM; a real K256 muon has ~13,500 rows (~53 primaries x 256), so expect minutes per variant
bash $R compare <RUN DIR printed by run>
```
Return `compare_*/t_ladder_summary.txt`, `t_ladder_summary.json`, `provenance.txt` and the `*/ablation.json` files.

## Limits
* 10 ns replay grid: threshold times finer than a few ns are interpolation, identical for all variants.
* The replay is ion-only with provisional Ar+ mobility; the prompt helper adds the electron step at birth time.
* Heed primary positive ions are not in the archives (listed in the approximation register).
