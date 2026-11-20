# Stage C: controlled 100ns vs 10ns ion-current resolution test

**Scope:** same 2100V P3 noFW analytic *whole-plane* field, same frozen K512 muon ROOT archive, same positive-ion mobility surrogate and exactly the same 256 training-ion starting positions. No new Geant4/Heed/AvalancheMicroscopic runs, no ion-reservoir resampling and no stereo-strip/electronics calculation. Completed e+ion benchmark has an audited 0.05ns electron ROOT which can be added for first-100ns charge context.

## Why

The previous K512 same-reservoir 100ns-binned fast ion model yielded +87.6327 vs +87.6082 fC upper-plane integrated charges, but **first-100ns direct vs fast ion charges** were ~+9.81679 vs ~+7.83016 fC. Much of this error compensated by 200-400ns. It is scientifically incorrect to attribute all of it to *spatial* kernel interpolation: shifting piecewise-averaged 100ns single-ion impulse templates by fractional bins can smear the steep initial current. The reference discretization itself must also be checked.

## Decision gates

1. **Reference binning:** recompute K512 direct DriftLineRKF waveform at 10ns for 100us; charge-conservatively rebin ten adjacent fine bins into each original 100ns bin. Compare against the old 100ns direct reference. If not close, we have a reference-resolution issue and must not claim a fast kernel timing improvement.
2. **Same training sites:** rebuild all 256 old ion template positions at 10ns with direct RKF (do NOT resample avalanche ions or train PCA). Output a CSV comparing old 100ns template bins to new templates rebinned to 100ns, per original ion site.
3. **Same reservoir, same spatial interpolation:** reconstruct the K512 event from four-neighbor raw IDW at 10ns; keep the original exact weighted 3072 ion births and their t_birth_ns, and the same mirrored whole-plane symmetry.
4. **Early error:** compare first 100/200/400/1000ns integrated ion charges and all-window L2 on both planes for 100ns and 10ns. An optional electron ROOT adds microscopic electron charge from the same K512 muon for context.
5. **Interpret:** if first-100ns error collapses after proper fine integration, time discretization was dominant. If not, next isolate position-template coverage and ion creation-time discretization, without blindly raising template count or adding NN.

Important: the 10ns rebuild may produce slightly different Sensor numerical integration results even after rebinning. ROOT direct and kernel are **both ion-only**. We must not call any result VMM shaped; selected 100ns VMM3a peaking is not a boxcar integral.

## Commands — independent worktree, no interference with overnight 1000 muons

Never switch, pull, clean, rebuild or alter `~/na6p/cathode-overnight-2100` while its checksum-locked 1000-muon job is running. Create a separate worktree solely for this analysis.

~~~bash
cd ~/na6p/cathode-overnight-2100
git fetch origin feature/p3-2100v-10ns-ion-timing-20261120
git worktree add --detach ~/na6p/cathode-timing-2100 FETCH_HEAD
cd ~/na6p/cathode-timing-2100
export GARFIELD_INSTALL="$HOME/na6p/install/garfieldpp"
export LD_LIBRARY_PATH="$GARFIELD_INSTALL/lib:${LD_LIBRARY_PATH:-}"
export HEED_DATABASE="$GARFIELD_INSTALL/share/Heed/database"
cmake -S test/mwpc-avalanche -B test/mwpc-avalanche/build-timing \
  -DCMAKE_PREFIX_PATH="$GARFIELD_INSTALL"
cmake --build test/mwpc-avalanche/build-timing --target \
  mwpc_plane_response_from_ions_root \
  mwpc_rebuild_plane_templates_fine \
  mwpc_plane_muon_kernel_compare \
  mwpc_compare_plane_timing_100ns_10ns -j2
~~~

Set exact file paths:

~~~bash
BIN="$PWD/test/mwpc-avalanche/build-timing"
DATA="$HOME/na6p/cathode-ion-2100-wsl/test_runs/mwpc_stage_c/cathode_2100_laptop_v1"
IONS="$DATA/weighted_k512/muon_0000/ion_births.root"
MOBILITY="$HOME/na6p/src/garfieldpp/Data/IonMobility_Ar+_Ar.txt"
OLDTRAIN="$HOME/na6p/cathode-compact-kernel-2100/test_runs/mwpc_stage_c/polar_kernel_v3/train256_test48_basis.root"
DIRECT100="$HOME/na6p/cathode-plane-response-2100/test_runs/mwpc_stage_c/plane_response_pilot_v1/muon0_weighted512_direct"
FAST100="$HOME/na6p/cathode-compact-kernel-2100/test_runs/mwpc_stage_c/muon_plane_waveform_v1/jpsi_muon0_K512"
ELECTRONS="$HOME/na6p/cathode-compact-kernel-2100/test_runs/mwpc_stage_c/whole_muon_electron_plus_ion_v1/muon0_direct_electron_retry1.root"
OUT="$PWD/test_runs/mwpc_stage_c/plane_timing_dt10_v1"
mkdir -p "$OUT"
for P in "$IONS" "$MOBILITY" "$OLDTRAIN" \
  "${DIRECT100}_waveforms.csv" "${DIRECT100}_summary.csv" \
  "${FAST100}_waveforms.csv" "${FAST100}_summary.csv"; do
  test -s "$P" || { echo "MISSING $P"; exit 1; }
done
~~~

Step 1: direct weighted RKF at 10ns (same exact original sampled macroion population):

~~~bash
"$BIN/mwpc_plane_response_from_ions_root" \
  --ion-file "$IONS" --out-prefix "$OUT/muon0_K512_direct_dt10" \
  --mobility "$MOBILITY" --mode direct \
  --dt-ns 10 --window-us 100 --max-step-mm 0.05
~~~

Step 2: rebuild ONLY the same 256 original ion-response templates, new 10ns binning:

~~~bash
"$BIN/mwpc_rebuild_plane_templates_fine" \
  --input-training-root "$OLDTRAIN" \
  --mobility "$MOBILITY" --out-prefix "$OUT/templates256_dt10" \
  --dt-ns 10 --window-us 100 --max-step-mm 0.05
~~~

Step 3: use exactly the same K512 ions and 4-neighbor raw IDW template interpolation at 10ns:

~~~bash
"$BIN/mwpc_plane_muon_kernel_compare" \
  --ion-file "$IONS" \
  --training-root "$OUT/templates256_dt10_training.root" \
  --reference-waveforms "$OUT/muon0_K512_direct_dt10_waveforms.csv" \
  --reference-summary "$OUT/muon0_K512_direct_dt10_summary.csv" \
  --out-prefix "$OUT/muon0_K512_fast_dt10" \
  --neighbors 4 --dt-ns 10 --window-us 100
~~~

Step 4: numerical comparison and standalone ROOT+PNG early-pulse diagnostic:

~~~bash
EXTRA=()
if test -s "$ELECTRONS"; then EXTRA=(--electron-root "$ELECTRONS"); fi
"$BIN/mwpc_compare_plane_timing_100ns_10ns" \
  --direct-100ns "${DIRECT100}_waveforms.csv" \
  --event-100ns "${FAST100}_waveforms.csv" \
  --direct-10ns "$OUT/muon0_K512_direct_dt10_waveforms.csv" \
  --event-10ns "$OUT/muon0_K512_fast_dt10_waveforms.csv" \
  "${EXTRA[@]}" \
  --out-prefix "$OUT/early_timing_compare"
cat "$OUT/early_timing_compare_windows.csv"
cat "$OUT/early_timing_compare_metrics.csv"
ls -lh "$OUT/early_timing_compare_early_overlay.png"
~~~

Files are written to the NEW worktree only, old 100ns sources and the active 1000-muon campaign are READ ONLY. Each run refuses to overwrite an existing prefix.

### Measurements to report

- reference time-grid error: direct100 minus direct10->100 per cathode and integrated time window;
- fast temporal-response error: fast10 minus direct10->100, compared with original fast100 minus direct100;
- current L2 for 10ns and rebinned 100ns on each plane;
- optional complete-muon electron charge for first 100ns and full window (same e data used before, not recalculated);
- verified same original crossing, true ion-birth count, retained ion count, station/det, 100us window;
- provisional Ar+ in pure argon mobility, fixed-field, no strips/no electronics remain unchanged.
