# 2100V whole-plane compact ion-kernel pilot

Status: code committed, NOT YET compiled or executed in the user's WSL. Accuracy must be inspected in held-out metrics. NOT real stereo cathode strips, Fast4 model, electrons, electronics, or quantitative ion mobility.

Source: real Garfield++ ion-birth starting positions from the already sealed 2100V weighted K256 ROOT archive. The Garfield++ real P3 nine-wire electric field and the two WHOLE cathode weighting fields are used to simulate individual positive ions with DriftLineRKF. Birth times are deliberately reset to t=0 for this isolated impulse-response study (original birth time and primary are saved in the training bank).

Three methods evaluated on the same DISJOINT held-out ion starts:
1. RAW k-nearest ion-response templates, inverse-distance interpolation (IDW, same hemisphere).
2. PCA temporal basis learned on the training RKF waveforms, rank 2/4/8/16 when numerically meaningful. Predict held-out coefficients via spatial IDW. This is the actual fast model candidate.
3. PCA oracle: project the TRUE held-out RKF waveform onto the PCA basis. This diagnoses compression alone, and is NOT a deployable predictor.

Source sampling automatically restricts the pilot to the major ion hemisphere if the sample does not support independent training/test on both; prints actual coverage. This makes no false claim to validate the absent minority hemisphere.

The error metrics are relative L2 of BOTH current waveforms; peak absolute cumulative-charge difference measured in elementary charges; absolute final charge differences for both planes; and early cumulative charge at 1 us. CSV reports held-out mean/median/90th percentile/max. A process PASS means calculations completed, NOT that accuracy passed a scientific threshold.

## Commands (isolated WSL clone: never modify previous locked microscopic data)

~~~bash
cd ~/na6p
git clone --depth 1 \
  --branch feature/p3-2100v-compact-kernel-pilot-20261118 \
  https://github.com/StasBren/NA6PRoot.git cathode-compact-kernel-2100
cd cathode-compact-kernel-2100
export LD_LIBRARY_PATH="$HOME/na6p/install/garfieldpp/lib:$LD_LIBRARY_PATH"

cmake -S test/mwpc-avalanche -B test/mwpc-avalanche/build-kernel \
  -DCMAKE_PREFIX_PATH="$HOME/na6p/install/garfieldpp"
cmake --build test/mwpc-avalanche/build-kernel \
  --target mwpc_compact_plane_kernel_pilot -j2

BIN=./test/mwpc-avalanche/build-kernel/mwpc_compact_plane_kernel_pilot
MOBILITY="$HOME/na6p/src/garfieldpp/Data/IonMobility_Ar+_Ar.txt"
DATA="$HOME/na6p/cathode-ion-2100-wsl/test_runs/mwpc_stage_c/cathode_2100_laptop_v1"
OUT="$PWD/test_runs/mwpc_stage_c/compact_kernel_pilot_v1"
mkdir -p "$OUT"

"$BIN" \
  --ion-file "$DATA/weighted_k256/muon_0000/ion_births.root" \
  --mobility "$MOBILITY" \
  --out-prefix "$OUT/jpsi0_K256_96train_48test" \
  --train 96 --test 48 --neighbors 4 \
  --seed 210018 --dt-ns 100 --window-us 100 --max-step-mm 0.05

cat "$OUT/jpsi0_K256_96train_48test_summary.csv"
~~~

Outputs: <prefix>_summary.csv; <prefix>_metrics.csv (per ion); <prefix>_basis.root with ROOT histograms of eigenvalues, mean and PCA time shapes plus a compact TTree training_bank containing birth positions, original primary/time, basis coefficients, and training waveforms. They are small and do not duplicate microscopic ROOT archives. No overwrite of results is allowed.

Interpretation:
- Small oracle error, large interpolated error: work on spatial coordinates/adaptive sampling rather than increasing rank.
- Large oracle error: PCA rank/time representation inadequate.
- Raw IDW poor too: more local coverage, physics-informed templates, or conditional modeling needed.
- A good held-out ion result is only a pilot. Next compare entire muon waveforms and compare K256/K512, then independent channels/stations and full ROOT truth, and later stereo strip weighting fields.
- No grid/kernel is assumed correct for real Ar/CO2 positive ion drift time because Ar+ in pure Ar is a PROVISIONAL mobility surrogate.
- Existing frozen ROOT references and historical 371MB C2 numerical weighting directory are not modified or deleted.

## Fixed held-out learning curve (2026-11-18)

Original 96-train / 48-test pilot showed 2.41% raw-IDW mean waveform L2, 2.40% PCA-8 with IDW and 0.83% PCA-16 oracle. Interpolation rather than temporal rank is dominant. However the original holdout changes with --train and therefore must not be compared to a different training size.

New source selects held-out ions FIRST and keeps the test IDs fixed when --train is changed but archive, --test and --seed are unchanged. It freezes hemisphere policy using the maximum allowed training size 256. Rerun 96 in a fresh output prefix.

~~~bash
cd ~/na6p/cathode-compact-kernel-2100
git pull --ff-only
cmake --build test/mwpc-avalanche/build-kernel --target mwpc_compact_plane_kernel_pilot -j2
BIN=./test/mwpc-avalanche/build-kernel/mwpc_compact_plane_kernel_pilot
MOBILITY="$HOME/na6p/src/garfieldpp/Data/IonMobility_Ar+_Ar.txt"
DATA="$HOME/na6p/cathode-ion-2100-wsl/test_runs/mwpc_stage_c/cathode_2100_laptop_v1"
OUT="$PWD/test_runs/mwpc_stage_c/kernel_learning_curve_v2"
mkdir -p "$OUT"
for N in 96 192 256; do
  "$BIN" --ion-file "$DATA/weighted_k256/muon_0000/ion_births.root" \
    --mobility "$MOBILITY" --out-prefix "$OUT/train${N}_test48" \
    --train "$N" --test 48 --neighbors 4 --seed 210018 \
    --dt-ns 100 --window-us 100 --max-step-mm 0.05
done
for N in 96 192 256; do
  echo "=== train=$N ==="
  cat "$OUT/train${N}_test48_summary.csv"
done
diff -u \
  <(awk -F, 'NR>1{print $1}' "$OUT/train96_test48_metrics.csv" | sort -nu) \
  <(awk -F, 'NR>1{print $1}' "$OUT/train256_test48_metrics.csv" | sort -nu)
~~~

An empty diff confirms the same held-out source rows. If interpolation error shrinks sharply with N, source coverage matters. If it plateaus well above the PCA oracle, investigate physically motivated 2D interpolation or conditional parametric response. Further independent muon/channel holdouts and real mixture ion mobility remain necessary.
