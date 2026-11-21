# Gate 5: 12x24x4 phase-resolved deterministic polar single-ion operator

**Status 2026-11-21:** C++/Python code committed, new ROOT generation,
new Python smoke test, and real-data replay not yet run on WSL.

## Numerical evidence

On the originally examined development muon0301, the existing
12x24 uncompressed polar (ln r,theta) response bank gave
native 10ns two-plane event L2 **4.211960%** against direct RKF.
The old 256-site IDW4 benchmark gave **6.230371%**.

A separate exact-spatial audit recomputing every original weighted
ion at its actual t_birth matched the frozen Garfield RKF reference
within L2 = 1.154866e-16. The SAME exact ion positions drifted
at t=0 and shifted via 10ns fractional-bin linear interpolation
yielded 3.997625% full-waveform L2, with first100ns L2=4.293627%.

Using 10ns Sensor response templates sampled at t_birth phase
0,2.5,5,7.5 ns (and integer-shifted 0ns waveform for 10ns),
the exact-spatial reference gave **0.924936% full L2** and
**0.993818% first100ns L2**; full100us integrated charge errors
were numerical zero. Late (10–100us) phase-only L2 was 0.001733%.

These observations motivate **joint spatial plus temporal**
interpolation; no claim yet that it will reach 0.925% on full muon.

## Experiment

Modify the already validated deterministic polar C++ grid generator
to accept --phase-subdivisions {1,2,4,8}, default 1 for backward
compatibility. For Nphase=4 it simulates 12*24*4 = **1152**
separate Garfield DriftLineRKF ion starts, all at the same frozen
central-wire (r,theta) coordinates but with phases 0,2.5,5,7.5ns.

The script replay_phase_resolved_polar_bank.py uses identical:

- two whole-plane plus/minus Garfield electrode responses,
- original complete 4096-weighted-macroion source and weights,
- exact original muon0301 direct-RKF CSV and summary,
- fractional birth time, but represented by neighboring Garfield
  phase responses and an integer-bin shift,
- ordinary uncompressed bilinear (ln r,theta) spatial interpolation.

The phase0 templates are first compared directly against the existing
independently generated 12x24 t=0 grid. Require relative L2 < 1e-8;
otherwise STOP before studying muon accuracy.

For a birth t/10ns=m+f, let k=floor(4f), eta=4f-k.
Deposit (1-eta) in phase k, integer bin m, and eta in phase k+1
at m, except phase 4 -> phase0 at integer bin m+1.
Whole-event signal from FFT convolution of sparse phase/node moment
streams with exact stored phase-node current waveforms.

A separate synthetic --self-test must show the FFT scheme exactly
matches an independent scalar double-linear interpolation with
mixed positive and negative v and all 10ns fractional-shift rules.

## Most important error decomposition

For EXACT SAME muon reference R, predicted phase4 spatial grid P,
and phase4 replay evaluated at exact ion positions S:

    P-R = (P-S) + (S-R)

and the native waveform relative squared error satisfies

    ||P-R||²/||R||² = ||P-S||²/||R||²
                       + ||S-R||²/||R||²
                       + 2 <P-S,S-R>/||R||².

The cross term must be retained: direct subtraction or quadrature
of percentages is not a valid decomposition without orthogonality.
This controls spatial interpolation error separately from
interpolation between precomputed birth phases.

## WSL commands (new outputs only)

```bash
cd "$HOME/na6p/response-rank-pilot-20261121"
git pull --ff-only
source .venv/bin/activate

ROOT_INSTALL="$HOME/na6p/install/root"
GARFIELD_INSTALL="$HOME/na6p/install/garfieldpp"
source "$ROOT_INSTALL/bin/thisroot.sh"

PY=test/mwpc-stage-c/cathode-ion-reference/replay_phase_resolved_polar_bank.py
python -m py_compile "$PY"
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python "$PY" --self-test

cmake -S test/mwpc-avalanche \
  -B test/mwpc-avalanche/build-polar \
  -DROOT_DIR="$ROOT_INSTALL/cmake" \
  -DCMAKE_PREFIX_PATH="$ROOT_INSTALL;$GARFIELD_INSTALL"
cmake --build test/mwpc-avalanche/build-polar \
  --target mwpc_deterministic_polar_plane_bank -j2

BIN="$PWD/test/mwpc-avalanche/build-polar/mwpc_deterministic_polar_plane_bank"
MOBILITY="$HOME/na6p/src/garfieldpp/Data/IonMobility_Ar+_Ar.txt"
PILOT="$PWD/test_runs/mwpc_stage_c/polar_grid_2100_v1"
GRID="$PILOT/grid12x24_phase4_dt10.root"

test ! -e "$GRID" || { echo 'GRID EXISTS, STOP'; exit 1; }
"$BIN" --mobility "$MOBILITY" --out-root "$GRID" \
  --n-radius 12 --n-angle 24 --phase-subdivisions 4 \
  --clearance-min-um 0.00001 --radius-max-um 80 \
  --theta-min-deg 0.001 --dt-ns 10 --window-us 100 \
  --max-step-mm 0.05 > "$PILOT/grid12x24_phase4_generation.log" 2>&1
tail -n 5 "$PILOT/grid12x24_phase4_generation.log"
```

Then only if 1152 RKF ion starts finished without errors:

```bash
DATA="$HOME/na6p/cathode-overnight-2100/test_runs/mwpc_stage_c/overnight_2100_weighted_k256_v1/muon_0301/ion_births.root"
BASE="$HOME/na6p/cathode-independent-stereo-2100/test_runs/mwpc_stage_c/independent_muon_dt10_v1"
OLD="$PILOT/replay_muon0301_12x24_windowed_v2/complete_waveforms.csv"
PHASE="$PILOT/exact_spatial_phase4_muon0301_v1_comparison.csv"
OUT="$PILOT/replay_muon0301_12x24_phase4_v1"

for file in "$GRID" "$DATA" "$OLD" "$PHASE";do
  test -s "$file" || { echo "MISSING $file"; exit 1; }
done
test ! -e "$OUT" || { echo 'OUTPUT EXISTS, STOP'; exit 1; }
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python "$PY" \
  --grid-root "$GRID" \
  --original-grid-root "$PILOT/grid12x24_dt10.root" \
  --ion-file "$DATA" \
  --reference-waveforms "$BASE/muon_0301_direct_dt10_waveforms.csv" \
  --reference-summary "$BASE/muon_0301_direct_dt10_summary.csv" \
  --baseline-waveforms "$BASE/muon_0301_fast_dt10_waveforms.csv" \
  --old-polar-waveforms "$OLD" \
  --exact-spatial-phase-comparison "$PHASE" \
  --out-dir "$OUT"

cat "$OUT/event_metrics.csv"
cat "$OUT/window_metrics.csv"
cat "$OUT/metadata.json"
```

## Caveats

The output waveform error is for *two ideal whole cathode planes*
of the provisional 2100V drift model, given ALL sampled ion births.
It does not establish U/V stereo spatial resolution or compress
the microscopic avalanche representation. The two untouched
muons 0408,0446 are not evaluated and remain final holdouts.

Do not assume the temporal operator needs 4 separate full waveform
banks in a deployed C++ algorithm. An adaptive high-time-resolution
cumulative induced-charge representation may be much more compact;
we should investigate that only after this clean accuracy gate.
