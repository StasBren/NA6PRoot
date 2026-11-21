# 2100V whole-plane ion response — exact-spatial sub-bin phase control

**Status:** code committed; modified C++ target not yet rebuilt/executed in WSL.

## The decisive timing result

On verified frozen phi muon0301 (4096 K256 weighted ions, 886727 true
ion births), re-running actual-time Garfield++ RKF directly reproduced
the frozen whole-plane 10ns/100us reference at L2 = 1.15487e-16.

Using **identical ion positions** but drifting at time zero and
interpolating the 10ns histogram to the actual birth phase yields:

| Time range | Exact positions + naive 10ns phase shift | 12x24 uncompressed polar grid |
|---|---:|---:|
| 0–100ns | 4.29363% | 4.46966% |
| 100ns–1us | 0.34280% | 1.57521% |
| 1–10us | 0.00946% | 0.67997% |
| 10–100us | 0.01018% | 4.98490% |
| entire waveform | 3.99762% | 4.21196% |

Early 100ns charge bias of exact-spatial shifted vs direct is
-0.01129/-0.01106 fC; 100us charge difference ~2e-9 fC, negligible.
328 of 4096 source sampled ions arrive later than 100us in both
recalculated original-time and zero-time RKF runs.

The **3.9976% is a floor for the existing t=0/10ns linear time-shift
procedure**, NOT an intrinsic 10ns-bin data floor or a mathematically
rigorous lower bound for any spatial-grid approximation: different
error components can cancel. The long-time grid error is chiefly
spatial under this control, but do not use a quadrature subtraction
of error magnitudes unless noise/residual orthogonality is checked.

## Next independent check: phase-resolved ion response without spatial error

The original exact-spatial C++ audit now supports
`--phase-subdivisions N` with N = 1, 2, 4 or 8. Default N=1
preserves the previous algorithm, and self-checks exact agreement.

For N=4, it runs the same **real saved ions**, at actual source
birth time for the reference, and also with discrete birth phases
0, 2.5, 5.0, 7.5ns (10ns endpoint obtained by shifting the 0ns
template exactly one bin). Per-ion phase kernels are interpolated
linearly in the phase fraction, followed by the **integer** source
birth-time bin shift. Nothing spatial is interpolated.

The phase-aware replay needs extra Garfield RKF runs on the
already frozen ion positions, but does not rerun Geant4,
TrackHeed, or microscopic avalanches. It is a diagnostic, not
yet a deployable precomputed grid.

**All** 4096 original weighted ions are retained, and the
recomputed original-time waveform must still match the frozen
reference within L2 < 1e-5 before results are accepted.
Columns `phase_corrected_*` in the fresh CSV distinguish this
from the original method.

## WSL commands

```bash
cd "$HOME/na6p/response-rank-pilot-20261121"
git pull --ff-only
source .venv/bin/activate
ROOT_INSTALL="$HOME/na6p/install/root"
GARFIELD_INSTALL="$HOME/na6p/install/garfieldpp"
source "$ROOT_INSTALL/bin/thisroot.sh"
cmake -S test/mwpc-avalanche \
  -B test/mwpc-avalanche/build-polar \
  -DROOT_DIR="$ROOT_INSTALL/cmake" \
  -DCMAKE_PREFIX_PATH="$ROOT_INSTALL;$GARFIELD_INSTALL"
cmake --build test/mwpc-avalanche/build-polar \
  --target mwpc_exact_spatial_time_shift_audit -j2

BIN="$PWD/test/mwpc-avalanche/build-polar/mwpc_exact_spatial_time_shift_audit"
MOBILITY="$HOME/na6p/src/garfieldpp/Data/IonMobility_Ar+_Ar.txt"
DATA="$HOME/na6p/cathode-overnight-2100/test_runs/mwpc_stage_c/overnight_2100_weighted_k256_v1/muon_0301/ion_births.root"
BASE="$HOME/na6p/cathode-independent-stereo-2100/test_runs/mwpc_stage_c/independent_muon_dt10_v1"
OUT="$PWD/test_runs/mwpc_stage_c/polar_grid_2100_v1/exact_spatial_phase4_muon0301_v1"

# Use the already SHA-verified muon 0301, no new archive writes.
test ! -e "${OUT}_summary.csv" && test ! -e "${OUT}_comparison.csv" || {
  echo 'OUTPUT EXISTS: choose a fresh out prefix'; exit 1;
}
"$BIN" --ion-file "$DATA" --mobility "$MOBILITY" \
  --reference-waveforms "$BASE/muon_0301_direct_dt10_waveforms.csv" \
  --reference-summary "$BASE/muon_0301_direct_dt10_summary.csv" \
  --out-prefix "$OUT" --phase-subdivisions 4
cat "${OUT}_summary.csv"
```

If `phase_corrected_full_L2` falls dramatically below ~4%,
the correct next deployment candidate is a **small phase-resolved
template bank** on the existing 12x24 spatial grid (or an
equivalent fine-time representation). Continue denser angular
coverage for late arrivals as a separate problem.

If it does not fall, inspect the RKF Sensor time discretization
and phase interpolation versus a direct phase sweep before
designing any more kernels.

## Interpreting Gemini's new review

- Its recommendation for angular resolution near the long-tail
  separatrix is supported by the 10–100us convergence data.
- The recommended radial coordinate r^2-a^2 linearizes radial
  drift **only for field E~1/r and constant mobility**. In
  actual Garfield++ mobility tables mu(E) can vary, so treat
  this as a feature hypothesis, not an exact coordinate.
- The observed training 95th-percentile radius was ~24um,
  not a demonstrated 99th percentile. Some target starts
  extend to ~64um; a 25um cutoff would require valid
  explicit overflow nodes.
- Fine-pulse 10ns output error does not imply a particular
  VMM3a pulse timing or hit efficiency without independently
  validated shaper/electronics and real U/V weighting fields.

Keep muons 0408 and 0446 untouched until model selection is frozen.
