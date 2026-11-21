# Gate 6: exact-position cumulative Ramo charge from RKF path points

**Status:** Source committed but C++ target not yet compiled or run in WSL.

## Key numerical facts from development muon 0301

The original K256 weighted archive contains 4,096 saved real ion
positions from 20 Heed primaries, representing 886,727 original births.
The existing direct 10ns whole-plane RKF waveform (100us window)
was independently reproduced with relative L2 = 1.154866e-16.

| Reconstruction (same muon 0301) | Whole L2 | Early 0–100ns L2 | 10–100us L2 |
|---|---:|---:|---:|
| Original IDW4 bank | 6.23037% | 5.37754% | 45.4581% |
| 12x24 polar grid, old 10ns timing | 4.21196% | 4.46966% | 4.98490% |
| 12x24x4 phase-aware polar grid | 1.35451% | 1.28163% | 4.98718% |
| Exact birth positions, old shift | 3.99762% | 4.29363% | 0.01018% |
| Exact birth positions, phase4 | 0.92494% | 0.99382% | 0.00173% |

Phase4 polar versus original IDW4 improves full-event L2 by ~78.3%,
but 100us window charge errors remain -0.21776/-0.28483 fC,
almost exactly unchanged by four-phase interpolation. We now know
the phase4 polar error decomposition: relative spatial residual norm
0.85160%, relative timing residual norm 0.92494%, and nonzero
normalized squared cross-term 2.5397e-5. Quadrature subtraction
of percentages is not valid.

## Physics proposal

In static deterministic ion drift, the cumulative induced charge on
readout k at elapsed ion-drift time tau is

    Q_k(tau) = e [phi_k(X(tau))-phi_k(X(0))].

For an actual ion birth at t_birth and 10ns readout bin [n*10,(n+1)*10],
the mean current is obtained by exact **endpoint differencing**:

    I_k[n] = [Q_k(max(0, (n+1)*10 - t_birth))
             - Q_k(max(0,n*10 - t_birth))] / 10ns.

Evaluate phi_k using the exact original 2100V 9-wire analytic Garfield
whole-plane weighting potential. Ion positions X(tau) come from the
adaptive Garfield DriftLineRKF trajectory at its saved points.
We currently linearly interpolate the cumulative **charge** between
successive RKF point times. That approximation can still fail if
Garfield's adaptive trajectory steps are too coarse near the anode
or cathode. The new test quantitatively measures this uncertainty;
it does NOT assume the claimed <0.1% or <0.05% accuracy beforehand.

## What the new C++ target does

`mwpc_exact_spatial_cumulative_charge_audit`: on exactly the original
muon0301 weighted K256 archive, replay every single ion directly at
its original birth time (checking full-event numerical identity with
the original frozen direct RKF CSV), and independently drift the same
ion from time t=0. Record the adaptive RKF drift-line points
`GetDriftLinePoint`; convert each point into q*delta weighting potential;
interpolate this Q(tau) at exact birth-shifted 10ns bin endpoints.

Safety and integrity:

- **No** new Heed/AvalancheMicroscopic simulation or input replacement.
- Same ion positions, weights, field, provisional mobility, RKF step
  and original two whole-cathode electrodes as the frozen reference.
- Fail on mismatched original direct waveform (L2 must be <1e-5).
- Explicit original positive-ion cathode endpoint charge.
- Record number of RKF path points, late arrivals, and worst first
  200ns gap between adaptive path point times, so a too-coarse path
  is visible instead of being silently called exact.
- New files use a non-overwriting output prefix, preserving frozen data.

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
  --target mwpc_exact_spatial_cumulative_charge_audit -j2

BIN="$PWD/test/mwpc-avalanche/build-polar/mwpc_exact_spatial_cumulative_charge_audit"
MOBILITY="$HOME/na6p/src/garfieldpp/Data/IonMobility_Ar+_Ar.txt"
DATA="$HOME/na6p/cathode-overnight-2100/test_runs/mwpc_stage_c/overnight_2100_weighted_k256_v1/muon_0301/ion_births.root"
BASE="$HOME/na6p/cathode-independent-stereo-2100/test_runs/mwpc_stage_c/independent_muon_dt10_v1"
OUT="$PWD/test_runs/mwpc_stage_c/polar_grid_2100_v1/exact_spatial_cumulative_Q_muon0301_v1"

if [[ -e "${OUT}_comparison.csv" || -e "${OUT}_summary.csv" ]]; then
  echo 'Audit output exists, refusing overwrite'
else
  "$BIN" --ion-file "$DATA" --mobility "$MOBILITY" \
    --reference-waveforms "$BASE/muon_0301_direct_dt10_waveforms.csv" \
    --reference-summary "$BASE/muon_0301_direct_dt10_summary.csv" \
    --out-prefix "$OUT"
  cat "${OUT}_summary.csv"
fi
```

## Interpreting the result

- If exact-spatial cumulative Q(t) error is substantially below
  0.925%, use the same trajectories and weighting potentials to
  design a compact time-adaptive spatial bank (next separate task).
- If the error is high, compare the recorded near-wire RKF path gaps
  to the first 100ns response timescale. Improve interpolation along
  the trajectory, not the polar grid, first.
- Time/space errors are not independent; any future comparison must
  preserve the full residual cross term.
- Do NOT extrapolate 0301 results to new muons or to physical U/V
  stereo strips. Leave locked muons 0408/0446 untouched until
  the full method is frozen. The current model still takes real
  ion births as input; it has not compressed the avalanche itself.

## Gemini review status

Gemini's cumulative-Q suggestion is mathematically appropriate under
the deterministic no-diffusion conditions, but its numerical
accuracy and 4x memory improvement are hypotheses until verified.
Its claim that subtraction of squared aggregate errors isolates
spatial effects omits the measured cross term. Near a trajectory
separatrix, smoothness in (r,theta) or in a rescaled arrival time
must also be established rather than assumed.
