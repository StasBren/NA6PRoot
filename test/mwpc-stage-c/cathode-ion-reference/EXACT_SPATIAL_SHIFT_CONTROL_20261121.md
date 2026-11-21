# Exact-spatial fractional 10ns replay control (2026-11-21)

## Motivation

After deterministic polar RKF grid refinement (all on development
muon0301, original true K256 ion positions, two ideal full cathode
planes, 2100 V symmetric P3, 100us window):

| grid | native joint L2 | 0–100ns L2 | 10–100us L2 |
|---|---:|---:|---:|
| 6x12 | 5.16445% | 5.19055% | 11.03768% |
| 6x24 | 5.13948% | 5.21280% | 4.98555% |
| 12x12 | 4.30920% | 4.49213% | 11.03780% |
| 12x24 | 4.21196% | 4.46966% | 4.98490% |

Strong angular improvement in the late tail; stronger radial improvement
in the early bins. Previous source-driven IDW4 joint L2 6.23037%.
Both the polar-grid model and existing IDW4 use the **same** current
discrete-time event construction: single-ion response calculated at
birth t=0, shifted by floor(t_birth/10ns) plus fractional linear
interpolation into the same 10ns histogram.

**The remaining early error cannot be attributed entirely to spatial
interpolation before we quantify that time-shift approximation.**

## Exact-spatial test

New target `mwpc_exact_spatial_time_shift_audit`: For **each** real ion
in previously SHA-verified muon0301 K256 archive, run the original
`NA60Plane::Transport::DriftIon` at original `t_birth` and record
the 10ns whole-cathode waveform; run the **same ion position** with
`t_birth=0` and then apply the existing fractional 10ns time shift.
Same mobility file, RKF step and field as original direct reference.

The recomputed at-original-time signal must agree with the frozen
direct-RKF CSV at joint relative L2 < 1e-5, or the tool **refuses**
to interpret the result. It must never alter original ROOT.

Reported: exact-spatial shifted-vs-direct joint/native 10ns L2,
four time-window L2 metrics, first100ns and full100us charge errors,
late-ion counts, complete intermediate currents for independent audit.

**Scientific distinction:** The test simulates ion RKF twice, but never
runs a Geant4/TrackHeed/AvalancheMicroscopic avalanche. The t=0
calculation has exact spatial coordinates and is a *model-free control*
of the approximate time shift. It is not yet a deployable fast scheme.

## WSL

```bash
cd "$HOME/na6p/response-rank-pilot-20261121"
git pull --ff-only
source .venv/bin/activate
export ROOT_INSTALL="$HOME/na6p/install/root"
export GARFIELD_INSTALL="$HOME/na6p/install/garfieldpp"
source "$ROOT_INSTALL/bin/thisroot.sh"
cmake -S test/mwpc-avalanche -B test/mwpc-avalanche/build-polar \
  -DROOT_DIR="$ROOT_INSTALL/cmake" \
  -DCMAKE_PREFIX_PATH="$ROOT_INSTALL;$GARFIELD_INSTALL"
cmake --build test/mwpc-avalanche/build-polar \
  --target mwpc_exact_spatial_time_shift_audit -j2

BIN="$PWD/test/mwpc-avalanche/build-polar/mwpc_exact_spatial_time_shift_audit"
MOBILITY="$HOME/na6p/src/garfieldpp/Data/IonMobility_Ar+_Ar.txt"
DATA="$HOME/na6p/cathode-overnight-2100/test_runs/mwpc_stage_c/overnight_2100_weighted_k256_v1/muon_0301/ion_births.root"
BASE="$HOME/na6p/cathode-independent-stereo-2100/test_runs/mwpc_stage_c/independent_muon_dt10_v1"
OUT="$PWD/test_runs/mwpc_stage_c/polar_grid_2100_v1/exact_spatial_time_shift_muon0301_v1"
"$BIN" --ion-file "$DATA" --mobility "$MOBILITY" \
  --reference-waveforms "$BASE/muon_0301_direct_dt10_waveforms.csv" \
  --reference-summary "$BASE/muon_0301_direct_dt10_summary.csv" \
  --out-prefix "$OUT"
cat "${OUT}_summary.csv"
```

Requires original source integrity previously verified by SHA manifest.
New output prefixes must not exist; never overwrite. Send last
progress and summary if execution fails.

**Interpretation:** If exact-spatial time-shift L2 is already near
4%, further radial/angular grid refinement cannot solve the main
issue under the current 10ns shift rule. Improve event timing
translation/signal integration first. If exact-spatial timing error
is small, prioritize spatial interpolation refinement, with extra
near-wire radial nodes and near-separatrix angular nodes.

Two other independent muons 0408 and 0446 remain sealed.
