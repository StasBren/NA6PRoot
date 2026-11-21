# Gate 4: deterministic central-wire polar RKF operator (2026-11-21)

**Status:** GitHub implementation committed; WSL C++ compilation and runtime
not yet verified. An experimental numerical gate, NOT detector performance
evidence, not a declaration of sub-1% accuracy.

## Motivation: audited source vs target coverage

Previously frozen single-ion bank: 256 starts, six Heed-primary groups, all
positive v, angular support ~31.60°–151.11° relative to central wire.
Development muon 0301: 886727 real ion births represented by 4096
inverse-inclusion weighted samples from 20 primaries. Its folded
positive-hemisphere angle extends ~88.77°–179.995°.
Target weighted ions outside original marginal intrinsic r/angle ranges:
**71.001%**. Source-bank leave-one-out 95%-nearest-neighbor threshold
exceeded by **70.205%** of target weighted ions in intrinsic (ln r, θ).
About **27.5%** source-weighted ions have |sin θ| < .15, compared with
zero in the source bank.

All 256 training ion starts and 4096 target sampled ions have central-wire
nearest index zero, ruling out multiple-wire translations as the cause
for this particular muon. But the finite nine-wire patch does not justify
reusing central-wire templates for noncentral wires.

Source and target radial medians are 16.91 and 16.69 um, respectively;
angular support mismatch is the main immediately observed problem.

Separatrix Gaussian feature in code is bounded 0–1 and ridge features
are RMS normalized. Do not attribute prior 241% failure to a
divergent 1/sinθ feature or absent scaling.

## Controlled architecture

### Stage 4A: small **uncompressed** polar operator

A *dedicated C++ target* runs only deterministic Garfield++
DriftLineRKF single positive ions at fixed (r, θ), using the
original Plane2100Response.hh geometry, mobility, fixed w=0, RKF
step size, 10ns bins and 100us current readout.
It does not call AvalancheMicroscopic or TrackHeed, nor alter any
archived muon.

- 6 log-spaced **r** nodes between (15um + 0.00001um) and 80um.
- 12 cosine-spaced positive-v **θ** nodes from 0.001° to 179.999°,
  clustered near horizontal v=0 but never exactly on that separatrix.
- 72 deterministic single-ion starts total; deliberately small.
- Explicit RKF late-arrival recording, so a 100us recorded window
  is not misrepresented as complete endpoint charge.
- ROOT `polar_bank` entries contain each original single-ion
  current impulse, both whole planes, exactly 10000 bins per plane.
- Fail on failed transport or invalid geometry. No silent truncation
  or pretend continuation of incomplete grid.

The **Python replay** reads the complete bank and the already
SHA-verified K256 original muon0301 ion-birth ROOT, finds central-wire
(r, θ) for each macroion, maps negative-v to positive-v with
electrode exchange, linearly interpolates **uncompressed waveforms**
on the (ln r, θ) grid, weights by inverse-inclusion, and shifts using
the same fractional 10ns bin rule as our original IDW4 control.
Time convolution uses FFT, but no SVD or regression.

It refuses extrapolation, noncentral-wire targets and output
overwriting. It compares to the already recorded direct-RKF
100us two-plane signal for exactly the same weighted source ions.

### Stage 4B after analyzing pilot

Only if 4A builds, transports, and passes its synthetic interpolation/
timing smoke test: attempt a denser polar-grid convergence scan
(e.g. 12×24 and/or 16×32). Grid generation must be completed one
controlled job at a time on a suitable machine. First compare
**uncompressed** grid interpolation; add optional SVD rank tests
only *after* spatial interpolation is understood.

Do NOT claim that a 16×32 grid is enough or that this model is faster
than the existing IDW4 until both have timed runs on the same machine.

## WSL commands: one small local pilot, not a PBS campaign

```bash
cd "$HOME/na6p/response-rank-pilot-20261121"
git pull --ff-only
source .venv/bin/activate
export GARFIELD_INSTALL="$HOME/na6p/install/garfieldpp"
export LD_LIBRARY_PATH="$GARFIELD_INSTALL/lib:$HOME/na6p/install/root/lib:${LD_LIBRARY_PATH:-}"

PY=test/mwpc-stage-c/cathode-ion-reference/replay_deterministic_polar_bank.py
python -m py_compile "$PY"
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python "$PY" --self-test

# Compile ONLY the independent polar RKF bank tool.
cmake -S test/mwpc-avalanche -B test/mwpc-avalanche/build-polar   -DCMAKE_PREFIX_PATH="$GARFIELD_INSTALL"
cmake --build test/mwpc-avalanche/build-polar   --target mwpc_deterministic_polar_plane_bank -j2

MOBILITY="$HOME/na6p/src/garfieldpp/Data/IonMobility_Ar+_Ar.txt"
BIN="$PWD/test/mwpc-avalanche/build-polar/mwpc_deterministic_polar_plane_bank"
PILOT="$PWD/test_runs/mwpc_stage_c/polar_grid_2100_v1"
mkdir -p "$PILOT"
test -s "$MOBILITY" && test -x "$BIN" || {
  echo "MOBILITY OR BINARY MISSING"; exit 1;
}

"$BIN" --mobility "$MOBILITY"   --out-root "$PILOT/grid6x12_dt10.root"   --n-radius 6 --n-angle 12   --clearance-min-um 0.00001 --radius-max-um 80   --theta-min-deg 0.001   --dt-ns 10 --window-us 100 --max-step-mm 0.05
```

If the full grid succeeds, compare on **0301 only**:

```bash
BANK="$PILOT/grid6x12_dt10.root"
DATA="$HOME/na6p/cathode-overnight-2100/test_runs/mwpc_stage_c/overnight_2100_weighted_k256_v1/muon_0301/ion_births.root"
BASE="$HOME/na6p/cathode-independent-stereo-2100/test_runs/mwpc_stage_c/independent_muon_dt10_v1"
REF="$BASE/muon_0301_direct_dt10_waveforms.csv"
SUMMARY="$BASE/muon_0301_direct_dt10_summary.csv"
OLD="$BASE/muon_0301_fast_dt10_waveforms.csv"
for P in "$BANK" "$DATA" "$REF" "$SUMMARY" "$OLD";do
  test -s "$P" || { echo "MISSING $P"; exit 1; }
done
# Verify the original source archive remains SHA-sealed.
(cd "$(dirname "$DATA")" && test -s complete.txt &&
  sha256sum -c SHA256SUMS) || {
  echo "SOURCE ARCHIVE INTEGRITY FAILED"; exit 1;
}
OUT="$PILOT/replay_muon0301_uncompressed"
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python "$PY"   --grid-root "$BANK" --ion-file "$DATA"   --reference-waveforms "$REF" --reference-summary "$SUMMARY"   --baseline-waveforms "$OLD" --out-dir "$OUT"
cat "$OUT/event_metrics.csv"
cat "$OUT/metadata.json"
```

Run grid generator only once per fresh output path: if it fails, retain the
partial diagnostic and do not overwrite it. Report the last GRID line and
the failure, rather than trying dozens of arbitrary angles. Close-to-v=0
RKF can be expensive or fail stringent endpoint checks.

## Evaluation criteria

On the same 0301 K256 reservoir:
- Native **10ns** joint and per-plane current L2, against direct RKF;
- Existing IDW4 (6.2304%) **unchanged baseline**;
- First-100ns induced charge errors, 100us recorded charge errors;
- Fraction of macroions outside grid: must be zero for interpolation;
- Number of late grid ion arrivals and max transit time;
- Actual grid RKF and replay wall time separately.

Small 6×12 grid is a smoke/coverage test, not an acceptance test.
No hard-coded <1% success threshold is justified before actual strip
readout error budget is defined. For real deployment, reconfirm
independent muons 0408/0446 **only after architecture is frozen**.
