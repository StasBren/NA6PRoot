# Independent-muon physical-moment replay — read-only Gate 2 (2026-11-21)

**Status:** script committed; local WSL self-test and real-muon results PENDING.
This is NOT a discovery of an exact minimum avalanche descriptor and does
NOT compute actual U/V stereo strips.

## Why

Gate 1 found that the 256 original single-ion starts have 99.28% of the
uncentered response-bank squared L2 energy in the first singular mode,
but the held-out-ion SVD oracle error remains 0.865% median with 32 modes.
Physical birth-coordinate prediction: radial/dipole ~3.925% median,
quadrupole ~3.673%, existing 4-neighbor IDW ~2.417%. Gate 1 bank is
**256 positive-v starts only**. The SVD oracle uses true test waveforms:
it is NOT a fast model.

Now test whether **summed birth-coordinate moments** can reconstruct
a complete *independent* muon's two-plane time-dependent signal, starting
from the same true weighted K256 ion birth samples as direct RKF.

## Method

- Train physical single-ion response operator on all 256 original
  10 ns positive-v templates from source muon0.
- Read a separately generated SHA-audited whole-muon weighted K256 ROOT
  archive. Preserve exact per-primary populations and sampling weights.
- For each ion, reflect negative-v starts into the positive hemisphere
  and swap complete-cathode +/- output channels. Reflection is exact
  only for the **symmetric 3+3 mm nine-wire analytic whole-plane model**,
  not the Prototype-3 hardware or physical U/V strips.
- Deposit weighted physical-coordinate features in *birth-time bins*,
  split fractional 10 ns shifts exactly as the previous C++ kernel.
- Convolve grouped moment time streams with the fitted response basis
  by FFT, using no per-ion RKF; sum all primaries and both hemispheres.
- Compare 10000-bin predicted ion currents to an already-recorded direct
  ion RKF waveform for the exact same source archive. Optionally
  include the *previous* independent-muon IDW4 output as a baseline.
- Write one row per primary containing sampled angular/radial
  descriptors, plus the event metrics, waveforms, provenance.

This is a **given-ion-birth compression/replay** test. Generating
avalanche gains, clouds and moments directly from new Heed electrons
is a separate, later problem. Also, the training bank has 256
sampled sites only; any bulk extrapolation or finite nine-wire edge
effects can harm independent-muon performance.

## WSL commands: test one independent muon first

In your pre-existing isolated worktree:
```bash
cd "$HOME/na6p/response-rank-pilot-20261121"
git pull --ff-only
source .venv/bin/activate

SCRIPT=test/mwpc-stage-c/cathode-ion-reference/study_event_moment_replay.py
python -m py_compile "$SCRIPT"
python "$SCRIPT" --self-test

BANK="$HOME/na6p/cathode-timing-2100/test_runs/mwpc_stage_c/plane_timing_dt10_v1/templates256_dt10_training.root"
DATA="$HOME/na6p/cathode-overnight-2100/test_runs/mwpc_stage_c/overnight_2100_weighted_k256_v1/muon_0301/ion_births.root"
test -s "$BANK" && test -s "$DATA" || {
  echo "MISSING INPUT BANK OR SHA-AUDITED MUON; STOP"
  exit 1
}
# The archived original K256 input has complete.txt + SHA256SUMS.
D="$(dirname "$DATA")"
(cd "$D" && test -s complete.txt && sha256sum -c SHA256SUMS) || {
  echo "SOURCE MUON NOT VERIFIED: STOP"
  exit 1
}

mapfile -t MATCHES < <(find "$HOME/na6p" -type f   -name 'muon_0301_direct_dt10_waveforms.csv' -print)
if (( ${#MATCHES[@]} != 1 )); then
  printf 'Expected one existing direct RKF result; found %s:\n' "${#MATCHES[@]}"
  printf '%s\n' "${MATCHES[@]}"
  exit 1
fi
REF="${MATCHES[0]}"
STEM="${REF%_waveforms.csv}"
SUMMARY="${STEM}_summary.csv"
BASELINE="${STEM%_direct_dt10}_fast_dt10_waveforms.csv"
test -s "$SUMMARY" || { echo "MISSING $SUMMARY"; exit 1; }

OUT="$PWD/test_runs/mwpc_stage_c/moment_replay_muon0301_v1"
CMD=(python "$SCRIPT" --bank-root "$BANK" --ion-file "$DATA"
  --reference-waveforms "$REF" --reference-summary "$SUMMARY"
  --out-dir "$OUT")
if [[ -s "$BASELINE" ]]; then
  CMD+=(--baseline-waveforms "$BASELINE")
else
  echo "Existing IDW baseline not found. Will run moments-only."
fi

OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 "${CMD[@]}"
cat "$OUT/event_metrics.csv"
cat "$OUT/metadata.json"
cat "$OUT/per_primary_sampled_descriptors.csv"
```

Only after muon0301 is checked should we replay independent
0408 (J/psi, station 3) and 0446 (omega, station 2). Use a NEW
output directory for each event and the matching independently
recorded direct-RKF reference. If there are multiple matches for
the reference CSV, STOP and select the correct known existing
analysis directory manually; do not mix independent event data.

## Integrity / limitations

- Script refuses to overwrite output directory and checks original
  per-primary ROOT sample counts and exact weighted ion totals.
- Direct summary must identify the exact crossing, station, detector,
  true ion count, stored macroion count, 10 ns bins and 100 us window.
- For existing IDW baseline, the full direct waveform columns must
  match the same reference CSV before it is used.
- No new avalanche simulation, Garfield ion drift, or FEM.
- Comparisons use ONLY avalanche positive ions. No electron prompt
  current, original Heed primary ions, ion-species mixture mobility,
  VMM3a electronics, or actual U/V geometry.
- Some RKF ion trajectories exceed 100 us. Compare **windowed** signals;
  this is not the complete Ramo endpoint charge.
- Both K256 sampling error and physical-input systematics remain
  separate from the moment-response approximation error.
- The new per-primary descriptor CSV is computed from weighted sampled
  ion positions. It does not contain an exact microscopic cloud
  descriptor when K256 is a subsample.
- Passing a synthetic test proves bookkeeping only, not real accuracy.
- **Stop on failed preflight or mismatched metadata.** Do not submit
  new simulation jobs or change the 1000-muon recovery worktree.

Return the terminal self-test, event_metrics.csv, metadata.json, and
per_primary_sampled_descriptors.csv. Do not claim success or failure
of nine-parameter *avalanche generation* from these results.
