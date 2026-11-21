# Gate 3: regression stability and independent-muon test (2026-11-21)

**Code committed; real WSL execution NOT YET VERIFIED.** No new
Garfield avalanche, RKF ion tracks, FEM, electronics or U/V strip fields.

## What prompted this

On the actual independent original phi muon 0301 (886727 true avalanche
ions, 4096 weighted K256 ions in the archive; 20 Heed primary rows):
- previous same-reservoir IDW4 whole-plane ion-waveform L2 = 6.2304%;
- monopole 7.0923%, dipole 10.1929%, radial/dipole 6.7664%;
- quadrupole 12.2651%, separatrix gated features 486.77%.
- first 100 ns radial/dipole errors +0.0510/+0.0498 fC;
  native 10 ns time-resolved L2 remains 6.77%.

This does **not** prove angular quadrupole or separatrix physics is absent.
Higher-order regression may be unstable and extrapolating outside the
frozen bank. Current bank has exactly 256 starts and all v>0, drawn from
six nonzero primaries of one original muon.

## New controlled model selection

`study_regularized_moment_replay.py` leaves out *entire original
Heed-primary groups* within muon0's single-ion training bank, trains
each physical-coordinate basis on all other groups, and selects a
ridge penalty strictly by the held-out group mean relative L2.
It then refits on all 256 original site responses and predicts the
already frozen **independent** muon0301 (not involved in fitting or
hyperparameter selection). Zero-penalty pseudoinverse is included as
the original exact comparison. Regularization does **not** fix model
misspecification; it only tests conditioning/generalization of
coefficient fits.

The script writes:
- `cv_grid.csv`: all lambda values, leave-one-primary-out errors;
- `selected_models.csv`: selected penalty and per-feature extrapolation
  proxy on original target sample, including weighted ion fraction
  outside training ranges;
- `event_metrics.csv`: independent original event native 10 ns L2,
  first 100 ns and 100 us window charge errors;
- `first_100ns_waveforms.csv`: early current samples for inspection;
- `metadata.json`: input provenance and scientific limitations.

As before, **only positive avalanche ions** and *ideal complete cathode
planes* are modeled. K256 inverse-inclusion weights preserve primary
ion census but not the full birth-cloud moments exactly. Direct RKF
and the moment prediction use the SAME archived sample.

## Commands (WSL)

Copy and run in the already isolated checkout. Do not enter or modify
`~/na6p/cathode-overnight-2100` for software updates.

```bash
cd "$HOME/na6p/response-rank-pilot-20261121"
git pull --ff-only
source .venv/bin/activate
SCRIPT=test/mwpc-stage-c/cathode-ion-reference/study_regularized_moment_replay.py
python -m py_compile "$SCRIPT"
python "$SCRIPT" --self-test

BANK="$HOME/na6p/cathode-timing-2100/test_runs/mwpc_stage_c/plane_timing_dt10_v1/templates256_dt10_training.root"
ROOT="$HOME/na6p/cathode-overnight-2100/test_runs/mwpc_stage_c/overnight_2100_weighted_k256_v1/muon_0301/ion_births.root"
OLD="$HOME/na6p/cathode-independent-stereo-2100/test_runs/mwpc_stage_c/independent_muon_dt10_v1"
REF="$OLD/muon_0301_direct_dt10_waveforms.csv"
SUMMARY="$OLD/muon_0301_direct_dt10_summary.csv"
IDW="$OLD/muon_0301_fast_dt10_waveforms.csv"
for path in "$BANK" "$ROOT" "$REF" "$SUMMARY" "$IDW"; do
  test -s "$path" || { echo "MISSING $path"; exit 1; }
done
# Inputs remain read only, with the previously verified provenance.
OUT="$PWD/test_runs/mwpc_stage_c/moment_replay_muon0301_ridge_cv_v1"
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python "$SCRIPT" \
  --bank-root "$BANK" --ion-file "$ROOT" \
  --reference-waveforms "$REF" --reference-summary "$SUMMARY" \
  --baseline-waveforms "$IDW" --out-dir "$OUT"

cat "$OUT/selected_models.csv"
cat "$OUT/event_metrics.csv"
cat "$OUT/metadata.json"
```

Send the self-test output, `selected_models.csv`,
`event_metrics.csv`, and if possible the first several lines of
`cv_grid.csv`.

## Decision gate

If ridge reduces extreme extrapolation without sacrificing
held-out-source-primary CV accuracy, evaluate the same frozen fit
across the other two independent muons (0408, 0446).

If it fails, the data do not justify increasing the number of
monomials: switch to stratified/interpolated physics response kernels,
finer coverage of near-separatrix phases, or separately validated
field-line response compression. A small single-ion PCA oracle does
not establish the accuracy of a predictive low-dimensional *avalanche*
model or actual physical U/V strips.
