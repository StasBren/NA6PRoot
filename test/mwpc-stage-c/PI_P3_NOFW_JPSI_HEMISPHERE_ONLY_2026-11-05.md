# p3_noFW: J/psi 1000 + 10000 hemisphere-only study

Scope: ONLY exact full-Garfield ion-birth hemisphere asymmetry
A_hemi = (N_ion_up - N_ion_down)/(N_ion_up + N_ion_down).
Zero denominator is undefined. No 90-degree narrow-sector quantities.
No cathode induced charge, no ion drift, no electronics.

Data sources: original 1000 fully completed and separate extra 10000 in
test_runs/mwpc_stage_c/pi_p3_noFW_jpsi1000_native_v1/
test_runs/mwpc_stage_c/pi_p3_noFW_jpsi_extra10000_native_v1/

Only SHA-sealed completed 20-parent chunks are included. Any actively
running, incomplete or unsealed chunk is excluded. Native source input
SHA and every accessed CSV against its SHA-sealed receipt are checked.
All 1000 and any new parents are checked for exact unique Geant4
parent IDs and strict disjointness. Parent, not ion, is the statistical unit.

## Launch

~~~bash
cd ~/na6p/src/NA6PRoot-mwpc
git switch feature/mwpc-stage-c
git pull --ff-only
PY=test_runs/mwpc_stage_c/.venv_pi_slides/bin/python
SCRIPT=test/mwpc-stage-c/pi_p3_noFW_jpsi_hemisphere_analysis_v1.py
"$PY" -m py_compile "$SCRIPT"
"$PY" "$SCRIPT" selftest

# Now: exactly 1000 original full-Garfield muons
"$PY" -u "$SCRIPT" analyze --old-only

# During extra-10000 production: 1000 + SHA-sealed completed new parents
"$PY" -u "$SCRIPT" analyze

# Once all 10000 additional original parents are complete:
"$PY" -u "$SCRIPT" analyze --require-complete
~~~

Each execution creates a NEW immutable snapshot:
test_runs/mwpc_stage_c/pi_p3_noFW_jpsi_hemisphere_v1/
snapshot_XXXXX_old1000_newYYYYY/
An identical snapshot is never overwritten; do not issue identical
commands twice without additional completed chunks.

## Plot and CSV outputs

01 exact event-level A_hemi histogram.
02 A_hemi and |A_hemi| versus true signed/absolute Geant4 incidence angles.
03 A_hemi and |A_hemi| versus true entry u phase around nearest 3mm anode.
04 |A_hemi| vs Heed number and A_hemi vs Heed upper/lower side imbalance.
05 observed A^2 vs independent-flipped avalanche side-sign null prediction.
06 2D phase-vs-incidence map of mean |A_hemi|; low-stat cells masked.
exact_hemisphere_per_muon.csv: exact counts, parent IDs, Geant4 angles,
entry phase, TrackHeed side imbalance and event-level null diagnostics.
summary.json: physical counts, bootstrap confidence intervals, all SHA.

Important physical distinction: full Garfield ion-birth census is used
for A_hemi, never the 128-position-per-avalanche reservoir.
Anode phase uses the EXACT C++ station-dependent parity rule.
Independent-sign-flip calculation retains each avalanche's full
ion count and intrinsic side-imbalance magnitude, independently randomizes
its orientation; it is diagnostic ONLY, not evidence that orientations
actually are random nor a replacement for Garfield physical simulation.

## First 1000 results from user-uploaded exact per-event CSV

999 defined event asymmetries. Mean +0.014319, SD 0.370703.
Mean |A| 0.299817; fraction |A|>0.5 is 0.19319.
Pearson corr(|A|, N_Heed) = +0.20807.
Pearson corr(|A|, |theta_u|) = +0.03027.
No clear broad dependence of magnitude on absolute incidence angle.
Neither this nor the sign-flip null identifies a physical cause.

Physics benchmark and original 1000 campaign are unmodified.
