# P3 full-cloud Avalanche LUT v2 — executable first gate (2026-11-04)

## Scope

This is an **independent, append-only research prototype** on branch
\`feature/mwpc-stage-c\`; it does not replace the validated \`guarded_A20\`
count/sector LUT, modify old SHA-sealed physical runs, or claim cathode-ready
detector response.

**Physical source already completed:** 120/120 SHA-frozen native Garfield
chunks with 960 independent original D500 starting electrons (480 for each
of \`p3_noFW\`, \`p3_FW75\`). The exact census is **2,981,533 ionization
BIRTH coordinates** (\`u_mm,v_mm,w_mm,t_ns\`): noFW 284,272, FW75 2,697,261.
Every new archive record is the complete cloud from ONE physical microscopic
avalanche, not a random 128/512/4096-point subsample and not an ion trajectory.

## What was implemented

\`pi_p3_joint_fullcloud_bank_v1.py\` provides four commands.

- \`selftest\`: local toy tests of complete-cloud replay, w displacement, absolute
  TrackHeed time offset, true zero-ion record, physical metal veto, abstention
  when no nearby data, deterministic holdout, proper-score arithmetic.
- \`build\`: read-only SHA audits of all native input chunks and old analysis.
  Independently streams their full birth CSVs, verifies all exact C++ ion
  counts/moments/sectors, then writes a new **binary lossless float64 bank**
  using \`numpy\` memory-mappable \`.npy\` files, one per geometry, plus
  \`index.json\` with start, cell, sample, zero vs nonzero, offset and length,
  and a SHA256 \`manifest.json\`. Final results are written in a distinct
  \`test_runs/mwpc_stage_c/pi_p3_joint_fullcloud_bank_v1/\` folder.
  Build creates a separate \`.INCOMPLETE\` staging folder first, refuses
  any prior output, and renames only on complete success. **Do not delete
  an incomplete directory without diagnosing the failure.** Original
  source CSVs, analyses, and physics checkpoints are NEVER changed.
- \`audit\`: verifies all source/checkpoint SHA receipts, stored binary
  and index SHA, exact offset/census continuity and all 960 identities.
- \`evaluate\`: **spatially disjoint** SHA-deterministic ~20% of D500 cells
  held out PER geometry; two safeguarded sampling models are tested with
  *only the remaining cells*, using nearest admissible cell versus
  distance-weighted guarded neighbors. **Zero prediction coverage is
  reported as abstention, not a fabricated zero avalanche or score.**
  Conditional ion-count distribution scored by exact discrete CRPS,
  zero-ion events by Brier score; whole-cloud seven-dimensional joint
  features scored by a multivariate energy score: log(1+N), mean and SD
  birth radius, mean and q90 birth time, mean longitudinal w, top sector
  fraction. Feature scales learned on TRAINING cells only. Scores are
  averaged over independent D500 cells; paired comparison includes
  exploratory 2,000 cell-bootstrap replicates. This is NOT a full 4D
  point-process test, and uses a sparse non-flux source.

The replay API is \`Bank.sample(u0,v0,w0,t0,g,rng,mode,train)\`: returns
\`n_ion_births\` and \`births_uvwt_mm_ns\` together from the **same**
avalanche. It shifts only \`w\` by the requested starting-wire coordinate
and adds **TrackHeed t0** to stored full-microscopic birth times. It does
NOT add D500 drift time (which would double count primary-electron drift);
it does NOT shift the cloud in transverse \`u,v\`, because doing so would
distort the actual wire electric field response. This discrete empirical
replay approximates the new source position with the chosen archive start.

The guards forbid mixing the two geometries, v hemispheres, near-wire
regimes, opposite sides of the anode near the wire, and opposite sides
of field wires in the \`near_field\` regime. Candidates live within a
5x5 cell index neighborhood. Nearest-cell and guarded-neighbor baselines
in this *new bank* are NOT the same trained model as historical
\`guarded_A20\`; **no model replacement or superiority is inferred**.
For a position inside conductor metal the API returns a physical
zero-ion cloud (\`mode=metal_veto\`); where no valid source cloud exists
it returns \`None\`, deliberately avoiding unsupported extrapolation.

## Execute in user's WSL

~~~bash
cd ~/na6p/src/NA6PRoot-mwpc
git switch feature/mwpc-stage-c
git pull --ff-only
PY=test_runs/mwpc_stage_c/.venv_pi_slides/bin/python

"$PY" -m py_compile test/mwpc-stage-c/pi_p3_joint_fullcloud_bank_v1.py
"$PY" test/mwpc-stage-c/pi_p3_joint_fullcloud_bank_v1.py selftest

set -o pipefail
"$PY" -u test/mwpc-stage-c/pi_p3_joint_fullcloud_bank_v1.py build \
  | tee "$HOME/na6p/p3_joint_fullcloud_bank_build_v1.log"

"$PY" -u test/mwpc-stage-c/pi_p3_joint_fullcloud_bank_v1.py audit \
  | tee "$HOME/na6p/p3_joint_fullcloud_bank_audit_v1.log"

"$PY" -u test/mwpc-stage-c/pi_p3_joint_fullcloud_bank_v1.py evaluate \
  | tee "$HOME/na6p/p3_joint_fullcloud_bank_evaluate_v1.log"
~~~

The source code and its toy selftest have been locally syntax-checked
and toy-tested; **actual WSL \`build/audit/evaluate\` are not yet observed**.
Do not interpret these as completed or invent their results.

## Decision after the first WSL report

1. Check coverage: how often a valid archive source exists within the
   physical guards. Sparse training coverage is *not* evidence that
   a bank is accurate; many abstentions call for better source sampling.
2. Compare nearest versus guarded with held-out whole-cell paired
   CRPS/Brier/joint-feature energy scores and their confidence intervals.
   Do not tune model weights against these held-out cells and then report
   that same sample as blind generalization.
3. Before merging into the existing fast native-muon route, implement
   **proper training-data growth** (20/50/100 physical independent
   avalanches per selected cell) and compare with old protected A20 on a
   new independent Geant4→Heed→Garfield parent sample. Count accuracy
   and space-time shape both matter; retaining each complete cloud
   preserves correlations but does not prove conditional accuracy.
4. Ion DRIFT, charge induction by Shockley–Ramo and U/V cathode strips
   remain separate follow-on steps. Preserve the native complete clouds
   as ground truth regardless of future interpolation/compression choices.

Do not rebuild or modify the SHA-frozen native executable or the
original full-cloud source; this is a Python read-only derivative.
