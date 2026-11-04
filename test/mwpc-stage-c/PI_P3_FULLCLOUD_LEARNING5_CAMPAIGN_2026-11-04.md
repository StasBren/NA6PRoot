# P3 five-complete-clouds-per-cell learning curve: approved source, staged physical run

## Rationale (grounded in actual WSL results)

The original 960 complete Garfield avalanche clouds and their SHA bank passed
their full-ion census audit (noFW: 480 avalanches/284,272 ion births;
FW75: 480 avalanches/2,697,261 ion births).

On the 900-cell D500 original-start grid, geometric guarded source **availability**
is 8919/9000 (99.10%) in noFW and 8836/9000 (98.18%) in FW75 from ten
held-as-source D500 start positions per cell. This says *a candidate archive
cloud exists*; it does not certify exact Garfield conditional physics,
full four-dimensional birth distributions or detector signal precision.
The locally supported fractions were 34.56% / 28.69%, respectively.
Guarded neighbors improved external spatial-box and sparse within-box
held-outs, motivating improved per-cell sample depth rather than an
unvalidated new interpolation radius.

We therefore stage an independent physical sampling campaign with
**five genuinely new full microscopic Garfield avalanches per D500 cell
per geometry**; 900 cells x 5 starts x 2 p3 geometries =
**9,000 independent single-electron Garfield runs**, exact
`u_mm,v_mm,w_mm,t_ns` birth locations and zero-ion outcomes.
This is the *first learning-curve point*, not acceptance of a full
conditional model, and does not replace the validated count/sector
Guarded A20.

## Physical source/provenance constraints

- Inputs: exact original SHA-protected D500 random initial electron sites,
  **never** previously calculated drift endpoints or a Heed muon replay.
- Original A20/A100 used D500 permutation passes 0..99.
- Completed 960 full-cloud extended campaign used passes 100..139.
- New learning5 selects five **gas-only** positions per cell from reserved
  passes **140..219**. Every original 48-cloud pilot start is also
  explicitly excluded regardless of pass number. New selection forbids
  all original 960 starts.
- Passes 460..499 are set aside for *later*, independently simulated
  single-electron validation, not touched by training. Note the old
  48-cloud pilot sampled from a larger >=100 pass range: remove any
  overlapping old pilot identities from the eventual reserved validation,
  rather than assuming every pass460..499 site is untouched.
- Both P3 configurations are simulated separately. The source is uniform
  by spatial cell, **NOT physical incoming-muon flux**.
- Reuse the *already compiled and SHA-pinned* C++ target
  `mwpc_phase_b_joint_cloud_lut_v1 --lut-sidecar-full-births`; do NOT
  rebuild or change its gas, field, seed, physics code or original binary.
- New **only** output root:
  `test_runs/mwpc_stage_c/pi_p3_fullcloud_learning5_v1/`.
  Earlier complete-cloud checkpoints, A20, A100, D500, full-cloud binary
  bank and all external validation results remain immutable.
- Every new 8-start chunk records original exact coordinates,
  source trial records, all ion births, summary moments and endpoint
  reservoirs, physics log, hashes and final `complete.json` SHA seal.
  `status` tolerates one most-recent UNSEALED active writer directory
  but never treats it as finished or silently deletes it. No automatic
  retry of an interrupted partially written chunk.

## FIRST EXECUTION GATE (safe, one native 8-start block only)

Run in WSL with the original Garfield libraries and activated Python
virtual environment. First prove the source and native parser, then
run **just one** eight-start chunk. No background production yet.

~~~bash
cd ~/na6p/src/NA6PRoot-mwpc
git switch feature/mwpc-stage-c
git pull --ff-only
PY=test_runs/mwpc_stage_c/.venv_pi_slides/bin/python
SCRIPT=test/mwpc-stage-c/pi_p3_fullcloud_learning5_campaign_v1.py

"$PY" -m py_compile "$SCRIPT"
"$PY" "$SCRIPT" selftest
df -h "$HOME/na6p"
"$PY" -u "$SCRIPT" prepare
"$PY" -u "$SCRIPT" preflight

set -o pipefail
"$PY" -u "$SCRIPT" run --max-chunks 1 --hours 2 \
  | tee "$HOME/na6p/p3_fullcloud_learning5_first8.log"

"$PY" "$SCRIPT" status
~~~

Expected if actual native Garfield physics succeeds:

- `PASS FROZEN new full-cloud P3 learning5 inputs: 9000...`
- `PASS real C++ parser 8 unique unused D500 starts and frozen binary`
- `PASS SHA FULL BIRTHS 0 ions ...`
- `VERIFIED CHUNKS 1 / 1125` and `8 / 9000` true starts.

**These native commands have NOT been executed by the assistant.** A
Python-only synthetic `selftest` is not evidence that any of the 9000
physical Garfield runs completed. Share output of the first block and
`status` before authorizing a longer native run. Check free space;
full-birth CSV files can be much larger than the existing binary bank.

## AFTER first-chunk gate PASS (NOT yet an execution instruction)

Continue remaining 1,124 eight-start blocks with a single
background writer and a soft clock cutoff between chunks; reuse the
same prepared input and C++ binary. No second writer while the old one
is active. Recheck status SHA seals after any interruption. If a chunk
remains unsealed after writer exit, inspect its `physics.log` and raw
data; do not blindly delete or overwrite it.

After completion, compare conditional full-cloud sampling at
1/3/5 independent training events per cell with held-out Garfield
events, per geometry and region, including count CRPS, zero-Brier,
spatial/time joint scores, uncertainty over independent source cells
and parent muon-level count/sector baseline. Full 4D point-process
and U/V induced-charge validation are later gates, not performed here.

## Execution status

GitHub script and source protocol are staged.
**Physical 9000-start campaign: not yet run.**
