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


## 2026-11-04: User first WSL attempt stopped in Python source selection; FIXED

The original script called `old.sample_for(g,box,pass_index)` from
`pi_p3_extended_full_cloud_campaign_v1.py` with `pass_index` in
`140..219`. That OLD helper is correctly **restricted to 100..139**.
Thus `selftest` and `prepare` failed with `ValueError: Invalid
reserved D500 trial pass`. Later `preflight/run/status` raised
`FileNotFoundError` only because `prepare` never created
`input/exact_D500_p3_unseen.csv`. The D500 file, SHA-sealed Garfield
data and prior full-cloud bank are NOT implicated; first eight new native
avalanches were NOT started. User's WSL reported **647 GB free**.

**Code fix:** add local `sample_for_pass(g,box,p)`:
`(offset_for(g,box)+STRIDE*p)%500` with pass range 0..499.
Use it in `reserved`, `check_plan` and `selftest`. The latter
now verifies every 500-pass mapping is a full bijection, old 0..139
and planned 140..219 and reserved 460..499 have disjoint indices,
and agrees exactly with the OLD helper for 100..139. Keep the old
restricted helper unchanged; source identities are the same permutation.

**New safe first-gate retry:**

~~~bash
cd ~/na6p/src/NA6PRoot-mwpc
git pull --ff-only
PY=test_runs/mwpc_stage_c/.venv_pi_slides/bin/python
SCRIPT=test/mwpc-stage-c/pi_p3_fullcloud_learning5_campaign_v1.py
set -euo pipefail
"$PY" -m py_compile "$SCRIPT"
"$PY" "$SCRIPT" selftest

# Do NOT delete any existing campaign directory. The failed attempt
# should have left the root absent, because source selection failed
# before INPUT.mkdir().
CAMPAIGN=test_runs/mwpc_stage_c/pi_p3_fullcloud_learning5_v1
if [ -e "$CAMPAIGN" ]; then
  echo "STOP: existing campaign root; inspect before preparing"; exit 1
fi

"$PY" -u "$SCRIPT" prepare
"$PY" -u "$SCRIPT" preflight
"$PY" -u "$SCRIPT" run --max-chunks 1 --hours 2
"$PY" "$SCRIPT" status
~~~

Stop if *any* guard or command fails. No 1125-chunk background writer
before one SHA-sealed C++ pilot chunk succeeds. This Python patch was
committed on branch `feature/mwpc-stage-c`; real native WSL tests
are pending and must not be claimed as completed.


## 2026-11-04: Native WSL first 8-start gate PASS — continue remaining 1124 chunks

The user executed the corrected source selection in WSL. Actual terminal output:

- `PASS P3 learning5 disjoint original D500 pass ranges and metal guards`.
- `PASS FROZEN new full-cloud P3 learning5 inputs: 9000 D500 original starts`.
- For each P3 geometry: **900 spatial cells x five different original D500 starts**,
  or 4500 starts per geometry. 72,000 reserved candidate passes per geometry
  were SHA-checked; seven old-pilot IDs were excluded for noFW and six for FW75.
  No gas-start inside-metal candidates in this selection.
- `PASS real C++ parser 8 unique unused D500 starts and frozen binary`.
- `FULL MICROSCOPIC GARFIELD p3 ONLY chunk 0 avalanches 8`.
- `PASS SHA FULL BIRTHS 0 ions 31019 zero 0 wall_s 3.4`.
- `SAFE STOP 1 / 1125`; `VERIFIED CHUNKS 1 / 1125`.
- `PASS new learning5 P3 physics status: SHA-sealed avalanches 8 / 9000`;
  `EXACT ION BIRTH POINTS 31019`; four full avalanches in each P3 geometry;
  `NEXT UNUSED CHUNK 1`.

**This is only one real eight-start physical checkpoint, not 9000 finished
avalanches.** The first block's 3.4s wall clock is not a reliable full-run
runtime forecast; future blocks may have much larger gain/time/CSV volume.

Next command to **run only once** after verifying no other active physical
writer; resumes from frozen chunk 1 and preserves completed chunk 0:

~~~bash
cd ~/na6p/src/NA6PRoot-mwpc
git switch feature/mwpc-stage-c
git pull --ff-only
PY=test_runs/mwpc_stage_c/.venv_pi_slides/bin/python
SCRIPT=test/mwpc-stage-c/pi_p3_fullcloud_learning5_campaign_v1.py
"$PY" "$SCRIPT" status

# Must be exactly one physical writer; do not start a second job in parallel.
nohup nice -n 15 "$PY" -u "$SCRIPT" run --max-chunks 1124 --hours 18 \
  > "$HOME/na6p/p3_fullcloud_learning5_rest.log" 2>&1 < /dev/null &
echo "P3 learning5 physical writer PID $!"

# Read-only progress (safe while latest chunk is unsealed/in flight):
tail -n 25 "$HOME/na6p/p3_fullcloud_learning5_rest.log"
"$PY" "$SCRIPT" status
~~~

`run` skips and rechecks the exact same SHA-sealed chunk 0 and will not
overwrite frozen source or checkpoint output. It verifies C++ source and
binary against the recorded `preflight.json` before any new native physics.
The optional 18-hour soft budget stops only between finished chunks; if
reached, repeat the **same** `run` invocation only after the old writer
has definitely exited and no unsealed partial chunk remains.

After `VERIFIED CHUNKS 1125 / 1125` and no active writer, perform the
strict whole-campaign source/checkpoint audit:

~~~bash
"$PY" -u "$SCRIPT" audit \
  | tee "$HOME/na6p/p3_fullcloud_learning5_audit.log"
~~~

Do not assert final complete physics or begin training model until the
actual WSL whole-campaign audit reports 9000 SHA-sealed full Garfield
avalanche records. This new 9000-start campaign is a **first learning-curve
point**, not detector response validation.
