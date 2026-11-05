# p3_noFW — 10,000 additional real J/psi parent muons, full microscopic Garfield++

## Purpose

Add **10,000 new muon gas crossings** from **10,000 different generated
J/psi parent events**, all **disjoint** from the previous fully SHA-sealed
1000-muon p3_noFW benchmark, using real archived Geant4
muon chords, real TrackHeed clusters/electrons and all-seed Garfield
AvalancheMicroscopic. Previous 1000 and all D500/A20/full-cloud
archives are untouched. This is a separate independent production
campaign, **not** a 1000→10000 replacement (total becomes 11000).
Use exactly the existing built
\`mwpc_phase_b_native_p3_noFW_jpsi_bench_v1\` physical executable.
No Geant4 rerun, no C++/Garfield rebuild, no other geometries.

The original 1000 generated 54,912 real Garfield avalanches, 31,043,558
ion births, 0 capped, 2 unsuccessful. The benchmark's measured
4.997s/real muon predicts ~49,970s (13.88 hours) for 10k new real
muons on the SAME computer, but distribution tails and disk writing
cost can change that. Soft time budget checks only BETWEEN completed
20-parent blocks, not during an individual parent.

## Frozen independent parent source

The selection uses SHA-verified ORIGINAL saved Geant4 J/psi native
crossing exports, chunk3 onward. It verifies frozen source manifests
and excludes all 1000 prior parent keys \`(source,chunk,event)\`.
The original chunk2 old/new120 tests are not touched (chunk2 is excluded).
Choose one eligible crossing per generated parent event, SHA rank
independently of any subsequent Garfield outcomes, then SHA-rank
new unique parents. If chunk3 has enough parents, only it is used;
otherwise the script searches later already-exported chunks up to 40.
**Missing chunks are not invented or regenerated**: if fewer than
10k additional eligible parents exist in available original exports,
\`prepare\` fails BEFORE creating any new campaign root.
You can export additional previously simulated original ROOT chunks
and retry \`prepare\` (with no change to Geant4 physics).

As with the original 1000, the selection is one gas crossing per unique
J/psi generated parent; it is NOT the entire physical all-chamber
crossing flux, and should not be interpreted as full detector efficiency.

## Preparation, preflight, first 20

~~~bash
cd ~/na6p/src/NA6PRoot-mwpc
git switch feature/mwpc-stage-c
git pull --ff-only

PY=test_runs/mwpc_stage_c/.venv_pi_slides/bin/python
SCRIPT=test/mwpc-stage-c/pi_p3_noFW_jpsi_extra10000_native_v1.py

"$PY" -m py_compile "$SCRIPT"
"$PY" "$SCRIPT" selftest

# Disk capacity: 10k means ~10x the prior 1000's 128-point/avalanche CSV.
df -h ~/na6p
du -sh test_runs/mwpc_stage_c/pi_p3_noFW_jpsi1000_native_v1/native_chunks

# Creates NEW 10k frozen input; old 1000 used only for exclusion.
"$PY" -u "$SCRIPT" prepare --start-chunk 3 --end-chunk 40
"$PY" -u "$SCRIPT" preflight

set -o pipefail
# Compute just 20 actual independent new J/psi muons first.
"$PY" -u "$SCRIPT" run --max-chunks 1 --hours 2 \
  | tee ~/na6p/p3_noFW_extra10000_first20.log
"$PY" "$SCRIPT" status
~~~

Do not advance if preflight or the first physics chunk fails, produces
unexpected avalanche caps, or the disk is becoming full. The first
20 are not pilot-only: they belong to the frozen 10000 sample and will
be counted in the final result.

## Remaining 9980 in the background (only after first chunk PASS)

~~~bash
cd ~/na6p/src/NA6PRoot-mwpc
PY=test_runs/mwpc_stage_c/.venv_pi_slides/bin/python
SCRIPT=test/mwpc-stage-c/pi_p3_noFW_jpsi_extra10000_native_v1.py

# One native writer only; completes at most 499 NEW blocks of 20.
nohup nice -n 10 "$PY" -u "$SCRIPT" \
  run --max-chunks 499 --hours 20 \
  > ~/na6p/p3_noFW_extra10000_remaining.log 2>&1 < /dev/null &
echo "native Garfield runner PID $!"

# Progress from another shell:
"$PY" "$SCRIPT" status
tail -n 25 ~/na6p/p3_noFW_extra10000_remaining.log
pgrep -af 'pi_p3_noFW_jpsi_extra10000_native_v1.py'
~~~

Do not launch a second writer while one is running. If the budget
expires between 20-parent chunks, each complete chunk has an immutable
\`complete.json\` SHA of all original native data. Resume using the
same \`run\` command; it reaudits and skips the 20-parent sealed
chunks and continues from the first unused chunk.

**If a process dies DURING an unsealed 20-parent chunk**, a folder may
contain partial physics. The runner intentionally refuses to overwrite
it; examine process status and \`physics.log\` before any recovery.
Do not \`rm\`/\`mv\` original data. The previous 1000 are preserved.

## Completion and outputs

~~~bash
"$PY" "$SCRIPT" status
"$PY" "$SCRIPT" audit
~~~

Pass: \`PASS ALL 10000 NEW JPSI P3 NOFW MUONS SHA-SEALED\`
plus independent rehash audit.

Output in
\`test_runs/mwpc_stage_c/pi_p3_noFW_jpsi_extra10000_native_v1/\`:
- \`input/native_muon_crossings_input.csv\`: selected SHA-frozen 10k
  distinct original Geant4 muon gas chords and Heed RNG seeds.
- \`input/input_manifest.json\`: hashes of source original ROOT/CSV
  chunks and previous 1000 exclusion proof.
- \`native_chunks/chunk_0000..0499/\`: for each 20 true muons:
  native Garfield crossing, seed census, Heed primaries/clusters,
  per-seed and per-muon timings, 128-point-per-avalanche unbiased
  ion BIRTH reservoir, physics.log, source CSV and SHA complete.json.
- **All ionisation births counted exactly** per avalanche/whole muon,
  but only up to 128 \`(u,v,w,t)\` positions stored per avalanche;
  no induced cathode signal, no ion drift or strip response.

The original 1000 used ~5.65 million saved ion position samples.
This additional 10k will likely create order 50–70 million saved
position records, several GB of CSV. Check sufficient disk capacity
before an overnight run, ideally tens of GB free. This is a full
physical Garfield runtime campaign, not a full uncompressed ion-cloud
archive. No additional simulated muon parents are secretly selected
from the first 1000.
