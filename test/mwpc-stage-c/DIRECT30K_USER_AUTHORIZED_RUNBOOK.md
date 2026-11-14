# User-authorized direct production: 10,000 muons/channel without pilot1000

This is a **separate, version-controlled PBS orchestration branch**.
The scientific physics worktree and successful frozen campaign remain at
`599e29522edb6fcec9dfcca2904308cb3186edce`, with actual binary
SHA-256 `2cf92cb27d3a180cd26d3ab5417f57293473864f6d685ee43bacf931712aa8f1`.
Do not merge or pull this control branch into the physics worktree.

## Scientific waiver and limitations

On 2026-11-14 the user explicitly requested **direct 3x10,000**
(Jpsi, Phi, Omega) without the intermediate 1,000-muon resource/science gate.
The real 20-muon gate was passed with 0 caps, 0 failed returns, 1,169
microscopic avalanches, Garfield walltime 576.53 s, measured maximum RSS
183,544 KiB and end-to-end PBS walltime 38m25s.
The 1,000-muon gate is **not passed and must not be fabricated**.

Risk deliberately accepted: a 20-muon resource sample cannot bound rare
high-Heed tails; some 20-muon chunks could exceed 4GB RAM / 8h walltime.
Total full run CPU/walltime and I/O variance are uncertain. Each PBS worker
keeps the frozen native SHA and independently audits all per-wire electron
collection, ion births, source crossing IDs, full-Heed replay and output CSVs.

Each of the 3 channels schedules 500 distinct 20-muon chunk indices 0..499.
Jpsi chunk 0000 is already present and must be SHA-audited and skipped.
All other 1,499 chunks are new physics computations. 500 tasks/channel,
maximum 2 concurrent/ch; at most 6 microscopic jobs running concurrently.

## Prepare a SEPARATE clean control worktree (on wipp-an1)

```bash
cd "$HOME/NA6PRoot-mwpc"
git fetch origin \
 refs/heads/feature/p3-wiretruth-30k-direct-full-20261114:refs/remotes/origin/feature/p3-wiretruth-30k-direct-full-20261114
git worktree add -b wiretruth30k-direct-control \
 "$HOME/NA6PRoot-wiretruth30k-direct-control" \
 origin/feature/p3-wiretruth-30k-direct-full-20261114

WT="$HOME/NA6PRoot-wiretruth30k"
CTRL="$HOME/NA6PRoot-wiretruth30k-direct-control"
ROOT="/storage/cgrp/stanis/mwpc_penning057_wiretruth30k_v1"
DIRECT="$CTRL/test/mwpc-stage-c/submit_wiretruth30k_direct_v1.sh"
WORKER="$CTRL/test/mwpc-stage-c/pbs_wiretruth30k_direct_v1.sh"

git -C "$WT" rev-parse HEAD
git -C "$CTRL" log -1 --oneline
bash -n "$DIRECT"
bash -n "$WORKER"
ls -l "$ROOT/gates/gate_20.json" "$ROOT/locked/binary.json" \
 "$ROOT/locked/selection.json"
```

The existing physics worktree must remain on exactly
`599e29522edb6fcec9dfcca2904308cb3186edce`. Do not run
`git pull` on `$WT`. Do not remove the old 30k campaign source archives.

## Submit each FULL channel, explicitly

```bash
bash "$DIRECT" Jpsi
bash "$DIRECT" Phi
bash "$DIRECT" Omega
```

Save all three actual IDs. Each command submits a PBS array with
`-J 0-499%2` using the **direct** worker script, not the original
`submit_wiretruth30k_stage_v1.sh`. The direct worker has a real pilot20 gate
requirement and **does not require or invent** `gate_1000.json`.

The runner records source Git SHA and separate control branch Git SHA,
refuses to run after a change to either worktree, and will not resubmit an
already submitted channel from the ledger without manual review.

## Monitor actual calculations

```bash
/opt/pbs/bin/qstat -t -u "$USER"
tail -n 10 "$ROOT/submissions/qsub_ids.tsv"
python3 "$WT/test/mwpc-stage-c/pi_wiretruth30k_production_v1.py" \
 status --root "$ROOT"
# For one chunk:
ls -lah "$ROOT/results/jpsi/chunk_0000/"
# For failed chunks, inspect $ROOT/logs, $ROOT/attempts, and real array exit codes.
```

No gates should be run while the requested chunks are incomplete. After
all actual 1,500 results are complete and independently audited:

```bash
python3 "$WT/test/mwpc-stage-c/pi_wiretruth30k_production_v1.py" \
 gate --root "$ROOT" --stage full
python3 "$WT/test/mwpc-stage-c/pi_wiretruth30k_production_v1.py" \
 report --root "$ROOT" --old-root /storage/cgrp/stanis/mwpc_penning057
python3 "$WT/test/mwpc-stage-c/pi_wiretruth30k_production_v1.py" \
 export --root "$ROOT"
```

Resource projections based only on the real first 20 muons: Garfield-only
~240 CPU-hours; raw original per-chunk CSV ~24.3 GiB for 1,500 chunks.
The PBS pilot 20 took 38m25s overall; under 6 concurrent jobs,
idealized total walltime could be several days, and congestion/tails may
lengthen it. Do not claim a completed campaign until all checks pass.
