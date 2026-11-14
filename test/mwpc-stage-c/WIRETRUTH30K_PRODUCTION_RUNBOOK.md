# NA60+ MWPC — 30k true per-anode Penning057 reference

Production branch: `feature/p3-wiretruth-30k-production-20261114`.
Base: `feature/p3-wiretruth-1000-20261114` (C++ and CMake target unchanged).
Physics: real 3x10k Geant4 crossing inputs from the previous **frozen**
original2m Jpsi/Phi/Omega campaign; NEW isolated Garfield/TrackHeed simulation.

**Important:** committing production scripts is NOT evidence that the
Garfield++ target compiled or any PBS task ran. This README defines the exact
HPC gates and produces job IDs only when the human runs qsub on wipp-an1.

## 0. Verify HPC and installed old pilot, without changing any worktree

Connect F5 VPN, then from the laptop:

```bash
ssh stanis@wipp-an1
hostname
/opt/pbs/bin/qstat -u "$USER"
/usr/wipp/bin/qtop
ls -ld /storage/cgrp/stanis/mwpc_penning057/original2m_{jpsi,phi,omega}_10000_v1
find "$HOME/NA6PRoot-mwpc/test_runs/mwpc_stage_c/p3_noFW_penning057_wiretruth_1000_v1/chunks" \
  -name complete.json 2>/dev/null | wc -l
```

The last command indicates the number of original 1000-pilot receipts in
THAT path only; 50 receipts plus a successful pilot audit are needed to call
it finished. This is NOT the new 3x10k production path.

## 1. Isolated Git worktree on the HPC host

Do not switch/reset the pre-existing `~/NA6PRoot-mwpc` worktree, because it
may have valuable uncommitted modifications. Create a NEW worktree:

```bash
cd ~/NA6PRoot-mwpc
git fetch origin feature/p3-wiretruth-30k-production-20261114
git worktree add -b wiretruth30k-local \
  "$HOME/NA6PRoot-wiretruth30k" \
  origin/feature/p3-wiretruth-30k-production-20261114
WT="$HOME/NA6PRoot-wiretruth30k"
ROOT=/storage/cgrp/stanis/mwpc_penning057_wiretruth30k_v1
SCRIPT="$WT/test/mwpc-stage-c/pi_wiretruth30k_production_v1.py"
SUBMIT="$WT/test/mwpc-stage-c/submit_wiretruth30k_stage_v1.sh"
git -C "$WT" status --short
git -C "$WT" rev-parse HEAD
python3 -m py_compile "$SCRIPT"
python3 "$SCRIPT" selftest
bash -n "$WT/test/mwpc-stage-c/pbs_wiretruth30k_build_v1.sh"
bash -n "$WT/test/mwpc-stage-c/pbs_wiretruth30k_chunk_v1.sh"
bash -n "$SUBMIT"
```

Selftest covers station-0 +1.5 mm anode phase, station-1..5 zero phase,
negative global wire IDs, microscopic charge/ion census, and GNU time parser.

## 2. Existing C++ target BUILD on a PBS compute node

This uses exactly the same CVMFS Garfield environment and 1-core PBS
settings as the earlier original2m campaign. No new physics C++ code.
The build job stages all 1,500 original 20-muon input chunks, checks old
manifest SHA-256, locks source commit, binary and Garfield dynamic libraries,
ROOT and compiler versions. Original 30k results are read-only.

```bash
bash "$SUBMIT" build
cat "$ROOT/submissions/qsub_ids.tsv"
```

**Record the REAL PBS job ID printed by qsub**, monitor with `qstat`,
then wait for the build job to finish successfully. Only then:

```bash
python3 "$SCRIPT" runtime-check --root "$ROOT"
find "$ROOT/inputs" -name 'chunk_*.csv' | wc -l
cat "$ROOT/locked/binary.json"
```

If CVMFS libraries change, STOP. Do not mix results from changed physics
versions. Make a NEW campaign folder/manifest if rebuilding is necessary.

## 3. Twenty-muon scientific pilot: one PBS job

```bash
bash "$SUBMIT" pilot20
/opt/pbs/bin/qstat -t -u "$USER"
```

The worker carries out `--validate-input`, `--validate-heed`, full
`--wiretruth-p3-noFW` with 128 ion birth reservoir points per seed, and
an independent audit of collecting-wire endpoint counts, ion birth regions,
seed and muon charge conservation, zero-gain and failures, angular-sector
asymmetries, exact per-wire w/t population moments, global station parity,
uniform reservoir weights, and original TrackHeed replay. It runs through
PBS **not** on the login node.

After pilot job completion:

```bash
python3 "$SCRIPT" audit-one --root "$ROOT" --channel Jpsi --index 0
python3 "$SCRIPT" gate --root "$ROOT" --stage 20
cat "$ROOT/gates/gate_20.json"
```

Review the measured peak RSS, wall seconds, data volume, count of
`AvalancheElectron()` false returns, and `source_old_heed_replay` status.
Matching TrackHeed is tested against the old SHA-sealed primary CSVs at
the same original crossing ID, primary ID, cluster, coordinates and times.
If these do NOT match, all new primaries are still archived and must be
used for future Fast4/Fast5-T comparisons.

## 4. One-thousand-muon cross-channel pilot

Only after `gate_20.json` is successful. It covers 400 Jpsi, 300 Phi and
300 Omega from the same frozen 30k selection: Jpsi chunks 0..19, Phi/Omega
chunks 0..14. Completed Jpsi chunk 0 is NOT recalculated.

```bash
bash "$SUBMIT" pilot1000
/opt/pbs/bin/qstat -t -u "$USER"
python3 "$SCRIPT" status --root "$ROOT"
```

After all 50 chunks finish and pass SHA audits:

```bash
python3 "$SCRIPT" gate --root "$ROOT" --stage 1000
cat "$ROOT/gates/gate_1000.json"
```

The gate checks actual 4 GB per-worker memory and 8 h PBS wall limits
(with margins). It estimates total 1,500-chunk CPU-wall hours and raw
data size from actual `/usr/bin/time -v` measurements. Assess tail
latencies, Garfield failure flags, disk quota and PBS availability before
starting full production. The projection is not a promised completion time.

## 5. The remaining 29,000 muons, only after the 1000 gate

```bash
bash "$SUBMIT" full
/opt/pbs/bin/qstat -t -u "$USER"
/usr/wipp/bin/qtop
python3 "$SCRIPT" status --root "$ROOT"
```

Three independent arrays, limited to **2 simultaneous jobs per channel**:

- Jpsi: chunks 20..499 (480 x20)
- Phi: chunks 15..499 (485 x20)
- Omega: chunks 15..499 (485 x20)

The pilot 50 chunks are excluded, and the final 1,500 chunks are unique.
If interrupted: completed chunk receipts are SHA-checked/skipped.
Failed attempts are preserved under `attempts/` and never overwritten;
workers claim each chunk to avoid concurrent duplicates. The stage
submission ledger prevents accidental duplicate array submissions.
Use `WT_RESUBMIT=1` only after inspecting job and chunk state, never
to indiscriminately submit a second concurrent full array.

Monitor:

```bash
cat "$ROOT/submissions/qsub_ids.tsv"
ls -lt "$ROOT/logs" | head
find "$ROOT/attempts" -name worker.log | tail
/opt/pbs/bin/qstat -t -u "$USER"
python3 "$SCRIPT" status --root "$ROOT"
```

No claim lock should be removed without first confirming that its worker
is absent from qstat. Failed/incomplete staged outputs need manual
quarantine, not blind deletion of physics data.

## 6. Final physics audit, reports and compact export

```bash
python3 "$SCRIPT" gate --root "$ROOT" --stage full
python3 "$SCRIPT" report --root "$ROOT"
python3 "$SCRIPT" export --root "$ROOT"
du -sh "$ROOT/compact" "$ROOT/results"
ls -lh "$ROOT/compact/"
```

The compact folder includes:

- `{jpsi,phi,omega}_wire_seed.csv.gz` and `_wire_muon.csv.gz`.
- `_seed.csv.gz`, `_heed_primaries.csv.gz`, `_heed_clusters.csv.gz`,
  `_timing.csv.gz`, `_muon_timing.csv.gz`, `_crossing.csv.gz`.
- `_geant4_crossings.csv.gz`: exact full 3D Geant4 entry/exit from
  frozen 30k original input.
- `wire_sharing_by_angle.csv`: n=0,1,2,3+ wire census and mean
  second-wire charge fraction Q2/Qtotal by |theta_u| bins, separately
  for Jpsi, Phi, Omega, All.
- `previous_30k_sanity_comparison.csv`: old/new per-channel Heed,
  collected electron and ion birth totals. Different Heed streams should
  be reported, not forced to match.
- `run_manifest.json`: locked build/versions/physics, 1500 chunks,
  actual PBS IDs, SHA for nine native CSVs each, old Heed replay flags,
  strict Jpsi calibration and Phi/Omega validation separation.
- `export_sha256.json`: fingerprints of all 27 compact CSV.gz outputs.

The uncompressed original per-avalanche 128-point birth reservoirs are
preserved ONLY under the full HPC results folder, to avoid unnecessary
download; the 128-point sample has unbiased weights N_ion/N_samples for
each avalanche, but its adequacy for cathode-signal simulation is NOT
established. Before replacing full ion-cloud statistics by these samples,
measure convergence against a larger reservoir (e.g. dedicated 512-point
or all-birth small study). No Fast4/Fast5-T validation is performed by
these production scripts; that comparison requires replaying archived
new Heed primaries and strict holdout on Phi/Omega.

## 7. Download summaries to the local WSL machine

Exit the HPC SSH session and in WSL with F5 VPN connected:

```bash
mkdir -p ~/na6p/wiretruth30k_reference
rsync -avP \
  stanis@wipp-an1:/storage/cgrp/stanis/mwpc_penning057_wiretruth30k_v1/compact/ \
  ~/na6p/wiretruth30k_reference/
explorer.exe "$(wslpath -w ~/na6p/wiretruth30k_reference)"
```

Do not copy the huge original ion reservoirs unless needed.

**No PBS job IDs are known until `qsub` has actually been run on HPC.**
They are stored in `$ROOT/submissions/qsub_ids.tsv` and in the final
per-chunk scientific manifests, not invented here.
