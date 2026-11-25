# 5620-muon 1800 V Heed ROOT export: resumable PBS plan

**Do not repeat the single 8-hour PBS export**: 5638200 and 5639751 were
canceled after nearly zero CPU time with slow access to physical source files.
The v2 detailed log for 5639751 reached the first 54 KB source Heed cluster
CSV, after about 4 minutes resolving original source locks.

## New method

1. Reuse frozen `wiretruth1800_full.root` (5620 muons / 48 successful shards).
2. Spawn 48 small PBS shard workers, limited to 2 concurrent workers. They need
   Python standard library only, not the ROOT/Garfield CVMFS initialization.
3. Each worker verifies its original frozen chunk assignments and receipts;
   streams exactly the four required original physics CSVs **once** through
   SHA-256 onto PBS node-local scratch; reads the local copies; checks original
   Heed cluster / electron / Garfield seed consistency; seals one immutable
   shard with normalized clusters.tsv / primaries.tsv and SHA checksums.
   Persisting each complete shard prevents previous progress loss.
4. After *all* 48 are sealed and passed checksum and population QA, submit a
   separately gated ROOT reducer. It copies the already verified 29 MB ROOT
   into local scratch, concatenates 48 independently sealed TSVs, adds
   `heed_clusters` and `heed_primaries` TTrees with the existing compiled
   converter, and saves one enriched ROOT with SHA checksums.
5. Never modify original physical archives, the original 5620 muon ROOT,
   canceled v1/v2 export areas, or original frozen analysis binaries.

### Commands (WIPP)

```bash
cd ~/NA6PRoot-wiretruth1800-dashboard
git fetch origin feature/p3-1800v-5620-subset-dashboard-20261125
git merge --ff-only FETCH_HEAD
python3 -m py_compile test/mwpc-stage-c/wiretruth1800/export_1800v_heed_sharded.py
bash -n test/mwpc-stage-c/wiretruth1800/export_1800v_heed_sharded.sh
bash -n test/mwpc-stage-c/wiretruth1800/pbs_1800v_heed_shard.sh
bash -n test/mwpc-stage-c/wiretruth1800/pbs_1800v_heed_reduce.sh

export MWPC_1800_SUBSET="$PWD/test_runs/wiretruth1800_SUBSET_5620_20261125T101342Z"
export MWPC_1800_ENRICHED="$PWD/test_runs/wiretruth1800_HEED_5620_SHARDED_v3"

bash test/mwpc-stage-c/wiretruth1800/export_1800v_heed_sharded.sh prepare
bash test/mwpc-stage-c/wiretruth1800/export_1800v_heed_sharded.sh stage-submit
```

Check live progress with:

```bash
bash test/mwpc-stage-c/wiretruth1800/export_1800v_heed_sharded.sh status
tail -n 25 "$MWPC_1800_ENRICHED/pbs-logs/shard_array_0.progress"
```

**Do not run `final-submit` until all 48 shards show complete and the PBS
array is no longer active.** Then:

```bash
bash test/mwpc-stage-c/wiretruth1800/export_1800v_heed_sharded.sh final-submit
tail -n 25 "$MWPC_1800_ENRICHED/pbs-logs/final_root.progress"
bash test/mwpc-stage-c/wiretruth1800/export_1800v_heed_sharded.sh status
```

Final output (only if successfully completed):
`$MWPC_1800_ENRICHED/enriched/wiretruth1800_enriched_5620.root`.

## Physical limitations

- Original complete electron **drift paths** and individual electron
  endpoint coordinates were not saved.
- Saved Heed electron starts use local coordinates; per-anode electron
  endpoint times include avalanche dynamics, not pure drift.
- Original ionization birth reservoir is separate and is *not included*
  in this two-tree HEED enrichment.
- This is 5620 muons, not the 6460 from the original 55-shard plan.

These wrappers have not yet been executed on WIPP. Inspect errors/live logs
and avoid duplicate PBS array submissions.
