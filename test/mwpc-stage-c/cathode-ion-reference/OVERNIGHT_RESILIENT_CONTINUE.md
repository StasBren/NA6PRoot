# Continue all 1000 original muon indices despite per-event failures (2026-11-19)

The previous WSL run halted on muon 181 when the verifier correctly flagged two `successful=0` primaries. This is **not** a reason to terminate the entire simulation: do not change physics or relax the verifier. The old process is already running from index 182. The new `run_overnight_wsl_2100_k256_resilient.sh` is designed to be submitted at any time and **waits on the existing `.runner.lock`**. Once the old worker exits or stops, it takes the same physics lock and tries every original event index from 0000 to 0999. Already verified events are skipped, failed or unsealed directories are preserved and skipped without overwriting. Only missing entries are simulated.

Per-muon errors that do NOT end the run: nonzero microscopic simulator exit (including timeout or OOM), failed ROOT audit/uncapped-physics proof, missing/incomplete event artifacts, or checksum failure. Each newly failed event retains `FAILED`, physics and audit logs. A failed event never receives a `complete.txt`. Audit criteria are unchanged.

CAMPAIGN-wide failures still stop intentionally: changed Git checkout or simulation/checker binary, insufficient disk, unreadable frozen input selection, or inability to create a new output directory. These cannot safely be ignored. The default new runner has `MWPC_OVERNIGHT_MAX_HOURS=0` (no time cap), `MWPC_OVERNIGHT_MAX_EVENTS=1000`, per-event timeout 90 minutes, 8 GiB free-disk reserve and the same 12 GiB address-space limit. It makes **up to 1000 original-input attempts**, not a guarantee of 1000 verified results.

## Safe takeover WITHOUT switching branches or stopping the active worker

The current `~/na6p/cathode-overnight-2100` worktree is detached at the original fixed Git SHA `d6ec458`, and `config.lock` hashes that Git commit and the two science binaries. **Do not git switch, git merge or git pull in that worktree while the old worker is running.** The commands below fetch the new file but do NOT change HEAD or the physics binaries.

~~~bash
cd ~/na6p/cathode-overnight-2100
git fetch origin feature/p3-2100v-wsl-skip-failed-1000-20261119
RESILIENT=test/mwpc-stage-c/cathode-ion-reference/run_overnight_wsl_2100_k256_resilient.sh
git show FETCH_HEAD:"$RESILIENT" > "$RESILIENT"
bash -n "$RESILIENT"
echo "HEAD=$(git rev-parse --short HEAD) — MUST stay d6ec458"
LOG="$HOME/na6p/overnight_2100_k256_resilient.log"
nohup bash "$RESILIENT" run > "$LOG" 2>&1 < /dev/null &
echo "Resilient worker PID=$!"
tail -n 5 "$LOG"
~~~

`WAIT_FOR_RUNNING_WORKER` is expected while the current worker is still active, and **does not mean the corrected one is stuck**. After the old worker exits, the resilient runner prints `RESILIENT_START`, skips SHA-verified events and pre-existing failed/unsealed folders, and processes missing events up through 0999.

Do NOT start a third worker. The lock prevents concurrent event execution, but duplicates complicate monitoring. You may close the terminal with `nohup`, but keep the laptop awake and plugged in.

## Monitoring

~~~bash
cd ~/na6p/cathode-overnight-2100
RESILIENT=test/mwpc-stage-c/cathode-ion-reference/run_overnight_wsl_2100_k256_resilient.sh
bash "$RESILIENT" status
tail -n 25 "$HOME/na6p/overnight_2100_k256_resilient.log"
tail -n 15 "$HOME/na6p/overnight_2100_k256_resumed.log"
~~~

The status counters are `VERIFIED`, `FAILED`, `UNSEALED`, and `ABSENT`; they sum to 1000. The campaign has attempted the last index when `ABSENT=0` and the resilient process logs `RESILIENT_END` (not necessarily 1000 verified outputs). Retain muon 181's failed archive and any subsequent failed records; classify them in the analysis rather than treating them as good simulated muons.
