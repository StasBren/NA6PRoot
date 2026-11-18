# 2026-11-18 overnight run: frozen 1000 real muons, 2100 V uncapped, K256

**Purpose:** accumulate useful ion-birth ROOT archives for tomorrow's Stage C cathode response without re-generating Geant4 or attempting 1000 extremely large FULL ion archives. Source Geant4 1000 crossings: 334 Jpsi, 333 phi, 333 omega, all 6 stations. Same existing 2100 V uncapped Geant4/TrackHeed/Garfield++ microscopic target that previously produced 3 audited K256 reference events. No Fast4 approximation is used to generate this source; all avalanches are true microscopic Garfield++.

**Caveat:** This is a potentially long laptop run. It cannot be started from GitHub or ChatGPT; the user must launch it in their own WSL. Windows sleep/hibernation or shutting down WSL stops it. There is NO guarantee 1000 will finish by morning. The default is MAX 12 hours of starting new jobs and at most 1000 sequential muons; if an event runs beyond the remaining budget, it may take up to its individual 90-minute timeout. Successful events remain SHA-audited; a timeout, memory shortage or physics mismatch is an unsealed failure, never a valid truncated result.

Each event uses genuine 2100V uncapped Garfield++ microscopic electrons + ions, but archives only K256 randomly sampled ion birth positions/times *per Heed primary* plus exact true ion count, electron endpoint census, primary positions and relevant physical-wire counts. No ion trajectories/shaper/stereo signals are computed in this overnight campaign. K256 is a statistical approximation to the ion *distribution*, never a reduction of microscopic avalanche gain.

## One-time environment and branch

Use the EXISTING non-frozen checkout at `~/na6p/cathode-compact-kernel-2100`; no clone, no modification to frozen `~/na6p/cathode-ion-2100-wsl` or historical source.

~~~bash
cd ~/na6p/cathode-compact-kernel-2100
git status --short
git fetch origin feature/p3-2100v-wsl-overnight-1000-k256-20261118:refs/remotes/origin/feature/p3-2100v-wsl-overnight-1000-k256-20261118
git switch --track origin/feature/p3-2100v-wsl-overnight-1000-k256-20261118
export GARFIELD_INSTALL="$HOME/na6p/install/garfieldpp"
export HEED_DATABASE="$GARFIELD_INSTALL/share/Heed/database"
export LD_LIBRARY_PATH="$GARFIELD_INSTALL/lib:${LD_LIBRARY_PATH:-}"
bash test/mwpc-stage-c/cathode-ion-reference/run_overnight_wsl_2100_k256.sh prepare
~~~

**Only when PREPARED succeeds**, submit the serial run to *your WSL laptop*, disconnect from terminal if desired:

~~~bash
RUNNER=test/mwpc-stage-c/cathode-ion-reference/run_overnight_wsl_2100_k256.sh
LOG="$HOME/na6p/overnight_2100_k256.log"
nohup bash "$RUNNER" run > "$LOG" 2>&1 < /dev/null &
echo "OVERNIGHT_PID=$! LOG=$LOG"
~~~

Do NOT change branch, rebuild binaries or modify the frozen selection while that worker runs; it compares binary/selection/Git SHA before every event. Do NOT run a second copy of the worker; it uses a nonblocking `flock`.

## Read-only status tomorrow

~~~bash
cd ~/na6p/cathode-compact-kernel-2100
bash test/mwpc-stage-c/cathode-ion-reference/run_overnight_wsl_2100_k256.sh status
tail -n 35 "$HOME/na6p/overnight_2100_k256.log"
du -sh test_runs/mwpc_stage_c/overnight_2100_weighted_k256_v1
~~~

The runner refuses to overwrite any incomplete or invalid event directory. When restarting, previously complete SHA-verified events are skipped. A failed event stays in place with a `FAILED` marker and log for inspection. To stop the running worker manually: `kill <OVERNIGHT_PID>`; stop gracefully between events if possible. Do not delete outputs to restart without auditing.

### Safety constraints

- exactly one Garfield muon at a time; per-muon `timeout 90m`, `ulimit -v` 12 GiB virtual address space by default (not physical memory); 8 GiB free-disk reserve; checks at event boundaries; max 12h start-new-event budget. No physical avalanche cap.
- All events are saved under an isolated run folder with `input.csv`, `ion_births.root`, `microscopic_*` CSVs, `physics.log`, `root_audit.log`, `task_metadata.txt`, `SHA256SUMS`, and `complete.txt` only after verification.
- Even if all 1000 microscopic events complete, only 3 original events have independent full-ion truth at this point; sampling errors remain unknown for the 1000 population.
- Running a computationally intense fixed-field 2100V uncapped model is NOT evidence that the actual prototype can operate at that voltage. Space charge/streamers not simulated. Original 2026 Prototype-3 geometry is not the same as the frozen 3mm/pm3mm numerical reference.

## Tomorrow's starting checklist

1. Inspect verified completion count and any failed/unsealed event; do not silently exclude it from statistics.
2. ROOT/C++ summary: true ion births/gain and collected-electron multiplicities, 3 channels and six stations, muon angle/position dependence, outlier failures.
3. Recompute ion RKF waveforms at ~10ns early bins for one same-event K512 benchmark; retain the audited 0.05ns microscopic electron ROOT. Quantify early-time ion+electron response and discretization error without more PCA tuning.
4. Once early response is stable, apply a *documented* selectable VMM3a peaking transfer function (25/50/100/200ns) without mistaking the 25ns DAQ timestamp grouping for peaking or box integration.
5. Use the expanded K256 archives as input to the two-whole-plane model and eventually the 3D U/V strip weighting model, using the existing approximation/error register as the presentation reference.
