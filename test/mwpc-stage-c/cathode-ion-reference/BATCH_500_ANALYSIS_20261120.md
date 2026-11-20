# Stage C: first audited 500+ genuine 2100V muon sample

**Code status:** C++/ROOT source committed, not yet compiled in user's WSL. **Physics status:** ongoing resilient 1000-event full microscopic avalanche campaign; log reached completed `muon_0571` on 2026-11-20; this does not establish exactly 572 valid results because some earlier events failed verification. The run is still ACTIVE and must not be interrupted.

## What is available right now

- The frozen selection interleaves channels jpsi40-bore, phi40-bore and omega40-bore, and includes all six stations. First 500 slots already form a useful preliminary channel/station sample.
- For each **sealed and ROOT-audited** output directory: one complete Geant4/TrackHeed/uncapped Garfield++ microscopic muon, exact per-Heed-primary avalanche ion-birth and collected-electron counts, K256 weighted ion-birth coordinates and times, physical collecting-wire census.
- Only one prior event has a complete audited microscopic prompt-electron waveform; the overnight K256 archives do **not** contain whole-muon electron currents, individual ion drift, VMM3a shaping or real stereo strip signals.

## Read-only analysis tool

`summary_audited_2100_muons.cc` reads ONLY closed ROOT `muons/primaries/ions/clusters` trees for events with `complete.txt`, `SHA256SUMS`, a passed ROOT audit and an uncapped physics log. It checks exact frozen crossing identity, true ion count and sampled ion count per Heed primary, and never modifies source data. The tool records the missing/unsealed/failed slot indices explicitly in `_excluded.csv`. It does not interpret failed avalanche events as a physical zero gain. The effective gain proxy reported is `collected_electrons / number_of_Heed_primaries`; gain is `NA` for muons with no Heed primaries. There is NO correction here for known physics uncertainties or run-selection bias.

Outputs: `<prefix>_events.csv` for each verified crossing, `<prefix>_excluded.csv` with every skipped selection slot, `<prefix>_groups.csv` with counts and mean/RMS per channel and station, `<prefix>_hist.root` and `<prefix>_overview.png` (station coverage, effective-gain distribution, true avalanche ion births). Plot values are current uncalibrated simulation, not measured detector signals.

## Run while overnight campaign is still running

**DO NOT switch branches, pull, rebuild, or edit the active overnight checkout.** Use a second Git worktree to avoid any lock change or risk to running task.

~~~bash
cd ~/na6p/cathode-overnight-2100
git fetch origin feature/p3-2100v-audited-500-muon-summary-20261120
git worktree add --detach ~/na6p/cathode-analysis-2100 FETCH_HEAD
cd ~/na6p/cathode-analysis-2100
export LD_LIBRARY_PATH="$HOME/na6p/install/garfieldpp/lib:${LD_LIBRARY_PATH:-}"
cmake -S test/mwpc-avalanche -B test/mwpc-avalanche/build-analysis \
  -DCMAKE_PREFIX_PATH="$HOME/na6p/install/garfieldpp"
cmake --build test/mwpc-avalanche/build-analysis \
  --target mwpc_summary_audited_2100_muons -j2
~~~

Snapshot the first 572 selected event slots as of the latest reported DONE index 571. The currently running event 572 is excluded. If there are more recent verified events, they can be analyzed in a new non-overwriting prefix later.

~~~bash
DATA="$HOME/na6p/cathode-overnight-2100/test_runs/mwpc_stage_c/overnight_2100_weighted_k256_v1"
SELECT="$PWD/test/mwpc-stage-c/p3_wiretruth_1000_selection.csv"
BIN=./test/mwpc-avalanche/build-analysis/mwpc_summary_audited_2100_muons
OUT="$PWD/test_runs/mwpc_stage_c/batch_2100_statistics_v1"
mkdir -p "$OUT"
"$BIN" --data "$DATA" --selection "$SELECT" \
  --stop-index 572 --out-prefix "$OUT/first572_slots"
cat "$OUT/first572_slots_groups.csv"
head -n 15 "$OUT/first572_slots_excluded.csv"
ls -lh "$OUT/first572_slots_overview.png"
~~~

Read-only status check on the ORIGINAL runner remains available:

~~~bash
cd ~/na6p/cathode-overnight-2100
bash test/mwpc-stage-c/cathode-ion-reference/run_overnight_wsl_2100_k256.sh status
tail -n 20 ~/na6p/overnight_2100_k256_resilient.log
~~~

## What the next scientific steps are

1. Validate channel and station coverage and gain/ion-birth distribution across hundreds of actual muons; examine outliers and audit failures as *separate* categories.
2. Choose a modest representative event subset (low/median/high gain and both cathode-side populations) for direct RKF whole-plane waveforms and prompt-electron replays. Do not perform full per-ion RKF for all ~million births per event.
3. Independently, compute a 10ns direct RKF ion response for the already verified K512 benchmark muon; compare bin-integration vs 100ns fast template-shifting to quantify the known large first-bin current bias, then examine VMM3a peaking-time transfer function.
4. Adapt 2-plane model to actual U/V stereo strips with validated numerical weighting fields, and update the approximation/error register for the presentation.

**Important geometry caveat:** current frozen simulation is analytic nine-wire P3-noFW at 3mm pitch, ±3mm cathodes, 2100V and provisional ion mobility, not a full as-built prototype field. Data is simulation truth conditional on this model, not detector calibration.
