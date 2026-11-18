# WSL laptop pilot — same 2100 V uncapped cathode-ion birth archive as PBS

**This is an input-generation pilot, not an individual cathode-strip signal simulation.** The original C++ Garfield++/ROOT executable and validation/checkers are reused unchanged from `feature/p3-2100v-cathode-ionbirth-reference-20261117`. Only a separate WSL runner is added.

## Commands in WSL (not on the cluster)

```bash
git clone --depth 1 --branch feature/p3-2100v-cathode-ionbirth-wsl-20261118 \
  https://github.com/StasBren/NA6PRoot.git ~/na6p/cathode-ion-2100-wsl
cd ~/na6p/cathode-ion-2100-wsl
RUN=test/mwpc-stage-c/cathode-ion-reference/run_local_cathode_reference.sh
bash "$RUN" prepare
bash "$RUN" weighted 0 256
bash "$RUN" status
```

Only after a complete, SHA-verified weighted run, and inspection of its `resource_time.txt` and `root_audit.log`, consider the much heavier full archive for the **same muon**:

```bash
bash "$RUN" full 0
bash "$RUN" compare 0 256
```

To check reservoir-size convergence, run `bash "$RUN" weighted 0 512`, then `bash "$RUN" compare 0 512`; similarly for K=1024. All run the exact same original Geant4 muon and Garfield random seeds, with independent sampler RNG.

## Safeguards and output

- Separate build: `test/mwpc-avalanche/build-cathode-2100-wsl/`
- Separate output: `test_runs/mwpc_stage_c/cathode_2100_laptop_v1/`; nothing under `/storage/cgrp/stanis` is touched
- Default per-muon external timeout 3600 seconds. Change for a **new** invocation with `MWPC_LOCAL_TIMEOUT_SEC=7200 bash "$RUN" weighted 0 256` if necessary. Timeout/OOM means FAILED and does not produce a valid physics event
- Optional process virtual-memory ceiling: `MWPC_LOCAL_MAX_VMEM_KB=<kilobytes>`, but this protects the computer at the cost of aborting the event, never artificially capping the avalanche
- `DisableAvalancheSizeLimit()` remains in the Garfield++ source. No ion/electron multiplication cutoff is introduced
- Each event has `input.csv`, `ion_births.root`, microscopic per-seed and per-wire CSVs, `physics.log`, `root_audit.log`, `resource_time.txt`, `metadata.txt`, `SHA256SUMS`, and `complete.txt` only after successful validation
- Existing incomplete directories are not overwritten; choose a new `MWPC_LOCAL_OUTPUT_DIR` if you want a fresh run

## Interpretation

The `weighted` mode saves a uniform random reservoir of min(N,K) actual ion births per Heed primary, with inverse inclusion weights N/n. The `full` mode streams every birth to a compressed ROOT TTree and preserves a first-primary whole-cathode electron control. Neither mode has the necessary validated physical stereo-strip weighting field, calibrated ion species/mobility or electronics; those are separate response-simulation stages.

All QA/statistical calculations and plots in this workflow are C++/ROOT, not Python. A successfully audited archive is a **birth-input reference**, not a completed cathode response.

**Status:** WSL script syntax checked with `bash -n`; the actual Garfield++ compile, physics run and ROOT audit require the user's WSL environment and are not claimed as performed.
