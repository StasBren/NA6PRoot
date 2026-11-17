# NA60+/DiCE P3 noFW: 1000 true Geant4 muons at 2100 V (UNCAPPED)

Prepared 2026-11-17. **Code committed; no claim of a completed Garfield++/PBS execution.** Separate from the archived 1800 V data, the running 1800 V HPC campaign, and the lower/high-voltage scans.

## Physics and limits

- Input: exactly the same immutable 1000 Geant4 crossings as the wire-truth reference: 334 J/psi, 333 phi, 333 omega, all six stations, with frozen `heed_seed` and original ordering.
- P3 noFW: 3 mm wire pitch, 30 um anode diameter, cathodes at +/-3 mm, anodes +2100 V, cathodes 0 V, Ar/CO2 70/30, Penning transfer Ar r=0.57, B=0.
- Full `Garfield::AvalancheMicroscopic` microscopic electron collisions for every TrackHeed primary, true anode electron endpoints, exact per-wire ion-birth counts/angular moments.
- **NO Garfield avalanche-size truncation**: the new independent source calls `avalanche.DisableAvalancheSizeLimit()`, not `EnableAvalancheSizeLimit(N)`. The CLI must pass `2100 0`; the latter is a mode sentinel, NOT a limit.
- No artificial ion-generation cap is configured. The bounded ion-birth reservoir retains up to **128 uniformly sampled individual ion births per Heed primary**. The ionisation callback still counts **every** ion birth exactly and retains the exact per-wire census. The reservoir uses its own RNG (independent of the Garfield avalanche RNG) and stores `reservoir_weight = N_birth / N_sample` for unbiased per-seed reweighting. It does not change the microscopic physics; this is a *storage sampling limit*, not an avalanche-electron or avalanche-ion production limit.
- The separate **50,000 input-Heed-primary budget** and <=24 crossings per executable invocation are *fail-fast input checks*, not truncation of avalanche electrons/ions. PBS runs exactly one muon per invocation. A bizarre event exceeding 50,000 primaries will fail and must not be silently discarded.
- Garfield++ gain counters and endpoint storage have finite machine representation/memory. A truly arbitrarily large avalanche is impossible to simulate on finite hardware. PBS requests 32 GB RAM, 12h walltime per muon, plus a shorter subprocess timeout. Those safeguards **fail the job**; they do not return a fictitious capped gain.
- Model does **not** include space charge, avalanche saturation, streamer formation, photon-feedback discharge or hardware stability. Results at this voltage are a calculation within fixed-field microscopic physics, not operating recommendations.

## Reconnect and check OLD 1800 V farm run

```bash
ssh -4 stanis@wipp-an1.weizmann.ac.il
qstat -u stanis
qstat -x -u stanis 2>/dev/null | tail -50
find /storage/cgrp/stanis -maxdepth 4 -type d -mtime -5 2>/dev/null | tail -60
```

**A blank `qstat -u` does not prove that all old muons finished.** Find the actual frozen 1800-V campaign's job manifest and per-chunk completion checksums; compare completed original crossing identities (not merely the number of CSV files) against the expected input list. Never overwrite its outputs.

## NEW 2100 V OpenPBS campaign (cluster login)

Clone to a SEPARATE directory and run **pilot first**:

```bash
git clone --depth 1 --branch feature/p3-hv2100-uncapped-ion128-20261117 \
  https://github.com/StasBren/NA6PRoot.git ~/NA6PRoot-hv2100-ion128

cd ~/NA6PRoot-hv2100-ion128
bash test/mwpc-stage-c/hv2100-uncapped/manage_2100_pbs.sh prepare
bash test/mwpc-stage-c/hv2100-uncapped/manage_2100_pbs.sh pilot
```

If the cluster's Garfield/ROOT install is not auto-discovered, load the existing `module`/CMake environment or `export Garfield_DIR=/path/to/lib/cmake/Garfield` first. The manager invokes CMake exactly once and freezes the binary/input/Git hashes.

Pilot submits 5 single-muon array tasks, at most 2 concurrent (`-J 0-4%2`).

Monitor:

```bash
qstat -u stanis
bash ~/NA6PRoot-hv2100-ion128/test/mwpc-stage-c/hv2100-uncapped/manage_2100_pbs.sh status
```

**Only when all five pilot muons are SHA-verified**, and the per-seed gain/time/memory/counter diagnostics are physically and numerically sensible, submit the rest:

```bash
bash ~/NA6PRoot-hv2100-ion128/test/mwpc-stage-c/hv2100-uncapped/manage_2100_pbs.sh full
```

The `full` action schedules indices 5..999 as independent OpenPBS array tasks with max 8 concurrent tasks, and **refuses to run unless the first five are complete and valid**.

The default output root is

`/storage/cgrp/stanis/mwpc_penning057/hv2100_uncapped_ion128_1000_v2/`.

For a different, independent campaign set `MWPC_2100_OUTPUT_DIR` to a NEW directory **before** `prepare`; retain that value for all `pilot/full/status` invocations. Do not change the binary/source/inputs in the middle of a campaign.

Each `muon_NNNN/` contains:
- original two-line input selection;
- `physics.log` containing `GARFIELD_AVALANCHE_SIZE_LIMIT_DISABLED=1` and `ION_BIRTH_RESERVOIR_PER_PRIMARY=128`;
- microscopic crossing, primary-seed, per-wire electron endpoint and exact ion-birth census CSVs;
- `microscopic_ion_birth_sample.csv` with sampled `u_mm,v_mm,w_mm,t_ns`, total births and weight per primary, station-aware `u_global_mm,w_global_mm`, nearest physical wire ID and wire coordinate, and sample size. **128 per Heed primary, not 128 per muon.**
- the PBS worker verifies each primary has min(128, N_birth) coordinate rows with `reservoir_weight=N_birth/N_sample` before marking completion;
- a full per-file `SHA256SUMS`, task metadata and `complete.txt` written **only on fully checked success**.

Tasks with missing/inconsistent outputs stay unsealed. A re-submission cannot overwrite incomplete folders; inspect/quarantine the exact failed index explicitly before re-running. The `status` action validates all 1000 results independently and reports `VERIFIED_COMPLETE=N /1000`. A valid 1000/1000 is not yet a detector response validation.

## Scientific follow-up

1. Inspect direct gain distributions including zeros, `RMS(g)/mean(g)`, positive-only distributions, ion-count/endpoint conservation, multi-anode fraction and rare large avalanches. ROOT/C++ only, ATLAS-like style labelled `NA60+ Simulation`.
2. Compare same 1000 muon identities to reference 1800 V and higher-HV samples with consistent selections; check station, incidence angle and production channel.
3. The 128-ion weighted reservoir is useful for preliminary kernels and Fast4/Fast5 inputs, but is NOT an exact full ion cloud. Validate against a small direct full-ion-birth reference and compare convergence at K=128,256,512 before treating arbitrary cathode waveforms as accurate.
4. If 2100 V regularly produces impractically large avalanches, increasing computing resources does not solve missing space-charge physics. Preserve failed logs and reconsider model scope rather than adding an unreported cap.

## Compatibility / revision

This **ion128/v2** branch must be run in a separate checkout; it does not reuse the earlier `hv2100_uncapped_1000_v1` outputs. Both new results and executable are SHA-locked. Existing 1800 V, old high-voltage scans and PBS outputs remain untouched. **Static source checked; compilation/execution on the cluster still unverified.**
