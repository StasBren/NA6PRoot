# Independent 2100 V wire-truth pilot: preparation specification (2026-11-24)

Status: **PREPARATION ONLY**. This document does not authorize PBS submission or a physics binary change.

## Target
- 1000 Geant4 muons per channel: Jpsi, Phi, Omega; total 3000.
- Use the *same original frozen Geant4 crossing IDs* as the 1800 V campaign, chunks 0000 through 0049 inclusive (20 muons per chunk per channel), with SHA-256 input comparison.
- Geometry: p3_noFW, 3 mm anode pitch, 30 um anode diameter, cathodes at v=±3 mm, 9-wire analytic patch.
- Gas: Ar/CO2 70/30, Penning r=0.57, lambda=0 cm, T=293.15 K, P=760 Torr, B=0.
- Change anode HV only: 1800 V -> **2100 V**. Preserve other physics and source selection.
- User-reported preliminary mean: sum(ion_births)/sum(heed_primaries) ≈ 30000 at 2100 V. This is not a verified production result or a guarantee of cap safety.

## Provenance and isolation
- **Never** edit or use the 1800 V frozen source root, ROOT shards, binaries, or analysis worktree for this run.
- Use independent worktree and new output root; never overwrite a path that exists.
- The original Geant4 chunks are inputs only, with SHA verification; capture separate C++ source SHA, compiled binary SHA, dependency versions, exact input IDs, counters, PBS resources and job ID.
- Do not import or silently reuse results from other existing 2100 V branches: confirm physics equivalence and exact input identity first.
- Existing 1800 V C++ ROOT analysis requires a separate fix for legitimate in-anode failed starts. Do not recompile its frozen binary in place.

## Physics and audit gates
1. Audit candidate 2100 V implementation for *all* anode-voltage-dependent constants and locks; review actual Garfield++ per-seed `successful`, `ion_births`, `collected_e`, `capped`, `garfield_ne`, `garfield_ni`.
2. Smoke test 1 original chunk (20 muons) **per channel**, collecting wall time, peak RSS and per-seed distributions.
3. Reject any unexplained failed seed; geometrically verified initial starts inside 15 um anode radius can be classified separately, without falsifying original Garfield success flag. Check all other seed and ion-sector conservation counters.
4. Evaluate proportion and max of `capped`; high-gain tail at 2100 V must be assessed before 1000/channel. Mean gain by ratio-of-sums is not a cap test.
5. Only after passing physical and resource review submit exactly 50 chunks per channel, with bounded PBS concurrency and no duplicate submissions.
6. Analyze `G_ion = sum(ion_births) / sum(heed_primaries)` and `G_collected = sum(collected_e) / sum(heed_primaries)` with consistent denominators across matched 1800 V and 2100 V crossing IDs. Also report gain conditional on physical gas starts, separately if needed.
