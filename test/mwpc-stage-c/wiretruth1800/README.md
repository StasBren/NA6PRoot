# 1800 V real-Geant4 ~8,000-muon wiretruth distribution analysis

**Status 2026-11-23:** New C++/ROOT analysis code + file-integrity snapshot
script are committed in an isolated branch. Not yet compiled/tested on the
WIPP runtime, and no analysis result can be claimed before actual PBS output.
Zero new Garfield++ muons are simulated.

## Scientific source and identity

Only existing **completed and SHA-verified** production chunk results from

`/storage/cgrp/stanis/mwpc_penning057_wiretruth30k_v1/results/{jpsi,phi,omega}/chunk_XXXX`

are included. Original independent source crossing bytes remain frozen at

`/storage/cgrp/stanis/mwpc_penning057/original2m_{channel}_10000_v1/chunks/chunk_XXXX.csv`.

Every recorded chunk contains exactly 20 original Geant4 crossings.
The 2026-11-16 status was around 131 Jpsi, 136 Phi and 50 Omega chunks;
the Omega 50–134 PBS array was launched to finish about 135.
The exact number analyzed is **whatever successful, complete receipts exist
at the time of the SHA-frozen snapshot**. It is never assumed to be 8,000.

Original physics lock:
- production Git `599e29522edb6fcec9dfcca2904308cb3186edce`
- original binary SHA `2cf92cb27d3a180cd26d3ab5417f57293473864f6d685ee43bacf931712aa8f1`
- 1800 V, 9-wire analytic P3 noFW, 3 mm physical wire pitch,
  30 μm diameter, cathodes at v=±3 mm, Ar/CO₂ (70/30), Penning r=0.57
- station zero physical anode phase = +1.5 mm;
  all other stations physical anode phase = 0 mm.
- These numerical values do **not** describe the final as-built MWPC
  chamber geometry; distinguish them in all reports and future parameter scans.

## What gets analyzed

Full, **exact**, native Garfield++ counts (not ion-reservoir estimators):

- `microscopic_crossing.csv`: independent original Geant4 muon totals;
- `microscopic_seed.csv`: actual Heed-primary Garfield electron/ion census;
- `microscopic_wire_seed.csv`: individual primary × actual physical anode,
  exact collected-electron endpoints and nearest-ion-birth assignments;
- `microscopic_wire_muon.csv`: complete muon × actual physical anode
  counts, ion-hemisphere moments and station-aware anode indices.

The Python section of the Bash script is restricted to **file integrity,
selection and SHA256**, not physics/statistical analysis. Before adding
a chunk, it compares the frozen source input hash and the four required
CSV hash+size values to its `verified_report.json` and immutable production
selection/binary lock. Only clean completed records enter the manifest.

C++/ROOT independently tests:
1. each chunk has 20 original independent Geant4 muons;
2. physical station, nearest wire, global wire centre and patch `-4..4`;
3. every Garfield primary has `successful=1`, `capped=0`;
4. exact per-primary `(Q_e,N_ion)` = sum across its physical wires;
5. exact per-muon `(Q_e,N_ion)` = sum across primaries = sum across physical
   wires = source Garfield crossing;
6. positive-ion up+down+boundary = exact total ion births on wire;
7. positive-ion count on nearest birth wire is distinct from collected
   electron endpoint count on actual collecting wire.
   
**All original crossings are included as independent observations**, even
if a muon has zero collected charge. Events where an anode has ion births but
zero collected electrons are also retained.

## Physics distributions

The ROOT file has four TTrees:
- `muons`: Geant4 channel, det/station/PDG, entrance u,w and signed angles,
  physical nearest anode, total collected charge, total ion births, number
  of collecting anodes, leading/second wire ratios, nearest-wire fraction;
- `wires`: each occupied physical muon×wire, actual collected electron count,
  exact ion-birth count, global/local index and ion hemisphere asymmetry;
- `seeds`: primary-level collected charge and born ions and **its own**
  second-wire fraction, not the corresponding muon's fraction;
- `seed_wires`: each primary×physical wire counts, with primary total and
  station-aware wire association.
- Both `wires` and `seed_wires` preserve the exact
  native mean and RMS of ion-birth global w and t, and collected-electron
  arrival t (mm/ns). Physically empty populations have NaN moments,
  not fabricated zero times or positions.

Per-channel (Jpsi/Phi/Omega), per-channel×station (18 combinations) and
all-channel histograms:
- `log_muon_Q`, `log_muon_ions`;
- `log_wire_Q`, `log_wire_ions`, positive populations only;
- `n_wire` including zero collected charge;
- `log_seedQ`, positive-primary populations only;
- `second_wire_fraction` where `Q2/Qtotal` uses ranked **actual collected
  endpoint counts**, not necessarily adjacent wires;
- `nearest_wire_fraction` versus original entrance-nearest anode;
- `birth_hemi`, an exact per-wire ion-birth asymmetry;
- mean probability `Q2/Qtotal > 0.05` vs absolute impact phase
  and absolute signed-incidence angle;
- `q_fraction_by_wire` 2D: **mean** `Q_wire/Q_muon` versus signed impact
  phase and relative anode index, including absent wire charges as zeros;
- `ion_fraction_by_wire` 2D: independent exact nearest-birth wire
  fraction versus signed phase.

`channel_station_summary.csv` includes sample counts, zero/single/two/3+
collecting anode muons, second-wire >5% fraction, its naive binomial SE,
and total ion births/collected electrons by channel and station.
Plots are ROOT-made PNGs with honest `NA60+ Simulation` labels.

**No cathode strips, strip weighting field, VMM shaping, reconstruction
efficiency or detector threshold is simulated by this analysis.**
Ion birth assignment to wire is not ion drift endpoint assignment.

## Added per-anode gain vs ion-birth asymmetry (same ROOT analysis)

The exact native Garfield++ wire records contain **two independently
defined ion-birth angular asymmetries**. Both are computed for positive-ion
creation points assigned to the nearest physical anode by birth position,
**not** from collected-electron endpoints:

- `A_ion_hemi = (N_up - N_down)/(N_up + N_down)` for ion birth azimuth
  phi in (0,180) versus (-180,0). Exact boundary births are excluded.
- `A_ion_TB = (N_top-N_bottom)/(N_top+N_bottom)` for phi in [45,135)
  versus [-135,-45). Ion births outside those angular sectors are excluded.
  When the relevant denominator is zero, **A is NaN/undefined**, not 0.

The wire-level gain proxy is **`G_wire=Q_collected_wire/N_Heed,muon`**;
the muon gain proxy is `G_muon=Q_collected_muon/N_Heed,muon`.
For a single original Heed primary, its intrinsic wire contribution is
`G_primary_wire=Q_collected_primary_wire`. These are distinct
normalizations: a wire with zero collected electrons can still have
true ion births and asymmetry; log10(gain) is undefined for zero gain.

The same analysis now creates (for all pooled data, Jpsi, Phi, Omega;
all 18 channel×station groups are available as ROOT distributions):

- `asym_gain_vs_A_hemi_wire_all.png`: 2D density, one row per ion-bearing
  physical anode with positive collected charge, versus log10(G_wire).
- `asym_gain_vs_A_hemi_muon_all.png`: independent muon-level
  gain vs muon's exact aggregate ion-birth hemisphere A.
- `asym_gain_vs_A_hemi_primary_wire_all.png`: one Heed-primary×anode
  contribution, preserving native ion birth sectors.
- `asym_gain_vs_A_TB_wire_all.png`: the restricted 45°–135° top/bottom
  sector asymmetry version.
- `asym_muon_sharing_vs_A_hemi_all.png`: independent muon's Q2/Qtotal
  versus aggregate A_hemi.
- `asym_ion_births_vs_A_hemi_wire_all.png`: exact N_ion vs A irrespective
  of collected-electron endpoints.
- `asym_charge_fraction_vs_A_hemi_wire_all.png`: actual Q_wire/Q_muon
  including zero-collected but ion-bearing wires.
- `asym_A_hemi_vs_A_TB_all.png`: relationship of the two distinct azimuthal
  sectors. Do not assume they are independent.
- `asym_A_hemi_vs_local_wire_all.png`: local anode offsets -4..+4,
  centered on the physical entrance-nearest anode (station0 half-pitch).
- `asym_A_hemi_vs_signed_impact_all.png`,
  `asym_A_hemi_vs_theta_u_all.png`,
  `asym_A_hemi_vs_theta_w_all.png`: signed geometry dependence.
- `asym_zero_collected_ion_wires_A_hemi_all.png`: exactly the group with
  ion births but **Q_wire=0**. This group is intentionally absent from
  the logarithmic gain panel but remains in `wires` and QA denominators.
- `asym_mean_A_vs_local_wire_all.png` and signed impact/angle profiles:
  mean A conditioned on the explanatory variable, **unweighted per
  ion-bearing physical anode**. These profile error bars are exploratory:
  multiple wires from one muon are correlated and do NOT represent
  a full event-cluster bootstrap uncertainty.

For each channel and station, extra
`asym_gain_vs_A_hemi_wire_jpsi_station0.png` etc. and
`asym_A_hemi_vs_local_wire_jpsi_station0.png` are exported.

`asymmetry_denominators.csv` reports precisely how many ion-bearing
wire/primary-wire rows, defined A values, zero-electron cases and
positive-gain/defined-A pairs contribute in each group. The `wires`,
`seeds`, `seed_wires` and `muons` ROOT TTrees preserve the new
`ion_A_TB`, true sector counts and unlogged gain values.

No cathode strip current, gas ion drift velocity or physical
charge-sign response is inferred by these plots.

## WIPP commands (new isolated worktree, never change ongoing jobs)

Connect (no campus VPN if direct WIPP access works):

```bash
ssh -o ServerAliveInterval=60 -o ServerAliveCountMax=3 \
  stanis@wipp-an1.weizmann.ac.il
```

On `wipp-an1`:

```bash
git clone --depth 1 \
 --branch feature/p3-1800v-wiretruth-8k-analysis-20261123 \
 https://github.com/StasBren/NA6PRoot.git \
 "$HOME/NA6PRoot-wiretruth1800-analysis"

cd "$HOME/NA6PRoot-wiretruth1800-analysis"

RUNNER=test/mwpc-stage-c/wiretruth1800/run_wiretruth1800_analysis.sh
WORKER=test/mwpc-stage-c/wiretruth1800/pbs_wiretruth1800_analysis.sh

git log -1 --oneline
bash -n "$RUNNER"
bash -n "$WORKER"
bash "$RUNNER" build

# Set this ONCE in this shell, re-use it to retrieve outputs.
export MWPC_1800_ANALYSIS_OUT="$PWD/test_runs/wiretruth1800_snapshot_$(date -u +%Y%m%dT%H%M%SZ)"
echo "RESULT_DIR=$MWPC_1800_ANALYSIS_OUT"

# Smoke test: 1 SHA-verified chunk x 20 genuine Geant4 muons per
# physics channel (60 muons); no Garfield/Heed reruns.
bash "$RUNNER" pilot
```

First check the `w1800pil` PBS job (`N`, 1 CPU, 4 GB,
1 hour). Only after its `Exit_status=0` and the produced 60-muon
`asymmetry_denominators.csv` and PNGs pass checks, submit the **full
available ~8k snapshot** into a DIFFERENT unique folder:

```bash
export MWPC_1800_ANALYSIS_OUT="$PWD/test_runs/wiretruth1800_full_$(date -u +%Y%m%dT%H%M%SZ)"
unset MWPC_1800_CHUNKS_PER_CHANNEL
bash "$RUNNER" submit
```

PBS queue: `N` (used by the earlier MWPC jobs), configurable with
`MWPC_1800_PBS_QUEUE` if your account has a different permitted queue.
If PBS denies access, check `/opt/pbs/bin/qstat -Qf N` before trying again;
never alter production jobs or use the WIPP login node for full analysis.

The full PBS job `w1800ana` requests one CPU/4GB/4h; the unbounded
snapshot includes all complete SHA-verified original chunks at that time.
No held/ongoing physics jobs are altered.

PBS job `w1800ana`, `N`, one CPU/4GB/4h, no Garfield++ physics.
If the queue rejects the resource request, inspect `qstat -Qf`
before choosing a legal queue. For any SHA mismatch, source inconsistency
or job wall-time failure, the job returns nonzero; do not claim results.

```bash
/opt/pbs/bin/qstat -u "$USER"
tail -n 5 "$PWD/test_runs/pbs-logs/submissions.tsv"

# Once the job has Exit_status=0:
cat "$MWPC_1800_ANALYSIS_OUT/channel_station_summary.csv"
cat "$MWPC_1800_ANALYSIS_OUT/analysis_readme.txt"

# Check all generated products.
(cd "$MWPC_1800_ANALYSIS_OUT" && sha256sum -c SHA256SUMS)

# PDFs are NOT needed. PNG figures and ROOT TTree are produced.
ls -lh "$MWPC_1800_ANALYSIS_OUT/"*.png
ls -lh "$MWPC_1800_ANALYSIS_OUT/wiretruth1800_8k.root"
```

To inspect exact PBS exit status use
`/opt/pbs/bin/qstat -x -f <jobid> | grep -E 'job_state|Exit_status|resources_used'`.
To inspect stdout/stderr, list `$PWD/test_runs/pbs-logs` and use the real
PBS job ID (avoid placeholder paths). If file I/O proves too slow under
PBS, preserve logs and benchmark the source path before extending walltime
or weakening the SHA checks.

For the first 60-muon pilot, send its `PASS ... real muons=60` log,
`asymmetry_denominators.csv` and
`asym_gain_vs_A_hemi_wire_all.png` /
`asym_A_hemi_vs_local_wire_all.png`. If the smoke test fails, do NOT
submit the full ~8k array; return the exact error to ChatGPT for a code fix.

**Files to send back for physical interpretation**: job stdout with final
`PASS ... muons=` counts, `channel_station_summary.csv`, the
`all_per_wire_fraction_collected_electrons.png`,
`all_per_wire_fraction_ion_births.png`,
`hSecond.png`, `share_probability_vs_phase.png`,
`share_probability_vs_theta.png`, and `hLogQWire.png`.
Report each channel's number of independent muons and do not silently
normalize per-wire and per-muon histograms as if they sampled equal units.

## Comparison to 2100 V next

Only after seeing the real 1800 V per-wire results, use a separate
2100 V population with **matching source crossing IDs/geometry** where
possible. Do not compare an 1800 V 3mm-pitch simulation to an as-built 4mm
wire or different cathode geometry as if only the voltage changed.
For a voltage scan, use conditional gain/charge distributions,
shared-wire probabilities and statistical uncertainty at fixed
incident-angle/impact-position/Heed conditions; additionally preserve
zero-charge events. No 2100V resimulation is authorized here.

## Claude handoff

See `CLAUDE_HANDOFF_1800V_PER_WIRE.md`. Ask Claude to review the
ROOT-native distributions and sampling/normalization only after a WIPP
compile succeeds or when a specific physics discrepancy is identified;
do not run competing unverified physics simulations.
