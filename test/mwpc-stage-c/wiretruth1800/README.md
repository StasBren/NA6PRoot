# 1800 V real-Geant4 ~8,000-muon wiretruth distribution analysis

**Status 2026-11-23:** New C++/ROOT analysis code + file-integrity snapshot
script are committed in an isolated branch. Not yet compiled/tested on the
WIPP runtime, and no analysis result can be claimed before actual PBS output.
Zero new Garfield++ muons are simulated.

## Scientific source and identity

Only existing **completed and SHA-verified** production chunk results from

\`/storage/cgrp/stanis/mwpc_penning057_wiretruth30k_v1/results/{jpsi,phi,omega}/chunk_XXXX\`

are included. Original independent source crossing bytes remain frozen at

\`/storage/cgrp/stanis/mwpc_penning057/original2m_{channel}_10000_v1/chunks/chunk_XXXX.csv\`.

Every recorded chunk contains exactly 20 original Geant4 crossings.
The 2026-11-16 status was around 131 Jpsi, 136 Phi and 50 Omega chunks;
the Omega 50–134 PBS array was launched to finish about 135.
The exact number analyzed is **whatever successful, complete receipts exist
at the time of the SHA-frozen snapshot**. It is never assumed to be 8,000.

Original physics lock:
- production Git \`599e29522edb6fcec9dfcca2904308cb3186edce\`
- original binary SHA \`2cf92cb27d3a180cd26d3ab5417f57293473864f6d685ee43bacf931712aa8f1\`
- 1800 V, 9-wire analytic P3 noFW, 3 mm physical wire pitch,
  30 μm diameter, cathodes at v=±3 mm, Ar/CO₂ (70/30), Penning r=0.57
- station zero physical anode phase = +1.5 mm;
  all other stations physical anode phase = 0 mm.
- These numerical values do **not** describe the final as-built MWPC
  chamber geometry; distinguish them in all reports and future parameter scans.

## What gets analyzed

Full, **exact**, native Garfield++ counts (not ion-reservoir estimators):

- \`microscopic_crossing.csv\`: independent original Geant4 muon totals;
- \`microscopic_seed.csv\`: actual Heed-primary Garfield electron/ion census;
- \`microscopic_wire_seed.csv\`: individual primary × actual physical anode,
  exact collected-electron endpoints and nearest-ion-birth assignments;
- \`microscopic_wire_muon.csv\`: complete muon × actual physical anode
  counts, ion-hemisphere moments and station-aware anode indices.

The Python section of the Bash script is restricted to **file integrity,
selection and SHA256**, not physics/statistical analysis. Before adding
a chunk, it compares the frozen source input hash and the four required
CSV hash+size values to its \`verified_report.json\` and immutable production
selection/binary lock. Only clean completed records enter the manifest.

C++/ROOT independently tests:
1. each chunk has 20 original independent Geant4 muons;
2. physical station, nearest wire, global wire centre and patch \`-4..4\`;
3. every Garfield primary has \`successful=1\`, \`capped=0\`;
4. exact per-primary \`(Q_e,N_ion)\` = sum across its physical wires;
5. exact per-muon \`(Q_e,N_ion)\` = sum across primaries = sum across physical
   wires = source Garfield crossing;
6. positive-ion up+down+boundary = exact total ion births on wire;
7. positive-ion count on nearest birth wire is distinct from collected
   electron endpoint count on actual collecting wire.
   
**All original crossings are included as independent observations**, even
if a muon has zero collected charge. Events where an anode has ion births but
zero collected electrons are also retained.

## Physics distributions

The ROOT file has four TTrees:
- \`muons\`: Geant4 channel, det/station/PDG, entrance u,w and signed angles,
  physical nearest anode, total collected charge, total ion births, number
  of collecting anodes, leading/second wire ratios, nearest-wire fraction;
- \`wires\`: each occupied physical muon×wire, actual collected electron count,
  exact ion-birth count, global/local index and ion hemisphere asymmetry;
- \`seeds\`: primary-level collected charge and born ions and **its own**
  second-wire fraction, not the corresponding muon's fraction;
- \`seed_wires\`: each primary×physical wire counts, with primary total and
  station-aware wire association.

Per-channel (Jpsi/Phi/Omega), per-channel×station (18 combinations) and
all-channel histograms:
- \`log_muon_Q\`, \`log_muon_ions\`;
- \`log_wire_Q\`, \`log_wire_ions\`, positive populations only;
- \`n_wire\` including zero collected charge;
- \`log_seedQ\`, positive-primary populations only;
- \`second_wire_fraction\` where \`Q2/Qtotal\` uses ranked **actual collected
  endpoint counts**, not necessarily adjacent wires;
- \`nearest_wire_fraction\` versus original entrance-nearest anode;
- \`birth_hemi\`, an exact per-wire ion-birth asymmetry;
- mean probability \`Q2/Qtotal > 0.05\` vs absolute impact phase
  and absolute signed-incidence angle;
- \`q_fraction_by_wire\` 2D: **mean** \`Q_wire/Q_muon\` versus signed impact
  phase and relative anode index, including absent wire charges as zeros;
- \`ion_fraction_by_wire\` 2D: independent exact nearest-birth wire
  fraction versus signed phase.

\`channel_station_summary.csv\` includes sample counts, zero/single/two/3+
collecting anode muons, second-wire >5% fraction, its naive binomial SE,
and total ion births/collected electrons by channel and station.
Plots are ROOT-made PNGs with honest \`NA60+ Simulation\` labels.

**No cathode strips, strip weighting field, VMM shaping, reconstruction
efficiency or detector threshold is simulated by this analysis.**
Ion birth assignment to wire is not ion drift endpoint assignment.

## WIPP commands (new isolated worktree, never change ongoing jobs)

Connect (no campus VPN if direct WIPP access works):

\`\`\`bash
ssh -o ServerAliveInterval=60 -o ServerAliveCountMax=3 \
  stanis@wipp-an1.weizmann.ac.il
\`\`\`

On \`wipp-an1\`:

\`\`\`bash
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

# ROOT-only analysis and file SHA audit in an independent PBS job.
bash "$RUNNER" submit
\`\`\`

PBS job \`w1800ana\`, \`normE\`, one CPU/4GB/4h, no Garfield++ physics.
If the queue rejects the resource request, inspect \`qstat -Qf\`
before choosing a legal queue. For any SHA mismatch, source inconsistency
or job wall-time failure, the job returns nonzero; do not claim results.

\`\`\`bash
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
\`\`\`

To inspect exact PBS exit status use
\`/opt/pbs/bin/qstat -x -f <jobid> | grep -E 'job_state|Exit_status|resources_used'\`.
To inspect stdout/stderr, list \`$PWD/test_runs/pbs-logs\` and use the real
PBS job ID (avoid placeholder paths). If file I/O proves too slow under
PBS, preserve logs and benchmark the source path before extending walltime
or weakening the SHA checks.

**Files to send back for physical interpretation**: job stdout with final
\`PASS ... muons=\` counts, \`channel_station_summary.csv\`, the
\`all_per_wire_fraction_collected_electrons.png\`,
\`all_per_wire_fraction_ion_births.png\`,
\`hSecond.png\`, \`share_probability_vs_phase.png\`,
\`share_probability_vs_theta.png\`, and \`hLogQWire.png\`.
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

See \`CLAUDE_HANDOFF_1800V_PER_WIRE.md\`. Ask Claude to review the
ROOT-native distributions and sampling/normalization only after a WIPP
compile succeeds or when a specific physics discrepancy is identified;
do not run competing unverified physics simulations.
