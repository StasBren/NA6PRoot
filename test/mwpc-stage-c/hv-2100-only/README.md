# NA60+ MWPC — 2100 V working-point verification

User decision (2026-11-17): focus on **2100 V** at anode diameter **30 um**.
The user reports an approximate gain of **80,000**, which must be treated as a
preliminary number until the actual ROOT mean, fluctuation and cap audit pass.

This is a dedicated, isolated follow-up; it **does not** simulate 2150, 2200 or
2250 V. The earlier high-voltage scan was exploratory and remains unmodified.

- Full Garfield++ `AvalancheMicroscopic` with Ar/CO2 70:30,
  Penning probability r=0.57, P3 noFW, 3 mm pitch, 6 mm cathode gap, B=0.
- Identical frozen Geant4 crossings/TrackHeed input selection as the 1000-muon
  physical collecting-wire reference.
- Actual collecting-wire and nearest-wire ion birth counts saved separately.
- Avalanche size cap 1,000,000 by default, and runner fails on any cap.
- No space charge, streamer, discharge or photon-feedback validation.
- Signed muon ion-birth A is not induced cathode-strip charge.

## Commands

```bash
git clone --depth 1 --branch feature/p3-2100V-only-20261117 \
 https://github.com/StasBren/NA6PRoot.git ~/na6p/p3-2100-only

# smoke test first: 1 muon at 2100 V
bash ~/na6p/p3-2100-only/test/mwpc-stage-c/hv-2100-only/run_local_2100.sh 1

# only if cap/runtime check passes:
bash ~/na6p/p3-2100-only/test/mwpc-stage-c/hv-2100-only/run_local_2100.sh 5
bash ~/na6p/p3-2100-only/test/mwpc-stage-c/hv-2100-only/run_local_2100.sh 20
```

Previously sealed one-muon chunks are checksum verified before reuse. A changed
binary/config or a capped avalanche stops the run. Do not run concurrently with
the old scan if laptop CPU/memory is constrained.

Results:
`~/na6p/p3-2100-only/test_runs/mwpc_stage_c/p3_2100V_only_v1/analysis/`

ROOT/C++ produces `2100V_gain_statistics.csv`, `positive_gain_2100V.pdf/png`,
and `signed_asymmetry_2100V.pdf/png`. Reports include
mean per-primary gain (zeros included), gain conditional on positive collection,
RMS/mean at per-avalanche and per-muon levels, and bootstrap errors that resample
independent muons.

**Important:** 1 or 5 muons cannot establish precise fluctuation or asymmetry
statistics. Even 20 is only a pilot. Voltage is the *user-selected development
point*, not a verified stable detector operating point. Keep it configurable
for subsequent cathode response implementation.

Status: files created and source inspected; actual compilation and execution
must be confirmed from the user's WSL logs.
