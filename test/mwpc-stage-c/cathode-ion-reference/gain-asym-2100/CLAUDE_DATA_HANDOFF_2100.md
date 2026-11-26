# NA60+ 2100 V, 997 original Garfield muons — handoff to Claude

Source confirmed on WSL:
`$HOME/na6p/cathode-overnight-2100/test_runs/mwpc_stage_c/overnight_2100_weighted_k256_v1/`
and verified through `run_gain_asym_2100.sh locate`: 997 sealed+ROOT-audited, three other slots.

This packager is **read-only** and reuses the existing strict `per_anode_census_2100.py` checker. It does not rerun any microscopic simulation or ion drift and does not touch 1800 V or the user's original archives. It needs a completed
`test_runs/mwpc_stage_c/2100_gain_asym/population_*/` from `run_gain_asym_2100.sh population`, with `rootfiles.txt` and the verified V0 CSVs.

## On the user's WSL computer

From `~/na6p/wt-gain-asym-2100` after fetching this branch:

```bash
PACK=test/mwpc-stage-c/cathode-ion-reference/gain-asym-2100/claude_handoff_2100.sh
bash "$PACK" locate
bash "$PACK" check
bash "$PACK" export
```

The resulting file is **`~/NA60_2100V_997_for_Claude.zip`** (change with `MWPC_CLAUDE_2100_ZIP`; no overwrites).
Open the folder in Windows via `explorer.exe "$(wslpath -w "$HOME")"`.

Inside the ZIP:

| File | Statistical level / meaning |
| --- | --- |
| `muons.csv` | One Geant4 muon; exact ion and collected electron totals, angles, station and source |
| `wires.csv` | One muon × global anode; exact ion up/down/boundary counts, collected endpoints, pooled birth and collection timing |
| `seeds.csv` | One Heed primary; exact ion counts, zero-ion events, start timing and position offsets |
| `seed_wires.csv` | One primary × anode; exact positive-ion counts, exact hemispheric A, electron endpoint counts, formation and collection timing, exact longitudinal RMS |
| `k256_fragment_moments.csv` | Sampled weighted azimuthal moments and K256 sampling uncertainties per primary/anode |
| `v0_wire_moments.csv` | Existing aggregate K256 moment exports |
| `original_avalanche_status.csv` | Original Garfield successful/failed/capped and electron loss counters, initial seed coordinates |
| `original_heed_primaries.csv` | Original Heed primary records |
| `original_geant4_inputs.csv` | Original individual input crossing metadata |
| `README_for_Claude.txt`, `manifest.json` | Conventions, QA, provenance and precise limitations |

**Not a drop-in 1800 V ROOT clone.** It is a compact set of CSV **equivalents** of the four conceptual ROOT levels. The original 1800 `wiretruth1800_full.root` has different branch names and some source metadata. Ask Claude to map by column name; do not rename fields blindly or assume the 1800 schema. The 2100 V ROOT archive contains extra primary start/ion-birth information absent from the older 1800 V export. We can later make an exact cross-compatible ROOT exporter after inspecting the two schemas on the same computer.

Asymmetry `A_UD=(N_up-N_down)/(N_up+N_down)`: positive v vs negative v around the anode; boundary ion births excluded from A denominator only. Ion births use **nearest birth anode**, collected electron counts use the **true collecting anode**. The latter is not induced cathode strip signal or ion collection. `G_ion` and `G_collected` must use Heed-primary denominators, with global means computed as ratios of sums. K256 angular moments are statistically sampled; exact birth/collection counts and selected longitudinal/timing moments are from frozen Garfield counters. The sample includes prior inspected/validation events and is intended for descriptive distributions only.

Idealized numerical 2100 V reference: P3 no field wires; nine 30 μm anodes, 3 mm pitch, cathodes at ±3 mm; Ar/CO2 70:30, Penning 0.57, B=0. Not the as-built Prototype-3. Not mixed with the 1800 V scan, or the separate planned 3×1000 matched-muon campaign.
