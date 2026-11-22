# 1800 V vs 2100 V consistency audit (repository-based, 2026-11-22)

Sources read: `PI_500_GARFIELD_P3_README.md`, `PI_AVALANCHE_LUT20_LUT500_REUSE_DECISION.md`,
`PI_FAST_LUT_TRANSPORT_ONLY_README.md`, `c_transport_box_map.cxx`, `pi_audit_avalanche_lut_reuse.py`,
`P3_WIRETRUTH_1000_README.md` (branch `feature/p3-wiretruth-30k-direct-full-20261114`); the 2100 V Garfield source
`phaseB_native_muon_p3_noFW_penning057_cathode_ion_reference_v1.cxx`, `Plane2100Response.hh`,
`p3_export_full_spectral_paths.cc`, `COMPACT_KERNEL_PILOT.md`. The Fast4 and Fast5-T sources themselves are not in the
repository; their dependencies below come from the project documents.

## What each existing quantity was computed at

| Quantity | Source | Voltage | Used by |
|---|---|---|---|
| Electron collection probability, collecting wire, drift time, end-w, hit azimuth `phi_hit`, incidence angle | LUT500: `c_transport_box_map.cxx` (`DriftElectron`, attachment on, 0.1 eV start), 30 x 30 cells (0.1 x 0.2 mm) x 500 trials, raw `garfield_electrons.csv` | 1800 V (`--hv` default 1800, frozen config) | Fast4 primary-level transport; LUT500 asymmetry dashboards; AvalancheMC-vs-microscopic study |
| Gain, ion counts, hemisphere A, cloud extents, ion-creation delay (mean 27.09 ns, RMS 14.95 ns) | 30k-muon microscopic campaign | 1800 V | Fast4 avalanche attachment, Fast5-T donor, slide 9 and 10 numbers |
| Wire-truth validation (19.2% vs 19.66% second anode > 5%) | `P3_WIRETRUTH_1000` | 1800 V | Fast4 / Fast5-T validation |
| Ion births, exact per-wire census, electron collection times per wire, per-primary electron losses | K256 archives + `microscopic_wire_seed.csv` + `microscopic_seed.csv` | 2100 V | V0, ion-response pilots |
| Ion drift field and whole-plane / U/V weighting | `Plane2100Response.hh`, path exporter, spectral LUT | 2100 V field; Ar+ in Ar mobility (provisional, field-independent table) | 2100 V ion-response chain |

The 2100 V ion-response chain reads no 1800 V input. Fast4 and Fast5-T are 1800 V models end to end.

## What a 2100 V statistical model would need if it kept an explicit electron-transport stage

All five LUT500 outputs are field dependent: collection probability (attachment competes with drift time), drift time,
longitudinal diffusion (`w_end - w_start`), collecting wire near cell boundaries, and the hit azimuth that sets
the avalanche orientation. None can be transferred to 2100 V without a 2100 V drift calibration.

## Can E3 skip the transport stage?

Yes, in principle, and V0 now tests it. Each 2100 V archive records, for every real Heed primary, its start position
and its outcome: collected electrons per wire, ion births per wire, the mean endpoint time of the collected avalanche
electrons, ion-birth time moments, mean ion w, and the ion-cloud moments. Conditioning E3 directly on start position
marginalizes the electron drift, so LUT500 is not needed. `v0_summary.txt` section [8] and `v0_transport_map.csv`
give these 2100 V conditional statistics on a 6 x 12 grid (|du| 0.25 mm x v 0.5 mm), with counts per cell.

Two things the archives do not contain: individual electron drift lines, and the electron arrival azimuth at the wire.
E3 would therefore have to sample the ion-cloud moments directly given start position, not derive them from an
arrival azimuth. Whether that conditional is adequate is an empirical V0 question.

## Decision rule for a new 2100 V drift calibration

Rebuild transport only if one of these holds after the real V0 run:
1. a cell that matters (where real muons put primaries) has too few 2100 V primaries for the statistics E3 keeps;
   then run a targeted 2100 V microscopic sample in those cells only, not all 450,000 trajectories;
2. E3 must be parametric in voltage or geometry for design scans; then transport alone is insufficient anyway,
   because gain and cloud shape change too, and full avalanche samples are needed at each design point.

## Comparison columns

With `--lut500 .../p3p0_anode_only/garfield_electrons.csv` V0 bins LUT500 identically and prints it next to the
2100 V numbers. The observables differ (single-electron `DriftElectron` endpoint time at 1800 V vs mean endpoint time
of the collected avalanche electrons at 2100 V); the comparison shows the size of the voltage effect, it does not
calibrate one from the other.
