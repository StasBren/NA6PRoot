# Exact per-anode Garfield++ statistics — 2100 V archival data

The existing 2100 V Plot Studio now includes **Individual anodes**, with two distinct observation levels:

1. **Primary × physical anode**: one original `microscopic_wire_seed.csv` row per Heed electron and anode. Exact up/down/boundary positive-ion births, mean ion-birth time, time RMS, mean electron collection time, electron RMS and longitudinal ion spread.
2. **Muon × physical anode**: exact sum of ion-up/down/boundary and collecting-electron counts across contributing primaries in one archived muon, plus appropriately weighted time moments.

The original 997-muon source is NOT simulated again. Each saved `population_*/rootfiles.txt` or `pilot_*/rootfiles.txt` contains immutable ROOT paths; the dashboard reads the adjacent original `microscopic_wire_seed.csv` files directly. ROOT summary `v0_primaries.csv` provides the exact primary timestamp and independently audited birth/electron counts. ROOT `v0_muons.csv` provides exact muon totals and detector identity.

The per-anode input loader **aborts rather than drops data** if source exact per-primary or per-muon ion/electron sums disagree with the ROOT census, if up+down+boundary != ion births, or if the recorded A differs from recomputed A.

## Timing: three different quantities

For primary p on wire a, the exact **mean ion formation delay** is

```
tau_ion[p,a] = mean_t_birth_ns[p,a] - t_primary_ns[p]
```

This is the mean ion birth time after a Heed electron's creation. It is **not** first ion-birth/onset time, not electron collection time and not ion drift to the cathode.

The exact **within-avalanche birth-time spread**, `sigma_t_birth_ns`, is the RMS width of the ion production time distribution for a single (p,a) contribution.

The exact **mean anode electron collection delay** is

```
tau_e[p,a] = mean_t_collected_ns[p,a] - t_primary_ns[p]
```

defined only for anodes with collected avalanche electrons. The collected-electron signal is an anode endpoint count, **not** induced cathode-strip charge.

For muon m and wire a:

```
B = sum_p N_ion[p,a]

tau_ion[m,a] = sum_p(N_ion[p,a] * tau_ion[p,a]) / B

sigma_tau_ion[m,a]^2
    = sum_p[N_ion[p,a] * (sigma_t_birth[p,a]^2 + tau_ion[p,a]^2)] / B
      - tau_ion[m,a]^2
```

The group RMS includes within-avalanche timing fluctuations **and** differences between per-primary average delays; each ion is always referred to its own Heed primary timestamp.

## Hemisphere convention and physical-anode assignment

```
phi = atan2(v, u - u_nearest_anode)
up = births with v > 0 ; down = births with v < 0
A_UD = (N_up - N_down)/(N_up+N_down)
```

Boundary births v=0 do not enter A's denominator but remain in total ion births. For muon–anode aggregation, compute A from summed original up/down counts, **not** from `mean(A_primary)`.

The `global_wire` label identifies the wire nearest each **ion birth** for the ion counts; collected-electron counts refer to their **actual collecting wire** with the same global wire identifier. These counts can differ. A wire may collect electrons for primaries whose ions were born nearest another wire.

Use the dashboard **station and detector ID filters** when comparing individual global-wire IDs; do not silently pool disparate detector geometries.

## In the dashboard

Select `Individual anodes`, optionally restrict station, detector and global wire. Available figures include:
- Exact A_UD distributions for each (primary,wire) and (muon,wire)
- Primary-to-ion mean formation delay and within-avalanche sigma_t
- Primary-to-**collected electron** mean arrival delay and 2D comparison with ion birth delay
- Ion formation delay versus exact ion gain, |A_UD| versus gain, sigma_w versus gain
- Per-muon anode ion yield, collected-electron count, shares of parent muon charge, and mean formation delay
- Muon–anode A_UD versus ion yield and versus mean delay
- A filterable **per-anode summary table** with exact birth/electron totals, event counts, mean A and timing

The summary table computes equal-weight means over (muon,anode) records; its exact total counts are not equal-weight estimates. Plot uncertainty bars come from resampling independent muons. The user's display controls—axis color (black), editable titles and labels, legends with mean/RMS, logarithmic scales, exports—remain available.

## Commands

In WSL, from the existing worktree:

```bash
cd ~/na6p/wt-gain-asym-2100
git fetch origin feature/2100v-997-gain-hemisphere-root-20261124
git switch --detach FETCH_HEAD
DASH=test/mwpc-stage-c/cathode-ion-reference/gain-asym-2100/run_plot_studio_2100.sh
bash "$DASH" check
bash "$DASH" run
```

The `check` action now runs **four synthetic math/closure tests** before launching the dashboard. They exercise exact counts, aggregation, boundary births, and a deliberate inconsistent census. These are unit tests, not scientific results.

If `population_*` is not present, the latest already-generated `pilot_*` dataset is displayed with a preliminary warning. The existing C++/ROOT population analysis (read-only) is the only prerequisite to full 997-muon dashboard statistics. **No new Garfield simulation is required.**

### Deliberately not claimed

The per-anode plots do **not** show time-resolved ion drift, ion collection at cathodes, first avalanche ionization time, or actual induced U/V cathode-strip currents. They summarize microscopic ion-production and anode-electron-collection data only.

Detector model: 2100 V idealized P3 no-field-wire 9-anode reference, 3 mm anode pitch, 30 μm anode diameter, 3 mm each side to whole cathode plane, Ar/CO2 70:30 with Penning Ar r=0.57, B=0. This is not automatically the as-built prototype.
