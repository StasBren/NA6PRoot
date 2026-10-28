# Parent-event asymmetry: electron counters and p3/p4 LUT selector

This update addresses the screenshot labelled **99,968 Parents with crossing**.
That screenshot uses ONE J/ψ PARENT EVENT per station, not one muon gas
crossing. Never insert crossing-level P(A) from `full34` directly into
that parent-level graph.

The existing p3/LUT500 ROOT-style HTML template is extended with:
- a geometry selector: **3.0 mm · LUT500 (original)**,
  **3.0 mm · LUT20**, **4.0 mm · LUT20**;
- two electron cards: classified anode endpoints noFW and FW75;
- an explicit breakdown `N_up` and `N_down` for each geometry;
- dynamic p/LUT labels, 41-bin P(A), RMS, parent counts,
  DET map and chamber table.
The *electron* counters are `total_up+total_down`, summed across all
selected parent events. Exact azimuthal boundary hits are excluded,
and the counters measure primary-electron Garfield LUT endpoint positions,
not avalanche multiplication or cathode induced charge.

The original p3 LUT500 summary already has `total_up/total_down` per
station/DET. For p3 and p4 LUT20, the completed 100k-parent
`pi_lut20_p4_muon_crossing_asymmetry_v1/full34/crossings.csv.gz`
is regrouped by `(chunk,event)`. For each parent, **all its muon
gas crossings inside the requested station (or DET)** contribute
their upper/lower electron counts BEFORE one value of A is computed.
Do NOT average crossing A values or sum crossing histogram bins.

## Run on WSL

```bash
cd ~/na6p/src/NA6PRoot-mwpc
git switch feature/mwpc-stage-c
git pull --ff-only

bash test/mwpc-stage-c/run_pi_parent_asymmetry_geometry_dashboard.sh check
bash test/mwpc-stage-c/run_pi_parent_asymmetry_geometry_dashboard.sh status
bash test/mwpc-stage-c/run_pi_parent_asymmetry_geometry_dashboard.sh build
```

Open new **combined** dashboard:

`test_runs/mwpc_stage_c/pi_lut500_parent_asymmetry_root_v1/dashboard_geometry_lut_compare.html`

Do NOT expect new options to appear in a cached old
`pi_lut500_parent_asymmetry_root_v1/dashboard.html`. The older file
remains unchanged. `build` only writes the separate combined HTML
and its provenance sidecar.

Physical source files are read-only. Before regrouping, code verifies
both original manifest SHA256 sources. It also checks that:
- input is exactly 100k generated J/ψ parent slots and all 6 MS stations;
- p3 LUT500 and p3/p4 LUT20 have matching original Geant4 crossings
  and Heed primaries per station;
- LUT20 CSV contains each expected configuration exactly once per
  physical muon gas crossing, with conserved electron counts;
- sum of parentized electrons equals frozen LUT20 crossing summary;
- number of parent-with-crossing observations agrees with original
  p3 LUT500 parent result;
- sum of physical DET electron counts equals station electron count.

The standalone synthetic regression includes **two crossings belonging
to the SAME parent**, one with A=+1 and one with A=-1, verifying the
station's parent A=0 without substituting crossing moments.

The full34 reader streams from gzip, only retaining the current parent
in memory. No Garfield/Heed/Geant4 calls and no source recomputation.

If the combined HTML exists and styling must be regenerated:

```bash
bash test/mwpc-stage-c/run_pi_parent_asymmetry_geometry_dashboard.sh replace-html
```

Only that combined HTML is replaced; no physics or original dashboard
files are overwritten.

## Important interpretation

The selector intentionally shows different Garfield precision:
p3/LUT500 vs p3/LUT20 and p4/LUT20. Displayed curves in any given
plot ALWAYS use the same pitch and LUT trial count for FW/noFW,
never an overlay of LUT20 against LUT500 as though they were
equivalent errors.

A separate muon crossing P(A) dashboard remains available under
`pi_lut20_p4_muon_crossing_asymmetry_v1/full34/crossing_dashboard.html`
and, when complete, the LUT500 crossing-level dataset. Parent-level
and crossing-level histograms answer different questions.

All plots follow the agreed ATLAS-inspired presentation conventions
without ATLAS branding: white canvas, black axes and ticks, blue solid
noFW, red dashed FW75, `NA60+ Simulation — Work in progress`.
