# Stage C reference results

This directory is for **small, reviewed, reproducible summaries** that we want
to keep with the Stage-C code.

Do not commit raw Geant4 ROOT files or the full high-statistics trajectory
production here.  Those persistent inputs live in the separate
`na60-dice-mwpc-simulation` study tree and are read in place.

Expected committed products are things such as:

- a final C0 `incident_angle_summary.csv`;
- a compact metadata JSON identifying the input suite and cuts;
- later, small Stage-C response/LUT validation summaries.

Runtime output from exploratory runs belongs under:

```text
test_runs/mwpc_stage_c/
```

which is already ignored by the repository.
