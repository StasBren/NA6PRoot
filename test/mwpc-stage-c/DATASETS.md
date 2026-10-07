# Stage C input datasets

Stage C intentionally separates **persistent transport inputs** from the
microscopic-response code.

## High-statistics MNP33-bore production

Current reference suite:

```text
$HOME/na60-dice-mwpc-simulation/geometry/station-layout-study/full-setup/
  mnp33-bore-runs/20260912-202112-0wrxhdg_/
```

This production was made in the separate
`StasBren/na60-dice-mwpc-simulation` study repository.

It contains 200,000 parent dimuon events per source for:

- J/psi;
- Phi;
- Omega.

Each source is split into two independent 100k-event Geant4 chunks.

Relevant persistent files in each chunk include:

- `MCKine.root`;
- `LocalTrajectoryAudit.root`;
- request/manifest metadata.

The audit records transported states at nominal MS planes and at the dedicated
MNP33 reference planes

```text
z = 365, 430, 495 cm.
```

## Why the raw production is not copied here

The ROOT production is large, already persistent, and belongs to the
transport/layout study.  Copying it into NA6PRoot would create two sources of
truth.

Stage-C analyses therefore read it in place.

## What belongs in this repository

Commit:

- analysis code;
- compact configuration;
- small reviewed summaries and metadata needed to reproduce conclusions.

Do not commit:

- full ROOT chunks;
- large transient plots/CSV dumps;
- regenerated copies of the same transport production.

Runtime output belongs under:

```text
test_runs/mwpc_stage_c/
```

and small reviewed results belong under:

```text
test/mwpc-stage-c/results/
```

## C0 selections

The current C0 analysis distinguishes:

- `all`;
- `geometry_bore`;
- `useful_bore`;
- `useful_bore_working_area`.

The last selection is the one intended to answer whether the large downstream
angle tails still correspond to muons that remain inside a finite station
working area.
