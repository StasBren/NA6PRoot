# Penning r=0.57 — interactive Plot Studio for 3 × 10,000 real muons

The existing 1000-J/ψ Plot Studio has been extended to select four immutable
datasets: **J/ψ (10,000)**, **φ (10,000)**, **ω (10,000)** and
**All (30,000)**. The three individual channels each have 500 original
Geant4 20-muon chunks. The pooled sample is a descriptive 1:1:1 muon
mixture, not an estimate of production-weighted dimuon resonance yields.

## Physics

3 mm wire pitch, no field wires, 1800 V anode, cathode planes ±3 mm,
Ar/CO₂ 70:30, Penning r=0.57, full TrackHeed + Garfield++
AvalancheMicroscopic. **No new physics calculation** is performed by the
dashboard. **No cathode signal or positive-ion drift is simulated**.

For every muon-level population, zero-retained-Heed events are counted but
have undefined muon gain. Avalanche gain includes exact collected-electron
counts, including zeros. Exact per-seed upper/lower ion-birth asymmetry is
computed from full native integer counters. Angular, spatial and birth-delay
plots use up to 128 sampled ion births per avalanche weighted by the exact
number of ion births. The sample time span is a *lower bound* on the full
ion-birth time span.

The original 12 plots plus the 5 aggregate temporal-correlation views and
4 aggregate anode-wrapping views are exposed in 30k mode. In addition,
two **real-event inspector** pages are available through the user's local
read-only SSHFS mount: (1) the original Geant4 muon path and all TrackHeed
clusters side by side with an exact per-cluster/per-primary avalanche-size
heatmap, and (2) one genuine microscopic Garfield ion cloud viewed in 2D,
rotatable 3D including longitudinal w, and angular distribution around the
nearest anode. The old 1000-muon native raw-file schema and names differ
from the original-2M 30k campaign, so the latter has its own strict source
adapter. The 30k individual-event pages do not load all 45 GB and are not
part of the aggregate multi-channel overlay.

## Cluster: Prepare (three PBS jobs)

```bash
ssh stanis@wipp-an1
cd ~/NA6PRoot-mwpc
git pull --ff-only
python3 -c "import numpy; print('numpy',numpy.__version__)"
python3 -m py_compile test/mwpc-stage-c/pi_penning057_30k_plotstudio_cache_v1.py
ROOT=/storage/cgrp/stanis/mwpc_penning057
mkdir -p "$ROOT/plotstudio_30k_v1/logs"
PBS="$PWD/test/mwpc-stage-c/pbs_penning057_30k_plotstudio_prepare_v1.sh"
J=$(qsub -v CHANNEL=Jpsi -o "$ROOT/plotstudio_30k_v1/logs/" -e "$ROOT/plotstudio_30k_v1/logs/" "$PBS")
P=$(qsub -v CHANNEL=Phi -o "$ROOT/plotstudio_30k_v1/logs/" -e "$ROOT/plotstudio_30k_v1/logs/" "$PBS")
O=$(qsub -v CHANNEL=Omega -o "$ROOT/plotstudio_30k_v1/logs/" -e "$ROOT/plotstudio_30k_v1/logs/" "$PBS")
printf 'Jpsi %s\nPhi %s\nOmega %s\n' "$J" "$P" "$O"
C="$PWD/test/mwpc-stage-c/pbs_penning057_30k_plotstudio_combine_v1.sh"
qsub -W "depend=afterok:$J:$P:$O" -o "$ROOT/plotstudio_30k_v1/logs/" -e "$ROOT/plotstudio_30k_v1/logs/" "$C"
```

PBS is used for large Lustre I/O and NPZ compression, not the login node.
Preparation is strictly non-overwriting and accepts a completed recovery
chunk only when the original chunk did not finish. Every native CSV/report
is SHA-256 checked before use. Only new
`/storage/cgrp/stanis/mwpc_penning057/plotstudio_30k_v1` files are written.

```bash
qstat -u stanis
grep -E '^PASS|^FINISH|Traceback|ERROR'  /storage/cgrp/stanis/mwpc_penning057/plotstudio_30k_v1/logs/*.log
du -sh /storage/cgrp/stanis/mwpc_penning057/plotstudio_30k_v1
```

## Local WSL: Open the interactive Plot Studio

Copy *only the compact derived caches* after PBS success, not the original
heavy ion-birth CSVs:

```bash
cd ~/na6p/src/NA6PRoot-mwpc
git pull --ff-only
mkdir -p test_runs/mwpc_stage_c
scp -r stanis@wipp-an1:/storage/cgrp/stanis/mwpc_penning057/plotstudio_30k_v1   test_runs/mwpc_stage_c/
export MWPC_PENNING_30K_STUDIO_ROOT="$PWD/test_runs/mwpc_stage_c/plotstudio_30k_v1"
python3 -m streamlit run   test/mwpc-stage-c/pi_p3_noFW_penning057_interactive_dashboard_v1.py
```

Select **Geant4 muon dataset** at the top of the sidebar. Settings for
histogram bin width, angular bins, axis scales, colours, legend text/position,
early-vs-late windows and camera/HTML/CSV exports are inherited from the
original Plot Studio. All four selectable datasets use the same bin definitions.

### Real single-muon / single-avalanche details through SSHFS

With the Weizmann F5 VPN connected, mount the original **read-only** archive
in local WSL and launch the existing dashboard:

```bash
mkdir -p ~/na6p/cluster_mwpc
if ! mountpoint -q ~/na6p/cluster_mwpc; then
  sshfs -o ro,reconnect,ServerAliveInterval=15 \
    wipp-an1:/storage/cgrp/stanis/mwpc_penning057 \
    ~/na6p/cluster_mwpc
fi
cd ~/na6p/src/NA6PRoot-mwpc
source test_runs/mwpc_stage_c/.venv_pi_slides/bin/activate
export MWPC_PENNING_30K_STUDIO_ROOT="$PWD/test_runs/mwpc_stage_c/plotstudio_30k_v1"
export MWPC_PENNING_30K_RAW_ROOT="$HOME/na6p/cluster_mwpc"
python -m streamlit run test/mwpc-stage-c/pi_p3_noFW_penning057_interactive_dashboard_v1.py
```

Turn **off** channel overlay and choose **Plot family → Detailed real
events · 2 pages**. The two pages are **Real muon · ALL clusters +
avalanche-size heatmap** and **One avalanche · 2D / 3D w / angular cloud**.
Each viewer selects a genuine `jpsi`, `phi` or `omega` event and local
channel muon index `0..9999`. The native event is mapped to the original
20-muon chunk, including the three successful repaired chunks. Raw files are
SHA-256 verified against `outputs.sha256` and original Geant4 input chunk
against its frozen manifest.

- In the first page, left plot shows actual muon entry/exit and individual
  Heed cluster centres in their original u,v geometry, including anode
  lattice/cathodes; right heatmap shows per-Heed-electron **exact** Garfield
  collected gain (or true ion births) within each cluster. Empty heatmap cells
  mean no accepted Heed electron, and zero values mean a true zero-gain
  microscopic avalanche.
- In the second page choose a real Heed electron/avalanche. The XY panel
  shows sampled ion-creation positions, an angular histogram uses
  `atan2(v,u-u_anode)`, and rotatable 3D shows the full sampled longitudinal
  w coordinate. A seed retains at most 128 uniform original ion birth
  coordinates; all exact ion/gain/sector counters are preserved.
- These are ion **births**, not ion drift or cathode pulses. Reading an
  individual large native reservoir through VPN may have some latency.

```bash
python test/mwpc-stage-c/pi_penning057_individual_views_30k_v1.py
```

### Overlay / comparison of 1–3 independent channels

After `git pull --ff-only`, leave the existing
`MWPC_PENNING_30K_STUDIO_ROOT` pointing at the four locally prepared
subdirectories and run the usual Streamlit dashboard (using its existing
`.venv_pi_slides` environment). Check **Overlay separate J/ψ, φ, ω
channels** at the very top of the sidebar, then select **any one, two or all
three channels** in **Channels on the same graph (1–3)**. The sidebar also
has a **Plot family** selector: original distributions (12), early-birth
correlations (5), and **Anode position & ion wrapping (4)**. Set the latter
to find the four later-added anode studies, both in single-dataset and
multi-channel comparison modes. All 21 are available under **All aggregate
plots** by default.

- 1D histograms and correlations: overlaid on ONE Plotly graph, with fixed
  colours (J/ψ orange, φ blue, ω green), channel-specific legends and
  single-file HTML/camera/CSV comparison exports.
- 2D heatmaps: 1–3 side-by-side panels in ONE figure with the same numerical
  colour scale; overlaying coloured images would conceal spatial structure.
- Per-channel probability or density curves normalise **independently**
  using shared physical histogram edges. Raw counts are available but
  represent different retained-electron/ion populations, not cross sections.
  Gain zero is a separate [0,1) bin, and zero-Heed muons do not contribute
  fictitious gain entries.
- Gain-versus-asymmetry uses common gain bins; bootstrapped 95% intervals
  for wire wrapping resample independent real parent muons **within each
  channel**, not the 1.55 million starting electrons. The early/late
  comparison's quantile groups are independently calculated by channel.
- The original **single dataset** selector (including descriptive `All`
  30,000-muon mixture) is preserved when the overlay checkbox is off.
  The individual raw-muon timelines/3D clouds still require separate
  original files and are not enabled by comparison mode.

No reprocessing or recombination of the physics caches is required.
The comparison helper can be sanity-checked without Garfield or Streamlit:

```bash
python3 test/mwpc-stage-c/pi_penning057_channel_compare_v1.py
```

All source metadata and SHA checksums remain in each dataset's
`provenance.json` and `temporal_correlations_v1.json`. Grouped statistics
are shown for all four datasets, regardless of which is selected.

*Do not run the original 1000-muon `prepare` command for this campaign.*


## 13 November 2026 — presentation-oriented Plot Studio update

The original pre-computed caches remain unchanged. This update only
changes which views are offered in the dashboard, and changes presentation
layout and read-only statistical plots.

- Retained the ion-weighted aggregate ion-birth azimuth view; hid its
  equal-avalanche-weighted duplicate from the menu (the data remain in NPZ).
- Hid the aggregate (u,v) ion-birth density, mean gain vs (u0,v0) and
  ion-birth delay vs azimuth. Their source data are preserved.
- Added **Ion-birth azimuth vs starting u0** alongside the v0 map.
  Both use 24 bins of their starting coordinate and 72 bins of
  5 degrees around the anode. Rows are normalized independently.
  128-point fully sampled avalanches are selected by default, and
  optional Heed-start coordinate slicing is provided.
- Both start-coordinate heatmaps offer **Compare with electron drift
  LUT500**. This uses the unchanged
  `test_runs/mwpc_stage_c/pi_pitch_fieldwire_500_p3_v1/pilot/p3p0_anode_only/garfield_electrons.csv`,
  with 900 spatial cells x 500 independent microscopic drift trials.
  Before use, the CSV SHA-256 and frozen geometry manifest are verified.
  The LUT histogram includes only collected electron endpoints with
  valid `phi_hit_deg`. For every same-coordinate start bin, it compares
  the electron *arrival* angle with a full-Penning Garfield ion-*birth*
  azimuth histogram. Each is normalized separately, and the displayed
  median conditional total-variation distance uses only bins containing
  at least 40 electron endpoints and 40 ion-producing avalanches.
  **These datasets are independent microscopic histories**: they are not
  paired trajectories, and they need not have the same population after
  attachment/loss selection. Arrival phi is defined around the actual
  collecting anode, while birth phi in this archived dataset is defined
  around the original nearest anode. Some avalanches involve a different
  wire. If birth histograms are ion-weighted, an additional gain-weighting
  difference exists; disable that weighting for the closest per-primary
  comparison. Do not report median TV as an isolated physical discrepancy
  in intrinsic avalanche angle without accounting for these effects.
- Added the u0 heatmap to the three-channel aggregate comparison.
  In channel overlay mode the drift comparison checkbox is hidden so
  the three independent-channel comparison remains unambiguous.
- Real-event pages now stack full-width charts vertically: true Geant4
  muon+clusters above its per-electron gain heatmap; one avalanche
  (u,v) sample above its angular histogram, then rotatable 3D (u,v,w).
  The local anode cross-section defaults to a tighter 0.1-mm zoom,
  equal physical axes use `constrain="domain"`, and birth-time colours
  may be relative to the earliest retained ion birth or the initiating
  Heed electron.
- All original full microscopic tracks, gain/ion counts and 128-position
  spatial reservoirs are untouched. Plots remain scientific visualizations
  of **ion creation**, not positive-ion drift or cathode induction.

### Update only the four modified source files in a worktree with other work

Avoid a broad `git pull` if there are local edits. Run the following
commands only if none of these four files has uncommitted local changes:

```bash
cd ~/na6p/src/NA6PRoot-mwpc
git fetch origin feature/mwpc-stage-c
for p in \
  test/mwpc-stage-c/pi_p3_noFW_penning057_interactive_dashboard_v1.py \
  test/mwpc-stage-c/pi_p3_noFW_penning057_anode_wrapping_v1.py \
  test/mwpc-stage-c/pi_penning057_individual_views_30k_v1.py \
  test/mwpc-stage-c/pi_penning057_channel_compare_v1.py
do
  test -z "$(git status --porcelain -- "$p")" || {
    echo "LOCAL EDIT: inspect before replacing $p"; exit 1;
  }
  git show "FETCH_HEAD:$p" > "$p.fetched" &&
    mv "$p.fetched" "$p"
done

PY=test_runs/mwpc_stage_c/.venv_pi_slides/bin/python
"$PY" -m py_compile \
  test/mwpc-stage-c/pi_p3_noFW_penning057_interactive_dashboard_v1.py \
  test/mwpc-stage-c/pi_p3_noFW_penning057_anode_wrapping_v1.py \
  test/mwpc-stage-c/pi_penning057_individual_views_30k_v1.py \
  test/mwpc-stage-c/pi_penning057_channel_compare_v1.py

"$PY" test/mwpc-stage-c/pi_p3_noFW_penning057_interactive_dashboard_v1.py selftest
"$PY" test/mwpc-stage-c/pi_p3_noFW_penning057_anode_wrapping_v1.py
"$PY" test/mwpc-stage-c/pi_penning057_individual_views_30k_v1.py selftest
"$PY" test/mwpc-stage-c/pi_penning057_channel_compare_v1.py
```

The local drift-LUT comparison needs the original full 450k-row
microscopic drift CSV and `frozen_config.json`. It cannot run with
only the compact 30k Plot Studio ZIP. To open both detailed-event
pages, the existing read-only SSHFS mount of the native Garfield
chunks is still required.
