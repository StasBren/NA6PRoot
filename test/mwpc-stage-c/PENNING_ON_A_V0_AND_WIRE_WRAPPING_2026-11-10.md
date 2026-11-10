# Penning-ON exact A(v0) and physical avalanche wrapping views

Adds 5 scientific views to the existing 1000-real-muon Penning ON
Streamlit/Plotly MWPC Plot Studio. Reuses its already-prepared temporal
correlations cache, SHA-verified 1000 Geant4 J/psi muons, Garfield
electron collected count, exact ion-sector callbacks and at most 128
original sampled birth positions per avalanche.

Geometry is original p3_noFW: 3mm anode pitch, 30um anode diameter,
planes at v=+/-3mm, 1800V Ar/CO2 70:30 with Penning transfer r=0.57.
No source, existing analysis cache, ion csv, or running 30k campaign is
modified; no new all-1000 prepare is necessary.

## 1. Signed mean A vs original Heed starting v0

The true starting electron's v0 is relative to the anode-wire plane,
with v=0 at the wire centres. Exact per-avalanche A =
(Nion_up - Nion_down)/(Nion_up + Nion_down). The graph plots
the arithmetic equal-avalanche mean signed A, separately in 0.25mm
v0 bins. All avalanches with a defined A are included, even if
they collected zero electrons at anode or produced <128 ions.

The uncertainty band resamples 1000 complete independent real
Geant4 parent muons (350 bootstrap replicates); it does not
misrepresent microscopic electrons as independent muon parents.
Minimum bin N=40; insufficiently populated bins are blank.

## 2. Exact opposite-side ion fraction vs |v0|

fOpp=(1-sign(v0)*A)/2, relative to the original electron's
START hemisphere, not its unknown measured arrival side.
This is the exact opposite share among the counted upper/lower
ion births (boundary counts excluded). The plot excludes
|v0|<0.025mm because the initial side is ambiguous there.
Bootstrap and 0.25mm bins follow the preceding view.

## 3. Ion-birth azimuth conditioned on starting v0

Rows condition on original v0; columns show the 72 saved five-degree
phi bins relative to nearest anode wire. Each row is normalised to
unit mass, so shapes can be compared at different initial positions.
Choose ion-birth weighting (individual saved ion = true ions /
number sampled) or equal weight per avalanche.
Optionally restrict to 128-sample avalanches for uniform precision.
The azimuth is where ions were CREATED; it is NOT incoming electron
angle, ion drift path or an induced cathode signal.

## 4. Circular angular dispersion vs collected gain

R = absolute value of the average of exp(i*phi), measured over
sampled ion birth azimuths. Plot mean (1-R) vs collected gain,
with 128-sample avalanches and complete muon-level bootstrap.
R close to one means narrowly concentrated birth angles;
R close to zero means circular spread OR antipodal cancellation.
Therefore this statistic does not demonstrate continuous wrapping.

## 5. Single REAL microscopic avalanche around a wire

Select any real Geant4 parent among the 1000 and any original
TrackHeed electron seed. The default seed produces many ions
so its saved 128-point birth cloud is informative. The viewer
reads ONLY that muon's original sealed 20-parent chunk and checks
SHA of raw seed, Heed primary and ion-birth coordinate CSV files.
The actual collected count and ion count are cross-checked against
the pre-existing temporal correlations NPZ.

2D cross-section: u relative to that seed's original local nearest
3mm anode centre, v relative to the anode plane, with equal physical
axis scaling. The real anode is a 15um RADIUS grey circle. The
initial Heed electron start is a blue star (may be outside zoom).
Sampled ion-birth points are coloured by physical birth delay
relative to the initiating Heed electron. Choose the cross-section
view half-width from 0.05 to 3mm; only within-range points are
visible and there is no silent rescaling.

3D rotatable view: true sampled (delta_u,v,w) ion-birth points
and a real-scale 30um-diameter wire cylinder along w, with each
point coloured by its ion-birth time. No simulated later ion
drift, electron arrival surface, induced cathode current or
full unsampled microscopic avalanche is implied.

## Fetch the two changed/new scripts without Git merge conflicts

    set -e
    cd ~/na6p/src/NA6PRoot-mwpc
    PY=test_runs/mwpc_stage_c/.venv_pi_slides/bin/python
    DASH=test/mwpc-stage-c/pi_p3_noFW_penning057_interactive_dashboard_v1.py
    WRAP=test/mwpc-stage-c/pi_p3_noFW_penning057_anode_wrapping_v1.py
    git fetch origin feature/mwpc-stage-c
    mkdir -p ~/na6p/script_backups
    for f in "$DASH" "$WRAP"; do
      if [ -f "$f" ]; then
        cp -p "$f" ~/na6p/script_backups/"$(basename "$f").before_wrapping"
      fi
      git show "FETCH_HEAD:$f" > "${f}.new" && mv "${f}.new" "$f"
    done
    "$PY" -m py_compile "$DASH" "$WRAP"
    "$PY" "$WRAP"
    "$PY" "$DASH" selftest

Existing temporal script and both NPZ caches remain unchanged.

## Start the same dashboard

    "$PY" -m streamlit run "$DASH" --server.address localhost --server.port 8501 --theme.base light

Open http://localhost:8501 in Windows and select one of the 5 new
plot types listed under 'Edit the plot'.

The source was reviewed and the new Python module includes
self-tests without Garfield; full UI results must be verified
against local Penning-ON data in WSL before interpreting numbers.
