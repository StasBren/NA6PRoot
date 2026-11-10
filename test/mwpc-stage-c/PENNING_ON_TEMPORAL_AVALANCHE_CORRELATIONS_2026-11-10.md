# Penning-ON: avalanche timing, gain, asymmetry and cluster chronology

Uses existing 1000 TRUE, SHA-checked J/psi muons at p3_noFW.
Current 30,000-muon cluster campaign and all original files are untouched.

## Physics caveat

The saved microscopic code starts the Garfield avalanche at each original
TrackHeed electron starting (u,v,w,t), not at an independently detected
arrival at the anode. Arrival time and arrival angle were NOT measured.
An early-time variable t10 is therefore defined as the 10th percentile
of up to 128 uniformly sampled ion BIRTH times per avalanche relative
to the muon entry clock. This is NOT t_drift or anode arrival angle.
Default 128-position-only selection reduces quantile sample-size bias,
but early birth timing can still correlate intrinsically with growth
and gain, not only the preceding electron drift.

## New six editable views in MWPC Plot Studio

1. Mean collected gain vs early ion-birth time: exact gain, sampled time.
2. Mean exact |A_ion,hemi| vs early ion-birth time.
3. Two ion-weighted birth azimuth histograms for early/late avalanches.
   Select earliest/latest quartiles or custom 25 ns / 100 ns time windows.
4. Early ion-birth time vs original starting distance |v0|.
5. Mean gain vs both early birth time and starting |v0|.
6. One real muon timeline: exact Heed starting times by cluster and
   approximate early birth times with marker size = gain, colour = exact A.

Every view provides dynamic plot-specific physical interpretation.
These observations are ion-BIRTH source correlations, NOT actual
induced cathode signals, not ion drift or weighting-field calculation.

## Update two scripts WITHOUT losing local plot edits

    cd ~/na6p/src/NA6PRoot-mwpc
    git fetch origin feature/mwpc-stage-c
    PY=test_runs/mwpc_stage_c/.venv_pi_slides/bin/python
    DASH=test/mwpc-stage-c/pi_p3_noFW_penning057_interactive_dashboard_v1.py
    TIME=test/mwpc-stage-c/pi_p3_noFW_penning057_timing_correlations_v1.py
    mkdir -p ~/na6p/script_backups
    cp -p "$DASH" ~/na6p/script_backups/dashboard_before_temporal_update.py
    test ! -e "$TIME" || cp -p "$TIME" ~/na6p/script_backups/temporal_before_update.py
    git show "FETCH_HEAD:$DASH" > "${DASH}.new" && mv "${DASH}.new" "$DASH"
    git show "FETCH_HEAD:$TIME" > "${TIME}.new" && mv "${TIME}.new" "$TIME"
    "$PY" -m py_compile "$DASH" "$TIME"
    "$PY" "$TIME" selftest
    "$PY" "$DASH" selftest

## Build the new, separate compact timing cache ONCE (no Garfield)

    set -o pipefail
    "$PY" -u "$TIME" prepare | tee ~/na6p/penning1000_temporal_prepare.log

Original interactive_dashboard_v1/dashboard_data.npz remains unchanged.
New files temporal_correlations_v1.npz and temporal_correlations_v1.json
are created in the same folder, refusing overwrite and revalidating SHA.

## Open the same light-theme dashboard

    "$PY" -m streamlit run "$DASH" --server.address localhost --server.port 8501 --theme.base light

Open http://localhost:8501 in the Windows web browser.

The first tests use a full-size sample of all 1,000 muons, not a handpicked
single example. The one-muon timeline then illustrates the behaviour.
Exact electron arrival timing and velocity-angle would require a new
isolated, small Garfield diagnostic rather than falsifying these fields.

## Plot Studio update: custom time ranges, working legend stats, relevant options

Angular comparison (`t10_angles`) now defaults to **Custom time intervals (ns)**.
Earlier and later curves each have independent FROM and TO fields.
The selection is `From <= t10 < To`, and t10 is the 10th percentile
of sampled ion-birth creation times **relative to original muon entry**.
Windows can have different widths. Empty or underpopulated (<30)
groups display a warning rather than inventing a signal.
Quantile groups remain available for equal-population comparisons.
Each curve has its own editable text; optional legend statistics
refer to that group's sample `t10` **mean and population RMS** in ns.

Mean/RMS legend text is now actually constructed for relevant views:
- collected gain (inclusive and positive only): exact electron gains
- avalanche / muon asymmetry: exact signed A values
- physical ion-creation delay: weighted binned approximation, marked ≈
- per-avalanche sampled time span: sampled time spans in ns
- mean |A| vs gain: underlying |A| population
- time-correlation mean gain / mean |A|: eligible original observables
- early/late angular comparison: separately selected group t10 in ns.

No meaningless linear angular mean is shown for 0/±180-degree periodic
azimuths; the mean/RMS checkbox is omitted from those non-time angular views.

Sidebar: gain bin width appears only for gain histograms; angular
bin width only for angular histograms; asymmetry bin count only for
asymmetry histograms; time-bin grouping only for the ion-delay histogram;
heatmap colour-scale controls only for heatmaps; 128-point reliability
filter omitted from the one-muon timeline. 2D temporal heatmaps now
respect the logarithmic colour scale control too.

No new data preparation needed. Refresh BOTH updated Python scripts
using the file-specific git-show commands documented above, run both
selftests and restart Streamlit. Cache NPZ files stay unchanged.
