# Heed-only 5 mm versus 6 mm gas-gap scan

This isolated study tests whether shortening the MWPC gas path changes primary
ionisation statistics and the probability of very low electron counts.
It deliberately does NOT infer spatial resolution or readout efficiency yet.

## Starting point

Branch: feature/mwpc-heed-gap-scan (based on feature/mwpc-avalanche-visualization).

Reuses the existing Garfield++ / Heed executable:
- test/mwpc-avalanche/phaseA_heed_statistics.cxx
- target: mwpc_phase_a_heed_stats

Reference configuration for BOTH runs:
- 10 GeV/c negative muon, perpendicular to the chamber
- Ar/CO2 70:30 at 293.15 K and 760 Torr
- anode diameter 30 um, pitch 4 mm, +1800 V
- 5 mm symmetric (2.5 + 2.5 mm) versus 6 mm symmetric (3 + 3 mm)
- fixed local u = 1 mm; magnetic field is irrelevant for this Heed-only statistic
- 20,000 tracks per configuration by default

**Important geometry detail:** Existing Heed code initializes the muon 0.1 mm
inside the +v cathode rather than on the exact plane. Thus the nominal
straight track lengths inside the gas are approximately 4.9 and 5.9 mm.
The analysis uses these lengths in its simple scaling benchmark.

The HV and wire parameters keep the analytic cell construction consistent.
No avalanche calculation is done in this test; the Heed primary-electron
statistics are not a scan of gas gain.

## Environment and build (WSL)

From your normal NA6PRoot source directory:

    source ~/na6p/env.sh
    source ~/na6p/install/garfieldpp/share/Garfield/setupGarfield.sh
    cd ~/na6p/src/NA6PRoot-mwpc
    git fetch origin
    git switch feature/mwpc-heed-gap-scan
    # If that branch has not been checked out locally:
    # git switch -c feature/mwpc-heed-gap-scan --track origin/feature/mwpc-heed-gap-scan

Build the existing Heed statistics target:

    cd test/mwpc-avalanche
    cmake -S . -B build -DGarfield_DIR="$HOME/na6p/install/garfieldpp/lib/cmake/Garfield"
    cmake --build build -j4 --target mwpc_phase_a_heed_stats

For the figures, install matplotlib in your Python environment if needed:

    python3 -m pip install matplotlib

## Run

Quick check (500 tracks per thickness):

    bash run_phaseA_heed_gap_scan.sh 500

More useful production statistics (20,000 tracks per thickness):

    bash run_phaseA_heed_gap_scan.sh 20000

A specific output directory (overwrites same-named CSV/log/figures there):

    bash run_phaseA_heed_gap_scan.sh 20000 "$HOME/na6p/heed_gap_results"

If no output directory is given, each execution creates a time-stamped
subdirectory under test/mwpc-avalanche/runs/.

## Files per run

- heed_5mm_tracks.csv, heed_6mm_tracks.csv: raw per-muon statistics;
  includes clusters, conduction-electron count, and total energy transfer.
- heed_5mm_clusters.csv, heed_6mm_clusters.csv: raw per-cluster data.
- heed_5mm.log, heed_6mm.log: Heed console output.
- summary.csv: count, mean, sample variance, SD, CV, median, upper
  quantiles and mean energy transfer for both thicknesses.
- low_tail.csv: empirical probabilities P(N_e < n) for selected n.
- comparison.txt: observed ratios and a simple independent-increment
  benchmark, noting the effective lengths 4.9 and 5.9 mm.
- electron_count_histogram.png: overlaid electron count distributions,
  with counts above the displayed range accumulated into the last bin.
- low_electron_tail.png: probability of fewer than n primary electrons.
- relative_fluctuations.png: observed relative standard deviation and
  expected 1/sqrt(effective length) comparison.

The raw CSVs are the actual simulated data; they allow further checks on
the outliers and on sample stability without rerunning Heed.

## Interpretation and limitations

1. A change in the lower tail of primary-electron count is only a **hint**
   toward a possible threshold-efficiency effect. Actual efficiency needs
   drift/attachment, gain, induced-strip response, shaping, noise and threshold.
2. The high-electron tail from energetic delta electrons makes the sample
   variance and CV noisy. Repeat with more tracks and independent runs
   before deciding whether differences from the 1/sqrt(L) benchmark matter.
3. The experiment varies the **ionisation path length**, holding the muon
   perpendicular. It is not yet a comparison of the complete MWPC geometry.
4. The 5 mm and 6 mm gaps are symmetric about the wire plane. The existing
   Prototype-3-like configuration is asymmetric (2 + 4 mm), a separate
   design choice that should be compared after this controlled first scan.
5. This is not connected to real Geant4 -> Heed hit handoff; the muons
   are generated directly in a standalone Heed/Garfield cell.
