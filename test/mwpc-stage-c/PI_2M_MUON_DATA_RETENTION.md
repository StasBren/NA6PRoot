# NA60+ 6M-parent production: muon truth retention contract

Status: Geant4 campaign is FROZEN with 100k parents per source successfully
tested. New work must NOT change native Geant4 generator, original
simulator binary, installed field map, detector layout or production
validator; doing so would break comparability and frozen hashes.

## Mandatory source-of-truth for every 50k-parent source/chunk

Keep ALL these originals indefinitely until copied to verified backup:
- MCKine.root: retained-event TParticle truth, momentum, production
  vertex, particle PDG and retained mother links. Ancestry might be
  pruned and remapped; never assume a stored index is an original ID.
- HitsMuonSpec.root: event-indexed MuonSpec vectors; EVERY gas volume
  crossing, in/out position and in/out momentum, detector ID, remapped
  track ID, Edep, crossing time and status. The gasStepsMuonSpec tree
  holds finer steps with positive Edep only; it omits zero-Edep steps.
- MuonReferencePlanes.root: virtual/interpolated muon trajectory at many
  Z planes, upstream and downstream. NOT real sensitive-volume hits.
- geometry.root: actual material and chamber geometry; na6pLayout.ini.
- raw_workers: essential original worker-local Geant4 hit IDs, gas
  step IDs and trajectory IDs; these cannot necessarily be joined to
  merged MCKine IDs without explicit validated remapping.
- completed.json, command.json, simulation.log, validation.log,
  campaign.json: event counts, seeds, source/worker and chunk mapping,
  SHA hashes and actual time costs.

Never overwrite the first 2 completed chunks per source or any resulting
100k dashboard data. All 40 chunks/source belong to SAME frozen campaign.

## Per-muon analysis fields that can be recovered from retained originals

Global identity: source, chunk (1..40), event (0..49999 local), track ID,
gas hit index, detector ID and station. This is a composite key;
bare event numbers and bare track IDs are not globally unique.

Parent event: original mother resonance, production vertex (x,y,z,t,
where t is actually provided), PDG, four-momentum, event/generator
weight if supplied by the source, original mu+mu- pair correlation,
pair mass, pT, rapidity and opening angle. Retain ZERO-hit events
for later detection-efficiency denominators. MCKine may contain
pruned ancestry: only use actual available original truth.

Each chamber crossing: global x,y,z for gas entry and exit,
local Garfield u,v,w in mm, lab and local in/out momentum, muon signed
angles, Geant4 deposited energy and path length, TOF and boundary
statuses. Preserve ALL hits of the same muon and other species,
not just earliest hit or nearest wire. The existing
geant_muon_crossings.csv exports most of these, but ROOT remains
authoritative for full event truth and steps.

Tracking to changed readouts: keep all real reference-plane states
with original worker IDs, detector transforms, absorber history and
momentum before/after spectrometer. Use for alignment, magnetic bending,
material budget, scattering, reconstruction residuals, rates and
incidence angles when physical transport geometry is unchanged.

Environmental provenance: archive ACTUAL bytes of field maps,
transport cuts, original gas and materials, exact ROOT/Geant4/NA6PRoot
versions, executable/library SHA hashes, generator source and seeds.
A field-map checksum alone is not useful if the underlying file is
later changed or deleted. Preserve exact six Garfield LUT contents
and Heed track/cluster calibration files separately.

## What must be regenerated later, but does NOT require a new Geant4

Original gas crossing truth plus deterministic shared Poisson/Heed seed
can regenerate per-cluster positions and multiplicities, per-primary
electron birth sites and six-LUT collection/migration/arrival times.
This requires preserving exact RNG algorithm, seed, marks calibration,
Garfield LUT data and unchanged cross-source/chunk/event/hit keys.
The current merged 20-position-bin JSON is only an aggregate, NOT a
per-muon or per-electron data bank.

Future analogue prototypes additionally need Garfield gain/avalanche
statistics, the weighting field, Ramo signal, U/V strip/charge sharing,
pulse shaping, bandwidth, noise, thresholds, dead time and reconstruction.
They are NOT stored in current 6 LUT timing histograms. Real B-field
gas drift requires different Garfield physics (existing LUT B=0).

A sampled per-electron truth bank could be made later from existing
gas crossings without re-running Geant4, but saving individual
electron paths for the entire run is potentially very large and should
be optional. For parent-event bootstrap, retain the tuple
(source,chunk,event) and never count electrons as independent parents.

## When the expensive Geant4 run truly DOES need repeating

Changes to absorber, chamber gas volume, thickness, material or
physical placement change Geant4 transport or Edep and require
new physical simulation if precision matters. A different magnetic
field also changes charged-particle trajectories. Realistic
minimum-bias heavy-ion occupancy, centrality, accidental overlays,
and non-signal backgrounds are NOT contained in the forced Jpsi/Phi/
Omega signal samples; separate simulated inputs are needed.

Changes to pitch/wires/strip readout/electronics within a physically
equivalent material layout generally permit reusing Geant4 tracks
while rerunning only Heed/Garfield/response and reconstruction.

## Mandatory gate

Execute the checked script before launching the all-weekend run:
    cd ~/na6p/src/NA6PRoot-mwpc
    git pull --ff-only
    bash test/mwpc-stage-c/run_pi_weekend_2m.sh check

This now calls the independent ROOT muon truth audit for the first
six verified 50k-parent chunks: source+chunk identity, MCKine track
truth joins, actual MuonSpec hit objects, event counts, original worker
files, gas-step and reference-plane trees, and native CSV input SHA.
The ROOT audit samples 200 distributed events PER chunk; it does not
prove every future event. Every subsequent chunk still passes its
original full ROOT structural validator and SHA checkpoint.

The current weekend producer will save 120 Geant4 chunks, 120 CSV
exports and 720 Heed+six-LUT station shards if adequate storage and
awake runtime permit. Final merge of all 40 chunks/source and
parent-event statistical errors require a separate validated step
and do not require repeating the original Geant4.
