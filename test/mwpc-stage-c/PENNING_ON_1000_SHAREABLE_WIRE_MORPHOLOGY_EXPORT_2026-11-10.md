# Penning-ON 1000-muon: compact upload package for avalanche-response LUT design

This is a read-only extractor, not a physics re-run and not the ongoing
30,000-muon campaign. It verifies the SHA of original C++ Garfield outputs
for 50 completed chunks and exports the already-saved up-to-128 ion-birth
coordinates per avalanche into a compact analysis ZIP.

## Why this package

Do not build independent gain, phi, and longitudinal w-smearing LUTs before
checking their correlations. The export records (a) exact per-seed collected
electron gain and ion-birth count/upper-lower asymmetry, (b) sample-derived
per-avalanche circular angle and w shape, and (c) joint ion birth
probabilities P(phi,delta_w | starting v0, collected gain).

Longitudinal delta_w is birth_w - original Heed start_w. Consequently
it INCLUDES drift diffusion from the initial electron position and
any intrinsic avalanche longitudinal spread. An independently observed
electron arrival w/t/phi is NOT currently saved and would need a
separate limited calibration with the drift/multiplication stages recorded.

## Fetch ONLY the new exporter, preserving local changes

```bash
set -e
cd ~/na6p/src/NA6PRoot-mwpc
git fetch origin feature/mwpc-stage-c
EX=test/mwpc-stage-c/pi_p3_noFW_penning057_export_wire_model_v1.py
PY=test_runs/mwpc_stage_c/.venv_pi_slides/bin/python
mkdir -p ~/na6p/script_backups
if [ -f "$EX" ]; then
  cp -p "$EX" ~/na6p/script_backups/wire_export_before_update.py
fi
git show "FETCH_HEAD:$EX" > "${EX}.new"
mv "${EX}.new" "$EX"
"$PY" -m py_compile "$EX"
"$PY" "$EX" selftest
```

## Produce ZIP

```bash
set -o pipefail
"$PY" -u "$EX" export | tee ~/na6p/penning_on_wire_export.log
"$PY" "$EX" inspect
ZIP=test_runs/mwpc_stage_c/pi_p3_noFW_jpsi1000_penning057_v1/wire_model_export_v1/NA60_PENNING_ON_1000_WIRE_MODEL_REVIEW_V1.zip
ls -lh "$ZIP"
```

Copy into the Windows user's Downloads folder to upload here:

```bash
cp "$ZIP" /mnt/c/Users/stanm/Downloads/
```

Attach the ZIP file in the next chat turn. Its filename is
NA60_PENNING_ON_1000_WIRE_MODEL_REVIEW_V1.zip.

## ZIP contents

- README.txt: physical definitions, caveats, units and file layout.
- summary.json: exact all-1000 census, summary quantiles and 50 native
  chunk SHA references.
- per_avalanche_morphology.csv: every original starting electron, including
  zero-ion avalanches. Exact G, ion up/down counts, A, initial coordinates,
  cluster membership, physical timing. Sampled longitudinal w spread
  (RMS and 5-95% width), angular mean/resultant/90% shortest arc/largest
  empty arc/number of occupied 30-degree sectors, radial spread and times.
- joint_ion_source_histograms.npz: ion-weighted AND equal-avalanche
  4D source hist [v0, collected gain, birth_phi, birth_delta_w] plus
  [v0, birth_phi, ion-creation delay] and [v0, birth_phi, radial distance].
  The axis edges and units are embedded as arrays.
- representative_original_birth_clouds.csv: up to 8 original, uniformly
  sampled avalanche clouds from each 20-muon chunk, chosen across
  signed asymmetry, initial v0, gain, along-wire and azimuthal spread.
  Each cloud preserves its original 3D+time points.
- 01..06 PNGs: white-background first-look figures A(v0), angular
  distribution vs v0, w spread vs v0, phi 90%-arc vs collected gain
  and phi-by-longitudinal-displacement for both start hemispheres.
- file_sha256.json: integrity hashes of all packaged outputs.

## Physics limits

- C++ p3_noFW finite 3-mm anode patch; 30-um diameter, +/-3mm cathodes.
- Gas Ar/CO2 70:30, Penning transfer r=0.57, 1800V.
- Two differently weighted histograms answer DIFFERENT questions:
  ion-weighted represents where all CREATED ions are born; equal
  avalanche represents the shape of a typical ion-producing avalanche.
- All spatial/timing morphology values come from <=128 positions per
  avalanche, not a complete birth cloud. Require full-128 samples when
  comparing per-avalanche azimuth and w widths vs gain.
- A exact up/down is ion CREATION hemisphere, not induced cathode charge.
- The angle phi is relative to the nearest wire where EACH ion is born;
  rare cases with multiple ion-producing wires are counted/flagged.
- A circular resultant near zero may mean broad spread OR two opposite
  narrow lobes, so do not infer a continuous sheath from that alone.
- The raw original 1000 muons have one SHA-selected crossing per J/psi
  parent; they are not a flux-unbiased station-crossing distribution.

## Proposed next model

Compare P(phi,delta_w | v0,G) against factorisation
P(phi | v0,G)*P(delta_w | v0,G). If correlations remain, preserve
them in a joint template or mixture model. Additional small Garfield
calibration should record true arrival position (phi,w), arrival time
and early multiplication, to disentangle diffusion in the drift
from the intrinsic avalanche spreading around / along the wire.
