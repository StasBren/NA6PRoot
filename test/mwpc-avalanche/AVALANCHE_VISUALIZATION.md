# Microscopic avalanche GIF

This branch adds a presentation-oriented visualisation of one real
`Garfield::AvalancheMicroscopic` avalanche.

## Physics chain shown

```text
one low-energy seed electron
  -> Garfield / Magboltz microscopic gas collisions
  -> ionising collisions in the high-field region near the anode
  -> secondary electrons
  -> stochastic avalanche multiplication
  -> electron collection / attachment
```

The trajectory points and their times come from Garfield's microscopic
drift-line storage. The potential/electric-field background is sampled directly
from the same `ComponentAnalyticField` used for the avalanche.

This is separate from the Heed muon GIF: Heed is not needed here because this
animation begins with one already-created conduction electron.

## Build

```bash
source ~/na6p/env.sh
source ~/na6p/install/garfieldpp/share/Garfield/setupGarfield.sh

cd ~/na6p/src/NA6PRoot-mwpc/test/mwpc-avalanche

cmake -S . -B build \
  -DGarfield_DIR="$HOME/na6p/install/garfieldpp/lib/cmake/Garfield"

cmake --build build -j4
```

## Generate the two default GIFs

```bash
bash make_avalanche_demo_gifs.sh
```

This produces:

```text
avalanche_centered_contours.gif
avalanche_offaxis_contours.gif
```

The first starts one electron directly above the central wire. The second
starts it 0.8 mm off axis. Both use actual Garfield equipotential contours.

## Optional field-direction overlay

The C++ event does not need to be rerun. Re-render the already-generated
centered event:

```bash
python3 render_avalanche_gif.py \
  --prefix avalanche_centered \
  --output avalanche_centered_field.gif \
  --background both \
  --half-width-mm 1.25
```

Supported background modes are:

```text
contours
field
both
none
```

The arrows show the local electric-field direction; the contours show the
actual analytic potential returned by Garfield.

## Main files

```text
phaseA_avalanche_visualization.cxx
render_avalanche_gif.py
make_avalanche_demo_gifs.sh
```

The C++ program writes:

```text
<prefix>_meta.txt
<prefix>_field.csv
<prefix>_endpoints.csv
<prefix>_paths.csv
```

The GIF timing uses the microscopic Garfield time stored along the electron
drift lines.

## Important boundary

The absolute avalanche gain is still uncalibrated. Penning transfer has not
been tuned to prototype data. The GIF is therefore a microscopic visualisation
of the present Phase-A reference model, not a calibrated measurement of final
chamber gain.
