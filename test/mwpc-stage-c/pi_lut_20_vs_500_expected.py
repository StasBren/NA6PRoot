#!/usr/bin/env python3
"""Fast J/psi MS0 20-vs-500 LUT convergence with exact mean Heed occupancy.

Geant4 chord -> analytic expected primary-electron occupancy for each
0.1 x 0.2 mm microscopic Garfield spatial bin. Heed model expectation:
mean_clusters_per_mm * mean_primaries_per_cluster * true gas chord length.
Geant4 positions and physical finite-anode alignment are retained.

THIS IS AN EXPECTATION OVER HEED RANDOM BIRTHS, not the original frozen
event-by-event RNG Heed realization. It removes Heed birth shot noise and
expensive per-muon RNG/multiplicity expansion; compare with the separate
original script if exact Heed realizations are needed. Same native Geant4
events and both independent Garfield 20/500 response maps.
"""
from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path

import numpy as np

import pi_lut_20_vs_500_convergence as old
from pi_2m_native_six_geometry_fold import crossing_provenance
from pi_2m_zero_path_compatible_resume import coincident_zero_chord
from pi_six_geometry_native_timing import first_indices, n_anodes

OUT=old.REPO/"test_runs/mwpc_stage_c/pi_pitch_fieldwire_20_vs_500_expected_v1"


def segment_fractions(u0, v0, u1, v1, nwire):
    """Partition a straight gas segment at every exact LUT grid boundary.

    Return spatial box, physical nearest-anode+1, and fractional 3D chord
    length. This is a deterministic geometric integration, NOT a claim
    of deterministic electron drift in the nonuniform MWPC E-field.
    """
    if not all(math.isfinite(x) for x in (u0,v0,u1,v1)):
        raise ValueError("Nonfinite native gas crossing coordinates")
    if max(abs(v0),abs(v1))>3.05:
        raise ValueError("Gas path endpoints outside v +/-3.05mm")
    cuts=[0.,1.]
    # This origin is the left Voronoi border of physical anode 0.
    # Periodic repetition extends into finite-detector edge margins.
    for a,b,origin,step in (
        (u0,u1,-nwire*old.PITCH/2.,0.1),
        (v0,v1,-3.,0.2),
    ):
        if a==b:
            continue
        low,high=min(a,b),max(a,b)
        first=math.floor((low-origin)/step)+1
        last=math.ceil((high-origin)/step)-1
        if first<=last:
            boundaries=origin+np.arange(first,last+1,dtype=float)*step
            t=(boundaries-a)/(b-a)
            cuts.extend(t[(t>0.)&(t<1.)].tolist())
    ts=np.unique(np.asarray(cuts,dtype=float))
    portions=np.diff(ts)
    tmid=(ts[:-1]+ts[1:])*0.5
    u=u0+(u1-u0)*tmid
    v=v0+(v1-v0)*tmid
    nearest,cell=first_indices(u,old.PITCH,nwire)
    if np.any((nearest < -1)|(nearest > nwire)):
        raise ValueError("Geant4 gas chord outside nearest-anode margin")
    iu=np.clip(np.floor((cell+old.PITCH/2.)/0.1).astype(int),0,29)
    iv=np.clip(np.floor((v+3.)/0.2).astype(int),0,29)
    if not np.isclose(portions.sum(),1.,atol=1.e-12,rtol=0.):
        raise AssertionError("Analytic Geant4 gas chord fractions must sum to 1")
    return iu*30+iv,nearest+1,portions


def expected_occupancies(station,chunks,max_hits,lam,marks):
    """True native gas tracks, integrated mean Poisson-Heed response."""
    nwire=n_anodes(station,old.PITCH)
    table=np.zeros((old.NCELL,nwire+2),dtype=np.float64)
    rate_per_mm=float(lam)*float(np.mean(marks))
    if not (rate_per_mm>0 and math.isfinite(rate_per_mm)):
        raise ValueError("Unphysical expected Heed electron density")
    seen=selected=zero_paths=0
    born=0.
    csv_hashes=[]
    for chunk in range(1,chunks+1):
        path,provenance=crossing_provenance("Jpsi",chunk)
        csv_hashes.append(provenance["csv_sha256"])
        with path.open(newline="",encoding="utf-8") as stream:
            for row in csv.DictReader(stream):
                if int(row["station"])!=station or abs(int(row["pdg"]))!=13:
                    continue
                seen+=1
                if max_hits and selected>=max_hits:
                    continue
                length=float(row["path_length_mm"])
                if length==0.:
                    if not coincident_zero_chord(row):
                        raise ValueError("Invalid exact-zero native chord")
                    zero_paths+=1
                    selected+=1
                    continue
                if not (math.isfinite(length) and 0.<length<1000.):
                    raise ValueError(f"Invalid native gas path {length}")
                bin_idx,wire_idx,dt=segment_fractions(
                    float(row["u_in_mm"]),float(row["v_in_mm"]),
                    float(row["u_out_mm"]),float(row["v_out_mm"]),nwire)
                expectation=rate_per_mm*length
                np.add.at(table,(bin_idx,wire_idx),expectation*dt)
                born+=expectation
                selected+=1
                if selected%10000==0:\n                    print(f"  progress: {selected:,} MS{station} crossings, "\n                          f"{born:,.0f} expected Heed electrons",flush=True)\n        print(f"J/psi chunk {chunk}/{chunks}: MS{station} crossings={selected:,}, "
              f"analytic expected primaries={born:,.2f}",flush=True)
    if born<=0 or not np.isclose(table.sum(),born,rtol=1.e-10,atol=1.e-6):
        raise AssertionError("Mean-Heed spatial occupancy conservation failed")
    return table,{
        "station":station,"chunks":chunks,"source":"jpsi40-bore",
        "selected_muon_gas_crossings":selected,
        "all_muon_gas_crossings_seen":seen,
        "zero_length_crossings":zero_paths,
        "heed_primaries":born,"heed_occupancy_mode":"analytic_expectation",
        "heed_density_per_mm":rate_per_mm,
        "heed_occupancy_definition":"lambda_Heed * mean(multiplicity) * "
            "fractional length of each native straight Geant4 gas segment "
            "inside each 0.1x0.2mm periodic Garfield LUT box",
        "native_csv_sha256":csv_hashes,
        "physical_anodes":nwire,"geometrical_lut_cells":old.NCELL,
        "full_source_exposure":chunks==40 and not max_hits,
    }


def main():
    old.occupancies=expected_occupancies
    old.OUT=OUT
    old.main()
    # The shared renderer predates the analytic occupancy mode.
    # Correct its model label explicitly, without touching frozen code.
    out=OUT
    if "--out" in sys.argv:
        out=Path(sys.argv[sys.argv.index("--out")+1])
    for flag in sys.argv:
        if flag.startswith("--out="):
            out=Path(flag.split("=",1)[1])
    path=out/"compare.html"
    doc=path.read_text(encoding="utf-8")
    before="shared deterministic Heed births."
    if before not in doc:
        raise AssertionError("Shared HTML renderer label changed unexpectedly")
    doc=doc.replace(before,"analytically expected Heed electrons "
        "(integrated over native Geant4 gas chords; not seeded Heed draws).")
    doc=doc.replace(
        "<p>Orange = 20 microscopic trials/cell; blue = 500.",
        "<p>Heed model: exact conditional mean spatial occupancy; no "
        "per-event Heed shot noise. Both LUTs share the same Geant4 "
        "track-derived occupancy.</p><p>Orange = 20 microscopic "
        "trials/cell; blue = 500.")
    path.write_text(doc,encoding="utf-8")
    print("MODE: analytic expected Heed cell occupancy, "
          "not the original per-event Heed RNG realization",flush=True)


if __name__=="__main__":
    main()
