#!/usr/bin/env python3
"""Offline end-to-end fixture for the Stage-C station dashboard generator."""
from argparse import Namespace
import csv
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import numpy as np
from build_pi_station_dashboard import generate
from mwpc_checkerboard_projection import SOURCES, build_chambers
from pi_fixed_wire_occupancy import n_wires

def write_csv(path, fields, rows):
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

def main():
    geometry = build_chambers()
    with TemporaryDirectory(prefix="stage_c_dashboard_") as tmp:
        root = Path(tmp)
        crossings = root/"crossings"
        crossings.mkdir()
        modules, stations, electrons, geometric, arrays = [], [], [], [], {}
        for source in SOURCES:
            for st, group in enumerate(geometry):
                for c in group:
                    n = 2 if c.col == 0 else 0
                    modules.append(dict(source=source,scenario_ms2_z_cm=530,
                                        station=st,det_id=c.det_id,muon_crossings=n))
                    for j in range(n_wires(st)):
                        active = int(n == 2 and j == n_wires(st)//2)
                        electrons.append(dict(source=source,station=st,
                            det_id=c.det_id,wire_index=j,chamber_crossings=n,
                            n_ge1e=active,n_ge5e=0,n_le100ns=active))
                        geometric.append(dict(source=source,station=st,
                            det_id=c.det_id,wire_index=j,N_crossings=n,
                            n_geom=active))
                total = 2*sum(c.col == 0 for c in group)
                stations.append(dict(source=source,scenario_ms2_z_cm=530,
                    station=st,reference_states=total,
                    projected_muons_with_gas=total//2,
                    chamber_crossings_total=total))
                key = source+"_MS"+str(st)
                arrays[key] = np.array([[total,0],[0,0]])
                arrays[key+"_xedges"] = np.array([-5.,0.,5.])
                arrays[key+"_yedges"] = np.array([-5.,0.,5.])
        def write(path, rows):
            write_csv(path, tuple(rows[0]), rows)
        write(crossings/"modules.csv", modules)
        write(crossings/"stations.csv", stations)
        wf, gf, xyf = root/"wires.csv", root/"geometric.csv", root/"xy.npz"
        write(wf, electrons)
        write(gf, geometric)
        np.savez_compressed(xyf, **arrays)
        destination = root/"interactive.html"
        generate(Namespace(
            crossings=crossings, wire_counts=wf, geom_counts=gf, xy=xyf,
            template=Path(__file__).with_name("pi_station_dashboard_template.html"),
            out=destination))
        html = destination.read_text(encoding="utf-8")
        assert "__PI_STATION_DATA__" not in html
        marker = '<script id="stationData" type="application/json">'
        data = json.loads(html.split(marker,1)[1].split("</script>",1)[0])
        assert len(data["geometry"]) == 6
        assert len(data["chambers"][SOURCES[0]]) == 298
        assert len(data["chambers"][SOURCES[0]]["17"]["e1"]) == 77
        assert len(data["chambers"][SOURCES[0]]["52"]["e1"]) == 128
        print("PASS: 3 sources, 6 stations, 298 chambers and 109536 anodes")
        print("PASS: 18 XY histograms and standalone embedded HTML")

if __name__ == "__main__":
    main()
