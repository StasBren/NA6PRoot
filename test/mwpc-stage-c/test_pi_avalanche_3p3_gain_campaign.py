#!/usr/bin/env python3
"""Offline synthetic checks; zero Garfield physics and zero Google Drive writes."""
from __future__ import annotations
import csv
import json
import tempfile
from pathlib import Path
from unittest.mock import patch

import pi_avalanche_3p3_gain_campaign as av
import pi_publish_run_to_drive as pub

def run():
    assert [g[0] for g in av.CONFIGS] == [
        "p3_noFW","p3_FW75","p4_noFW","p4_FW75"]
    assert {(g[1],g[2]) for g in av.CONFIGS} == {
        (3.0,0.0),(3.0,75.0),(4.0,0.0),(4.0,75.0)}
    assert av.COMMON["gap_minus_mm"]==av.COMMON["gap_plus_mm"]==3.0
    with tempfile.TemporaryDirectory(prefix="pi_four_geometry_smoke_") as td:
        path=Path(td)
        csvfile=path/"gain.csv"
        with csvfile.open("w",newline="") as f:
            writer=csv.writer(f)
            writer.writerow(av.COLS)
            writer.writerow([0,10,9,7,1,2,9.0])
            writer.writerow([1,0,0,0,0,1,0])
            writer.writerow([2,20,19,13,1,4,10.0])
        a=av.analyze(csvfile,3)
        assert a["single_seed_events"]==3
        assert a["zero_collected_events"]==1
        assert a["effective_collected_gain_mean"]==20/3
        assert a["gain_p50"]==7 and a["gain_p90"]==13
        csvfile.write_text(csvfile.read_text().replace("2,20,19,13","2,20,19,-13"))
        try:
            av.analyze(csvfile,3)
        except ValueError:pass
        else:raise AssertionError("Negative avalanche counts were accepted")
        src=path/"physics"
        src.mkdir()
        (src/"summary.json").write_text('{"status":"complete"}')
        (src/"manifest.json").write_text('{"schema":"fixture"}')
        (src/"gain_results.csv").write_text("event,mean\n1,3\n")
        (src/"crossings.csv.gz").write_bytes(b"never-upload-heavy-raw")
        (src/"heavy.html").write_bytes(b"also-excluded")
        with patch.object(pub,"BUNDLE_ROOT",path/"stage"),\
             patch.object(pub,"REPO",av.REPO):
            d=pub.prepare(src,"fixture_test_3p3")
            assert (d/"summary.json").is_file()
            assert (d/"gain_results.csv").is_file()
            assert not (d/"crossings.csv.gz").exists()
            assert not (d/"heavy.html").exists()
            assert (d/"publication_manifest.json").is_file()
            assert pub.prepare(src,"fixture_test_3p3")==d
            (src/"summary.json").write_text('{"status":"altered"}')
            try:
                pub.prepare(src,"fixture_test_3p3")
            except ValueError:pass
            else:raise AssertionError("Mutable run was accepted after staging")
    print("PASS four p3/p4 ±3mm geometries, single-seed gain stats, "
          "negative-count rejection, Drive summary allowlist and immutable bundles.")

if __name__=="__main__":run()
