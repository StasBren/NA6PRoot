#!/usr/bin/env python3
"""Simulation-FREE regression of 100k benchmark stage in frozen 2M suite.

Creates three synthetic PRE-COMPLETED checkpoints with dummy files.
Checks that the real runner skips completed chunks, does not launch
na6psim, and NEVER publishes a false full-campaign-completed marker.
"""
from __future__ import annotations
import json
import tempfile
from pathlib import Path
from types import SimpleNamespace

from run_pi_geant4_campaign import (
    BASE_SEED, CHANNELS, ESSENTIAL, simulate,
)


def main():
    with tempfile.TemporaryDirectory(prefix="pi_2m_campaign_test_") as t:
        suite=Path(t)
        for ch in CHANNELS:
            chunk=suite/ch.lower()/"chunk0001"
            chunk.mkdir(parents=True)
            (chunk/"completed.json").write_text(json.dumps({
                "channel":ch,"events":50_000,"seed":BASE_SEED[ch],
                "chunk":1}),encoding="utf-8")
            for name in ESSENTIAL:
                (chunk/name).write_bytes(b"synthetic presence only")
        # Full physical campaign target remains frozen at 2M per channel.
        args=SimpleNamespace(events_per_source=2_000_000,
                             chunk_size=50_000,
                             max_chunks_per_source=1,
                             hours=0.)
        simulate(args,suite,{
            "parallel_runner":"/NONEXISTENT/na6psim_parallel",
            "root":"/NONEXISTENT/root",
        })
        assert not (suite/"campaign_complete.json").exists()
        assert not any(suite.rglob("chunk0002"))
        for ch in CHANNELS:
            assert (suite/ch.lower()/"chunk0001"/"completed.json").exists()
    print("PASS: frozen 2M/channel benchmark cap does NOT run Geant4, "
          "preserves existing chunks and cannot mark full campaign complete")


if __name__=="__main__":
    main()
