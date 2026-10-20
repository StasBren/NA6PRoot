#!/usr/bin/env python3
"""Synthetic check: measured chunk timings only; mean/RMS and errors."""
import json
from pathlib import Path
from tempfile import TemporaryDirectory

from pi_stagec_performance_report import analyse, stats


def main():
    zero = stats([])
    assert zero["n"] == 0 and zero["mean"] is None
    one = stats([2])
    assert one["mean"] == 2 and one["rms"] == 0
    assert one["mean_error"] is None
    two = stats([1, 3])
    assert two["n"] == 2 and two["mean"] == 2 and two["rms"] == 1
    assert two["mean_error"] is not None and two["mean_error"] > 0
    assert two["rms_error"] is not None and two["rms_error"] > 0
    assert stats([1, 3]) == two  # seeded bootstrap reproducibility

    with TemporaryDirectory(prefix="pi_wall_runtime_") as tmp:
        suite = Path(tmp)
        (suite/"campaign.json").write_text(
            json.dumps({"frozen_inputs": {
                "channels": ["Jpsi"], "events_per_source": 6,
                "chunk_size": 2}}))
        channel = suite/"jpsi"
        for i in range(1, 4):
            p = channel/f"chunk{i:04d}"
            p.mkdir(parents=True)
            row = {"channel": "Jpsi", "chunk": i, "events": 2}
            if i != 3:  # legacy completed checkpoint without a clock
                t = float(i)
                row["wall_seconds"] = {"geant4": t,
                                       "root_validation": 0.5*t,
                                       "chunk_total": 2*t}
                row["timing_clock"] = \
                    "time.perf_counter, elapsed wall seconds"
            (p/"completed.json").write_text(json.dumps(row))
        report = analyse(suite)
        group = report["channels"]["Jpsi"]
        assert group["timed_chunks"] == 2
        assert group["legacy_chunks_without_timings"] == 1
        assert group["completed_events"] == 6
        assert group["stages"]["geant4"]["mean"] == 0.75
        assert group["stages"]["geant4"]["rms"] == 0.25
        (channel/"chunk0002"/"completed.json").unlink()
        group = analyse(suite)["channels"]["Jpsi"]
        assert group["completed_events"] == 4
        assert group["stages"]["geant4"]["mean_error"] is None
    print("PASS: true-only timings, historical missing markers, "
          "mean/RMS and seeded chunk-bootstrap errors")


if __name__ == "__main__":
    main()
