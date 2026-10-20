#!/usr/bin/env python3
"""Synthetic test: LUT per-box mean/RMS/error and lookup throughput."""
import numpy as np
from pi_centered_lut_quality_benchmark import assess, microbenchmark


def main():
    pool = np.empty((600, 8, 3), dtype=float)
    pool[:, :, 0] = 1
    pool[:, :, 1] = 0
    pool[:, :, 2] = np.array([10., 20., 30., 40.,
                               10., 20., 30., 40.])
    pool[0, 7, 0] = 0
    pool[0, 7, 2] = np.nan
    pool[1, 7, 2] = np.nan  # collected, unknown time
    rows = assess(pool, bootstrap=20)
    assert len(rows) == 600
    assert rows[0]["n_collected"] == 7
    assert rows[0]["n_collected_with_known_time"] == 7
    assert rows[1]["n_collected"] == 8
    assert rows[1]["n_collected_unknown_time"] == 1
    assert rows[2]["n_collected"] == 8
    assert rows[2]["mean_drift_ns"] == 25
    assert rows[2]["rms_drift_ns"] == np.sqrt(125)
    assert rows[2]["mean_stat_error_ns"] > 0
    assert rows[2]["rms_stat_error_ns"] > 0
    perf = microbenchmark(pool, n_electrons=3000, batch=1000)
    assert perf["n_synthetic_positions"] == 3000
    assert perf["batches"] == 3
    assert perf["lut_nanoseconds_per_electron"]["n"] == 3
    print("PASS: per-box collected/unknown counters, mean/RMS errors "
          "and vectorized lookup timing separation")


if __name__ == "__main__":
    main()
