#!/usr/bin/env python3
"""Command-line front end of the cathode response engine (one persistent process per batch).

  engine_cli.py info  KERNEL_DIR                         kernel identity (sha256) + transport / weighting ids
  engine_cli.py run   --kernel KERNEL_DIR --kernel-id SHA --moments-dir DIR --out OUT [--model K2A] [--impl vector]
                      [--format npz|replay] [--repeat 1]
        Evaluates every DIR/*.moments.csv (written by `p3_birth_moments --batch`) with one kernel opened once.
        Writes OUT/<event>.npz (ion, electron, labels, t_center_ns, meta) or replay-format directories, and
        OUT/engine_timing.json: process start-up, kernel verify/load, per-event moments-read and evaluate times
        (cold = first event, warm = the rest), peak RSS.
  engine_cli.py bench --kernel KERNEL_DIR --moments FILE [--scales 1,2,4,8,16] [--repeat 10]
        Latency and memory scaling with the number of clouds (clouds replicated in time: same physics per cloud),
        vector vs loop implementation.
Normal operation needs no development-muon index, run directory or S0 file.
"""
import argparse
import json
import resource
import sys
import time
import tracemalloc
from pathlib import Path

T_START = time.perf_counter()
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np  # noqa: E402
import cathode_engine as ce  # noqa: E402

T_IMPORT = time.perf_counter() - T_START


def rss_mib():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024


def cmd_info(a):
    k = ce.KernelTable(a.kernel, verify=True, preload=False)
    print(json.dumps(k.describe(), indent=1))


def cmd_run(a):
    k = ce.KernelTable(a.kernel, expected_id=a.kernel_id, verify=True, preload=True)
    eng = ce.CathodeResponseEngine(k, model=a.model, half_strips=a.half_strips, window_ns=a.window_ns)
    a.out.mkdir(parents=True, exist_ok=False)
    files = sorted(Path(a.moments_dir).glob("*.moments.csv"))
    if not files:
        raise SystemExit("no *.moments.csv in " + str(a.moments_dir))
    ev = []
    for i, f in enumerate(files):
        name = f.name[:-len(".moments.csv")]
        t0 = time.perf_counter(); m = ce.read_moments_csv(f); t1 = time.perf_counter()
        best = None
        for _ in range(a.repeat):
            s = time.perf_counter(); r = eng.evaluate(m, a.impl); e = time.perf_counter() - s
            best = e if best is None else min(best, e)
        t2 = time.perf_counter()
        if a.format == "npz":
            np.savez(a.out / f"{name}.npz", ion=r.ion, electron=r.electron, labels=np.array(r.labels),
                     t_center_ns=r.t_center_ns, meta=json.dumps(r.meta))
        else:
            r.write_replay_format(a.out, name)
        ev.append({"event": name, "clouds": m.n, "read_moments_s": t1 - t0, "evaluate_s": best,
                   "write_s": time.perf_counter() - t2, "warm": i > 0, "peak_rss_mib_so_far": rss_mib()})
        print(f"{name}: {m.n} clouds, evaluate {best * 1e3:.2f} ms{' (cold: first event)' if i == 0 else ''}")
    warm = [e["evaluate_s"] for e in ev if e["warm"]]
    rep = {"tool": "engine_cli run", "model": a.model, "impl": a.impl, "kernel": k.describe(), "import_s": T_IMPORT,
           "events": ev, "n_events": len(ev), "cold_first_event_evaluate_s": ev[0]["evaluate_s"],
           "warm_evaluate_s_median": float(np.median(warm)) if warm else None,
           "process_wall_s": time.perf_counter() - T_START, "peak_rss_mib": rss_mib(),
           "throughput_events_per_s_warm": (len(warm) / sum(warm)) if warm else None}
    (a.out / "engine_timing.json").write_text(json.dumps(rep, indent=1) + "\n")
    wm = rep["warm_evaluate_s_median"]
    print(f"PASS engine run: {len(ev)} events; kernel verify {k.verify_s:.2f} s, load {k.load_s:.2f} s; "
          f"warm evaluate median {'n/a' if wm is None else f'{wm * 1e3:.3f} ms'}")


def cmd_bench(a):
    k = ce.KernelTable(a.kernel, verify=False, preload=True)
    eng = ce.CathodeResponseEngine(k, model=a.model)
    m0 = ce.read_moments_csv(a.moments)
    rows = []
    for sc in [int(x) for x in a.scales.split(",")]:
        cols = {c: np.tile(v, sc) for c, v in m0.cols.items()}
        cols["t_ns"] = cols["t_ns"] + np.repeat(np.arange(sc), m0.n) * 1e-3   # distinct clouds, same physics
        m = ce.Moments(cols, m0.base_U, m0.base_V)
        res = {}
        for impl in ("vector", "compiled", "loop"):
            if impl == "compiled":
                try:
                    ce.CathodeResponseEngine.compiled_library()
                except OSError:
                    continue
            if impl == "loop" and sc > a.loop_max_scale:
                continue
            eng.evaluate(m, impl)
            tracemalloc.start(); eng.evaluate(m, impl); peak = tracemalloc.get_traced_memory()[1]; tracemalloc.stop()
            ts = []
            for _ in range(a.repeat if impl == "vector" else max(1, a.repeat // 5)):
                s = time.perf_counter(); eng.evaluate(m, impl); ts.append(time.perf_counter() - s)
            res[impl] = {"median_ms": 1e3 * float(np.median(ts)), "min_ms": 1e3 * float(np.min(ts)), "temp_peak_mib": peak / 2**20}
        rows.append({"clouds": m.n, **res})
        print(f"clouds {m.n:6d}: " + "; ".join(f"{i} median {v['median_ms']:.2f} ms, temp {v['temp_peak_mib']:.1f} MiB" for i, v in res.items()))
    out = {"kernel": k.describe(), "model": a.model, "rows": rows, "process_peak_rss_mib": rss_mib()}
    if a.out:
        Path(a.out).write_text(json.dumps(out, indent=1) + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    i = sub.add_parser("info"); i.add_argument("kernel", type=Path)
    r = sub.add_parser("run")
    r.add_argument("--kernel", type=Path, required=True); r.add_argument("--kernel-id", required=True)
    r.add_argument("--moments-dir", type=Path, required=True); r.add_argument("--out", type=Path, required=True)
    r.add_argument("--model", default="K2A", choices=list(ce.MODELS)); r.add_argument("--impl", default="vector", choices=["vector", "compiled", "loop"])
    r.add_argument("--format", default="npz", choices=["npz", "replay"]); r.add_argument("--repeat", type=int, default=1)
    r.add_argument("--half-strips", type=int, default=5); r.add_argument("--window-ns", type=float, default=7000.0)
    b = sub.add_parser("bench")
    b.add_argument("--kernel", type=Path, required=True); b.add_argument("--moments", type=Path, required=True)
    b.add_argument("--model", default="K2A", choices=list(ce.MODELS)); b.add_argument("--scales", default="1,2,4,8,16")
    b.add_argument("--repeat", type=int, default=10); b.add_argument("--loop-max-scale", type=int, default=4)
    b.add_argument("--out", type=Path)
    a = ap.parse_args()
    {"info": cmd_info, "run": cmd_run, "bench": cmd_bench}[a.cmd](a)


if __name__ == "__main__":
    main()
