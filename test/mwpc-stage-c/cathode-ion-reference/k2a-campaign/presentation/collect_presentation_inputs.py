#!/usr/bin/env python3
"""Collect the small, read-only inputs for the K2A signal-response presentation into ONE new upload folder.

  collect_presentation_inputs.py NEW_UPLOAD_DIR --analysis ANALYSIS_DIR [--stress DIR] [--queue QUEUE.tsv]
        [--pf PREFLIGHT_DIR] [--val-root DIR] [--kernel KERNEL_DIR] [--green GREEN_DIR]
        [--examples auto|SLOT,SLOT,...] [--extra-examples SLOT,...] [--max-file-mb 25] [--max-example-mb 60]

Reads only; writes only inside NEW_UPLOAD_DIR (which must not exist). Never copies ion paths (vertices.f64), numpy
tables (*.npy), ROOT files, group exports or anything of a sealed slot.

  report     ANALYSIS_DIR/{validation_report.md,.json, per_muon.csv, combined_from.txt, provenance.txt, slots.txt}
  per slot   for every slot of the report (its run from combined_from.txt, else ANALYSIS_DIR): outcome.json, steps.tsv,
             timing.tsv, cmp_*.json/.txt, shift_assignment.json, label_union.json, shaped_*/per_muon_Tp*.jsonl and
             *.log, prod/compare_moments.md/.log, prod/births.log, prod/birth_moments.csv.json, K2A_prod/result.json,
             K2A_imposed/result.json, far_dir.txt; from the S0 far dir: groups/groups.json and
             S0groups/g_*/timing_attempts.json. Plus each run's provenance.txt and slots.txt.
  stress     --stress DIR: its report files and the same per-slot files (kept in stress/, separate provenance)
  queue      --queue FILE (the overnight queue), copied as is
  preflight  --pf: inventory/slots.csv WITHOUT sealed rows, domain/per_event.csv, and every non-sealed moments file
             moments_nonsealed/**/muon_*.moments.csv(+ .event.json, .frame.json): population features for the join
  probes     --val-root (default ~/na6p/k2a_validation): wire_probe, kernel_probe, closure_sladder, closure_0637,
             domain_scan (md/json/csv/txt only, no waveforms)
  kernel     --kernel / --green: kernel.json, axes.json, and a small slice of K.npy (K0, Kc1, KA2 at the central D node
             of every wire and radius, all tau; npz, < 1 MB) for the kernel figure
  examples   waveforms of a few slots: shaped_layout/muon_M/* (window-aligned S0_T1 / S0_T0 / K2A / K2Aimp, csv) and
             the full-length S0_T1 replay + prompt (to 100 us). auto = the accuracy-set muon with the median and the one
             with the largest model U shaped L2 at Tp 100, plus 0637 if present; --extra-examples adds slots.
Writes MANIFEST.tsv (source, destination, bytes, sha256), SKIPPED.tsv (with reason) and SUMMARY.txt.
"""
import argparse
import csv
import hashlib
import json
import os
import shutil
import statistics
import sys
from pathlib import Path

NEVER = (".npy", ".root", ".f64", ".f32", ".bin")
PER_SLOT = ["outcome.json", "steps.tsv", "timing.tsv", "shift_assignment.json", "far_dir.txt",
            "prod/compare_moments.md", "prod/compare_moments.log", "prod/births.log", "prod/birth_moments.csv.json",
            "K2A_prod/result.json", "K2A_imposed/result.json", "label_union.log", "time_grid.json"]
PER_SLOT_GLOB = ["cmp_*.json", "cmp_*.txt", "shaped_*/per_muon_Tp*.jsonl", "shaped_*/*.json", "shaped_*.log",
                 "shaped_layout/muon_*/label_union.json", "*.json"]
REPORT = ["validation_report.md", "validation_report.json", "per_muon.csv", "combined_from.txt", "provenance.txt",
          "slots.txt", "selection.csv", "criteria_freeze_check.txt"]


class Collector:
    def __init__(self, out, max_file):
        self.out, self.max_file, self.man, self.skip = out, max_file, [], []

    def copy(self, src, dst_rel, cap=None):
        src = Path(src)
        if not src.is_file():
            return False
        if src.suffix in NEVER or "vertices" in src.name:
            self.skip.append((str(src), "never copied (ion paths / binary tables)")); return False
        n = src.stat().st_size
        if n > (cap or self.max_file):
            self.skip.append((str(src), f"{n / 1e6:.1f} MB above the size cap")); return False
        dst = self.out / dst_rel
        if dst.exists():
            return True
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)
        self.man.append((str(src), str(dst_rel), n, hashlib.sha256(dst.read_bytes()).hexdigest()))
        return True

    def tree(self, src, dst_rel, pats, cap=None):
        src = Path(src)
        for p in pats:
            for f in sorted(src.glob(p)):
                if f.is_file():
                    self.copy(f, Path(dst_rel) / f.relative_to(src), cap)


def runs_of(analysis):
    cf = analysis / "combined_from.txt"
    return [Path(l.strip()) for l in cf.read_text().splitlines() if l.strip()] if cf.is_file() else [analysis]


def slot_home(runs):
    home = {}
    for r in runs:
        st = r / "slots.txt"
        for s in (st.read_text().split() if st.is_file() else []):
            home.setdefault(f"{int(s):04d}", r)
    return home


def sealed_slots(pf):
    p = pf / "inventory" / "slots.csv" if pf else None
    if not p or not p.is_file():
        return set()
    return {r["slot"] for r in csv.DictReader(open(p, newline="")) if "sealed" in r.get("role", "")}


def collect_run(C, analysis, tag, sealed):
    for f in REPORT:
        C.copy(analysis / f, Path(tag) / "report" / f)
    runs = runs_of(analysis)
    home = slot_home(runs)
    bad = sorted(set(home) & sealed)
    if bad:
        raise SystemExit(f"refusing: sealed slots {bad} appear in {analysis}")
    for i, r in enumerate(runs):
        for f in ("provenance.txt", "slots.txt", "selection.csv", "criteria_freeze_check.txt"):
            C.copy(r / f, Path(tag) / "runs" / f"{i:03d}_{r.name}" / f)
    for s, r in sorted(home.items()):
        V = r / "muons" / f"muon_{s}"
        d = Path(tag) / "muons" / f"muon_{s}"
        for f in PER_SLOT:
            C.copy(V / f, d / f)
        C.tree(V, d, PER_SLOT_GLOB)
        fd = V / "far_dir.txt"
        if fd.is_file():
            E = Path(fd.read_text().strip())
            C.copy(E / "groups" / "groups.json", d / "far" / "groups.json")
            C.tree(E / "S0groups", d / "far" / "S0groups", ["g_*/timing_attempts.json"])
    return home


def pick_examples(analysis, home):
    pm = analysis / "per_muon.csv"
    if not pm.is_file():
        return []
    rows = [r for r in csv.DictReader(open(pm, newline="")) if r.get("accuracy_v1_primary") in ("PASS", "FAIL")]
    vals = []
    for r in rows:
        try:
            vals.append((float(r["model_U_shaped_L2_Tp100"]), r["slot"]))
        except (KeyError, ValueError):
            pass
    vals.sort()
    pick = []
    if vals:
        med = statistics.median(v for v, _ in vals)
        pick += [min(vals, key=lambda x: abs(x[0] - med))[1], vals[-1][1]]
    if "0637" in home:
        pick.append("0637")
    return list(dict.fromkeys(pick))


def collect_examples(C, home, slots, cap):
    done = []
    for s in slots:
        r = home.get(s)
        if r is None:
            C.skip.append((s, "example slot not in the report")); continue
        V = r / "muons" / f"muon_{s}"
        lay = V / "shaped_layout" / f"muon_{s}"
        files = [f for f in sorted(lay.rglob("*.csv"))] if lay.is_dir() else []
        fd = V / "far_dir.txt"
        if fd.is_file():
            E = Path(fd.read_text().strip())
            files += [E / "S0_local" / "replay" / "waveforms.csv", E / "S0_local" / "prompt" / "prompt_waveforms.csv"]
        files += [V / "K2A_prod" / "replay" / "waveforms.csv", V / "K2A_prod" / "prompt" / "prompt_waveforms.csv"]
        tot = sum(f.stat().st_size for f in files if f.is_file())
        if tot > cap:
            C.skip.append((s, f"example waveforms {tot / 1e6:.1f} MB above --max-example-mb")); continue
        for f in files:
            if f.is_file():
                rel = f.relative_to(lay) if lay in f.parents else (
                    Path("S0T1_full") / f.relative_to(f.parents[1]) if "S0_local" in f.parts else
                    Path("K2A_full") / f.relative_to(f.parents[1]))
                C.copy(f, Path("examples") / f"muon_{s}" / rel, cap)
        done.append(s)
    return done


def kernel_slice(C, kd, gd):
    for f in ("kernel.json", "kernel_meta.json"):
        C.copy(kd / f, Path("kernel") / f)
    if gd:
        C.copy(gd / "axes.json", Path("kernel") / "green_axes.json")
    try:
        import numpy as np
        K = np.load(kd / "K.npy", mmap_mode="r")            # (basis, side, wire, r, D, tau)
        meta = json.loads((kd / "kernel.json").read_text())
        names = meta.get("basis") or meta.get("names") or []
        idx = [i for i, n in enumerate(names) if n in ("K0", "Kc1", "KA2")] or list(range(min(3, K.shape[0])))
        jd = K.shape[4] // 2
        sl = np.asarray(K[idx][:, :, :, :, jd, :], dtype=np.float32)
        dst = C.out / "kernel" / "K_slice_centralD.npz"
        np.savez_compressed(dst, K=sl, basis=np.array([names[i] if names else str(i) for i in idx]), D_index=jd,
                            shape_full=np.array(K.shape))
        C.man.append((str(kd / "K.npy") + f" [slice basis {idx}, D index {jd}]", "kernel/K_slice_centralD.npz",
                      dst.stat().st_size, hashlib.sha256(dst.read_bytes()).hexdigest()))
    except Exception as e:  # noqa: BLE001 - optional item, reported
        C.skip.append((str(kd / "K.npy"), f"kernel slice not written: {e}"))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out", type=Path); ap.add_argument("--analysis", type=Path, required=True)
    ap.add_argument("--stress", type=Path); ap.add_argument("--queue", type=Path); ap.add_argument("--pf", type=Path)
    ap.add_argument("--val-root", type=Path, default=Path.home() / "na6p" / "k2a_validation")
    ap.add_argument("--kernel", type=Path); ap.add_argument("--green", type=Path)
    ap.add_argument("--examples", default="auto"); ap.add_argument("--extra-examples", default="")
    ap.add_argument("--max-file-mb", type=float, default=25); ap.add_argument("--max-example-mb", type=float, default=60)
    a = ap.parse_args()
    if a.out.exists():
        raise SystemExit(f"refusing: {a.out} exists")
    if not (a.analysis / "per_muon.csv").is_file():
        raise SystemExit(f"no per_muon.csv in {a.analysis}")
    a.out.mkdir(parents=True)
    C = Collector(a.out, a.max_file_mb * 1e6)
    sealed = sealed_slots(a.pf)
    home = collect_run(C, a.analysis, "main", sealed)
    shome = collect_run(C, a.stress, "stress", sealed) if a.stress and (a.stress / "per_muon.csv").is_file() else {}
    if a.stress and not shome:
        C.skip.append((str(a.stress), "no per_muon.csv: stress run not collected"))
    if a.queue:
        C.copy(a.queue, Path("queue") / a.queue.name) or C.skip.append((str(a.queue), "queue file missing"))
    if a.pf:
        inv = a.pf / "inventory" / "slots.csv"
        if inv.is_file():
            rows = [r for r in csv.DictReader(open(inv, newline="")) if "sealed" not in r.get("role", "")]
            dst = a.out / "preflight" / "inventory_slots_nonsealed.csv"; dst.parent.mkdir(parents=True, exist_ok=True)
            with open(dst, "w", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=list(rows[0]), lineterminator="\n"); w.writeheader(); w.writerows(rows)
            C.man.append((str(inv) + " [sealed rows removed]", "preflight/inventory_slots_nonsealed.csv",
                          dst.stat().st_size, hashlib.sha256(dst.read_bytes()).hexdigest()))
        C.copy(a.pf / "domain" / "per_event.csv", Path("preflight") / "domain_per_event.csv")
        mn = a.pf / "moments_nonsealed"
        for f in sorted(mn.rglob("muon_*.moments.csv*")) if mn.is_dir() else []:
            s = f.name.split("muon_")[1][:4]
            if s in sealed:
                C.skip.append((str(f), "sealed slot")); continue
            if f.name.endswith((".moments.csv", ".event.json", ".frame.json")):
                C.copy(f, Path("preflight") / "moments" / f.name)
    for d in ("wire_probe", "kernel_probe", "closure_sladder", "closure_0637", "domain_scan"):
        C.tree(a.val_root / d, Path("probes") / d, ["*.md", "*.json", "*.csv", "*.txt", "*/*.md", "*/*.json"])
    if a.kernel:
        kernel_slice(C, a.kernel, a.green)
    ex = pick_examples(a.analysis, home) if a.examples == "auto" else [f"{int(x):04d}" for x in a.examples.split(",") if x.strip()]
    ex += [f"{int(x):04d}" for x in a.extra_examples.split(",") if x.strip() and f"{int(x):04d}" not in ex]
    allhome = {**shome, **home}
    if set(ex) & sealed:
        raise SystemExit("refusing: an example slot is sealed")
    got = collect_examples(C, allhome, ex, a.max_example_mb * 1e6)
    with open(a.out / "MANIFEST.tsv", "w") as fh:
        fh.write("source\tdestination\tbytes\tsha256\n")
        fh.writelines(f"{s}\t{d}\t{n}\t{h}\n" for s, d, n, h in C.man)
    with open(a.out / "SKIPPED.tsv", "w") as fh:
        fh.write("item\treason\n"); fh.writelines(f"{s}\t{r}\n" for s, r in C.skip)
    tot = sum(m[2] for m in C.man)
    summ = [f"upload folder {a.out}", f"files {len(C.man)}, {tot / 1e6:.1f} MB; skipped {len(C.skip)} (SKIPPED.tsv)",
            f"report slots {len(home)} (main), {len(shome)} (stress)", f"examples {got}",
            f"sealed slots excluded by the preflight inventory: {len(sealed)} (none copied)"]
    (a.out / "SUMMARY.txt").write_text("\n".join(summ) + "\n")
    print("\n".join(summ))


if __name__ == "__main__":
    sys.exit(main())
