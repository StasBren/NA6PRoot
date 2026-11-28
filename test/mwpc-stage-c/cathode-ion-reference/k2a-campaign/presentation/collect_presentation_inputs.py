#!/usr/bin/env python3
"""Collect the small inputs for the K2A signal-response presentation into ONE new upload folder, read-only.

  collect_presentation_inputs.py NEW_UPLOAD_DIR --analysis REPORT_DIR --pf PREFLIGHT_DIR
        [--also-report OTHER_REPORT_DIR ...] [--stress DIR] [--queue QUEUE.tsv] [--val-root DIR]
        [--runs-glob GLOB,...] [--kernel KERNEL_DIR] [--green GREEN_DIR]
        [--examples auto|none|SLOT,...] [--extra-examples SLOT,...] [--max-file-mb 25] [--max-example-mb 60]

READ-ONLY. Every input is only opened for reading; the only directory created is NEW_UPLOAD_DIR, which must not exist
and must not lie inside any input directory (so a running campaign is never written to). Never copied: ion paths
(vertices.f64), numpy tables (*.npy, except a separately written small kernel slice), ROOT files, group exports, and
anything of a sealed slot. All checks run before NEW_UPLOAD_DIR is created where possible.

--analysis   the aggregate report the presentation uses (k2a_validation_report.py output: per_muon.csv,
             validation_report.md/.json, combined_from.txt, provenance.txt, slots.txt). Any report works, e.g. a NEW
             aggregate of the running queue; the report directory itself is never modified. --also-report copies the
             report files (only) of further reports, e.g. the frozen 153-event snapshot, under reports/<name>/.
per slot     for every slot of --analysis (run directory from combined_from.txt / slots.txt, else --runs-glob under
             --val-root, default overnight_800/muon_*): outcome.json, steps.tsv, timing.tsv, cmp_*.json/.txt,
             shift_assignment.json, shaped_*/per_muon_Tp*.jsonl, prod/compare_moments.md/.log, prod/births.log,
             prod/birth_moments.csv.json, K2A_prod/result.json, label_union.json and other top-level json; from the
             S0 far dir groups/groups.json and S0groups/g_*/timing_attempts.json; each run's provenance.txt, slots.txt
--pf         REQUIRED. inventory/slots.csv (must be readable, with slot and role columns and at least one sealed row,
             so that sealed exclusion is verifiable; copied WITHOUT sealed rows), domain/per_event.csv, and every
             non-sealed p3_birth_moments product under moments_nonsealed/ with its relative path preserved:
             <id>.moments.csv, <id>.moments.csv.event.json, <id>.moments.csv.frame.json, <id>.moments.csv.json,
             <id>.error.txt, batch_timing.json; also <id>.event.json / <id>.frame.json if a converter wrote those names.
             <id> = muon_NNNN or NNNN.
--stress     another report (e.g. the high-multiplicity muons), collected the same way under stress/ (own provenance)
--kernel     kernel.json, and a small slice of K.npy written as kernel/K_slice_centralD.npz: bases K0, Kc1, KA2 at the
             D node nearest 0, every side, wire, radius and tau, with the axes. Generated ONLY if kernel.json lists those
             basis names and its green_axes (wires, r_um, D_mm, tau_ns) match K.npy's shape; otherwise reported, not
             substituted. --green adds the Green table's axes.json.
examples     window-aligned waveforms (shaped_layout/muon_M/*) plus full-length S0_T1 and K2A replay/prompt for a few
             slots: auto = accuracy-set muons with the median and the largest model U shaped L2 (Tp 100) and 0637.
Outputs: MANIFEST.tsv (source, destination, bytes, sha256), SKIPPED.tsv, CHECKS.tsv (check, status, detail),
MISSING_ESSENTIAL.tsv, SUMMARY.txt. Exit 0 = COMPLETE, 2 = INCOMPLETE (essential data missing: the folder is written
and every gap is listed), 1 = refused (nothing written).
"""
import argparse
import csv
import hashlib
import json
import os
import re
import shutil
import statistics
import sys
from collections import Counter, defaultdict
from pathlib import Path

NEVER = (".npy", ".root", ".f64", ".f32", ".bin")
PER_SLOT = ["outcome.json", "steps.tsv", "timing.tsv", "shift_assignment.json", "far_dir.txt",
            "prod/compare_moments.md", "prod/compare_moments.log", "prod/births.log", "prod/birth_moments.csv.json",
            "K2A_prod/result.json", "K2A_imposed/result.json", "label_union.log"]
PER_SLOT_GLOB = ["cmp_*.json", "cmp_*.txt", "shaped_*/per_muon_Tp*.jsonl", "shaped_*/*.json", "shaped_*.log",
                 "shaped_layout/muon_*/label_union.json", "*.json"]
REPORT = ["validation_report.md", "validation_report.json", "per_muon.csv", "combined_from.txt", "provenance.txt",
          "slots.txt", "selection.csv", "criteria_freeze_check.txt"]
COMPLETED = {"VALID", "VALID_T1_ONLY", "OUT_OF_DOMAIN", "S0_FAILED", "FAIL"}
PF_NAME = re.compile(r"^(?P<id>.+?)(?P<kind>\.moments\.csv\.event\.json|\.moments\.csv\.frame\.json|\.moments\.csv\.json|"
                     r"\.moments\.csv|\.event\.json|\.frame\.json|\.error\.txt)$")
KIND = {".moments.csv": "moments", ".moments.csv.event.json": "event", ".event.json": "event",
        ".moments.csv.frame.json": "frame", ".frame.json": "frame", ".moments.csv.json": "timing", ".error.txt": "error"}
SLICE_BASES = ("K0", "Kc1", "KA2")


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def slot_of(eid):
    m = re.fullmatch(r"(?:muon_)?(\d{4})", eid)
    return m.group(1) if m else None


class Collector:
    def __init__(self, out, max_file):
        self.out, self.max_file = out, max_file
        self.man, self.skip, self.checks, self.missing = [], [], [], []

    def check(self, name, ok, detail="", level="FAIL"):
        self.checks.append((name, "PASS" if ok else level, detail))
        return ok

    def need(self, slot, item, why):
        self.missing.append((slot, item, why))

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
        shutil.copyfile(src, dst)                       # content only: no metadata written back anywhere
        self.man.append((str(src), str(dst_rel), n, sha(dst)))
        return True

    def tree(self, src, dst_rel, pats, cap=None):
        src = Path(src)
        for p in pats:
            for f in sorted(src.glob(p)):
                if f.is_file():
                    self.copy(f, Path(dst_rel) / f.relative_to(src), cap)

    def wrote(self, path, label):
        self.man.append((label, str(path.relative_to(self.out)), path.stat().st_size, sha(path)))


def read_inventory(pf):
    inv = pf / "inventory" / "slots.csv"
    if not inv.is_file():
        raise SystemExit(f"refusing: preflight inventory {inv} not found")
    rows = list(csv.DictReader(open(inv, newline="")))
    if not rows or not {"slot", "role"} <= set(rows[0]):
        raise SystemExit(f"refusing: {inv} is empty or lacks slot/role columns")
    sealed = {r["slot"] for r in rows if "sealed" in r["role"]}
    if not sealed:
        raise SystemExit(f"refusing: {inv} has no sealed rows, so sealed exclusion cannot be verified")
    return inv, rows, sealed


def report_rows(d):
    p = d / "per_muon.csv"
    if not p.is_file():
        raise SystemExit(f"refusing: no per_muon.csv in {d}")
    rows = list(csv.DictReader(open(p, newline="")))
    for r in rows:
        r["slot"] = f"{int(r['slot']):04d}"
    return rows


def run_dirs(report, rows, val_root, globs, C, tag):
    cf = report / "combined_from.txt"
    runs = [Path(l.strip()) for l in cf.read_text().splitlines() if l.strip()] if cf.is_file() else [report]
    home = {}
    for r in runs:
        st = r / "slots.txt"
        for s in (st.read_text().split() if st.is_file() else []):
            home.setdefault(f"{int(s):04d}", r)
    for row in rows:
        s = row["slot"]
        if s in home and (home[s] / "muons" / f"muon_{s}").is_dir():
            continue
        hit = [d.parents[1] for g in globs for d in sorted((val_root / g).parent.glob(Path(g).name + f"/muons/muon_{s}"))]
        if hit:
            home[s] = hit[0]
            if len(hit) > 1:
                C.check(f"{tag}: run dir of {s}", False, f"{len(hit)} candidates, using {hit[0]}", "WARN")
    return runs, home


def collect_report(C, report, tag, sealed, inv_slots, val_root, globs, pf_have):
    rows = report_rows(report)
    bad = sorted({r["slot"] for r in rows} & sealed)
    if bad:
        raise SystemExit(f"refusing: sealed slots {bad} in {report}")
    for f in REPORT:
        C.copy(report / f, Path(tag) / "report" / f)
    runs, home = run_dirs(report, rows, val_root, globs, C, tag)
    for i, r in enumerate(sorted({*runs, *home.values()}, key=str)):
        for f in ("provenance.txt", "slots.txt", "selection.csv", "criteria_freeze_check.txt"):
            C.copy(r / f, Path(tag) / "runs" / f"{i:03d}_{r.name}" / f)
    prov = defaultdict(set)
    for r in {*home.values()}:
        p = r / "provenance.txt"
        kv = dict(l.split("=", 1) for l in p.read_text().splitlines() if "=" in l) if p.is_file() else {}
        for k in ("git_commit", "kernel", "strip_eval_sha256", "criteria_v1_sha256", "births_bin_sha256"):
            prov[k].add(kv.get(k, "MISSING"))
    for k, v in prov.items():
        C.check(f"{tag}: one {k} across runs", len(v) == 1 and "MISSING" not in v, "; ".join(sorted(v))[:300],
                "FAIL" if k == "kernel" else "WARN")         # one kernel is essential; harness differences are flagged
    cnt = Counter(r.get("execution", "") for r in rows)
    acc = Counter(r.get("accuracy_v1_primary", "") for r in rows)
    for r in rows:
        s, ex = r["slot"], r.get("execution", "")
        d = Path(tag) / "muons" / f"muon_{s}"
        if s not in inv_slots:
            C.need(s, "preflight inventory row", "slot not in the non-sealed inventory")
        if ex not in COMPLETED:
            C.check(f"{tag}: {s} completed", False, f"execution '{ex}' (in progress or inconclusive): not required", "WARN")
            continue
        for kind in ("moments", "event", "frame"):
            if kind not in pf_have.get(s, set()):
                C.need(s, f"preflight {kind} file", "not found under moments_nonsealed")
        V = home.get(s, Path("/nonexistent")) / "muons" / f"muon_{s}"
        if not V.is_dir():
            C.need(s, "per-slot run directory", "not found via combined_from.txt, slots.txt or --runs-glob"); continue
        if not (home[s] / "provenance.txt").is_file():
            C.need(s, "run provenance.txt", str(home[s]))
        for f in PER_SLOT:
            C.copy(V / f, d / f)
        C.tree(V, d, PER_SLOT_GLOB)
        for f in ("outcome.json", "steps.tsv", "timing.tsv"):
            if not (V / f).is_file():
                C.need(s, f, str(V))
        if r.get("accuracy_v1_primary") in ("PASS", "FAIL"):
            for f in ("cmp_S0T1_vs_K2A.json", "shaped_refS0T1_ionelec/per_muon_Tp100.jsonl"):
                if not (V / f).is_file():
                    C.need(s, f, "accuracy-set slot")
        if ex == "FAIL" and (V / "steps.tsv").is_file():
            fails = [l.split("\t")[0] for l in (V / "steps.tsv").read_text().splitlines() if "\tFAIL\t" in l]
            C.check(f"{tag}: {s} failing steps", False, ", ".join(fails) or "none listed", "INFO")
        fd = V / "far_dir.txt"
        if fd.is_file():
            E = Path(fd.read_text().strip())
            C.copy(E / "groups" / "groups.json", d / "far" / "groups.json")
            C.tree(E / "S0groups", d / "far" / "S0groups", ["g_*/timing_attempts.json"])
    return rows, home, cnt, acc


def collect_preflight(C, pf, inv, inv_rows, sealed):
    keep = [r for r in inv_rows if r["slot"] not in sealed]
    dst = C.out / "preflight" / "inventory_slots_nonsealed.csv"; dst.parent.mkdir(parents=True, exist_ok=True)
    with open(dst, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(inv_rows[0]), lineterminator="\n"); w.writeheader(); w.writerows(keep)
    C.wrote(dst, f"{inv} [{len(inv_rows) - len(keep)} sealed rows removed]")
    C.check("preflight inventory read", True, f"{len(inv_rows)} rows, {len(sealed)} sealed (excluded), {len(keep)} kept")
    C.copy(pf / "domain" / "per_event.csv", Path("preflight") / "domain_per_event.csv")
    have, unparsed, nfile = defaultdict(set), [], 0
    mn = pf / "moments_nonsealed"
    if not C.check("preflight moments_nonsealed present", mn.is_dir(), str(mn)):
        return have
    for root, _, files in os.walk(mn):
        for fn in sorted(files):
            src = Path(root) / fn
            rel = Path("preflight") / "moments_nonsealed" / src.relative_to(mn)
            if fn == "batch_timing.json":
                C.copy(src, rel); continue
            m = PF_NAME.match(fn)
            if not m:
                continue
            s = slot_of(m.group("id"))
            if s is None:
                unparsed.append(fn); continue
            if s in sealed:
                C.skip.append((str(src), "sealed slot")); continue
            if C.copy(src, rel):
                have[s].add(KIND[m.group("kind")]); nfile += 1
    n_ok = sum(1 for s in have if {"moments", "event", "frame"} <= have[s])
    C.check("preflight moments collected", n_ok > 0, f"{nfile} files; {n_ok} slots with moments+event+frame; "
            f"{sum(1 for s in have if 'error' in have[s])} with error.txt; unparsed names {len(unparsed)}")
    C.check("no sealed slot among collected preflight files", not (set(have) & sealed), "")
    return have


def kernel_slice(C, kd, gd):
    C.copy(kd / "kernel.json", Path("kernel") / "kernel.json")
    if gd:
        C.copy(gd / "axes.json", Path("kernel") / "green_axes.json")
    why = None
    try:
        import numpy as np
        meta = json.loads((kd / "kernel.json").read_text())
        names, ax = meta.get("basis"), meta.get("green_axes") or {}
        K = np.load(kd / "K.npy", mmap_mode="r")
        if not isinstance(names, list):
            why = "kernel.json has no 'basis' list"
        elif any(b not in names for b in SLICE_BASES):
            why = f"basis {names} lacks one of {SLICE_BASES}"
        elif any(k not in ax for k in ("wires", "r_um", "D_mm", "tau_ns")):
            why = "kernel.json green_axes lacks wires / r_um / D_mm / tau_ns"
        else:
            want = (len(names), 2, len(ax["wires"]), len(ax["r_um"]), len(ax["D_mm"]), len(ax["tau_ns"]))
            if tuple(K.shape) != want:
                why = f"K.npy shape {tuple(K.shape)} != (basis, side, wires, r, D, tau) {want} from kernel.json"
    except Exception as e:  # noqa: BLE001 - reported
        why = f"{type(e).__name__}: {e}"
    if why:
        C.check("kernel slice", False, f"NOT generated: {why}", "WARN"); return
    idx = [names.index(b) for b in SLICE_BASES]
    jd = int(np.argmin(np.abs(np.asarray(ax["D_mm"], dtype=float))))     # the D node nearest 0 (cloud on the strip)
    dst = C.out / "kernel" / "K_slice_centralD.npz"
    np.savez_compressed(dst, K=np.asarray(K[idx][:, :, :, :, jd, :], dtype=np.float32), basis=np.array(SLICE_BASES),
                        sides=np.array(["U", "V"]), wires=np.array(ax["wires"]), r_um=np.array(ax["r_um"]),
                        tau_ns=np.array(ax["tau_ns"]), D_index=jd, D_mm=float(ax["D_mm"][jd]))
    C.wrote(dst, f"{kd / 'K.npy'} [bases {SLICE_BASES} (indices {idx}), D index {jd}]")
    C.check("kernel slice", True, f"bases {SLICE_BASES} at indices {idx}; D = {ax['D_mm'][jd]} mm; shape verified")


def pick_examples(rows, home):
    vals = []
    for r in rows:
        if r.get("accuracy_v1_primary") in ("PASS", "FAIL"):
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
        if s not in home:
            C.skip.append((s, "example slot not in the report")); continue
        V = home[s] / "muons" / f"muon_{s}"
        lay = V / "shaped_layout" / f"muon_{s}"
        files = [(f, Path(f.relative_to(lay))) for f in sorted(lay.rglob("*.csv"))] if lay.is_dir() else []
        fd = V / "far_dir.txt"
        if fd.is_file():
            E = Path(fd.read_text().strip())
            files += [(E / "S0_local" / x / y, Path("S0T1_full") / x / y)
                      for x, y in (("replay", "waveforms.csv"), ("prompt", "prompt_waveforms.csv"))]
        files += [(V / "K2A_prod" / x / y, Path("K2A_full") / x / y)
                  for x, y in (("replay", "waveforms.csv"), ("prompt", "prompt_waveforms.csv"))]
        tot = sum(f.stat().st_size for f, _ in files if f.is_file())
        if tot > cap:
            C.skip.append((s, f"example waveforms {tot / 1e6:.1f} MB above --max-example-mb")); continue
        for f, rel in files:
            C.copy(f, Path("examples") / f"muon_{s}" / rel, cap)
        done.append(s)
    return done


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out", type=Path); ap.add_argument("--analysis", type=Path, required=True)
    ap.add_argument("--pf", type=Path, required=True); ap.add_argument("--also-report", type=Path, nargs="*", default=[])
    ap.add_argument("--stress", type=Path); ap.add_argument("--queue", type=Path)
    ap.add_argument("--val-root", type=Path, default=Path.home() / "na6p" / "k2a_validation")
    ap.add_argument("--runs-glob", default="overnight_800/muon_*")
    ap.add_argument("--kernel", type=Path); ap.add_argument("--green", type=Path)
    ap.add_argument("--examples", default="auto"); ap.add_argument("--extra-examples", default="")
    ap.add_argument("--max-file-mb", type=float, default=25); ap.add_argument("--max-example-mb", type=float, default=60)
    a = ap.parse_args()
    out = a.out.expanduser().resolve()
    if out.exists():
        raise SystemExit(f"refusing: {out} exists")
    inputs = [a.analysis, a.pf, a.val_root, a.stress, a.kernel, a.green, *a.also_report, a.queue.parent if a.queue else None]
    for d in [Path(x).expanduser().resolve() for x in inputs if x]:
        if out == d or d in out.parents:
            raise SystemExit(f"refusing: the upload folder {out} lies inside the input {d}")
    inv, inv_rows, sealed = read_inventory(a.pf)                 # checks before anything is written
    rows = report_rows(a.analysis)
    if {r["slot"] for r in rows} & sealed:
        raise SystemExit("refusing: the report contains sealed slots")
    ex_req = [] if a.examples in ("auto", "none") else [f"{int(x):04d}" for x in a.examples.split(",") if x.strip()]
    ex_req += [f"{int(x):04d}" for x in a.extra_examples.split(",") if x.strip()]
    if set(ex_req) & sealed:
        raise SystemExit("refusing: an example slot is sealed")
    out.mkdir(parents=True)
    C = Collector(out, a.max_file_mb * 1e6)
    inv_slots = {r["slot"] for r in inv_rows if r["slot"] not in sealed}
    have = collect_preflight(C, a.pf, inv, inv_rows, sealed)
    globs = [g for g in a.runs_glob.split(",") if g.strip()]
    rows, home, cnt, acc = collect_report(C, a.analysis, "main", sealed, inv_slots, a.val_root, globs, have)
    extra = []
    for d in a.also_report:
        r2 = report_rows(d)
        if {r["slot"] for r in r2} & sealed:
            C.check(f"also-report {d.name}", False, "contains sealed slots: not copied"); continue
        for f in REPORT:
            C.copy(d / f, Path("reports") / d.name / f)
        extra.append(f"{d.name}: {len(r2)} slots, per_muon.csv sha256 {sha(d / 'per_muon.csv')[:16]}")
    shome = {}
    if a.stress:
        if (a.stress / "per_muon.csv").is_file():
            _, shome, _, _ = collect_report(C, a.stress, "stress", sealed, inv_slots, a.val_root, globs, have)
        else:
            C.check("stress report", False, f"no per_muon.csv in {a.stress}", "WARN")
    if a.queue:
        C.check("queue file", C.copy(a.queue, Path("queue") / a.queue.name), str(a.queue), "WARN")
    for d in ("wire_probe", "kernel_probe", "closure_sladder", "closure_0637", "domain_scan"):
        C.tree(a.val_root / d, Path("probes") / d, ["*.md", "*.json", "*.csv", "*.txt", "*/*.md", "*/*.json"])
    if a.kernel:
        kernel_slice(C, a.kernel, a.green)
    ex = [] if a.examples == "none" else (pick_examples(rows, home) if a.examples == "auto" else [])
    ex = list(dict.fromkeys(ex + ex_req))
    got = collect_examples(C, {**shome, **home}, ex, a.max_example_mb * 1e6)
    for name, data, head in (("MANIFEST.tsv", C.man, "source\tdestination\tbytes\tsha256"),
                             ("SKIPPED.tsv", C.skip, "item\treason"), ("CHECKS.tsv", C.checks, "check\tstatus\tdetail"),
                             ("MISSING_ESSENTIAL.tsv", C.missing, "slot\titem\treason")):
        with open(out / name, "w") as fh:
            fh.write(head + "\n"); fh.writelines("\t".join(map(str, t)) + "\n" for t in data)
    fails = [c for c in C.checks if c[1] == "FAIL"]
    status = "COMPLETE" if not C.missing and not fails else "INCOMPLETE"
    summ = [f"status {status}  ({len(C.missing)} essential items missing, {len(fails)} failed checks: see "
            "MISSING_ESSENTIAL.tsv, CHECKS.tsv)",
            f"upload folder {out}",
            f"report {a.analysis}: {len(rows)} slots; per_muon.csv sha256 {sha(a.analysis / 'per_muon.csv')[:16]}",
            "  execution " + ", ".join(f"{k or '-'} {v}" for k, v in sorted(cnt.items())),
            "  accuracy_v1_primary " + ", ".join(f"{k or '-'} {v}" for k, v in sorted(acc.items())),
            *(["also-report " + x for x in extra]),
            f"stress slots {len(shome)}; examples {got}",
            f"sealed slots in the inventory: {len(sealed)}, none copied",
            f"files {len(C.man)}, {sum(m[2] for m in C.man) / 1e6:.1f} MB; skipped {len(C.skip)}"]
    (out / "SUMMARY.txt").write_text("\n".join(summ) + "\n")
    print("\n".join(summ))
    return 0 if status == "COMPLETE" else 2


if __name__ == "__main__":
    sys.exit(main())
