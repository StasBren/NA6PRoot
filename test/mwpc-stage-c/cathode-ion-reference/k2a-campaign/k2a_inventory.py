#!/usr/bin/env python3
"""Phase 1a, READ-ONLY: check the 1000-slot 2100 V K256 archive against the existing audited inventory.

  k2a_inventory.py --archive ARCHIVE --selection SELECTION_DIR --out NEW_DIR [--sha-jobs N | --no-sha]
                   [--reference-run RUN]

ARCHIVE        overnight_2100_weighted_k256_v1 (muon_0000 ... muon_0999); nothing in it is written or moved
SELECTION_DIR  phase_space_2100_signal40_v1: archive_inventory.csv (997 rows), excluded_slots.txt (3),
               development20.csv, sealed_validation20.csv, previously_inspected_not_sealed.txt, FROZEN_SOURCE_SHA256.txt
--reference-run  the validated S-ladder run: its muon_* with S0 references are marked role dev_reference

Checks, per slot 0000..0999:
  seal files present (complete.txt, SHA256SUMS, root_audit.log, ion_births.root, microscopic_wire_muon.csv), ROOT size,
  full `sha256sum -c SHA256SUMS` (default; --no-sha skips it and says so), membership in archive_inventory.csv and
  excluded_slots.txt; the two lists must partition 0..999 and agree with the seal files.
  Selection files: FROZEN_SOURCE_SHA256.txt re-verified. Sealed set: ONLY the index column of sealed_validation20.csv is
  read (no other column, no selected40.csv).
Roles: sealed_validation > dev_reference (S0 exists) > development > previously_inspected > pool.
Outputs (NEW_DIR, refuses to overwrite): slots.csv, excluded_slots_detail.md, events_nonsealed.txt, events_sealed.txt,
events_development.txt (batch lists "muon_NNNN PATH/ion_births.root"), inventory_check.json. Exit 1 if inconsistent.
"""
import argparse
import csv
import hashlib
import json
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

SEAL = ("complete.txt", "SHA256SUMS", "root_audit.log", "ion_births.root", "microscopic_wire_muon.csv")


def index_col(path):
    with open(path, newline="") as f:
        return [int(r["index"]) for r in csv.DictReader(f) if r.get("index", "").strip()]


def lines(path):
    return [int(x) for x in Path(path).read_text().split() if x.strip()] if Path(path).exists() else None


def sha256(p, chunk=1 << 24):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(chunk), b""):
            h.update(b)
    return h.hexdigest()


def sha_check(d):
    if not (d / "SHA256SUMS").is_file():
        return "no_SHA256SUMS"
    r = subprocess.run(["sha256sum", "-c", "--quiet", "SHA256SUMS"], cwd=d, capture_output=True, text=True)
    return "ok" if r.returncode == 0 else "FAIL: " + (r.stdout + r.stderr).strip().replace("\n", " | ")[:300]


def describe_dir(d, maxfiles=40):
    if not d.exists():
        return ["directory absent"]
    out = []
    for p in sorted(d.iterdir())[:maxfiles]:
        out.append(f"- `{p.name}` {'dir' if p.is_dir() else str(p.stat().st_size) + ' B'}")
        if p.is_file() and p.stat().st_size <= 4096 and p.suffix in (".txt", ".log", ".err", ".out", ".json", ""):
            try:
                head = p.read_text(errors="replace").splitlines()[:6]
            except OSError:
                head = []
            out += ["      " + h[:160] for h in head]
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--archive", type=Path, required=True); ap.add_argument("--selection", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True); ap.add_argument("--reference-run", type=Path)
    ap.add_argument("--sha-jobs", type=int, default=2); ap.add_argument("--no-sha", action="store_true")
    ap.add_argument("--slots", type=int, default=1000)
    a = ap.parse_args()
    if a.out.exists():
        raise SystemExit(f"refusing: {a.out} exists")
    S = a.selection
    problems = []
    # frozen selection files
    frozen = {}
    fs = S / "FROZEN_SOURCE_SHA256.txt"
    if fs.exists():
        for l in fs.read_text().splitlines():
            if l.strip():
                h, name = l.split(None, 1); name = name.strip().lstrip("*")
                frozen[name] = (h == sha256(S / name)) if (S / name).exists() else False
        if not all(frozen.values()):
            problems.append("selection files differ from FROZEN_SOURCE_SHA256.txt: " + ",".join(k for k, v in frozen.items() if not v))
    else:
        problems.append("FROZEN_SOURCE_SHA256.txt missing")
    with open(S / "archive_inventory.csv", newline="") as f:
        inv = {int(r["index"]): r for r in csv.DictReader(f)}
    excl = lines(S / "excluded_slots.txt") or []
    dev = set(index_col(S / "development20.csv"))
    sealed = set(index_col(S / "sealed_validation20.csv"))          # index column only
    prev = set(lines(S / "previously_inspected_not_sealed.txt") or [])
    ref = set()
    if a.reference_run:
        ref = {int(p.name[5:]) for p in a.reference_run.glob("muon_*") if (p / "S0" / "replay" / "waveforms.csv").exists()}
    if dev & sealed:
        problems.append(f"development and sealed sets overlap: {sorted(dev & sealed)}")
    if ref & sealed:
        problems.append(f"reference-run muons in the sealed set: {sorted(ref & sealed)}")
    if not ref <= dev:
        problems.append(f"reference-run muons not in development20: {sorted(ref - dev)}")
    slots = range(a.slots)
    if set(inv) | set(excl) != set(slots) or set(inv) & set(excl):
        problems.append("archive_inventory.csv and excluded_slots.txt do not partition the slots")
    a.out.mkdir(parents=True)
    t0 = time.time()
    rows = []
    for i in slots:
        d = a.archive / f"muon_{i:04d}"
        pres = {f: (d / f).is_file() for f in SEAL}
        role = ("sealed_validation" if i in sealed else "dev_reference" if i in ref else "development" if i in dev
                else "previously_inspected" if i in prev else "pool")
        rows.append({"slot": f"{i:04d}", "dir": d.is_dir(), **{f"has_{k.split('.')[0]}": v for k, v in pres.items()},
                     "root_bytes": (d / "ion_births.root").stat().st_size if pres["ion_births.root"] else 0,
                     "in_inventory": i in inv, "in_excluded_list": i in excl, "role": role,
                     "channel": inv[i]["channel"] if i in inv else "", "station": inv[i]["station"] if i in inv else "",
                     "sha_check": "not_run"})
    sealed_ok = [r for r in rows if all(r[f"has_{k.split('.')[0]}"] for k in SEAL)]
    if not a.no_sha:
        with ThreadPoolExecutor(max(1, a.sha_jobs)) as ex:
            for r, res in zip(sealed_ok, ex.map(lambda r: sha_check(a.archive / f"muon_{r['slot']}"), sealed_ok)):
                r["sha_check"] = res
    for r in rows:
        allseal = all(r[f"has_{k.split('.')[0]}"] for k in SEAL)
        shaok = r["sha_check"] in ("ok", "not_run")
        if r["in_inventory"] and allseal and shaok:
            r["status"] = "USABLE"
        elif r["in_excluded_list"] and not (allseal and r["sha_check"] == "ok"):
            r["status"] = "EXCLUDED"
        else:
            r["status"] = "INCONSISTENT"
            problems.append(f"slot {r['slot']}: inventory={r['in_inventory']} excluded={r['in_excluded_list']} "
                            f"seal={allseal} sha={r['sha_check']}")
    with open(a.out / "slots.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n")
        w.writeheader(); w.writerows(rows)
    L = ["# Excluded slots (read-only description; nothing was regenerated)", ""]
    for r in rows:
        if r["status"] != "USABLE":
            d = a.archive / f"muon_{r['slot']}"
            missing = [k for k in SEAL if not r[f"has_{k.split('.')[0]}"]]
            L += [f"## muon_{r['slot']}  status {r['status']}", "",
                  f"missing seal files: {', '.join(missing) or 'none'}; sha: {r['sha_check']}; listed in excluded_slots.txt: "
                  f"{r['in_excluded_list']}", ""] + describe_dir(d) + [""]
    (a.out / "excluded_slots_detail.md").write_text("\n".join(L) + "\n")
    for name, sel in (("events_nonsealed.txt", lambda r: r["role"] != "sealed_validation"),
                      ("events_sealed.txt", lambda r: r["role"] == "sealed_validation"),
                      ("events_development.txt", lambda r: r["role"] in ("development", "dev_reference"))):
        with open(a.out / name, "w") as f:
            for r in rows:
                if r["status"] == "USABLE" and sel(r):
                    f.write(f"muon_{r['slot']} {a.archive / ('muon_' + r['slot']) / 'ion_births.root'}\n")
    usable = [r for r in rows if r["status"] == "USABLE"]
    count = lambda key: {k: sum(1 for r in usable if r[key] == k) for k in sorted({r[key] for r in usable})}
    summ = {"tool": "k2a_inventory.py", "archive": str(a.archive), "selection": str(S), "slots": a.slots,
            "usable": len(usable), "excluded": [r["slot"] for r in rows if r["status"] == "EXCLUDED"],
            "inconsistent": [r["slot"] for r in rows if r["status"] == "INCONSISTENT"],
            "sha_verified": None if a.no_sha else sum(1 for r in usable if r["sha_check"] == "ok"),
            "sha_mode": "not run (--no-sha)" if a.no_sha else "sha256sum -c SHA256SUMS on every sealed slot",
            "by_role": count("role"), "by_channel": count("channel"),
            "selection_frozen_sha_ok": frozen, "inventory_sha256": sha256(S / "archive_inventory.csv"),
            "n_development": len(dev), "n_sealed": len(sealed), "n_reference_run": len(ref),
            "reference_run_muons": sorted(f"{i:04d}" for i in ref),
            "total_root_bytes": sum(r["root_bytes"] for r in usable), "wall_s": time.time() - t0, "problems": problems,
            "pass": not problems}
    (a.out / "inventory_check.json").write_text(json.dumps(summ, indent=1) + "\n")
    print(json.dumps({k: summ[k] for k in ("usable", "excluded", "inconsistent", "sha_verified", "by_role", "by_channel",
                                           "total_root_bytes", "wall_s")}, indent=1))
    for p in problems:
        print("PROBLEM", p)
    print(("PASS" if not problems else "FAIL") + f" inventory: {len(usable)} usable, excluded {summ['excluded']}")
    return 0 if not problems else 1


if __name__ == "__main__":
    raise SystemExit(main())
