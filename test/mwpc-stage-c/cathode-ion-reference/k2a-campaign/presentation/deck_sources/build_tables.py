#!/usr/bin/env python3
"""Join the collected validation data into analysis tables (read-only on the upload)."""
import glob
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

U = Path("/home/claude/pres/data/upload_20261128_143502")
OUT = Path("/home/claude/pres/an")


def load_outcomes(tag):
    rows = []
    for f in sorted(glob.glob(str(U / tag / "muons" / "muon_*" / "outcome.json"))):
        o = json.load(open(f))
        s = o["muon"]
        r = {"slot": s, "outcome": o["outcome"], "acc": (o.get("accuracy") or {}).get("status")}
        dom = o.get("k2a_domain") or {}
        r["k2a_clouds"] = dom.get("clouds"); r["k2a_groups"] = dom.get("groups")
        for comp, pref in (("S0T1_vs_K2A", "m"), ("S0T0_vs_K2A", "tot"), ("S0T0_vs_S0T1", "tim")):
            c = (o.get("comparisons") or {}).get(comp) or {}
            for k, v in c.items():
                if isinstance(v, (int, float)) and not isinstance(v, bool):
                    r[f"{pref}:{k}"] = v
        sh = o.get("shaped") or {}
        for run, var, pref in (("refS0T1_ionelec", "K2A", "sh"), ("refS0T1_ion", "K2A", "shion"),
                               ("refS0T0_ionelec", "K2A", "shT0"), ("refS0T1_ionelec", "S0T0", "shtim")):
            m = ((sh.get(run) or {}).get("metrics") or {}).get(var) or {}
            for tp, d in m.items():
                for k, v in d.items():
                    if isinstance(v, (int, float)) and not isinstance(v, bool):
                        r[f"{pref}{tp}:{k}"] = v
        rt = o.get("runtime") or {}
        for side in ("S0", "K2A"):
            for k, v in (rt.get(side) or {}).items():
                if isinstance(v, (int, float)) and not isinstance(v, bool):
                    r[f"rt{side}:{k}"] = v
        for k, v in (rt.get("speedup_measured") or {}).items():
            r[f"sp:{k}"] = v
        dr = (o.get("criteria_v1_dry_run") or {}).get("primary_vs_S0T1") or {}
        r["exemptions"] = ";".join(str(e) for e in dr.get("exemptions") or [])
        r["n_checks"] = dr.get("n_checks")
        trn = (o.get("truncation") or {}).get("planes") or {}
        for p in ("U", "V"):
            r[f"trunc_{p}_Qout_fC"] = (trn.get(p) or {}).get("Q_outside_window_fC")
        rows.append(r)
    return pd.DataFrame(rows)


def population():
    dom = pd.read_csv(U / "preflight" / "domain_per_event.csv")
    dom["slot"] = dom["event"].str.replace("muon_", "", regex=False)
    inv = pd.read_csv(U / "preflight" / "inventory_slots_nonsealed.csv", dtype={"slot": str})
    ev = []
    for f in glob.glob(str(U / "preflight" / "moments_nonsealed" / "*" / "muon_*.moments.csv.event.json")):
        e = json.load(open(f)); e["slot"] = os.path.basename(f)[5:9]; ev.append(e)
    ev = pd.DataFrame(ev)[["slot", "ion_rows_saved", "path_length_mm", "momentum_GeV"]]
    pop = inv[["slot", "role", "channel", "station", "root_bytes"]].merge(dom.drop(columns=["channel", "station", "role"]), on="slot", how="left")
    pop = pop.merge(ev, on="slot", how="left")
    # per-slot moments aggregates (all non-sealed events)
    agg = []
    for f in glob.glob(str(U / "preflight" / "moments_nonsealed" / "*" / "muon_*.moments.csv")):
        m = pd.read_csv(f)
        if m.empty:
            continue
        q = m["N_fC"].values
        agg.append({"slot": os.path.basename(f)[5:9], "m_clouds": len(m), "m_charge_fC": q.sum(),
                    "m_mean_abs_A": float((q * m["A"].abs()).sum() / q.sum()), "m_primaries": m["primary"].nunique(),
                    "m_wires": m["wire"].nunique()})
    pop = pop.merge(pd.DataFrame(agg), on="slot", how="left")
    q = pd.read_csv(U / "queue" / "queue.tsv", sep="\t", header=None, names=["slot", "q_rows", "q_station", "q_channel"],
                    dtype={"slot": str})
    q["queue_pos"] = np.arange(len(q))
    pop = pop.merge(q, on="slot", how="left")
    return pop


def main():
    pm = pd.read_csv(U / "main" / "report" / "per_muon.csv", dtype={"slot": str})
    oc = load_outcomes("main")
    pop = population()
    t = pm.merge(oc, on="slot", how="left").merge(pop.drop(columns=[c for c in ("channel", "station") if False]), on="slot",
                                                   how="left", suffixes=("", "_pop"))
    t["station_i"] = t["station_pop"].astype("Int64")
    t["channel_i"] = t["channel_pop"]
    acc = t[(t.accuracy_v1_primary == "PASS") | (t.accuracy_v1_primary == "FAIL")].copy()
    snap = pd.read_csv(U / "reports" / "analysis_20261128_115131" / "per_muon.csv", dtype={"slot": str})
    stress_pm = pd.read_csv(U / "stress" / "report" / "per_muon.csv", dtype={"slot": str})
    soc = load_outcomes("stress")
    stress = stress_pm.merge(soc, on="slot").merge(pop, on="slot", how="left", suffixes=("", "_pop"))
    t.to_pickle(OUT / "all.pkl"); acc.to_pickle(OUT / "acc.pkl"); pop.to_pickle(OUT / "pop.pkl")
    snap.to_pickle(OUT / "snap.pkl"); stress.to_pickle(OUT / "stress.pkl")
    print("all", t.shape, "acc", acc.shape, "pop", pop.shape, "snap", snap.shape, "stress", stress.shape)
    print("pop roles", pop.role.value_counts().to_dict())
    print("acc stations", acc.station_i.value_counts().sort_index().to_dict())
    print("acc channels", acc.channel_i.value_counts().to_dict())
    print("missing pop join in acc:", acc.station_i.isna().sum(), " missing sh:", acc["sh100:U_L2_rel"].isna().sum())


if __name__ == "__main__":
    main()
