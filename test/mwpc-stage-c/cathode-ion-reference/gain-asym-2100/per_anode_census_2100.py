#!/usr/bin/env python3
"""Read-only exact per-anode Garfield++ census from the 2100 V source archive.

Hierarchy:
 - one row per (muon, Heed primary, physical/nearest-birth anode)
 - one row per (muon, physical/nearest-birth anode)
No K256 estimates, ion drift, induction, or re-simulation are used.
"""
from __future__ import annotations

from pathlib import Path
import re
import numpy as np
import pandas as pd

REQUIRED = (
    "primary_index", "wire_index_global", "ion_births", "collected_e",
    "ion_up", "ion_down", "ion_boundary", "A_ion_hemi",
    "mean_t_birth_ns", "sigma_t_birth_ns",
    "mean_t_collected_ns", "sigma_t_collected_ns",
    "sigma_w_birth_mm",
)
NUMERIC = (*REQUIRED,)
LABEL = re.compile(r"^muon_(\d+)$")


def require(ok: bool, reason: str):
    if not bool(ok):
        raise ValueError("2100 V original anode census: " + reason)


def verified_join(raw: pd.DataFrame, pr: pd.DataFrame, mu: pd.DataFrame):
    """Strict join against already ROOT-audited V0 primary and muon counters."""
    ids = ("muon_slot", "primary")
    actual = raw.groupby(list(ids), as_index=False)[["ion_births", "collected_e"]].sum()
    expected = pr[["muon_slot", "primary", "ion_births_exact", "collected_electrons"]].copy()
    checked = expected.merge(actual, on=list(ids), how="left",
                             validate="one_to_one")
    checked[["ion_births", "collected_e"]] = checked[
        ["ion_births", "collected_e"]].fillna(0)
    require(np.array_equal(checked["ion_births"].to_numpy(dtype=np.int64),
                           checked["ion_births_exact"].to_numpy(dtype=np.int64)),
            "per-primary exact ion-birth sum disagrees with archived ROOT count")
    require(np.array_equal(checked["collected_e"].to_numpy(dtype=np.int64),
                           checked["collected_electrons"].to_numpy(dtype=np.int64)),
            "per-primary physical collection sum disagrees with archived ROOT count")
    event = (actual.groupby("muon_slot", as_index=False)[["ion_births","collected_e"]]
             .sum().rename(columns={"collected_e":"collected_e_anode"}))
    ev = mu[["muon_slot","ion_births_exact","collected_electrons"]].merge(
        event,on="muon_slot",how="left",validate="one_to_one")
    ev[["ion_births","collected_e_anode"]] = ev[
        ["ion_births","collected_e_anode"]].fillna(0)
    require(np.array_equal(ev["ion_births"].to_numpy(np.int64),
                           ev["ion_births_exact"].to_numpy(np.int64)),
            "per-muon ion-birth total disagrees with original ROOT archive")
    require(np.array_equal(ev["collected_e_anode"].to_numpy(np.int64),
                           ev["collected_electrons"].to_numpy(np.int64)),
            "per-muon collected-electron total disagrees with original ROOT archive")


def read_per_anode(run_dir: Path, primaries: pd.DataFrame,
                   muons: pd.DataFrame):
    """Returns (fragment, per_muon_anode) dataframes after strict census closure.

    The original immutable wire-census inputs are located using rootfiles.txt.
    Selected5 and population997 each retain their own exact input-root list.
    """
    rootlist = Path(run_dir) / "rootfiles.txt"
    require(rootlist.is_file(), f"missing source manifest: {rootlist}")
    paths = [Path(x.strip()) for x in rootlist.read_text().splitlines() if x.strip()]
    require(len(paths) == len(set(paths)), "duplicate source ROOT path")
    require(len(paths) == len(muons), "rootfiles.txt event count differs from V0 output")
    frames = []
    for root in paths:
        match = LABEL.fullmatch(root.parent.name)
        require(bool(match), f"unrecognized muon source path: {root}")
        slot = int(match.group(1))
        p = root.parent / "microscopic_wire_seed.csv"
        require(p.is_file(), f"missing immutable original wire census: {p}")
        df = pd.read_csv(p, low_memory=False)
        missing=set(REQUIRED)-set(df.columns)
        require(not missing, f"missing wire columns in {p}: {sorted(missing)}")
        if df.empty:
            continue
        df = df.loc[:,list(REQUIRED)].copy()
        df.insert(0,"muon_slot",slot)
        frames.append(df)
    require(frames,"all original per-wire census CSVs are empty")
    raw=pd.concat(frames,ignore_index=True)
    raw=raw.rename(columns={"primary_index":"primary",
                           "wire_index_global":"global_wire",
                           "ion_births":"ion_births",
                           "collected_e":"collected_e"})
    # Require whole-number original microscopic counts, not rounded estimates.
    for c in ("muon_slot","primary","global_wire","ion_births",
              "collected_e","ion_up","ion_down","ion_boundary"):
        raw[c]=pd.to_numeric(raw[c],errors="raise")
        require(bool(raw[c].notna().all()),f"null exact counter {c}")
        require(bool(np.isfinite(raw[c].to_numpy(dtype=float)).all()),
                f"nonfinite exact counter {c}")
        require(bool(np.all(raw[c].to_numpy(dtype=float)==
                            np.floor(raw[c].to_numpy(dtype=float)))),
                f"fractional original integer counter {c}")
        raw[c]=raw[c].astype("int64")
    require(not raw.duplicated(["muon_slot","primary","global_wire"]).any(),
            "duplicate primary-anode records")
    for c in ("ion_births","collected_e","ion_up","ion_down","ion_boundary"):
        require(bool((raw[c]>=0).all()),f"negative original counter {c}")
    require(bool(((raw["ion_up"]+raw["ion_down"]+raw["ion_boundary"])
                  ==raw["ion_births"]).all()),
            "up+down+boundary birth counters fail exact closure")
    denom=raw["ion_up"]+raw["ion_down"]
    raw["A_UD_exact"]=np.divide(
        (raw["ion_up"]-raw["ion_down"]).to_numpy(dtype=float),
        denom.to_numpy(dtype=float),
        out=np.full(len(raw),np.nan),where=denom.to_numpy()>0)
    a_source=pd.to_numeric(raw["A_ion_hemi"],errors="coerce")
    mask=denom>0
    require(bool(np.allclose(raw.loc[mask,"A_UD_exact"],
                             a_source.loc[mask],rtol=1.e-5,atol=1.e-5)),
            "hemispheric A differs from source A_ion_hemi")
    # Source may encode undefined A as a blank/NaN or zero on a per-row basis.
    # Reject any finite nonzero A if there were no up/down births.
    require(bool((a_source.loc[~mask].isna() |
                  (a_source.loc[~mask].abs()<1e-12)).all()),
            "source A values defined where up+down=0")
    verified_join(raw,primaries,muons)
    t0=primaries[["muon_slot","primary","t_primary_ns",
                  "ion_births_exact","collected_electrons"]].copy()
    require(not t0.duplicated(["muon_slot","primary"]).any(),
            "duplicate V0 primary ID")
    raw=raw.merge(t0,on=["muon_slot","primary"],how="left",
                  validate="many_to_one")
    require(bool(raw["t_primary_ns"].notna().all()),
            "some anode fragments do not match an original Heed start")
    raw=raw.merge(muons[["muon_slot","channel","station","det_id",
                          "theta_u_deg","n_primaries","ion_births_exact",
                          "collected_electrons"]]
                  .rename(columns={"ion_births_exact":"muon_ion_births",
                                   "collected_electrons":"muon_collected_e"}),
                  on="muon_slot",how="left",validate="many_to_one")
    # mean_t_birth is the mean timestamp of avalanche-ion creation, not
    # the earliest ion timestamp and NOT the later ion-drift arrival time.
    for c in ("mean_t_birth_ns","sigma_t_birth_ns",
              "mean_t_collected_ns","sigma_t_collected_ns",
              "sigma_w_birth_mm"):
        raw[c]=pd.to_numeric(raw[c],errors="coerce")
    raw["ion_birth_delay_ns"]=raw["mean_t_birth_ns"]-raw["t_primary_ns"]
    raw.loc[raw["ion_births"]==0,"ion_birth_delay_ns"]=np.nan
    raw["electron_collection_delay_ns"]=raw["mean_t_collected_ns"]-raw["t_primary_ns"]
    raw.loc[raw["collected_e"]==0,"electron_collection_delay_ns"]=np.nan
    require(bool(np.isfinite(raw.loc[raw["ion_births"]>0,"ion_birth_delay_ns"]).all()),
            "missing positive-ion mean-birth time on positive-ion source row")
    require(bool(np.isfinite(raw.loc[raw["ion_births"]>0,"sigma_t_birth_ns"]).all()),
            "missing within-fragment ion-birth sigma_t")
    require(bool(np.isfinite(raw.loc[raw["collected_e"]>0,
                                     "electron_collection_delay_ns"]).all()),
            "missing electron collection time where electrons were collected")
    raw["absA_UD_exact"]=raw["A_UD_exact"].abs()
    raw["G_ion_primary"]=raw["ion_births_exact"]
    raw["ion_share_muon"]=np.divide(raw["ion_births"],raw["muon_ion_births"],
          out=np.full(len(raw),np.nan),where=raw["muon_ion_births"]>0)
    raw["electron_share_muon"]=np.divide(raw["collected_e"],raw["muon_collected_e"],
          out=np.full(len(raw),np.nan),where=raw["muon_collected_e"]>0)
    raw["G_ion_anode_contribution"]=np.divide(
          raw["ion_births"],raw["n_primaries"],
          out=np.full(len(raw),np.nan),where=raw["n_primaries"]>0)

    # Exact aggregation by one original muon and global physical wire ID.
    # Primary origins differ, so variance below concerns (ion birth time -
    # its OWN primary start), including both intra- and inter-primary spread.
    raw["birth_weighted_delay"]=raw["ion_births"]*raw["ion_birth_delay_ns"].fillna(0.)
    raw["birth_second_moment"]=raw["ion_births"]*(
         raw["ion_birth_delay_ns"].fillna(0.)**2+
         raw["sigma_t_birth_ns"].fillna(0.)**2)
    raw["collection_weighted_delay"]=raw["collected_e"]*(
         raw["electron_collection_delay_ns"].fillna(0.))
    raw["collection_second_moment"]=raw["collected_e"]*(
         raw["electron_collection_delay_ns"].fillna(0.)**2+
         raw["sigma_t_collected_ns"].fillna(0.)**2)
    grouped=raw.groupby(["muon_slot","global_wire"],as_index=False).agg(
         ion_births=("ion_births","sum"),
         collected_e=("collected_e","sum"),
         ion_up=("ion_up","sum"),ion_down=("ion_down","sum"),
         ion_boundary=("ion_boundary","sum"),
         n_primaries_at_anode=("primary","nunique"),
         birth_weighted_delay=("birth_weighted_delay","sum"),
         birth_second_moment=("birth_second_moment","sum"),
         collection_weighted_delay=("collection_weighted_delay","sum"),
         collection_second_moment=("collection_second_moment","sum"))
    total_ud=(grouped["ion_up"]+grouped["ion_down"]).to_numpy(float)
    grouped["A_UD_exact"]=np.divide(
         (grouped["ion_up"]-grouped["ion_down"]).to_numpy(float),
         total_ud,out=np.full(len(grouped),np.nan),where=total_ud>0)
    grouped["absA_UD_exact"]=grouped["A_UD_exact"].abs()
    n_ion=grouped["ion_births"].to_numpy(dtype=float)
    n_col=grouped["collected_e"].to_numpy(dtype=float)
    grouped["ion_birth_delay_ns"]=np.divide(
         grouped["birth_weighted_delay"],n_ion,
         out=np.full(len(grouped),np.nan),where=n_ion>0)
    grouped["electron_collection_delay_ns"]=np.divide(
         grouped["collection_weighted_delay"],n_col,
         out=np.full(len(grouped),np.nan),where=n_col>0)
    m_ion=grouped["ion_birth_delay_ns"].fillna(0.).to_numpy(float)
    m_col=grouped["electron_collection_delay_ns"].fillna(0.).to_numpy(float)
    second_ion=np.divide(grouped["birth_second_moment"],n_ion,
         out=np.full(len(grouped),np.nan),where=n_ion>0)
    second_col=np.divide(grouped["collection_second_moment"],n_col,
         out=np.full(len(grouped),np.nan),where=n_col>0)
    grouped["sigma_t_birth_relative_ns"]=np.sqrt(
         np.maximum(0.,second_ion-m_ion*m_ion))
    grouped["sigma_t_collected_relative_ns"]=np.sqrt(
         np.maximum(0.,second_col-m_col*m_col))
    grouped=grouped.merge(
         muons[["muon_slot","channel","station","det_id","theta_u_deg",
                "n_primaries","ion_births_exact","collected_electrons"]],
         on="muon_slot",how="left",validate="many_to_one")
    grouped["ion_share_muon"]=np.divide(grouped["ion_births"],grouped["ion_births_exact"],
         out=np.full(len(grouped),np.nan),where=grouped["ion_births_exact"]>0)
    grouped["electron_share_muon"]=np.divide(grouped["collected_e"],
          grouped["collected_electrons"],out=np.full(len(grouped),np.nan),
          where=grouped["collected_electrons"]>0)
    grouped["G_ion_anode_contribution"]=np.divide(
         grouped["ion_births"],grouped["n_primaries"],
         out=np.full(len(grouped),np.nan),where=grouped["n_primaries"]>0)
    grouped["n_primary_fragments"]=grouped["n_primaries_at_anode"]
    # These records are exact origin-corrected statistics for the immutable
    # Garfield dataset, not cathode charge or ion arrival-time estimates.
    return raw,grouped
