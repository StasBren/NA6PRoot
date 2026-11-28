#!/usr/bin/env python3
"""Numbers quoted on the slides, computed from the joined tables (accuracy set = 192 VALID, non-OOD events)."""
import json

import numpy as np
import pandas as pd

A = pd.read_pickle("/home/claude/pres/an/acc.pkl")
T = pd.read_pickle("/home/claude/pres/an/all.pkl")
P = pd.read_pickle("/home/claude/pres/an/pop.pkl")
S = pd.read_pickle("/home/claude/pres/an/snap.pkl")
ST = pd.read_pickle("/home/claude/pres/an/stress.pkl")


def d(x):
    x = pd.to_numeric(x, errors="coerce").dropna().values
    return {"n": int(len(x)), "mean": float(x.mean()), "rms": float(np.sqrt((x ** 2).mean())), "std": float(x.std(ddof=1)) if len(x) > 1 else 0.0,
            "median": float(np.median(x)), "median_abs": float(np.median(np.abs(x))), "p90_abs": float(np.quantile(np.abs(x), .9)),
            "max_abs": float(np.abs(x).max()), "min": float(x.min()), "max": float(x.max())}


st = {"n_slots": len(T), "execution": T.execution.value_counts().to_dict(),
      "accuracy": T.accuracy_v1_primary.value_counts().to_dict(), "n_acc": len(A),
      "acc_roles": A.role.value_counts().to_dict(), "acc_station": {str(k): int(v) for k, v in A.station_i.value_counts().sort_index().items()},
      "acc_channel": A.channel_i.value_counts().to_dict(), "exemptions": A.loc[A.exemptions.astype(str) != "", ["slot", "exemptions"]].values.tolist()}
keys = {
    "sh100_U_L2": "sh100:U_L2_rel", "sh100_V_L2": "sh100:V_L2_rel", "sh200_U_L2": "sh200:U_L2_rel", "sh200_V_L2": "sh200:V_L2_rel",
    "shion100_U_L2": "shion100:U_L2_rel", "shion100_V_L2": "shion100:V_L2_rel",
    "unsh_ie_U": "m:ion_plus_electron_step.U.unshaped_L2_rel", "unsh_ie_V": "m:ion_plus_electron_step.V.unshaped_L2_rel",
    "unsh_ion_U": "m:ion_only.U.unshaped_L2_rel", "unsh_ion_V": "m:ion_only.V.unshaped_L2_rel",
    "Q_U": "m:ion_plus_electron_step.U.window_charge_rel_diff", "Q_V": "m:ion_plus_electron_step.V.window_charge_rel_diff",
    "UV": "sh100:UV_charge_ratio_rel", "du_100": "sh100:u_mm", "dw_100": "sh100:w_mm", "du_200": "sh200:u_mm", "dw_200": "sh200:w_mm",
    "du_ion100": "shion100:u_mm", "dw_ion100": "shion100:w_mm",
    "dt10_U": "sh100:U_t_f0.1_ns", "dt10_V": "sh100:V_t_f0.1_ns", "dt10_U200": "sh200:U_t_f0.1_ns", "dt10_V200": "sh200:V_t_f0.1_ns",
    "lead_U": "sh100:U_lead_peak_rel", "lead_V": "sh100:V_lead_peak_rel", "sec_U": "sh100:U_second_over_lead", "sec_V": "sh100:V_second_over_lead",
    "tim_sh100_U": "shtim100:U_L2_rel", "tim_unsh_ion_U": "tim:ion_only.U.unshaped_L2_rel",
    "tot_sh100_U": "shT0100:U_L2_rel",
}
st["dist"] = {k: d(A[v]) for k, v in keys.items()}
# decisions
dec = {}
for tp in ("100", "200"):
    for p in ("U", "V"):
        dec[f"lead_changed_{p}{tp}"] = int(A[f"sh{tp}:{p}_lead_index_changed"].sum())
        for f in ("0.05", "0.1"):
            dec[f"flips_{p}{tp}_f{f}"] = int(A[f"sh{tp}:{p}_hitflips_f{f}"].sum())
st["decisions"] = dec
# runtime over the accuracy set
rt = {"S0_T0_total": d(A["rtS0:total_T0_single_pass"]), "S0_T1_total": d(A["rtS0:total_T1_single_pass"]),
      "S0_rkf": d(A["rtS0:rkf_export_archive_read_plus_ion_drift"]), "S0_replay_T1": d(A["rtS0:weighting_replay_T1_chosen_attempt_sum"]),
      "S0_replay_T0": d(A["rtS0:weighting_replay_T0_chosen_attempt_sum"]), "S0_estep_T0": d(A["rtS0:electron_step_T0"]),
      "K2A_io": d(A["rtK2A:archive_io_read_validate_min"]), "K2A_mom": d(A["rtK2A:moments_construction_min"]),
      "K2A_eval": d(A["rtK2A:response_eval_warm_median"]), "K2A_cold_open": d(A["rtK2A:kernel_open_verify_once"]),
      "sp_eval_T0": d(A["sp:S0_T0/cathode_response_eval"]), "sp_mom_T0": d(A["sp:S0_T0/response_from_births_incl_moments"]),
      "sp_e2e_T0": d(A["sp:S0_T0/end_to_end_from_archive_incl_io"]), "sp_cold_T0": d(A["sp:S0_T0/single_event_cold_incl_kernel_open"]),
      "sp_eval_T1": d(A["sp:S0_T1/cathode_response_eval"]), "sp_e2e_T1": d(A["sp:S0_T1/end_to_end_from_archive_incl_io"]),
      "negative_substeps": {c: A.loc[A[c] < 0, "slot"].tolist() for c in A.columns if c.startswith("rtS0:") and (A[c] < 0).any()},
      "rkf_max_slot": A.loc[A["rtS0:rkf_export_archive_read_plus_ion_drift"].idxmax(), "slot"]}
st["runtime"] = rt
# complexity
st["clouds"] = d(A["k2a_clouds"]); st["charge"] = d(A["m_charge_fC"]); st["rows"] = d(A["ion_rows_saved"])
pp = P[P.m_clouds.notna()]
st["pop_n_with_moments"] = int(len(pp)); st["pop_rows"] = d(pp["ion_rows_saved"]); st["pop_clouds"] = d(pp["m_clouds"])
st["pop_charge"] = d(pp["m_charge_fC"])
q = P[P.queue_pos.notna()]
st["queue"] = {"n": int(len(q)), "processed_in_report": int(q.slot.isin(T.slot).sum()),
               "rows_processed_max": float(q[q.slot.isin(T.slot)].q_rows.max()), "rows_pending_median": float(q[~q.slot.isin(T.slot)].q_rows.median())}
# fraction of population with rows <= max processed
st["pop_frac_rows_le_maxacc"] = float((pp.ion_rows_saved <= A.ion_rows_saved.max()).mean())
# snapshot
sacc = S[S.accuracy_v1_primary == "PASS"]
st["snapshot"] = {"n": len(S), "execution": S.execution.value_counts().to_dict(), "acc": S.accuracy_v1_primary.value_counts().to_dict(),
                  "sh100_U_L2_median": float(sacc.model_U_shaped_L2_Tp100.median()), "sh100_V_L2_median": float(sacc.model_V_shaped_L2_Tp100.median()),
                  "sh100_U_L2_max": float(sacc.model_U_shaped_L2_Tp100.max()), "sh100_V_L2_max": float(sacc.model_V_shaped_L2_Tp100.max()),
                  "u_max_um": float(sacc.model_u_um_Tp100.abs().max()), "w_max_um": float(sacc.model_w_um_Tp100.abs().max()),
                  "sp_eval_median": float(sacc.speedup_cathode_response_eval_T0.median()), "sp_e2e_median": float(sacc.speedup_end_to_end_from_archive_T0.median())}
st["new_from_per_muon_csv"] = {"sh100_U_L2_median": float(A.model_U_shaped_L2_Tp100.median()), "sh100_V_L2_median": float(A.model_V_shaped_L2_Tp100.median()),
                               "u_max_um": float(A.model_u_um_Tp100.abs().max()), "w_max_um": float(A.model_w_um_Tp100.abs().max())}
st["stress"] = ST[["slot", "execution", "accuracy_v1_primary", "model_U_shaped_L2_Tp100", "model_V_shaped_L2_Tp100", "k2a_clouds",
                   "ion_rows_saved", "S0_T0_response_s", "K2A_response_eval_s"]].to_dict("records")
st["truncation_Qout_fC"] = {p: d(A[f"trunc_{p}_Qout_fC"]) for p in ("U", "V")}
json.dump(st, open("/home/claude/pres/an/stats.json", "w"), indent=1, default=str)
for k in ("execution", "accuracy", "acc_roles", "exemptions", "decisions", "queue", "pop_frac_rows_le_maxacc", "snapshot", "new_from_per_muon_csv", "stress"):
    print(k, st[k])
for k in ("sh100_U_L2", "sh100_V_L2", "shion100_U_L2", "unsh_ie_U", "unsh_ion_U", "du_100", "dw_100", "du_ion100", "dt10_U", "Q_U", "UV", "lead_U", "tim_sh100_U"):
    x = st["dist"][k]; print(f"{k:14s} mean {x['mean']:.3g} rms {x['rms']:.3g} std {x['std']:.3g} med {x['median']:.3g} p90|| {x['p90_abs']:.3g} max|| {x['max_abs']:.3g}")
for k, x in rt.items():
    if isinstance(x, dict) and "median" in x:
        print(f"{k:14s} median {x['median']:.4g} min {x['min']:.4g} max {x['max']:.4g}")
print("neg", rt["negative_substeps"], "rkfmax", rt["rkf_max_slot"])
print("clouds", st["clouds"]["min"], st["clouds"]["median"], st["clouds"]["max"], "rows acc", st["rows"]["min"], st["rows"]["max"], "pop rows median", st["pop_rows"]["median"], "pop clouds med", st["pop_clouds"]["median"], st["pop_clouds"]["max"])
