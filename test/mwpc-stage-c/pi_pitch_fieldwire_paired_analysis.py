#!/usr/bin/env python3
"""Read-only paired Garfield analysis of PI p=3/3.5/4 mm, +/- 75um field wires.

Requires the completed frozen FULL-PITCH microscopic Garfield pilot, not
the narrow smoke results. Within a pitch, map each independent Garfield
trial from anode-only to field-wire-on by exact (box_id,sample,random_seed)
and verified physical birth coordinates. This exposes selection bias between
conditional mean drift times and true within-electron transport changes.

Statistical uncertainty: fixed-stratum bootstrap over separately seeded
microscopic single-electron trials WITHIN EACH u/v box (N=20 / cell).
The same resampled indices are applied to both electrode configurations.
This represents LUT Monte-Carlo error for the specified UNIFORM BOX
population; does NOT represent uncertainty across GEANT muon parent events.

No new Garfield/Geant4/Heed simulation is performed.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import Counter
from pathlib import Path

import numpy as np

from pi_pitch_fieldwire_pilot import DU, NV, PITCHES, tag

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
DEFAULT_INPUT = REPO / "test_runs/mwpc_stage_c/pi_pitch_fieldwire_pilot_v1/pilot"
DEFAULT_OUTPUT = REPO / "test_runs/mwpc_stage_c/pi_pitch_fieldwire_paired_v1"
GATES = (50., 100.)
FIELDS = ("collected", "field_wire_collected", "drift_time_ns",
          "incidence_normal_deg", "status", "u_start_mm", "v_start_mm",
          "random_seed", "box_id", "sample", "wire_index")


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_map(path, pitch, samples):
    boxes = round(pitch / DU) * NV
    required = set(FIELDS) | {"iu", "iv"}
    table = {}
    status = Counter()
    if not path.is_file():
        raise FileNotFoundError(path)
    with path.open(newline="", encoding="utf-8") as file:
        reader = csv.DictReader(file)
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(f"Missing columns in {path}: {required-set(reader.fieldnames or [])}")
        for row in reader:
            ib, ix = int(row["box_id"]), int(row["sample"])
            iu, iv = int(row["iu"]), int(row["iv"])
            if (ib != iu * NV + iv or not (0 <= ib < boxes)
                    or not (0 <= iv < NV) or not (0 <= ix < samples)):
                raise ValueError(f"Bad u/v/sample grid: {path}")
            k = ib, ix
            if k in table:
                raise ValueError(f"Duplicate Garfield trial {k}")
            collected = int(row["collected"])
            fw = int(row["field_wire_collected"])
            if collected not in (0, 1) or fw not in (0, 1) or collected+fw > 1:
                raise ValueError("Mutually exclusive electrode collection broken")
            time = float(row["drift_time_ns"])
            angle = float(row["incidence_normal_deg"])
            if collected and math.isfinite(time) and time < 0:
                raise ValueError("Negative collected drift time")
            status[int(row["status"])] += 1
            table[k] = {
                "collected": bool(collected),
                "field": bool(fw),
                "time": time if collected and math.isfinite(time)
                          and time >= 0 else math.nan,
                "angle": angle if collected and math.isfinite(angle)
                           and 0 <= angle <= 180 else math.nan,
                "seed": int(row["random_seed"]),
                "u": float(row["u_start_mm"]),
                "v": float(row["v_start_mm"]),
                "status": int(row["status"]),
            }
    if len(table) != boxes*samples:
        raise ValueError(f"Incomplete map: {len(table)} / {boxes*samples}")
    return table, dict(status)


def match_maps(pitch, input_dir, field_diameter=75.):
    bare_dir = input_dir / tag(pitch, 0.)
    field_dir = input_dir / tag(pitch, field_diameter)
    a_config = json.loads((bare_dir/"frozen_config.json").read_text())
    b_config = json.loads((field_dir/"frozen_config.json").read_text())
    if a_config["pitch_mm"] != pitch or b_config["pitch_mm"] != pitch:
        raise ValueError("Wrong pitch config")
    if a_config["field_wire_diam_um"] != 0. or \
            b_config["field_wire_diam_um"] != field_diameter:
        raise ValueError("Wrong field-wire diameter")
    identity = ("gap_minus_mm", "gap_plus_mm", "anode_voltage_v",
                "anode_diameter_um", "gas", "magnetic_field_t",
                "n_u", "n_v", "du_mm", "dv_mm", "samples_per_box",
                "seed", "cpp_sha256", "binary_sha256")
    for key in identity:
        if a_config[key] != b_config[key]:
            raise ValueError(f"Unpaired physical controls: {key}")
    samples = int(a_config["samples_per_box"])
    if samples < 2 or a_config.get("mode") != "pilot" or \
            b_config.get("mode") != "pilot":
        raise ValueError("Paired analysis requires full-pitch pilot maps")
    a_path = bare_dir/"garfield_electrons.csv"
    b_path = field_dir/"garfield_electrons.csv"
    bare, sa = read_map(a_path, pitch, samples)
    field, sb = read_map(b_path, pitch, samples)
    if bare.keys() != field.keys():
        raise ValueError("Garfield trial pairs have inconsistent IDs")
    nbox = round(pitch/DU)*NV
    shape = (nbox, samples)
    keys = [(ib, j) for ib in range(nbox) for j in range(samples)]
    a_col = np.empty(shape, dtype=bool)
    b_col = np.empty(shape, dtype=bool)
    a_time = np.full(shape, np.nan)
    b_time = np.full(shape, np.nan)
    a_angle = np.full(shape, np.nan)
    b_angle = np.full(shape, np.nan)
    nfield = 0
    transitions = Counter()
    for (ib,j) in keys:
        a, b = bare[ib,j], field[ib,j]
        if a["seed"] != b["seed"] or \
           not math.isclose(a["u"], b["u"], rel_tol=0., abs_tol=1.e-8) or \
           not math.isclose(a["v"], b["v"], rel_tol=0., abs_tol=1.e-8):
            raise ValueError(f"Initial Garfield electron not paired at {ib}/{j}")
        a_col[ib,j] = a["collected"]
        b_col[ib,j] = b["collected"]
        a_time[ib,j] = a["time"]
        b_time[ib,j] = b["time"]
        a_angle[ib,j] = a["angle"]
        b_angle[ib,j] = b["angle"]
        nfield += b["field"]
        transitions[(int(a["collected"]), int(b["collected"]))] += 1
    if sum(transitions.values()) != a_col.size:
        raise ValueError("Paired Garfield transition accounting failed")
    return {
        "pitch": pitch, "samples": samples,
        "bare": a_col, "field": b_col,
        "bare_time": a_time, "field_time": b_time,
        "bare_angle": a_angle, "field_angle": b_angle,
        "transition_counts": {f"{a}->{b}": transitions[a,b]
                              for a in (0,1) for b in (0,1)},
        "status_bare": sa, "status_field": sb,
        "field_wire_ends": nfield,
        "input_sha256": {str(a_path):sha256(a_path),
                         str(b_path):sha256(b_path)},
    }


def moments(values):
    good = np.asarray(values, dtype=float)
    good = good[np.isfinite(good)]
    if not len(good):
        return math.nan, math.nan
    m = float(good.mean())
    return m, float(np.sqrt(np.mean((good-m)**2)))


def metrics(a, b, ta, tb, aa, ab):
    both = a & b
    diff = tb[both] - ta[both]
    diff = diff[np.isfinite(diff)]
    a_time = ta[a & np.isfinite(ta)]
    b_time = tb[b & np.isfinite(tb)]
    a_mean, a_rms = moments(a_time)
    b_mean, b_rms = moments(b_time)
    ai = aa[a & np.isfinite(aa)]
    bi = ab[b & np.isfinite(ab)]
    angle_a, angle_rms_a = moments(ai)
    angle_b, angle_rms_b = moments(bi)
    return {
        "anode_fraction_no": float(a.mean()),
        "anode_fraction_field": float(b.mean()),
        "delta_anode_fraction_pp": float(100*(b.mean()-a.mean())),
        "mean_time_no_ns": a_mean,
        "mean_time_field_ns": b_mean,
        "delta_conditional_mean_time_ns": b_mean-a_mean,
        "rms_time_no_ns": a_rms,
        "rms_time_field_ns": b_rms,
        "delta_conditional_rms_time_ns": b_rms-a_rms,
        "n_both_collected": int(np.count_nonzero(both)),
        "paired_mean_time_delta_both_ns":
            float(diff.mean()) if len(diff) else math.nan,
        "paired_rms_time_delta_both_ns":
            moments(diff)[1],
        "n_paired_valid_times": int(len(diff)),
        "angle_mean_no_deg": angle_a,
        "angle_mean_field_deg": angle_b,
        "delta_conditional_mean_angle_deg": angle_b-angle_a,
        "angle_rms_no_deg": angle_rms_a,
        "angle_rms_field_deg": angle_rms_b,
        "delta_conditional_rms_angle_deg": angle_rms_b-angle_rms_a,
        "anode_le50ns_fraction_no":
            float(np.mean(a & np.isfinite(ta) & (ta<=50))),
        "anode_le50ns_fraction_field":
            float(np.mean(b & np.isfinite(tb) & (tb<=50))),
        "anode_le100ns_fraction_no":
            float(np.mean(a & np.isfinite(ta) & (ta<=100))),
        "anode_le100ns_fraction_field":
            float(np.mean(b & np.isfinite(tb) & (tb<=100))),
    }


def analyse(pair, replicas=400, seed=260020):
    a, b = pair["bare"], pair["field"]
    ta, tb = pair["bare_time"], pair["field_time"]
    aa, ab = pair["bare_angle"], pair["field_angle"]
    nbox, samples = a.shape
    point = metrics(a,b,ta,tb,aa,ab)
    rng = np.random.default_rng(seed)
    row_index = np.arange(nbox)[:,None]
    variations = {k:[] for k in point if k not in
                  ("n_both_collected","n_paired_valid_times")}
    # Preserve exactly 20 Garfield trials per spatial cell in each bootstrap.
    for _ in range(replicas):
        chosen = rng.integers(0, samples, (nbox,samples))
        mm = metrics(a[row_index,chosen], b[row_index,chosen],
                     ta[row_index,chosen], tb[row_index,chosen],
                     aa[row_index,chosen], ab[row_index,chosen])
        for key in variations:
            variations[key].append(mm[key])
    result = {
        "pitch_mm": pair["pitch"],
        "n_spatial_cells": nbox,
        "trials_per_cell": samples,
        "n_electron_pairs": int(a.size),
        "field_wire_endpoints": pair["field_wire_ends"],
        "transitions": pair["transition_counts"],
        "raw_endpoint_status_no_field": pair["status_bare"],
        "raw_endpoint_status_with_field": pair["status_field"],
        "input_sha256": pair["input_sha256"],
        "point_estimate": point,
        "bootstrap_standard_errors": {
            k: float(np.std([v for v in values if math.isfinite(v)],
                            ddof=1)) if sum(math.isfinite(v) for v in values)>1
            else None for k,values in variations.items()
        },
        "uncertainty_unit": (
            "Single Garfield electron, resampled WITHIN fixed initial (u,v) "
            "cell and paired across field/no-field configurations. "
            "No parent-Geant4-event uncertainty or interpolation systematic."),
    }
    return result


def spatial_profiles(pair):
    n_u = round(pair["pitch"]/DU)
    result = []
    for iu in range(n_u):
        sl = slice(iu*NV, (iu+1)*NV)
        point = metrics(*(pair[k][sl] for k in (
            "bare", "field", "bare_time", "field_time",
            "bare_angle", "field_angle")))
        result.append({"pitch_mm":pair["pitch"],
                       "u_center_mm": -pair["pitch"]/2+(iu+.5)*DU,
                       "u_over_pitch":(-pair["pitch"]/2+(iu+.5)*DU)/pair["pitch"],
                       "u_bin":iu,
                       **point})
    return result


def run(input_dir, out, diameter=75., replicas=400):
    out.mkdir(parents=True, exist_ok=True)
    studies = []
    profile = []
    for pitch in PITCHES:
        pair = match_maps(pitch, input_dir, diameter)
        studies.append(analyse(pair, replicas=replicas, seed=260020+int(pitch*10)))
        profile.extend(spatial_profiles(pair))
    table = out/"pi_pitch_fieldwire_paired_u_profiles.csv"
    fields = list(profile[0])
    with table.open("w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=fields)
        writer.writeheader()
        writer.writerows(profile)
    output = {
        "scope": "EXISTING full-pitch 126k Garfield single-electron trials",
        "field_wire_diameter_um": diameter,
        "drift_model": "B=0, Ar/CO2 70:30, 3+3mm, +1800V anodes",
        "statistics": "Arithmetic mean, centered RMS; paired stratified "
                      "bootstrap Monte-Carlo errors across common initial "
                      "electron positions, NOT errors on correlated "
                      "source-parent GEANT muon events.",
        "warning": "Cross-pitch electron births cover different physical "
                   "u intervals; results are UNIFORM over respective unit "
                   "cells, NOT matched native Geant4/Heed tracks.",
        "field_wire_hits_zero_not_proof": "Electrons prefer +1800V anodes; "
                     "noncollected endpoint statuses still need inspection.",
        "groups":studies,
    }
    json_path=out/"pi_pitch_fieldwire_paired_summary.json"
    json_path.write_text(json.dumps(output, indent=2, allow_nan=False)+"\n",
                         encoding="utf-8")
    print("PASS: paired same-start Gaussian-bootstrap Garfield physics study")
    print("Output:",json_path)
    print("By initial u:",table)
    for row in studies:
        m=row["point_estimate"];e=row["bootstrap_standard_errors"]
        print(f"p={row['pitch_mm']:.1f}mm, "
              f"delta collection={m['delta_anode_fraction_pp']:.3f} +/- "
              f"{e['delta_anode_fraction_pp']:.3f} percentage points, "
              f"delta mean all collected={m['delta_conditional_mean_time_ns']:.3f}"
              f" +/- {e['delta_conditional_mean_time_ns']:.3f} ns, "
              f"delta mean SAME collected electrons="
              f"{m['paired_mean_time_delta_both_ns']:.3f} +/- "
              f"{e['paired_mean_time_delta_both_ns']:.3f} ns "
              f"(N={m['n_paired_valid_times']})")
        print("  collection transition counts:",row["transitions"])
        print("  status without field:",row["raw_endpoint_status_no_field"])
        print("  status with field   :",row["raw_endpoint_status_with_field"])
    return output


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input",type=Path,default=DEFAULT_INPUT)
    ap.add_argument("--out",type=Path,default=DEFAULT_OUTPUT)
    ap.add_argument("--field-wire-diam-um",type=float,default=75.)
    ap.add_argument("--bootstrap",type=int,default=400)
    args=ap.parse_args()
    if args.bootstrap<2:
        ap.error("--bootstrap must be >=2")
    run(args.input,args.out,args.field_wire_diam_um,args.bootstrap)


if __name__=="__main__":
    main()
