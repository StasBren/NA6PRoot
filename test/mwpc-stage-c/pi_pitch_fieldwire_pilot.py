#!/usr/bin/env python3
"""PI wire-pitch/grounded-field-wire microscopic Garfield parameter pilot.

This is a NEW, independent Garfield-only study, NOT the old p=4 mm LUT.
For every p = 3,3.5,4 mm use a centered v=-3..+3 mm chamber and a
zero-volt field wire halfway between each pair of +1800 V anodes (optional).
The simulated electron birth sites are uniformly sampled within each box
of one pitch; they are NOT conditioned on native GEANT muon angle yet.

Usage (after building mwpc_stage_c_transport_box_map):
  python3 pi_pitch_fieldwire_pilot.py --mode smoke --samples 2 \
      --field-wire-diam-um 75
  python3 pi_pitch_fieldwire_pilot.py --mode pilot --samples 20 \
      --field-wire-diam-um 75

WARNING: 75 um is the CHOSEN study parameter, not a verified
NA60+ mechanical field-wire diameter. All 6 configurations are independent
and have isolated, SHA-pinned resumable output directories.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import subprocess
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
ROOT = REPO / "test_runs/mwpc_stage_c/pi_pitch_fieldwire_pilot_v1"
CPP = HERE / "c_transport_box_map.cxx"
DEFAULT_BINARY = HERE / "build/mwpc_stage_c_transport_box_map"
PITCHES = (3.0, 3.5, 4.0)
DU, DV = 0.1, 0.2
NV = 30
GAPS = (3., 3.)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def same_physics_manifest(old, current):
    """Source-only reporting changes do not alter microscopic Garfield data."""
    return ({k: v for k, v in old.items() if k != "runner_sha256"} ==
            {k: v for k, v in current.items() if k != "runner_sha256"})


def tag(pitch, fw):
    return f"p{pitch:.1f}".replace(".", "p") + (
        f"_field{fw:.1f}um".replace(".", "p") if fw else "_anode_only")


def expected_rows(pitch, samples, mode):
    nu = round(pitch / DU)
    return (1 if mode == "smoke" else nu) * NV * samples


def validate_column(path, iu, pitch, samples):
    n_u = round(pitch/DU)
    if not path.is_file():
        raise ValueError(f"Missing partial map column: {path}")
    seen = set()
    n = 0
    with path.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        required = {
            "iu", "iv", "box_id", "sample", "u_center_mm",
            "v_center_mm", "collected", "drift_time_ns",
            "field_wire_collected", "field_wire_index",
            "phi_hit_deg", "incidence_normal_deg",
        }
        if not required.issubset(reader.fieldnames or ()):
            raise ValueError("Missing pitch/field-wire Garfield columns")
        for row in reader:
            ib, j, k, sample = (int(row[key]) for key in (
                "iu", "iv", "box_id", "sample"))
            if ib != iu or not (0 <= ib < n_u) or not (0 <= j < NV) \
                or k != ib*NV+j or not (0 <= sample < samples):
                raise ValueError(f"Wrong pitch-grid bin in {path}")
            if not (math.isclose(float(row["u_center_mm"]),
                                 -pitch/2+(ib+.5)*DU, abs_tol=1.e-7) and
                    math.isclose(float(row["v_center_mm"]),
                                 -3+(j+.5)*DV, abs_tol=1.e-7)):
                raise ValueError("Garfield pitch-grid coordinate mismatch")
            c, fw = int(row["collected"]), int(row["field_wire_collected"])
            if c not in (0, 1) or fw not in (0, 1) or c+fw > 1:
                raise ValueError("Invalid anode/field-wire classification")
            key = j, sample
            if key in seen:
                raise ValueError("Repeated Garfield trial")
            seen.add(key)
            n += 1
    if n != NV*samples or len(seen) != n:
        raise ValueError(f"Incomplete Garfield column {iu}: {n}/{NV*samples}")
    return n


def load_column_timing(dest):
    """Read a measured Garfield runtime, never infer it from file mtimes."""
    path = dest.with_suffix(".timing.json")
    if not path.is_file():
        return None
    value = json.loads(path.read_text(encoding="utf-8"))
    if (value.get("csv_sha256") != digest(dest) or
            value.get("timing_kind") != "monotonic_wall_seconds" or
            not math.isfinite(float(value["seconds"])) or
            float(value["seconds"]) <= 0):
        raise ValueError(f"Inconsistent or invalid measured timing: {path}")
    return float(value["seconds"])


def save_column_timing(dest, seconds, n_electrons):
    if not math.isfinite(seconds) or seconds <= 0 or n_electrons <= 0:
        raise ValueError("Invalid Garfield measured runtime")
    info = {
        "timing_kind": "monotonic_wall_seconds",
        "seconds": float(seconds),
        "n_electrons": int(n_electrons),
        "csv_sha256": digest(dest),
        "measurement": "subprocess.run Garfield including process startup",
    }
    sidecar = dest.with_suffix(".timing.json")
    temp = sidecar.with_name(sidecar.name + ".tmp")
    temp.write_text(json.dumps(info, indent=2)+"\n", encoding="utf-8")
    os.replace(temp, sidecar)


def mean_rms(values, seed, reps=300):
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x)]
    n = len(x)
    if n == 0:
        return dict(n=0, mean=None, rms=None,
                    mean_error=None, rms_error=None)
    mean = float(x.mean())
    rms = float(np.sqrt(np.mean((x-mean)**2)))
    if n < 2:
        return dict(n=n, mean=mean, rms=rms,
                    mean_error=None, rms_error=None)
    rng = np.random.default_rng(seed)
    bm = np.empty(reps)
    br = np.empty(reps)
    for k in range(reps):
        sample = x[rng.integers(0, n, size=n)]
        bm[k] = sample.mean()
        br[k] = np.sqrt(np.mean((sample-bm[k])**2))
    return dict(n=n, mean=mean, rms=rms,
                mean_error=float(bm.std(ddof=1)),
                rms_error=float(br.std(ddof=1)))


def summarise(csv_path, pitch, fw, mode, samples):
    times, incidence = [], []
    total = collected = grounded = time_unknown = 0
    with csv_path.open(newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            total += 1
            c = int(row["collected"])
            g = int(row["field_wire_collected"])
            if c+g > 1:
                raise ValueError("Electron collected on anode and field wire")
            if g:
                if fw <= 0. or not (-4 <= int(row["field_wire_index"]) < 4):
                    raise ValueError("Field-wire hit absent from declared geometry")
            collected += c
            grounded += g
            if c:
                dt = float(row["drift_time_ns"])
                if math.isfinite(dt) and dt >= 0:
                    times.append(dt)
                else:
                    time_unknown += 1
                a = float(row["incidence_normal_deg"])
                if math.isfinite(a):
                    if not (0 <= a <= 180):
                        raise ValueError("Unphysical incidence angle")
                    incidence.append(a)
    if total != expected_rows(pitch, samples, mode):
        raise ValueError("Incomplete pitch/field-wire total exposure")
    result = dict(
        pitch_mm=pitch, field_wire_diameter_um=fw,
        n_electrons=total, n_anode_collected=collected,
        n_field_wire_collected=grounded,
        n_not_on_wires=total-collected-grounded,
        n_anode_collected_unknown_time=time_unknown,
        anode_collection_fraction=collected/total,
        field_wire_collection_fraction=grounded/total,
        drift_time_ns=mean_rms(times, 24100+int(10*pitch)+int(fw)),
        electron_incidence_normal_deg=mean_rms(
            incidence, 34100+int(10*pitch)+int(fw)),
        analysis_scope=("one central u column ONLY: smoke, not "
                        "representative of full unit cell" if mode == "smoke"
                        else "uniformly populated Garfield box grid "
                        "across one complete wire pitch"),
        statistical_errors=(
            "Bootstrap over independently seeded Garfield electron trials, "
            "conditional on reaching an anode for time/angle. No Geant4 "
            "muon-event uncertainty or model systematics."),
    )
    return result


def config(pitch, fw, mode, samples, seed, binary, output):
    nu = round(pitch/DU)
    if not math.isclose(nu*DU, pitch, abs_tol=1.e-9):
        raise ValueError("Pitch must tile DU exactly")
    directory = output / tag(pitch, fw)
    directory.mkdir(parents=True, exist_ok=True)
    fingerprint = {
        "study": "PI_GARFIELD_PITCH_FIELDWIRE_V1",
        "pitch_mm": pitch, "field_wire_diam_um": fw,
        "field_wire_voltage_v": 0,
        "field_wire_offsets_mm": "+/- p/2 between anodes",
        "gap_minus_mm": 3., "gap_plus_mm": 3.,
        "anode_voltage_v": 1800., "anode_diameter_um": 30.,
        "gas": "Ar/CO2 70:30, 293.15 K, 760 Torr",
        "magnetic_field_t": 0., "n_u": nu, "n_v": NV,
        "du_mm": DU, "dv_mm": DV,
        "samples_per_box": samples, "seed": seed,
        "mode": mode, "cpp_sha256": digest(CPP),
        "binary_sha256": digest(binary),
        "runner_sha256": digest(Path(__file__)),
        "collection_interpretation":
            "Geometric endpoint-near-electrode Garfield classification; "
            "field-wire collection is NOT anode collection.",
    }
    manifest = directory / "frozen_config.json"
    columns = directory / "columns"
    columns.mkdir(exist_ok=True)
    if manifest.is_file():
        old = json.loads(manifest.read_text(encoding="utf-8"))
        # Report-only script edits may change the runner hash, while
        # CPP/binary/grid/HV/gas/seeds define frozen Garfield physics.
        if not same_physics_manifest(old, fingerprint):
            if any(columns.glob("iu_*.csv")):
                # Protect completed microscopic Garfield physics data.
                raise ValueError(
                    f"Incompatible frozen config with completed columns: {directory}. "
                    "Keep these results; choose a NEW --out namespace.")
            # An ELF loader failure (exit 127) may have left an empty run.
            # Keep its original manifest rather than asking users to rm it.
            archive = directory / (
                "frozen_config.pre_retry_"
                + hashlib.sha256(manifest.read_bytes()).hexdigest()[:12]
                + ".json")
            if archive.exists():
                raise ValueError(f"Previous incompatible manifest already archived: {archive}")
            manifest.rename(archive)
            print(f"Archived failed launch settings (NO completed columns): {archive}",
                  flush=True)
    if not manifest.is_file():
        manifest.write_text(json.dumps(fingerprint, indent=2)+"\n",
                            encoding="utf-8")
    iu_values = ([nu//2] if mode == "smoke" else list(range(nu)))
    observations = []
    unknown_column_times = 0
    for iu in iu_values:
        dest = columns / f"iu_{iu:03d}.csv"
        if dest.exists():
            validate_column(dest, iu, pitch, samples)
            elapsed_saved = load_column_timing(dest)
            if elapsed_saved is None:
                unknown_column_times += 1
            else:
                observations.append({"iu": iu, "wall_seconds": elapsed_saved,
                                     "n_electrons": NV*samples,
                                     "reused_checkpoint": True})
            print(f"SKIP validated {tag(pitch,fw)} / u-column {iu}")
            continue
        tmp = columns / f"iu_{iu:03d}.csv.tmp"
        log = columns / f"iu_{iu:03d}.log"
        tmp.unlink(missing_ok=True)
        cmd = [
            str(binary), "--gap-minus-mm", "3",
            "--gap-plus-mm", "3", "--wire-pitch-mm", str(pitch),
            "--wire-diam-um", "30",
            "--field-wire-diam-um", str(fw),
            "--hv", "1800", "--b-tesla", "0",
            "--half-wires", "4", "--u-min-mm", str(-pitch/2),
            "--u-max-mm", str(pitch/2), "--v-min-mm", "-3",
            "--v-max-mm", "3", "--du-mm", str(DU),
            "--dv-mm", str(DV), "--electrons-per-box", str(samples),
            "--base-seed", str(seed), "--iu-min", str(iu),
            "--iu-max", str(iu), "--output", str(tmp),
        ]
        print(f"RUN {tag(pitch,fw)} u-column {iu}/{nu-1}: "
              f"{NV*samples} Garfield drift trials", flush=True)
        t0 = time.perf_counter()
        with log.open("w", encoding="utf-8") as stream:
            done = subprocess.run(cmd, cwd=REPO, stdout=stream,
                                  stderr=subprocess.STDOUT, check=False)
        elapsed = time.perf_counter() - t0
        if done.returncode:
            tail = "\n".join(log.read_text(
                encoding="utf-8", errors="replace").splitlines()[-45:])
            raise RuntimeError(
                f"Garfield failed with exit {done.returncode}: {log}\n"
                f"Last 45 log lines:\n{tail}\n"
                "Exit 127 commonly indicates a missing runtime shared library; "
                "check the ELF ldd preflight printed by the shell runner.")
        validate_column(tmp, iu, pitch, samples)
        os.replace(tmp, dest)
        save_column_timing(dest, elapsed, NV*samples)
        observations.append({"iu": iu, "wall_seconds": elapsed,
                             "n_electrons": NV*samples,
                             "reused_checkpoint": False})
        print(f"TIMING {tag(pitch,fw)} / u-column {iu}: "
              f"{elapsed:.2f} s ({elapsed/(NV*samples):.4f} s/electron)",
              flush=True)
    combined = directory / "garfield_electrons.csv"
    temporary = directory / "garfield_electrons.csv.tmp"
    with temporary.open("w", newline="", encoding="utf-8") as o:
        header = None
        for iu in iu_values:
            with (columns / f"iu_{iu:03d}.csv").open(
                    newline="", encoding="utf-8") as f:
                row_header = f.readline()
                if header is None:
                    header = row_header
                    o.write(header)
                elif row_header != header:
                    raise ValueError("Different Garfield CSV schemas")
                for line in f:
                    o.write(line)
    os.replace(temporary, combined)
    result = summarise(combined, pitch, fw, mode, samples)
    result["analysis_runner_sha256"] = digest(Path(__file__))
    result["column_wall_times"] = observations
    result["columns_missing_measured_wall_time"] = unknown_column_times
    result["garfield_wall_seconds_per_electron"] = mean_rms(
        [r["wall_seconds"] / r["n_electrons"] for r in observations],
        250020 + int(10*pitch) + int(fw)) if observations else None
    result["timing_includes_process_startup"] = True
    result["raw_data"] = str(combined)
    (directory / "summary.json").write_text(json.dumps(
        result, indent=2)+"\n", encoding="utf-8")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=("smoke", "pilot"), default="smoke")
    parser.add_argument("--samples", type=int, default=None)
    parser.add_argument("--field-wire-diam-um", type=float, default=75.,
                        help="Study default 75 um; may override explicitly")
    parser.add_argument("--seed", type=int, default=260020)
    parser.add_argument("--only-pitch", type=float, choices=PITCHES,
                        default=None,
                        help="Run exactly one pitch with both field-wire settings. "
                             "Default retains the original six geometries.")
    parser.add_argument("--binary", type=Path, default=DEFAULT_BINARY)
    parser.add_argument("--out", type=Path, default=ROOT)
    args = parser.parse_args()
    if args.samples is None:
        args.samples = 2 if args.mode == "smoke" else 20
    if args.samples < 1 or args.seed < 0:
        parser.error("samples must be positive; seed nonnegative")
    if not (0 < args.field_wire_diam_um < 1000):
        parser.error("field wire diameter must be positive and <1000 um")
    if not args.binary.is_file():
        parser.error("Garfield executable absent; build Stage-C CMake target")
    output = args.out / args.mode
    output.mkdir(parents=True, exist_ok=True)
    results = []
    pitches=(args.only_pitch,) if args.only_pitch is not None else PITCHES
    for pitch in pitches:
        for fw in (0., args.field_wire_diam_um):
            results.append(config(pitch, fw, args.mode, args.samples,
                                  args.seed, args.binary.resolve(), output))
    summary = output / ("two_pitch_configurations.json"
                        if args.only_pitch is not None
                        else "all_six_configurations.json")
    summary.write_text(json.dumps({
        "scope": "NEW Garfield-only one-cell B=0 field-wire pilot",
        "mode": args.mode,
        "field_diameter_provisional_um": args.field_wire_diam_um,
        "results": results,
        "not": "Not Geant4 2M run, not conditional on muon entrance angle, "
               "not finite chamber edge model or verified hardware diameter",
    }, indent=2)+"\n", encoding="utf-8")
    print(f"PASS: all {len(results)} independent configurations validated:", summary)
    for row in results:
        t = row["drift_time_ns"]
        print(f"p={row['pitch_mm']:.1f}mm, field={row['field_wire_diameter_um']:.1f}um "
              f"n={row['n_electrons']}, collected={row['n_anode_collected']}, "
              f"field_collected={row['n_field_wire_collected']}; "
              f"time mean={t['mean']} +/- {t['mean_error']} ns; "
              f"RMS={t['rms']} +/- {t['rms_error']} ns")


if __name__ == "__main__":
    main()
