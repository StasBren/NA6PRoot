#!/usr/bin/env python3

import configparser
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time


def take_option_value(args, index):
    arg = args[index]
    if "=" in arg and arg.startswith("--"):
        return arg.split("=", 1)[1], index + 1
    if index + 1 >= len(args):
        raise SystemExit(f"missing value for {arg}")
    return args[index + 1], index + 2


def parse_bool(value):
    return str(value).strip().lower() not in ("0", "false", "no", "off")


def parse_passthrough(args):
    nevents = 1
    config_strings = []
    load_ini = None
    do_digitization = True
    digitize_only = False
    rnd_seed = None
    base = []
    i = 0
    while i < len(args):
        arg = args[i]
        key = arg.split("=", 1)[0] if arg.startswith("--") else arg
        if key in ("-n", "--nevents"):
            val, i = take_option_value(args, i)
            nevents = int(val)
        elif arg.startswith("-n") and len(arg) > 2:
            nevents = int(arg[2:])
            i += 1
        elif key == "--configKeyValues":
            val, i = take_option_value(args, i)
            config_strings.append(val)
        elif key == "--load-ini":
            val, ni = take_option_value(args, i)
            load_ini = val
            base.extend([arg, val] if "=" not in arg else [arg])
            i = ni
        elif key in ("--event-offset", "--skip-events"):
            _, i = take_option_value(args, i)
        elif key in ("-r", "--rnd-seed"):
            val, i = take_option_value(args, i)
            rnd_seed = int(val)
        elif arg.startswith("-r") and len(arg) > 2:
            rnd_seed = int(arg[2:])
            i += 1
        elif key in ("--doDigitization", "-dig"):
            val, i = take_option_value(args, i)
            do_digitization = parse_bool(val)
        elif key == "--digitize-only":
            digitize_only = True
            i += 1
        else:
            base.append(arg)
            i += 1
    return nevents, config_strings, load_ini, do_digitization, digitize_only, rnd_seed, base

def split_config(config_strings):
    output_dir = None
    kept = []
    for conf in config_strings:
        for token in conf.split(";"):
            token = token.strip()
            if not token:
                continue
            if token.startswith("keyval.output_dir="):
                output_dir = token.split("=", 1)[1]
            else:
                kept.append(token)
    return output_dir, kept


def read_ini_output_dir(path):
    if not path:
        return None
    cp = configparser.ConfigParser()
    cp.optionxform = str
    try:
        cp.read(path)
    except configparser.Error:
        return None
    if cp.has_section("keyval") and cp.has_option("keyval", "output_dir"):
        return cp.get("keyval", "output_dir")
    return None


def normalize_output_dir(path):
    if path in (None, "", "none"):
        return Path(".").resolve()
    return Path(os.path.expandvars(os.path.expanduser(path))).resolve()


def chunk_sizes(nevents, n_parallel):
    n_workers = min(nevents, n_parallel)
    base = nevents // n_workers
    rem = nevents % n_workers
    chunks = []
    offset = 0
    for worker in range(n_workers):
        n = base + (1 if worker < rem else 0)
        chunks.append((worker, offset, n))
        offset += n
    return chunks


def build_config(kept_tokens, output_dir):
    tokens = list(kept_tokens)
    tokens.append(f"keyval.output_dir={output_dir}")
    return ";".join(tokens)


def run_workers(worker_cmds):
    procs = []
    for worker, cmd, log_path in worker_cmds:
        log_path.parent.mkdir(parents=True, exist_ok=True)
        print("running", " ".join(map(str, cmd)))
        log = open(log_path, "w")
        procs.append((worker, cmd, log_path, log, subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT)))

    failed = []
    while procs:
        alive = []
        for worker, cmd, log_path, log, proc in procs:
            ret = proc.poll()
            if ret is None:
                alive.append((worker, cmd, log_path, log, proc))
                continue
            log.close()
            if ret != 0:
                failed.append((worker, ret, log_path))
        procs = alive
        if procs:
            time.sleep(1)

    if failed:
        msg = "\n".join(f"worker {w} failed with code {ret}, log: {log}" for w, ret, log in failed)
        raise SystemExit(msg)


def make_worker_commands(na6psim, base_args, kept_config, tmp_root, chunks, rnd_seed):
    commands = []
    for worker, offset, nevents in chunks:
        worker_dir = tmp_root / f"worker{worker:03d}"
        cmd = [
            na6psim,
            *base_args,
            "-n", str(nevents),
            "--event-offset", str(offset),
            "--skip-events", str(offset),
            "--doDigitization", "false",
        ]
        if rnd_seed is not None and rnd_seed >= 0:
            cmd.extend(["-r", str(rnd_seed + worker)])
        elif rnd_seed is not None:
            cmd.extend(["-r", str(rnd_seed)])
        cmd.extend(["--configKeyValues", build_config(kept_config, worker_dir)])
        commands.append((worker, cmd, worker_dir / "na6psim.log"))
    return commands


def run_final_digitization(na6psim, base_args, kept_config, final_out, rnd_seed):
    cmd = make_digitize_only_command(na6psim, base_args, kept_config, final_out, rnd_seed)
    print("digitizing", " ".join(map(str, cmd)))
    subprocess.run(cmd, check=True)


def print_plan(worker_cmds, tmp_root, final_out, do_digitization):
    print(f"will run {len(worker_cmds)} processes:")
    for _, cmd, log_path in worker_cmds:
        print("  " + " ".join(map(str, cmd)))
        print(f"    log: {log_path}")
    action = "merge results and digitize" if do_digitization else "merge results"
    print(f"with worker output to {tmp_root}, then {action} to {final_out}")


def print_digitize_only_plan(na6psim, base_args, kept_config, final_out, rnd_seed):
    cmd = make_digitize_only_command(na6psim, base_args, kept_config, final_out, rnd_seed)
    print("will run 1 process:")
    print("  " + " ".join(map(str, cmd)))
    print(f"with digitization output to {final_out}")


def make_digitize_only_command(na6psim, base_args, kept_config, final_out, rnd_seed):
    cmd = [na6psim, *base_args, "--digitize-only"]
    if rnd_seed is not None:
        cmd.extend(["-r", str(rnd_seed)])
    cmd.extend(["--configKeyValues", build_config(kept_config, final_out)])
    return cmd

def merge_outputs(tmp_root, final_out, chunks):
    hadd = shutil.which("hadd")
    if not hadd:
        raise SystemExit("hadd was not found in PATH; source ROOT environment before running na6psim_parallel")

    final_out.mkdir(parents=True, exist_ok=True)
    worker_dirs = [tmp_root / f"worker{worker:03d}" for worker, _, _ in chunks]
    root_names = sorted({p.name for d in worker_dirs for p in d.glob("*.root") if p.name != "geometry.root"})

    for name in root_names:
        inputs = [d / name for d in worker_dirs if (d / name).exists()]
        if not inputs:
            continue
        out = final_out / name
        if len(inputs) == 1:
            shutil.copy2(inputs[0], out)
        else:
            subprocess.run([hadd, "-f", str(out), *map(str, inputs)], check=True)

    for name in ("geometry.root", "na6pLayout.ini"):
        src = worker_dirs[0] / name
        if src.exists():
            shutil.copy2(src, final_out / name)


def split_wrapper_args(argv):
    workers = 1
    keep_worker_output = False
    passthrough = []
    i = 0
    while i < len(argv):
        arg = argv[i]
        if arg in ("-h", "--help"):
            print("usage: na6psim_parallel [--workers N] [--keep-worker-output] [na6psim arguments]")
            print()
            print("Runs na6psim in event chunks and merges worker ROOT outputs.")
            print("All arguments except --workers and --keep-worker-output are passed to na6psim.")
            raise SystemExit(0)
        if arg == "--":
            passthrough.extend(argv[i + 1:])
            break
        if arg in ("--workers", "--n_parallel"):
            if i + 1 >= len(argv):
                raise SystemExit(f"missing value for {arg}")
            workers = int(argv[i + 1])
            i += 2
            continue
        if arg.startswith("--workers=") or arg.startswith("--n_parallel="):
            workers = int(arg.split("=", 1)[1])
            i += 1
            continue
        if arg == "--keep-worker-output":
            keep_worker_output = True
            i += 1
            continue
        passthrough.append(arg)
        i += 1
    if workers < 1:
        raise SystemExit("--workers must be positive")
    return workers, keep_worker_output, passthrough

def main():
    workers, keep_worker_output, passthrough = split_wrapper_args(sys.argv[1:])
    nevents, config_strings, load_ini, do_digitization, digitize_only, rnd_seed, base_args = parse_passthrough(passthrough)

    out_from_args, kept_config = split_config(config_strings)
    final_out = normalize_output_dir(out_from_args or read_ini_output_dir(load_ini))
    tmp_parent = final_out.parent if final_out.name else Path(".").resolve()
    tmp_parent.mkdir(parents=True, exist_ok=True)
    tmp_root = Path(tempfile.mkdtemp(prefix=f".{final_out.name or 'na6psim'}_workers_", dir=tmp_parent))

    na6psim = shutil.which("na6psim")
    if not na6psim:
        raise SystemExit("na6psim was not found in PATH")

    if digitize_only:
        print_digitize_only_plan(na6psim, base_args, kept_config, final_out, rnd_seed)
        shutil.rmtree(tmp_root, ignore_errors=True)
        run_final_digitization(na6psim, base_args, kept_config, final_out, rnd_seed)
        return

    if nevents < 0:
        raise SystemExit("parallel mode needs an explicit non-negative -n/--nevents value")
    if nevents == 0:
        raise SystemExit("nothing to run: -n/--nevents is 0")

    chunks = chunk_sizes(nevents, workers)
    worker_cmds = make_worker_commands(na6psim, base_args, kept_config, tmp_root, chunks, rnd_seed)
    print_plan(worker_cmds, tmp_root, final_out, do_digitization)
    try:
        run_workers(worker_cmds)
        merge_outputs(tmp_root, final_out, chunks)
        if do_digitization:
            run_final_digitization(na6psim, base_args, kept_config, final_out, rnd_seed)
    finally:
        if keep_worker_output:
            print(f"kept worker output: {tmp_root}")
        else:
            shutil.rmtree(tmp_root, ignore_errors=True)


if __name__ == "__main__":
    main()
