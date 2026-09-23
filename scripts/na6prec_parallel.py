#!/usr/bin/env python3

import configparser
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time
import re


def take_option_value(args, index):
    arg = args[index]
    if "=" in arg and arg.startswith("--"):
        return arg.split("=", 1)[1], index + 1
    if index + 1 >= len(args):
        raise SystemExit(f"missing value for {arg}")
    return args[index + 1], index + 2


def parse_bool(value):
    return str(value).strip().lower() not in ("0", "false", "no", "off")


def chunk_sizes(nevents, n_parallel, first=0):
    n_workers = min(nevents, n_parallel)
    base = nevents // n_workers
    rem = nevents % n_workers
    chunks = []
    offset = first
    for worker in range(n_workers):
        n = base + (1 if worker < rem else 0)
        chunks.append((worker, offset, n))
        offset += n
    return chunks

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

def normalize_output_dir(path, default):
    if path in (None, "", "none"):
        return default.resolve()
    return Path(os.path.expandvars(os.path.expanduser(path))).resolve()


def split_config(values):
    output, kept = None, []
    for value in values:
        for token in value.split(";"):
            token = token.strip()
            if not token:
                continue
            if token.startswith("keyval.output_dir="):
                output = token.split("=", 1)[1]
            else:
                kept.append(token)
    return output, kept


def build_config(kept_tokens, output_dir):
    tokens = list(kept_tokens)
    tokens.append(f"keyval.output_dir={output_dir}")
    return ";".join(tokens)


def input_event_count(input_dir):
    candidates = (("MCKine.root", "mckine"),
                  ("HitsVerTel.root", "hitsVerTel"),
                  ("DigitsVerTel.root", "digitsVerTel"),
                  ("HitsMuonSpec.root", "hitsMuonSpec"),
                  ("ClustersVerTel.root", "clustersVerTel"),
                  ("ClustersMuonSpec.root", "clustersMuonSpec"))
    for filename, tree in candidates:
        path = input_dir / filename
        if not path.is_file():
            continue
        code = (f'TFile f("{path}"); auto* t = f.Get<TTree>("{tree}"); '
                f'if (t) std::cout << "NA6PREC_EVENT_COUNT=" << t->GetEntries() << std::endl;')
        result = subprocess.run(["root", "-l", "-b", "-q", "-e", code],
                                text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        match = re.search(r"NA6PREC_EVENT_COUNT=(\d+)", result.stdout)
        if result.returncode == 0 and match:
            return int(match.group(1)), path, tree
    raise SystemExit("could not determine the event count; pass -n/--nevents explicitly")


def link_inputs(source, destination, excluded):
    """Symlink shared input files into a worker directory, skipping worker outputs."""
    destination.mkdir(parents=True, exist_ok=True)
    for item in source.iterdir():
        if item.name in excluded or not item.is_file():
            continue
        target = destination / item.name
        if not target.exists() and not target.is_symlink():
            target.symlink_to(item.resolve())


def run_workers(commands):
    active, failed = [], []
    for worker, command, cwd, log_path in commands:
        cwd.mkdir(parents=True, exist_ok=True)
        print("running", " ".join(map(str, command)))
        log = open(log_path, "w")
        active.append((worker, log_path, log, subprocess.Popen(command, cwd=cwd, stdout=log, stderr=subprocess.STDOUT)))
    while active:
        remaining = []
        for worker, log_path, log, process in active:
            status = process.poll()
            if status is None:
                remaining.append((worker, log_path, log, process))
            else:
                log.close()
                if status:
                    failed.append((worker, status, log_path))
        active = remaining
        if active:
            time.sleep(1)
    if failed:
        raise SystemExit("\n".join(f"worker {worker} failed with code {status}, log: {log}" for worker, status, log in failed))


def merge_outputs(worker_dirs, output, names):
    hadd = shutil.which("hadd")
    if not hadd:
        raise SystemExit("hadd was not found in PATH; source the ROOT environment first")
    output.mkdir(parents=True, exist_ok=True)
    for name in names:
        inputs = [directory / name for directory in worker_dirs if (directory / name).is_file()]
        if not inputs:
            continue
        target = output / name
        temporary = output / f".{name}.merging"
        subprocess.run([hadd, "-f", str(temporary), *map(str, inputs)], check=True)
        temporary.replace(target)


def parse_passthrough(args):
    nevents = None
    first, last = 0, None
    config_values, base = [], []
    load_ini = None
    flags = {"--doHitsToRecPoints": True, "--doDigitsToRecPoints": False,
             "--doTrackletVertex": True, "--doVTTracking": True,
             "--doMSTracking": True, "--doMatching": True}
    stage_aliases = {"-hitcl": "--doHitsToRecPoints",
                     "-cl": "--doDigitsToRecPoints",
                     "-vert": "--doTrackletVertex",
                     "-vt": "--doVTTracking",
                     "-ms": "--doMSTracking",
                     "-mt": "--doMatching"}
    stage_options = set(flags) | set(stage_aliases) | {"--firstevent", "--lastevent", "-f", "-l"}
    i = 0
    while i < len(args):
        arg = args[i]
        key = arg.split("=", 1)[0] if arg.startswith("--") else arg
        if key in ("--nevents", "-n"):
            value, i = take_option_value(args, i)
            nevents = int(value)
        elif arg.startswith("-n") and len(arg) > 2:
            nevents = int(arg[2:])
            i += 1
        elif key in stage_options:
            value, i = take_option_value(args, i)
            if key in ("--firstevent", "-f"):
                first = int(value)
            elif key in ("--lastevent", "-l"):
                last = int(value)
            else:
                flags[stage_aliases.get(key, key)] = parse_bool(value)
        elif key == "--configKeyValues":
            value, i = take_option_value(args, i)
            config_values.append(value)
        elif key in ("--load-ini", "--load-recoparam"):
            value, i = take_option_value(args, i)
            value = str(Path(value).resolve())
            if key == "--load-ini":
                load_ini = value
            base.extend([key, value])
        elif key in ("--geometry", "-g"):
            value, i = take_option_value(args, i)
            base.extend([key, str(Path(value).resolve())])
        else:
            base.append(arg)
            i += 1
    if nevents is not None and nevents <= 0:
        raise SystemExit("--nevents must be positive")
    if first < 0:
        first = 0
    if nevents is not None and (last is None or last < 0):
        last = first + nevents - 1
    if last is not None and last < first:
        raise SystemExit("--lastevent must not precede --firstevent")
    return nevents, first, last, config_values, load_ini, flags, base


def main():
    workers, keep, input_dir = 1, False, Path.cwd()
    passthrough, i = [], 0
    while i < len(sys.argv[1:]):
        arg = sys.argv[1:][i]
        if arg in ("-h", "--help"):
            print("usage: na6prec_parallel [--workers N] [--keep-worker-output] [na6prec arguments]")
            print()
            print("Runs na6prec in event chunks, merges stage outputs, then runs matching.")
            print("The event count is read from the input ROOT files unless -n/--nevents is provided.")
            raise SystemExit(0)
        if arg == "--":
            passthrough.extend(sys.argv[1:][i + 1:])
            break
        if arg in ("--workers", "--n_parallel", "--input-dir"):
            if i + 1 >= len(sys.argv[1:]):
                raise SystemExit(f"missing value for {arg}")
            value = sys.argv[1:][i + 1]
            if arg == "--input-dir": input_dir = Path(value).resolve()
            else: workers = int(value)
            i += 2
        elif arg == "--keep-worker-output":
            keep = True; i += 1
        elif arg.startswith("--workers=") or arg.startswith("--n_parallel="):
            workers = int(arg.split("=", 1)[1]); i += 1
        else:
            passthrough.append(arg); i += 1
    if workers < 1 or not input_dir.is_dir():
        raise SystemExit("--workers must be positive and --input-dir must exist")
    nevents, first, last, configs, load_ini, flags, base = parse_passthrough(passthrough)
    if nevents is None:
        total, source, tree = input_event_count(input_dir)
        last = total - 1 if last is None or last < 0 else min(last, total - 1)
        nevents = last - first + 1
        print(f"using {nevents} events from {source.name}:{tree}")
    else:
        nevents = min(nevents, last - first + 1)
    if nevents <= 0:
        raise SystemExit("the selected event range is empty")
    output_from_args, kept = split_config(configs)
    final = normalize_output_dir(output_from_args or read_ini_output_dir(load_ini), Path.cwd())
    final.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(prefix=f".{final.name or 'na6prec'}_workers_", dir=final.parent))
    executable = shutil.which("na6prec")
    if not executable:
        raise SystemExit("na6prec was not found in PATH")
    event_chunks = chunk_sizes(nevents, workers, first)
    completed = False
    try:
        if flags["--doHitsToRecPoints"] or flags["--doDigitsToRecPoints"]:
            commands, dirs = [], []
            for worker, offset, count in event_chunks:
                directory = temporary / "clusters" / f"worker{worker:03d}"; dirs.append(directory)
                link_inputs(input_dir, directory, {"ClustersVerTel.root", "ClustersMuonSpec.root", "TracksVerTel.root", "TracksMuonSpec.root", "TracksMatching.root", "VerticesVerTel.root"})
                command = [executable, *base, "--firstevent", str(offset), "--lastevent", str(offset + count - 1),
                           "--doHitsToRecPoints", str(flags["--doHitsToRecPoints"]).lower(), "--doDigitsToRecPoints", str(flags["--doDigitsToRecPoints"]).lower(),
                           "--doTrackletVertex", "false", "--doVTTracking", "false", "--doMSTracking", "false", "--doMatching", "false",
                           "--configKeyValues", build_config(kept, directory)]
                commands.append((worker, command, directory, directory / "na6prec.log"))
            print(f"clusterizing {nevents} events with {len(commands)} workers")
            run_workers(commands); merge_outputs(dirs, final, ("ClustersVerTel.root", "ClustersMuonSpec.root"))
        link_inputs(input_dir, final, {"TracksVerTel.root", "TracksMuonSpec.root", "TracksMatching.root", "VerticesVerTel.root"})

        if flags["--doTrackletVertex"] or flags["--doVTTracking"] or flags["--doMSTracking"]:
            commands, dirs = [], []
            for worker, offset, count in event_chunks:
                directory = temporary / "tracking" / f"worker{worker:03d}"; dirs.append(directory)
                link_inputs(final, directory, {"TracksVerTel.root", "TracksMuonSpec.root", "TracksMatching.root", "VerticesVerTel.root"})
                link_inputs(input_dir, directory, {"TracksVerTel.root", "TracksMuonSpec.root", "TracksMatching.root", "VerticesVerTel.root"})
                command = [executable, *base, "--firstevent", str(offset), "--lastevent", str(offset + count - 1),
                           "--doHitsToRecPoints", "false", "--doDigitsToRecPoints", "false",
                           "--doTrackletVertex", str(flags["--doTrackletVertex"]).lower(), "--doVTTracking", str(flags["--doVTTracking"]).lower(),
                           "--doMSTracking", str(flags["--doMSTracking"]).lower(), "--doMatching", "false",
                           "--configKeyValues", build_config(kept, directory)]
                commands.append((worker, command, directory, directory / "na6prec.log"))
            print(f"tracking {nevents} events with {len(commands)} workers")
            run_workers(commands); merge_outputs(dirs, final, ("TracksVerTel.root", "TracksMuonSpec.root", "VerticesVerTel.root"))
            
        if flags["--doMatching"]:
            command = [executable, *base, "--firstevent", str(first), "--lastevent", str(first + nevents - 1),
                       "--doHitsToRecPoints", "false", "--doDigitsToRecPoints", "false", "--doTrackletVertex", "false",
                       "--doVTTracking", "false", "--doMSTracking", "false", "--doMatching", "true",
                       "--configKeyValues", build_config(kept, final)]
            print("matching merged tracks")
            subprocess.run(command, cwd=final, check=True)
        completed = True
    finally:
        if keep or not completed: print(f"kept worker output: {temporary}")
        else: shutil.rmtree(temporary, ignore_errors=True)


if __name__ == "__main__":
    main()
