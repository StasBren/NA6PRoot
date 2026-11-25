#!/usr/bin/env python3
"""Fallback for '/usr/bin/time -v -o OUT cmd...' where GNU time is not installed: runs cmd, writes the child's peak RSS
in GNU time's wording ('Maximum resident set size (kbytes): N') and the wall time to OUT, returns the child's exit code.
  rss_run.py OUT cmd [args...]"""
import resource
import subprocess
import sys
import time

out, cmd = sys.argv[1], sys.argv[2:]
t0 = time.monotonic()
rc = subprocess.call(cmd)
wall = time.monotonic() - t0
ru = resource.getrusage(resource.RUSAGE_CHILDREN)
with open(out, "w") as f:
    f.write(f"\tCommand being timed: \"{' '.join(cmd)}\"\n\tMaximum resident set size (kbytes): {ru.ru_maxrss}\n"
            f"\tElapsed (wall clock) time (seconds): {wall:.3f}\n\tExit status: {rc}\n\t(measured by rss_run.py: getrusage RUSAGE_CHILDREN)\n")
sys.exit(rc if rc >= 0 else 128 - rc)
