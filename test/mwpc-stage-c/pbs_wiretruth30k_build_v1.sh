#!/usr/bin/env bash
#PBS -N wt30build
#PBS -q N
#PBS -m n
#PBS -l select=1:ncpus=1:mem=4gb
#PBS -l walltime=02:00:00
#PBS -l io=0.1

# Build existing dedicated Garfield CMake target on compute node,
# stage IDENTICAL old 30k original Geant4 inputs, freeze binary/runtime.
set -Eeuo pipefail
REPO=${WT_REPO:-"$HOME/NA6PRoot-wiretruth30k"}
ROOT=${WT_ROOT:-/storage/cgrp/stanis/mwpc_penning057_wiretruth30k_v1}
PROGRAM="$REPO/test/mwpc-stage-c/pi_wiretruth30k_production_v1.py"
BUILD="$REPO/test/mwpc-avalanche/build-wiretruth30k-cvmfs"
TARGET=mwpc_p3_noFW_penning057_wiretruth_v1

export GARFIELD_INSTALL=/cvmfs/sft-nightlies.cern.ch/lcg/nightlies/dev3/Mon/Garfield++/HEAD/x86_64-el9-gcc13-opt
source "$GARFIELD_INSTALL/Garfield++-env.sh"
echo "BUILD_JOB=${PBS_JOBID:-manual} NODE=$(hostname)"
echo "GIT_COMMIT=$(git -C "$REPO" rev-parse HEAD)"
echo "GARFIELD_INSTALL=$GARFIELD_INSTALL"
echo "ROOT_VERSION=$(root-config --version)"
echo "GCC_VERSION=$(g++ -dumpfullversion)"
cmake -S "$REPO/test/mwpc-avalanche" -B "$BUILD" -DCMAKE_BUILD_TYPE=Release
cmake --build "$BUILD" --target "$TARGET" -j1
python3 -m py_compile "$PROGRAM"
python3 "$PROGRAM" prepare --root "$ROOT" \
  --old-root /storage/cgrp/stanis/mwpc_penning057
python3 "$PROGRAM" lock-binary --root "$ROOT"
python3 "$PROGRAM" runtime-check --root "$ROOT"
echo "PASS BUILD LOCKED, no expensive Garfield muons run"
