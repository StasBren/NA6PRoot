#!/usr/bin/env bash
# Safely rebuild current feature/mwpc-stage-c code with the optional
# layout-independent muon plane audit using an EXISTING CMake build.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd -P)"
if ! command -v cmake >/dev/null 2>&1; then
  echo "ERROR: cmake not in PATH" >&2
  exit 2
fi

CANDIDATES=()
if [[ -n "${NA6PBUILDDIR:-}" ]]; then CANDIDATES+=("$NA6PBUILDDIR"); fi
CANDIDATES+=(
  "$REPO/build"
  "$HOME/na6p/build"
  "$HOME/na6p/build/NA6PRoot"
  "$HOME/na6p/build/NA6PRoot-mwpc"
  "$HOME/na6p/src/NA6PRoot-mwpc/build"
)
MATCHES=()
for dir in "${CANDIDATES[@]}"; do
  [[ -f "$dir/CMakeCache.txt" ]] || continue
  home="$(sed -n 's/^CMAKE_HOME_DIRECTORY:INTERNAL=//p' "$dir/CMakeCache.txt" | head -1)"
  [[ -d "$home" ]] || continue
  canonical="$(cd "$home" && pwd -P)"
  [[ "$canonical" == "$REPO" ]] || continue
  realdir="$(cd "$dir" && pwd -P)"
  already=false
  for old in "${MATCHES[@]}"; do
    if [[ "$old" == "$realdir" ]]; then already=true; fi
  done
  if [[ "$already" == false ]]; then MATCHES+=("$realdir"); fi
done
if [[ "${#MATCHES[@]}" != 1 ]]; then
  echo "ERROR: expected exactly ONE existing CMake build tree for $REPO" >&2
  echo "Found ${#MATCHES[@]} matching build directories:" >&2
  for p in "${MATCHES[@]}"; do echo "  $p" >&2; done
  echo "Provide the correct existing NA6PBUILDDIR, not a new guessed tree." >&2
  exit 3
fi
BUILD="${MATCHES[0]}"
PREFIX="$(sed -n 's/^CMAKE_INSTALL_PREFIX:PATH=//p' "$BUILD/CMakeCache.txt" | head -1)"
if [[ -z "$PREFIX" ]]; then
 echo "ERROR: Cannot establish install prefix in CMake cache" >&2
 exit 4
fi
echo "Branch source : $REPO"
echo "Build tree    : $BUILD"
echo "Install prefix: $PREFIX"
echo "Reconfiguring to include the new simulator reference-audit header..."
cmake -S "$REPO" -B "$BUILD"
echo "Building and installing updated Geant4 simulator and simLib..."
cmake --build "$BUILD" --target install --parallel 4
# Success concerns the target built and installed from THIS checkout.
# A login shell may still point to the OLD NA6PRoot installation, but
# run_pi_geant4_campaign.sh now activates THIS prefix automatically.
if [[ ! -x "$PREFIX/bin/na6psim" ||
      ! -x "$PREFIX/bin/na6psim_parallel" ||
      ! -f "$PREFIX/lib/libsimLib.so" ]]; then
 echo "ERROR: Stage-C installed files incomplete at $PREFIX" >&2
 exit 5
fi
echo "PASS: rebuilt Stage-C na6psim and libsimLib installed at $PREFIX"
echo "The PI campaign wrapper will select this installation automatically."
