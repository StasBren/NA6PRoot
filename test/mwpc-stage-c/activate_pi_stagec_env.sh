#!/usr/bin/env bash
# Source this file from a Stage-C runner. Select the installed simulator
# built from THIS checkout, without changing the caller's login shell.
#
# Both na6psim_parallel (PATH) and its worker na6psim (PATH), as well as
# the ELF loader (LD_LIBRARY_PATH), MUST resolve to the SAME install.
#
# Caller must define REPO as absolute realpath of this checkout.
if [[ -z "${REPO:-}" || ! -d "$REPO" ]]; then
  echo "ERROR: source activation with REPO set to the checkout" >&2
  return 2
fi
REPO="$(cd "$REPO" && pwd -P)"
_build_candidates=()
if [[ -n "${NA6PBUILDDIR:-}" ]]; then
  _build_candidates+=("$NA6PBUILDDIR")
fi
_build_candidates+=(
  "$REPO/build"
  "$HOME/na6p/build"
  "$HOME/na6p/build/NA6PRoot"
  "$HOME/na6p/build/NA6PRoot-mwpc"
  "$HOME/na6p/src/NA6PRoot-mwpc/build"
)
_matching_builds=()
for _dir in "${_build_candidates[@]}"; do
  [[ -f "$_dir/CMakeCache.txt" ]] || continue
  _source="$(sed -n 's/^CMAKE_HOME_DIRECTORY:INTERNAL=//p' "$_dir/CMakeCache.txt" | head -1)"
  [[ -d "$_source" ]] || continue
  _source="$(cd "$_source" && pwd -P)"
  [[ "$_source" == "$REPO" ]] || continue
  _dir="$(cd "$_dir" && pwd -P)"
  _found=0
  for _old in "${_matching_builds[@]}"; do
    [[ "$_old" == "$_dir" ]] && _found=1
  done
  if [[ $_found -eq 0 ]]; then
    _matching_builds+=("$_dir")
  fi
done
if [[ ${#_matching_builds[@]} -ne 1 ]]; then
  echo "ERROR: expected ONE CMake build belonging to $REPO, found ${#_matching_builds[@]}" >&2
  echo "Set NA6PBUILDDIR to the existing matching build tree." >&2
  return 3
fi
_stagec_build="${_matching_builds[0]}"
_stagec_prefix="$(sed -n 's/^CMAKE_INSTALL_PREFIX:PATH=//p' "$_stagec_build/CMakeCache.txt" | head -1)"
if [[ -z "$_stagec_prefix" || ! -x "$_stagec_prefix/bin/na6psim" ||
      ! -x "$_stagec_prefix/bin/na6psim_parallel" ||
      ! -f "$_stagec_prefix/lib/libsimLib.so" ]]; then
  echo "ERROR: Stage-C build exists but its installed simulator, wrapper or simLib is missing." >&2
  echo "Build: $_stagec_build; prefix: $_stagec_prefix" >&2
  echo "Run bash test/mwpc-stage-c/rebuild_pi_geant4_stagec.sh first." >&2
  return 4
fi
_stagec_prefix="$(cd "$_stagec_prefix" && pwd -P)"
export PATH="$_stagec_prefix/bin:$PATH"
export LD_LIBRARY_PATH="$_stagec_prefix/lib:${LD_LIBRARY_PATH:-}"
export ROOT_INCLUDE_PATH="$_stagec_prefix/include:${ROOT_INCLUDE_PATH:-}"
export NA6P_STAGEC_INSTALL_PREFIX="$_stagec_prefix"
# Match the canonical installed init.sh.in project-root variable; field,
# geometry and any external ROOT macros must see the same installation.
export NA6PROOT_ROOT="$_stagec_prefix"
hash -r
_active_sim="$(command -v na6psim)"
_active_parallel="$(command -v na6psim_parallel)"
if [[ "$(readlink -f "$_active_sim")" != "$(readlink -f "$_stagec_prefix/bin/na6psim")" ||
      "$(readlink -f "$_active_parallel")" != "$(readlink -f "$_stagec_prefix/bin/na6psim_parallel")" ]]; then
  echo "ERROR: cannot activate this source checkout's native tools" >&2
  return 5
fi
echo "Stage-C installed simulator selected:"
echo "  na6psim:          $_active_sim"
echo "  na6psim_parallel: $_active_parallel"
echo "  libraries:        $_stagec_prefix/lib"
unset _build_candidates _matching_builds _dir _source _found _old
unset _stagec_build _stagec_prefix _active_sim _active_parallel
