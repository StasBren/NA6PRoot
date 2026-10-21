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
# A Garfield-only caller needs NA6PRoot, not the ROOT merging executable.
# Opt in for native Geant4 ROOT/hadd production and run full loader checks.
if [[ "${NA6P_STAGEC_REQUIRE_ROOT_HADD:-0}" == "1" ]]; then
# ROOT is installed separately from NA6PRoot. Resolve ROOT binaries AND
# ROOT shared libraries from one matching installation before any Geant4 run.
_stagec_root_candidates=()
if [[ -n "${NA6P_ROOT_PREFIX:-}" ]]; then
  _stagec_root_candidates+=("$NA6P_ROOT_PREFIX")
fi
if [[ -n "${ROOTSYS:-}" ]]; then
  _stagec_root_candidates+=("$ROOTSYS")
fi
_stagec_root_candidates+=("$HOME/na6p/install/root")
if command -v root-config >/dev/null 2>&1; then
  _stagec_root_cfg="$(root-config --prefix 2>/dev/null || true)"
  [[ -z "$_stagec_root_cfg" ]] || _stagec_root_candidates+=("$_stagec_root_cfg")
fi
_stagec_root_prefix=""
_stagec_root_libdir=""
for _candidate in "${_stagec_root_candidates[@]}"; do
  [[ -n "$_candidate" ]] || continue
  for _lib in "$_candidate/lib" "$_candidate/lib64"; do
    if [[ -x "$_candidate/bin/root" && -x "$_candidate/bin/hadd" &&
          -f "$_lib/libCore.so" && -f "$_lib/libRIO.so" &&
          -f "$_lib/libTree.so" ]]; then
      _stagec_root_prefix="$(cd "$_candidate" && pwd -P)"
      _stagec_root_libdir="$(cd "$_lib" && pwd -P)"
      break 2
    fi
  done
done
if [[ -z "$_stagec_root_prefix" ]]; then
  echo "ERROR: cannot resolve BOTH root and hadd from the SAME ROOT install." >&2
  echo "Checked ROOT installation directories:" >&2
  printf '  %s\n' "${_stagec_root_candidates[@]}" >&2
  echo "Use NA6P_ROOT_PREFIX=/absolute/path/to/root/install if needed." >&2
  echo "Do not start Geant4 until bin/root, bin/hadd and libraries agree." >&2
  return 6
fi
export ROOTSYS="$_stagec_root_prefix"
export PATH="$_stagec_root_prefix/bin:$PATH"
export LD_LIBRARY_PATH="$_stagec_root_libdir:${LD_LIBRARY_PATH:-}"
export ROOT_INCLUDE_PATH="$_stagec_root_prefix/include:${ROOT_INCLUDE_PATH:-}"
hash -r
_active_root_real="$(readlink -f "$(command -v root)")"
_active_hadd_real="$(readlink -f "$(command -v hadd)")"
_expected_root_real="$(readlink -f "$_stagec_root_prefix/bin/root")"
_expected_hadd_real="$(readlink -f "$_stagec_root_prefix/bin/hadd")"
if [[ "$_active_root_real" != "$_expected_root_real" ]] ||
   [[ "$_active_hadd_real" != "$_expected_hadd_real" ]]; then
  echo "ERROR: root/hadd resolve to different ROOT installations" >&2
  return 7
fi
echo "ROOT runtime selected:"
echo "  root:              $(command -v root)"
echo "  hadd:              $(command -v hadd)"
echo "  libraries:         $_stagec_root_libdir"
echo "  ROOTSYS:           $_stagec_root_prefix"
for _program in "$_stagec_root_prefix/bin/hadd" "$_stagec_prefix/bin/na6psim"; do
  _unresolved="$(ldd "$_program" 2>/dev/null | grep 'not found' || true)"
  if [[ -n "$_unresolved" ]]; then
    echo "ERROR: unresolved dynamic libraries for $_program:" >&2
    echo "$_unresolved" >&2
    return 8
  fi
done
unset _stagec_root_candidates _stagec_root_cfg _candidate _lib
unset _stagec_root_prefix _stagec_root_libdir _program _unresolved
unset _active_root_real _active_hadd_real _expected_root_real _expected_hadd_real

fi

echo "Stage-C installed simulator selected:"
echo "  na6psim:          $_active_sim"
echo "  na6psim_parallel: $_active_parallel"
echo "  libraries:        $_stagec_prefix/lib"
unset _build_candidates _matching_builds _dir _source _found _old
unset _stagec_build _stagec_prefix _active_sim _active_parallel
