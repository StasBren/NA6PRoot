#!/usr/bin/env bash
# Controlled cathode strip-geometry study: ONE alternative readout geometry on the five development muons of an
# existing S-ladder run. Weighting side only: archived trajectories, Green-grid paths, avalanche, electronics unchanged.
# READ-ONLY on archives, the frozen LUT/replay/exporter and the baseline run; every product goes to new directories.
#
#   bash run_geometry_study.sh selftest WORK [ALT_LUT GEOM]    SYNTHETIC tests (restrip, C++ core, + toy pipeline if LUT given)
#   bash run_geometry_study.sh overlay-identity                rebuild the BASELINE LUT through the overlay; compare with the
#                                                              frozen MWPC_SPECTRAL_LUT array by array (no-op proof, ~4 min)
#   bash run_geometry_study.sh check GEOM                      geometry id, transport-coupled constants, frozen module sha
#   bash run_geometry_study.sh lut GEOM                        LUT for GEOM (frozen builder, default settings) + independent
#                                                              numerical checks (geom_lut_check) of it AND of the frozen LUT
#   bash run_geometry_study.sh prepare BASE TAG GEOM LUT IDX.. per muon: strip frame from the archive (p3_birth_moments,
#                                                              default geometry) -> restriped T1 export -> S0 (exact replay)
#                                                              and electron step with LUT
#   bash run_geometry_study.sh green BASE TAG GEOM LUT         Green table from BASE/green_paths_TAG (REUSED, not re-drifted)
#   bash run_geometry_study.sh closure BASE TAG GEOM           G_TAG vs S0 per muon (s_synthesize G) + shaped comparison
#   bash run_geometry_study.sh kernel BASE TAG GEOM            K1/K2/K1A/K2A kernels + export-path moments + apply
#   bash run_geometry_study.sh births BASE TAG GEOM IDX..      birth-record moments with GEOM vs export path; K2A/K2 waveforms
#   bash run_geometry_study.sh engine BASE TAG GEOM            batch reader + response engine (compiled/vector) with GEOM
#   bash run_geometry_study.sh report BASE TAG GEOM            accuracy ladder vs S0 and vs G + geometry effect (base vs alt)
#   bash run_geometry_study.sh hashes BASE TAG GEOM            geometry/weighting/transport/kernel identities -> manifest
# BASE = baseline S-ladder run (muon_*/T1_export, S0, green_paths_TAG, green_TAG, kernel_TAG, *kern_TAG).
# GEOM = geometry json (geometries/*.json). Output: BASE/geometry_<name>_<id12>/ (a new run directory).
# Environment: as run_s_ladder.sh (MWPC_NUMW_DIR, archive and selection paths), plus MWPC_BIRTHS_BIN (p3_birth_moments
# built from THIS commit: it needs --stereo-tan), NA60_KERNEL_EVAL_LIB for `engine`, MWPC_GEOM_OUT for `lut` (default
# $REPO/test_runs/mwpc_stage_c/geometry_study_2100). No sealed muon, no PBS, no avalanche, no Geant4.
set -Eeo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SA="$(cd "$HERE/../spatial-ablation" && pwd)"; TA="$(cd "$HERE/../timing-ablation" && pwd)"
REPO="$(cd "$HERE/../../../.." && pwd)"
set -u
umask 077
ARCHIVE="${MWPC_ARCHIVE_2100_K256:-$HOME/na6p/cathode-overnight-2100/test_runs/mwpc_stage_c/overnight_2100_weighted_k256_v1}"
SELECTION="${MWPC_SIGNAL40_OUTPUT:-$HOME/NA6PRoot-signal40-2100/test_runs/mwpc_stage_c/phase_space_2100_signal40_v1}"
DEV="$SELECTION/development20.csv"; SEALED="$SELECTION/sealed_validation20.csv"
GOUT="${MWPC_GEOM_OUT:-$REPO/test_runs/mwpc_stage_c/geometry_study_2100}"
HS="${MWPC_HALF_STRIPS:-5}"; DSUB="${MWPC_S_DSUB:-10}"
BASEGEOM="$HERE/geometries/baseline.json"
fail(){ echo "GEOMETRY STUDY FAILED: $*" >&2; exit 2; }
need(){ [[ -n "${!1:-}" ]] || fail "set $1"; [[ -e "${!1}" ]] || fail "$1=${!1} does not exist"; }
idx(){ awk -F, 'NR==1{for(i=1;i<=NF;i++)if($i=="index")c=i; if(!c)exit 3; next} NF{printf "%04d\n",$c}' "$1"; }
dev_list(){ [[ -s "$DEV" && -s "$SEALED" ]] || fail "selection lists missing"; comm -23 <(idx "$DEV" | sort) <(idx "$SEALED" | sort); }
verified(){ local d="$ARCHIVE/muon_$1"; [[ -s "$d/ion_births.root" && -s "$d/SHA256SUMS" ]] && (cd "$d" && sha256sum -c SHA256SUMS >/dev/null 2>&1); }
gval(){ python3 -c "import sys; sys.path.insert(0,'$HERE'); import weighting_geometry as wg; g=wg.Geometry.load('$1'); print($2)"; }
gdir(){ echo "$1/geometry_$(gval "$2" 'g.name')_$(gval "$2" 'g.geometry_id[:12]')"; }
GR(){ python3 "$HERE/geom_run.py" "$GEOM" "$MWPC_NUMW_DIR" "$@"; }
nonbaseline(){ [[ $(gval "$GEOM" 'g.is_baseline') == False ]] || fail "GEOM is the baseline geometry: nothing to study"; }

case "${1:-}" in
  selftest)
    W="${2:?WORK dir (new)}"; [[ ! -e "$W" ]] || fail "$W exists"; mkdir -p "$W"
    python3 "$HERE/tests/test_restrip.py" "$W/restrip" "$HERE/geometries/stereo015.json"
    bash "$HERE/tests/test_core_geometry.sh" "$W/core" "$HERE/geometries/stereo015.json"
    if [[ $# -ge 4 ]]; then need MWPC_NUMW_DIR; need MWPC_SPECTRAL_LUT
      bash "$HERE/tests/test_pipeline_geometry.sh" "$W/pipeline" "$MWPC_NUMW_DIR" "$MWPC_SPECTRAL_LUT" "$3" "$4"; fi;;
  overlay-identity)
    need MWPC_NUMW_DIR; need MWPC_SPECTRAL_LUT; mkdir -p "$GOUT/luts"; L="$GOUT/luts/overlay_identity_$(date +%Y%m%d_%H%M%S).npz"
    GEOM="$BASEGEOM"; GR lut-build "$L" | grep -v retained
    python3 - "$MWPC_SPECTRAL_LUT" "$L" <<'PY'
import json, sys
import numpy as np
a, b = (np.load(p) for p in sys.argv[1:])
ok = all(a[k].shape == b[k].shape and a[k].tobytes() == b[k].tobytes() for k in ("near", "bulk"))
ma, mb = json.loads(a["meta"].item()), json.loads(b["meta"].item())
print({k: float(np.max(abs(a[k] - b[k]))) for k in ("near", "bulk")}, "meta equal:", ma == mb)
print(("PASS overlay with the baseline geometry rebuilds the frozen LUT bit for bit" if ok and ma == mb else
       "NOT BIT-IDENTICAL: the frozen LUT was built with another numpy/scipy or settings; see max |diff| above"))
PY
    ;;
  check)
    GEOM="${2:?GEOM}"; need MWPC_NUMW_DIR
    GR describe;;
  lut)
    GEOM="${2:?GEOM}"; need MWPC_NUMW_DIR; need MWPC_SPECTRAL_LUT; nonbaseline
    mkdir -p "$GOUT/luts"; N="$(gval "$GEOM" 'g.name')_$(gval "$GEOM" 'g.geometry_id[:12]')"; L="$GOUT/luts/lut_$N.npz"
    [[ ! -e "$L" ]] || fail "$L exists"
    t0=$(date +%s); GR lut-build "$L" | tee "$GOUT/luts/lut_$N.build.log"; echo "lut_build_s=$(( $(date +%s) - t0 ))" >> "$GOUT/luts/lut_$N.build.log"
    python3 "$HERE/geom_lut_check.py" "$GEOM" "$MWPC_NUMW_DIR" "$L" "$GOUT/luts/lut_$N.check.json" | grep -v retained | tee "$GOUT/luts/lut_$N.check.txt"
    [[ -e "$GOUT/luts/frozen_baseline.check.json" ]] || python3 "$HERE/geom_lut_check.py" "$BASEGEOM" "$MWPC_NUMW_DIR" "$MWPC_SPECTRAL_LUT" \
        "$GOUT/luts/frozen_baseline.check.json" | grep -v retained | tee "$GOUT/luts/frozen_baseline.check.txt"
    grep -q '^PASS' "$GOUT/luts/lut_$N.check.txt" || fail "LUT check"
    echo "LUT $L  (pass it to prepare/green)";;
  prepare)
    BASE="${2:?BASE run}"; TAG="${3:?TAG}"; GEOM="${4:?GEOM}"; LUT="${5:?LUT}"; shift 5; need MWPC_NUMW_DIR; need MWPC_BIRTHS_BIN; nonbaseline
    [[ $# -ge 1 ]] || fail "list muon indices"
    G="$(gdir "$BASE" "$GEOM")"; mkdir -p "$G"; cp -n "$GEOM" "$G/geometry.json"
    for M in "$@"; do
      M=$(printf %04d "$((10#$M))"); R0="$BASE/muon_$M"; R="$G/muon_$M"
      dev_list | grep -qx "$M" || fail "$M is not a non-sealed development muon"; verified "$M" || fail "archive muon_$M not verified"
      [[ -s "$R0/T1_export/vertices.f64" && -s "$R0/S0/replay/waveforms.csv" ]] || fail "$R0 is not a prepared baseline muon"
      [[ ! -e "$R" ]] || fail "$R exists"; mkdir -p "$R/frame"
      "$MWPC_BIRTHS_BIN" "$ARCHIVE/muon_$M/ion_births.root" "$R/frame/default_geometry.moments.csv" > "$R/frame/births.log"
      python3 "$HERE/geom_restrip_export.py" "$R0/T1_export" "$R/frame/default_geometry.moments.csv.frame.json" "$GEOM" "$R/T1_export" | tee "$R/restrip.log"
      t0=$(date +%s)
      GR replay "$LUT" "$R/T1_export" "$R/S0/replay" --half-strips "$HS" > "$R/S0.replay.log" 2>&1 || fail "S0 replay $M ($R/S0.replay.log)"
      echo "S0_replay_s=$(( $(date +%s) - t0 ))" >> "$R/timing.txt"
      GR prompt "$LUT" "$R/T1_export" "$R/S0/prompt" --replay-dir "$MWPC_NUMW_DIR" --half-strips "$HS" > "$R/S0.prompt.log" 2>&1 || fail "electron step $M"
      echo "PREPARED $R"
    done;;
  green)
    BASE="${2:?BASE run}"; TAG="${3:?TAG}"; GEOM="${4:?GEOM}"; LUT="${5:?LUT}"; need MWPC_NUMW_DIR; nonbaseline
    G="$(gdir "$BASE" "$GEOM")"; GP="$BASE/green_paths_$TAG"; GD="$G/green_$TAG"
    [[ -s "$GP/vertices.f64" ]] || fail "$GP missing"; [[ ! -e "$GD" ]] || fail "$GD exists"
    python3 "$SA/s_grid_spec.py" "$G"/muon_*/T1_export --half-strips "$HS" --strip-pitch-mm "$(gval "$GEOM" 'g.strip_pitch_mm')" | tee "$G/grid_spec_$TAG.txt"
    DM=$(sed -n 's/.*--d-max-mm \([0-9.]*\)$/\1/p' "$G/grid_spec_$TAG.txt"); [[ -n "$DM" ]] || fail "grid spec parse"
    # the reused grid must still cover every ion birth of the restriped exports (transport-side spec is unchanged)
    W0=$(sed -n 's/^--wires \([^ ]*\) --r-um.*/\1/p' "$BASE/grid_spec_$TAG.txt"); W1=$(sed -n 's/^--wires \([^ ]*\) --r-um.*/\1/p' "$G/grid_spec_$TAG.txt")
    R0=$(sed -n 's/.*--r-um \([^ ]*\)$/\1/p' "$BASE/grid_spec_$TAG.txt"); R1=$(sed -n 's/.*--r-um \([^ ]*\)$/\1/p' "$G/grid_spec_$TAG.txt")
    [[ "$W0" == "$W1" && "$R0" == "$R1" ]] || fail "wires/radii of the restriped exports differ from the baseline grid spec: transport changed?"
    t0=$(date +%s)
    GR green "$LUT" "$GP" "$GD" --replay-dir "$MWPC_NUMW_DIR" --d-max-mm "$DM" --d-sub "$DSUB" | tee "$G/green_table_$TAG.log"
    echo "green_table_${TAG}_s=$(( $(date +%s) - t0 ))" >> "$G/timing.txt"
    python3 - "$BASE/green_$TAG/axes.json" "$GD/axes.json" <<'PY'
import json, sys
a, b = (json.load(open(p)) for p in sys.argv[1:])
assert a["grid_vertices_sha256"] == b["grid_vertices_sha256"] and a["grid_ions_sha256"] == b["grid_ions_sha256"], "transport id differs"
assert a["lut_sha256"] != b["lut_sha256"], "same LUT as the baseline"
print(f"PASS transport_id unchanged {b['grid_vertices_sha256'][:16]}; weighting_id {a['lut_sha256'][:12]} -> {b['lut_sha256'][:12]}; "
      f"stereo {a['stereo_tan']} -> {b['stereo_tan']}, pitch {a['pitch_mm']} -> {b['pitch_mm']}")
PY
    ;;
  closure)
    BASE="${2:?BASE run}"; TAG="${3:?TAG}"; GEOM="${4:?GEOM}"; G="$(gdir "$BASE" "$GEOM")"; GD="$G/green_$TAG"; S=$(date +%H%M%S)
    [[ -s "$GD/axes.json" ]] || fail "run green first"
    for R in "$G"/muon_*; do [[ -e "$R/G_$TAG" ]] || python3 "$SA/s_synthesize.py" "$GD" "$R/T1_export" "$R/S0/replay/waveforms.csv" "$R" --variants G --suffix "_$TAG" | tee -a "$R/synth_$TAG.log"; done
    GA=(--stereo-tan "$(gval "$GEOM" 'g.stereo_tan')" --strip-pitch-mm "$(gval "$GEOM" 'g.strip_pitch_mm')")
    python3 "$TA/t_shaped_compare.py" "$G" --variants "S0,G_$TAG" --tp 100,200 --no-prompt --plots "${GA[@]}" --out "$G/closure_1_ion_only_${TAG}_$S"
    python3 "$TA/t_shaped_compare.py" "$G" --variants "S0,G_$TAG" --tp 100,200 --plots "${GA[@]}" --out "$G/closure_2_ion_plus_electron_step_${TAG}_$S"
    echo "return: $G/closure_*_${TAG}_$S/t_ladder_summary.txt";;
  kernel)
    BASE="${2:?BASE run}"; TAG="${3:?TAG}"; GEOM="${4:?GEOM}"; G="$(gdir "$BASE" "$GEOM")"; GD="$G/green_$TAG"; KD="$G/kernel_$TAG"
    [[ -s "$GD/axes.json" ]] || fail "run green first"
    [[ -e "$KD" ]] || python3 "$SA/s_kernel.py" build "$GD" "$KD" | tee "$G/kernel_build_$TAG.log"
    for R in "$G"/muon_*; do
      [[ -e "$R/cloud_moments.csv" ]] || python3 "$SA/s_kernel.py" moments "$R/T1_export" "$R/cloud_moments.csv"
      for m in K1 K2 K1A K2A; do
        [[ -e "$R/${m}kern_$TAG" ]] || python3 "$SA/s_kernel.py" apply "$KD" "$R/cloud_moments.csv" "$R/S0/replay/waveforms.csv" "$R" --model "$m" --suffix "_$TAG" --repeat 5 | tee -a "$R/kernel_$TAG.log"
      done
    done;;
  births)
    BASE="${2:?BASE run}"; TAG="${3:?TAG}"; GEOM="${4:?GEOM}"; shift 4; need MWPC_BIRTHS_BIN; G="$(gdir "$BASE" "$GEOM")"; KD="$G/kernel_$TAG"
    [[ -s "$KD/kernel.json" ]] || fail "run kernel first"
    GA=(--stereo-tan "$(gval "$GEOM" 'g.stereo_tan')" --strip-pitch-mm "$(gval "$GEOM" 'g.strip_pitch_mm')")
    for M in "$@"; do
      M=$(printf %04d "$((10#$M))"); R="$G/muon_$M"; dev_list | grep -qx "$M" || fail "$M is sealed or unknown"
      B="$R/births_$(date +%Y%m%d_%H%M%S)"; mkdir -p "$B"
      "$MWPC_BIRTHS_BIN" "$ARCHIVE/muon_$M/ion_births.root" "$B/birth_moments.csv" --repeat 5 "${GA[@]}" | tee "$B/births.log"
      grep -q '^PASS birth moments' "$B/births.log" || fail "p3_birth_moments $M"
      python3 "$SA/s_compare_moments.py" "$R/cloud_moments.csv" "$B/birth_moments.csv" --report "$B/compare_moments.md" > /dev/null \
        || { cat "$B/compare_moments.md"; fail "birth moments of $M differ from the restriped export path"; }
      for m in K2A K2; do
        python3 "$SA/s_kernel.py" apply "$KD" "$B/birth_moments.csv" "$R/S0/replay/waveforms.csv" "$B" --model "$m" --suffix "_${TAG}_births" --repeat 5
        python3 "$SA/s_compare_waveforms.py" "$R/${m}kern_$TAG" "$B/${m}kern_${TAG}_births" --report "$B/compare_waveforms_$m.md" > /dev/null \
          || { cat "$B/compare_waveforms_$m.md"; fail "$m waveforms of $M differ"; }
      done
      echo "PASS muon $M ($B)"
    done;;
  engine)
    BASE="${2:?BASE run}"; TAG="${3:?TAG}"; GEOM="${4:?GEOM}"; need MWPC_BIRTHS_BIN; need NA60_KERNEL_EVAL_LIB
    G="$(gdir "$BASE" "$GEOM")"; KD="$G/kernel_$TAG"; E="$G/engine_${TAG}_$(date +%Y%m%d_%H%M%S)"; mkdir -p "$E/compare"; : > "$E/events.txt"
    for R in "$G"/muon_*; do M=${R##*/muon_}; dev_list | grep -qx "$M" || fail "$M sealed"; echo "muon_$M $ARCHIVE/muon_$M/ion_births.root" >> "$E/events.txt"; done
    "$MWPC_BIRTHS_BIN" --batch "$E/events.txt" "$E/moments" --stereo-tan "$(gval "$GEOM" 'g.stereo_tan')" --strip-pitch-mm "$(gval "$GEOM" 'g.strip_pitch_mm')" | tee "$E/births_batch.log"
    grep -q '^PASS birth moments batch' "$E/births_batch.log" || fail "batch birth reader"
    EC="$SA/engine/engine_cli.py"; KID=$(python3 "$EC" info "$KD" | python3 -c "import json,sys; print(json.load(sys.stdin)['kernel_id'])"); echo "kernel_id=$KID" > "$E/kernel_id.txt"
    for impl in compiled vector; do
      python3 "$EC" run --kernel "$KD" --kernel-id "$KID" --moments-dir "$E/moments" --out "$E/replay_$impl" --impl "$impl" --format replay
      python3 "$EC" run --kernel "$KD" --kernel-id "$KID" --moments-dir "$E/moments" --out "$E/bench_$impl" --impl "$impl" --format npz --repeat 5
    done
    for line in $(cut -d' ' -f1 "$E/events.txt"); do
      python3 "$SA/s_compare_moments.py" "$G/$line/cloud_moments.csv" "$E/moments/$line.moments.csv" --report "$E/compare/${line}_moments.md" > /dev/null || fail "moments $line"
      for impl in compiled vector; do
        python3 "$SA/s_compare_waveforms.py" "$G/$line/K2Akern_$TAG" "$E/replay_$impl/$line" --report "$E/compare/${line}_wave_$impl.md" > /dev/null || fail "waveforms $line $impl"
      done
    done
    python3 "$SA/s_engine_report.py" "$E" --out "$E/engine_report.md" > /dev/null || true
    echo "PASS engine with $(gval "$GEOM" 'g.name'): $E";;
  report)
    BASE="${2:?BASE run}"; TAG="${3:?TAG}"; GEOM="${4:?GEOM}"; G="$(gdir "$BASE" "$GEOM")"; S="$G/report_${TAG}_$(date +%Y%m%d_%H%M%S)"
    GA=(--stereo-tan "$(gval "$GEOM" 'g.stereo_tan')" --strip-pitch-mm "$(gval "$GEOM" 'g.strip_pitch_mm')")
    V="S0,G_$TAG,K1kern_$TAG,K2kern_$TAG,K1Akern_$TAG,K2Akern_$TAG"; VG="G_$TAG,K1kern_$TAG,K2kern_$TAG,K1Akern_$TAG,K2Akern_$TAG"
    for sig in ion_only electron_step; do
      NP=(); [[ $sig == ion_only ]] && NP=(--no-prompt)
      python3 "$TA/t_shaped_compare.py" "$G" --variants "$V" --tp 100,200 "${NP[@]}" "${GA[@]}" --out "$S/vs_S0_$sig"
      python3 "$TA/t_shaped_compare.py" "$G" --variants "$VG" --tp 100,200 "${NP[@]}" "${GA[@]}" --out "$S/vs_G_$sig"
      python3 "$SA/s_ladder_report.py" --vs-s0 "$S/vs_S0_$sig" --vs-g "$S/vs_G_$sig" --tag "$TAG" --out "$S/accuracy_$sig.md"
      python3 "$HERE/geom_physics_compare.py" --base "$BASE" --alt "$G" --base-geom "$BASEGEOM" --alt-geom "$GEOM" \
        --variants "S0,K2Akern_$TAG" "${NP[@]}" --out "$S/geometry_effect_$sig" > /dev/null
    done
    echo "return: $S/accuracy_*.md and $S/geometry_effect_*/summary.md";;
  hashes)
    BASE="${2:?BASE run}"; TAG="${3:?TAG}"; GEOM="${4:?GEOM}"; G="$(gdir "$BASE" "$GEOM")"
    python3 "$HERE/geom_manifest.py" "$BASE" "$G" "$TAG" --out "$G/identities_$(date +%Y%m%d_%H%M%S).json";;
  *) sed -n 2,26p "$0"; exit 1;;
esac
