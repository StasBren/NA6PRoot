#!/usr/bin/env bash
#PBS -N mwpc_2100_raw
#PBS -l select=1:ncpus=1:mem=32gb
#PBS -l walltime=12:00:00
#PBS -j oe
# One *complete* true-Geant4 muon per OpenPBS array element. No avalanche size cap.
# Resource failures are FAILURES, not valid gain measurements.
set -Eeuo pipefail
umask 027

: "${PBS_ARRAY_INDEX:?This file must run as a PBS array task}"
: "${MWPC_2100_REPO:?Pass MWPC_2100_REPO with qsub -v}"
: "${MWPC_2100_OUT:?Pass MWPC_2100_OUT with qsub -v}"
: "${MWPC_2100_SELECTION_SHA:?Pass MWPC_2100_SELECTION_SHA with qsub -v}"
: "${MWPC_2100_BINARY_SHA:?Pass MWPC_2100_BINARY_SHA with qsub -v}"
: "${MWPC_2100_GIT_SHA:?Pass MWPC_2100_GIT_SHA with qsub -v}"

idx="${PBS_ARRAY_INDEX}"
[[ "$idx" =~ ^[0-9]+$ ]] && (( idx >= 0 && idx < 1000 )) || {
  echo "Invalid array index $idx (need 0..999)" >&2; exit 2;
}
REPO="${MWPC_2100_REPO}"
OUT="${MWPC_2100_OUT}"
SELECTION="$REPO/test/mwpc-stage-c/p3_wiretruth_1000_selection.csv"
EXE="$REPO/test/mwpc-avalanche/build-hv2100-uncapped/mwpc_p3_noFW_penning057_2100_uncapped_v1"
printf -v tag 'muon_%04d' "$idx"
DIR="$OUT/$tag"

# PBS can start with a minimal environment; load only our known environment setup.
if [[ -f "$HOME/na6p/env.sh" ]]; then source "$HOME/na6p/env.sh"; fi
if [[ -f "$HOME/na6p/install/garfieldpp/share/Garfield/setupGarfield.sh" ]]; then
  source "$HOME/na6p/install/garfieldpp/share/Garfield/setupGarfield.sh"
fi

[[ -s "$SELECTION" && -x "$EXE" ]] || {
  echo "Missing source selection or Garfield executable" >&2; exit 3;
}
[[ "$(sha256sum "$SELECTION" | awk '{print $1}')" == "$MWPC_2100_SELECTION_SHA" ]] || {
  echo "SELECTION SHA changed since submission" >&2;exit 3;
}
[[ "$(sha256sum "$EXE" | awk '{print $1}')" == "$MWPC_2100_BINARY_SHA" ]] || {
  echo "BINARY SHA changed since submission" >&2;exit 3;
}
[[ "$(git -C "$REPO" rev-parse HEAD)" == "$MWPC_2100_GIT_SHA" ]] || {
  echo "GIT revision changed since submission" >&2;exit 3;
}
[[ -s "$OUT/config.lock" ]] || {
  echo "Missing immutable run lock" >&2;exit 3;
}

if [[ -e "$DIR/complete.txt" ]];then
  (cd "$DIR" && sha256sum -c SHA256SUMS >/dev/null) || {
    echo "Completed output checksum mismatch: $DIR" >&2; exit 3;
  }
  echo "SKIP already complete and SHA verified: $tag";exit 0
fi
# No deletion or rewriting of partial results: quarantine explicitly before resubmitting.
if ! mkdir "$DIR" 2>/dev/null;then
  echo "Refusing to overwrite incomplete or concurrent result directory: $DIR" >&2
  exit 4
fi

# Index 0 => CSV line 2, index 999 => line 1001. Original crossings/seeds unchanged.
awk -v target="$((idx+2))" 'NR==1 || NR==target' "$SELECTION" > "$DIR/input.csv"
[[ "$(wc -l < "$DIR/input.csv")" -eq 2 ]] || {
  echo "Failed to extract exactly one Geant4 crossing" >&2;exit 4;
}
{
  printf 'index=%s\nPBS_JOBID=%s\nstart_utc=%s\n' "$idx" "${PBS_JOBID:-unknown}" "$(date -u -Is)"
  printf 'git_sha=%s\nselection_sha256=%s\nbinary_sha256=%s\n' \
      "$MWPC_2100_GIT_SHA" "$MWPC_2100_SELECTION_SHA" "$MWPC_2100_BINARY_SHA"
  printf 'physics=P3_noFW,ArCO2_70_30,Penning_r0.57,B0,HV2100\n'
  printf 'Garfield_avalanche_size_limit=DISABLED\n'
} > "$DIR/task_metadata.txt"

# One uncapped full microscopic avalanche per true TrackHeed primary electron.
# The trailing zero is NOT a cutoff; this independent executable insists on zero.
# timeout exits NONZERO and does not seal the result. Resource walltime / RAM
# are for machine protection, never for truncating avalanche gains.
if ! timeout --signal=TERM --kill-after=30s 11h30m \
     "$EXE" --wiretruth-hv-p3-noFW "$DIR/input.csv" "$DIR/microscopic" 2100 0 \
     > "$DIR/physics.log" 2>&1; then
  echo "FAILED (OOM / timeout / Garfield): $DIR" >&2
  tail -n 30 "$DIR/physics.log" >&2
  exit 5
fi

for kind in crossing seed wire_muon wire_seed ion_birth_sample;do
  [[ -s "$DIR/microscopic_${kind}.csv" ]] || {
    echo "Missing microscopic_${kind}.csv" >&2;exit 5;
  }
done
[[ "$(wc -l < "$DIR/microscopic_crossing.csv")" -eq 2 ]] || {
  echo "Expected exactly 1 completed muon crossing" >&2;exit 5;
}
grep -q '^GARFIELD_AVALANCHE_SIZE_LIMIT_DISABLED=1$' "$DIR/physics.log" || {
  echo "No uncapped Garfield startup proof" >&2;exit 5;
}
# Verify per-seed bookkeeping; this does not exclude legitimate zero-gain seeds.
awk -F, '
  NR==1 {
    for(i=1;i<=NF;i++){
      if($i=="capped") cap=i;
      if($i=="successful") ok=i;
    }
    if(cap==0 || ok==0) exit 2;
  }
  NR>1 { if($(cap)+0!=0 || $(ok)+0!=1) bad++; }
  END { if(NR<2 || bad>0) exit 1; }
' "$DIR/microscopic_seed.csv" || {
  echo "Capped or failed primary avalanche: $DIR" >&2;exit 5;
}

(cd "$DIR" && sha256sum input.csv physics.log task_metadata.txt microscopic_*.csv > SHA256SUMS)
printf 'DONE HV=2100 index=%s muons=1 no_avalanche_size_limit=1\n' "$idx" > "$DIR/complete.txt"
echo "PASS 2100V uncapped: $tag"
