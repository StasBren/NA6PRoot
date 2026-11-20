#!/usr/bin/env bash
# Read ONLY seal-verified independent muons. No microscopic regeneration.
# 10ns single-ion response bank: SAME 256 first-muon training sites.
set -Eeuo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"
BIN="$REPO/test/mwpc-avalanche/build-timing"
RUN="${MWPC_OVERNIGHT_OUT:-$HOME/na6p/cathode-overnight-2100/test_runs/mwpc_stage_c/overnight_2100_weighted_k256_v1}"
TRAIN="${MWPC_TIMING_TRAIN_ROOT:-$HOME/na6p/cathode-timing-2100/test_runs/mwpc_stage_c/plane_timing_dt10_v1/templates256_dt10_training.root}"
MOBILITY="${MWPC_MOBILITY:-$HOME/na6p/src/garfieldpp/Data/IonMobility_Ar+_Ar.txt}"
OUT="$REPO/test_runs/mwpc_stage_c/${MWPC_UNSEEN_OUTPUT_TAG:-unseen_muons_10ns_v1}"
fail(){ echo "UNSEEN TEST STOP: $*" >&2; exit 2; }
[[ -x "$BIN/mwpc_plane_response_from_ions_root" &&
   -x "$BIN/mwpc_plane_muon_kernel_compare" ]] || fail "build required ion-response targets first"
[[ -s "$TRAIN" && -s "$MOBILITY" ]] || fail "training bank or ion mobility missing"
[[ -d "$RUN" ]] || fail "missing microscopic muon data directory"
mkdir -p "$OUT"
RESULTS="$OUT/independent_muon_10ns_metrics.csv"
[[ ! -e "$RESULTS" ]] || fail "refuse overwrite $RESULTS; change MWPC_UNSEEN_OUTPUT_TAG"
echo "index,source,station,plane,window_ns,direct_charge_fC,fast_charge_fC,error_fC,error_pct_of_direct,full_waveform_relative_L2" > "$RESULTS"
for spec in '0435:jpsi40-bore:4' '0511:phi40-bore:3' '0101:omega40-bore:4'; do
  IFS=: read -r idx source station <<<"$spec"
  dir="$RUN/muon_$idx"
  [[ -s "$dir/complete.txt" && -s "$dir/SHA256SUMS" &&
     -s "$dir/ion_births.root" ]] || fail "unsealed muon_$idx"
  (cd "$dir" && sha256sum -c SHA256SUMS >/dev/null) ||
    fail "checksum mismatch muon_$idx"
  grep -qx 'PASS CATHODE ROOT PHYSICS AUDIT' "$dir/root_audit.log" ||
    fail "no valid ROOT audit muon_$idx"
  grep -qx 'GARFIELD_AVALANCHE_SIZE_LIMIT_DISABLED=1' "$dir/physics.log" ||
    fail "uncapped Garfield flag missing"
  d="$OUT/muon_$idx"
  mkdir -p "$d"
  direct="$d/direct_K256_dt10"
  fast="$d/fast_K256_dt10"
  [[ ! -e "${direct}_summary.csv" &&
     ! -e "${fast}_summary.csv" ]] || fail "do not overwrite output muon_$idx"
  echo "BEGIN independent muon_$idx"
  "$BIN/mwpc_plane_response_from_ions_root" \
    --ion-file "$dir/ion_births.root" \
    --out-prefix "$direct" --mobility "$MOBILITY" --mode direct \
    --dt-ns 10 --window-us 100 --max-step-mm 0.05 \
    >"$d/direct.log" 2>&1 || {
      tail -n 15 "$d/direct.log";fail "direct RKF muon_$idx"; }
  "$BIN/mwpc_plane_muon_kernel_compare" \
    --ion-file "$dir/ion_births.root" \
    --training-root "$TRAIN" \
    --reference-waveforms "${direct}_waveforms.csv" \
    --reference-summary "${direct}_summary.csv" \
    --out-prefix "$fast" --neighbors 4 --dt-ns 10 --window-us 100 \
    >"$d/fast.log" 2>&1 || {
      tail -n 15 "$d/fast.log";fail "kernel muon_$idx"; }
  awk -F, -v idx="$idx" -v source="$source" -v station="$station" '
    BEGIN{dt=10; m[1]=100;m[2]=200;m[3]=400;m[4]=1000;m[5]=100000}
    NR==1 {
      if ($1!="time_center_ns" || $2!="plus_kernel_fC_per_ns" ||
          $4!="plus_direct_fC_per_ns")exit 3;
      next
    }
    {
      n++;
      for(s=1;s<=2;s++){
        k=(s==1?$2:$3)+0; r=(s==1?$4:$5)+0;
        delta=k-r; E[s]+=delta*delta; D[s]+=r*r;
        for(j=1;j<=5;j++)if(n*dt<=m[j]){
          Qk[s,j]+=k*dt;Qr[s,j]+=r*dt
        }
      }
    }
    END{
      if(n!=10000)exit 5;
      for(s=1;s<=2;s++){
        plane=(s==1?"plus":"minus");
        l2=sqrt(E[s]/(D[s]>0?D[s]:1e-60));
        for(j=1;j<=5;j++){
          delta=Qk[s,j]-Qr[s,j];
          pct=(Qr[s,j]==0?"NA":sprintf("%.9g",100*delta/Qr[s,j]));
          printf("%s,%s,%s,%s,%d,%.15g,%.15g,%.15g,%s,%.12g\n",
            idx,source,station,plane,m[j],Qr[s,j],Qk[s,j],delta,pct,l2)
        }
      }
    }
  ' "${fast}_waveforms.csv" >> "$RESULTS" ||
    fail "cannot summarize waveforms muon_$idx"
  echo "DONE independent muon_$idx"
done
echo "PASS three independent 10ns whole-plane ion-response comparisons"
echo "SUMMARY $RESULTS"
echo "LIMIT: no stereo strips, no VMM, no independent K256 reservoir variance test"
