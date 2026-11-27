# Shared helpers of run_k2a_validation.sh and k2a_validation_k2a.sh (sourced, not executed).
# Needs: V (muon output dir, steps.tsv), TF (timing.tsv), HS, INL, RC, TA.
vstep(){ printf '%s\t%s\t%s\t%s\n' "$1" "$2" "$3" "${4:-}" >> "$V/steps.tsv"; echo "  [$2] $1: $3"; }
tm(){ # tm LABEL CMD...: CMD unchanged; appends label, wall_s, peak_rss_kib (GNU time only), exit to $TF
  local label="$1"; shift
  local rc=0 t0 t1 rss=""
  if [[ $(type -t "$1") != function && -x /usr/bin/time ]]; then
    local tf; tf=$(mktemp)
    /usr/bin/time -f '%e %M' -o "$tf" "$@" || rc=$?
    read -r t1 rss < <(tail -1 "$tf"); rm -f "$tf"
  else
    t0=$(date +%s.%N); "$@" || rc=$?; t1=$(date +%s.%N)
    t1=$(python3 -c "print(f'{$t1 - $t0:.3f}')")
  fi
  printf '%s\t%s\t%s\t%s\n' "$label" "$t1" "$rss" "$rc" >> "$TF"
  return $rc
}
stripcmp(){ # NAME REF TEST TOL
  local rc=0
  python3 "$INL/k2a_strip_compare.py" "$2" "$3" "$V/$1.json" --tol "$4" --name "$1" > "$V/$1.log" 2>&1 || rc=$?
  case $rc in 0) vstep "$1" PASS "$(tail -1 "$V/$1.log")" "$V/$1.json";; 3) vstep "$1" INCONCLUSIVE "$(tail -1 "$V/$1.log")" "$V/$1.json";;
              *) vstep "$1" FAIL "$(tail -1 "$V/$1.log")" "$V/$1.json";; esac
}
physcmp(){ # NAME REF CAND
  local rc=0
  python3 "$RC/k2a_recenter_compare.py" "$2" "$3" --out "$V/cmp_$1.json" --kind k2a_vs_s0 --half-strips "$HS" --label "$1" > "$V/cmp_$1.txt" 2>&1 || rc=$?
  local st; st=$(python3 -c "import json,sys; print(json.load(open(sys.argv[1]))['status'])" "$V/cmp_$1.json" 2>/dev/null || echo FAIL)
  [[ $rc -le 3 ]] || st=FAIL
  vstep "cmp_$1" "$st" "$(head -c 300 "$V/cmp_$1.txt" | tr '\t\n' '  ')" "$V/cmp_$1.json"
}
shaped(){ # NAME VARIANTS [--no-prompt]  (layout V/shaped_layout/muon_M/<variant>)
  local n="$1" vars="$2"; shift 2
  local rc=0
  python3 "$TA/t_shaped_compare.py" "$V/shaped_layout" --variants "$vars" --out "$V/shaped_$n" --bootstrap 0 "$@" > "$V/shaped_$n.log" 2>&1 || rc=$?
  if [[ $rc -eq 0 ]]; then vstep "shaped_$n" OK "variants $vars ${*:-ion+electron}" "$V/shaped_$n/per_muon_Tp100.jsonl"
  else vstep "shaped_$n" FAIL "t_shaped_compare exit $rc: $(tail -1 "$V/shaped_$n.log" | head -c 250)" "$V/shaped_$n.log"; fi
}
evstep(){ # NAME RC LOG -> OK / BLOCKED / FAIL for a k2a_strip_eval run
  case $2 in
    0) vstep "$1" OK "$(tail -1 "$3" | head -c 300)" "$V/$1/result.json";;
    4) vstep "$1" BLOCKED "$(tail -1 "$3" | head -c 300)" "$V/$1/result.json";;
    *) vstep "$1" FAIL "k2a_strip_eval exit $2: $(tail -1 "$3" | head -c 250)" "$3";;
  esac
}
