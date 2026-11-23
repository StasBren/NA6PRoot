#!/usr/bin/env bash
# Read-only diagnosis of event 0878. No Garfield reruns or archive writes.
set -Eeuo pipefail
if [[ -f "$HOME/na6p/env.sh" ]]; then source "$HOME/na6p/env.sh"; fi
if [[ -f "$HOME/na6p/install/garfieldpp/share/Garfield/setupGarfield.sh" ]]; then
  source "$HOME/na6p/install/garfieldpp/share/Garfield/setupGarfield.sh"
fi
archive="$HOME/na6p/cathode-overnight-2100/test_runs/mwpc_stage_c/overnight_2100_weighted_k256_v1"
if [[ -n "$(printenv MWPC_ARCHIVE_2100_K256 || true)" ]]; then archive="$MWPC_ARCHIVE_2100_K256"; fi
d="$archive/muon_0878"
[[ -d "$d" ]] || { echo "Missing: $d" >&2; exit 2; }
echo "=== Read-only 2100V muon_0878 diagnostic ==="
echo "Archive: $d"
echo
echo "=== File sizes ==="
for name in ion_births.root microscopic_wire_muon.csv microscopic_wire_seed.csv microscopic_seed.csv SHA256SUMS complete.txt root_audit.log; do
  if [[ -e "$d/$name" ]]; then stat -c '%n  bytes=%s' "$d/$name"; else echo "MISSING $name"; fi
done
echo
echo "=== Manifest format ==="
if [[ -f "$d/SHA256SUMS" ]]; then
  LC_ALL=C sed -n '1,6l' "$d/SHA256SUMS"
  if LC_ALL=C grep -Eq '^[0-9a-fA-F]{64} [ *]' "$d/SHA256SUMS"; then
    (cd "$d" && sha256sum -c SHA256SUMS) || echo "WARNING: checksum check failed"
  else
    echo "SHA256SUMS is not standard sha256sum -c format; do not infer corruption."
  fi
fi
echo
echo "=== CSV schemas and independent sums ==="
python3 - "$d" <<'PY'
import csv,sys
from pathlib import Path
d=Path(sys.argv[1])
for name in ("microscopic_wire_muon.csv","microscopic_wire_seed.csv","microscopic_seed.csv"):
    print("FILE",name)
    p=d/name
    if not p.is_file():
        print("  MISSING"); continue
    try:
        with p.open(encoding="utf-8-sig",newline="") as f:
            rd=csv.reader(f)
            print("  header:",repr(next(rd,None)))
            print("  data 1:",repr(next(rd,None)))
            print("  data 2:",repr(next(rd,None)))
        with p.open(encoding="utf-8-sig",newline="") as f:
            rd=csv.DictReader(f)
            need={"primary_index","collected_e","ion_births"}
            names=set(rd.fieldnames or [])
            if not need.issubset(names):
                print("  Missing expected columns:",sorted(need-names))
                continue
            n=0; ions=0; collected=0; levels=set()
            for row in rd:
                if not row or not any(row.values()): continue
                n+=1; ions+=int(row["ion_births"]); collected+=int(row["collected_e"])
                levels.add(int(row["primary_index"]))
            print(f"  SUM rows={n} collected={collected} ion_births={ions}")
            print("  primary_index sample:",sorted(levels)[:15])
    except Exception as exc:
        print("  PARSE ERROR:",type(exc).__name__,str(exc))
PY
echo
echo "=== ROOT primary sums ==="
if command -v root >/dev/null 2>&1; then
  tmp="$(mktemp -d /tmp/na60_0878_diag.XXXXXXXX)"
  trap 'rm -rf "$tmp"' EXIT
  export MWPC_AUDIT_ROOT_FILE="$d/ion_births.root"
  cat > "$tmp/na60_0878_census.C" <<'CPP'
#include <TFile.h>
#include <TTree.h>
#include <TLeaf.h>
#include <cstdlib>
#include <iostream>
#include <memory>
void na60_0878_census() {
 const char* path=std::getenv("MWPC_AUDIT_ROOT_FILE");
 if (!path) { std::cerr<<"No ROOT env path\n"; return; }
 std::unique_ptr<TFile> f(TFile::Open(path,"READ"));
 if (!f || f->IsZombie()) { std::cerr<<"Cannot read ROOT\n"; return; }
 auto* m=dynamic_cast<TTree*>(f->Get("muons"));
 auto* p=dynamic_cast<TTree*>(f->Get("primaries"));
 if (!m || !p) { std::cerr<<"Missing ROOT muons/primaries\n"; return; }
 std::cout<<"ROOT muons="<<m->GetEntries()<<" primaries="<<p->GetEntries()<<"\n";
 Long64_t ni=0, ne=0, totalI=0, totalE=0;
 for(const char* bname:{"ion_births_exact","collected_electrons"}) {
   auto* b=p->GetBranch(bname);
   if (!b) { std::cerr<<"Missing "<<bname<<"\n"; return; }
   auto* leaf=b->GetLeaf(bname);
   std::cout<<"branch "<<bname<<" type="<<(leaf?leaf->GetTypeName():"unknown")<<"\n";
 }
 if(p->SetBranchAddress("ion_births_exact",&ni)<0 ||
    p->SetBranchAddress("collected_electrons",&ne)<0) {
    std::cerr<<"Incompatible census branch types\n";return;
 }
 for(Long64_t j=0;j<p->GetEntries();++j){
   p->GetEntry(j);
   totalI+=ni; totalE+=ne;
 }
 std::cout<<"ROOT SUM collected="<<totalE<<" ion_births="<<totalI<<"\n";
}
CPP
  root -l -b -q "$tmp/na60_0878_census.C" || echo "WARNING: ROOT diagnostic failed"
else echo "ROOT executable not in environment"; fi
echo
echo "=== Existing root_audit.log tail ==="
if [[ -f "$d/root_audit.log" ]]; then tail -n 15 "$d/root_audit.log"; fi
echo
echo "=== complete.txt ==="
if [[ -f "$d/complete.txt" ]]; then sed -n '1,10p' "$d/complete.txt"; fi
echo "DONE (read-only). Do not delete or skip event 878 yet."
