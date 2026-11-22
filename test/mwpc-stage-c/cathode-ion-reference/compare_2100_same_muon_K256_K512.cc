// ROOT-only cross-run identical microscopic-physics fingerprint of K256 vs K512.
// Sampling size may change. Source Geant4/Heed inputs, Garfield seeds, exact
// ion callbacks and collected-electron endpoints MUST NOT change.
#include <TFile.h>
#include <TTree.h>
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>
namespace {
void Require(bool b,const std::string& s) {
  if(!b)throw std::runtime_error("2100V paired K256/K512 fingerprint: "+s);
}
struct Muon {
  std::string id,source;
  int st=0,det=0,heed=0,k=0;
};
struct Primary {
  int index=0,cluster=0,seed=0,success=0;
  long long births=0,collected=0,endpoints=0,saved=0;
  double u=0,v=0,w=0,t=0;
};
struct Ref {
  Muon m;
  std::vector<Primary> prim;
};
Ref Load(const char* path){
  std::unique_ptr<TFile> f(TFile::Open(path,"READ"));
  Require(f&&!f->IsZombie(),std::string("missing ROOT ")+path);
  auto* m=dynamic_cast<TTree*>(f->Get("muons"));
  auto* p=dynamic_cast<TTree*>(f->Get("primaries"));
  Require(m&&p&&m->GetEntries()==1,"missing ROOT source trees");
  Ref R;std::string* id=nullptr,*ch=nullptr,*mode=nullptr;
  m->SetBranchAddress("crossing_id",&id);
  m->SetBranchAddress("source_channel",&ch);
  m->SetBranchAddress("storage_mode",&mode);
  m->SetBranchAddress("station",&R.m.st);
  m->SetBranchAddress("det_id",&R.m.det);
  m->SetBranchAddress("heed_seed",&R.m.heed);
  m->SetBranchAddress("reservoir_per_primary",&R.m.k);
  m->GetEntry(0);
  Require(id&&ch&&mode&&*mode=="weighted","wrong weighted reference mode");
  R.m.id=*id;R.m.source=*ch;
  Primary a;
  p->SetBranchAddress("primary_index",&a.index);
  p->SetBranchAddress("cluster_index",&a.cluster);
  p->SetBranchAddress("garfield_seed",&a.seed);
  p->SetBranchAddress("successful",&a.success);
  p->SetBranchAddress("ion_births_exact",&a.births);
  p->SetBranchAddress("ion_births_saved",&a.saved);
  p->SetBranchAddress("collected_electrons",&a.collected);
  p->SetBranchAddress("electron_endpoints",&a.endpoints);
  p->SetBranchAddress("u_patch_cm",&a.u);
  p->SetBranchAddress("v_absolute_cm",&a.v);
  p->SetBranchAddress("w_relative_cm",&a.w);
  p->SetBranchAddress("t_primary_ns",&a.t);
  for(Long64_t i=0;i<p->GetEntries();++i){
    p->GetEntry(i);
    Require(a.index==i&&a.success==1&&a.births>=0&&a.collected>=0,
      "invalid original successful primary");
    Require(a.saved==std::min<long long>(a.births,R.m.k),
      "sample size not consistent with original ion birth count");
    R.prim.push_back(a);
  }
  return R;
}
void Equal(double x,double y,const std::string& key){
  Require(std::isfinite(x)&&std::isfinite(y) &&
          std::abs(x-y)<=1e-10*std::max(1.,std::max(std::abs(x),std::abs(y))),
          "changed physical Heed start "+key);
}
}
int main(int argc,char**argv){
 try{
  Require(argc==3,"usage compare_2100_same_muon_K256_K512 K256.root K512.root");
  const auto a=Load(argv[1]),b=Load(argv[2]);
  Require(a.m.k==256&&b.m.k==512,
      "must compare original weighted K256 vs new K512");
  Require(a.m.id==b.m.id&&a.m.source==b.m.source&&
          a.m.st==b.m.st&&a.m.det==b.m.det&&a.m.heed==b.m.heed,
          "different original Geant4 crossing, station or Heed RNG seed");
  Require(a.prim.size()==b.prim.size(),
          "different accepted TrackHeed primary population");
  long long ions=0,charge=0;
  for(size_t i=0;i<a.prim.size();++i){
    const auto& x=a.prim[i];const auto& y=b.prim[i];
    Require(x.index==y.index&&x.cluster==y.cluster&&x.seed==y.seed&&
            x.births==y.births&&x.collected==y.collected&&
            x.endpoints==y.endpoints&&x.success==y.success,
            "different uncapped microscopic avalanche at primary "+
                std::to_string(i));
    Equal(x.u,y.u,"u");Equal(x.v,y.v,"v");
    Equal(x.w,y.w,"w");Equal(x.t,y.t,"time");
    ions+=x.births;charge+=x.collected;
  }
  std::cout<<"PASS same 2100V Geant4/TrackHeed/Garfield physics K256 vs K512"
           <<" crossing="<<a.m.id<<" primaries="<<a.prim.size()
           <<" exact_births="<<ions<<" exact_collected_e="<<charge
           <<" only independently sampled ion coordinates/weights differ\n";
  return 0;
 }catch(const std::exception& e){
  std::cerr<<e.what()<<"\n";return 2;
 }
}
