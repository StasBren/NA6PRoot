// ROOT/C++ unbiased reservoir convergence against SAME-MUON complete callbacks.
// Compares conditional primary-by-primary moments and local collecting-wire
// birth associations, reporting exact finite-population sampling uncertainty.
// Input: full/ion_births.root weighted/ion_births.root (same frozen muon).
// This does not test physical stereo-strip waveforms or ion transport.
#include <TFile.h>
#include <TTree.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <string>
namespace {
void Must(bool ok,const std::string& msg){
 if(!ok)throw std::runtime_error("RESERVOIR CONVERGENCE: "+msg);
}
struct Summary {
 long long n=0;
 double weight=0;
 std::array<double,4> sum{},sum2{};
 std::array<double,9> nByWire{};
 void Add(const std::array<double,4>& x,int wire,double wt) {
   Must(wire>=-4&&wire<=4,"nearest local wire invalid");
   ++n;weight+=wt;
   for(int k=0;k<4;++k) {
      sum[k]+=wt*x[k];sum2[k]+=wt*x[k]*x[k];
   }
   nByWire[wire+4]+=wt;
 }
};
struct Loaded{
 std::string id,mode;
 int K=0;
 std::map<int,long long> totals,keep;
 std::map<int,Summary> stats;
 long long total=0;
};
Loaded Read(const char* path){
 std::unique_ptr<TFile> file(TFile::Open(path,"READ"));
 Must(file&&!file->IsZombie(),"could not open ROOT file "+std::string(path));
 auto* muons=dynamic_cast<TTree*>(file->Get("muons"));
 auto* prim=dynamic_cast<TTree*>(file->Get("primaries"));
 auto* ions=dynamic_cast<TTree*>(file->Get("ions"));
 Must(muons&&prim&&ions&&muons->GetEntries()==1,"missing fundamental ROOT trees");
 std::string *id=nullptr,*mode=nullptr;int K=0;
 muons->SetBranchAddress("crossing_id",&id);
 muons->SetBranchAddress("storage_mode",&mode);
 muons->SetBranchAddress("reservoir_per_primary",&K);
 muons->GetEntry(0);
 Loaded L;L.id=*id;L.mode=*mode;L.K=K;
 int pi=0;long long n=0,saved=0;
 prim->SetBranchAddress("primary_index",&pi);
 prim->SetBranchAddress("ion_births_exact",&n);
 prim->SetBranchAddress("ion_births_saved",&saved);
 for(Long64_t i=0;i<prim->GetEntries();++i) {
   prim->GetEntry(i);
   Must(n>=0&&saved>=0&&saved<=n,"invalid exact Heed primary census");
   Must(L.totals.emplace(pi,n).second,"duplicate Heed primary");
   L.keep[pi]=saved;L.total+=n;
 }
 int iw=0,local=0;
 double u=0,v=0,w=0,t=0,weight=0;
 ions->SetBranchAddress("primary_index",&iw);
 ions->SetBranchAddress("local_wire_index",&local);
 ions->SetBranchAddress("u_patch_cm",&u);
 ions->SetBranchAddress("v_absolute_cm",&v);
 ions->SetBranchAddress("w_relative_cm",&w);
 ions->SetBranchAddress("t_birth_ns",&t);
 ions->SetBranchAddress("inverse_inclusion_weight",&weight);
 for(Long64_t i=0;i<ions->GetEntries();++i){
   ions->GetEntry(i);
   Must(L.totals.count(iw),"birth refers to unknown primary");
   L.stats[iw].Add({u,v,w,t},local,weight);
 }
 for(const auto& [p,N]:L.totals){
   const auto& s=L.stats[p];
   Must(s.n==L.keep[p],"ROOT row count and parent primary count mismatch");
   Must(std::abs(s.weight-double(N))<1e-7*std::max(1.,double(N)),
        "Horvitz-Thompson birth census does not reconstruct exact count");
 }
 return L;
}
}
int main(int argc,char**argv){
 try{
   Must(argc==3,"usage compare_sampling_to_full_root full.root sampled.root");
   const auto full=Read(argv[1]), sampled=Read(argv[2]);
   Must(full.id==sampled.id&&full.mode=="full"&&sampled.mode=="weighted",
        "must compare same crossing, full versus weighted");
   Must(full.K==0&&(sampled.K==256||sampled.K==512||sampled.K==1024),
        "wrong source archive sampling design");
   Must(full.totals==sampled.totals,"Garfield callback census changed between modes");
   Must(full.total>0,"no ion statistics for convergence test");
   double fpcDiag=0;
   const char* dim[]={"u_patch_cm","v_absolute_cm","w_relative_cm","t_birth_ns"};
   std::cout<<"NA60+ Simulation - ROOT/C++ unbiased ion reservoir convergence\n"
            <<"crossing="<<full.id<<" reservoir_K="<<sampled.K
            <<" total_exact_births="<<full.total<<"\n";
   for(int j=0;j<4;++j){
     double exact=0,estimated=0,variance=0;
     for(const auto& [p,N]:full.totals){
       if(N==0)continue;
       const auto& F=full.stats.at(p),&S=sampled.stats.at(p);
       const long long n=S.n;
       Must(F.n==N&&n==std::min<long long>(N,sampled.K),
            "parent conditional sample size wrong");
       const double mu=F.sum[j]/N;
       const double m2=F.sum2[j]/N;
       const double sampleVar=std::max(0.,m2-mu*mu);
       const double trueVariance=N>1?sampleVar*N/(N-1.):0.;
       // Variance of Horvitz-Thompson global moment estimator, fixed
       // exact n=K per primary, sampled without replacement.
       variance+=double(N)*double(N)*(1.-double(n)/N)*
                 trueVariance/n;
       exact+=F.sum[j];estimated+=S.sum[j];
     }
     exact/=full.total;estimated/=full.total;
     const double se=std::sqrt(variance)/full.total;
     const double z=se>0?(estimated-exact)/se:0.;
     std::cout<<dim[j]<<": full="<<exact<<" sampled="<<estimated
              <<" delta="<<(estimated-exact)
              <<" exact_design_SE="<<se<<" z="<<z<<"\n";
     if(se==0)Must(std::abs(estimated-exact)<1e-9,
                  "zero-variance conditional mean differs");
     fpcDiag+=z*z;
   }
   for(int wire=-4;wire<=4;++wire){
     double exact=0,estimated=0,var=0;
     for(const auto& [p,N]:full.totals){
       if(N==0)continue;
       const double fullWire=full.stats.at(p).nByWire[wire+4];
       const auto& S=sampled.stats.at(p);
       const double probability=fullWire/N;
       const double sampleVar=N>1?
           probability*(1.-probability)*N/(N-1.):0.;
       var+=double(N)*N*(1.-double(S.n)/N)*sampleVar/S.n;
       exact+=fullWire;estimated+=S.nByWire[wire+4];
     }
     const double se=std::sqrt(var),z=se>0?(estimated-exact)/se:0.;
     std::cout<<"nearest_local_wire="<<wire
              <<" true_births="<<exact<<" HT_births="<<estimated
              <<" SE="<<se<<" z="<<z<<"\n";
   }
   // Finite-population z scores are diagnostics, not a fabricated proof of
   // waveform accuracy; small-K tails and boundary conditional bins need
   // independent replicates and future weighted Ramo-current validation.
   std::cout<<"PASS exact counts/source identity; reported conditional ROOT "
               "sampling deviation and finite-population uncertainty\n";
   return 0;
 }catch(const std::exception& e){
   std::cerr<<e.what()<<"\n";return 2;
 }
}
