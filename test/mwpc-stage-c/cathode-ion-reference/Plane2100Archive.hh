// Read the frozen one-Geant4-muon ROOT ion-birth archive without materializing
// millions of ion records. Preserve exact primary populations/weights.
#pragma once
#include <TFile.h>
#include <TTree.h>
#include <algorithm>
#include <cmath>
#include <memory>
#include <string>
#include <vector>
#include "Plane2100Response.hh"

namespace NA60Plane {
class Archive {
  std::unique_ptr<TFile> f_;
  TTree* ions_=nullptr;
  int primary_=0;
  double u_=0,v_=0,w_=0,t_=0,wt_=0;
  std::vector<long long> exact_,saved_;
  std::vector<long long> seen_;
  std::vector<double> summedWeights_;
public:
  std::string id,source,mode;
  int station=0,det=0,K=0;
  long long nTrue=0,nSaved=0;
  Archive(const std::string& name) {
    f_.reset(TFile::Open(name.c_str(),"READ"));
    Require(f_&&!f_->IsZombie(),"cannot open ROOT archive "+name);
    auto* mu=dynamic_cast<TTree*>(f_->Get("muons"));
    auto* prim=dynamic_cast<TTree*>(f_->Get("primaries"));
    ions_=dynamic_cast<TTree*>(f_->Get("ions"));
    Require(mu&&prim&&ions_&&mu->GetEntries()==1,"missing muons/primaries/ions tree");
    std::string *a=nullptr,*b=nullptr,*c=nullptr;
    mu->SetBranchAddress("crossing_id",&a);
    mu->SetBranchAddress("source_channel",&b);
    mu->SetBranchAddress("storage_mode",&c);
    mu->SetBranchAddress("station",&station);
    mu->SetBranchAddress("det_id",&det);
    mu->SetBranchAddress("reservoir_per_primary",&K);
    mu->GetEntry(0);
    id=*a;source=*b;mode=*c;
    Require(mode=="weighted"||mode=="full","unknown source storage mode");
    int p=0,ok=0;long long N=0,n=0;
    prim->SetBranchAddress("primary_index",&p);
    prim->SetBranchAddress("ion_births_exact",&N);
    prim->SetBranchAddress("ion_births_saved",&n);
    prim->SetBranchAddress("successful",&ok);
    for(Long64_t j=0;j<prim->GetEntries();++j) {
      prim->GetEntry(j);
      Require(p==(int)j && ok==1 && N>=0 && n>=0 && n<=N,
              "invalid/failed microscopic primary");
      if(mode=="full")Require(N==n,"FULL ion file is incomplete");
      else Require(K>0&&n==std::min<long long>(N,K),
                   "weighted reservoir length disagrees with true population");
      exact_.push_back(N);saved_.push_back(n);
      nTrue+=N;nSaved+=n;
    }
    seen_.assign(exact_.size(),0);
    summedWeights_.assign(exact_.size(),0.);
    Require(ions_->GetEntries()==nSaved,"ion ROOT row count disagrees with primaries");
    ions_->SetBranchAddress("primary_index",&primary_);
    ions_->SetBranchAddress("u_patch_cm",&u_);
    ions_->SetBranchAddress("v_absolute_cm",&v_);
    ions_->SetBranchAddress("w_relative_cm",&w_);
    ions_->SetBranchAddress("t_birth_ns",&t_);
    ions_->SetBranchAddress("inverse_inclusion_weight",&wt_);
  }
  long long Entries()const{return nSaved;}
  struct Birth {double u,v,w,t,weight;int primary;};
  Birth Get(long long i) {
    Require(i>=0&&i<nSaved,"ion row index outside archive");
    ions_->GetEntry(i);
    Require(primary_>=0&&primary_<(int)exact_.size(),
            "ion row has invalid parent Heed primary");
    const auto idx=(std::size_t)primary_;
    Require(std::isfinite(wt_)&&wt_>0.,"invalid statistical weight");
    const double expected=double(exact_[idx])/double(saved_[idx]);
    Require(std::abs(wt_-expected)<1.e-8*std::max(1.,expected),
            "source sample ion weight does not match original census");
    ++seen_[idx];summedWeights_[idx]+=wt_;
    return {u_,v_,w_,t_,wt_,primary_};
  }
  void CheckComplete()const{
    for(std::size_t i=0;i<exact_.size();++i) {
      Require(seen_[i]==saved_[i],"not all sampled ion rows processed");
      Require(std::abs(summedWeights_[i]-double(exact_[i])) <
              1.e-6*std::max(1.,double(exact_[i])),
              "weighted ion census not conserved");
    }
  }
};
}
