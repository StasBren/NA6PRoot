// Stage C: rebuild the SAME preselected 256 single-ion RKF templates on a
// finer grid. No new avalanche simulations and no new training-site sampling.
// Compare each 100ns stored template against a 10ns RKF recalculation
// integrated back into the ORIGINAL 100ns bins.
#include "Plane2100Response.hh"
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TNamed.h>
#include <TBranch.h>
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>
namespace {
std::string Arg(int n,char**v,const std::string&key) {
  for(int i=1;i+1<n;++i)if(key==v[i])return v[i+1];
  return "";
}
double Num(int n,char**v,const std::string&key,double def){
  const std::string s=Arg(n,v,key);
  return s.empty()?def:std::stod(s);
}
}
int main(int argc,char**argv) {
 try {
  const auto existing=Arg(argc,argv,"--input-training-root");
  const auto mobility=Arg(argc,argv,"--mobility");
  const auto prefix=Arg(argc,argv,"--out-prefix");
  const double dt=Num(argc,argv,"--dt-ns",10.);
  const double windowUs=Num(argc,argv,"--window-us",100.);
  const double maxStepMm=Num(argc,argv,"--max-step-mm",0.05);
  NA60Plane::Require(!existing.empty()&&!mobility.empty()&&!prefix.empty(),
                     "usage: --input-training-root OLD.root --mobility FILE "
                     "--out-prefix PREFIX --dt-ns 10 --window-us 100");
  const int bins=(int)std::lround(windowUs*1000./dt);
  NA60Plane::Require(dt>0&&bins>=2&&bins<=10000&&
       std::abs(bins*dt-windowUs*1000.)<1.e-7,
       "invalid fine time grid");
  const int ratio=(int)std::lround(100./dt);
  NA60Plane::Require(ratio>=1&&std::abs(ratio*dt-100.)<1.e-9,
       "new step must subdivide 100ns coarse step exactly");
  const std::string rootPath=prefix+"_training.root";
  const std::string checkPath=prefix+"_template_check.csv";
  NA60Plane::Require(!std::filesystem::exists(rootPath)&&
                     !std::filesystem::exists(checkPath),
                     "refusing to overwrite fine template results");

  std::unique_ptr<TFile> old(TFile::Open(existing.c_str(),"READ"));
  NA60Plane::Require(old&&!old->IsZombie(),"cannot read original training ROOT");
  TTree* source=dynamic_cast<TTree*>(old->Get("training_bank"));
  auto* originalMean=dynamic_cast<TH1D*>(old->Get("average_plus"));
  auto* originalArchive=dynamic_cast<TNamed*>(old->Get("input_ion_archive"));
  NA60Plane::Require(source&&originalMean&&originalArchive&&
       source->GetEntries()>0 && source->GetEntries()<=256,
       "missing 1..256 original RKF training templates");
  const int oldBins=originalMean->GetNbinsX();
  NA60Plane::Require(oldBins*ratio==bins &&
       std::abs(originalMean->GetXaxis()->GetXmax()-windowUs*1000.)<1.e-7,
       "old/new training time windows not equal");
  double u=0,v=0,w=0,t=0;
  int row=0,primary=0;
  std::vector<double>* oldWave=nullptr;
  for(const auto& name:{"u_cm","v_cm","w_cm","original_birth_t_ns",
      "original_ion_row","original_primary","waveform_plus_then_minus"}) {
    NA60Plane::Require(source->GetBranch(name)!=nullptr,
                       "original training bank lacks branch "+std::string(name));
  }
  source->SetBranchAddress("u_cm",&u);
  source->SetBranchAddress("v_cm",&v);
  source->SetBranchAddress("w_cm",&w);
  source->SetBranchAddress("original_birth_t_ns",&t);
  source->SetBranchAddress("original_ion_row",&row);
  source->SetBranchAddress("original_primary",&primary);
  source->SetBranchAddress("waveform_plus_then_minus",&oldWave);

  NA60Plane::Transport rkf(mobility,dt,bins,maxStepMm);
  std::ofstream check(checkPath);
  NA60Plane::Require(bool(check),"cannot create template consistency CSV");
  check<<std::setprecision(15)
       <<"template_index,original_ion_row,u_cm,v_cm,"
         "plus_old_Q_fC,plus_rebinned_10ns_Q_fC,"
         "minus_old_Q_fC,minus_rebinned_10ns_Q_fC,"
         "plus_first100_old_fC,plus_first100_new_fC,"
         "minus_first100_old_fC,minus_first100_new_fC,"
         "plus_max_100ns_current_diff_fC_per_ns,"
         "minus_max_100ns_current_diff_fC_per_ns\n";

  TFile file(rootPath.c_str(),"RECREATE",
             "Stage C identical 256 single-ion RKF training sites at fine dt",4);
  NA60Plane::Require(!file.IsZombie(),"could not create fine ROOT");
  TNamed inputArchive("input_ion_archive",originalArchive->GetTitle());
  TNamed sourceBank("original_training_root",existing.c_str());
  TNamed note("scope",
    "Exactly SAME original training positions and per-ion RKF drift; "
    "fine temporal recomputation only. Whole cathode planes, positive ions, "
    "Ar+ in Ar provisional mobility for Ar/CO2. No new ion births, "
    "no Fast4 avalanche model, no stereo strips, no electronics.");
  inputArchive.Write();sourceBank.Write();note.Write();

  TTree bank("training_bank","Original same-ion positions with fine time RKF");
  std::vector<double> wave;
  double uOut=0,vOut=0,wOut=0,tOut=0;
  int rowOut=0,primaryOut=0;
  bank.Branch("original_ion_row",&rowOut);
  bank.Branch("original_primary",&primaryOut);
  bank.Branch("u_cm",&uOut);bank.Branch("v_cm",&vOut);
  bank.Branch("w_cm",&wOut);bank.Branch("original_birth_t_ns",&tOut);
  bank.Branch("waveform_plus_then_minus",&wave);

  TH1D meanPlus("average_plus",
       "Same-position mean fine ion-current response;time since birth (ns);I (fC/ns)",
       bins,0.,bins*dt);
  TH1D meanMinus("average_minus",
       "Same-position mean fine ion-current response;time since birth (ns);I (fC/ns)",
       bins,0.,bins*dt);
  std::vector<double> meanP(bins,0.),meanM(bins,0.);
  for(Long64_t index=0;index<source->GetEntries();++index) {
    NA60Plane::Require(source->GetEntry(index)>0&&oldWave,
                       "cannot load original training row");
    NA60Plane::Require((int)oldWave->size()==2*oldBins,
                       "original training response has unexpected length");
    rkf.CheckStart(u,v,w,0.);
    rkf.Drift(u,v,w,0.);
    wave.assign(2*bins,0.);
    double oldQ[2]{},newQ[2]{},old100[2]{},new100[2]{},maxI[2]{};
    for(int c=0;c<2;++c) {
      for(int k=0;k<bins;++k) {
        const double current=rkf.Current(c,k);
        wave[c*bins+k]=current;
        if(c==0)meanP[k]+=current/source->GetEntries();
        else meanM[k]+=current/source->GetEntries();
      }
      for(int coarse=0;coarse<oldBins;++coarse) {
        double sum=0.;
        for(int m=0;m<ratio;++m)
          sum+=wave[c*bins+coarse*ratio+m]/ratio;
        const double previous=(*oldWave)[c*oldBins+coarse];
        oldQ[c]+=100.*previous;
        newQ[c]+=100.*sum;
        maxI[c]=std::max(maxI[c],std::abs(sum-previous));
        if(coarse==0){old100[c]=100.*previous;new100[c]=100.*sum;}
      }
    }
    check<<index<<","<<row<<","<<u<<","<<v<<","
         <<oldQ[0]<<","<<newQ[0]<<","<<oldQ[1]<<","<<newQ[1]<<","
         <<old100[0]<<","<<new100[0]<<","<<old100[1]<<","<<new100[1]<<","
         <<maxI[0]<<","<<maxI[1]<<"\n";
    uOut=u;vOut=v;wOut=w;tOut=t;rowOut=row;primaryOut=primary;
    bank.Fill();
    if((index+1)%32==0)
      std::cout<<"rebuilt_same_position_RKF="<<index+1<<"/"
               <<source->GetEntries()<<"\n";
  }
  for(int k=0;k<bins;++k) {
    meanPlus.SetBinContent(k+1,meanP[k]);
    meanMinus.SetBinContent(k+1,meanM[k]);
  }
  file.cd();bank.Write();meanPlus.Write();meanMinus.Write();
  file.Close();check.close();
  std::cout<<"PASS fine single-ion templates rebuilt WITHOUT moving any "
    "training position, n="<<source->GetEntries()<<" dt_ns="<<dt
    <<" bins="<<bins<<" old_grid_100ns_subdivision="<<ratio<<"\n"
    <<"Results "<<rootPath<<" and "<<checkPath<<"\n"
    <<"NOTE: numerical construction complete; investigate old-vs-new first-100ns "
    "differences before interpreting event timing.\n";
  return 0;
 }catch(const std::exception&e){
  std::cerr<<"FINE TEMPLATE REBUILD FAILED: "<<e.what()<<"\n";
  return 2;
 }
}
