// Standalone audited FULL-MUON prompt-electron current replay (two WHOLE planes).
// Input = one frozen Geant4 + Garfield K512 microscopic ion-birth ROOT archive.
// During microscopic replay EVERY Heed primary must exactly match the archived
// Garfield seed, ion callback census, collected electron count, and start.
// The microphysics archive is read-only. No positive-ion-current approximation.
#pragma once
#include <TFile.h>
#include <TTree.h>
#include <TNamed.h>
#include <TBranch.h>
#include <TH1D.h>
#include "Garfield/Sensor.hh"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace NA60Electron {
inline void Must(bool condition,const std::string&msg) {
  if(!condition)throw std::runtime_error("FULL MUON ELECTRON: "+msg);
}
struct Expected {
  int primary=-1,seed=0,successful=0;
  long long births=0,collected=0;
  double u=0,v=0,w=0,t=0;
};
class Archive {
  std::unique_ptr<TFile> truth_,out_;
  std::vector<Expected> expected_;
  std::vector<double> plus_,minus_;
  std::string expectedId_,sourcePath_,outPath_,crossingId_;
  int sourceStation_=0,sourceDetector_=0,sourceHeedSeed_=0;
  long long referenceBirths_=0,referenceCollected_=0;
  int seen_=0;double maxEndTime_=0;
  static constexpr double dt_=0.05;
  static constexpr int bins_=4000;
public:
  Archive(const std::string&reference,const std::string&output):
      sourcePath_(reference),outPath_(output),plus_(bins_,0.),minus_(bins_,0.) {
    Must(!std::filesystem::exists(output),"output already exists "+output);
    truth_.reset(TFile::Open(reference.c_str(),"READ"));
    Must(truth_&&!truth_->IsZombie(),"cannot open archived microscopic ROOT "+reference);
    auto* mu=dynamic_cast<TTree*>(truth_->Get("muons"));
    auto* prim=dynamic_cast<TTree*>(truth_->Get("primaries"));
    Must(mu&&prim&&mu->GetEntries()==1,"missing or non-unique source muon/primaries tree");
    std::string* id=nullptr;
    mu->SetBranchAddress("crossing_id",&id);
    mu->SetBranchAddress("station",&sourceStation_);
    mu->SetBranchAddress("det_id",&sourceDetector_);
    mu->SetBranchAddress("heed_seed",&sourceHeedSeed_);
    Must(mu->GetEntry(0)>0&&id,"source muon entry incomplete");
    expectedId_=*id;
    Expected row;
    prim->SetBranchAddress("primary_index",&row.primary);
    prim->SetBranchAddress("garfield_seed",&row.seed);
    prim->SetBranchAddress("successful",&row.successful);
    prim->SetBranchAddress("ion_births_exact",&row.births);
    prim->SetBranchAddress("collected_electrons",&row.collected);
    prim->SetBranchAddress("u_patch_cm",&row.u);
    prim->SetBranchAddress("v_absolute_cm",&row.v);
    prim->SetBranchAddress("w_relative_cm",&row.w);
    prim->SetBranchAddress("t_primary_ns",&row.t);
    for(Long64_t i=0;i<prim->GetEntries();++i) {
      Must(prim->GetEntry(i)>0&&row.primary==(int)i&&row.successful==1,
           "source contains missing or failed Heed primary");
      Must(row.births>=0&&row.collected>=0,"negative microscopic census");
      expected_.push_back(row);
      referenceBirths_+=row.births;referenceCollected_+=row.collected;
    }
    Must(!expected_.empty(),"source has no successful Heed primaries");
    out_.reset(TFile::Open(output.c_str(),"RECREATE","NA60+ direct electron plane response",4));
    Must(out_&&!out_->IsZombie(),"cannot create output ROOT "+output);
  }
  void BeginMuon(const std::string&id,int station,int detector,int heedSeed,int nPrimaries) {
    Must(id==expectedId_ && station==sourceStation_ &&
         detector==sourceDetector_ && heedSeed==sourceHeedSeed_ &&
         nPrimaries==(int)expected_.size(),
         "replayed Geant4 crossing or Heed primary count does not match frozen source");
    crossingId_=id;
  }
  void BeginPrimary(int index,int seed,double u,double v,double w,double t) {
    Must(index==seen_&&index>=0&&index<(int)expected_.size(),
         "wrong or duplicate Heed primary index");
    const auto&e=expected_[index];
    Must(seed==e.seed && std::abs(u-e.u)<1.e-9 &&
         std::abs(v-e.v)<1.e-9 && std::abs(w-e.w)<1.e-9 &&
         std::abs(t-e.t)<1.e-9,
         "replayed microscopic seed differs from frozen ion-birth reference");
  }
  void EndPrimary(int index,long long births,long long collected,bool success,
                  double maxElectronEndpointT,Garfield::Sensor&sensor) {
    Must(index==seen_&&index<(int)expected_.size(),"wrong finished primary index");
    const auto&e=expected_[index];
    Must(success && births==e.births && collected==e.collected,
         "microscopic ion/electron census changed when enabling electron currents");
    Must(std::isfinite(maxElectronEndpointT)&&maxElectronEndpointT>=0 &&
         maxElectronEndpointT<dt_*bins_,
         "electron endpoint extends beyond 200 ns recording window");
    maxEndTime_=std::max(maxEndTime_,maxElectronEndpointT);
    for(int bin=0;bin<bins_;++bin) {
      const double plus=sensor.GetElectronSignal("cathode_plus",bin);
      const double minus=sensor.GetElectronSignal("cathode_minus",bin);
      Must(std::isfinite(plus)&&std::isfinite(minus),
           "nonfinite Garfield direct electron current");
      plus_[bin]+=plus;minus_[bin]+=minus;
    }
    ++seen_;
    std::cout<<"AUDITED direct electron primary "<<seen_<<"/"
             <<expected_.size()<<" births="<<births
             <<" collected_electrons="<<collected<<"\n";
  }
  void Finish() {
    Must(seen_==(int)expected_.size() && crossingId_==expectedId_,
         "not every Heed primary contributed an electron current");
    out_->cd();
    TNamed provenance("source_microscopic_ion_archive",sourcePath_.c_str());
    TNamed limitations("scope",
      "Whole P3 noFW cathode planes ONLY; direct microscopic Garfield++ electron "
      "current all archived primaries. No ion signal in this ROOT file, no strips, "
      "no electronics/shaping. Same Geant4 crossing and exact per-primary "
      "Garfield ion/electron census verified against read-only archive. "
      "200 ns electron window; 0.05 ns bins.");
    provenance.Write(); limitations.Write();
    TTree meta("muon","Audited single complete Geant4 muon");
    std::string id=crossingId_;
    int station=sourceStation_,det_id=sourceDetector_,heedSeed=sourceHeedSeed_;
    int nPrimaries=(int)expected_.size();
    long long births=referenceBirths_,collected=referenceCollected_;
    double dt=dt_,maxEnd=maxEndTime_,window=dt_*bins_;
    meta.Branch("crossing_id",&id);
    meta.Branch("station",&station);
    meta.Branch("det_id",&det_id);
    meta.Branch("heed_seed",&heedSeed);
    meta.Branch("heed_primaries",&nPrimaries);
    meta.Branch("exact_ion_births",&births);
    meta.Branch("collected_electrons",&collected);
    meta.Branch("dt_ns",&dt);
    meta.Branch("time_window_ns",&window);
    meta.Branch("max_electron_endpoint_ns",&maxEnd);
    meta.Fill();meta.Write();
    TTree signals("electron_signal","FULL MUON direct Garfield electron induced currents");
    int bin=0;double time=0,plus=0,minus=0;
    signals.Branch("time_bin",&bin);
    signals.Branch("time_center_ns",&time);
    signals.Branch("plus_electron_fC_per_ns",&plus);
    signals.Branch("minus_electron_fC_per_ns",&minus);
    double qp=0,qm=0;
    for(bin=0;bin<bins_;++bin) {
      time=(bin+.5)*dt_;
      plus=plus_[bin];minus=minus_[bin];
      signals.Fill();qp+=dt_*plus;qm+=dt_*minus;
    }
    signals.Write();
    TH1D hp("electron_plus","NA60+ direct muon electron signal;time (ns);I_{e,+} (fC/ns)",
             bins_,0,dt_*bins_);
    TH1D hm("electron_minus","NA60+ direct muon electron signal;time (ns);I_{e,-} (fC/ns)",
             bins_,0,dt_*bins_);
    for(int j=0;j<bins_;++j){hp.SetBinContent(j+1,plus_[j]);hm.SetBinContent(j+1,minus_[j]);}
    hp.Write();hm.Write();
    out_->Write();out_->Close();
    std::cout<<"PASS audited COMPLETE-MUON direct ELECTRON-only plane currents."
             <<" crossing="<<crossingId_<<" primaries="<<seen_
             <<" ions_checked="<<referenceBirths_
             <<" collected_electrons_checked="<<referenceCollected_
             <<" plus_electron_Q_fC="<<qp<<" minus_electron_Q_fC="<<qm
             <<" electron_max_endpoint_ns="<<maxEndTime_
             <<" ROOT="<<outPath_<<"\n";
  }
};
} // namespace NA60Electron
