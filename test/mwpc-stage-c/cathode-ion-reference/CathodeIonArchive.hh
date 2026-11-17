// Isolated 2100 V cathode-response INPUT reference (not a strip-current model).
// ROOT trees: immutable Geant4 muon, Heed cluster/primary, actual ion births,
// electron endpoints, and a one-primary direct whole-cathode electron signal.
// Each ion belongs to a concrete original Heed primary via primary_index.
// No physical positive-ion species is inferred from the ionisation callback.
#pragma once
#include <TFile.h>
#include <TTree.h>
#include <array>
#include <cmath>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace CathodeRef {
class Archive {
  std::unique_ptr<TFile> f_;
  TTree *muons_=nullptr,*clusters_=nullptr,*primaries_=nullptr,
        *ions_=nullptr,*ends_=nullptr,*electronSignal_=nullptr;
  std::string mode_,source_,crossingId_;
  int reservoir_=0,station_=0,det_=0,heedSeed_=0;
  int chunk_=0,event_=0,hit_=0,track_=0,pdg_=0;
  double in_[3]{},out_[3]{},mom_[3]{};
  double thetaU_=0,thetaW_=0,path_=0,momentum_=0,edep_=0;
  double pitch_=0.3,phaseMm_=0,phaseCm_=0,centerGlobalCm_=0;
  int baseWire_=0;
  int clusterIndex_=0,primaryIndex_=0,garfieldSeed_=0,success_=0;
  double primaryU_=0,primaryV_=0,primaryW_=0,primaryTime_=0;
  double primaryEnergy_=0;
  long long ionTotal_=0, ionSaved_=0, electronCollected_=0,
            electronEndpoints_=0;
  double inversePi_=0;
  // ion row: patch position retained alongside true physical position.
  long long ordinal_=0;
  int ionLocalWire_=0,ionGlobalWire_=0;
  double ionU_=0,ionV_=0,ionW_=0,ionT_=0;
  double globalU_=0,globalW_=0,wireCenterMm_=0,ionWeight_=0,ionPi_=0;
  // Direct whole-plane electron-signal control for first Heed primary only.
  int endpointStatus_=0;
  double eX0_=0,eY0_=0,eZ0_=0,eT0_=0,eE0_=0,
         eX1_=0,eY1_=0,eZ1_=0,eT1_=0,eE1_=0;
  int electronBin_=0;
  double electronTimeNs_=0,electronPlus_=0,electronMinus_=0;
  int acceptedClusters_=0,acceptedElectrons_=0;
  double clusterU_=0,clusterV_=0,clusterW_=0,clusterT_=0;
  void I(TTree* t,const char* key,int& v){
    t->Branch(key,&v,(std::string(key)+"/I").c_str());
  }
  void L(TTree* t,const char* key,long long& v){
    t->Branch(key,&v,(std::string(key)+"/L").c_str());
  }
  void D(TTree* t,const char* key,double& v){
    t->Branch(key,&v,(std::string(key)+"/D").c_str());
  }
  static void Ensure(bool b,const std::string& what) {
    if(!b)throw std::runtime_error("Cathode ion archive: "+what);
  }
 public:
  Archive(const std::string& path,const std::string& mode,int reservoir)
      : mode_(mode),reservoir_(reservoir) {
    Ensure((mode=="weighted" && (reservoir==256||reservoir==512||reservoir==1024)) ||
           (mode=="full"&&reservoir==0),"invalid mode/reservoir");
    f_.reset(TFile::Open(path.c_str(),"RECREATE","Uncapped cathode ion truth",4));
    Ensure(bool(f_)&&!f_->IsZombie(),"ROOT file creation failed: "+path);
    f_->SetCompressionLevel(4);
    muons_=new TTree("muons","Geant4 immutable chamber-local coordinates and momentum");
    clusters_=new TTree("clusters","TrackHeed accepted cluster provenance");
    primaries_=new TTree("primaries","Original Heed seeds with exact callback birth census");
    ions_=new TTree("ions","True Garfield ionisation callback coordinates/times; no ion species identification");
    ends_=new TTree("electron_endpoints","First-Heed-primary full electron endpoints (full mode only)");
    electronSignal_=new TTree("electron_control","Direct Garfield whole-cathode electron-only current; NOT stereo strips");
    for(TTree* t:{muons_,clusters_,primaries_,ions_,ends_,electronSignal_})
      t->SetDirectory(f_.get());
    ions_->SetAutoFlush(-16LL*1024*1024); // limit ROOT buffered basket memory
    ends_->SetAutoFlush(-8LL*1024*1024);
    muons_->Branch("crossing_id",&crossingId_);
    muons_->Branch("source_channel",&source_);
    muons_->Branch("storage_mode",&mode_);
    I(muons_,"reservoir_per_primary",reservoir_);
    I(muons_,"source_chunk",chunk_);I(muons_,"event",event_);
    I(muons_,"hit_index",hit_);I(muons_,"track_id",track_);
    I(muons_,"station",station_);I(muons_,"det_id",det_);
    I(muons_,"pdg",pdg_);I(muons_,"heed_seed",heedSeed_);
    I(muons_,"base_global_wire",baseWire_);
    D(muons_,"pitch_cm",pitch_);D(muons_,"station_phase_mm",phaseMm_);
    D(muons_,"patch_shift_cm",centerGlobalCm_);
    D(muons_,"theta_u_deg",thetaU_);D(muons_,"theta_w_deg",thetaW_);
    D(muons_,"path_length_mm",path_);D(muons_,"momentum_GeV",momentum_);
    D(muons_,"edep_eV",edep_);
    for(int j=0;j<3;++j) {
      const char* k[3]={"u","v","w"};
      D(muons_,(std::string(k[j])+"_in_mm").c_str(),in_[j]);
      D(muons_,(std::string(k[j])+"_out_mm").c_str(),out_[j]);
      D(muons_,(std::string("p")+k[j]+"_in_GeV").c_str(),mom_[j]);
    }
    I(clusters_,"cluster_index",clusterIndex_);
    D(clusters_,"u_relative_cm",clusterU_);D(clusters_,"v_absolute_cm",clusterV_);
    D(clusters_,"w_relative_cm",clusterW_);D(clusters_,"t_ns",clusterT_);
    I(clusters_,"n_original",acceptedClusters_);
    I(clusters_,"n_accepted",acceptedElectrons_);
    I(primaries_,"primary_index",primaryIndex_);
    I(primaries_,"cluster_index",clusterIndex_);
    I(primaries_,"garfield_seed",garfieldSeed_);
    I(primaries_,"successful",success_);
    D(primaries_,"u_patch_cm",primaryU_);
    D(primaries_,"v_absolute_cm",primaryV_);
    D(primaries_,"w_relative_cm",primaryW_);
    D(primaries_,"t_primary_ns",primaryTime_);
    D(primaries_,"start_energy_eV",primaryEnergy_);
    L(primaries_,"ion_births_exact",ionTotal_);
    L(primaries_,"ion_births_saved",ionSaved_);
    L(primaries_,"collected_electrons",electronCollected_);
    L(primaries_,"electron_endpoints",electronEndpoints_);
    D(primaries_,"inverse_inclusion_probability",inversePi_);
    I(ions_,"primary_index",primaryIndex_);I(ions_,"cluster_index",clusterIndex_);
    I(ions_,"garfield_seed",garfieldSeed_);
    L(ions_,"callback_ordinal",ordinal_);
    I(ions_,"local_wire_index",ionLocalWire_);
    I(ions_,"global_wire_index",ionGlobalWire_);
    D(ions_,"wire_u_global_mm",wireCenterMm_);
    D(ions_,"u_patch_cm",ionU_);D(ions_,"v_absolute_cm",ionV_);
    D(ions_,"w_relative_cm",ionW_);D(ions_,"t_birth_ns",ionT_);
    D(ions_,"u_chamber_cm",globalU_);
    D(ions_,"w_chamber_cm",globalW_);
    D(ions_,"inclusion_probability",ionPi_);
    D(ions_,"inverse_inclusion_weight",ionWeight_);
    I(ends_,"primary_index",primaryIndex_);
    I(ends_,"endpoint_status",endpointStatus_);
    D(ends_,"u0_patch_cm",eX0_);D(ends_,"v0_cm",eY0_);
    D(ends_,"w0_relative_cm",eZ0_);D(ends_,"t0_ns",eT0_);
    D(ends_,"energy0_eV",eE0_);
    D(ends_,"u1_patch_cm",eX1_);D(ends_,"v1_cm",eY1_);
    D(ends_,"w1_relative_cm",eZ1_);D(ends_,"t1_ns",eT1_);
    D(ends_,"energy1_eV",eE1_);
    I(electronSignal_,"primary_index",primaryIndex_);
    I(electronSignal_,"time_bin",electronBin_);
    D(electronSignal_,"time_center_ns",electronTimeNs_);
    D(electronSignal_,"whole_cathode_plus_electron_fC_per_ns",electronPlus_);
    D(electronSignal_,"whole_cathode_minus_electron_fC_per_ns",electronMinus_);
  }
  template<class Crossing>
  void Muon(const Crossing& h,const std::string& id,double patchPhaseCm,
            double physicalPitchMm,int baseGlobalWire) {
    crossingId_=id;source_=h.source;chunk_=h.chunk;event_=h.event;
    hit_=h.hit;track_=h.track;station_=h.station;det_=h.det;
    pdg_=h.pdg;heedSeed_=static_cast<int>(h.heedSeed);
    thetaU_=h.thetaU;thetaW_=h.thetaW;path_=h.length;
    momentum_=h.momentum;edep_=h.edep;pitch_=physicalPitchMm*.1;
    phaseCm_=patchPhaseCm;baseWire_=baseGlobalWire;
    phaseMm_=(station_==0?1.5:0.);
    // Patch x=0 is the real anode nearest the physical Geant4 entrance.
    centerGlobalCm_=.1*h.a[0]-patchPhaseCm;
    const double expected=.1*(phaseMm_+baseWire_*physicalPitchMm);
    Ensure(std::abs(centerGlobalCm_-expected)<1.e-8,
           "station-aware physical wire offset/phase mismatch");
    for(int j=0;j<3;++j){in_[j]=h.a[j];out_[j]=h.b[j];mom_[j]=h.p[j];}
    muons_->Fill();
  }
  template<class Cluster>
  void ClusterEntry(const Cluster& c,int index) {
    clusterIndex_=index;clusterU_=c.u;clusterV_=c.v;
    clusterW_=c.w;clusterT_=c.t;
    acceptedClusters_=c.nOriginal;acceptedElectrons_=c.nKept;
    clusters_->Fill();
  }
  template<class Primary>
  void BeginPrimary(const Primary& p,int index,int garfieldSeed) {
    primaryIndex_=index;clusterIndex_=p.clusterIndex;garfieldSeed_=garfieldSeed;
    primaryU_=p.u+phaseCm_;primaryV_=p.v;primaryW_=p.w;
    primaryTime_=p.t;primaryEnergy_=0.1;
  }
  bool Full()const{return mode_=="full";}
  int Reservoir()const{return reservoir_;}
  void Ion(double u,double v,double w,double t,long long ordinal,
           double inclusionProbability,double weight) {
    Ensure(std::isfinite(u)&&std::isfinite(v)&&std::isfinite(w)&&
           std::isfinite(t)&&inclusionProbability>0.&&inclusionProbability<=1.,
           "invalid physical ion birth coordinates/inclusion probability");
    ionU_=u;ionV_=v;ionW_=w;ionT_=t;ordinal_=ordinal;
    ionLocalWire_=int(std::llround(u/pitch_));
    Ensure(ionLocalWire_>=-4&&ionLocalWire_<=4,
           "callback birth outside nine-wire patch");
    ionGlobalWire_=baseWire_+ionLocalWire_;
    wireCenterMm_=phaseMm_+ionGlobalWire_*10.*pitch_;
    globalU_=centerGlobalCm_+u;globalW_=.1*in_[2]+w;
    Ensure(std::abs(10.*globalU_-wireCenterMm_-
                10.*(u-ionLocalWire_*pitch_))<1.e-8,
           "local-to-physical position mapping mismatch");
    ionPi_=inclusionProbability;ionWeight_=weight;
    ions_->Fill();
  }
  void PrimaryResult(long long nIons,long long nSaved,long long collected,
                     long long endpoints,bool successful) {
    Ensure(nIons>=0&&nSaved>=0&&nSaved<=nIons&&collected>=0,
           "invalid primary census");
    Ensure((mode_=="full" && nSaved==nIons) ||
           (mode_=="weighted" && nSaved==std::min<long long>(nIons,reservoir_)),
           "recorded ion count disagrees with callback census");
    ionTotal_=nIons;ionSaved_=nSaved;electronCollected_=collected;
    electronEndpoints_=endpoints;success_=successful?1:0;
    inversePi_=nSaved?double(nIons)/double(nSaved):0.;
    primaries_->Fill();
  }
  void Endpoint(double x0,double y0,double z0,double t0,double e0,
                double x1,double y1,double z1,double t1,double e1,int status) {
    eX0_=x0;eY0_=y0;eZ0_=z0;eT0_=t0;eE0_=e0;
    eX1_=x1;eY1_=y1;eZ1_=z1;eT1_=t1;eE1_=e1;
    endpointStatus_=status;ends_->Fill();
  }
  void ElectronBin(int bin,double t,double plus,double minus) {
    electronBin_=bin;electronTimeNs_=t;
    electronPlus_=plus;electronMinus_=minus;electronSignal_->Fill();
  }
  void Finish() {
    Ensure(muons_->GetEntries()==1,"detailed archive must contain exactly one muon");
    f_->cd();
    muons_->Write();clusters_->Write();primaries_->Write();
    ions_->Write();ends_->Write();electronSignal_->Write();
    f_->Write();f_->Close();
  }
};
} // namespace CathodeRef
