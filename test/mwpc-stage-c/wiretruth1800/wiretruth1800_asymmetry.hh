// NA60+ Simulation -- additional exact 1800V Garfield++ per-anode ion-birth asymmetry observables.
// This is NOT a cathode-strip signal model. ALL asymmetries are calculated from
// exact ionisation callbacks assigned to the nearest physical anode.
// Physics and plots: C++/ROOT only; zero-collecting ion-bearing wires retained.
#pragma once
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TFile.h>
#include <TColor.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <limits>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
namespace MWPC1800Asym {
inline double Undefined(){return std::numeric_limits<double>::quiet_NaN();}
inline double SignedAsym(long long pos,long long neg){
  if(pos<0||neg<0)throw std::runtime_error("Negative native ion sector census");
  return pos+neg>0?double(pos-neg)/double(pos+neg):Undefined();
}
// A_hemi = (N(phi>0) - N(phi<0)) / (N(phi>0)+N(phi<0)), excluding boundary phi=0, +/-180.
// A_TB = (N(phi in [45,135)) - N(phi in [-135,-45))) / (same total), excluding side.
// Unlike A_hemi, A_TB excludes angles near +/-u. Both use BIRTH ion azimuth,
// not positive ion drift endpoints nor collected-electron polarity.
class Plots {
public:
  TH2D *wireGainA=nullptr,*seedWireGainA=nullptr,*muonGainA=nullptr,
       *muonShareA=nullptr,
       *ionYieldA=nullptr,*wireFracA=nullptr,*aHemiTB=nullptr,
       *aByLocal=nullptr,*aByImpact=nullptr,*aByThetaU=nullptr,
       *aByThetaW=nullptr,*tbGainA=nullptr;
  TH1D *zeroQIonA=nullptr,*aHemi=nullptr,*aTB=nullptr;
  TProfile *meanWireLogGainA=nullptr,*meanSeedLogGainA=nullptr,
           *meanAByLocal=nullptr,*meanAByImpact=nullptr,
           *meanAByThetaU=nullptr,*meanAByThetaW=nullptr;
  long long muonsWithDefinedA=0,muonsWithGainAndA=0;
  long long ionWireRows=0,wireWithA=0,wireWithTB=0,
            wireWithZeroQ=0,wireWithGainAndA=0;
  long long seedWireRows=0,seedWireWithA=0,seedWireWithZeroQ=0,
            seedWireWithGainAndA=0;
  explicit Plots(const std::string& scope){
    auto name=[&](const char*x){return scope+"_asym_"+x;};
    wireGainA=new TH2D(name("wire_gain_vs_A_hemi").c_str(),
       ";A_{ion,hemi} per anode;log_{10}(Q_{wire}/N_{Heed,muon})",
       64,-1,1,100,-4.5,7.5);
    seedWireGainA=new TH2D(name("primary_wire_gain_vs_A_hemi").c_str(),
       ";A_{ion,hemi} per primary and wire;log_{10}(Q_{primary,wire})",
       64,-1,1,100,-.5,7.5);
    muonGainA=new TH2D(name("muon_gain_vs_A_hemi").c_str(),
       ";A_{ion,hemi} per muon;log_{10}(Q_{muon}/N_{Heed,muon})",
       64,-1,1,100,-4.5,7.5);
    muonShareA=new TH2D(name("muon_second_wire_fraction_vs_A_hemi").c_str(),
       ";A_{ion,hemi} per muon;Q_{2}/Q_{muon}",
       64,-1,1,50,0,.5);
    ionYieldA=new TH2D(name("wire_ion_births_vs_A_hemi").c_str(),
       ";A_{ion,hemi} per anode;log_{10}(N_{ion birth,wire})",
       64,-1,1,100,-.5,9.5);
    wireFracA=new TH2D(name("wire_charge_fraction_vs_A_hemi").c_str(),
       ";A_{ion,hemi} per anode;Q_{wire}/Q_{muon}",
       64,-1,1,60,0,1.0);
    aHemiTB=new TH2D(name("A_hemi_vs_A_TB").c_str(),
       ";A_{ion,hemi};A_{ion,TB}",64,-1,1,64,-1,1);
    aByLocal=new TH2D(name("A_hemi_vs_local_wire").c_str(),
       ";Local physical anode relative to entrance-nearest wire;A_{ion,hemi}",
       9,-4.5,4.5,64,-1,1);
    aByImpact=new TH2D(name("A_hemi_vs_signed_impact").c_str(),
       ";Signed entrance impact u-u_{nearest wire} [mm];A_{ion,hemi}",
       50,-1.5,1.5,64,-1,1);
    aByThetaU=new TH2D(name("A_hemi_vs_theta_u").c_str(),
       ";Signed Geant4 #theta_{u} [deg];A_{ion,hemi}",
       72,-90,90,64,-1,1);
    aByThetaW=new TH2D(name("A_hemi_vs_theta_w").c_str(),
       ";Signed Geant4 #theta_{w} [deg];A_{ion,hemi}",
       72,-90,90,64,-1,1);
    tbGainA=new TH2D(name("wire_gain_vs_A_TB").c_str(),
       ";A_{ion,TB} per anode;log_{10}(Q_{wire}/N_{Heed,muon})",
       64,-1,1,100,-4.5,7.5);
    zeroQIonA=new TH1D(name("ion_A_zero_collected_wire").c_str(),
       ";A_{ion,hemi};Ion-bearing anodes with Q_{wire}=0",
       64,-1,1);
    aHemi=new TH1D(name("A_hemi").c_str(),
       ";A_{ion,hemi};Ion-bearing muon-anode rows",64,-1,1);
    aTB=new TH1D(name("A_TB").c_str(),
       ";A_{ion,TB};Ion-bearing muon-anode rows",64,-1,1);
    meanWireLogGainA=new TProfile(name("mean_log_wire_gain_vs_A").c_str(),
       ";A_{ion,hemi};Mean log_{10}(Q_{wire}/N_{Heed,muon})",
       32,-1,1);
    meanSeedLogGainA=new TProfile(name("mean_log_seed_wire_gain_vs_A").c_str(),
       ";A_{ion,hemi};Mean log_{10}(Q_{primary,wire})",
       32,-1,1);
    meanAByLocal=new TProfile(name("mean_A_vs_local_wire").c_str(),
       ";Local anode index relative to entrance-nearest;Mean A_{ion,hemi}",
       9,-4.5,4.5);
    meanAByImpact=new TProfile(name("mean_A_vs_impact").c_str(),
       ";Signed Geant4 impact relative to nearest wire [mm];Mean A_{ion,hemi}",
       30,-1.5,1.5);
    meanAByThetaU=new TProfile(name("mean_A_vs_theta_u").c_str(),
       ";Signed Geant4 #theta_{u} [deg];Mean A_{ion,hemi}",
       36,-90,90);
    meanAByThetaW=new TProfile(name("mean_A_vs_theta_w").c_str(),
       ";Signed Geant4 #theta_{w} [deg];Mean A_{ion,hemi}",
       36,-90,90);
    for(TH1* h:All()) {
      h->SetDirectory(nullptr);
      h->Sumw2();
    }
  }
  std::array<TH1*,21> All() const{
    return {wireGainA,seedWireGainA,muonGainA,muonShareA,ionYieldA,wireFracA,
      aHemiTB,aByLocal,aByImpact,aByThetaU,aByThetaW,tbGainA,
      zeroQIonA,aHemi,aTB,meanWireLogGainA,meanSeedLogGainA,
      meanAByLocal,meanAByImpact,meanAByThetaU,meanAByThetaW};
  }
  void Write()const{
    for(TH1* h:All())if(h)h->Write();
  }
  void FillMuon(long long ionUp,long long ionDown,long long collected,int nHeed,
                double fSecond){
    const double a=SignedAsym(ionUp,ionDown);
    if(!std::isfinite(a))return;
    ++muonsWithDefinedA;
    if(collected>0 && std::isfinite(fSecond))muonShareA->Fill(a,fSecond);
    // Undefined if no accepted Heed primaries, even if there were ions.
    if(collected>0&&nHeed>0){
      muonGainA->Fill(a,std::log10(double(collected)/nHeed));
      ++muonsWithGainAndA;
    }
  }
  void FillWire(long long ionUp,long long ionDown,
                long long top,long long bottom,
                long long q,long long ionBorn,long long qMuon,int nHeed,
                int localWire,double impactMm,double signedThetaU,double signedThetaW){
    if(ionBorn<=0)return;
    ++ionWireRows;
    if(q==0)++wireWithZeroQ;
    const double a=SignedAsym(ionUp,ionDown);
    const double tb=SignedAsym(top,bottom);
    if(std::isfinite(a)){
      ++wireWithA;
      aHemi->Fill(a);
      aByLocal->Fill(localWire,a);
      aByImpact->Fill(impactMm,a);
      aByThetaU->Fill(signedThetaU,a);
      aByThetaW->Fill(signedThetaW,a);
      meanAByLocal->Fill(localWire,a);
      meanAByImpact->Fill(impactMm,a);
      meanAByThetaU->Fill(signedThetaU,a);
      meanAByThetaW->Fill(signedThetaW,a);
      ionYieldA->Fill(a,std::log10(double(ionBorn)));
      if(qMuon>0)wireFracA->Fill(a,double(q)/qMuon);
      if(q==0)zeroQIonA->Fill(a);
      if(q>0&&nHeed>0){
        const double lg=std::log10(double(q)/nHeed);
        wireGainA->Fill(a,lg);
        meanWireLogGainA->Fill(a,lg);
        ++wireWithGainAndA;
      }
    }
    if(std::isfinite(tb)){
      ++wireWithTB;
      aTB->Fill(tb);
      if(q>0&&nHeed>0)tbGainA->Fill(tb,std::log10(double(q)/nHeed));
      if(std::isfinite(a))aHemiTB->Fill(a,tb);
    }
  }
  void FillPrimaryWire(long long up,long long down,long long q,long long ions){
    if(ions<=0)return;
    ++seedWireRows;
    if(q==0)++seedWireWithZeroQ;
    const double a=SignedAsym(up,down);
    if(std::isfinite(a)){
      ++seedWireWithA;
      if(q>0){
        const double lg=std::log10(double(q)); // one Heed electron: per-wire gain is Q
        seedWireGainA->Fill(a,lg);
        meanSeedLogGainA->Fill(a,lg);
        ++seedWireWithGainAndA;
      }
    }
  }
};
inline void DrawOne(const std::filesystem::path&dir,const std::string&ch,
                    TH2D* h,const std::string&which,const std::string&description){
  if(h->GetEntries()==0)return;
  TCanvas canvas(("c_asym_"+ch+"_"+which).c_str(),"",980,710);
  canvas.SetLeftMargin(.14);canvas.SetRightMargin(.16);canvas.SetBottomMargin(.13);
  gStyle->SetOptStat(0);
  h->Draw("COLZ");
  TLatex label;label.SetNDC();label.SetTextFont(42);
  label.SetTextSize(.036);label.DrawLatex(.16,.92,"NA60+ Simulation");
  label.SetTextSize(.025);label.DrawLatex(.16,.875,description.c_str());
  canvas.SaveAs((dir/("asym_"+which+"_"+ch+".png")).string().c_str());
}
inline void DrawOneD(const std::filesystem::path&dir,const std::string&ch,
                     TH1D* h,const std::string&which,const std::string&desc){
  if(h->GetEntries()==0)return;
  TCanvas canvas(("c_asym_1d_"+ch+"_"+which).c_str(),"",960,680);
  canvas.SetLeftMargin(.14);canvas.SetBottomMargin(.14);
  gStyle->SetOptStat(0);
  h->SetLineColor(kBlue+1);h->SetLineWidth(2);
  h->Draw("HIST");
  TLatex label;label.SetNDC();label.SetTextFont(42);
  label.SetTextSize(.036);label.DrawLatex(.16,.92,"NA60+ Simulation");
  label.SetTextSize(.024);label.DrawLatex(.16,.875,desc.c_str());
  canvas.SaveAs((dir/("asym_"+which+"_"+ch+".png")).string().c_str());
}
inline void DrawProfile(const std::filesystem::path&dir,const std::string&ch,
                        TProfile* h,const std::string&which){
  if(h->GetEntries()==0)return;
  TCanvas canvas(("c_asym_profile_"+ch+"_"+which).c_str(),"",970,700);
  canvas.SetLeftMargin(.14);canvas.SetBottomMargin(.14);
  gStyle->SetOptStat(0);
  h->SetLineColor(kBlue+1);h->SetMarkerColor(kBlue+1);
  h->SetMarkerStyle(20);h->SetMinimum(-1.05);h->SetMaximum(1.05);
  h->Draw("E1");
  TLatex label;label.SetNDC();label.SetTextFont(42);
  label.SetTextSize(.036);label.DrawLatex(.16,.92,"NA60+ Simulation");
  label.SetTextSize(.023);
  label.DrawLatex(.16,.875,"Ion-birth asymmetry, one entry per physical muon-anode (correlated)");
  canvas.SaveAs((dir/("asym_"+which+"_"+ch+".png")).string().c_str());
}
inline void Draw(const std::filesystem::path&out,const std::map<std::string,std::unique_ptr<Plots>>&byScope){
  // Separate files are intentionally made only for pooled/physics-channel
  // populations. All 18 station-conditioned source histograms are in ROOT.
  for(const std::string&ch:{"all","jpsi","phi","omega"}){
    Plots& p=*byScope.at(ch);
    const std::string desc="1800 V P3 noFW, exact Garfield ion births + electron endpoints; "+ch;
    DrawOne(out,ch,p.wireGainA,"gain_vs_A_hemi_wire",desc);
    DrawOne(out,ch,p.seedWireGainA,"gain_vs_A_hemi_primary_wire",desc);
    DrawOne(out,ch,p.muonGainA,"gain_vs_A_hemi_muon",desc);
    DrawOne(out,ch,p.muonShareA,"muon_sharing_vs_A_hemi",desc);
    DrawOne(out,ch,p.tbGainA,"gain_vs_A_TB_wire",desc);
    DrawOne(out,ch,p.ionYieldA,"ion_births_vs_A_hemi_wire",desc);
    DrawOne(out,ch,p.wireFracA,"charge_fraction_vs_A_hemi_wire",desc);
    DrawOne(out,ch,p.aHemiTB,"A_hemi_vs_A_TB",desc);
    DrawOne(out,ch,p.aByLocal,"A_hemi_vs_local_wire",desc);
    DrawOne(out,ch,p.aByImpact,"A_hemi_vs_signed_impact",desc);
    DrawOne(out,ch,p.aByThetaU,"A_hemi_vs_theta_u",desc);
    DrawOne(out,ch,p.aByThetaW,"A_hemi_vs_theta_w",desc);
    DrawOneD(out,ch,p.zeroQIonA,"zero_collected_ion_wires_A_hemi",
             "N_{ion birth,wire}>0, Q_{collected,wire}=0; "+ch);
    DrawOneD(out,ch,p.aHemi,"ion_A_hemi_distribution",
             "Ion-bearing physical anodes, unweighted per anode; "+ch);
    DrawOneD(out,ch,p.aTB,"ion_A_TB_distribution",
             "Ion-bearing physical anodes, top/bottom sectors; "+ch);
    DrawProfile(out,ch,p.meanAByLocal,"mean_A_vs_local_wire");
    DrawProfile(out,ch,p.meanAByImpact,"mean_A_vs_signed_impact");
    DrawProfile(out,ch,p.meanAByThetaU,"mean_A_vs_signed_theta_u");
    DrawProfile(out,ch,p.meanAByThetaW,"mean_A_vs_signed_theta_w");
  }
  // Every channel x station also receives its own wire-resolved physics
  // diagnostic; the full distributions/profiles remain in the ROOT file.
  for(const std::string&ch:{"jpsi","phi","omega"})for(int station=0;station<6;++station){
    const auto key=ch+"_station"+std::to_string(station);
    Plots& p=*byScope.at(key);
    const auto desc="1800 V P3 noFW; "+ch+" station "+std::to_string(station);
    DrawOne(out,key,p.wireGainA,"gain_vs_A_hemi_wire",desc);
    DrawOne(out,key,p.aByLocal,"A_hemi_vs_local_wire",desc);
  }
}
} // namespace MWPC1800Asym
