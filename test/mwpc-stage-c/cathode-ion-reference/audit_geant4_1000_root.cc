// C++/ROOT-only audit of the EXACT frozen 1000 true Geant4 crossings.
// Geometry-only audit: Heed counts/charge sharing/gain require actual Garfield.
// Writes station-/channel-conditioned ROOT histograms and ATLAS-like PNGs.
// Usage: audit_geant4_1000_root original1000.csv output/audit.root
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TROOT.h>
#include <array>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
namespace {
std::vector<std::string> Split(const std::string& s) {
  std::stringstream ss(s);std::string w;std::vector<std::string> out;
  while(std::getline(ss,w,',')) {
    if(!w.empty()&&w.back()=='\r')w.pop_back();
    out.push_back(w);
  }
  return out;
}
void Must(bool b,const std::string& why){
  if(!b)throw std::runtime_error("INPUT ROOT AUDIT: "+why);
}
std::string Normalized(const std::string& s) {
  if(s.find("jpsi")!=std::string::npos||s.find("Jpsi")!=std::string::npos)
    return "jpsi";
  if(s.find("phi")!=std::string::npos||s.find("Phi")!=std::string::npos)
    return "phi";
  if(s.find("omega")!=std::string::npos||s.find("Omega")!=std::string::npos)
    return "omega";
  throw std::runtime_error("Unknown source channel: "+s);
}
}
int main(int argc,char**argv){
 try {
  Must(argc==3,"usage audit_geant4_1000_root input.csv output.root");
  std::ifstream in(argv[1]);Must(bool(in),"missing original 1000 selection");
  std::string line;std::getline(in,line);
  const std::string expected=
    "source,chunk,event,hit_index,track_id,station,det_id,pdg,"
    "u_in_mm,v_in_mm,w_in_mm,u_out_mm,v_out_mm,w_out_mm,"
    "theta_u_in_deg,theta_w_in_deg,path_length_mm,edep_eV,"
    "pu_in_GeV,pv_in_GeV,pw_in_GeV,momentum_GeV,heed_seed";
  if(!line.empty()&&line.back()=='\r')line.pop_back();
  Must(line==expected,"Geant4 schema does not match native Garfield source");
  std::filesystem::create_directories(std::filesystem::path(argv[2]).parent_path());
  std::unique_ptr<TFile> file(TFile::Open(argv[2],"RECREATE"));
  Must(file&&!file->IsZombie(),"cannot create ROOT histograms");
  int index=0,station=0,det=0,heedSeed=0,signU=0,signW=0;
  double thetaU=0,thetaW=0,u=0,v=0,w=0,du=0,path=0,momentum=0,
         uToEdge=0, wToEdge=0;
  std::string channel;
  auto* tree=new TTree("crossings","Exact original 1000 Geant4 inputs; no simulated Garfield observables");
  tree->Branch("index",&index);tree->Branch("channel",&channel);
  tree->Branch("station",&station);tree->Branch("det_id",&det);
  tree->Branch("heed_seed",&heedSeed);
  tree->Branch("theta_u_deg",&thetaU);tree->Branch("theta_w_deg",&thetaW);
  tree->Branch("theta_u_sign",&signU);tree->Branch("theta_w_sign",&signW);
  tree->Branch("u_in_mm",&u);tree->Branch("v_in_mm",&v);tree->Branch("w_in_mm",&w);
  tree->Branch("impact_phase_mm",&du);
  tree->Branch("edge_u_distance_mm",&uToEdge);
  tree->Branch("edge_w_distance_mm",&wToEdge);
  tree->Branch("path_length_mm",&path);tree->Branch("momentum_GeV",&momentum);
  std::map<std::string,std::array<int,6>> counts;
  std::map<std::string,TH1D*> hist;
  for(auto ch:{"jpsi","phi","omega"}){
    for(int st=0;st<6;++st) {
      const std::string stem=std::string(ch)+"_st"+std::to_string(st);
      auto add=[&](const std::string& obs,int bins,double lo,double hi){
        auto* h=new TH1D((stem+"_"+obs).c_str(),
            (std::string(ch)+" station "+std::to_string(st)+";"+obs+";Events").c_str(),
            bins,lo,hi);
        h->Sumw2();h->SetDirectory(file.get());
        hist[stem+"_"+obs]=h;
      };
      add("theta_u_deg",100,-90,90);
      add("theta_w_deg",100,-90,90);
      add("impact_phase_mm",60,-1.5,1.5);
      add("w_in_mm",100,-1000,1000);
      add("momentum_GeV",80,0,80);
      add("path_mm",100,0,100);
      add("u_edge_mm",70,-20,350);
      add("w_edge_mm",70,-20,600);
    }
  }
  while(std::getline(in,line)){
    if(line.empty())continue;
    auto f=Split(line);
    Must(f.size()==23,"exact Geant4 row has wrong number of columns");
    channel=Normalized(f[0]);station=std::stoi(f[5]);det=std::stoi(f[6]);
    Must(station>=0&&station<6,"invalid station");
    u=std::stod(f[8]);v=std::stod(f[9]);w=std::stod(f[10]);
    thetaU=std::stod(f[14]);thetaW=std::stod(f[15]);
    path=std::stod(f[16]);momentum=std::stod(f[21]);
    heedSeed=std::stoi(f[22]);
    const double pitch=3.;
    const double width=station==0?307.2:512.;
    const int n=int(std::floor(width/pitch+.5));
    const double phase=n%2==0?1.5:0.;
    du=u-(std::round((u-phase)/pitch)*pitch+phase);
    uToEdge=width*.5-std::abs(u);
    // W half-length is not frozen here. Use an unbounded proxy; DO NOT
    // advertise this as true chamber-edge clearance.
    wToEdge=-1.;
    signU=thetaU>0?1:(thetaU<0?-1:0);
    signW=thetaW>0?1:(thetaW<0?-1:0);
    Must(std::isfinite(du)&&std::abs(du)<=1.50000001,
         "physical wire-centered phase inconsistent");
    const std::string stem=channel+"_st"+std::to_string(station);
    hist.at(stem+"_theta_u_deg")->Fill(thetaU);
    hist.at(stem+"_theta_w_deg")->Fill(thetaW);
    hist.at(stem+"_impact_phase_mm")->Fill(du);
    hist.at(stem+"_w_in_mm")->Fill(w);
    hist.at(stem+"_momentum_GeV")->Fill(momentum);
    hist.at(stem+"_path_mm")->Fill(path);
    hist.at(stem+"_u_edge_mm")->Fill(uToEdge);
    counts[channel][station]++;
    tree->Fill();index++;
  }
  Must(index==1000,"frozen Geant4 selection not exactly 1000 rows");
  Must(std::accumulate(counts["jpsi"].begin(),counts["jpsi"].end(),0)==334 &&
       std::accumulate(counts["phi"].begin(),counts["phi"].end(),0)==333 &&
       std::accumulate(counts["omega"].begin(),counts["omega"].end(),0)==333,
       "frozen physics-channel 334/333/333 counts drifted");
  for(auto ch:{"jpsi","phi","omega"}){
    std::cout<<"SOURCE "<<ch<<" station counts:";
    for(int st=0;st<6;++st)std::cout<<" "<<counts[ch][st];
    std::cout<<"\n";
  }
  file->cd();tree->Write();
  for(const auto& [key,h]:hist)h->Write();
  // ROOT-only ATLAS-like visual diagnostics, with honest simulation label.
  gROOT->SetBatch(true);gStyle->SetOptStat(0);
  for(auto ch:{"jpsi","phi","omega"}){
    TCanvas canvas(("c_"+std::string(ch)).c_str(),"",900,650);
    canvas.SetLeftMargin(.13);canvas.SetBottomMargin(.12);
    TLegend legend(.68,.56,.88,.86);legend.SetBorderSize(0);
    bool first=true;
    for(int st=0;st<6;++st){
      const std::string key=std::string(ch)+"_st"+std::to_string(st)+"_theta_u_deg";
      TH1D* h=hist.at(key);
      h->SetLineColor(1+st);h->SetLineWidth(2);
      h->SetTitle(";Signed #theta_{u} [deg];Geant4 muon crossings");
      h->Draw(first?"HIST":"HIST SAME");first=false;
      legend.AddEntry(h,("Station "+std::to_string(st)).c_str(),"l");
    }
    legend.Draw();
    TLatex label;label.SetNDC(true);label.SetTextSize(.037);
    label.DrawLatex(.17,.89,"NA60+ Simulation");
    label.SetTextSize(.027);
    label.DrawLatex(.17,.84,(std::string(ch)+", unprocessed Geant4 input").c_str());
    const auto png=std::filesystem::path(argv[2]).parent_path()/
      (std::string(ch)+"_theta_u_stations.png");
    canvas.SaveAs(png.string().c_str());
  }
  file->Write();file->Close();
  std::cout<<"PASS ROOT/C++ frozen Geant4 1000 channel/station geometry audit, "
              "NO inferred Heed/gain/multi-wire statistics\n";
  return 0;
 }catch(const std::exception& e){
  std::cerr<<e.what()<<"\n";return 2;
 }
}
