// C++/ROOT geometry-stratified 252-event reference selection without running
// Garfield. Inclusion probabilities are exact CONDITIONAL ON input strata.
// Subsequent gain/multiwire stress sample MUST use actual completed 2100 V
// Garfield ROOT/CSV data; geometry-only selection never invents charge.
#include <TFile.h>
#include <TTree.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>
namespace{
void Require(bool ok,const std::string& text){
 if(!ok)throw std::runtime_error("STRATIFIED 1000 INPUT: "+text);
}
struct Row{
 int index=0,station=0,det=0,sign=0,boundary=0;
 double du=0,thetaU=0,thetaW=0,w=0,length=0,momentum=0;
 std::string channel;
};
struct Stratum{
 std::vector<Row> events;
 int sample=0;
 double ideal=0;
};
}
int main(int argc,char**argv) {
 try{
  Require(argc==3,"usage select_detailed_1000_root audit.root output/selection.root");
  std::unique_ptr<TFile> in(TFile::Open(argv[1],"READ"));
  Require(in&&!in->IsZombie(),"unreadable source audit.root");
  auto* t=dynamic_cast<TTree*>(in->Get("crossings"));
  Require(t&&t->GetEntries()==1000,"input audit not exactly 1000 frozen Geant4 muons");
  Row row;std::string* channel=nullptr;
  t->SetBranchAddress("index",&row.index);
  t->SetBranchAddress("channel",&channel);
  t->SetBranchAddress("station",&row.station);
  t->SetBranchAddress("det_id",&row.det);
  t->SetBranchAddress("theta_u_deg",&row.thetaU);
  t->SetBranchAddress("theta_w_deg",&row.thetaW);
  t->SetBranchAddress("impact_phase_mm",&row.du);
  t->SetBranchAddress("w_in_mm",&row.w);
  t->SetBranchAddress("path_length_mm",&row.length);
  t->SetBranchAddress("momentum_GeV",&row.momentum);
  using Key=std::tuple<std::string,int,int,int>;
  std::map<Key,Stratum> groups;
  std::vector<Row> all;
  for(Long64_t i=0;i<t->GetEntries();++i){
    t->GetEntry(i);row.channel=*channel;
    Require(row.index==i&&row.station>=0&&row.station<=5,
            "incorrect original frozen Geant4 index or station");
    // Zero signed angle is kept in non-negative stratum; original signed
    // angle remains unchanged in event metadata, no perturbation.
    row.sign=row.thetaU<0?-1:+1;
    row.boundary=std::abs(row.du)>=1.0?1:0;
    groups[{row.channel,row.station,row.sign,row.boundary}].events.push_back(row);
    all.push_back(row);
  }
  const int goal=252;
  int allocated=0;
  for(auto& [key,g]:groups){
    const int N=g.events.size();
    g.ideal=double(goal)*N/all.size();
    // Cover every existing station/channel/sign/phase stratum at least once.
    g.sample=std::min(N,std::max(1,int(std::floor(g.ideal))));
    allocated+=g.sample;
  }
  Require(allocated<=goal,"more occupied strata than geometry reference budget");
  while(allocated<goal) {
    auto best=groups.end();
    double benefit=-1e100;
    for(auto it=groups.begin();it!=groups.end();++it){
      auto& g=it->second;
      if(g.sample>=(int)g.events.size())continue;
      const double delta=g.ideal-g.sample;
      if(delta>benefit){best=it;benefit=delta;}
    }
    Require(best!=groups.end(),"can't allocate requested 252 events");
    ++best->second.sample;++allocated;
  }
  std::filesystem::path path=argv[2];
  std::filesystem::create_directories(path.parent_path());
  std::unique_ptr<TFile> out(TFile::Open(argv[2],"RECREATE","NA60+ input selection"));
  Require(out&&!out->IsZombie(),"cannot write selection ROOT");
  TTree chosen("weighted_selection",
               "Stratified random SRSWOR 252 Geant4 original indices; no gain bias");
  int index=0,station=0,det=0,sign=0,boundary=0,nCell=0,nDraw=0;
  double pi=0,weight=0,thetaU=0,thetaW=0,du=0;
  std::string cname;
  chosen.Branch("index",&index);
  chosen.Branch("channel",&cname);
  chosen.Branch("station",&station);chosen.Branch("det_id",&det);
  chosen.Branch("theta_u_sign",&sign);
  chosen.Branch("anode_boundary",&boundary);
  chosen.Branch("cell_population",&nCell);chosen.Branch("cell_draws",&nDraw);
  chosen.Branch("inclusion_probability",&pi);chosen.Branch("inverse_inclusion_weight",&weight);
  chosen.Branch("theta_u_deg",&thetaU);chosen.Branch("theta_w_deg",&thetaW);
  chosen.Branch("impact_phase_mm",&du);
  std::ofstream csv(path.parent_path()/"weighted_252_original_indices.csv");
  Require(bool(csv),"cannot write weighted selection indices");
  csv<<"index,channel,station,det_id,theta_u_deg,theta_w_deg,impact_phase_mm,"
        "signed_theta_group,anode_boundary,N_stratum,n_stratum,selection_pi,weight\n";
  std::mt19937_64 random(0xA160C21002026ULL);
  std::map<int,int> picked;
  for(auto& [key,g]:groups){
    std::shuffle(g.events.begin(),g.events.end(),random);
    nCell=g.events.size();nDraw=g.sample;
    pi=double(nDraw)/nCell;weight=1./pi;
    for(int j=0;j<nDraw;++j){
      const auto& r=g.events[j];
      index=r.index;cname=r.channel;station=r.station;det=r.det;
      sign=r.sign;boundary=r.boundary;thetaU=r.thetaU;thetaW=r.thetaW;du=r.du;
      Require(picked.emplace(index,1).second,"sampled duplicate frozen index");
      chosen.Fill();
      csv<<index<<","<<cname<<","<<station<<","<<det<<","
         <<thetaU<<","<<thetaW<<","<<du<<","<<sign<<","<<boundary<<","
         <<nCell<<","<<nDraw<<","<<pi<<","<<weight<<"\n";
    }
  }
  Require(chosen.GetEntries()==goal,"missing 252 stratified weighted inputs");
  // Geometry stress sample is deterministic/non-representative; DO NOT use
  // these to fit population distributions or pretend known sampling pi.
  std::sort(all.begin(),all.end(),[](const Row& a,const Row& b){
    const double aScore=3*std::abs(a.du)/1.5+
       .04*(std::abs(a.thetaU)+std::abs(a.thetaW))+
       .002*std::abs(a.w);
    const double bScore=3*std::abs(b.du)/1.5+
       .04*(std::abs(b.thetaU)+std::abs(b.thetaW))+
       .002*std::abs(b.w);
    return aScore==bScore?a.index<b.index:aScore>bScore;
  });
  std::ofstream stress(path.parent_path()/"geometry_stress_24_indices.csv");
  Require(bool(stress),"cannot write 24-case stress test indices");
  stress<<"index,channel,station,theta_u_deg,theta_w_deg,impact_phase_mm,"
           "label,probability\n";
  for(int i=0;i<24;++i){
    const auto& r=all[i];
    stress<<r.index<<","<<r.channel<<","<<r.station<<","<<r.thetaU<<","
          <<r.thetaW<<","<<r.du<<",DELIBERATELY_OVERSAMPLED_GEOMETRY,NA\n";
  }
  out->cd();chosen.Write();out->Write();out->Close();
  std::cout<<"PASS ROOT/C++ weighted stratified Geant4 selection n=252 "
           <<"source=1000, conditional within-stratum pi=n_h/N_h; "
           <<"separate geometry stress=24 (non-representative)\n";
  std::cout<<"HEED multiplicities, charge sharing, and high-gain tails "
           <<"NOT known until actual uncapped Garfield muons are computed\n";
  return 0;
 }catch(const std::exception& e){
  std::cerr<<e.what()<<"\n";return 2;
 }
}
