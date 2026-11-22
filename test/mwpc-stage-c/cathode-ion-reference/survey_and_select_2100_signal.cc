// NA60+ 2100V reference: SELECT from existing K256 real-Geant4 ROOT archives.
// 20 Jpsi development + 10 Phi and 10 Omega SEALED validation. No new Garfield.
// All 997 verified records, not 40 purposive events, are statistical reference.
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
namespace fs=std::filesystem;
namespace {
const std::set<int> previouslyInspected={301,408,446,435,511,101};
void Need(bool pass,const std::string&why){
  if(!pass)throw std::runtime_error("2100V archive selection: "+why);
}
std::vector<std::string> Split(std::string s){
  if(!s.empty()&&s.back()=='\r')s.pop_back();
  std::stringstream ss(s);std::string x;std::vector<std::string> out;
  while(std::getline(ss,x,','))out.push_back(x);
  return out;
}
std::map<std::string,size_t> Header(const std::string& s) {
  std::map<std::string,size_t> out;const auto v=Split(s);
  for(size_t i=0;i<v.size();++i)Need(out.emplace(v[i],i).second,
                                  "duplicate CSV header "+v[i]);
  return out;
}
std::string Cell(const std::vector<std::string>&row,
                 const std::map<std::string,size_t>&m,const std::string&key) {
  const auto it=m.find(key);
  Need(it!=m.end()&&it->second<row.size(),"missing column "+key);
  return row[it->second];
}
std::string Channel(std::string s){
  std::transform(s.begin(),s.end(),s.begin(),
       [](unsigned char c){return std::tolower(c);});
  if(s.find("jpsi")!=std::string::npos)return "jpsi";
  if(s.find("omega")!=std::string::npos)return "omega";
  if(s.find("phi")!=std::string::npos)return "phi";
  throw std::runtime_error("unknown Geant4 source channel: "+s);
}
struct Event{
  int index=-1,st=-1,det=-1,heedSeed=0,nPrim=0;
  int nWires=0,phaseBin=0,signU=0,signW=0;
  bool share5=false,tail=false,wExtreme=false;
  std::string channel,id;
  double u=0,w=0,thetaU=0,thetaW=0,phase=0,path=0,p=0;
  double gain=0,secondFrac=0,ions=0,charge=0;
};
Event ReadOne(const fs::path&dir,int idx){
  std::unique_ptr<TFile> file(TFile::Open((dir/"ion_births.root").string().c_str(),"READ"));
  Need(file&&!file->IsZombie(),"unreadable ROOT "+dir.string());
  auto* m=dynamic_cast<TTree*>(file->Get("muons"));
  auto* p=dynamic_cast<TTree*>(file->Get("primaries"));
  Need(m&&p&&m->GetEntries()==1,"bad ROOT muon/primaries "+dir.string());
  Event e;e.index=idx;
  std::string *source=nullptr,*identity=nullptr;
  m->SetBranchAddress("source_channel",&source);
  m->SetBranchAddress("crossing_id",&identity);
  m->SetBranchAddress("station",&e.st);
  m->SetBranchAddress("det_id",&e.det);
  m->SetBranchAddress("heed_seed",&e.heedSeed);
  m->SetBranchAddress("u_in_mm",&e.u);
  m->SetBranchAddress("w_in_mm",&e.w);
  m->SetBranchAddress("theta_u_deg",&e.thetaU);
  m->SetBranchAddress("theta_w_deg",&e.thetaW);
  m->SetBranchAddress("path_length_mm",&e.path);
  m->SetBranchAddress("momentum_GeV",&e.p);
  m->GetEntry(0);
  Need(source&&identity&&!identity->empty()&&e.st>=0&&e.st<6,
       "bad source/station/identity for slot "+std::to_string(idx));
  e.channel=Channel(*source);e.id=*identity;
  const int nAnodes=int(std::floor((e.st==0?307.2:512.)/3.+.5));
  const double phaseMm=nAnodes%2==0?1.5:0.;
  e.phase=e.u-(std::round((e.u-phaseMm)/3.)*3.+phaseMm);
  Need(std::abs(e.phase)<=1.50000001,"bad 3mm physical wire phase");
  const double d=std::abs(e.phase);
  e.phaseBin=d<.5?0:d<1.?1:2;
  e.signU=(e.thetaU>0)-(e.thetaU<0);
  e.signW=(e.thetaW>0)-(e.thetaW<0);
  e.nPrim=int(p->GetEntries());
  long long ion=0,collected=0;
  p->SetBranchAddress("ion_births_exact",&ion);
  p->SetBranchAddress("collected_electrons",&collected);
  for(Long64_t k=0;k<p->GetEntries();++k){
    p->GetEntry(k);Need(ion>=0&&collected>=0,"bad primary census");
    e.ions+=ion;e.charge+=collected;
  }
  e.gain=e.nPrim?e.charge/e.nPrim:0.;
  std::ifstream wire(dir/"microscopic_wire_muon.csv");
  Need(bool(wire),"no physical-wire per-muon truth "+dir.string());
  std::string line;Need(bool(std::getline(wire,line)),"no wire CSV header");
  const auto hdr=Header(line);
  long long summedIon=0,summedCharge=0;
  std::vector<long long> positives;
  while(std::getline(wire,line)){
    if(line.empty())continue;const auto row=Split(line);
    Need(std::stoi(Cell(row,hdr,"primary_index"))==-1,
         "incorrect wire_muon aggregation level");
    const auto ce=std::stoll(Cell(row,hdr,"collected_e"));
    const auto ib=std::stoll(Cell(row,hdr,"ion_births"));
    Need(ce>=0&&ib>=0,"negative per-wire census");
    summedCharge+=ce;summedIon+=ib;
    if(ce>0)positives.push_back(ce);
  }
  Need(summedCharge==static_cast<long long>(e.charge) &&
       summedIon==static_cast<long long>(e.ions),
       "ROOT/physical-wire charge or ion counts disagree at "+std::to_string(idx));
  std::sort(positives.begin(),positives.end(),std::greater<long long>());
  e.nWires=int(positives.size());
  e.secondFrac=positives.size()>1&&summedCharge>0?
    double(positives[1])/summedCharge:0.;
  e.share5=e.secondFrac>.05;
  return e;
}
uint64_t Tie(int index){
  uint64_t x=uint64_t(index)^0x9e3779b97f4a7c15ULL;
  x^=x>>30;x*=0xbf58476d1ce4e5b9ULL;
  x^=x>>27;x*=0x94d049bb133111ebULL;
  return x^(x>>31);
}
// Purposive phase-space diversity, NOT probability sampling.
// The complete archival pool alone is used for population statistics.
std::vector<const Event*> Choose(const std::vector<Event>&pool,
                        const std::string&channel,int count){
  std::vector<const Event*> options,chosen;
  for(const auto& e:pool)if(e.channel==channel &&
     !previouslyInspected.count(e.index))options.push_back(&e);
  Need((int)options.size()>=count,"insufficient "+channel+" archives");
  std::map<std::string,int> usedBins;std::set<int> usedIdx;
  while((int)chosen.size()<count){
    const Event*best=nullptr;double bestScore=-1e100;
    for(auto* e:options){
      if(usedIdx.count(e->index))continue;
      const auto novelty=[&](const std::string&name,double w){
        return w/(1.+usedBins[name]);
      };
      double s=0;
      s+=novelty("st"+std::to_string(e->st),35.);
      s+=novelty("det"+std::to_string(e->st)+":"+std::to_string(e->det),6.);
      s+=novelty("phase"+std::to_string(e->phaseBin),20.);
      s+=novelty("share"+std::to_string(e->share5),18.);
      s+=novelty("signU"+std::to_string(e->signU),10.);
      s+=novelty("signW"+std::to_string(e->signW),10.);
      s+=novelty("tail"+std::to_string(e->tail),12.);
      s+=novelty("wExtreme"+std::to_string(e->wExtreme),7.);
      s+=novelty("thetaHigh"+std::to_string(std::abs(e->thetaU)>15),6.);
      s+=novelty("pathHigh"+std::to_string(e->path>6.),5.);
      s+=novelty("zero"+std::to_string(e->charge==0),4.);
      s+=double(Tie(e->index)%1000)*1e-9;
      if(s>bestScore){best=e;bestScore=s;}
    }
    Need(best!=nullptr,"could not complete phase-space selection");
    chosen.push_back(best);usedIdx.insert(best->index);
    auto inc=[&](const std::string&name){usedBins[name]++;};
    inc("st"+std::to_string(best->st));
    inc("det"+std::to_string(best->st)+":"+std::to_string(best->det));
    inc("phase"+std::to_string(best->phaseBin));
    inc("share"+std::to_string(best->share5));
    inc("signU"+std::to_string(best->signU));
    inc("signW"+std::to_string(best->signW));
    inc("tail"+std::to_string(best->tail));
    inc("wExtreme"+std::to_string(best->wExtreme));
    inc("thetaHigh"+std::to_string(std::abs(best->thetaU)>15));
    inc("pathHigh"+std::to_string(best->path>6.));
    inc("zero"+std::to_string(best->charge==0));
  }
  return chosen;
}
}
int main(int argc,char**argv){
 try{
  Need(argc==3,"usage survey_and_select_2100_signal VERIFIED_K256_ARCHIVE OUTPUT_DIR");
  const fs::path archiveRoot=fs::absolute(argv[1]),output=fs::absolute(argv[2]);
  Need(fs::is_directory(archiveRoot),"missing 2100V K256 ROOT archive");
  Need(!fs::exists(output),"output exists: choose new frozen output path");
  fs::create_directories(output);
  std::vector<Event> pool;std::vector<int> excluded;
  for(int index=0;index<1000;++index){
    std::ostringstream dir;dir<<"muon_"<<std::setw(4)<<std::setfill('0')<<index;
    const fs::path p=archiveRoot/dir.str();
    if(!fs::is_regular_file(p/"complete.txt")||
       !fs::is_regular_file(p/"SHA256SUMS")||
       !fs::is_regular_file(p/"root_audit.log")||
       !fs::is_regular_file(p/"ion_births.root")){
      excluded.push_back(index);continue;
    }
    pool.push_back(ReadOne(p,index));
  }
  Need(pool.size()>=900,"fewer than 900 sealed-file candidates; audit source archive");
  std::map<std::string,int> counts;
  std::map<std::string,std::vector<double>> gains,absW;
  for(const auto& e:pool){
    counts[e.channel]++;gains[e.channel].push_back(e.gain);
    absW[e.channel].push_back(std::abs(e.w));
  }
  for(auto& [ch,v]:gains)std::sort(v.begin(),v.end());
  for(auto& [ch,v]:absW)std::sort(v.begin(),v.end());
  for(auto& e:pool){
    const auto&g=gains.at(e.channel);const auto&w=absW.at(e.channel);
    e.tail=e.gain>=g[std::min(g.size()-1,size_t(.9*g.size()))];
    e.wExtreme=std::abs(e.w)>=w[std::min(w.size()-1,size_t(.8*w.size()))];
  }
  Need(counts["jpsi"]>=326&&counts["phi"]>=325&&counts["omega"]>=325,
       "too many missing source-channel files");
  const auto dev=Choose(pool,"jpsi",20);
  const auto holdPhi=Choose(pool,"phi",10);
  const auto holdOmega=Choose(pool,"omega",10);
  std::set<int> unique;
  for(auto group:{dev,holdPhi,holdOmega})
    for(const Event*e:group)Need(unique.insert(e->index).second,
                               "duplicate dev/validation index");
  Need(unique.size()==40,"not exactly 40 distinct reference muons");
  std::ofstream csv(output/"selected40.csv");
  Need(bool(csv),"cannot write selection CSV");
  csv<<"role,index,channel,station,det_id,theta_u_deg,theta_w_deg,"
       "u_in_mm,w_in_mm,anode_impact_phase_mm,phase_class,"
       "path_length_mm,momentum_GeV,heed_primaries,ion_births_exact,"
       "collected_electrons,mean_gain_per_primary,positive_collecting_wires,"
       "second_wire_fraction,multi_wire_gt5pct,gain_top10pct,w_extreme_proxy,"
       "crossing_id,selection_design\n";
  std::unique_ptr<TFile> root(TFile::Open((output/"survey.root").string().c_str(),"RECREATE"));
  Need(root&&!root->IsZombie(),"cannot write survey.root");
  int index=0,station=0,det=0,shared=0,phaseClass=0;
  double phase=0,thetaU=0,thetaW=0,gain=0,second=0;
  std::string channel,role;
  TTree selected("selected40","20 Jpsi development and 20 sealed Phi/Omega validation: purposefully stratified");
  selected.SetDirectory(nullptr);
  selected.Branch("index",&index);selected.Branch("channel",&channel);
  selected.Branch("role",&role);selected.Branch("station",&station);
  selected.Branch("det_id",&det);selected.Branch("phase_mm",&phase);
  selected.Branch("theta_u_deg",&thetaU);selected.Branch("theta_w_deg",&thetaW);
  selected.Branch("mean_gain",&gain);selected.Branch("second_wire_fraction",&second);
  selected.Branch("multianode_gt5pct",&shared);selected.Branch("phase_class",&phaseClass);
  for(auto [label,group]:{std::pair<std::string,std::vector<const Event*>>{"development",dev},
                          {"sealed_validation",holdPhi},{"sealed_validation",holdOmega}}){
    for(auto*e:group){
      index=e->index;station=e->st;det=e->det;
      phase=e->phase;thetaU=e->thetaU;thetaW=e->thetaW;gain=e->gain;
      second=e->secondFrac;shared=e->share5;phaseClass=e->phaseBin;
      channel=e->channel;role=label;
      selected.Fill();
      csv<<std::setprecision(15)<<role<<","<<index<<","<<channel<<","<<station<<","<<det
         <<","<<thetaU<<","<<thetaW<<","<<e->u<<","<<e->w<<","<<phase
         <<","<<phaseClass<<","<<e->path<<","<<e->p<<","<<e->nPrim<<","
         <<e->ions<<","<<e->charge<<","<<gain<<","<<e->nWires<<","
         <<second<<","<<int(e->share5)<<","<<int(e->tail)<<","
         <<int(e->wExtreme)<<","<<e->id
         <<",PURPOSIVE_MODEL_DIAGNOSTICS_NOT_POPULATION_WEIGHTED\n";
    }
  }
  Need(selected.GetEntries()==40,"corrupt selected 40 ROOT tree");
  std::ofstream inventory(output/"archive_inventory.csv");
  Need(bool(inventory),"cannot write 997 archive inventory");
  inventory<<"index,channel,station,det_id,theta_u_deg,theta_w_deg,"
            "anode_phase_mm,path_mm,momentum_GeV,heed_primaries,"
            "ion_births_exact,collected_e,mean_gain,n_positive_wires,"
            "second_wire_fraction,multiwire_gt5pct\n";
  TH1D hPhase("h_phase",";Impact relative to nearest physical anode [mm];Verified muons",60,-1.5,1.5);
  TH1D hSecond("h_second",";Second-wire fraction of collected charge;Verified muons",100,0,1);
  TH1D hTheta("h_theta",";Signed #theta_{u} [deg];Verified muons",100,-90,90);
  TH1D hGain("h_gain",";Collected electrons / Heed primary;Verified muons",100,0,200000);
  hPhase.SetDirectory(nullptr);hSecond.SetDirectory(nullptr);
  hTheta.SetDirectory(nullptr);hGain.SetDirectory(nullptr);
  for(const auto&e:pool){
    inventory<<std::setprecision(15)<<e.index<<","<<e.channel<<","<<e.st<<","<<e.det
             <<","<<e.thetaU<<","<<e.thetaW<<","<<e.phase<<","<<e.path
             <<","<<e.p<<","<<e.nPrim<<","<<e.ions<<","<<e.charge<<","
             <<e.gain<<","<<e.nWires<<","<<e.secondFrac<<","<<int(e.share5)<<"\n";
    hPhase.Fill(e.phase);hSecond.Fill(e.secondFrac);
    hTheta.Fill(e.thetaU);hGain.Fill(e.gain);
  }
  root->cd();selected.Write();hPhase.Write();hSecond.Write();
  hTheta.Write();hGain.Write();root->Write();root->Close();
  std::ofstream missing(output/"excluded_slots.txt");
  for(int i:excluded)missing<<std::setfill('0')<<std::setw(4)<<i<<"\n";
  std::ofstream previous(output/"previously_inspected_not_sealed.txt");
  for(int i:previouslyInspected)previous<<std::setw(4)<<std::setfill('0')<<i<<"\n";
  std::ofstream scope(output/"selection_scope.txt");
  scope<<"NA60+ Simulation: 2100 V, P3 noFW, K256 source; archived original Geant4 muons\n"
       <<"archive_root="<<archiveRoot<<"\n"
       <<"candidates_with_seal_files="<<pool.size()<<" of 1000 (this scan does NOT rehash all files)\n"
       <<"channel_counts=jpsi:"<<counts["jpsi"]<<",phi:"<<counts["phi"]
       <<",omega:"<<counts["omega"]<<"\n"
       <<"selection=20 Jpsi development; 10 Phi + 10 Omega sealed validation\n"
       <<"previously_seen_indices_excluded=0301,0408,0446,0435,0511,0101\n"
       <<"population_statistics=all available verified muons, NOT selected40\n"
       <<"sample=purposive coverage of station,det,impact phase,signed angles,"
         "multiwire>5pct,top10pct gain,large abs(w),gas path\n"
       <<"edge_w=UNKNOWN, large abs(w) is ONLY a proxy\n"
       <<"ion_birth_archive=weighted K256, not full ion positions\n"
       <<"simulation_jobs_submitted=0\n";
  std::cout<<"PASS survey ROOT/C++ source pool "<<pool.size()<<"/1000; "
           <<"Jpsi="<<counts["jpsi"]<<" Phi="<<counts["phi"]
           <<" Omega="<<counts["omega"]<<"\n"
           <<"PASS 20 Jpsi development / 10 Phi + 10 Omega sealed validation "
           <<"from existing 2100V K256 archive (NO new Garfield jobs)\n"
           <<"OUTPUT "<<output<<"\n";
  return 0;
 }catch(const std::exception&e){
   std::cerr<<"SURVEY FAILED: "<<e.what()<<"\n";return 2;
 }
}
