// NA60+ Simulation: per-wire collected electron and ion birth distributions
// from actual verified P3-noFW 1800 V 3x~2700 original Geant4 muons.
// Physics/statistics done exclusively in C++/ROOT. No Garfield rerun.
// Input manifest is prepared by a separate SHA256-only file-integrity script.
// Output ROOT TTree+histograms, channel/station summary and ATLAS-like plots.
// No strip induced-current or detector efficiency claims are made here.
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TMath.h>
#include <algorithm>
#include <cctype>
#include <array>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>
namespace fs=std::filesystem;
namespace {
void Need(bool ok,const std::string& why) {
  if(!ok)throw std::runtime_error("1800 V WIRETRUTH AUDIT: "+why);
}
std::vector<std::string> Split(std::string line,char delim=',') {
  if(!line.empty()&&line.back()=='\r')line.pop_back();
  std::vector<std::string> out;std::stringstream ss(line);std::string token;
  while(std::getline(ss,token,delim))out.push_back(token);
  return out;
}
using Header=std::map<std::string,std::size_t>;
Header Columns(const std::string& header) {
  Header h;auto tokens=Split(header);
  for(std::size_t i=0;i<tokens.size();++i)
    Need(h.emplace(tokens[i],i).second,"duplicate header "+tokens[i]);
  return h;
}
std::string Get(const std::vector<std::string>&r,const Header&h,
                const std::string&key) {
  const auto i=h.find(key);
  Need(i!=h.end()&&i->second<r.size(),"missing required field "+key);
  return r[i->second];
}
double Val(const std::vector<std::string>&r,const Header&h,const std::string&key) {
  const auto s=Get(r,h,key);Need(!s.empty(),"missing numeric "+key);
  const double v=std::stod(s);Need(std::isfinite(v),"nonfinite "+key);return v;
}
long long Count(const std::vector<std::string>&r,const Header&h,const std::string&key) {
  const auto s=Get(r,h,key);Need(!s.empty(),"missing count "+key);
  const auto n=std::stoll(s);Need(n>=0,"negative count "+key);return n;
}
std::string Lower(std::string s) {
  std::transform(s.begin(),s.end(),s.begin(),
                 [](unsigned char x){return std::tolower(x);});
  return s;
}
std::string Canon(const std::string&src) {
  const auto s=Lower(src);
  if(s.find("jpsi")!=std::string::npos)return "jpsi";
  if(s.find("omega")!=std::string::npos)return "omega";
  if(s.find("phi")!=std::string::npos)return "phi";
  throw std::runtime_error("Unknown source channel "+src);
}
std::string Identity(const std::vector<std::string>&r,const Header&h) {
  const std::array<std::string,7> f={
    "source","chunk","event","hit_index","track_id","station","det_id"
  };
  std::string id=Get(r,h,"source");
  // Garfield's crossing ID formats the six numeric fields using
  // std::to_string(integer). The original selection may have zero padding.
  for(std::size_t i=1;i<f.size();++i){
    const auto raw=Get(r,h,f[i]);
    std::size_t end=0;
    const auto value=std::stoll(raw,&end);
    Need(end==raw.size(),"invalid original integer crossing ID");
    id+="|"+std::to_string(value);
  }
  return id;
}
struct Muon {
  std::string id,channel;
  int station=0,det=0,pdg=0,base=0,heed=0;
  double u=0,w=0,thetaU=0,thetaW=0,path=0,momentum=0,impact=0;
  long long q[9]={},ions[9]={},up[9]={},down[9]={},boundary[9]={};
  long long nExpected=-1,qExpected=-1,nHeed=-1,capped=0,failed=0;
  bool crossingSeen=false;
  std::set<int> seenWire;
};
struct Seed {
  long long q=0,n=0;int nWireQ=0;
  int garfield=0;
};
using SeedKey=std::pair<std::string,int>;
using WireKey=std::tuple<std::string,int,int>;
struct Counts {
  long long nMuon=0,nZero=0,nSingle=0,nDouble=0,nMany=0,nShare5=0,
            nWireQPositive=0,nWireIonOnly=0,nHeed=0,nCharge=0,nIons=0;
  double sumQ2=0;
};
struct Hist {
  TH1D *hNWire,*hLogQ,*hLogIon,*hLogQWire,*hLogIonWire,*hSecond,*hNear,*hA,*hSeedLogQ,*hOffset;
  TH2D *hPhaseSecond,*hThetaSecond;
  TProfile *pMultiPhase,*pMultiTheta;
  TProfile2D *pQShare,*pIonShare;
  explicit Hist(const std::string&tag) {
    auto name=[&](const char*s){return tag+"_"+s;};
    hNWire=new TH1D(name("n_wire").c_str(),";Anodes with collected electrons per muon;Muons",11,-.5,10.5);
    hLogQ=new TH1D(name("log_muon_Q").c_str(),";log_{10}(collected electrons / muon);Muons with Q > 0",95,0,9.5);
    hLogIon=new TH1D(name("log_muon_ions").c_str(),";log_{10}(ion births / muon);Muons with N > 0",95,0,9.5);
    hLogQWire=new TH1D(name("log_wire_Q").c_str(),";log_{10}(collected electrons / wire);Positive-charge wire rows",95,0,9.5);
    hLogIonWire=new TH1D(name("log_wire_ions").c_str(),";log_{10}(ion births / wire);Positive-ion wire rows",95,0,9.5);
    hSecond=new TH1D(name("second_wire_fraction").c_str(),";Q_{2}/Q_{total};Muons with Q > 0",50,0,.5);
    hNear=new TH1D(name("nearest_wire_fraction").c_str(),";Q_{nearest Geant4 entry wire}/Q_{total};Muons with Q > 0",55,0,1.1);
    hA=new TH1D(name("birth_hemi").c_str(),";(N_{up}-N_{down})/(N_{up}+N_{down});Wire ion groups",50,-1,1);
    hSeedLogQ=new TH1D(name("log_seedQ").c_str(),";log_{10}(collected electrons / Heed primary);Positive-primary rows",95,0,9.5);
    hOffset=new TH1D(name("wire_offset").c_str(),";Collecting wire index - Geant4 nearest wire;Collected-electron fraction",9,-4.5,4.5);
    hPhaseSecond=new TH2D(name("phase_q2").c_str(),";|u_{entry}-u_{nearest wire}| [mm];Q_{2}/Q_{total}",30,0,1.5,50,0,.5);
    hThetaSecond=new TH2D(name("theta_q2").c_str(),";Signed #theta_{u} [deg];Q_{2}/Q_{total}",72,-90,90,50,0,.5);
    pMultiPhase=new TProfile(name("multi_phase").c_str(),";|u_{entry}-u_{nearest wire}| [mm];P(Q_{2}/Q_{total} > 5%)",15,0,1.5);
    pMultiTheta=new TProfile(name("multi_theta").c_str(),";|#theta_{u}| [deg];P(Q_{2}/Q_{total} > 5%)",30,0,90);
    pQShare=new TProfile2D(name("q_fraction_by_wire").c_str(),
      ";Signed impact relative to nearest wire [mm];Collecting-wire offset;Mean Q_{wire}/Q_{muon}",30,-1.5,1.5,9,-4.5,4.5);
    pIonShare=new TProfile2D(name("ion_fraction_by_wire").c_str(),
      ";Signed impact relative to nearest wire [mm];Ion-birth nearest-wire offset;Mean N_{ion,wire}/N_{ion,muon}",30,-1.5,1.5,9,-4.5,4.5);
    const std::array<TH1*,16> all={hNWire,hLogQ,hLogIon,hLogQWire,
      hLogIonWire,hSecond,hNear,hA,hSeedLogQ,hOffset,
      hPhaseSecond,hThetaSecond,pMultiPhase,pMultiTheta,pQShare,pIonShare};
    for(TH1* h:all) {
      h->Sumw2();
      // Keep histograms alive for ATLAS-like PNG rendering after TFile closes.
      h->SetDirectory(nullptr);
    }
  }
  void Write()const {
    const std::array<TH1*,16> all={hNWire,hLogQ,hLogIon,hLogQWire,
      hLogIonWire,hSecond,hNear,hA,hSeedLogQ,hOffset,
      hPhaseSecond,hThetaSecond,pMultiPhase,pMultiTheta,pQShare,pIonShare};
    for(TH1* h:all)h->Write();
  }
};
struct TRow {
  std::string crossing_id,channel;
  int station=0,det_id=0,pdg=0,primary_index=0,
      global_wire=0,local_wire=0,nearest_wire=0,n_active=0,n_active_primary=0;
  double u_entry_mm=0,w_entry_mm=0,theta_u_deg=0,theta_w_deg=0,
         impact_mm=0,wire_u_global_mm=0,q2_fraction=0,q2_over_q1=0,
         nearest_q_fraction=0,ion_A=0,mean_t_birth_ns=0,mean_t_collected_ns=0;
  Long64_t q=0,ions=0,heed_primaries=0,q_total=0,ion_total=0;
};
void Branches(TTree& t,TRow&r){
  t.Branch("crossing_id",&r.crossing_id);
  t.Branch("channel",&r.channel);
  t.Branch("station",&r.station);t.Branch("det_id",&r.det_id);t.Branch("pdg",&r.pdg);
  t.Branch("primary_index",&r.primary_index);
  t.Branch("global_wire",&r.global_wire);t.Branch("local_wire",&r.local_wire);
  t.Branch("nearest_wire",&r.nearest_wire);
  t.Branch("n_active",&r.n_active);t.Branch("n_active_primary",&r.n_active_primary);
  t.Branch("u_entry_mm",&r.u_entry_mm);t.Branch("w_entry_mm",&r.w_entry_mm);
  t.Branch("theta_u_deg",&r.theta_u_deg);t.Branch("theta_w_deg",&r.theta_w_deg);
  t.Branch("impact_mm",&r.impact_mm);
  t.Branch("wire_u_global_mm",&r.wire_u_global_mm);
  t.Branch("q2_fraction",&r.q2_fraction);t.Branch("q2_over_q1",&r.q2_over_q1);
  t.Branch("nearest_q_fraction",&r.nearest_q_fraction);
  t.Branch("ion_A",&r.ion_A);
  t.Branch("mean_t_birth_ns",&r.mean_t_birth_ns);
  t.Branch("mean_t_collected_ns",&r.mean_t_collected_ns);
  t.Branch("q",&r.q);t.Branch("ions",&r.ions);
  t.Branch("heed_primaries",&r.heed_primaries);
  t.Branch("q_total",&r.q_total);t.Branch("ion_total",&r.ion_total);
}
struct Session {
  std::map<std::string,std::unique_ptr<Hist>> hist;
  std::map<std::string,Counts> counts;
  Session(){
    for(const auto&channel:{"all","jpsi","phi","omega"}){
      hist.emplace(channel,std::make_unique<Hist>(channel));
      counts[channel]=Counts{};
    }
    for(const auto&channel:{"jpsi","phi","omega"})
      for(int st=0;st<6;++st){
        const std::string key=std::string(channel)+"_station"+std::to_string(st);
        hist.emplace(key,std::make_unique<Hist>(key));
        counts[key]=Counts{};
      }
  }
  std::vector<std::string> Keys(const Muon&m) {
    return {"all",m.channel,m.channel+"_station"+std::to_string(m.station)};
  }
};
void Drawing(const fs::path&out,Session&s){
  gROOT->SetBatch(true);
  gStyle->SetOptStat(0);
  const std::array<std::string,3> cs={"jpsi","phi","omega"};
  const std::array<int,3> colors={kBlack,kBlue+1,kRed+1};
  struct Panel{std::string field,title;};
  const std::vector<Panel> panels={
    {"hLogQWire","Collected charge per physical anode (Q > 0)"},
    {"hLogIonWire","Ion births assigned to nearest physical anode (N > 0)"},
    {"hLogQ","Collected charge per real Geant4 muon (Q > 0)"},
    {"hNWire","Collected-electron wire multiplicity"},
    {"hSecond","Second-wire fraction of collected charge"},
    {"hNear","Fraction collected on entrance-nearest wire"},
    {"hA","Per-wire positive-ion hemisphere asymmetry"}
  };
  for(const auto& panel:panels){
    TCanvas canvas(("c_"+panel.field).c_str(),"",960,690);
    canvas.SetLeftMargin(.13);canvas.SetBottomMargin(.13);
    TLegend legend(.70,.63,.89,.84);legend.SetBorderSize(0);legend.SetFillStyle(0);
    bool first=true;
    double globalMax=0;
    auto get=[&](const std::string&ch)->TH1D*{
      Hist&h=*s.hist.at(ch);
      if(panel.field=="hLogQWire")return h.hLogQWire;
      if(panel.field=="hLogIonWire")return h.hLogIonWire;
      if(panel.field=="hLogQ")return h.hLogQ;
      if(panel.field=="hNWire")return h.hNWire;
      if(panel.field=="hSecond")return h.hSecond;
      if(panel.field=="hNear")return h.hNear;
      return h.hA;
    };
    for(const auto& ch:cs)globalMax=std::max(globalMax,get(ch)->GetMaximum());
    for(size_t i=0;i<cs.size();++i){
      TH1D*h=get(cs[i]);
      h->SetLineColor(colors[i]);h->SetLineWidth(2);
      h->SetMaximum(std::max(1.,globalMax*1.25));
      h->Draw(first?"HIST":"HIST SAME");
      first=false;legend.AddEntry(h,cs[i].c_str(),"l");
    }
    legend.Draw();
    TLatex label;label.SetNDC();label.SetTextFont(42);
    label.SetTextSize(.038);label.DrawLatex(.16,.90,"NA60+ Simulation");
    label.SetTextSize(.024);label.DrawLatex(.16,.86,"1800 V, P3 noFW, Geant4 + Heed + Garfield++");
    const fs::path filename=out/(panel.field+".png");
    canvas.SaveAs(filename.string().c_str());
  }
  for(const std::string& kind:{"phase","theta"}){
    TCanvas canvas(("c_share_"+kind).c_str(),"",960,690);
    canvas.SetLeftMargin(.13);canvas.SetBottomMargin(.13);
    TLegend legend(.71,.69,.89,.85);legend.SetBorderSize(0);
    bool first=true;
    for(size_t i=0;i<cs.size();++i){
      TProfile* p=kind=="phase"?
        s.hist.at(cs[i])->pMultiPhase:s.hist.at(cs[i])->pMultiTheta;
      p->SetLineColor(colors[i]);p->SetMarkerColor(colors[i]);
      p->SetMarkerStyle(20+i);p->SetMinimum(0);p->SetMaximum(1);
      p->Draw(first?"E1":"E1 SAME");first=false;
      legend.AddEntry(p,cs[i].c_str(),"lp");
    }
    legend.Draw();
    TLatex label;label.SetNDC();label.SetTextFont(42);
    label.SetTextSize(.038);label.DrawLatex(.16,.90,"NA60+ Simulation");
    label.SetTextSize(.024);label.DrawLatex(.16,.86,"Measured two-wire share > 5% of muon charge");
    canvas.SaveAs((out/("share_probability_vs_"+kind+".png")).string().c_str());
  }
  for(const auto&channel:{"all","jpsi","phi","omega"}){
    for(const std::string&kind:{"collected_electrons","ion_births"}){
      TCanvas canvas(("c_2d_"+std::string(channel)+"_"+kind).c_str(),"",1000,710);
      canvas.SetLeftMargin(.13);canvas.SetRightMargin(.16);canvas.SetBottomMargin(.13);
      TProfile2D* h=kind=="collected_electrons"?
        s.hist.at(channel)->pQShare:s.hist.at(channel)->pIonShare;
      h->SetMinimum(0);h->SetMaximum(1);h->Draw("COLZ");
      TLatex label;label.SetNDC();label.SetTextFont(42);
      label.SetTextSize(.038);label.DrawLatex(.17,.92,"NA60+ Simulation");
      canvas.SaveAs((out/(std::string(channel)+"_per_wire_fraction_"+kind+".png")).string().c_str());
    }
  }
}
}
int main(int argc,char**argv) {
  try {
    Need(argc==3,"usage: wiretruth1800_distributions chunk_manifest.tsv output_directory");
    const fs::path manifest=fs::absolute(argv[1]),out=fs::absolute(argv[2]);
    Need(fs::is_regular_file(manifest),"missing frozen SHA-verified manifest");
    Need(!fs::exists(out),"output directory already exists: refusing overwrite");
    fs::create_directories(out);
    std::unique_ptr<TFile> root(TFile::Open((out/"wiretruth1800_8k.root").string().c_str(),"RECREATE"));
    Need(root&&!root->IsZombie(),"ROOT output cannot be opened");
    TTree muons("muons","Exact collected-anode endpoints per Geant4 crossing");
    TTree wires("wires","Positive electron and/or positive ion per physical anode");
    TTree seeds("seeds","Per-Heed-primary collected electrons, ion birth census and sharing");
    TTree seedWires("seed_wires","Per-Heed-primary exact physical-wire collected electrons and ion birth census");
    // Stack-owned TTrees must not be deleted by TFile::Close().
    // Their entries are explicitly written to the ROOT output below.
    muons.SetDirectory(nullptr);wires.SetDirectory(nullptr);
    seeds.SetDirectory(nullptr);seedWires.SetDirectory(nullptr);
    TRow mrow,wrow,srow,swrow;
    Branches(muons,mrow);Branches(wires,wrow);Branches(seeds,srow);
    Branches(seedWires,swrow);
    Session study;
    std::ifstream mf(manifest);Need(bool(mf),"cannot open snapshot manifest");
    std::string line;std::getline(mf,line);
    Need(line=="channel\tchunk\toriginal_input\twire_muon\twire_seed\tcrossing\tseed",
         "wrong snapshot manifest header");
    std::set<std::pair<std::string,int>> usedChunks;
    std::set<std::string> usedCrossings;
    int nChunks=0;
    while(std::getline(mf,line)){
      if(line.empty())continue;
      const auto cols=Split(line,'\t');Need(cols.size()==7,"malformed manifest record");
      const std::string channel=Lower(cols[0]);
      Need(channel=="jpsi"||channel=="phi"||channel=="omega",
           "invalid frozen physics channel");
      const int index=std::stoi(cols[1]);
      Need(index>=0&&index<500&&usedChunks.emplace(channel,index).second,
           "duplicate/invalid source chunk");
      ++nChunks;
      std::map<std::string,Muon> events;
      {
        std::ifstream in(cols[2]);Need(bool(in),"missing frozen Geant4 input "+cols[2]);
        std::string header;Need(bool(std::getline(in,header)),"empty original input");
        const auto h=Columns(header);
        int n=0;
        while(std::getline(in,line)){
          if(line.empty()||line=="\r")continue;
          const auto row=Split(line);const auto id=Identity(row,h);
          Muon m;m.id=id;m.channel=Canon(Get(row,h,"source"));
          Need(m.channel==channel,"frozen physics channel mismatch");
          m.station=int(Val(row,h,"station"));m.det=int(Val(row,h,"det_id"));
          m.pdg=int(Val(row,h,"pdg"));m.heed=int(Val(row,h,"heed_seed"));
          m.u=Val(row,h,"u_in_mm");m.w=Val(row,h,"w_in_mm");
          m.thetaU=Val(row,h,"theta_u_in_deg");
          m.thetaW=Val(row,h,"theta_w_in_deg");
          m.path=Val(row,h,"path_length_mm");
          m.momentum=Val(row,h,"momentum_GeV");
          Need(m.station>=0&&m.station<=5&&std::abs(m.pdg)==13,
               "invalid original physical muon/station");
          const double phaseMm=m.station==0?1.5:0.;
          m.base=int(std::llround((m.u-phaseMm)/3.));
          m.impact=m.u-(phaseMm+3.*m.base);
          Need(std::abs(m.impact)<=1.500000001,"wrong 3 mm station-dependent phase");
          Need(usedCrossings.insert(id).second,"duplicate independent original crossing");
          Need(events.emplace(id,m).second,"duplicate crossing in single original input");
          ++n;
        }
        Need(n==20,"original source chunk does not contain exactly 20 muons");
      }
      {
        std::ifstream in(cols[5]);Need(bool(in),"missing Garfield crossing totals");
        std::string header;Need(bool(std::getline(in,header)),"empty crossing totals");
        auto h=Columns(header);int n=0;
        while(std::getline(in,line)){
          if(line.empty())continue;
          auto r=Split(line);
          const auto cid=Get(r,h,"source|chunk|event|hit_index|track_id|station|det_id");
          auto it=events.find(cid);Need(it!=events.end(),"unknown crossing in physics output");
          auto&m=it->second;Need(!m.crossingSeen,"duplicate crossing physics row");
          m.qExpected=Count(r,h,"collected_avalanche_electrons");
          m.nExpected=Count(r,h,"ion_births");
          m.nHeed=Count(r,h,"heed_primaries");
          m.capped=Count(r,h,"capped_seeds");
          m.failed=Count(r,h,"failed_seeds");
          Need(m.capped==0&&m.failed==0,"capped/failed avalanche in completed chunk");
          m.crossingSeen=true;++n;
        }
        Need(n==20,"Garfield crossing output missing muons");
      }
      std::map<WireKey,std::array<long long,2>> bySeedWire;
      std::map<std::pair<std::string,int>,std::array<long long,2>> byMuonWire;
      std::map<SeedKey,Seed> originalSeeds;
      {
        std::ifstream in(cols[6]);Need(bool(in),"missing Garfield individual seeds");
        std::string header;Need(bool(std::getline(in,header)),"empty seed CSV");
        auto h=Columns(header);
        while(std::getline(in,line)){
          if(line.empty())continue;auto r=Split(line);
          auto cid=Get(r,h,"crossing_id");Need(events.count(cid),"unknown seed crossing");
          int pi=int(Val(r,h,"primary_index"));
          Need(pi>=0,"invalid primary index");
          Need(Count(r,h,"capped")==0 &&
               Count(r,h,"successful")==1,"failed/capped primary must be excluded");
          Seed s;s.q=Count(r,h,"collected_e");s.n=Count(r,h,"ion_births");
          s.garfield=int(Val(r,h,"garfield_seed"));
          Need(originalSeeds.emplace(SeedKey{cid,pi},s).second,
               "duplicate original Heed primary");
        }
      }
      {
        std::ifstream in(cols[4]);Need(bool(in),"missing per-primary/physical-wire CSV");
        std::string header;Need(bool(std::getline(in,header)),"empty wire_seed CSV");
        auto h=Columns(header);
        while(std::getline(in,line)){
          if(line.empty())continue;auto r=Split(line);
          const auto cid=Get(r,h,"crossing_id");
          auto it=events.find(cid);Need(it!=events.end(),"unknown per-seed wire crossing");
          const auto&m=it->second;
          const int pi=int(Val(r,h,"primary_index")),k=int(Val(r,h,"wire_index_global"));
          const long long q=Count(r,h,"collected_e"),n=Count(r,h,"ion_births");
          Need(pi>=0&&originalSeeds.count({cid,pi}),"wire seed missing original primary");
          Need(std::abs(k-m.base)<=4 && (q>0||n>0),
               "invalid physical wire or empty per-seed wire row");
          const double wireU=Val(r,h,"wire_u_global_mm");
          Need(std::abs(wireU-(m.station==0?1.5:0.)-3.*k)<1e-7,
               "real wire centre / station phase mismatch");
          Need(bySeedWire.emplace(WireKey{cid,pi,k},std::array<long long,2>{q,n}).second,
               "duplicate per-primary physical wire");
          byMuonWire[{cid,k}][0]+=q;
          byMuonWire[{cid,k}][1]+=n;
        }
      }
      std::map<SeedKey,std::array<long long,2>> seedSums;
      for(const auto&[key,v]:bySeedWire){
        const SeedKey p{std::get<0>(key),std::get<1>(key)};
        seedSums[p][0]+=v[0];seedSums[p][1]+=v[1];
        if(v[0])originalSeeds.at(p).nWireQ++;
      }
      for(const auto&[key,s]:originalSeeds){
        const auto v=seedSums[key];
        Need(v[0]==s.q&&v[1]==s.n,"per-primary wire/seed endpoint-ion sum mismatch");
      }
      {
        std::ifstream in(cols[3]);Need(bool(in),"missing per-muon physical wire CSV");
        std::string header;Need(bool(std::getline(in,header)),"empty wire_muon CSV");
        auto h=Columns(header);
        while(std::getline(in,line)){
          if(line.empty())continue;auto r=Split(line);
          const auto cid=Get(r,h,"crossing_id");
          auto it=events.find(cid);Need(it!=events.end(),"unknown wire-muon crossing");
          auto&m=it->second;
          Need(int(Val(r,h,"station"))==m.station,"wire station does not match G4");
          Need(int(Val(r,h,"primary_index"))==-1&&
               int(Val(r,h,"cluster_index"))==-1,
               "wire_muon row is not per-muon aggregation");
          int k=int(Val(r,h,"wire_index_global")),off=k-m.base;
          Need(std::abs(off)<=4,"physical anode outside 9-wire patch");
          Need(m.seenWire.insert(k).second,"duplicate wire within one muon");
          const double physicalU=Val(r,h,"wire_u_global_mm");
          Need(std::abs(physicalU-(m.station==0?1.5:0.)-3.*k)<1e-6,
               "global physical anode index mismatch");
          const auto q=Count(r,h,"collected_e"),n=Count(r,h,"ion_births");
          Need(q>0||n>0,"unexpected empty wire record");
          const long long up=Count(r,h,"ion_up");
          const long long down=Count(r,h,"ion_down");
          const long long boundary=Count(r,h,"ion_boundary");
          Need(up+down+boundary==n,"signed positive-ion birth population inconsistent");
          const std::size_t slot=std::size_t(off+4);
          m.q[slot]=q;m.ions[slot]=n;
          m.up[slot]=up;m.down[slot]=down;m.boundary[slot]=boundary;
          const auto& fromSeeds=byMuonWire.at({cid,k});
          Need(fromSeeds[0]==q&&fromSeeds[1]==n,
               "wire_muon != independently summed Heed-primary wire truth");
        }
      }
      for(auto&[cid,m]:events){
        Need(m.crossingSeen,"missing original Garfield crossing");
        const auto qsum=std::accumulate(std::begin(m.q),std::end(m.q),0LL);
        const auto nsum=std::accumulate(std::begin(m.ions),std::end(m.ions),0LL);
        Need(qsum==m.qExpected&&nsum==m.nExpected,
             "per-muon real physical wire electron/ion census mismatch");
        long long qSeeds=0,nSeeds=0;int nPrimary=0;
        for(const auto&[key,s]:originalSeeds)
          if(key.first==cid){qSeeds+=s.q;nSeeds+=s.n;++nPrimary;}
        Need(qsum==qSeeds&&nsum==nSeeds&&nPrimary==m.nHeed,
             "seed/muon real charge and Heed-primary census mismatch");
        std::vector<std::pair<long long,int>> positive;
        for(int i=0;i<9;++i)if(m.q[i])positive.push_back({m.q[i],i});
        std::sort(positive.begin(),positive.end(),std::greater<>());
        const int nw=int(positive.size());
        const double f2=nw>=2?double(positive[1].first)/qsum:0.;
        const double f21=nw>=2?double(positive[1].first)/positive[0].first:0.;
        const bool sharing=f2>.05;
        const double near=qsum?double(m.q[4])/qsum:std::numeric_limits<double>::quiet_NaN();
        TRow base;base.crossing_id=cid;base.channel=channel;
        base.station=m.station;base.det_id=m.det;base.pdg=m.pdg;
        base.nearest_wire=m.base;base.n_active=nw;
        base.u_entry_mm=m.u;base.w_entry_mm=m.w;
        base.theta_u_deg=m.thetaU;base.theta_w_deg=m.thetaW;
        base.impact_mm=m.impact;base.q_total=qsum;
        base.ion_total=nsum;base.heed_primaries=m.nHeed;
        base.q2_fraction=qsum?f2:std::numeric_limits<double>::quiet_NaN();
        base.q2_over_q1=qsum?f21:std::numeric_limits<double>::quiet_NaN();
        base.nearest_q_fraction=near;
        mrow=base;mrow.q=qsum;mrow.ions=nsum;muons.Fill();
        for(auto&key:study.Keys(m)){
          auto&c=study.counts.at(key);auto&h=*study.hist.at(key);
          ++c.nMuon;c.nCharge+=qsum;c.nIons+=nsum;c.nHeed+=m.nHeed;
          c.nZero+=(nw==0);c.nSingle+=(nw==1);c.nDouble+=(nw==2);c.nMany+=(nw>=3);
          c.nShare5+=sharing;
          h.hNWire->Fill(nw);
          if(qsum){h.hLogQ->Fill(std::log10(double(qsum)));h.hSecond->Fill(f2);
            h.hNear->Fill(near);h.hPhaseSecond->Fill(std::abs(m.impact),f2);
            h.hThetaSecond->Fill(m.thetaU,f2);h.pMultiPhase->Fill(std::abs(m.impact),sharing);
            h.pMultiTheta->Fill(std::abs(m.thetaU),sharing);c.sumQ2+=f2;
          }
          if(nsum)h.hLogIon->Fill(std::log10(double(nsum)));
          // Zero-charge wires MUST participate in the conditional mean:
          // one independent Geant4 muon contributes all nine fractions.
          for(int i=0;i<9;++i){
            if(qsum)h.pQShare->Fill(m.impact,i-4,double(m.q[i])/qsum);
            if(nsum)h.pIonShare->Fill(m.impact,i-4,double(m.ions[i])/nsum);
          }
        }
        for(int i=0;i<9;++i){
          if(!m.q[i]&&!m.ions[i])continue;
          const int k=m.base+i-4;wrow=base;
          wrow.global_wire=k;wrow.local_wire=i-4;
          wrow.wire_u_global_mm=(m.station==0?1.5:0.)+3.*k;
          wrow.q=m.q[i];wrow.ions=m.ions[i];
          wrow.ion_A=(m.up[i]+m.down[i])?
            double(m.up[i]-m.down[i])/(m.up[i]+m.down[i]):
            std::numeric_limits<double>::quiet_NaN();
          wires.Fill();
          for(auto&key:study.Keys(m)){
            auto&c=study.counts.at(key);auto&h=*study.hist.at(key);
            if(wrow.q){++c.nWireQPositive;h.hLogQWire->Fill(std::log10(double(wrow.q)));
              h.hOffset->Fill(i-4,double(wrow.q)/std::max(1LL,qsum));
            }
            if(wrow.ions)h.hLogIonWire->Fill(std::log10(double(wrow.ions)));
            if(wrow.ions&&!wrow.q)++c.nWireIonOnly;
            if(std::isfinite(wrow.ion_A))h.hA->Fill(wrow.ion_A);
          }
        }
        for(const auto&[key,s]:originalSeeds){
          if(key.first!=cid)continue;
          std::vector<std::pair<long long,int>> rankedPrimary;
          long long primaryNearest=0;
          auto it=bySeedWire.lower_bound(WireKey{cid,key.second,std::numeric_limits<int>::min()});
          const auto begin=it;
          for(;it!=bySeedWire.end()&&std::get<0>(it->first)==cid&&
                       std::get<1>(it->first)==key.second;++it){
            const int wire=std::get<2>(it->first);
            const long long q=it->second[0];
            if(q)rankedPrimary.push_back({q,wire});
            if(wire==m.base)primaryNearest=q;
          }
          std::sort(rankedPrimary.begin(),rankedPrimary.end(),std::greater<>());
          const double fp2=rankedPrimary.size()>1?
            double(rankedPrimary[1].first)/s.q:0.;
          const double fp21=rankedPrimary.size()>1?
            double(rankedPrimary[1].first)/rankedPrimary[0].first:0.;
          srow=base;
          srow.primary_index=key.second;
          srow.q=s.q;srow.ions=s.n;
          srow.q_total=s.q;srow.ion_total=s.n;
          srow.n_active=s.nWireQ;srow.n_active_primary=s.nWireQ;
          srow.q2_fraction=s.q?fp2:std::numeric_limits<double>::quiet_NaN();
          srow.q2_over_q1=s.q?fp21:std::numeric_limits<double>::quiet_NaN();
          srow.nearest_q_fraction=s.q?
            double(primaryNearest)/s.q:std::numeric_limits<double>::quiet_NaN();
          seeds.Fill();
          for(auto wIt=begin;wIt!=it;++wIt) {
            const int wire=std::get<2>(wIt->first);
            swrow=srow;
            swrow.global_wire=wire;swrow.local_wire=wire-m.base;
            swrow.wire_u_global_mm=(m.station==0?1.5:0.)+3.*wire;
            swrow.q=wIt->second[0];swrow.ions=wIt->second[1];
            seedWires.Fill();
          }
          if(s.q)
            for(auto&group:study.Keys(m))
              study.hist.at(group)->hSeedLogQ->Fill(std::log10(double(s.q)));
        }
      }
    }
    Need(nChunks>0,"no completed 1800V physics chunks in manifest");
    root->cd();muons.Write();wires.Write();seeds.Write();seedWires.Write();
    for(auto&[key,h]:study.hist)h->Write();
    root->Write();root->Close();
    std::ofstream summary(out/"channel_station_summary.csv");
    Need(bool(summary),"cannot write per-channel/station independent-muon summary");
    summary<<"scope,n_muons,n_zero_charge,n_single_anode,n_two_anode,n_3plus_anodes,"
             "n_5pct_multianode,frac_multianode_gt5pct,binomial_se,"
             "mean_q2_over_qtotal_positives,total_collected_e,total_ion_births,"
             "total_heed_primaries,n_positive_charge_wire_rows,"
             "n_ion_but_zero_e_wire_rows\n";
    for(const auto&[key,c]:study.counts){
      const double f=c.nMuon?double(c.nShare5)/c.nMuon:0.;
      const double se=c.nMuon?std::sqrt(f*(1-f)/c.nMuon):0.;
      const double meanQ2=c.nMuon>c.nZero?
         c.sumQ2/(c.nMuon-c.nZero):0.;
      Need(c.nMuon==c.nZero+c.nSingle+c.nDouble+c.nMany,
           "invalid independent-muon multiplicity count");
      summary<<std::setprecision(16)<<key<<","<<c.nMuon<<","<<c.nZero<<","
             <<c.nSingle<<","<<c.nDouble<<","<<c.nMany<<","<<c.nShare5
             <<","<<f<<","<<se<<","<<meanQ2<<","<<c.nCharge<<","
             <<c.nIons<<","<<c.nHeed<<","<<c.nWireQPositive<<","
             <<c.nWireIonOnly<<"\n";
    }
    summary.close();
    Drawing(out,study);
    std::ofstream status(out/"analysis_readme.txt");
    status<<"NA60+ Simulation: verified 1800 V, P3 noFW, ArCO2 70:30, Penning 0.57\n"
          <<"Frozen original Geant4 crossings, TrackHeed and microscopic Garfield++\n"
          <<"Source chunks="<<nChunks<<" original muons="<<study.counts.at("all").nMuon<<"\n"
          <<"This is a partial 30k production snapshot, NOT a completed 30k gate.\n"
          <<"Collected electron endpoint counts assigned to actual anode.\n"
          <<"Ion births assigned independently to nearest anode at ionisation.\n"
          <<"Zero electrons with nonzero ion births are retained.\n"
          <<"Q2 = second-largest positive collecting-wire count; f2=Q2/sum(all Q).\n"
          <<"All zero-gain muons included in multiplicity and sharing denominator.\n"
          <<"Nearest wire uses original Geant4 entrance, 3mm pitch, station0 +1.5mm.\n"
          <<"Uncertainties in summary are naive independent-muon binomial SE;\n"
          <<"no systematic error or correlated-channel uncertainty estimated.\n"
          <<"Muon wire/seed/parent and complete original counts independently agree.\n"
          <<"Only SHA-frozen completed chunks enter; no cap/failure in accepted events.\n"
          <<"No strip weighting field, ion transport, electronics, or detection threshold.\n"
          <<"WARNING: 1800 V historical code labels P3 but numerical geometry is\n"
          <<"3mm pitch and symmetric +/-3mm cathodes; not the as-built chamber.\n";
    std::cout<<"PASS 1800V WIRETRUTH C++/ROOT real muons="
             <<study.counts.at("all").nMuon
             <<" chunks="<<nChunks<<" Jpsi="<<study.counts.at("jpsi").nMuon
             <<" Phi="<<study.counts.at("phi").nMuon
             <<" Omega="<<study.counts.at("omega").nMuon
             <<" zero_charge="<<study.counts.at("all").nZero
             <<" Q2_fraction_over_5pct="<<study.counts.at("all").nShare5
             <<" output="<<out<<"\n";
    return 0;
  }catch(const std::exception&e){
    std::cerr<<"FAILED: "<<e.what()<<"\n";return 2;
  }
}
