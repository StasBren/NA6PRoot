// NA60+ Stage C: READ-ONLY statistical summary of completed 2100V K256
// genuine-Geant4 / Heed / uncapped Garfield++ muon archives.
// Designed to run while a separate WSL overnight writer continues working.
// Only completed/sealed + ROOT-audited folders contribute to distributions.
// Missing/failed/unsealed selection slots remain in a separate exclusion CSV.
// This is NOT a cathode waveform, signal-kernel, strip, or VMM simulation.
#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TTree.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace fs=std::filesystem;
namespace {
void Check(bool b,const std::string&what) {
 if(!b)throw std::runtime_error("AUDITED MUON SUMMARY: "+what);
}
std::string Arg(int n,char**v,const std::string&key) {
 for(int j=1;j+1<n;++j)if(key==v[j])return v[j+1];
 return "";
}
std::vector<std::string> Split(const std::string&s) {
 std::istringstream ss(s);std::string x;std::vector<std::string> out;
 while(std::getline(ss,x,',')) {
   if(!x.empty()&&x.back()=='\r')x.pop_back();
   out.push_back(x);
 }
 return out;
}
std::map<std::string,int> Cols(const std::vector<std::string>&header) {
 std::map<std::string,int> out;
 for(int i=0;i<(int)header.size();++i)
   Check(out.emplace(header[i],i).second,"repeated CSV header");
 return out;
}
std::string Field(const std::vector<std::string>&row,
                 const std::map<std::string,int>&names,
                 const std::string&name) {
 const auto it=names.find(name);
 Check(it!=names.end()&&it->second<(int)row.size(),
       "missing selection column: "+name);
 return row[it->second];
}
bool Contains(const fs::path&path,const std::string&needle) {
 std::ifstream in(path);if(!in)return false;
 std::string row;
 while(std::getline(in,row))
   if(row.find(needle)!=std::string::npos)return true;
 return false;
}
int Channel(const std::string&s) {
 if(s=="jpsi40-bore")return 0;
 if(s=="phi40-bore")return 1;
 if(s=="omega40-bore")return 2;
 return -1;
}
const char* labels[]={"J/#psi","#phi","#omega"};
const char* plain[]={"jpsi40-bore","phi40-bore","omega40-bore"};
struct Selected {
 std::string id,channel;
 int station=-1,det=-1;
};
std::vector<Selected> ReadSelection(const std::string&path) {
 std::ifstream in(path);
 Check(bool(in),"selection file not readable "+path);
 std::string header,line;
 Check(bool(std::getline(in,header)),"empty frozen selection");
 const auto names=Cols(Split(header));
 std::vector<Selected> out;
 while(std::getline(in,line)) {
   if(line.empty())continue;
   const auto v=Split(line);
   Selected r;
   r.channel=Field(v,names,"source");
   r.station=std::stoi(Field(v,names,"station"));
   r.det=std::stoi(Field(v,names,"det_id"));
   r.id=r.channel+"|"+Field(v,names,"chunk")+"|"+Field(v,names,"event")+
        "|"+Field(v,names,"hit_index")+"|"+Field(v,names,"track_id")+
        "|"+Field(v,names,"station")+"|"+Field(v,names,"det_id");
   Check(Channel(r.channel)>=0&&r.station>=0&&r.station<6,
         "unexpected selected channel or station");
   out.push_back(std::move(r));
 }
 Check(out.size()==1000,"frozen selection is not 1000 crossings");
 return out;
}
struct Event {
 std::string id,channel;
 int station=0,det=0,primaries=0,clusters=0;
 long long ions=0,saved=0,collected=0;
 double angleU=0,angleW=0,pathMm=0,momentum=0;
 double gain=std::numeric_limits<double>::quiet_NaN();
};
Event Read(const fs::path&rootPath,const Selected&expect) {
 std::unique_ptr<TFile> file(TFile::Open(rootPath.c_str(),"READ"));
 Check(file&&!file->IsZombie(),"cannot read sealed ROOT "+rootPath.string());
 auto* mu=dynamic_cast<TTree*>(file->Get("muons"));
 auto* pr=dynamic_cast<TTree*>(file->Get("primaries"));
 auto* ion=dynamic_cast<TTree*>(file->Get("ions"));
 auto* cl=dynamic_cast<TTree*>(file->Get("clusters"));
 Check(mu&&pr&&ion&&cl&&mu->GetEntries()==1,
       "sealed ROOT missing expected trees or single muon");
 int reservoir=0;
 std::string *id=nullptr,*source=nullptr,*mode=nullptr;
 Event ev;
 mu->SetBranchAddress("crossing_id",&id);
 mu->SetBranchAddress("source_channel",&source);
 mu->SetBranchAddress("storage_mode",&mode);
 mu->SetBranchAddress("station",&ev.station);
 mu->SetBranchAddress("det_id",&ev.det);
 mu->SetBranchAddress("reservoir_per_primary",&reservoir);
 mu->SetBranchAddress("theta_u_deg",&ev.angleU);
 mu->SetBranchAddress("theta_w_deg",&ev.angleW);
 mu->SetBranchAddress("path_length_mm",&ev.pathMm);
 mu->SetBranchAddress("momentum_GeV",&ev.momentum);
 Check(mu->GetEntry(0)>0&&id&&source&&mode,
       "sealed ROOT missing muon metadata");
 ev.id=*id;ev.channel=*source;
 Check(ev.id==expect.id&&ev.channel==expect.channel&&
       ev.station==expect.station&&ev.det==expect.det &&
       *mode=="weighted"&&reservoir==256,
       "archive identity/physics mode differs from frozen selection");
 ev.clusters=(int)cl->GetEntries();
 int primary=0,successful=0;long long births=0,saved=0,collected=0;
 pr->SetBranchAddress("primary_index",&primary);
 pr->SetBranchAddress("successful",&successful);
 pr->SetBranchAddress("ion_births_exact",&births);
 pr->SetBranchAddress("ion_births_saved",&saved);
 pr->SetBranchAddress("collected_electrons",&collected);
 ev.primaries=(int)pr->GetEntries();
 for(int i=0;i<ev.primaries;++i) {
   Check(pr->GetEntry(i)>0&&primary==i&&successful==1&&
         births>=0&&saved==std::min<long long>(births,256)&&collected>=0,
         "unexpected failed primary or invalid ion census");
   ev.ions+=births;ev.saved+=saved;ev.collected+=collected;
 }
 Check(ev.saved==ion->GetEntries(),"saved ion count differs from ROOT rows");
 Check(std::isfinite(ev.pathMm)&&ev.pathMm>0 &&
       std::isfinite(ev.momentum)&&ev.momentum>=0,
       "invalid muon kinematics");
 if(ev.primaries>0)ev.gain=double(ev.collected)/ev.primaries;
 return ev;
}
struct Group {
 long long count=0,primaries=0,ions=0,collected=0,zeroPrimary=0;
 double gainSum=0,gain2=0;long long nGains=0;
 void Add(const Event& e) {
   ++count;primaries+=e.primaries;ions+=e.ions;collected+=e.collected;
   if(e.primaries==0)++zeroPrimary;
   else {gainSum+=e.gain;gain2+=e.gain*e.gain;++nGains;}
 }
};
} // namespace

int main(int argc,char**argv) {
 try{
 const auto folder=Arg(argc,argv,"--data");
 const auto selection=Arg(argc,argv,"--selection");
 const auto prefix=Arg(argc,argv,"--out-prefix");
 const auto stopString=Arg(argc,argv,"--stop-index");
 Check(!folder.empty()&&!selection.empty()&&!prefix.empty()&&
       !stopString.empty(),
       "usage: --data OVERNIGHT_RUN_DIR --selection SELECTION_CSV "
       "--stop-index EXCLUSIVE_UP_TO_1000 --out-prefix PREFIX");
 const int stop=std::stoi(stopString);
 Check(stop>=1&&stop<=1000,"stop index must be 1..1000 exclusive");
 const auto selected=ReadSelection(selection);
 const fs::path base(folder);
 Check(fs::is_directory(base),"overnight directory absent");
 for(const char*ext:{"_events.csv","_excluded.csv","_groups.csv","_hist.root","_overview.png"})
   Check(!fs::exists(prefix+ext),"refusing to overwrite "+prefix+ext);
 std::ofstream events(prefix+"_events.csv");
 std::ofstream excluded(prefix+"_excluded.csv");
 std::ofstream groups(prefix+"_groups.csv");
 Check(bool(events)&&bool(excluded)&&bool(groups),
       "cannot create result CSV; create output directory first");
 events<<std::setprecision(15)
       <<"selection_index,crossing_id,channel,station,det_id,heed_clusters,"
         "heed_primaries,true_ion_births,retained_K256_ion_rows,"
         "collected_electrons,effective_gain_per_primary,"
         "theta_u_deg,theta_w_deg,path_mm,momentum_GeV\n";
 excluded<<"selection_index,crossing_id,channel,station,reason\n";
 std::array<std::array<Group,6>,3> by{};
 std::array<Group,3> channels{};
 long long verified=0,unsealed=0,missing=0,invalid=0;
 TH1D logGain[3]={
  TH1D("jpsi_loggain","J/#psi;log_{10}(1 + collected e^{-}/Heed primary);Events",100,0,8),
  TH1D("phi_loggain","#phi;log_{10}(1 + collected e^{-}/Heed primary);Events",100,0,8),
  TH1D("omega_loggain","#omega;log_{10}(1 + collected e^{-}/Heed primary);Events",100,0,8)};
 TH1D logBirth[3]={
  TH1D("jpsi_logbirth","J/#psi;log_{10}(1 + avalanche ion births);Events",100,0,9),
  TH1D("phi_logbirth","#phi;log_{10}(1 + avalanche ion births);Events",100,0,9),
  TH1D("omega_logbirth","#omega;log_{10}(1 + avalanche ion births);Events",100,0,9)};
 TH2D counts("accepted_by_station_channel",
            "Audited 2100V event coverage;station (0 to 5);channel",
            6,-.5,5.5,3,-.5,2.5);
 for(int c=0;c<3;++c)counts.GetYaxis()->SetBinLabel(c+1,labels[c]);
 for(int idx=0;idx<stop;++idx){
   const auto& sel=selected[idx];
   std::ostringstream tag;tag<<"muon_"<<std::setfill('0')<<std::setw(4)<<idx;
   const fs::path dir=base/tag.str();
   std::string reason;
   if(!fs::exists(dir)) {reason="ABSENT";++missing;}
   else if(!fs::exists(dir/"complete.txt")) {
     reason=fs::exists(dir/"FAILED")?"FAILED_UNSEALED":"UNSEALED_OR_ACTIVE";
     ++unsealed;
   }else if(!fs::is_regular_file(dir/"ion_births.root") ||
            !fs::exists(dir/"SHA256SUMS") ||
            !Contains(dir/"root_audit.log","PASS CATHODE ROOT PHYSICS AUDIT") ||
            !Contains(dir/"physics.log","GARFIELD_AVALANCHE_SIZE_LIMIT_DISABLED=1")) {
     reason="INVALID_SEAL_OR_AUDIT";++invalid;
   }else{
     try {
       const auto e=Read(dir/"ion_births.root",sel);
       const int c=Channel(e.channel);
       Check(c>=0,"unexpected channel");
       ++verified;
       by[c][e.station].Add(e);channels[c].Add(e);
       counts.Fill(e.station,c);
       if(e.primaries>0)logGain[c].Fill(std::log10(1+e.gain));
       logBirth[c].Fill(std::log10(1+double(e.ions)));
       events<<idx<<","<<e.id<<","<<e.channel<<","<<e.station<<","<<e.det
             <<","<<e.clusters<<","<<e.primaries<<","<<e.ions<<","<<e.saved
             <<","<<e.collected<<",";
       if(e.primaries>0)events<<e.gain;else events<<"NA";
       events<<","<<e.angleU<<","<<e.angleW<<","<<e.pathMm
             <<","<<e.momentum<<"\n";
     }catch(const std::exception&e){
       reason="SEALED_ROOT_INCONSISTENT";
       ++invalid;
       std::cerr<<"EXCLUDED index "<<idx<<": "<<e.what()<<"\n";
     }
   }
   if(!reason.empty())excluded<<idx<<","<<sel.id<<","<<sel.channel
                              <<","<<sel.station<<","<<reason<<"\n";
 }
 groups<<"channel,station,verified_muons,zero_Heed_primary_muons,"
         "sum_Heed_primaries,sum_avalanche_ion_births,sum_collected_electrons,"
         "mean_event_gain_per_primary,rms_event_gain_per_primary\n";
 for(int c=0;c<3;++c){
   for(int st=0;st<7;++st){
     const Group&g=st==6?channels[c]:by[c][st];
     groups<<plain[c]<<","<<(st==6?"ALL":std::to_string(st))
           <<","<<g.count<<","<<g.zeroPrimary<<","<<g.primaries
           <<","<<g.ions<<","<<g.collected<<",";
     if(g.nGains) {
       double mean=g.gainSum/g.nGains;
       double rms=std::sqrt(std::max(0.,g.gain2/g.nGains-mean*mean));
       groups<<std::setprecision(14)<<mean<<","<<rms<<"\n";
     }else groups<<"NA,NA\n";
   }
 }
 events.close();excluded.close();groups.close();

 gROOT->SetBatch(kTRUE);gStyle->SetOptStat(0);
 TFile root((prefix+"_hist.root").c_str(),"RECREATE",
            "NA60+ audited 2100V K256 muon statistics",4);
 Check(!root.IsZombie(),"cannot write ROOT histograms");
 counts.Write();
 for(int c=0;c<3;++c){logGain[c].Write();logBirth[c].Write();}
 TCanvas canvas("overview","NA60+ audited muons",1500,550);
 canvas.Divide(3,1);
 canvas.cd(1);
 counts.SetMinimum(0);counts.Draw("COLZ TEXT");
 canvas.cd(2);
 logGain[0].SetLineColor(kBlack);
 logGain[1].SetLineColor(kBlue+1);
 logGain[2].SetLineColor(kRed+1);
 logGain[0].SetLineWidth(2);logGain[1].SetLineWidth(2);logGain[2].SetLineWidth(2);
 logGain[0].SetTitle("Collected-electron gain;log_{10}(1+G_{eff});muons");
 logGain[0].Draw("HIST");logGain[1].Draw("HIST SAME");logGain[2].Draw("HIST SAME");
 TLegend gainLegend(.64,.68,.9,.88);gainLegend.SetFillStyle(0);
 for(int c=0;c<3;++c)gainLegend.AddEntry(&logGain[c],labels[c],"l");
 gainLegend.DrawClone();
 canvas.cd(3);
 for(int c=0;c<3;++c){
   logBirth[c].SetLineColor(c==0?kBlack:(c==1?kBlue+1:kRed+1));
   logBirth[c].SetLineWidth(2);
 }
 logBirth[0].SetTitle("Avalanche ion births;log_{10}(1+N_{ion});muons");
 logBirth[0].Draw("HIST");
 logBirth[1].Draw("HIST SAME");
 logBirth[2].Draw("HIST SAME");
 canvas.cd(0);canvas.Modified();canvas.Update();
 root.cd();canvas.Write();
 canvas.SaveAs((prefix+"_overview.png").c_str());
 root.Close();
 std::cout<<"PASS read-only audited muon statistics, selection slots [0,"
          <<stop<<") verified="<<verified<<" unsealed="<<unsealed
          <<" absent="<<missing<<" invalid="<<invalid<<"\n";
 for(int c=0;c<3;++c)std::cout<<"channel="<<plain[c]
        <<" verified="<<channels[c].count
        <<" primaries="<<channels[c].primaries
        <<" ion_births="<<channels[c].ions
        <<" collected_electrons="<<channels[c].collected<<"\n";
 std::cout<<"Files: "<<prefix<<"_events.csv, _excluded.csv, _groups.csv,"
          <<" _hist.root, _overview.png\n";
 std::cout<<"NOTE no full ion trajectories, no cathode currents, no electronics. "
          <<"The sample excludes failed/unsealed events explicitly.\n";
 return 0;
 }catch(const std::exception&e) {
 std::cerr<<e.what()<<"\n";return 2;
 }
}
