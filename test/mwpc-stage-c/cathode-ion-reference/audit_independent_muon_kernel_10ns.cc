// Stage C independent-muon generalization check for the FIXED 256-template
// response kernel. Three genuinely different audited muons; not template
// training muon0. Read-only microscopic ROOT, direct RKF and IDW outputs.
// No electronics, new avalanches, new ion sampling, or readout strips.
#include <TFile.h>
#include <TTree.h>
#include <TNamed.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TStyle.h>
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
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs=std::filesystem;
namespace {
void Need(bool b,const std::string&s){
 if(!b)throw std::runtime_error("INDEPENDENT MUON KERNEL AUDIT: "+s);
}
std::string Arg(int argc,char**argv,const std::string&k){
 for(int i=1;i+1<argc;++i)if(k==argv[i])return argv[i+1];
 return "";
}
std::vector<std::string> Split(const std::string&s) {
 std::vector<std::string> v;std::string a;std::stringstream ss(s);
 while(std::getline(ss,a,',')) {
  if(!a.empty()&&a.back()=='\r')a.pop_back();
  v.push_back(a);
 }
 return v;
}
std::map<std::string,std::string> Summary(const std::string&path) {
 std::ifstream f(path);std::string h,row;
 Need(bool(f)&&bool(std::getline(f,h))&&bool(std::getline(f,row)),
      "missing or empty summary "+path);
 auto keys=Split(h),vals=Split(row);
 Need(keys.size()==vals.size(),"summary column mismatch "+path);
 std::map<std::string,std::string> out;
 for(std::size_t i=0;i<keys.size();++i)
  Need(out.emplace(keys[i],vals[i]).second,"duplicate summary key");
 return out;
}
std::string V(const std::map<std::string,std::string>&m,const std::string&key) {
 const auto it=m.find(key);
 Need(it!=m.end(),"missing summary key "+key);
 return it->second;
}
std::array<std::vector<double>,4> ReadFast(const std::string&file) {
 std::ifstream f(file);std::string h,row;
 Need(bool(f)&&bool(std::getline(f,h)),
      "missing event waveform "+file);
 Need(h=="time_center_ns,plus_kernel_fC_per_ns,minus_kernel_fC_per_ns,"
         "plus_direct_fC_per_ns,minus_direct_fC_per_ns,"
         "plus_kernel_cumulative_fC,minus_kernel_cumulative_fC,"
         "plus_direct_cumulative_fC,minus_direct_cumulative_fC",
      "wrong fast waveform CSV schema "+file);
 std::array<std::vector<double>,4> signal;
 for(auto&v:signal)v.reserve(10000);
 for(int bin=0;bin<10000;++bin) {
  Need(bool(std::getline(f,row)),"truncated fast waveform "+file);
  const auto parts=Split(row);
  Need(parts.size()==9,"malformed 10ns waveform row");
  Need(std::abs(std::stod(parts[0])-(bin+.5)*10.)<1.e-7,
       "10ns timestamp mismatch");
  for(int c=0;c<4;++c) {
   const double y=std::stod(parts[1+c]);
   Need(std::isfinite(y),"nonfinite response current");
   signal[c].push_back(y);
  }
 }
 Need(!std::getline(f,row),"extra waveform row beyond 100us");
 return signal;
}
double Charge(const std::vector<double>&a,int ns){
 Need(ns%10==0&&ns>=0&&ns<=100000,"unsupported ns window");
 return 10.*std::accumulate(a.begin(),a.begin()+ns/10,0.);
}
double RelL2(const std::vector<double>&a,const std::vector<double>&b) {
 Need(a.size()==b.size(),"L2 shape mismatch");
 double e=0,r=0;
 for(std::size_t j=0;j<a.size();++j) {
  const double d=a[j]-b[j];e+=d*d;r+=b[j]*b[j];
 }
 return std::sqrt(e/std::max(r,1.e-50));
}
double RebinnedL2(const std::vector<double>&a,const std::vector<double>&b){
 Need(a.size()==10000&&b.size()==10000,"expected 10000 raw bins");
 std::vector<double> ar(1000),br(1000);
 for(int j=0;j<1000;++j)
  for(int k=0;k<10;++k){
   ar[j]+=a[10*j+k]/10.;
   br[j]+=b[10*j+k]/10.;
  }
 return RelL2(ar,br);
}
const std::array<int,3> events={301,408,446};
struct Metadata {
 std::string crossing,channel;
 int station=0,primaries=0;
 long long exact=0,stored=0,collected=0;
};
Metadata Meta(const fs::path&path){
 std::unique_ptr<TFile> f(TFile::Open(path.c_str(),"READ"));
 Need(f&&!f->IsZombie(),"cannot open audited muon archive "+path.string());
 TTree* mu=dynamic_cast<TTree*>(f->Get("muons"));
 TTree* prim=dynamic_cast<TTree*>(f->Get("primaries"));
 TTree* ions=dynamic_cast<TTree*>(f->Get("ions"));
 Need(mu&&prim&&ions&&mu->GetEntries()==1,
      "missing frozen muon census");
 Metadata m;
 std::string *id=nullptr,*channel=nullptr;
 mu->SetBranchAddress("crossing_id",&id);
 mu->SetBranchAddress("source_channel",&channel);
 mu->SetBranchAddress("station",&m.station);
 Need(mu->GetEntry(0)>0&&id&&channel,"invalid muon metadata");
 m.crossing=*id;m.channel=*channel;
 long long exact=0,saved=0,elec=0;int successful=0;
 prim->SetBranchAddress("ion_births_exact",&exact);
 prim->SetBranchAddress("ion_births_saved",&saved);
 prim->SetBranchAddress("collected_electrons",&elec);
 prim->SetBranchAddress("successful",&successful);
 m.primaries=(int)prim->GetEntries();
 for(int i=0;i<m.primaries;++i) {
  Need(prim->GetEntry(i)>0&&successful==1&&
       exact>=0&&saved==std::min<long long>(exact,256),
       "archived primary failed physics audit");
  m.exact+=exact;m.stored+=saved;m.collected+=elec;
 }
 Need(m.stored==ions->GetEntries(),"K256 sample row census inconsistent");
 return m;
}
}
int main(int argc,char**argv){
 try{
  const auto data=Arg(argc,argv,"--data");
  const auto results=Arg(argc,argv,"--results");
  const auto training=Arg(argc,argv,"--training-root");
  const auto prefix=Arg(argc,argv,"--out-prefix");
  Need(!data.empty()&&!results.empty()&&!training.empty()&&!prefix.empty(),
       "usage: --data OVERNIGHT_DIR --results ANALYSIS_DIR "
       "--training-root 10NS_KERNEL.root --out-prefix PREFIX");
  std::unique_ptr<TFile> train(TFile::Open(training.c_str(),"READ"));
  Need(train&&!train->IsZombie(),"cannot read fixed training ROOT");
  auto* bank=dynamic_cast<TTree*>(train->Get("training_bank"));
  auto* original=dynamic_cast<TNamed*>(train->Get("input_ion_archive"));
  Need(bank&&original&&bank->GetEntries()==256,
       "require unchanged 256-template training archive");
  const auto oldSource=std::string(original->GetTitle());
  Need(oldSource.find("muon_0000")!=std::string::npos ||
       oldSource.find("muon0")!=std::string::npos,
       "training source not identifiable as original muon0");
  for(auto ext:{"_events.csv","_windows.csv","_diagnostic.root","_early.png"})
    Need(!fs::exists(prefix+ext),"refusing to overwrite result "+prefix+ext);
  std::ofstream output(prefix+"_events.csv"),windows(prefix+"_windows.csv");
  Need(bool(output)&&bool(windows),"cannot write analysis csv");
  output<<std::setprecision(16)
    <<"index,channel,station,crossing_id,n_Heed_primaries,true_ions,"
      "saved_K256_ions,collected_electrons,effective_gain_per_Heed,"
      "plus_native10_L2,minus_native10_L2,plus_rebinned100_L2,"
      "minus_rebinned100_L2,plus_final_direct_fC,plus_final_fast_fC,"
      "minus_final_direct_fC,minus_final_fast_fC\n";
  windows<<std::setprecision(16)
    <<"index,channel,station,window_ns,plane,direct_ion_fC,"
      "fast_ion_fC,difference_fC,relative_difference_percent\n";
  gROOT->SetBatch(kTRUE);gStyle->SetOptStat(0);
  TFile fout((prefix+"_diagnostic.root").c_str(),"RECREATE",
       "Independent 3 muons: fixed K256 avalanche and 256 RKF templates");
  Need(!fout.IsZombie(),"cannot create diagnostics ROOT");
  TNamed limitations("limitations",
    "Independent muons selected from audited 1000 selection: indices "
    "301(phi),408(jpsi),446(omega). Uncapped microscopic avalanche, 2100V. "
    "K256 per-primary sampled ion births, weighted direct RKF vs "
    "unchanged 256-position/4-neighbor single-ion template bank "
    "at 10ns time bins. Two complete ideal cathode planes only, "
    "no electron current or electronics, no physical strips, "
    "provisional Ar+ in Ar ion mobility.");
  limitations.Write();
  TCanvas canvas("independent_muon_10ns","Independent 10ns muon kernel",1600,1100);
  canvas.Divide(2,3);
  std::vector<std::unique_ptr<TH1D>> hists;
  double totalError2=0,totalRef2=0;
  for(std::size_t i=0;i<events.size();++i){
   const int idx=events[i];
   std::ostringstream tag;tag<<"muon_"<<std::setfill('0')<<std::setw(4)<<idx;
   const fs::path archive=fs::path(data)/tag.str()/"ion_births.root";
   const auto m=Meta(archive);
   const std::string base=std::string(results)+"/"+tag.str();
   const auto direct=Summary(base+"_direct_dt10_summary.csv");
   const auto fast=Summary(base+"_fast_dt10_summary.csv");
   for(const auto& sm:{direct,fast}) {
    Need(V(sm,"crossing_id")==m.crossing &&
         std::stoi(V(sm,"station"))==m.station &&
         std::stoll(V(sm,"true_ion_births"))==m.exact &&
         std::abs(std::stod(V(sm,"dt_ns"))-10.)<1.e-10 &&
         std::abs(std::stod(V(sm,"window_us"))-100.)<1.e-10,
         "different physical muon or incorrect time grid in summary");
   }
   Need(std::stoll(V(direct,"retained_ion_rows"))==m.stored &&
        std::stoll(V(fast,"retained_macroions"))==m.stored &&
        std::stoi(V(fast,"training_templates"))==512 &&
        std::stoi(V(fast,"neighbors"))==4,
        "different reservoir or spatial interpolation configuration");
   const auto w=ReadFast(base+"_fast_dt10_waveforms.csv");
   for(int side=0;side<2;++side) {
    const auto& kernel=w[side];const auto& ref=w[side+2];
    const double L2=RelL2(kernel,ref);
    const double rebL2=RebinnedL2(kernel,ref);
    if(side==0){
      output<<idx<<","<<m.channel<<","<<m.station<<","<<m.crossing
            <<","<<m.primaries<<","<<m.exact<<","<<m.stored
            <<","<<m.collected<<","
            <<(m.primaries?double(m.collected)/m.primaries:0.)
            <<","<<L2<<",";
    }else{
      // Output row written in one pass below.
    }
    for(int window:{50,100,200,400,1000,100000}) {
      const double q0=Charge(ref,window),q1=Charge(kernel,window);
      windows<<idx<<","<<m.channel<<","<<m.station<<","
             <<window<<","<<(side?"minus":"plus")<<","
             <<q0<<","<<q1<<","<<q1-q0<<","
             <<(std::abs(q0)>1.e-15?100.*(q1-q0)/q0:
                 std::numeric_limits<double>::quiet_NaN())<<"\n";
    }
    for(int bin=0;bin<10000;++bin){
      const double d=kernel[bin]-ref[bin];
      totalError2+=d*d;totalRef2+=ref[bin]*ref[bin];
    }
    auto h0=std::make_unique<TH1D>(
       (tag.str()+"_direct_"+std::to_string(side)).c_str(),
       (m.channel+" station "+std::to_string(m.station)+
        (side?" lower":" upper")+";time since muon (ns);ion current (fC/ns)").c_str(),
       100,0.,1000.);
    auto h1=std::make_unique<TH1D>(
       (tag.str()+"_fast_"+std::to_string(side)).c_str(),
       h0->GetTitle(),100,0.,1000.);
    h0->SetDirectory(nullptr);h1->SetDirectory(nullptr);
    for(int bin=0;bin<100;++bin){
      h0->SetBinContent(bin+1,ref[bin]);h1->SetBinContent(bin+1,kernel[bin]);
    }
    h0->SetLineColor(kBlack);h0->SetLineWidth(2);
    h1->SetLineColor(kBlue+1);h1->SetLineWidth(2);
    canvas.cd(i*2+side+1);
    h0->Draw("HIST");h1->Draw("HIST SAME");
    TLegend leg(.52,.72,.89,.88);
    leg.SetBorderSize(0);leg.SetFillStyle(0);
    leg.AddEntry(h0.get(),"Direct weighted RKF","l");
    leg.AddEntry(h1.get(),"Fixed fast template kernel","l");
    leg.DrawClone();
    fout.cd();h0->Write();h1->Write();
    hists.push_back(std::move(h0));hists.push_back(std::move(h1));
   }
   // Event one-row metrics in fixed upper/lower order.
   const double upL2=RelL2(w[0],w[2]),dnL2=RelL2(w[1],w[3]);
   const double upReb=RebinnedL2(w[0],w[2]),dnReb=RebinnedL2(w[1],w[3]);
   // Write complete row below; no partial records.
   output<<dnL2<<","<<upReb<<","<<dnReb<<","
         <<Charge(w[2],100000)<<","<<Charge(w[0],100000)<<","
         <<Charge(w[3],100000)<<","<<Charge(w[1],100000)<<"\n";
   std::cout<<std::setprecision(7)<<"independent "<<tag.str()
       <<" source="<<m.channel<<" gain="
       <<(m.primaries?double(m.collected)/m.primaries:0.)
       <<" upper_10ns_L2="<<upL2<<" lower_10ns_L2="<<dnL2
       <<" upper_100ns_rebinned_L2="<<upReb
       <<" lower_100ns_rebinned_L2="<<dnReb<<"\n";
  }
  output.close();windows.close();
  canvas.cd(0);canvas.Modified();canvas.Update();
  fout.cd();canvas.Write();canvas.SaveAs((prefix+"_early.png").c_str());fout.Close();
  std::cout<<"PASS independent same-K256 10ns kernel check, three genuinely "
    "distinct 2100V muons, combined_10ns_waveform_L2="
    <<std::sqrt(totalError2/std::max(totalRef2,1.e-50))
    <<" (not validation of real strips, ions in Ar/CO2, or electronics)\n";
  return 0;
 }catch(const std::exception&e){
  std::cerr<<e.what()<<"\n";return 2;
 }
}
