// Same-muon 100ns -> 10ns comparison, two ideal WHOLE cathodes, ion only.
// Uses FOUR existing waveform files: direct RKF and IDW fast ion for each dt.
// Time rebin is CHARGE CONSERVING (sum 10 consecutive 10ns current bins).
// Optional microscopic electron ROOT adds integrated electron context, not shaper.
// Does not run Garfield++, read/write overnight muon campaign, or calibrate VMM.
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TNamed.h>
#include <TROOT.h>
#include <TStyle.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
namespace {
using Sides=std::array<std::vector<double>,2>;
void Check(bool b,const std::string&s) {
 if(!b)throw std::runtime_error("PLANE TIMING AUDIT: "+s);
}
std::string Arg(int argc,char**argv,const std::string&k) {
 for(int i=1;i+1<argc;++i)if(k==argv[i])return argv[i+1];
 return "";
}
std::vector<std::string> Split(const std::string&s) {
 std::vector<std::string> cols;std::string x;std::istringstream in(s);
 while(std::getline(in,x,','))cols.push_back(x);
 return cols;
}
std::map<std::string,std::string> LoadSummary(const std::string&path) {
 std::ifstream f(path);std::string head,line;
 Check(bool(f)&&bool(std::getline(f,head))&&bool(std::getline(f,line)),
       "cannot read summary "+path);
 auto keys=Split(head),values=Split(line);
 Check(keys.size()==values.size(),"summary column count mismatch");
 std::map<std::string,std::string> out;
 for(std::size_t i=0;i<keys.size();++i)
   Check(out.emplace(keys[i],values[i]).second,"duplicate summary column");
 return out;
}
std::string Value(const std::map<std::string,std::string>&m,
                  const std::string&k) {
 const auto it=m.find(k);
 Check(it!=m.end(),"summary missing "+k);
 return it->second;
}
std::string SummaryOf(const std::string&waveform) {
 const std::string suffix="_waveforms.csv";
 Check(waveform.size()>suffix.size()&&
    waveform.substr(waveform.size()-suffix.size())==suffix,
    "expected file ending _waveforms.csv: "+waveform);
 return waveform.substr(0,waveform.size()-suffix.size())+"_summary.csv";
}
Sides ReadWave(const std::string&file,bool event,double dt,int bins) {
 std::ifstream f(file);std::string line;
 Check(bool(f)&&bool(std::getline(f,line)),"missing waveform "+file);
 const std::string directHead=
   "time_center_ns,ion_plus_fC_per_ns,ion_minus_fC_per_ns,"
   "cumulative_plus_fC,cumulative_minus_fC";
 const std::string eventHead=
   "time_center_ns,plus_kernel_fC_per_ns,minus_kernel_fC_per_ns,"
   "plus_direct_fC_per_ns,minus_direct_fC_per_ns,"
   "plus_kernel_cumulative_fC,minus_kernel_cumulative_fC,"
   "plus_direct_cumulative_fC,minus_direct_cumulative_fC";
 Check(line==(event?eventHead:directHead),"unexpected waveform CSV schema "+file);
 Sides current;
 for(auto &v:current)v.reserve(bins);
 for(int i=0;i<bins;++i) {
   Check(bool(std::getline(f,line)),"waveform ended early "+file);
   const auto c=Split(line);
   Check(c.size()==(event?9u:5u),"wrong waveform column count");
   const double t=std::stod(c[0]);
   Check(std::isfinite(t)&&std::abs(t-(i+0.5)*dt)<1.e-7*dt,
         "waveform time bin does not match expected grid");
   for(int side=0;side<2;++side) {
     const double value=std::stod(c[side+1]);
     Check(std::isfinite(value),"nonfinite current in waveform "+file);
     current[side].push_back(value);
   }
 }
 Check(!std::getline(f,line),"unexpected extra waveform bins "+file);
 return current;
}
Sides ReadEventDirect(const std::string&file,double dt,int bins) {
 std::ifstream f(file);std::string line;
 Check(bool(f)&&bool(std::getline(f,line)),"cannot read event waveforms "+file);
 Sides current;
 for(auto &v:current)v.reserve(bins);
 for(int i=0;i<bins;++i) {
   Check(bool(std::getline(f,line)),"event waveforms truncated");
   auto c=Split(line);
   Check(c.size()==9u,"event waveforms 9 columns required");
   Check(std::abs(std::stod(c[0])-(i+0.5)*dt)<1.e-7*dt,
         "event waveforms time grid mismatch");
   for(int side=0;side<2;++side)
     current[side].push_back(std::stod(c[side+3]));
 }
 Check(!std::getline(f,line),"unexpected extra event waveform bins");
 return current;
}
Sides Rebin(const Sides&fine,int ratio) {
 Check(ratio>=1&&fine[0].size()==fine[1].size() &&
       fine[0].size()%(std::size_t)ratio==0,
       "cannot rebin fine waveform");
 Sides coarse;
 for(int s=0;s<2;++s){
   const auto n=fine[s].size()/ratio;
   coarse[s].assign(n,0.);
   for(std::size_t j=0;j<n;++j)
     for(int k=0;k<ratio;++k)
       coarse[s][j]+=fine[s][j*ratio+k]/ratio;
 }
 return coarse;
}
double SumCharge(const std::vector<double>&v,int n,double dt) {
 Check(n>=0&&n<=(int)v.size(),"integral extends beyond current window");
 return std::accumulate(v.begin(),v.begin()+n,0.)*dt;
}
double RelativeL2(const std::vector<double>&a,const std::vector<double>&b){
 Check(a.size()==b.size(),"L2 vector lengths differ");
 double e=0,r=0;
 for(std::size_t i=0;i<a.size();++i){
   e+=(a[i]-b[i])*(a[i]-b[i]);
   r+=b[i]*b[i];
 }
 return std::sqrt(e/std::max(r,1.e-60));
}
struct Electron {
 std::string crossing;
 Sides currents;
 bool present=false;
};
Electron ReadElectron(const std::string&path){
 Electron e;
 if(path.empty())return e;
 std::unique_ptr<TFile> file(TFile::Open(path.c_str(),"READ"));
 Check(file&&!file->IsZombie(),"cannot read complete electron ROOT "+path);
 auto* meta=dynamic_cast<TTree*>(file->Get("muon"));
 auto* signals=dynamic_cast<TTree*>(file->Get("electron_signal"));
 Check(meta&&signals&&meta->GetEntries()==1,"electron ROOT missing audited signal");
 std::string* id=nullptr;double dt=0,window=0;
 meta->SetBranchAddress("crossing_id",&id);
 meta->SetBranchAddress("dt_ns",&dt);
 meta->SetBranchAddress("time_window_ns",&window);
 Check(meta->GetEntry(0)>0&&id&&
       std::abs(dt-0.05)<1.e-10&&std::abs(window-200.)<1.e-10,
       "electron ROOT incompatible with archived 0.05ns / 200ns bins");
 e.crossing=*id;
 Check(signals->GetEntries()==4000,"electron signal has wrong number of bins");
 int index=0;double t=0,plus=0,minus=0;
 signals->SetBranchAddress("time_bin",&index);
 signals->SetBranchAddress("time_center_ns",&t);
 signals->SetBranchAddress("plus_electron_fC_per_ns",&plus);
 signals->SetBranchAddress("minus_electron_fC_per_ns",&minus);
 for(int i=0;i<4000;++i) {
   Check(signals->GetEntry(i)>0&&index==i &&
       std::abs(t-(i+.5)*.05)<1.e-8&&
       std::isfinite(plus)&&std::isfinite(minus),
       "malformed electron signal");
   e.currents[0].push_back(plus);
   e.currents[1].push_back(minus);
 }
 e.present=true;
 return e;
}
double ElectronCharge(const Electron&e,int side,int windowNs) {
 if(!e.present)return 0.;
 const int n=std::min(windowNs*20,4000);
 return SumCharge(e.currents[side],n,.05);
}
void PlotEarly(const std::string&prefix,const Sides&oldDirect,
     const Sides&oldFast,const Sides&newDirect,const Sides&newFast,
     TFile&root) {
 gROOT->SetBatch(kTRUE);gStyle->SetOptStat(0);
 TCanvas canvas("early_100ns_10ns",
                "NA60+ 2100V ion response timing validation",1450,900);
 canvas.Divide(2,2);
 std::vector<std::unique_ptr<TH1D>> stored;
 stored.reserve(16);
 for(int side=0;side<2;++side)for(int kind=0;kind<2;++kind) {
   const int count=kind?100:50; // first 1000ns cumulative / 500ns current
   auto make=[&](const std::string&name) {
     auto hist=std::make_unique<TH1D>(
       (name+"_"+std::to_string(side)+"_"+std::to_string(kind)).c_str(),
       ((side?"Lower":"Upper")+std::string(" plane;time (ns);")+
        (kind?"cumulative induced charge (fC)":"ion current (fC/ns)")).c_str(),
       count,0,count*10.);
     hist->SetDirectory(nullptr);
     return hist;
   };
   auto fineDirect=make("fine_direct"),fineFast=make("fine_fast");
   auto coarseDirect=make("coarse_direct"),coarseFast=make("coarse_fast");
   double q0=0,q1=0,q2=0,q3=0;
   for(int i=0;i<count;++i) {
     const int coarse=i/10;
     const double rNew=newDirect[side][i],fNew=newFast[side][i];
     const double rOld=oldDirect[side][coarse],fOld=oldFast[side][coarse];
     q0+=rNew*10.;q1+=fNew*10.;q2+=rOld*10.;q3+=fOld*10.;
     fineDirect->SetBinContent(i+1,kind?q0:rNew);
     fineFast->SetBinContent(i+1,kind?q1:fNew);
     coarseDirect->SetBinContent(i+1,kind?q2:rOld);
     coarseFast->SetBinContent(i+1,kind?q3:fOld);
   }
   fineDirect->SetLineColor(kBlack);fineDirect->SetLineWidth(2);
   fineFast->SetLineColor(kBlue+1);fineFast->SetLineWidth(2);
   coarseDirect->SetLineColor(kOrange+7);coarseDirect->SetLineStyle(2);
   coarseFast->SetLineColor(kGreen+2);coarseFast->SetLineStyle(2);
   double low=0,high=0;
   for(auto *h:{fineDirect.get(),fineFast.get(),
                coarseDirect.get(),coarseFast.get()}){
     low=std::min(low,h->GetMinimum());
     high=std::max(high,h->GetMaximum());
   }
   Check(std::isfinite(low)&&std::isfinite(high)&&low<high,
         "early timing plot has invalid values");
   fineDirect->SetMinimum(low-.12*(high-low));
   fineDirect->SetMaximum(high+.18*(high-low));
   canvas.cd(kind*2+side+1);
   fineDirect->Draw("HIST");fineFast->Draw("HIST SAME");
   coarseDirect->Draw("HIST SAME");coarseFast->Draw("HIST SAME");
   TLegend leg(.45,.65,.9,.89);
   leg.SetFillStyle(0);leg.SetBorderSize(0);
   leg.AddEntry(fineDirect.get(),"Direct RKF, 10 ns","l");
   leg.AddEntry(fineFast.get(),"Fast kernel, 10 ns","l");
   leg.AddEntry(coarseDirect.get(),"Direct RKF, 100 ns","l");
   leg.AddEntry(coarseFast.get(),"Fast kernel, 100 ns","l");
   leg.DrawClone();
   root.cd();
   fineDirect->Write();fineFast->Write();
   coarseDirect->Write();coarseFast->Write();
   stored.push_back(std::move(fineDirect));
   stored.push_back(std::move(fineFast));
   stored.push_back(std::move(coarseDirect));
   stored.push_back(std::move(coarseFast));
 }
 canvas.cd(0);canvas.Modified();canvas.Update();
 root.cd();canvas.Write();
 canvas.SaveAs((prefix+"_early_overlay.png").c_str());
}
} // namespace
int main(int argc,char**argv){
 try {
  const auto oldDirectPath=Arg(argc,argv,"--direct-100ns");
  const auto oldEventPath=Arg(argc,argv,"--event-100ns");
  const auto newDirectPath=Arg(argc,argv,"--direct-10ns");
  const auto newEventPath=Arg(argc,argv,"--event-10ns");
  const auto electronPath=Arg(argc,argv,"--electron-root");
  const auto prefix=Arg(argc,argv,"--out-prefix");
  Check(!oldDirectPath.empty()&&!oldEventPath.empty()&&
        !newDirectPath.empty()&&!newEventPath.empty()&&!prefix.empty(),
        "usage: --direct-100ns FILE --event-100ns FILE --direct-10ns FILE "
        "--event-10ns FILE --out-prefix PREFIX [--electron-root ROOT]");
  for(const auto&ext:{"_windows.csv","_metrics.csv","_early.root","_early_overlay.png"})
   Check(!std::filesystem::exists(prefix+ext),
         "refusing to overwrite existing timing comparison "+prefix+ext);
  const auto oldD=LoadSummary(SummaryOf(oldDirectPath));
  const auto oldF=LoadSummary(SummaryOf(oldEventPath));
  const auto fineD=LoadSummary(SummaryOf(newDirectPath));
  const auto fineF=LoadSummary(SummaryOf(newEventPath));
  const auto refId=Value(oldD,"crossing_id");
  const auto refIons=Value(oldD,"true_ion_births");
  const auto refSaved=Value(oldD,"retained_ion_rows");
  const auto refStation=Value(oldD,"station");
  const auto refDet=Value(oldD,"det_id");
  Check(Value(oldF,"crossing_id")==refId &&
        Value(fineD,"crossing_id")==refId &&
        Value(fineF,"crossing_id")==refId &&
        Value(oldF,"true_ion_births")==refIons &&
        Value(fineD,"true_ion_births")==refIons &&
        Value(fineF,"true_ion_births")==refIons &&
        Value(oldF,"retained_macroions")==refSaved &&
        Value(fineD,"retained_ion_rows")==refSaved &&
        Value(fineF,"retained_macroions")==refSaved,
        "input waveforms do not represent the SAME microscopic muon/reservoir");
  for(const auto&report:{oldD,oldF,fineD,fineF})
    Check(Value(report,"station")==refStation&&
          Value(report,"det_id")==refDet&&
          std::abs(std::stod(Value(report,"window_us"))-100.)<1.e-8,
          "mismatched station, detector or 100us window");
  Check(std::abs(std::stod(Value(oldD,"dt_ns"))-100.)<1.e-9&&
        std::abs(std::stod(Value(oldF,"dt_ns"))-100.)<1.e-9&&
        std::abs(std::stod(Value(fineD,"dt_ns"))-10.)<1.e-9&&
        std::abs(std::stod(Value(fineF,"dt_ns"))-10.)<1.e-9,
        "unexpected timing grid in summary");
  Check(Value(oldD,"archive_mode")=="weighted" &&
        Value(fineD,"archive_mode")=="weighted"&&
        Value(oldD,"response_mode")=="direct" &&
        Value(fineD,"response_mode")=="direct",
        "direct reference must be weighted RKF, not endpoint charge shortcut");

  auto rd100=ReadWave(oldDirectPath,false,100.,1000);
  auto fk100=ReadWave(oldEventPath,true,100.,1000);
  auto rd10=ReadWave(newDirectPath,false,10.,10000);
  auto fk10=ReadWave(newEventPath,true,10.,10000);
  auto eventOldRef=ReadEventDirect(oldEventPath,100.,1000);
  auto eventFineRef=ReadEventDirect(newEventPath,10.,10000);
  for(int side=0;side<2;++side)for(int bin=0;bin<1000;++bin)
    Check(std::abs(eventOldRef[side][bin]-rd100[side][bin]) <
      1.e-10*std::max(1.e-6,std::abs(rd100[side][bin])),
      "100ns event/direct RKF reference waveform disagrees");
  for(int side=0;side<2;++side)for(int bin=0;bin<10000;++bin)
    Check(std::abs(eventFineRef[side][bin]-rd10[side][bin]) <
      1.e-10*std::max(1.e-6,std::abs(rd10[side][bin])),
      "10ns event/direct RKF reference waveform disagrees");
  auto rd10to100=Rebin(rd10,10),fk10to100=Rebin(fk10,10);
  const auto electron=ReadElectron(electronPath);
  Check(!electron.present||electron.crossing==refId,
        "microscopic electrons do not match the muon crossing");

  std::ofstream windows(prefix+"_windows.csv");
  Check(bool(windows),"cannot write window comparisons");
  windows<<std::setprecision(16)
    <<"window_ns,plane,direct100_charge_fC,direct10_rebinned_charge_fC,"
      "kernel100_charge_fC,kernel10_rebinned_charge_fC,"
      "RKF_grid_difference_fC,old_kernel_minus_old_RKF_fC,"
      "fine_kernel_minus_fine_RKF_fC,"
      "microscopic_electron_charge_fC,"
      "fine_direct_electron_plus_ion_charge_fC,"
      "fine_kernel_electron_plus_ion_charge_fC\n";
  for(int window:{100,200,400,1000,10000,100000}) {
    const int N=window/100;
    for(int side=0;side<2;++side) {
      const double a=SumCharge(rd100[side],N,100.);
      const double b=SumCharge(rd10to100[side],N,100.);
      const double c=SumCharge(fk100[side],N,100.);
      const double d=SumCharge(fk10to100[side],N,100.);
      const double e=ElectronCharge(electron,side,window);
      windows<<window<<","<<(side?"minus":"plus")
             <<","<<a<<","<<b<<","<<c<<","<<d
             <<","<<a-b<<","<<c-a<<","<<d-b<<","
             <<e<<","<<b+e<<","<<d+e<<"\n";
    }
  }
  windows.close();
  std::ofstream metrics(prefix+"_metrics.csv");
  Check(bool(metrics),"cannot write L2 comparison");
  metrics<<std::setprecision(16)
    <<"plane,L2_old_kernel_vs_old_RKF_100ns,"
      "L2_fine_kernel_vs_fine_RKF_10ns,"
      "L2_fine_kernel_rebinned_vs_fine_RKF_rebinned_100ns,"
      "L2_old_RKF_vs_rebinned_new_RKF_100ns,"
      "L2_old_kernel_vs_new_RKF_rebinned_100ns,"
      "L2_fine_kernel_rebinned_vs_old_RKF_100ns\n";
  for(int side=0;side<2;++side){
    metrics<<(side?"minus":"plus")<<","
       <<RelativeL2(fk100[side],rd100[side])<<","
       <<RelativeL2(fk10[side],rd10[side])<<","
       <<RelativeL2(fk10to100[side],rd10to100[side])<<","
       <<RelativeL2(rd100[side],rd10to100[side])<<","
       <<RelativeL2(fk100[side],rd10to100[side])<<","
       <<RelativeL2(fk10to100[side],rd100[side])<<"\n";
  }
  metrics.close();

  TFile root((prefix+"_early.root").c_str(),"RECREATE",
      "NA60+ whole-plane 100ns vs 10ns ion timing audit",4);
  Check(!root.IsZombie(),"cannot write audit ROOT plots");
  TNamed scope("scope",
   "Same original K512 microscopic ion reservoir, identical 256 RKF "
   "training positions. Compares 100ns and 10ns same-muon direct and "
   "fast nearest-template ion-only currents. Optional full-muon direct "
   "electron ROOT adds integrated current for context. No VMM response; "
   "same fixed ideal nine-wire whole cathode geometry, provisional mobility. "
   "Individual errors: reference grid, spatial interpolation, temporal "
   "resampling remain to be interpreted with explicit tables.");
  scope.Write();PlotEarly(prefix,rd100,fk100,rd10,fk10,root);root.Close();
  std::cout<<"PASS controlled same-reservoir 100ns vs 10ns comparison for "
      <<refId<<" (not a claim of accuracy)\n"
      <<"See "<<prefix<<"_windows.csv and _metrics.csv, plus _early_overlay.png\n"
      <<"Check 100ns direct vs rebinned 10ns direct FIRST; "
      <<"then check fast first100/200ns charge error.\n";
  return 0;
 }catch(const std::exception&e){
  std::cerr<<"PLANE TIMING COMPARISON FAILED: "<<e.what()<<"\n";return 2;
 }
}
