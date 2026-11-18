// NA60+ Stage C: assemble complete muon e+ion currents on two WHOLE planes.
// Electron current: direct Garfield replay of ALL Heed primaries, 0.05 ns bins.
// Ion current: original SAME-event K512 100-ns direct RKF and kernel waveforms.
// No electronics shaper is applied. 100ns bins are NOT enough to resolve the
// true VMM3a peaking-time response, especially its earliest current samples.
#include <TFile.h>
#include <TTree.h>
#include <TNamed.h>
#include <TH1D.h>
#include <TCanvas.h>
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
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
namespace {
void Must(bool b,const std::string&msg) {
 if(!b)throw std::runtime_error("MUON e+ion ASSEMBLY: "+msg);
}
std::string Arg(int argc,char**argv,const std::string&key){
 for(int i=1;i+1<argc;++i)if(key==argv[i])return argv[i+1];
 return "";
}
std::vector<std::string> Split(const std::string&s) {
 std::vector<std::string> v;std::istringstream ss(s);std::string x;
 while(std::getline(ss,x,','))v.push_back(x);
 return v;
}
struct Meta {
 std::string id;
 int station=0,det=0,primaries=0;
 long long nIon=0,nElectron=0;
 double dt=0,window=0,maxEndpoint=0;
};
struct Sample {double t=0,ekp=0,ekm=0,dp=0,dm=0,kp=0,km=0;};
Meta ReadElectrons(const std::string&path,std::vector<Sample>&bins,double dt) {
 std::unique_ptr<TFile> root(TFile::Open(path.c_str(),"READ"));
 Must(root&&!root->IsZombie(),"cannot open complete electron response "+path);
 auto* meta=dynamic_cast<TTree*>(root->Get("muon"));
 auto* sig=dynamic_cast<TTree*>(root->Get("electron_signal"));
 Must(meta&&sig&&meta->GetEntries()==1,
      "direct electron ROOT is missing audited whole muon");
 Meta m;
 std::string* id=nullptr;
 meta->SetBranchAddress("crossing_id",&id);
 meta->SetBranchAddress("station",&m.station);
 meta->SetBranchAddress("det_id",&m.det);
 meta->SetBranchAddress("heed_primaries",&m.primaries);
 meta->SetBranchAddress("exact_ion_births",&m.nIon);
 meta->SetBranchAddress("collected_electrons",&m.nElectron);
 meta->SetBranchAddress("dt_ns",&m.dt);
 meta->SetBranchAddress("time_window_ns",&m.window);
 meta->SetBranchAddress("max_electron_endpoint_ns",&m.maxEndpoint);
 Must(meta->GetEntry(0)>0&&id,"incomplete electron metadata");
 m.id=*id;
 Must(m.primaries>0&&m.nIon>0&&m.nElectron>0 &&
      std::abs(m.dt-0.05)<1.e-10 && m.window==200.,
      "invalid full-electron metadata or unexpected fine time grid");
 int timeBin=0;double t=0,plus=0,minus=0;
 sig->SetBranchAddress("time_bin",&timeBin);
 sig->SetBranchAddress("time_center_ns",&t);
 sig->SetBranchAddress("plus_electron_fC_per_ns",&plus);
 sig->SetBranchAddress("minus_electron_fC_per_ns",&minus);
 Must(sig->GetEntries()==4000,"direct-electron 200ns signal is incomplete");
 for(Long64_t i=0;i<sig->GetEntries();++i) {
   Must(sig->GetEntry(i)>0&&timeBin==i &&
        std::abs(t-(i+0.5)*m.dt)<1.e-8 &&
        std::isfinite(plus)&&std::isfinite(minus),
        "malformed fine-bin direct-electron waveform");
   const int coarse=(int)std::floor(t/dt);
   Must(coarse>=0&&coarse<(int)bins.size(),
        "electron waveform extends beyond coarse ion window");
   bins[coarse].ekp+=plus*m.dt/dt;
   bins[coarse].ekm+=minus*m.dt/dt;
 }
 return m;
}
void VerifyIonIdentity(const std::string&path,const Meta&m,int expectedBins,double dt) {
 std::ifstream in(path);Must(bool(in),"cannot read ion-only summary "+path);
 std::string header,record;
 Must(bool(std::getline(in,header))&&bool(std::getline(in,record)),
      "empty ion-only muon summary");
 auto keys=Split(header),vals=Split(record);
 Must(keys.size()==vals.size(),"ion summary columns mismatch");
 auto get=[&](const std::string&key){
   auto it=std::find(keys.begin(),keys.end(),key);
   Must(it!=keys.end(),"required ion summary key missing: "+key);
   return vals[std::distance(keys.begin(),it)];
 };
 Must(get("crossing_id")==m.id &&
      std::stoi(get("station"))==m.station &&
      std::stoi(get("det_id"))==m.det &&
      std::stoll(get("true_ion_births"))==m.nIon &&
      std::stoll(get("retained_macroions"))>0 &&
      std::abs(std::stod(get("dt_ns"))-dt)<1.e-9*dt &&
      std::abs(std::stod(get("window_us"))*1000.-expectedBins*dt)<1.e-7,
      "ion/electron results refer to different crossing, population or time grid");
}
void ReadIons(const std::string&path,std::vector<Sample>&data,double dt){
 std::ifstream in(path);Must(bool(in),"cannot read ion-only event CSV "+path);
 std::string header,line;
 Must(bool(std::getline(in,header)),"empty ion-only event CSV");
 const std::string expected=
   "time_center_ns,plus_kernel_fC_per_ns,minus_kernel_fC_per_ns,"
   "plus_direct_fC_per_ns,minus_direct_fC_per_ns,"
   "plus_kernel_cumulative_fC,minus_kernel_cumulative_fC,"
   "plus_direct_cumulative_fC,minus_direct_cumulative_fC";
 Must(header==expected,"unexpected ion event waveform schema");
 for(int i=0;i<(int)data.size();++i){
   Must(bool(std::getline(in,line)),"ion waveform truncated");
   const auto cols=Split(line);Must(cols.size()==9,"ion waveform has wrong column count");
   auto&x=data[i];
   x.t=std::stod(cols[0]);x.kp=std::stod(cols[1]);x.km=std::stod(cols[2]);
   x.dp=std::stod(cols[3]);x.dm=std::stod(cols[4]);
   Must(std::abs(x.t-(i+0.5)*dt)<1.e-7*dt &&
      std::isfinite(x.dp)&&std::isfinite(x.dm) &&
      std::isfinite(x.kp)&&std::isfinite(x.km),
      "ion waveform has nonfinite currents or mismatched grid");
 }
 Must(!std::getline(in,line),"ion waveform has extra bins");
}
double Charge(const std::vector<Sample>&rows,double dt,int channel,int n) {
 double total=0;
 for(int j=0;j<n;++j){
   const auto&r=rows[j];
   if(channel==0)total+=dt*r.ekp;
   else if(channel==1)total+=dt*r.ekm;
   else if(channel==2)total+=dt*r.dp;
   else if(channel==3)total+=dt*r.dm;
   else if(channel==4)total+=dt*r.kp;
   else total+=dt*r.km;
 }
 return total;
}
void MakePlot(const std::vector<Sample>&data,double dt,
              const std::string&prefix,TFile&out) {
 gROOT->SetBatch(kTRUE);gStyle->SetOptStat(0);
 TCanvas canvas("electron_plus_ion_planes",
                "NA60+ complete muon electron and ion cathode currents",1350,940);
 canvas.Divide(2,2);
 // ROOT pads hold pointers to objects passed to TH1::Draw().
 // Prior stack-local histograms were destroyed before canvas.SaveAs(),
 // leaving only cloned legends. Keep all plot histograms alive until print.
 std::vector<std::unique_ptr<TH1D>> drawn;
 drawn.reserve(16);
 for(int c=0;c<2;++c)for(int zoom=0;zoom<2;++zoom){
   const int n=zoom?std::min((int)data.size(),20):(int)data.size();
   const std::string suffix=std::string(c?"minus":"plus")+
                            (zoom?"_first2us":"_full");
   const std::string caption=std::string(c?"Lower":"Upper")+
             " whole-plane induced current"+
             (zoom?" (first 2 #mus)":" (100 #mus)");
   auto electron=std::make_unique<TH1D>(("electron_"+suffix).c_str(),
             (caption+";t (ns);current (fC/ns)").c_str(),n,0,n*dt);
   auto ion=std::make_unique<TH1D>(("ion_"+suffix).c_str(),"",n,0,n*dt);
   auto sum=std::make_unique<TH1D>(("total_"+suffix).c_str(),
             (caption+";t (ns);current (fC/ns)").c_str(),n,0,n*dt);
   auto sumKernel=std::make_unique<TH1D>(("kernel_total_"+suffix).c_str(),"",n,0,n*dt);
   for(TH1D* h:{electron.get(),ion.get(),sum.get(),sumKernel.get()})
     h->SetDirectory(nullptr);
   for(int j=0;j<n;++j){
     const auto&p=data[j];
     const double e=c?p.ekm:p.ekp;
     const double i=c?p.dm:p.dp;
     const double k=c?p.km:p.kp;
     electron->SetBinContent(j+1,e);
     ion->SetBinContent(j+1,i);
     sum->SetBinContent(j+1,e+i);
     sumKernel->SetBinContent(j+1,e+k);
   }
   electron->SetLineColor(kRed+1);electron->SetLineWidth(2);
   ion->SetLineColor(kBlue+1);ion->SetLineWidth(2);
   sum->SetLineColor(kBlack);sum->SetLineWidth(2);
   sumKernel->SetLineColor(kGreen+2);sumKernel->SetLineWidth(2);
   sumKernel->SetLineStyle(2);
   const double low=std::min({0.,electron->GetMinimum(),
                              ion->GetMinimum(),sum->GetMinimum(),
                              sumKernel->GetMinimum()});
   const double high=std::max({0.,electron->GetMaximum(),
                               ion->GetMaximum(),sum->GetMaximum(),
                               sumKernel->GetMaximum()});
   Must(std::isfinite(low)&&std::isfinite(high)&&high>low,
        "plot has empty or invalid signal ranges");
   sum->SetMinimum(low-.10*std::max(1.e-10,high-low));
   sum->SetMaximum(high+.16*std::max(1.e-10,high-low));
   canvas.cd(zoom*2+c+1);
   sum->Draw("HIST");
   ion->Draw("HIST SAME");
   electron->Draw("HIST SAME");
   sumKernel->Draw("HIST SAME");
   TLegend leg(.43,.69,.91,.89);
   leg.SetFillStyle(0);leg.SetBorderSize(0);leg.SetTextSize(.032);
   leg.AddEntry(sum.get(),"Total: electron + direct RKF ion","l");
   leg.AddEntry(sumKernel.get(),"Total: electron + fast ion kernel","l");
   leg.AddEntry(electron.get(),"Microscopic prompt electrons","l");
   leg.AddEntry(ion.get(),"Direct weighted ion response","l");
   leg.DrawClone();
   out.cd();electron->Write();ion->Write();sum->Write();sumKernel->Write();
   drawn.push_back(std::move(electron));
   drawn.push_back(std::move(ion));
   drawn.push_back(std::move(sum));
   drawn.push_back(std::move(sumKernel));
 }
 canvas.cd(0);
 canvas.Modified();
 canvas.Update();
 out.cd();canvas.Write();
 canvas.SaveAs((prefix+"_electron_ion_planes.png").c_str());
}
}
int main(int argc,char**argv){
 try{
  const std::string electron=Arg(argc,argv,"--electron-root");
  const std::string ion=Arg(argc,argv,"--ion-waveforms");
  const std::string ionSummary=Arg(argc,argv,"--ion-summary");
  const std::string prefix=Arg(argc,argv,"--out-prefix");
  Must(!electron.empty()&&!ion.empty()&&!ionSummary.empty()&&!prefix.empty(),
    "usage: --electron-root FILE --ion-waveforms CSV --ion-summary CSV --out-prefix PREFIX");
  constexpr double dt=100.;
  constexpr int nBins=1000;
  for(const std::string&ext:{"_summary.csv","_waveforms.csv",
                              "_response.root","_electron_ion_planes.png"})
   Must(!std::filesystem::exists(prefix+ext),"output exists: "+prefix+ext);
  std::vector<Sample> rows(nBins);
  const auto meta=ReadElectrons(electron,rows,dt);
  VerifyIonIdentity(ionSummary,meta,nBins,dt);
  ReadIons(ion,rows,dt);
  double eP=Charge(rows,dt,0,nBins),eM=Charge(rows,dt,1,nBins);
  double ionP=Charge(rows,dt,2,nBins),ionM=Charge(rows,dt,3,nBins);
  double kP=Charge(rows,dt,4,nBins),kM=Charge(rows,dt,5,nBins);
  const double first100ePlus=Charge(rows,dt,0,1);
  const double first100eMinus=Charge(rows,dt,1,1);
  const double first100ionPlus=Charge(rows,dt,2,1);
  const double first100ionMinus=Charge(rows,dt,3,1);
  std::ofstream csv(prefix+"_waveforms.csv");
  Must(bool(csv),"cannot save combined waveform CSV");
  csv<<std::setprecision(16)
     <<"time_center_ns,plus_direct_ion_fC_per_ns,minus_direct_ion_fC_per_ns,"
       "plus_electron_fC_per_ns,minus_electron_fC_per_ns,"
       "plus_total_direct_fC_per_ns,minus_total_direct_fC_per_ns,"
       "plus_total_kernel_fC_per_ns,minus_total_kernel_fC_per_ns\n";
  for(auto&p:rows)
   csv<<p.t<<","<<p.dp<<","<<p.dm<<","<<p.ekp<<","<<p.ekm<<","
      <<p.dp+p.ekp<<","<<p.dm+p.ekm<<","
      <<p.kp+p.ekp<<","<<p.km+p.ekm<<"\n";
  csv.close();
  std::ofstream summary(prefix+"_summary.csv");
  Must(bool(summary),"cannot save summary");
  summary<<std::setprecision(16)
    <<"crossing_id,station,det_id,Heed_primaries,true_ion_births,"
    "collected_electrons,plus_electron_Q_fC,minus_electron_Q_fC,"
    "plus_direct_ion_Q_fC,minus_direct_ion_Q_fC,"
    "plus_total_direct_Q_fC,minus_total_direct_Q_fC,"
    "plus_total_kernel_Q_fC,minus_total_kernel_Q_fC,"
    "plus_first100ns_electron_Q_fC,minus_first100ns_electron_Q_fC,"
    "plus_first100ns_direct_ion_Q_fC,minus_first100ns_direct_ion_Q_fC,"
    "electron_max_endpoint_ns,coarse_dt_ns\n";
  summary<<meta.id<<","<<meta.station<<","<<meta.det<<","<<meta.primaries
    <<","<<meta.nIon<<","<<meta.nElectron<<","
    <<eP<<","<<eM<<","<<ionP<<","<<ionM<<","
    <<eP+ionP<<","<<eM+ionM<<","<<eP+kP<<","<<eM+kM<<","
    <<first100ePlus<<","<<first100eMinus<<","
    <<first100ionPlus<<","<<first100ionMinus<<","
    <<meta.maxEndpoint<<","<<dt<<"\n";
  summary.close();
  TFile out((prefix+"_response.root").c_str(),"RECREATE",
            "NA60+ complete muon e+ion whole-plane signal",4);
  Must(!out.IsZombie(),"cannot save ROOT response");
  TNamed provenance("source_paths",
    ("direct_electron_ROOT="+electron+"\nion_only_CSV="+ion+
     "\nion_only_summary="+ionSummary).c_str());
  provenance.Write();
  TNamed scope("scope",
    "Two WHOLE cathode planes, NOT stereo strips. Complete muon: direct "
    "Garfield microscopic electrons, K512 weighted RKF ions, fast-ion kernel "
    "alternative. 100ns ion waveform bins; 0.05ns electron bins averaged to "
    "100ns for comparison. NO shaper/noise/threshold. 100ns is NOT a "
    "VMM3a box integration window; early ion waveform needs finer timing "
    "before 25/50/100/200ns frontend peaking can be simulated.");
  scope.Write();
  MakePlot(rows,dt,prefix,out);
  out.Close();
  std::cout<<std::setprecision(10)
    <<"PASS WHOLE-MUON direct microscopic electron + weighted ion plane currents\n"
    <<"crossing="<<meta.id<<" primaries="<<meta.primaries
    <<" electron_end_ns_max="<<meta.maxEndpoint<<"\n"
    <<"plus_Q_e_fC="<<eP<<" plus_Q_ion_fC="<<ionP
    <<" total_Q_fC="<<eP+ionP<<"\n"
    <<"minus_Q_e_fC="<<eM<<" minus_Q_ion_fC="<<ionM
    <<" total_Q_fC="<<eM+ionM<<"\n"
    <<"first100ns_plus_Q_e_fC="<<first100ePlus
    <<" first100ns_plus_Q_ion_fC="<<first100ionPlus<<"\n"
    <<"first100ns_minus_Q_e_fC="<<first100eMinus
    <<" first100ns_minus_Q_ion_fC="<<first100ionMinus<<"\n"
    <<"No VMM shaping modeled; early 100ns ion error persists.\n";
  return 0;
 } catch(const std::exception&e) {
  std::cerr<<"ELECTRON+ION ASSEMBLY FAILED: "<<e.what()<<"\n";
  return 2;
 }
}
