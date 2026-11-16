// ROOT/C++ ONLY. Garfield++ voltage and wire-diameter validation.
// Uses complete muons as independent units for bootstrap uncertainty.
// No Gaussian gain fits. ATLAS-inspired plots marked "NA60+ Simulation".
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TRandom3.h"
#include "TError.h"
#include "TString.h"
#include <algorithm>
#include <array>
#include <cmath>
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

namespace HvDia {
struct Reader {
  std::ifstream is;
  std::map<std::string,int> key;
  std::string source;
  explicit Reader(std::string name):is(name),source(name){
    if(!is)throw std::runtime_error("No CSV: "+source);
    std::string line;if(!std::getline(is,line))
      throw std::runtime_error("Empty CSV: "+source);
    auto items=Split(line);
    for(std::size_t i=0;i<items.size();++i)key[items[i]]=int(i);
  }
  static std::vector<std::string> Split(const std::string& line){
    std::stringstream ss(line);std::string x;std::vector<std::string> out;
    while(std::getline(ss,x,',')){
      if(!x.empty()&&x.back()=='\r')x.pop_back();
      out.push_back(x);
    }
    return out;
  }
  bool Next(std::vector<std::string>& out){
    std::string s;if(!std::getline(is,s))return false;
    out=Split(s);return true;
  }
  std::string Get(const std::vector<std::string>& row,const std::string& k)const{
    auto it=key.find(k);
    if(it==key.end())throw std::runtime_error(source+": absent column "+k);
    if(it->second>=int(row.size()))throw std::runtime_error("short row in "+source);
    return row[it->second];
  }
  double Num(const std::vector<std::string>& row,const std::string& k)const{
    return std::stod(Get(row,k));
  }
};
struct Moments{
  double sum=0.,sum2=0.;
  long long n=0;
  void Add(double x){if(!std::isfinite(x))return;sum+=x;sum2+=x*x;n++;}
  void Merge(const Moments& b){sum+=b.sum;sum2+=b.sum2;n+=b.n;}
  double Mean()const{return n?sum/n:0.;}
  double RMS()const{return n?std::sqrt(std::max(0.,sum2/n-Mean()*Mean())):0.;}
  double CV()const{return Mean()>0?RMS()/Mean():0.;}
};
struct Muon{
  std::string id;
  double gmu=std::numeric_limits<double>::quiet_NaN();
  double A=std::numeric_limits<double>::quiet_NaN();
  int nHeed=0;
  Moments all,positive;
};
struct Sample{
  int voltage=0,diameter=0;
  std::vector<Muon> muons;
  std::vector<double> allG,positiveG,muG,As;
  Moments all,positive,muon;
  double zero=0.;
  double eMeanAll=0.,eCvAll=0.,eCvPositive=0.,eCvMuon=0.,eMeanMuon=0.;
  double eZero=0.;
};
double SD(const std::vector<double>& v){
  if(v.size()<2)return 0.;
  double s=0.;for(double x:v)s+=x;
  const double m=s/v.size();double sq=0.;
  for(double x:v)sq+=(x-m)*(x-m);
  return std::sqrt(sq/(v.size()-1));
}
Sample Read(const std::string& dir,int voltage,int diameter,bool firstFive=false) {
  Sample o;o.voltage=voltage;o.diameter=diameter;
  std::string prefix=dir+"/V_"+std::to_string(voltage)+"_D"+
    std::to_string(diameter)+"/";
  if(firstFive)prefix+="chunk_0000/micro_";
  else prefix+="merged_";
  std::map<std::string,int> id;
  {
    Reader f(prefix+"crossing.csv");
    std::vector<std::string> row;
    const std::string crossingKey="source|chunk|event|hit_index|track_id|station|det_id";
    while(f.Next(row)){
      Muon m;
      m.id=f.Get(row,crossingKey);
      m.gmu=f.Num(row,"G_eff_per_Heed_primary");
      m.A=f.Num(row,"A_ion_hemi");
      m.nHeed=int(f.Num(row,"heed_primaries"));
      if(id.count(m.id))throw std::runtime_error("Duplicated crossing "+m.id);
      id[m.id]=int(o.muons.size());
      o.muons.push_back(m);
      if(std::isfinite(m.gmu)){o.muon.Add(m.gmu);o.muG.push_back(m.gmu);}
      if(std::isfinite(m.A))o.As.push_back(m.A);
    }
  }
  {
    Reader f(prefix+"seed.csv");
    std::vector<std::string> row;
    while(f.Next(row)){
      std::string key=f.Get(row,"crossing_id");
      auto it=id.find(key);
      if(it==id.end())throw std::runtime_error("Orphan seed row "+key);
      if(f.Num(row,"capped")!=0 || f.Num(row,"successful")!=1)
        throw std::runtime_error("Capped/failed avalanche in "+prefix);
      const double gain=f.Num(row,"collected_e");
      if(!std::isfinite(gain)||gain<0)throw std::runtime_error("Bad gain");
      auto& m=o.muons[it->second];
      m.all.Add(gain);o.all.Add(gain);o.allG.push_back(gain);
      if(gain>0){
        m.positive.Add(gain);o.positive.Add(gain);o.positiveG.push_back(gain);
      }else o.zero+=1.;
    }
  }
  for(const auto& m:o.muons)
    if(m.all.n!=m.nHeed)
      throw std::runtime_error("Mismatch Heed primary vs seeds: "+m.id);
  if(o.muons.empty())throw std::runtime_error("No muons");
  return o;
}
void Bootstrap(Sample& s,int B=500){
  if(s.muons.size()<2)return;
  TRandom3 rng(1923+s.voltage+100*s.diameter);
  std::vector<double> mAll,cvAll,cvPos,cvMuon,mMuon,zf;
  for(int b=0;b<B;b++){
    Moments all,pos,mu;
    long long zero=0;
    for(std::size_t j=0;j<s.muons.size();j++){
      const auto& event=s.muons[rng.Integer(UInt_t(s.muons.size()))];
      all.Merge(event.all);
      pos.Merge(event.positive);
      mu.Add(event.gmu);
      zero+=event.all.n-event.positive.n;
    }
    if(all.n>0&&pos.n>0&&mu.n>0){
      mAll.push_back(all.Mean());cvAll.push_back(all.CV());
      cvPos.push_back(pos.CV());cvMuon.push_back(mu.CV());
      mMuon.push_back(mu.Mean());zf.push_back(double(zero)/all.n);
    }
  }
  s.eMeanAll=SD(mAll);s.eCvAll=SD(cvAll);s.eCvPositive=SD(cvPos);
  s.eCvMuon=SD(cvMuon);s.eMeanMuon=SD(mMuon);s.eZero=SD(zf);
}
void Style(){
  gROOT->SetBatch(kTRUE);gStyle->SetOptStat(0);
  gStyle->SetFrameLineWidth(1);gStyle->SetLineWidth(2);
  gStyle->SetLabelFont(42,"XYZ");gStyle->SetTitleFont(42,"XYZ");
  gStyle->SetLabelSize(.04,"XYZ");gStyle->SetTitleSize(.046,"XYZ");
  gStyle->SetPadLeftMargin(.16);gStyle->SetPadBottomMargin(.14);
  gStyle->SetPadTopMargin(.07);gStyle->SetPadRightMargin(.05);
  gStyle->SetErrorX(0);
}
void Brand(const char* desc){
  TLatex x;x.SetNDC();x.SetTextFont(42);x.SetTextSize(.048);
  x.DrawLatex(.21,.90,"#bf{NA60+} #it{Simulation}");
  x.SetTextSize(.030);x.DrawLatex(.21,.845,desc);
}
void Save(TCanvas& c,const std::string& folder,const char* name){
  c.SaveAs((folder+"/"+name+".pdf").c_str());
  c.SaveAs((folder+"/"+name+".png").c_str());
}
void Summary(const std::array<Sample,3>& s,const Sample& wire50,const Sample& wire30_5,const std::string& out){
  std::ofstream f(out+"/analysis/gain_summary.csv");
  f<<"voltage_V,anode_diameter_um,n_muons,n_primary,mean_gain_all,bootstrap_error_mean,"
   <<"rms_over_mean_all,bootstrap_error_cv_all,mean_gain_positive,rms_over_mean_positive,"
   <<"bootstrap_error_cv_pos,zero_fraction,bootstrap_error_zero,mean_muon_gain,"
   <<"bootstrap_error_muon_mean,rms_over_mean_muon,bootstrap_error_cv_muon\n";
  f<<std::setprecision(12);
  const auto row=[&](const Sample& a){
    f<<a.voltage<<","<<a.diameter<<","<<a.muons.size()<<","<<a.all.n<<","
     <<a.all.Mean()<<","<<a.eMeanAll<<","<<a.all.CV()<<","<<a.eCvAll<<","
     <<a.positive.Mean()<<","<<a.positive.CV()<<","<<a.eCvPositive<<","
     <<a.zero/a.all.n<<","<<a.eZero<<","<<a.muon.Mean()<<","
     <<a.eMeanMuon<<","<<a.muon.CV()<<","<<a.eCvMuon<<"\n";
  };
  for(const auto& a:s)row(a);
  row(wire30_5);row(wire50);
}
void GainVoltage(const std::array<Sample,3>& a,const std::string& out){
  TCanvas c("c_gain","Gain voltage",900,650);
  TGraphErrors g(3);
  double high=0.;
  for(int j=0;j<3;j++){
    g.SetPoint(j,a[j].voltage,a[j].all.Mean());
    g.SetPointError(j,0.,a[j].eMeanAll);
    high=std::max(high,a[j].all.Mean()+a[j].eMeanAll);
  }
  g.SetTitle(";Anode high voltage [V];Mean collected gain [electrons / primary]");
  g.SetMarkerStyle(20);g.SetMarkerSize(1.25);g.SetLineWidth(2);
  g.SetMarkerColor(kBlack);g.SetLineColor(kBlack);
  g.SetMinimum(0.);g.SetMaximum(high*1.3);g.Draw("APL");
  g.GetXaxis()->SetLimits(1970,2130);
  Brand("Anode diameter = 30 #mum; all Heed primaries incl. zero response");
  Save(c,out+"/analysis","mean_gain_vs_voltage_30um");
}
void CVVoltage(const std::array<Sample,3>& a,const std::string& out){
  TCanvas c("c_cv","RMS / mean",900,650);
  TGraphErrors gAll(3),gPos(3),gMuon(3);
  double hi=0.;
  for(int i=0;i<3;i++){
    const double v=a[i].voltage;
    gAll.SetPoint(i,v,a[i].all.CV());gAll.SetPointError(i,0,a[i].eCvAll);
    gPos.SetPoint(i,v,a[i].positive.CV());gPos.SetPointError(i,0,a[i].eCvPositive);
    gMuon.SetPoint(i,v,a[i].muon.CV());gMuon.SetPointError(i,0,a[i].eCvMuon);
    hi=std::max(hi,std::max({a[i].all.CV()+a[i].eCvAll,a[i].positive.CV()+a[i].eCvPositive,a[i].muon.CV()+a[i].eCvMuon}));
  }
  gAll.SetTitle(";Anode high voltage [V];RMS / mean");
  gAll.SetLineColor(kBlack);gAll.SetMarkerColor(kBlack);
  gPos.SetLineColor(kBlue+1);gPos.SetMarkerColor(kBlue+1);
  gMuon.SetLineColor(kRed+1);gMuon.SetMarkerColor(kRed+1);
  for(auto* g:{&gAll,&gPos,&gMuon}){
    g->SetLineWidth(2);g->SetMarkerStyle(20);g->SetMarkerSize(1.25);
  }
  gAll.SetMinimum(0);gAll.SetMaximum(hi>0?hi*1.24:1.);
  gAll.Draw("APL");gAll.GetXaxis()->SetLimits(1970,2130);
  gPos.Draw("PL SAME");gMuon.Draw("PL SAME");
  TLegend leg(.48,.68,.92,.85);leg.SetBorderSize(0);leg.SetFillStyle(0);
  leg.AddEntry(&gAll,"Per primary (all, including zero)","lp");
  leg.AddEntry(&gPos,"Per primary (positive only)","lp");
  leg.AddEntry(&gMuon,"Per muon effective gain","lp");leg.Draw();
  Brand("30 #mum anode; 500 muon-level bootstrap resamples");
  Save(c,out+"/analysis","gain_fluctuations_vs_voltage_30um");
}
void GainDistributions(const std::array<Sample,3>& a,const std::string& out){
  TCanvas c("c_dist","Gain distributions",900,650);
  TLegend leg(.68,.70,.90,.88);leg.SetBorderSize(0);leg.SetFillStyle(0);
  const std::array<int,3> colors{{kBlack,kBlue+1,kRed+1}};
  std::array<std::unique_ptr<TH1D>,3> h;
  double highest=0.;
  for(int i=0;i<3;i++){
    const std::string name="positive_gain_"+std::to_string(i);
    h[i]=std::make_unique<TH1D>(name.c_str(),"",100,0,6.2);
    h[i]->SetDirectory(nullptr);h[i]->SetLineColor(colors[i]);h[i]->SetLineWidth(2);
    for(double g:a[i].positiveG)if(g>0)h[i]->Fill(std::log10(g));
    if(h[i]->Integral("width")>0)h[i]->Scale(1./h[i]->Integral("width"));
    highest=std::max(highest,h[i]->GetMaximum());
    leg.AddEntry(h[i].get(),(std::to_string(a[i].voltage)+" V").c_str(),"l");
  }
  h[0]->SetTitle(";log_{10}(collected electrons per primary);Probability density");
  h[0]->SetMaximum(highest*1.28);h[0]->Draw("HIST");
  for(int i=1;i<3;i++)h[i]->Draw("HIST SAME");
  leg.Draw();Brand("30 #mum anode; positive gain only; P(g=0) separate");
  Save(c,out+"/analysis","positive_gain_distributions_30um");
}
void CompareDiameter(const Sample& d30,const Sample& d50,const std::string& out){
  TCanvas c("c_diam","Diameter control, 5 common muons",900,650);
  TGraphErrors g(2);
  g.SetPoint(0,30,d30.all.Mean());g.SetPointError(0,0,d30.eMeanAll);
  g.SetPoint(1,50,d50.all.Mean());g.SetPointError(1,0,d50.eMeanAll);
  g.SetTitle(";Anode wire diameter [#mum];Mean collected gain [electrons / primary]");
  g.SetMarkerStyle(20);g.SetMarkerSize(1.45);g.SetLineWidth(2);
  g.SetMarkerColor(kBlue+1);g.SetLineColor(kBlue+1);
  const double maxY=std::max(d30.all.Mean()+d30.eMeanAll,d50.all.Mean()+d50.eMeanAll);
  g.SetMinimum(0);g.SetMaximum(maxY*1.32);g.Draw("APL");
  g.GetXaxis()->SetLimits(25,55);
  Brand("2000 V; same 5 muons and Heed seeds; exploratory");
  Save(c,out+"/analysis","diameter_2000V_5muon_comparison");
}
} // namespace HvDia

void Analyze_HVdiam_ROOT(const char* workDir) {
  using namespace HvDia;
  try{
    Style();
    const std::string out=workDir;
    gSystem->mkdir((out+"/analysis").c_str(),true);
    std::array<Sample,3> s{{
      Read(out,2000,30),Read(out,2050,30),Read(out,2100,30)
    }};
    for(auto& a:s)Bootstrap(a);
    auto d50=Read(out,2000,50);Bootstrap(d50);
    auto d30_5=Read(out,2000,30,true);Bootstrap(d30_5);
    Summary(s,d50,d30_5,out);
    GainVoltage(s,out);CVVoltage(s,out);GainDistributions(s,out);
    CompareDiameter(d30_5,d50,out);
    for(const auto& a:s)
      std::cout<<"HV "<<a.voltage<<" V, D "<<a.diameter<<" um: n_muons="
               <<a.muons.size()<<" n_primaries="<<a.all.n
               <<" gain="<<a.all.Mean()<<" RMS/mean="<<a.all.CV()<<std::endl;
    std::cout<<"At 2000V, SAME first five muons: diameter30 gain="
             <<d30_5.all.Mean()<<" diameter50 gain="<<d50.all.Mean()<<"\n";
    std::cout<<"ROOT results: "<<out<<"/analysis\n";
  }catch(const std::exception& e){
    ::Error("Analyze_HVdiam_ROOT","%s",e.what());
    throw;
  }
}
