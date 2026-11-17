// ROOT/C++ only, five-point truly uncapped Garfield++ avalanche census.
// Whole-muon bootstrap, no Gaussian fit, no incomplete-muon conditioning.
#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TRandom3.h"
#include "TSystem.h"
#include "TError.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace UncappedHV {
const std::array<int,5> voltages{{2000,2050,2100,2150,2200}};
const std::array<int,5> colors{{kBlack,kRed+1,kBlue+1,kGreen+2,kMagenta+2}};
std::vector<std::string> Split(const std::string& s){
  std::stringstream ss(s);std::string x;std::vector<std::string> out;
  while(std::getline(ss,x,',')) {
    if(!x.empty()&&x.back()=='\r')x.pop_back();
    out.push_back(x);
  }
  return out;
}
struct CSV{
  std::ifstream in;
  std::map<std::string,int> col;
  std::string path;
  explicit CSV(const std::string& p):in(p),path(p){
    if(!in)throw std::runtime_error("Cannot read "+p);
    std::string s;
    if(!std::getline(in,s))throw std::runtime_error("Empty CSV "+p);
    auto a=Split(s);
    for(int i=0;i<int(a.size());i++)col[a[i]]=i;
  }
  bool Next(std::vector<std::string>& row){
    std::string s;if(!std::getline(in,s))return false;
    row=Split(s);return true;
  }
  std::string Get(const std::vector<std::string>& row,const std::string& key)const{
    const auto i=col.find(key);
    if(i==col.end())throw std::runtime_error("Missing field "+key+" in "+path);
    if(i->second>=int(row.size()))throw std::runtime_error("Short CSV row "+path);
    return row[i->second];
  }
  double Num(const std::vector<std::string>& row,const std::string& key)const {
    const auto s=Get(row,key);
    if(s.empty())return std::numeric_limits<double>::quiet_NaN();
    return std::stod(s);
  }
};
struct Event {
  std::string id;
  int nPrimary=0;
  double Gmu=std::numeric_limits<double>::quiet_NaN();
  double A=std::numeric_limits<double>::quiet_NaN();
  std::vector<double> gain,ions,endpoints,garfieldNE,garfieldNI;
};
struct Data {
  int V=0;
  std::vector<Event> muons;
  std::vector<double> gains,positive,ions,positiveIons,endpoints,garfieldNE,garfieldNI,gmu,A;
  double zero=0,meanGainErr=0,cvGainErr=0,cvPositiveErr=0,cvMuonErr=0;
  double meanIonErr=0,cvIonErr=0,meanMuonErr=0;
};
double Mean(const std::vector<double>& x){
  long double sum=0.;int n=0;
  for(const auto v:x)if(std::isfinite(v)){sum+=v;n++;}
  return n?static_cast<double>(sum/n):0.;
}
double CV(const std::vector<double>& x){
  const double m=Mean(x);
  if(m<=0)return 0;
  long double s=0;int n=0;
  for(const double v:x)if(std::isfinite(v)){s+=(v-m)*(v-m);n++;}
  return n?std::sqrt(static_cast<double>(s/n))/m:0.;
}
double SD(const std::vector<double>& x){
  if(x.size()<2)return 0;
  const double m=Mean(x);long double s=0;
  for(double v:x)s+=(v-m)*(v-m);
  return std::sqrt(static_cast<double>(s/(x.size()-1)));
}
Data Read(const std::string& dir,int v){
  Data d;d.V=v;
  const std::string pref=dir+"/HV_"+std::to_string(v)+"/merged_";
  std::map<std::string,std::size_t> lookup;
  {
    CSV f(pref+"crossing.csv");std::vector<std::string> row;
    const std::string idKey="source|chunk|event|hit_index|track_id|station|det_id";
    while(f.Next(row)){
      Event e;
      e.id=f.Get(row,idKey);
      if(lookup.count(e.id))throw std::runtime_error("Duplicate muon ID "+e.id);
      e.nPrimary=static_cast<int>(f.Num(row,"heed_primaries"));
      e.Gmu=f.Num(row,"G_eff_per_Heed_primary");
      e.A=f.Num(row,"A_ion_hemi");
      lookup[e.id]=d.muons.size();
      if(std::isfinite(e.Gmu))d.gmu.push_back(e.Gmu);
      if(std::isfinite(e.A))d.A.push_back(e.A);
      d.muons.push_back(std::move(e));
    }
  }
  {
    CSV f(pref+"seed.csv");std::vector<std::string> row;
    while(f.Next(row)){
      const auto id=f.Get(row,"crossing_id");
      const auto it=lookup.find(id);
      if(it==lookup.end())throw std::runtime_error("Orphan avalanche "+id);
      if(f.Num(row,"successful")!=1 || f.Num(row,"capped")!=0)
        throw std::runtime_error("Failed/capped avalanche "+id);
      const double gain=f.Num(row,"collected_e");
      const double ion=f.Num(row,"ion_births");
      const double nend=f.Num(row,"electron_endpoints");
      const double ni=f.Num(row,"garfield_ni");
      const double ne=f.Num(row,"garfield_ne");
      if(!std::isfinite(gain)||!std::isfinite(ion)||!std::isfinite(nend)||
         gain<0||ion<0||nend<0)
        throw std::runtime_error("Invalid avalanche counts "+id);
      auto& e=d.muons[it->second];
      e.gain.push_back(gain);e.ions.push_back(ion);
      e.endpoints.push_back(nend);e.garfieldNE.push_back(ne);e.garfieldNI.push_back(ni);
      d.gains.push_back(gain);d.ions.push_back(ion);
      d.endpoints.push_back(nend);d.garfieldNI.push_back(ni);d.garfieldNE.push_back(ne);
      if(gain>0)d.positive.push_back(gain);else d.zero+=1;
      if(ion>0)d.positiveIons.push_back(ion);
    }
  }
  for(const auto& e:d.muons)
    if(e.nPrimary!=int(e.gain.size()) || e.ions.size()!=e.gain.size())
      throw std::runtime_error("Primary-count mismatch: "+e.id);
  if(d.muons.empty()||d.gains.empty())throw std::runtime_error("Empty voltage "+std::to_string(v));
  return d;
}
void Bootstrap(Data& d,int B=500){
  const int n=int(d.muons.size());if(n<2)return;
  TRandom3 rng(1221+d.V);
  std::vector<double> g,cvg,cvpos,cvion,mi,cvm,mg;
  for(int k=0;k<B;k++){
    std::vector<double> gains,positive,ions,gm;
    for(int j=0;j<n;j++){
      const auto& e=d.muons[rng.Integer(n)];
      gains.insert(gains.end(),e.gain.begin(),e.gain.end());
      ions.insert(ions.end(),e.ions.begin(),e.ions.end());
      for(double v:e.gain)if(v>0)positive.push_back(v);
      if(std::isfinite(e.Gmu))gm.push_back(e.Gmu);
    }
    if(!gains.empty()&&!positive.empty()&&!gm.empty()){
      g.push_back(Mean(gains));cvg.push_back(CV(gains));
      cvpos.push_back(CV(positive));
      mi.push_back(Mean(ions));cvion.push_back(CV(ions));
      mg.push_back(Mean(gm));cvm.push_back(CV(gm));
    }
  }
  d.meanGainErr=SD(g);d.cvGainErr=SD(cvg);d.cvPositiveErr=SD(cvpos);
  d.meanIonErr=SD(mi);d.cvIonErr=SD(cvion);
  d.meanMuonErr=SD(mg);d.cvMuonErr=SD(cvm);
}
double Max(const std::vector<double>& x){
  return x.empty()?0:*std::max_element(x.begin(),x.end());
}
double Tail(const std::vector<double>& x,double threshold){
  if(x.empty())return 0;
  return double(std::count_if(x.begin(),x.end(),[&](double v){return v>threshold;}))/x.size();
}
void WriteSummary(const std::array<Data,5>& data,const std::string& path){
  std::ofstream f(path+"/analysis/uncapped_avalanche_summary.csv");
  if(!f)throw std::runtime_error("Cannot write ROOT summary");
  f<<"voltage_V,muons,primaries,zero_fraction,mean_collected_gain,bootstrap_error_mean_gain,"
     "rms_over_mean_all,bootstrap_error_cv_all,mean_positive_gain,"
     "rms_over_mean_positive,bootstrap_error_cv_positive,mean_muon_effective_gain,"
     "bootstrap_error_muon_mean,rms_over_mean_muon,bootstrap_error_cv_muon,"
     "mean_ion_births,bootstrap_error_mean_ions,rms_over_mean_ion_births,"
     "bootstrap_error_cv_ions,max_ion_births,max_electron_endpoints,"
     "max_collected_gain,mean_garfield_ni,mean_garfield_ne,"
     "fraction_ion_births_gt_1e6,fraction_ion_births_gt_1e7\n";
  f<<std::setprecision(15);
  for(const auto& d:data) {
    f<<d.V<<","<<d.muons.size()<<","<<d.gains.size()<<","<<d.zero/d.gains.size()<<","
      <<Mean(d.gains)<<","<<d.meanGainErr<<","<<CV(d.gains)<<","<<d.cvGainErr<<","
      <<Mean(d.positive)<<","<<CV(d.positive)<<","<<d.cvPositiveErr<<","
      <<Mean(d.gmu)<<","<<d.meanMuonErr<<","<<CV(d.gmu)<<","<<d.cvMuonErr<<","
      <<Mean(d.ions)<<","<<d.meanIonErr<<","<<CV(d.ions)<<","<<d.cvIonErr<<","
      <<Max(d.ions)<<","<<Max(d.endpoints)<<","<<Max(d.gains)<<","
      <<Mean(d.garfieldNI)<<","<<Mean(d.garfieldNE)<<","
      <<Tail(d.ions,1e6)<<","<<Tail(d.ions,1e7)<<"\n";
  }
}
void Style(){
  gROOT->SetBatch(kTRUE);gStyle->SetOptStat(0);
  gStyle->SetCanvasColor(kWhite);gStyle->SetPadColor(kWhite);
  gStyle->SetLabelFont(42,"XYZ");gStyle->SetTitleFont(42,"XYZ");
  gStyle->SetTitleSize(.047,"XYZ");gStyle->SetLabelSize(.042,"XYZ");
  gStyle->SetPadLeftMargin(.17);gStyle->SetPadRightMargin(.05);
  gStyle->SetPadTopMargin(.08);gStyle->SetPadBottomMargin(.13);
  gStyle->SetErrorX(0);
}
void Brand(const std::string& description){
  TLatex latex;latex.SetNDC();latex.SetTextFont(42);latex.SetTextSize(.044);
  latex.DrawLatex(.19,.91,"#bf{NA60+} #it{Simulation}");
  latex.SetTextSize(.027);latex.DrawLatex(.19,.86,description.c_str());
}
void Save(TCanvas& c,const std::string& dir,const std::string& name){
  c.SaveAs((dir+"/analysis/"+name+".png").c_str());
  c.SaveAs((dir+"/analysis/"+name+".pdf").c_str());
}
void Trend(const std::array<Data,5>& data,const std::string& dir,bool ions){
  TCanvas c(ions?"c_ions_voltage":"c_gain_voltage","Uncapped gain and ions",900,640);
  TGraphErrors g(5);
  double ymax=0.;
  for(int i=0;i<5;i++){
    const auto& d=data[i];
    const double m=ions?Mean(d.ions):Mean(d.gains);
    const double err=ions?d.meanIonErr:d.meanGainErr;
    g.SetPoint(i,d.V,m);g.SetPointError(i,0.,err);
    ymax=std::max(ymax,m+err);
  }
  g.SetMarkerStyle(20);g.SetMarkerSize(1.3);g.SetLineWidth(2);
  g.SetMarkerColor(kBlue+1);g.SetLineColor(kBlue+1);
  g.SetTitle(ions?";Anode voltage [V];Mean ion births / Heed primary":
                  ";Anode voltage [V];Mean collected electrons / Heed primary");
  g.SetMinimum(0.);g.SetMaximum(ymax>0?1.28*ymax:1.);
  g.Draw("APL");g.GetXaxis()->SetLimits(1970,2230);
  Brand("P3 noFW, 30 #mum wire, Penning r=0.57; NO size cutoff");
  Save(c,dir,ions?"uncapped_mean_ion_births_vs_voltage":"uncapped_mean_gain_vs_voltage");
}
void Fluctuations(const std::array<Data,5>& data,const std::string& dir){
  TCanvas c("c_uncapped_cv","RMS/mean versus voltage",900,640);
  TGraphErrors ga(5),gp(5),gi(5);
  double high=0;
  for(int i=0;i<5;i++){
    const auto& d=data[i];ga.SetPoint(i,d.V,CV(d.gains));
    ga.SetPointError(i,0.,d.cvGainErr);
    gp.SetPoint(i,d.V,CV(d.positive));gp.SetPointError(i,0.,d.cvPositiveErr);
    gi.SetPoint(i,d.V,CV(d.ions));gi.SetPointError(i,0.,d.cvIonErr);
    high=std::max(high,std::max({CV(d.gains)+d.cvGainErr,CV(d.positive)+d.cvPositiveErr,CV(d.ions)+d.cvIonErr}));
  }
  for(auto* g:{&ga,&gp,&gi}){g->SetMarkerStyle(20);g->SetMarkerSize(1.1);g->SetLineWidth(2);}
  ga.SetLineColor(kBlack);ga.SetMarkerColor(kBlack);
  gp.SetLineColor(kBlue+1);gp.SetMarkerColor(kBlue+1);
  gi.SetLineColor(kRed+1);gi.SetMarkerColor(kRed+1);
  ga.SetTitle(";Anode voltage [V];RMS / mean (per primary)");
  ga.SetMinimum(0.);ga.SetMaximum(high>0?1.25*high:1.);ga.Draw("APL");
  ga.GetXaxis()->SetLimits(1970,2230);
  gp.Draw("PL SAME");gi.Draw("PL SAME");
  TLegend leg(.54,.68,.93,.85);leg.SetBorderSize(0);leg.SetFillStyle(0);
  leg.AddEntry(&ga,"Collected e, including zeros","lp");
  leg.AddEntry(&gp,"Collected e, positive only","lp");
  leg.AddEntry(&gi,"Ion births, including zeros","lp");leg.Draw();
  Brand("NO avalanche-size limit; bootstrap resamples whole muons");
  Save(c,dir,"uncapped_gain_fluctuations_vs_voltage");
}
void Distributions(const std::array<Data,5>& data,const std::string& dir,bool ions){
  TCanvas c(ions?"c_ion_hist":"c_gain_hist","Distributions",900,640);
  std::array<std::unique_ptr<TH1D>,5> hist;
  TLegend leg(.69,.66,.9,.88);leg.SetFillStyle(0);leg.SetBorderSize(0);
  double high=0.;
  for(int i=0;i<5;i++){
    const auto& vec=ions?data[i].positiveIons:data[i].positive;
    const std::string name=(ions?"ions_":"gain_")+std::to_string(data[i].V);
    hist[i]=std::make_unique<TH1D>(name.c_str(),"",180,0.,9.);
    hist[i]->SetDirectory(nullptr);hist[i]->SetLineColor(colors[i]);hist[i]->SetLineWidth(2);
    for(double x:vec)if(x>0)hist[i]->Fill(std::log10(x));
    if(hist[i]->Integral("width")>0)hist[i]->Scale(1./hist[i]->Integral("width"));
    high=std::max(high,hist[i]->GetMaximum());
    leg.AddEntry(hist[i].get(),(std::to_string(data[i].V)+" V").c_str(),"l");
  }
  hist[0]->SetTitle(ions?";log_{10}(ionization births per primary);Probability density":
                          ";log_{10}(collected electrons per primary);Probability density");
  hist[0]->SetMaximum(high>0?1.25*high:1.);hist[0]->Draw("HIST");
  for(int i=1;i<5;i++)hist[i]->Draw("HIST SAME");
  leg.Draw();Brand("Full Garfield++ avalanches; nonzero values; NO internal size limit");
  Save(c,dir,ions?"uncapped_ion_birth_distributions":"uncapped_collected_gain_distributions");
}
} // namespace UncappedHV

void AnalyzeUncappedHV_ROOT(const char* folder){
  try {
    using namespace UncappedHV;
    Style();const std::string dir(folder);
    gSystem->mkdir((dir+"/analysis").c_str(),kTRUE);
    std::array<Data,5> data;
    for(int i=0;i<5;i++){
      data[i]=Read(dir,voltages[i]);Bootstrap(data[i]);
      std::cout<<"UNCAPPED HV "<<data[i].V<<"V: "
               <<data[i].muons.size()<<" muons, "<<data[i].gains.size()
               <<" primaries; meanGain="<<Mean(data[i].gains)
               <<", maxIonBirths="<<Max(data[i].ions)
               <<", maxEndpoints="<<Max(data[i].endpoints)
               <<", RMS/mean="<<CV(data[i].gains)<<"\n";
    }
    WriteSummary(data,dir);
    Trend(data,dir,false);Trend(data,dir,true);
    Fluctuations(data,dir);
    Distributions(data,dir,false);Distributions(data,dir,true);
    std::cout<<"COMPLETE ROOT uncapped analysis: "<<dir<<"/analysis/\n";
  }catch(const std::exception& e){
    ::Error("AnalyzeUncappedHV_ROOT","%s",e.what());
    throw;
  }
}
