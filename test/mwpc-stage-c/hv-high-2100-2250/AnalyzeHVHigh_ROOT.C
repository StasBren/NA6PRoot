// 2100-2250V experimental HIGH-HV extension. ROOT/C++ only.
// Interpret results ONLY if cap / numerical / space-charge diagnostics pass.
// ROOT/C++ only: mean gain, microscopic gain fluctuations, muon-effective
// gain fluctuations, and signed asymmetry versus voltage.
// All statistical error estimates use bootstrap resampling WHOLE muons.
// No Gaussian fits; all seed-level and event-level distributions preserved.
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TH1D.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TRandom3.h"
#include "TROOT.h"
#include "TAxis.h"
#include "TMath.h"
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
#include <numeric>
#include <stdexcept>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

namespace HvAna {
constexpr std::array<int,4> VS{{2100,2150,2200,2250}};
constexpr std::array<int,4> COL{{kBlack,kRed+1,kBlue+1,kGreen+2}};

std::vector<std::string> Cells(const std::string& line) {
  std::vector<std::string> a;std::stringstream s(line);std::string x;
  while(std::getline(s,x,',')) a.push_back(x);
  if(!line.empty()&&line.back()==',')a.push_back("");
  return a;
}
struct CSV {
  std::ifstream in;
  std::unordered_map<std::string,std::size_t> col;
  std::string filename;
  explicit CSV(const std::string& path):in(path),filename(path) {
    if(!in)throw std::runtime_error("Unable to read "+path);
    std::string head;
    if(!std::getline(in,head))throw std::runtime_error("Empty CSV "+path);
    if(!head.empty()&&head.back()=='\r')head.pop_back();
    const auto fields=Cells(head);
    for(std::size_t i=0;i<fields.size();++i)col[fields[i]]=i;
  }
  bool Next(std::vector<std::string>& values) {
    std::string line;
    if(!std::getline(in,line))return false;
    if(!line.empty()&&line.back()=='\r')line.pop_back();
    values=Cells(line);return true;
  }
  std::string Get(const std::vector<std::string>& v,
                  const std::string& key)const{
    auto it=col.find(key);
    if(it==col.end())throw std::runtime_error(filename+": missing column "+key);
    return it->second<v.size()?v[it->second]:"";
  }
  double Num(const std::vector<std::string>& v,const std::string& key)const {
    const auto s=Get(v,key);
    return s.empty()?std::numeric_limits<double>::quiet_NaN():std::stod(s);
  }
};
double Mean(const std::vector<double>& x) {
  double s=0;int n=0;for(double v:x)if(std::isfinite(v)){s+=v;++n;}
  return n?s/n:std::numeric_limits<double>::quiet_NaN();
}
double RMS(const std::vector<double>& x){
  const double m=Mean(x);if(!std::isfinite(m))return m;
  double s=0;int n=0;for(double v:x)if(std::isfinite(v)){s+=(v-m)*(v-m);++n;}
  return n?std::sqrt(std::max(0.,s/n)):0.;
}
double Sd(const std::vector<double>& x) {
  if(x.size()<2)return 0;
  const double m=Mean(x);double s=0;
  for(double v:x)s+=(v-m)*(v-m);
  return std::sqrt(std::max(0.,s/(x.size()-1)));
}
struct Muon {
  std::string id;
  double gmu=std::numeric_limits<double>::quiet_NaN();
  double Amu=std::numeric_limits<double>::quiet_NaN();
  int nHeed=0;
  std::vector<double> gain; // all primary avalanches, including zero
};
struct Moments {
  double s=0,s2=0;std::size_t n=0;
  void Add(double x){s+=x;s2+=x*x;++n;}
  void Merge(const Moments& b){s+=b.s;s2+=b.s2;n+=b.n;}
  double Mean() const{return n?s/n:0.;}
  double RMS()const{return n?std::sqrt(std::max(0.,s2/n-Mean()*Mean())):0.;}
  double Relative()const{return Mean()>0?RMS()/Mean():0.;}
};
struct VSample {
  int V=0;
  std::vector<Muon> muons;
  std::vector<double> positiveGain,allGain,gmu,A;
  std::unique_ptr<TH1D> hPos,hGmu,hA;
  std::size_t capped=0;
  double zeros=0.;
  double seedMean=0.,posMean=0.,gmuMean=0.;
  double rAll=0.,rPos=0.,rMuon=0.;
  double eAll=0.,ePos=0.,eMuon=0.,eMean=0.;
  double eGmuMean=0.;
};
VSample Load(const std::string& root,int V) {
  VSample a;a.V=V;
  const std::string prefix=root+"/HV_"+std::to_string(V)+"/merged_";
  std::unordered_map<std::string,std::size_t> index;
  {
    CSV f(prefix+"crossing.csv");std::vector<std::string> row;
    const std::string idfield="source|chunk|event|hit_index|track_id|station|det_id";
    while(f.Next(row)) {
      Muon m;
      m.id=f.Get(row,idfield);
      m.nHeed=static_cast<int>(f.Num(row,"heed_primaries"));
      m.gmu=f.Num(row,"G_eff_per_Heed_primary");
      m.Amu=f.Num(row,"A_ion_hemi");
      if(index.count(m.id))throw std::runtime_error("Duplicate crossing "+m.id);
      index[m.id]=a.muons.size();
      if(std::isfinite(m.gmu))a.gmu.push_back(m.gmu);
      if(std::isfinite(m.Amu))a.A.push_back(m.Amu);
      a.muons.push_back(std::move(m));
    }
  }
  {
    CSV f(prefix+"seed.csv");std::vector<std::string> row;
    while(f.Next(row)) {
      const std::string id=f.Get(row,"crossing_id");
      auto it=index.find(id);
      if(it==index.end())throw std::runtime_error("Seed references missing muon "+id);
      const double g=f.Num(row,"collected_e"),cap=f.Num(row,"capped");
      if(!std::isfinite(g)||g<0)throw std::runtime_error("Invalid seed gain");
      a.muons[it->second].gain.push_back(g);
      a.allGain.push_back(g);
      if(g>0)a.positiveGain.push_back(g);
      else a.zeros+=1;
      if(cap!=0)a.capped++;
    }
  }
  for(const auto& m:a.muons)if(static_cast<int>(m.gain.size())!=m.nHeed)
    throw std::runtime_error("Heed seed count mismatch in "+m.id);
  if(a.capped>0)throw std::runtime_error("Capped avalanches in "+std::to_string(V));
  a.seedMean=Mean(a.allGain);
  a.posMean=Mean(a.positiveGain);
  a.gmuMean=Mean(a.gmu);
  a.rAll=RMS(a.allGain)/a.seedMean;
  a.rPos=RMS(a.positiveGain)/a.posMean;
  a.rMuon=RMS(a.gmu)/a.gmuMean;
  return a;
}
void Bootstrap(VSample& sample,int reps=400){
  TRandom3 rng(8241+sample.V);
  std::vector<double> meanAll,ratAll,ratPos,ratMuon,meanMuon;
  const std::size_t N=sample.muons.size();
  if(N<2)return;
  for(int b=0;b<reps;++b) {
    Moments all,pos,gm;
    for(std::size_t j=0;j<N;++j) {
      const Muon& m=sample.muons[rng.Integer(static_cast<UInt_t>(N))];
      if(std::isfinite(m.gmu))gm.Add(m.gmu);
      for(const double g:m.gain) {
        all.Add(g);
        if(g>0)pos.Add(g);
      }
    }
    if(all.n&&pos.n&&gm.n&&all.Mean()>0&&pos.Mean()>0&&gm.Mean()>0) {
      meanAll.push_back(all.Mean());ratAll.push_back(all.Relative());
      ratPos.push_back(pos.Relative());ratMuon.push_back(gm.Relative());
      meanMuon.push_back(gm.Mean());
    }
  }
  sample.eMean=Sd(meanAll);sample.eAll=Sd(ratAll);
  sample.ePos=Sd(ratPos);sample.eMuon=Sd(ratMuon);
  sample.eGmuMean=Sd(meanMuon);
}
void Style(){
  gROOT->SetBatch(kTRUE);gStyle->SetOptStat(0);
  gStyle->SetCanvasColor(kWhite);gStyle->SetPadColor(kWhite);
  gStyle->SetFrameFillColor(kWhite);gStyle->SetFrameLineColor(kBlack);
  gStyle->SetTitleFont(42,"XYZ");gStyle->SetLabelFont(42,"XYZ");
  gStyle->SetTitleSize(.048,"XYZ");gStyle->SetLabelSize(.043,"XYZ");
  gStyle->SetPadLeftMargin(.15);gStyle->SetPadBottomMargin(.14);
  gStyle->SetPadRightMargin(.04);gStyle->SetPadTopMargin(.09);
  gStyle->SetErrorX(0.);
}
void Label(const std::string& note) {
  TLatex t;t.SetNDC(true);t.SetTextFont(42);t.SetTextSize(.044);
  t.DrawLatex(.19,.87,"#bf{NA60+} #it{Simulation}");
  t.SetTextSize(.029);t.DrawLatex(.19,.82,note.c_str());
}
void Save(TCanvas& c,const std::string& folder,const std::string& name){
  c.SaveAs((folder+"/"+name+".png").c_str());
  c.SaveAs((folder+"/"+name+".pdf").c_str());
}
void GraphFrame(const char* title,const double lo,const double hi) {
  // noop: plotting uses TGraph axis directly
}
void DrawGraph(VSample* data[],const std::string& root,
               bool isMean) {
  const std::string folder=root+"/analysis";
  TCanvas c("c_voltage","Voltage trend",900,650);
  TGraphErrors g1(4),g2(4),g3(4);
  for(int i=0;i<4;++i){
    const auto& a=*data[i];
    g1.SetPoint(i,a.V,isMean?a.seedMean:a.rAll);
    g1.SetPointError(i,0,isMean?a.eMean:a.eAll);
    g2.SetPoint(i,a.V,isMean?a.gmuMean:a.rPos);
    g2.SetPointError(i,0,isMean?a.eGmuMean:a.ePos);
    g3.SetPoint(i,a.V,a.rMuon);g3.SetPointError(i,0,a.eMuon);
  }
  for(auto* g:{&g1,&g2,&g3}){
    g->SetMarkerStyle(20);g->SetMarkerSize(1.25);g->SetLineWidth(2);
  }
  g1.SetMarkerColor(kBlack);g1.SetLineColor(kBlack);
  g2.SetMarkerColor(kBlue+1);g2.SetLineColor(kBlue+1);
  g3.SetMarkerColor(kRed+1);g3.SetLineColor(kRed+1);
  g1.SetTitle(isMean?";Anode high voltage [V];Mean collected gain [e^{-}]":";Anode high voltage [V];RMS / mean");
  double maxY=0.;
  for(int i=0;i<4;++i) {
    const auto& a=*data[i];
    if(isMean)maxY=std::max(maxY,std::max(a.seedMean+a.eMean,a.gmuMean+a.eGmuMean));
    else maxY=std::max(maxY,std::max({a.rAll+a.eAll,a.rPos+a.ePos,a.rMuon+a.eMuon}));
  }
  g1.SetMinimum(0.);g1.SetMaximum(maxY>0?1.2*maxY:1.);
  g1.Draw("AP");g1.GetXaxis()->SetLimits(2070,2280);
  g1.GetYaxis()->SetTitleOffset(1.25);
  g2.Draw("P SAME");
  if(!isMean)g3.Draw("P SAME");
  TLegend legend(.57,.72,.90,.88);legend.SetBorderSize(0);legend.SetFillStyle(0);
  legend.AddEntry(&g1,isMean?"Per-primary mean gain":"Per-primary, zeros included","lp");
  legend.AddEntry(&g2,isMean?"Per-muon effective gain":"Per-primary, positive only","lp");
  if(!isMean)legend.AddEntry(&g3,"Per-muon effective gain","lp");
  legend.Draw();
  Label("P3 noFW, Ar/CO_{2} 70:30, Penning r=0.57");
  Save(c,folder,isMean?"mean_gain_vs_HV":"gain_RMS_over_mean_vs_HV");
}
std::unique_ptr<TH1D> MakeHist(const char* stem,int v,
                               const std::vector<double>& x,
                               bool loggain){
  const std::string n=std::string(stem)+std::to_string(v);
  const double xmin=loggain?0.:-1.0,xmax=loggain?6.7:1.0;
  auto h=std::make_unique<TH1D>(n.c_str(),"",loggain?134:40,xmin,xmax);
  h->SetDirectory(nullptr);h->SetLineWidth(2);
  for(const double g:x) {
    if(!std::isfinite(g))continue;
    if(loggain){if(g>0)h->Fill(std::log10(g));}
    else h->Fill(g);
  }
  const double integral=h->Integral("width");
  if(integral>0)h->Scale(1./integral);
  return h;
}
void DrawOverlay(VSample* samples[],const std::string& root,
                 bool isGain){
  TCanvas c("c_hist","Voltage distributions",900,650);
  TLegend legend(.65,.67,.90,.88);legend.SetFillStyle(0);legend.SetBorderSize(0);
  double ymax=0.;
  std::array<std::unique_ptr<TH1D>,4> hs;
  for(int i=0;i<4;++i){
    hs[i]=MakeHist(isGain?"gain":"asym",VS[i],
                   isGain?samples[i]->positiveGain:samples[i]->A,isGain);
    hs[i]->SetLineColor(COL[i]);
    ymax=std::max(ymax,hs[i]->GetMaximum());
    legend.AddEntry(hs[i].get(),(std::to_string(VS[i])+" V").c_str(),"l");
  }
  hs[0]->SetTitle(isGain?";log_{10}(collected gain per primary);Probability density":";Signed muon ion-birth asymmetry A_{#mu};Probability density");
  hs[0]->SetMaximum(ymax*1.3);
  hs[0]->GetYaxis()->SetTitleOffset(1.35);
  hs[0]->Draw("HIST");
  for(int i=1;i<4;++i)hs[i]->Draw("HIST SAME");
  legend.Draw();
  Label(isGain?"Only g > 0; zero response reported separately":"Signed event asymmetry; no Gaussian fit");
  Save(c,root+"/analysis",isGain?"positive_gain_distributions":"signed_asymmetry_distributions");
}
void WriteCSV(const std::array<VSample,4>& data,const std::string& root){
  std::ofstream o(root+"/analysis/gain_voltage_summary.csv");
  o<<"voltage_V,n_muons,n_avalanches,n_positive,zero_fraction,mean_gain_all,mean_gain_all_bootstrap_error,"
      "rms_over_mean_all,bootstrap_err_all,mean_gain_positive,rms_over_mean_positive,"
      "bootstrap_err_positive,mean_muon_effective_gain,bootstrap_err_muon_mean,"
      "rms_over_mean_muon,bootstrap_err_muon\n";
  o<<std::setprecision(12);
  for(const auto& a:data) {
    o<<a.V<<","<<a.muons.size()<<","<<a.allGain.size()<<","<<a.positiveGain.size()<<","
     <<a.zeros/a.allGain.size()<<","<<a.seedMean<<","<<a.eMean<<","
     <<a.rAll<<","<<a.eAll<<","<<a.posMean<<","<<a.rPos<<","<<a.ePos<<","
     <<a.gmuMean<<","<<a.eGmuMean<<","<<a.rMuon<<","<<a.eMuon<<"\n";
  }
}
}
void AnalyzeHVHigh_ROOT(const char* outputDir) {
  using namespace HvAna;
  try{
    Style();
    const std::string root=outputDir;
    gSystem->mkdir((root+"/analysis").c_str(),true);
    std::array<VSample,4> data;
    for(int i=0;i<4;++i){
      data[i]=Load(root,VS[i]);Bootstrap(data[i]);
      std::cout<<"HV="<<VS[i]<<" muons="<<data[i].muons.size()
               <<" primaries="<<data[i].allGain.size()
               <<" meanGain="<<data[i].seedMean
               <<" RMS/mean(all)="<<data[i].rAll
               <<" RMS/mean(pos)="<<data[i].rPos
               <<" RMS/mean(muon)="<<data[i].rMuon<<"\n";
    }
    VSample* ptrs[4] = {&data[0],&data[1],&data[2],&data[3]};
    WriteCSV(data,root);
    DrawGraph(ptrs,root,true);DrawGraph(ptrs,root,false);
    DrawOverlay(ptrs,root,true);DrawOverlay(ptrs,root,false);
    std::cout<<"ROOT HIGH-HV analysis finished: "<<root<<"/analysis/\n";
  }catch(const std::exception& e){
    ::Error("AnalyzeHVHigh_ROOT","%s",e.what());
    throw;
  }
}
