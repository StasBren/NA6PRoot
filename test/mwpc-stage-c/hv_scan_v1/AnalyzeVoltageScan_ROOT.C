// ROOT-only physics analysis and ATLAS-style plots for the 1000-muon HV scan.
// The statistical unit for error bars is the independent muon, NOT individual
// ion births or individual avalanches. Bootstrap resamples entire muon records.
// root -l -b -q 'AnalyzeVoltageScan_ROOT.C("/path/to/hv_scan_wiretruth_1000_v1")'
#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TColor.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TH1D.h>
#include <TFile.h>
#include <TRandom3.h>
#include <TSystem.h>
#include <TString.h>
#include <TLine.h>
#include <TAxis.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <cstdint>
#include <limits>

namespace NA60Voltage {
using Row=std::map<std::string,std::string>;
const std::vector<int> kVoltages={1700,1750,1800,1850,1900};
std::vector<std::string> split(const std::string& text,char delimiter=','){
  std::vector<std::string> out;std::stringstream ss(text);std::string val;
  while(std::getline(ss,val,delimiter))out.push_back(val);
  if(!text.empty()&&text.back()==delimiter)out.push_back("");
  return out;
}
std::vector<Row> csv(const std::string& path){
  std::ifstream is(path);if(!is)throw std::runtime_error("Cannot open "+path);
  std::string line;if(!std::getline(is,line))throw std::runtime_error("Empty "+path);
  auto names=split(line);std::vector<Row> records;
  while(std::getline(is,line)){
    if(line.empty())continue;
    auto cols=split(line);
    if(cols.size()!=names.size())throw std::runtime_error("Bad CSV column count: "+path);
    Row r;for(std::size_t i=0;i<names.size();++i)r[names[i]]=cols[i];
    records.push_back(std::move(r));
  }
  return records;
}
std::string field(const Row& r,const std::string& name){
  const auto it=r.find(name);
  if(it==r.end())throw std::runtime_error("Missing field: "+name);
  return it->second;
}
double number(const Row& r,const std::string& name,double missing=std::numeric_limits<double>::quiet_NaN()){
  auto s=field(r,name);return s.empty()?missing:std::stod(s);
}
std::string source(const std::string& id){
  if(id.rfind("jpsi",0)==0)return "Jpsi";
  if(id.rfind("phi",0)==0)return "Phi";
  if(id.rfind("omega",0)==0)return "Omega";
  return "Other";
}
struct Wire {double q=0,ion=0,up=0,down=0;};
struct Muon {
  std::string id,ch;
  int station=-1,nseed=0;
  double gain=std::numeric_limits<double>::quiet_NaN();
  double A=std::numeric_limits<double>::quiet_NaN();
  double s=0,s2=0,sp=0,sp2=0;
  int npos=0,zeros=0;
  bool multi=false,opposite=false,hasOpp=false;
  double weightedLocalAbs=std::numeric_limits<double>::quiet_NaN();
  std::vector<Wire> wires;
};
struct ScanPoint {int V=0,chunks=0;std::vector<Muon> muons;};
ScanPoint read(const std::string& output,int V){
  ScanPoint o;o.V=V;
  std::map<std::string,Muon> m;
  for(int i=0;i<50;++i){
    const std::string dir=output+Form("/voltages/V%d/chunk_%04d",V,i);
    std::ifstream pass(dir+"/ROOT_AUDIT_PASS.txt");
    if(!pass.good())continue;
    ++o.chunks;
    const std::string p=dir+"/microscopic_";
    for(const auto& r:csv(p+"crossing.csv")){
      std::string id=field(r,"source|chunk|event|hit_index|track_id|station|det_id");
      if(m.count(id))throw std::runtime_error("Duplicate crossing id at V="+std::to_string(V));
      Muon x;x.id=id;x.ch=source(id);x.station=int(number(r,"station"));
      x.nseed=int(number(r,"heed_primaries"));
      x.gain=number(r,"G_eff_per_Heed_primary");
      x.A=number(r,"A_ion_hemi");
      m[id]=x;
    }
    for(const auto& r:csv(p+"seed.csv")){
      const std::string id=field(r,"crossing_id");
      auto it=m.find(id);if(it==m.end())throw std::runtime_error("seed orphan");
      const double gain=number(r,"collected_e");
      it->second.s+=gain;it->second.s2+=gain*gain;
      if(gain>0){++it->second.npos;it->second.sp+=gain;it->second.sp2+=gain*gain;}
      else ++it->second.zeros;
    }
    for(const auto& r:csv(p+"wire_muon.csv")){
      const std::string id=field(r,"crossing_id");
      auto it=m.find(id);if(it==m.end())throw std::runtime_error("wire orphan");
      Wire w;w.q=number(r,"collected_e");w.ion=number(r,"ion_births");
      w.up=number(r,"ion_up");w.down=number(r,"ion_down");
      it->second.wires.push_back(w);
    }
  }
  if(o.chunks==0)return o;
  o.muons.reserve(m.size());
  for(auto& kv:m){
    auto& x=kv.second;
    if(x.nseed!=x.npos+x.zeros)throw std::runtime_error("seed count mismatch "+x.id);
    double totCharge=0,ions=0,localAbs=0;
    for(const auto& w:x.wires){totCharge+=w.q;ions+=w.ion;}
    std::sort(x.wires.begin(),x.wires.end(),[](const Wire& a,const Wire& b){return a.q>b.q;});
    x.multi=(x.wires.size()>1 && totCharge>0 && x.wires[1].q>0.05*totCharge);
    for(const auto& w:x.wires){
      if(w.up+w.down>0)localAbs+=w.ion*std::abs((w.up-w.down)/(w.up+w.down));
    }
    if(ions>0)x.weightedLocalAbs=localAbs/ions;
    if(x.multi){
      const auto& a=x.wires[0];const auto& b=x.wires[1];
      if(a.up+a.down>0 && b.up+b.down>0){
        const double A1=(a.up-a.down)/(a.up+a.down),A2=(b.up-b.down)/(b.up+b.down);
        x.opposite=A1*A2<0.;x.hasOpp=true;
      }
    }
    o.muons.push_back(x);
  }
  if((int)o.muons.size()!=20*o.chunks)throw std::runtime_error("Unexpected muon count at HV "+std::to_string(V));
  return o;
}
struct Measures {
  double gain=0,cvMuon=0,cvSeed=0,cvPos=0,zero=0,asymMean=0,asymRms=0,asymAbs=0,multi=0,localAbs=0,opposite=0;
};
Measures measure(const std::vector<Muon>& m,const std::vector<std::size_t>& indices){
  Measures x;double eventS=0,eventS2=0,nEvent=0;
  double seedS=0,seedS2=0,ns=0,seedPosS=0,seedPosS2=0,nPos=0,zeros=0;
  double aS=0,aS2=0,aAbs=0,nA=0,multi=0,localS=0,nLocal=0,opp=0,nOpp=0;
  for(auto ix:indices){
    const auto& z=m.at(ix);
    if(std::isfinite(z.gain)){eventS+=z.gain;eventS2+=z.gain*z.gain;++nEvent;}
    seedS+=z.s;seedS2+=z.s2;ns+=z.nseed;
    seedPosS+=z.sp;seedPosS2+=z.sp2;nPos+=z.npos;zeros+=z.zeros;
    if(std::isfinite(z.A)){aS+=z.A;aS2+=z.A*z.A;aAbs+=std::abs(z.A);++nA;}
    multi+=z.multi?1:0;
    if(z.multi && std::isfinite(z.weightedLocalAbs)){localS+=z.weightedLocalAbs;++nLocal;}
    if(z.hasOpp){nOpp++;if(z.opposite)opp++;}
  }
  auto cv=[](double s,double s2,double n){
    if(n<=0||s==0)return std::numeric_limits<double>::quiet_NaN();
    const double mean=s/n;return std::sqrt(std::max(0.,s2/n-mean*mean))/mean;
  };
  x.gain=nEvent?eventS/nEvent:0;
  x.cvMuon=cv(eventS,eventS2,nEvent);
  x.cvSeed=cv(seedS,seedS2,ns);
  x.cvPos=cv(seedPosS,seedPosS2,nPos);
  x.zero=ns?zeros/ns:0;
  x.asymMean=nA?aS/nA:0;
  x.asymRms=nA?std::sqrt(std::max(0.,aS2/nA-x.asymMean*x.asymMean)):0;
  x.asymAbs=nA?aAbs/nA:0;
  x.multi=indices.empty()?0:multi/indices.size();
  x.localAbs=nLocal?localS/nLocal:0;
  x.opposite=nOpp?opp/nOpp:0;
  return x;
}
Measures calc(const std::vector<Muon>& m){
  std::vector<std::size_t> id(m.size());for(std::size_t i=0;i<m.size();++i)id[i]=i;
  return measure(m,id);
}
Measures bootstrapUncertainty(const std::vector<Muon>& m,unsigned seed=4562026){
  if(m.size()<2)return Measures{};
  TRandom3 rng(seed);const unsigned B=400;
  std::vector<Measures> samples; samples.reserve(B);
  for(unsigned b=0;b<B;++b){
    std::vector<std::size_t> take(m.size());
    for(auto& i:take)i=rng.Integer(m.size());
    samples.push_back(measure(m,take));
  }
  Measures mean{},sum2{};
#define ACC(key) {double sm=0,ss=0;for(const auto& v:samples){sm+=v.key;ss+=v.key*v.key;}mean.key=sm/B;sum2.key=std::sqrt(std::max(0.,ss/B-mean.key*mean.key));}
  ACC(gain);ACC(cvMuon);ACC(cvSeed);ACC(cvPos);ACC(zero);ACC(asymMean);ACC(asymRms);ACC(asymAbs);ACC(multi);ACC(localAbs);ACC(opposite);
#undef ACC
  return sum2;
}
void atlas(TPad* p,const char* subtitle){
  p->cd();TLatex t;t.SetNDC();t.SetTextFont(72);t.SetTextSize(.045);t.DrawLatex(.13,.88,"ATLAS");
  t.SetTextFont(42);t.SetTextSize(.040);t.DrawLatex(.28,.88,"Internal");
  t.SetTextSize(.028);t.DrawLatex(.13,.82,subtitle);
}
void style(){
  gROOT->SetBatch(kTRUE);gStyle->SetOptStat(0);gStyle->SetTitleFont(42,"XYZ");
  gStyle->SetLabelFont(42,"XYZ");gStyle->SetPadTickX(1);gStyle->SetPadTickY(1);
  gStyle->SetPadLeftMargin(.12);gStyle->SetPadBottomMargin(.13);
}
TGraphErrors* graph(const std::vector<ScanPoint>& scans,
                     const std::vector<Measures>& values,
                     const std::vector<Measures>& errors,
                     const char* what,Color_t color){
  auto g=new TGraphErrors(scans.size());
  for(std::size_t i=0;i<scans.size();++i){
    const auto& v=values[i];const auto& e=errors[i];
    double y=0,dy=0;
    const std::string w=what;
#define SELECT(key) if(w==#key){y=v.key;dy=e.key;}
    SELECT(gain) else SELECT(cvMuon) else SELECT(cvSeed) else SELECT(cvPos)
    else SELECT(zero) else SELECT(asymMean) else SELECT(asymRms)
    else SELECT(asymAbs) else SELECT(multi) else SELECT(localAbs)
    else SELECT(opposite)
#undef SELECT
    g->SetPoint(i,scans[i].V,y);
    g->SetPointError(i,0,dy);
  }
  g->SetMarkerStyle(20);g->SetMarkerSize(1.25);g->SetMarkerColor(color);
  g->SetLineColor(color);g->SetLineWidth(2);
  return g;
}
void plotOne(const std::string& out,const std::vector<ScanPoint>& scans,
             const std::vector<Measures>& values,const std::vector<Measures>& errors,
             const char* variable,const char* title,const char* unit,
             const char* dest){
  if(scans.empty())return;
  auto c=new TCanvas(Form("c%s",dest),title,900,620);
  auto gr=graph(scans,values,errors,variable,kBlack);
  gr->SetTitle(Form(";%s;%s","Anode voltage [V]",unit));
  gr->Draw("AP");gr->GetXaxis()->SetLimits(1670,1930);
  gr->GetYaxis()->SetTitleOffset(1.15);
  gr->GetYaxis()->SetRangeUser(0,gr->GetYaxis()->GetXmax()*1.18);
  atlas((TPad*)gPad,"1000 Geant4 muons / voltage; bootstrap over muons");
  c->Print((out+"/"+dest+".pdf").c_str());
  c->Print((out+"/"+dest+".png").c_str());delete c;
}
void plotDistributions(const std::string& out,const std::vector<ScanPoint>& scans,
                       TFile& f){
  if(scans.empty())return;
  const int colors[]={kBlack,kAzure+1,kGreen+2,kOrange+7,kRed+1};
  double maxG=5000;
  for(const auto& s:scans)for(const auto& x:s.muons)if(std::isfinite(x.gain))maxG=std::max(maxG,x.gain);
  maxG=std::ceil(maxG/1000.)*1000.;
  auto c=new TCanvas("cDistributions","Full distributions: fixed 1000 Geant4 muons",1200,520);
  c->Divide(2,1);
  TLegend leg(.62,.64,.86,.87);leg.SetBorderSize(0);leg.SetFillStyle(0);
  std::vector<TH1D*> gains,asym;
  for(std::size_t i=0;i<scans.size();++i){
    const auto& s=scans[i];
    auto hg=new TH1D(Form("hGmu_V%d",s.V),";Effective gain G_{#mu} [e / Heed electron];Probability density",80,0,maxG*1.05);
    auto ha=new TH1D(Form("hA_V%d",s.V),";Signed ion-birth asymmetry A_{#mu};Probability density",35,-1.,1.);
    hg->Sumw2();ha->Sumw2();
    for(const auto& mu:s.muons){
      if(std::isfinite(mu.gain))hg->Fill(mu.gain);
      if(std::isfinite(mu.A))ha->Fill(mu.A);
    }
    if(hg->Integral("width")>0)hg->Scale(1./hg->Integral("width"));
    if(ha->Integral("width")>0)ha->Scale(1./ha->Integral("width"));
    for(auto* h:{hg,ha}){h->SetLineColor(colors[i%5]);h->SetLineWidth(2);h->SetStats(0);}
    gains.push_back(hg);asym.push_back(ha);
    leg.AddEntry(hg,Form("%d V (N=%lu)",s.V,(unsigned long)s.muons.size()),"l");
    f.cd();hg->Write();ha->Write();
  }
  c->cd(1);gPad->SetLeftMargin(.14);
  double maxY=0;for(auto* h:gains)maxY=std::max(maxY,h->GetMaximum());
  for(std::size_t i=0;i<gains.size();++i){gains[i]->SetMaximum(maxY*1.50);gains[i]->Draw(i?"HIST SAME":"HIST");}
  leg.Draw();atlas((TPad*)gPad,"Garfield++ exact microscopic avalanches");
  c->cd(2);gPad->SetLeftMargin(.14);maxY=0;for(auto* h:asym)maxY=std::max(maxY,h->GetMaximum());
  for(std::size_t i=0;i<asym.size();++i){asym[i]->SetMaximum(maxY*1.50);asym[i]->Draw(i?"HIST SAME":"HIST");}
  atlas((TPad*)gPad,"Signed A_{#mu}; no Gaussian fit");
  c->Print((out+"/voltage_distributions.pdf").c_str());
  c->Print((out+"/voltage_distributions.png").c_str());
  delete c;
}
} // namespace NA60Voltage

void AnalyzeVoltageScan_ROOT(const char* dir){
  using namespace NA60Voltage;
  try{
    style();
    const std::string base=dir;
    const std::string output=base+"/reports";
    gSystem->mkdir(output.c_str(),kTRUE);
    std::vector<ScanPoint> scans;
    for(const auto V:kVoltages){
      auto x=read(base,V);
      if(x.chunks==0)continue;
      std::cout<<"ROOT analyse "<<V<<" V: "<<x.chunks<<" / 50 chunks, "<<x.muons.size()<<" muons\n";
      scans.push_back(std::move(x));
    }
    if(scans.empty())throw std::runtime_error("No ROOT-audited chunks to analyse");
    std::vector<Measures> vv,ee;
    std::ofstream result(output+"/voltage_summary_ROOT.csv");
    result<<"voltage_V,muons,sealed_chunks,mean_Gmu,SE_mean_Gmu,CV_Gmu,SE_CV_Gmu,"
             "CV_seed_all,SE_CV_seed_all,CV_seed_positive,SE_CV_seed_positive,"
             "zero_fraction,SE_zero_fraction,mean_A_mu,SE_mean_A_mu,RMS_A_mu,SE_RMS_A_mu,"
             "mean_abs_A_mu,SE_mean_abs_A_mu,multiwire_gt5_fraction,SE_multiwire,"
             "mean_local_absA_multi,SE_local_absA_multi,opposite_sign_fraction,SE_opposite_sign\n";
    for(const auto& s:scans){
      const auto x=calc(s.muons);
      const auto e=bootstrapUncertainty(s.muons,20261116+s.V);
      vv.push_back(x);ee.push_back(e);
      result<<std::setprecision(9)<<s.V<<","<<s.muons.size()<<","<<s.chunks<<","<<x.gain<<","<<e.gain<<","<<x.cvMuon<<","<<e.cvMuon<<","<<x.cvSeed<<","<<e.cvSeed<<","<<x.cvPos<<","<<e.cvPos<<","<<x.zero<<","<<e.zero<<","<<x.asymMean<<","<<e.asymMean<<","<<x.asymRms<<","<<e.asymRms<<","<<x.asymAbs<<","<<e.asymAbs<<","<<x.multi<<","<<e.multi<<","<<x.localAbs<<","<<e.localAbs<<","<<x.opposite<<","<<e.opposite<<"\n";
    }
    result.close();
    std::ofstream byStation(output+"/voltage_station_summary_ROOT.csv");
    byStation<<"voltage_V,station,muons,mean_Gmu,CV_Gmu,RMS_A_mu,multiwire_gt5_fraction\n";
    for(const auto& s:scans){
      for(int st=0;st<6;++st){
        std::vector<Muon> sub;
        for(const auto& mu:s.muons)if(mu.station==st)sub.push_back(mu);
        if(sub.empty())continue;
        const auto a=calc(sub);
        byStation<<s.V<<","<<st<<","<<sub.size()<<","<<a.gain<<","<<a.cvMuon<<","<<a.asymRms<<","<<a.multi<<"\n";
      }
    }
    byStation.close();
    TFile rootFile((output+"/voltage_distributions.root").c_str(),"RECREATE");
    plotDistributions(output,scans,rootFile);
    rootFile.Close();
    plotOne(output,scans,vv,ee,"gain","Mean effective gain","Mean G_{#mu} [e / Heed primary]","mean_gain_vs_HV");
    plotOne(output,scans,vv,ee,"cvSeed","Avalanche gain fluctuation","RMS(g) / mean(g), all seeds","seed_gain_RMS_over_mean_vs_HV");
    plotOne(output,scans,vv,ee,"cvMuon","Muon gain fluctuation","RMS(G_{#mu}) / mean(G_{#mu})","muon_gain_RMS_over_mean_vs_HV");
    plotOne(output,scans,vv,ee,"zero","Zero response","Fraction with g = 0","zero_gain_fraction_vs_HV");
    plotOne(output,scans,vv,ee,"asymRms","Asymmetry spread","RMS of signed A_{#mu}","muon_asymmetry_RMS_vs_HV");
    plotOne(output,scans,vv,ee,"multi","Charge sharing",">5% charge on second wire","multiwire_fraction_vs_HV");
    std::cout<<"SUCCESS: ROOT plots and summary tables written to "<<output<<"\n";
    for(const auto& s:scans)if(s.chunks!=50)
      std::cout<<"WARNING: V="<<s.V<<" has only "<<s.chunks<<" chunks; PILOT statistics only\n";
  }catch(const std::exception& e){
    std::cerr<<"ROOT VOLTAGE ANALYSIS FAILED: "<<e.what()<<"\n";
    gSystem->Exit(11);
  }
}
