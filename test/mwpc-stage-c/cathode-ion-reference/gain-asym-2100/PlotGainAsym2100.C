// NA60+ 2100 V exact up/down ion-BIRTH asymmetry plots.
// Input: trusted v0_archive_stats CSV exports, not sampled K256 angular counts.
// No microscopic transport, simulation, or 1800 V comparison.
#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH1F.h>
#include <TGraphErrors.h>
#include <TLatex.h>
#include <TNamed.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TAxis.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
constexpr int NB=24, BOOT=250;
constexpr double XMAX=8.4;
struct Point { int slot; double ion, A; };
struct Bin { double sum=0.,ab=0.; int n=0; };
using Bins=std::array<Bin,NB>;
void check(bool x,const std::string& m){if(!x)throw std::runtime_error("gain/asym: "+m);}
std::vector<std::string> csv(const std::string& s){
 std::stringstream f(s);std::string t;std::vector<std::string> v;
 while(std::getline(f,t,',')){if(!t.empty()&&t.back()=='\r')t.pop_back();v.push_back(t);}
 if(!s.empty()&&s.back()==',')v.push_back("");return v;
}
double num(const std::string& s) {
 try {return s.empty()?NAN:std::stod(s);}catch(...){return NAN;}
}
std::vector<Point> rows(const std::string& file,const std::string& gain,const std::string& a){
 std::ifstream f(file);check(bool(f),"missing "+file);
 std::string s;check(bool(std::getline(f,s)),"empty "+file);
 auto cols=csv(s);std::map<std::string,int> idx;
 for(size_t i=0;i<cols.size();++i)idx[cols[i]]=int(i);
 for(auto name:{"muon_slot"})check(idx.count(name)!=0,std::string("missing ")+name);
 check(idx.count(gain)!=0 && idx.count(a)!=0,"missing "+gain+" or "+a+" in "+file);
 std::vector<Point> v;
 while(std::getline(f,s)){if(s.empty())continue;auto c=csv(s);
  check(c.size()>=cols.size(),"short row in "+file);
  v.push_back({std::stoi(c[idx["muon_slot"]]),num(c[idx[gain]]),num(c[idx[a]])});
 }return v;
}
void style(){gROOT->SetBatch(kTRUE);gStyle->SetOptStat(0);gStyle->SetOptTitle(0);
 gStyle->SetPadTickX(1);gStyle->SetPadTickY(1);gStyle->SetPadGridX(0);gStyle->SetPadGridY(0);
 gStyle->SetLabelFont(42,"XYZ");gStyle->SetTitleFont(42,"XYZ");
 gStyle->SetLabelSize(.043,"XYZ");gStyle->SetTitleSize(.050,"XYZ");
 gStyle->SetTitleOffset(1.20,"X");gStyle->SetTitleOffset(1.35,"Y");
 gStyle->SetFrameLineWidth(2);
}
void label(){TLatex t;t.SetNDC();t.SetTextFont(42);t.SetTextSize(.034);
 t.DrawLatex(.16,.92,"#bf{NA60+} #it{Simulation}   #scale[0.83]{P3 noFW, 2100 V}");
}
void saveHist(TH1* h,const std::string& path,const char* opt="HIST",bool twoD=false){
 TCanvas c("c_hist","",920,690);c.SetTicks(1,1);c.SetLeftMargin(.14);
 c.SetRightMargin(twoD?.17:.06);c.SetBottomMargin(.15);c.SetTopMargin(.10);
 h->SetStats(false);h->SetLineWidth(2);h->SetLineColor(kBlack);
 h->Draw(opt);label();c.SaveAs((path+".png").c_str());c.SaveAs((path+".pdf").c_str());
}
TGraphErrors makeGraph(const std::map<int,Bins>& mu,bool absolute,
                       const std::string& kind,std::ostream& summary){
 TGraphErrors g;g.SetName((kind+(absolute?"_absolute":"_signed")).c_str());
 g.SetMarkerStyle(20);g.SetMarkerColor(kBlue+1);g.SetLineColor(kBlue+1);
 g.SetMarkerSize(.85);
 std::vector<const Bins*> units;for(const auto& kv:mu)units.push_back(&kv.second);
 if(units.empty())return g;
 std::mt19937_64 rng(absolute?441128:441129);
 std::uniform_int_distribution<size_t> choice(0,units.size()-1);
 for(int i=0;i<NB;++i){
  double s=0.;int n=0,individual=0;
  for(auto* unit:units){const auto b=(*unit)[i];s+=absolute?b.ab:b.sum;n+=b.n;individual+=b.n>0;}
  if(n<20 || individual<4)continue;
  double se=0.;std::vector<double> boot;boot.reserve(BOOT);
  for(int j=0;j<BOOT;++j){double sm=0.;int nm=0;
   for(size_t k=0;k<units.size();++k){const auto b=(*units[choice(rng)])[i];
    sm+=absolute?b.ab:b.sum;nm+=b.n;}
   if(nm>0)boot.push_back(sm/nm);
  }
  if(boot.size()>1){double ave=std::accumulate(boot.begin(),boot.end(),0.)/boot.size();
   for(double v:boot)se+=(v-ave)*(v-ave);
   se=std::sqrt(se/(boot.size()-1));}
  double loglo=XMAX*i/NB,loghi=XMAX*(i+1)/NB;
  double x=std::pow(10.,(loglo+loghi)/2)-1.;
  int p=g.GetN();g.SetPoint(p,x,s/n);g.SetPointError(p,0.,se);
  summary<<kind<<','<<(absolute?"abs":"signed")<<','<<loglo<<','<<loghi<<','
         <<n<<','<<individual<<','<<std::setprecision(10)<<s/n<<','<<se<<'\n';
 }return g;
}
void saveGraph(TGraphErrors& g,const std::string& out,bool absolute,const char* subject){
 if(g.GetN()==0)return;
 TCanvas c("c_graph","",920,690);c.SetTicks(1,1);c.SetLogx();
 c.SetLeftMargin(.14);c.SetRightMargin(.06);c.SetBottomMargin(.15);c.SetTopMargin(.10);
 std::string axis=std::string(";Ionization gain ")+subject+";"+
  (absolute?"Mean |A_{UD}|":"Mean A_{UD}");
 TH1F* f=c.DrawFrame(1,absolute?0:-1.05,1e8,1.05,axis.c_str());
 f->GetXaxis()->SetMoreLogLabels();f->GetXaxis()->SetNoExponent();
 g.Draw("P SAME");label();c.SaveAs((out+".png").c_str());c.SaveAs((out+".pdf").c_str());
}
}

void PlotGainAsym2100(const char* csvDirectory,const char* outDirectory){
 check(csvDirectory&&outDirectory,"two directory arguments are required");
 style();std::string in=csvDirectory,out=outDirectory;
 check(gSystem->mkdir(out.c_str(),true)==0||!gSystem->AccessPathName(out.c_str()),"cannot make output");
 auto pr=rows(in+"/v0_primaries.csv","ion_births_exact","exact_A");
 auto mu=rows(in+"/v0_muons.csv","ion_births_exact","A_exact");
 auto np=rows(in+"/v0_muons.csv","n_primaries","A_exact");
 check(!pr.empty()&&!mu.empty()&&mu.size()==np.size(),"missing population entries");
 TH1D hP("primary_hemispheric_A",";A_{UD};Heed primary avalanches",60,-1.,1.);
 TH1D hM("muon_hemispheric_A",";A_{UD};Muons",60,-1.,1.);
 TH2D h2P("primary_A_vs_ion_births",";log_{10}(1+ion births per primary);A_{UD}",72,0,XMAX,60,-1,1);
 TH2D h2M("muon_A_vs_ion_gain",";log_{10}(1+G_{ion,muon});A_{UD}",72,0,XMAX,60,-1,1);
 TH1D hG("primary_ion_births",";log_{10}(1+ion births per primary);Heed primaries",72,0,XMAX);
 std::map<int,Bins> primaryBins,muonBins;
 long long nPositive=0,nZero=0,nUndef=0,nMuValid=0;double allIon=0,allHeed=0;
 for(const auto& r:pr) {
  if(!std::isfinite(r.ion)||r.ion<0)throw std::runtime_error("negative/invalid primary ion count");
  hG.Fill(std::log10(1.+r.ion));
  if(r.ion==0){++nZero;continue;}
  if(!std::isfinite(r.A)){++nUndef;continue;}
  ++nPositive;hP.Fill(r.A);h2P.Fill(std::log10(1.+r.ion),r.A);
  int b=int(std::log10(1.+r.ion)/XMAX*NB);
  if(b>=0&&b<NB){auto& x=primaryBins[r.slot][b];x.sum+=r.A;x.ab+=std::abs(r.A);++x.n;}
 }
 for(size_t i=0;i<mu.size();++i) {
  check(mu[i].slot==np[i].slot,"muon CSV join failed");
  double n=np[i].ion,nIon=mu[i].ion;
  check(std::isfinite(n)&&n>=0&&std::isfinite(nIon)&&nIon>=0,"invalid muon counts");
  allHeed+=n;allIon+=nIon;
  if(n<=0||!std::isfinite(mu[i].A))continue;
  ++nMuValid;double gain=nIon/n,a=mu[i].A;
  hM.Fill(a);h2M.Fill(std::log10(1.+gain),a);
  int b=int(std::log10(1.+gain)/XMAX*NB);
  if(b>=0&&b<NB){auto& x=muonBins[mu[i].slot][b];x.sum+=a;x.ab+=std::abs(a);++x.n;}
 }
 check(nPositive>0&&nMuValid>0&&allHeed>0,"no exact Garfield hemisphere census available");
 std::ofstream binned(out+"/gain_asym_profiles.csv");
 binned<<"unit,observable,log10G_lo,log10G_hi,n_entries,n_muons,mean,muon_bootstrap_sd\n";
 auto pS=makeGraph(primaryBins,false,"primary",binned);
 auto pA=makeGraph(primaryBins,true,"primary",binned);
 auto mS=makeGraph(muonBins,false,"muon",binned);
 auto mA=makeGraph(muonBins,true,"muon",binned);
 TFile root((out+"/gain_asym_2100.root").c_str(),"CREATE");
 check(!root.IsZombie(),"cannot create new ROOT output");
 hP.Write();hM.Write();h2P.Write();h2M.Write();hG.Write();
 pS.Write();pA.Write();mS.Write();mA.Write();
 TNamed scope("scope","Exact microscopic A_ion_hemi (v>0 vs v<0), boundary excluded from A denominator. Ion births only, not ion drift or collected cathode charge. Includes full 997 descriptive population: not an independent model validation split.");
 scope.Write();root.Close();
 saveHist(&hP,out+"/01_primary_A_distribution");
 saveHist(&hM,out+"/02_muon_A_distribution");
 saveHist(&hG,out+"/03_ion_births_per_primary");
 saveHist(&h2P,out+"/04_primary_A_vs_ion_births_2D","COLZ",true);
 saveHist(&h2M,out+"/05_muon_A_vs_ion_gain_2D","COLZ",true);
 saveGraph(pS,out+"/06_primary_mean_A_vs_gain",false,"per primary");
 saveGraph(pA,out+"/07_primary_mean_absA_vs_gain",true,"per primary");
 saveGraph(mS,out+"/08_muon_mean_A_vs_gain",false,"per muon");
 saveGraph(mA,out+"/09_muon_mean_absA_vs_gain",true,"per muon");
 std::ofstream summary(out+"/summary.txt");
 summary<<std::setprecision(12)<<"NA60+ 2100 V original Garfield microscopic BIRTH counts\n"
        <<"muons="<<mu.size()<<" primaries="<<pr.size()<<" positive_A_defined="<<nPositive
        <<" zero_ion="<<nZero<<" positive_A_undefined="<<nUndef
        <<" valid_muon_A="<<nMuValid<<"\n"
        <<"global_mean_ion_gain_ratio_of_sums="<<allIon/allHeed<<"\n"
        <<"Errors: muon-block bootstrap, 250 replicates; >=20 rows and >=4 muons/bin.\n"
        <<"DESCRIPTIVE ONLY (includes sealed development/validation holdouts).\n";
 std::cout<<"PASS 2100 V exact ion-birth hemispheric gain plots: "<<out<<std::endl;
}
