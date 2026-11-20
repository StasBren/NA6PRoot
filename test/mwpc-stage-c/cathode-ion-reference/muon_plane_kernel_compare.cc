// NA60+ Stage C: WHOLE-MUON two-COMPLETE-plane ion current via response templates.
// The input ion births are from a real Geant4/Heed/microscopic Garfield muon.
// Reference is the SAME archived weighted ions evaluated with direct ion RKF.
// No new avalanche simulation, NO prompt electrons, NO stereo strips/electronics.
// A nearest-neighbor kernel is an approximation, not an exact ion drift.
//
// The output comparison separates kernel interpolation error (same reservoir)
// from the unmeasured variance due to ion-birth sampling (different reservoir).
#include "Plane2100Archive.hh"
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
#include <chrono>
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
#include <utility>
#include <vector>

namespace {
using Clock=std::chrono::steady_clock;
using Birth=NA60Plane::Archive::Birth;
struct Template {
 double u=0,v=0;
 std::vector<double> wave;
};
struct Wave {
 std::array<std::vector<double>,2> cur;
 int n=0;
 double dt=0;
};
std::string Arg(int argc,char**argv,const std::string&opt) {
 for(int i=1;i+1<argc;++i)if(opt==argv[i])return argv[i+1];
 return "";
}
double Number(int argc,char**argv,const std::string&opt,double def) {
 const auto s=Arg(argc,argv,opt);
 return s.empty()?def:std::stod(s);
}
double Dist(const Birth&x,const Template&y) {
 // Exactly the Cartesian+radial distance used by the existing raw-IDW
 // benchmark, NOT a new choice optimized against the target muon.
 const double du=(x.u-y.u)/.003,dv=(x.v-y.v)/.003;
 const double r0=std::hypot(x.u,x.v),r1=std::hypot(y.u,y.v);
 NA60Plane::Require(r0>0&&r1>0,"invalid zero radius");
 const double dr=std::log(r0/r1)/.4;
 return du*du+dv*dv+dr*dr;
}
std::vector<Template> ReadBank(const std::string&path,int bins,double dt,
                                std::string&origin) {
 std::unique_ptr<TFile> f(TFile::Open(path.c_str(),"READ"));
 NA60Plane::Require(f&&!f->IsZombie(),"cannot read kernel training bank");
 auto* bank=dynamic_cast<TTree*>(f->Get("training_bank"));
 auto* ax=dynamic_cast<TH1D*>(f->Get("average_plus"));
 auto* meta=dynamic_cast<TNamed*>(f->Get("input_ion_archive"));
 NA60Plane::Require(bank&&ax&&meta&&bank->GetEntries()>0,
                    "missing original compact-pilot training metadata");
 NA60Plane::Require(ax->GetNbinsX()==bins &&
   std::abs(ax->GetXaxis()->GetXmax()-dt*bins)<1.e-7*dt*bins,
   "training waveform time grid differs from requested event waveform");
 origin=meta->GetTitle();
 double u=0,v=0;
 std::vector<double>* vec=nullptr;
 bank->SetBranchAddress("u_cm",&u);
 bank->SetBranchAddress("v_cm",&v);
 bank->SetBranchAddress("waveform_plus_then_minus",&vec);
 std::vector<Template> out;
 out.reserve(2*bank->GetEntries());
 int positive=0,negative=0;
 for(Long64_t j=0;j<bank->GetEntries();++j) {
   NA60Plane::Require(bank->GetEntry(j)>0&&vec,
                      "training bank corrupted or missing row");
   NA60Plane::Require(std::isfinite(u)&&std::isfinite(v)&&v!=0.&&
       vec->size()==2u*(unsigned)bins,"invalid training impulse response");
   for(double x:*vec)NA60Plane::Require(std::isfinite(x),
                       "nonfinite training current");
   out.push_back({u,v,*vec});
   if(v>0)++positive;else ++negative;
 }
 // The P3 noFW analytic geometry has identical grounded cathode planes
 // at v=+-0.3 cm and no magnetic field. Reflection v -> -v interchanges
 // plus and minus weighting electrodes without changing current sign:
 // K_plus(u,-v,t)=K_minus(u,+v,t), and vice versa.
 // Extend a SINGLE-HEMISPHERE training bank by exact model symmetry
 // rather than inventing a statistical kernel or discarding negative ions.
 // This is NOT valid for asymmetric cathode gaps, B fields or real strips.
 if(positive==0 || negative==0) {
   const auto original=out.size();
   for(std::size_t i=0;i<original;++i) {
     Template reflected=out[i];
     reflected.v=-reflected.v;
     for(int j=0;j<bins;++j) {
       reflected.wave[j]=out[i].wave[bins+j];
       reflected.wave[bins+j]=out[i].wave[j];
     }
     out.push_back(std::move(reflected));
   }
   std::cout<<"MODEL SYMMETRY: expanded "<<original
            <<" original training templates by "<<original
            <<" v-reflected/electrode-swapped templates; "
            <<"exact for symmetric P3 noFW two-WHOLE-plane model ONLY.\n";
 }
 NA60Plane::Require(out.size()>=4,"too few response templates");
 return out;
}
std::vector<std::pair<int,double>> Nearest(const Birth&b,
      const std::vector<Template>&bank,int k) {
 std::vector<std::pair<double,int>> distances;
 distances.reserve(bank.size());
 for(int j=0;j<(int)bank.size();++j) {
   if((bank[j].v>0)!=(b.v>0))continue;
   distances.emplace_back(Dist(b,bank[j]),j);
 }
 NA60Plane::Require((int)distances.size()>=k,
   "kernel training bank lacks enough ions on target cathode hemisphere");
 std::partial_sort(distances.begin(),distances.begin()+k,distances.end());
 std::vector<std::pair<int,double>> out;
 if(distances[0].first<1.e-16) {
   out.emplace_back(distances[0].second,1.);
   return out;
 }
 double total=0;
 for(int j=0;j<k;++j) {
   const double weight=1./std::max(1.e-12,distances[j].first);
   out.emplace_back(distances[j].second,weight);
   total+=weight;
 }
 for(auto&item:out)item.second/=total;
 return out;
}
Wave EmptyWave(int bins,double dt) {
 Wave w;w.n=bins;w.dt=dt;
 for(auto&x:w.cur)x.assign(bins,0.);
 return w;
}
// Match the direct-RKF SUMMARY to the exact ROOT crossing, reservoir and
// time grid. The waveform CSV alone does not contain event identity.
void ValidateReferenceSummary(const std::string&path,
                 const NA60Plane::Archive&input,double dt,double windowUs) {
 std::ifstream f(path);
 NA60Plane::Require(bool(f),"cannot read direct RKF summary "+path);
 std::string header,line,field;
 NA60Plane::Require(bool(std::getline(f,header)) &&
                    bool(std::getline(f,line)),
                    "empty direct RKF summary");
 const std::string expected=
   "crossing_id,channel,station,det_id,archive_mode,response_mode,"
   "true_ion_births,retained_ion_rows,dt_ns,window_us,max_ion_arrival_ns,"
   "plus_endpoint_side_estimate_fC,minus_endpoint_side_estimate_fC,"
   "plus_RKF_window_fC,minus_RKF_window_fC,"
   "plus_RKF_ramo_fC,minus_RKF_ramo_fC";
 const std::string lateSuffix=
   ",late_ion_rows,late_ion_weighted_population,"
   "plus_unrecorded_tail_fC,minus_unrecorded_tail_fC,"
   "max_late_arrival_ns";
 NA60Plane::Require(header==expected||header==expected+lateSuffix,
    "unexpected direct RKF summary schema");
 std::istringstream csv(line);
 std::vector<std::string> cells;
 while(std::getline(csv,field,','))cells.push_back(field);
 NA60Plane::Require(cells.size()==(header==expected?17u:22u),
                    "direct RKF summary record column count wrong");
 NA60Plane::Require(cells[0]==input.id &&
   cells[2]==std::to_string(input.station) &&
   cells[3]==std::to_string(input.det) &&
   cells[4]=="weighted"&&cells[5]=="direct"&&
   std::stoll(cells[6])==input.nTrue &&
   std::stoll(cells[7])==input.nSaved &&
   std::abs(std::stod(cells[8])-dt)<1.e-10*dt &&
   std::abs(std::stod(cells[9])-windowUs)<1.e-10*windowUs,
   "reference belongs to a DIFFERENT muon/reservoir/time grid");
}
Wave ReadRef(const std::string&path,int bins,double dt) {
 std::ifstream f(path);
 NA60Plane::Require(bool(f),"could not read direct-RKF waveform CSV: "+path);
 std::string line;
 NA60Plane::Require(bool(std::getline(f,line)) &&
     line=="time_center_ns,ion_plus_fC_per_ns,ion_minus_fC_per_ns,"
           "cumulative_plus_fC,cumulative_minus_fC",
     "unexpected direct RKF CSV schema");
 Wave w=EmptyWave(bins,dt);
 for(int j=0;j<bins;++j) {
   NA60Plane::Require(bool(std::getline(f,line)),
                      "reference ends before requested time window");
   std::replace(line.begin(),line.end(),',',' ');
   std::istringstream s(line);
   double t=0,plus=0,minus=0,cp=0,cm=0;
   NA60Plane::Require(bool(s>>t>>plus>>minus>>cp>>cm),
                      "malformed reference current row");
   NA60Plane::Require(std::abs(t-(j+.5)*dt)<1.e-6*dt &&
        std::isfinite(plus)&&std::isfinite(minus),
        "reference and kernel time grids do not match");
   w.cur[0][j]=plus;w.cur[1][j]=minus;
 }
 while(std::getline(f,line))NA60Plane::Require(line.empty(),
                       "reference contains extra current rows");
 return w;
}
struct Loss {double plus=0,minus=0;};
void Deposit(Wave&out,const std::vector<double>&impulse,int bins,
             const Birth&b,double coefficient,Loss&missing) {
 NA60Plane::Require(b.t>=0&&std::isfinite(b.t),
                    "invalid real microscopic ion creation time");
 const double shift=b.t/out.dt;
 const int offset=(int)std::floor(shift);
 const double frac=shift-offset;
 NA60Plane::Require(frac>=0&&frac<1.&&offset>=0,
                    "invalid fractional time shift");
 for(int c=0;c<2;++c) {
   double lost=0.;
   for(int j=0;j<bins;++j) {
     const double i=coefficient*impulse[c*bins+j];
     if(i==0)continue;
     const int first=j+offset, second=first+1;
     if(first<bins)out.cur[c][first]+=(1.-frac)*i;
     else lost+=(1.-frac)*i*out.dt;
     if(second<bins)out.cur[c][second]+=frac*i;
     else lost+=frac*i*out.dt;
   }
   if(c==0)missing.plus+=lost;else missing.minus+=lost;
 }
}
std::array<double,2> Charges(const Wave&w) {
 std::array<double,2> q{};
 for(int c=0;c<2;++c)
   for(double x:w.cur[c])q[c]+=w.dt*x;
 return q;
}
struct Metric {
 double relL2=0,plusRelL2=0,minusRelL2=0;
 double maxCumFc=0,early1usFc=0,early10usFc=0;
 std::array<double,2> peakTRefUs{},peakTFastUs{};
};
Metric Compare(const Wave&fast,const Wave&ref) {
 NA60Plane::Require(fast.n==ref.n&&fast.dt==ref.dt,"reference grid mismatch");
 Metric m;
 double totalError2=0,totalRef2=0;
 for(int c=0;c<2;++c) {
   double se=0,sr=0,cumulative=0;
   int peakRef=0,peakFast=0;
   for(int j=0;j<ref.n;++j) {
     const double y=fast.cur[c][j],r=ref.cur[c][j];
     const double d=y-r;
     se+=d*d;sr+=r*r;
     totalError2+=d*d;totalRef2+=r*r;
     cumulative+=d*ref.dt;
     m.maxCumFc=std::max(m.maxCumFc,std::abs(cumulative));
     if((j+1)*ref.dt<=1000.+1.e-10)
       m.early1usFc=std::max(m.early1usFc,std::abs(cumulative));
     if((j+1)*ref.dt<=10000.+1.e-10)
       m.early10usFc=std::max(m.early10usFc,std::abs(cumulative));
     if(std::abs(r)>std::abs(ref.cur[c][peakRef]))peakRef=j;
     if(std::abs(y)>std::abs(fast.cur[c][peakFast]))peakFast=j;
   }
   if(c==0)m.plusRelL2=std::sqrt(se/std::max(sr,1.e-60));
   else m.minusRelL2=std::sqrt(se/std::max(sr,1.e-60));
   m.peakTRefUs[c]=(peakRef+.5)*ref.dt/1000.;
   m.peakTFastUs[c]=(peakFast+.5)*ref.dt/1000.;
 }
 m.relL2=std::sqrt(totalError2/std::max(totalRef2,1.e-60));
 return m;
}
void Plot(const Wave&fast,const Wave&ref,const std::string&prefix,
          TFile&root) {
 gROOT->SetBatch(kTRUE);
 gStyle->SetOptStat(0);
 TCanvas canvas("muon_planes","NA60+ simulated muon cathode planes",1400,900);
 canvas.Divide(2,2);
 std::vector<std::unique_ptr<TH1D>> hist;
 hist.reserve(8);
 for(int c=0;c<2;++c) {
   const std::string side=c==0?"Upper (+v)":"Lower (-v)";
   for(int kind=0;kind<2;++kind) {
     std::vector<double> cumRef(ref.n,0.),cumFast(ref.n,0.);
     if(kind==1) {
       for(int j=0;j<ref.n;++j) {
         cumRef[j]=(j?cumRef[j-1]:0.)+ref.dt*ref.cur[c][j];
         cumFast[j]=(j?cumFast[j-1]:0.)+fast.dt*fast.cur[c][j];
       }
     }
     const std::string id=std::to_string(c)+"_"+std::to_string(kind);
     const std::string label=side+(kind==0?" ion current":" cumulative ion charge");
     const std::string axes=";t since window start (#mus);"+
        std::string(kind==0?"current (fC/ns)":"induced charge (fC)");
     auto hr=std::make_unique<TH1D>(("ref_"+id).c_str(),
                        (label+axes).c_str(),ref.n,0,ref.n*ref.dt/1000.);
     auto hf=std::make_unique<TH1D>(("kernel_"+id).c_str(),
                        (label+axes).c_str(),ref.n,0,ref.n*ref.dt/1000.);
     for(int j=0;j<ref.n;++j) {
       hr->SetBinContent(j+1,kind==0?ref.cur[c][j]:cumRef[j]);
       hf->SetBinContent(j+1,kind==0?fast.cur[c][j]:cumFast[j]);
     }
     hr->SetLineColor(kBlack);
     hf->SetLineColor(kBlue+1);
     hf->SetLineStyle(2);
     hr->SetLineWidth(2);hf->SetLineWidth(2);
     const double low=std::min(hr->GetMinimum(),hf->GetMinimum());
     const double high=std::max(hr->GetMaximum(),hf->GetMaximum());
     const double span=std::max(1.e-10,high-low);
     hr->SetMinimum(low-.12*span);
     hr->SetMaximum(high+.12*span);
     canvas.cd(kind*2+c+1);
     hr->Draw("HIST");
     hf->Draw("HIST SAME");
     TLegend leg(.48,.74,.9,.88);
     leg.SetBorderSize(0);leg.SetFillStyle(0);
     leg.AddEntry(hr.get(),"Weighted RKF direct","l");
     leg.AddEntry(hf.get(),"Fast kernel (IDW)","l");
     leg.DrawClone();
     hist.push_back(std::move(hr));hist.push_back(std::move(hf));
   }
 }
 canvas.cd(0);
 canvas.Update();
 root.cd();
 for(auto&h:hist)h->Write();
 canvas.Write("muon_plane_overlay");
 canvas.SaveAs((prefix+"_muon_overlay.png").c_str());
}
} // namespace
int main(int argc,char**argv) {
 try{
  const std::string ions=Arg(argc,argv,"--ion-file");
  const std::string training=Arg(argc,argv,"--training-root");
  const std::string direct=Arg(argc,argv,"--reference-waveforms");
  const std::string directSummary=Arg(argc,argv,"--reference-summary");
  const std::string prefix=Arg(argc,argv,"--out-prefix");
  const int neighbors=(int)Number(argc,argv,"--neighbors",4);
  const double dt=Number(argc,argv,"--dt-ns",100.);
  const double tmaxUs=Number(argc,argv,"--window-us",100.);
  NA60Plane::Require(!ions.empty()&&!training.empty()&&!direct.empty()&&
      !directSummary.empty()&&!prefix.empty(),
    "usage: --ion-file weighted.root --training-root basis.root "
    "--reference-waveforms direct_waveforms.csv "
    "--reference-summary direct_summary.csv --out-prefix PREFIX");
  NA60Plane::Require(neighbors>=1&&neighbors<=8&&dt>0&&tmaxUs>0&&tmaxUs<=500,
                     "invalid event kernel grid/neighbor config");
  const int bins=(int)std::ceil(1000.*tmaxUs/dt);
  NA60Plane::Require(bins>=2&&bins<=10000,"event time bins must be 2..10000");
  for(const std::string& ext:{"_summary.csv","_waveforms.csv",
                              "_response.root","_muon_overlay.png"}) {
    NA60Plane::Require(!std::filesystem::exists(prefix+ext),
                "refusing to overwrite existing result "+prefix+ext);
  }
  const auto start=Clock::now();
  std::string trainingOrigin;
  const auto bank=ReadBank(training,bins,dt,trainingOrigin);
  const auto loaded=Clock::now();
  NA60Plane::Archive input(ions);
  NA60Plane::Require(input.mode=="weighted","require weighted primary-accurate ROOT archive");
  ValidateReferenceSummary(directSummary,input,dt,tmaxUs);
  Wave kernel=EmptyWave(bins,dt);
  Loss truncated;
  for(long long j=0;j<input.Entries();++j) {
    const auto birth=input.Get(j);
    NA60Plane::Require(std::isfinite(birth.u)&&std::isfinite(birth.v)&&
                      std::isfinite(birth.t)&&birth.v!=0,
                      "invalid/unsupported ion birth");
    const auto near=Nearest(birth,bank,neighbors);
    // IDW interpolates actual precomputed Garfield current templates;
    // macroion weight is EXACT original per-primary Horvitz-Thompson weight.
    for(const auto &[index,weight]:near)
      Deposit(kernel,bank[index].wave,bins,birth,birth.weight*weight,truncated);
  }
  input.CheckComplete();
  const auto constructed=Clock::now();
  const auto reference=ReadRef(direct,bins,dt);
  const auto metric=Compare(kernel,reference);
  const auto q=Charges(kernel), qr=Charges(reference);
  const auto compared=Clock::now();

  std::ofstream summary(prefix+"_summary.csv");
  NA60Plane::Require(bool(summary),"cannot create summary");
  summary<<std::setprecision(15);
  summary<<"crossing_id,station,det_id,true_ion_births,retained_macroions,"
            "training_templates,neighbors,dt_ns,window_us,"
            "plus_direct_fC,plus_kernel_fC,plus_error_fC,"
            "minus_direct_fC,minus_kernel_fC,minus_error_fC,"
            "combined_L2,plus_L2,minus_L2,max_cumulative_error_fC,"
            "max_early_1us_cumulative_error_fC,"
            "max_early_10us_cumulative_error_fC,"
            "plus_peak_time_direct_us,plus_peak_time_kernel_us,"
            "minus_peak_time_direct_us,minus_peak_time_kernel_us,"
            "truncated_kernel_plus_fC,truncated_kernel_minus_fC,"
            "bank_loading_seconds,kernel_construction_seconds,comparison_seconds\n";
  summary<<input.id<<","<<input.station<<","<<input.det<<","
         <<input.nTrue<<","<<input.nSaved<<","<<bank.size()<<","
         <<neighbors<<","<<dt<<","<<tmaxUs<<","
         <<qr[0]<<","<<q[0]<<","<<(q[0]-qr[0])<<","
         <<qr[1]<<","<<q[1]<<","<<(q[1]-qr[1])<<","
         <<metric.relL2<<","<<metric.plusRelL2<<","
         <<metric.minusRelL2<<","<<metric.maxCumFc<<","
         <<metric.early1usFc<<","<<metric.early10usFc<<","
         <<metric.peakTRefUs[0]<<","<<metric.peakTFastUs[0]<<","
         <<metric.peakTRefUs[1]<<","<<metric.peakTFastUs[1]<<","
         <<truncated.plus<<","<<truncated.minus<<","
         <<std::chrono::duration<double>(loaded-start).count()<<","
         <<std::chrono::duration<double>(constructed-loaded).count()<<","
         <<std::chrono::duration<double>(compared-constructed).count()<<"\n";
  summary.close();
  std::ofstream csv(prefix+"_waveforms.csv");
  NA60Plane::Require(bool(csv),"cannot create waveform CSV");
  csv<<std::setprecision(15)
     <<"time_center_ns,plus_kernel_fC_per_ns,minus_kernel_fC_per_ns,"
       "plus_direct_fC_per_ns,minus_direct_fC_per_ns,"
       "plus_kernel_cumulative_fC,minus_kernel_cumulative_fC,"
       "plus_direct_cumulative_fC,minus_direct_cumulative_fC\n";
  double ck[2]{},cr[2]{};
  for(int j=0;j<bins;++j) {
    for(int c=0;c<2;++c) {
      ck[c]+=kernel.cur[c][j]*dt;
      cr[c]+=reference.cur[c][j]*dt;
    }
    csv<<(j+.5)*dt<<","<<kernel.cur[0][j]<<","<<kernel.cur[1][j]
       <<","<<reference.cur[0][j]<<","<<reference.cur[1][j]
       <<","<<ck[0]<<","<<ck[1]<<","<<cr[0]<<","<<cr[1]<<"\n";
  }
  csv.close();

  TFile root((prefix+"_response.root").c_str(),"RECREATE",
              "NA60+ muon two-plane ion-only signal",4);
  NA60Plane::Require(!root.IsZombie(),"cannot create ROOT muon signal output");
  TNamed provenance("provenance",
    ("event_ions="+ions+"\ntrained_ion_archive="+trainingOrigin+
     "\ntraining_basis="+training+"\nweighted_RKF_reference="+direct+
     "\nweighted_RKF_summary="+directSummary).c_str());
  provenance.Write();
  TNamed limits("limits",
   "This is a COMPLETE-MUON ION-ONLY current for both WHOLE cathode planes "
   "using weighted Garfield ion births and raw nearest response templates. "
   "For single-hemisphere banks, reflect v and swap upper/lower currents "
   "using the exact symmetry of the ideal P3 noFW analytic model. "
   "No prompt electrons, physical U/V strips or electronics. Ion transport "
   "training mobility: provisional Ar+ in pure Ar for Ar/CO2 70:30. "
   "Training and evaluation may share the same underlying muon. "
   "Accuracy vs direct RKF is a KERNEL-only comparison, NOT reservoir or physical error.");
  limits.Write();
  Plot(kernel,reference,prefix,root);
  root.Close();
  std::cout<<std::setprecision(8)<<"PASS complete-muon ION-ONLY whole-plane kernel test\n"
    <<"crossing="<<input.id<<" true_births="<<input.nTrue
    <<" weighted_birth_rows="<<input.nSaved<<"\n"
    <<"plus_RKF="<<qr[0]<<" fC plus_kernel="<<q[0]<<" fC\n"
    <<"minus_RKF="<<qr[1]<<" fC minus_kernel="<<q[1]<<" fC\n"
    <<"relative_L2 plus="<<metric.plusRelL2
    <<" minus="<<metric.minusRelL2
    <<" combined="<<metric.relL2<<"\n"
    <<"kernel_construction_s="
    <<std::chrono::duration<double>(constructed-loaded).count()<<"\n"
    <<"IMPORTANT: PASS = completed calculation, not accepted physics accuracy.\n"
    <<"NO electron currents, NO stereo-strip readout, NO electronics.\n";
  return 0;
 }catch(const std::exception& e) {
  std::cerr<<"MUON PLANE KERNEL FAILURE: "<<e.what()<<"\n";
  return 2;
 }
}
