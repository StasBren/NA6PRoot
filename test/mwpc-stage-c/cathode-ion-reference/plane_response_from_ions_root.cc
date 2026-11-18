// NA60+ 2100V microscopic ion-birth archive -> WHOLE-cathode response.
// direct : fully drift EACH SAVED macroion (weighted archives only).
// charge : fast initial-v endpoint-side HYPOTHESIS for charge only;
//          no claimed drift/waveform from FULL archives in this mode.
#include "Plane2100Archive.hh"
#include <TFile.h>
#include <TH1D.h>
#include <TNamed.h>
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
std::string Arg(int n,char**v,const std::string& key) {
 for(int j=1;j+1<n;++j)if(v[j]==key)return v[j+1];
 return "";
}
double Number(int n,char**v,const std::string& key,double def) {
 const auto x=Arg(n,v,key);return x.empty()?def:std::stod(x);
}
}

int main(int argc,char**argv) {
 try {
  const std::string source=Arg(argc,argv,"--ion-file");
  const std::string prefix=Arg(argc,argv,"--out-prefix");
  const std::string mobility=Arg(argc,argv,"--mobility");
  const std::string mode=Arg(argc,argv,"--mode");
  const double dt=Number(argc,argv,"--dt-ns",100.);
  const double windowUs=Number(argc,argv,"--window-us",100.);
  const double maxStepMm=Number(argc,argv,"--max-step-mm",0.05);
  const long long maxDirect=(long long)Number(argc,argv,"--max-direct-ions",20000.);
  NA60Plane::Require(!source.empty()&&!prefix.empty()&&!mobility.empty(),
    "usage: --ion-file ROOT --out-prefix PREFIX --mobility FILE --mode direct|charge");
  NA60Plane::Require(mode=="direct"||mode=="charge","mode must be direct or charge");
  NA60Plane::Require(dt>0.&&windowUs>0.&&windowUs<=500.&&maxDirect>0,
    "invalid time window or direct-event limit");
  const double n=std::ceil(1000.*windowUs/dt);
  NA60Plane::Require(n>=2.&&n<=10000.,"waveform must contain 2..10000 bins");
  const int bins=(int)n;
  NA60Plane::Archive input(source);
  NA60Plane::Require(mode!="direct"||input.mode=="weighted",
    "do not silently drift a FULL archive; use charge mode for FULL");
  NA60Plane::Require(mode!="direct"||input.Entries()<=maxDirect,
    "direct-ion processing cost exceeds user limit; NO truncation performed");
  for(const char* extension:{"_summary.csv","_waveforms.csv","_response.root"}) {
    NA60Plane::Require(!std::filesystem::exists(prefix+extension),
      "output already exists; choose a fresh prefix");
  }
  NA60Plane::Transport tr(mobility,dt,bins,maxStepMm);
  std::array<double,2> qPred{},qDirect{},qRamo{};
  std::vector<double> iPlus(bins,0.),iMinus(bins,0.);
  double lastArrivalNs=0.;
  long long processed=0;
  for(long long j=0;j<input.Entries();++j) {
    const auto b=input.Get(j);
    tr.CheckStart(b.u,b.v,b.w,b.t);
    const auto pred=tr.PredictedEndpointCharge(b.u,b.v,b.w);
    for(int k=0;k<2;++k)qPred[k]+=b.weight*pred[k];
    if(mode=="direct") {
      auto result=tr.Drift(b.u,b.v,b.w,b.t);
      lastArrivalNs=std::max(lastArrivalNs,result.endT);
      for(int k=0;k<2;++k) {
        qDirect[k]+=b.weight*result.integralQ[k];
        qRamo[k]+=b.weight*result.endpointQ[k];
      }
      for(int ib=0;ib<bins;++ib) {
        iPlus[ib]+=b.weight*tr.Current(0,ib);
        iMinus[ib]+=b.weight*tr.Current(1,ib);
      }
    }
    ++processed;
    if(processed%1000==0)
      std::cout<<"processed_real_ion_birth_rows="<<processed
               <<"/"<<input.Entries()<<"\n";
  }
  input.CheckComplete();
  NA60Plane::Require(processed==input.nSaved,
      "partial or missing ion source rows");

  if(mode=="direct") {
    for(int k=0;k<2;++k) {
      const double tol=std::max(0.005*std::abs(qRamo[k]),
                         0.0005*NA60Plane::eFc*input.nTrue);
      NA60Plane::Require(std::abs(qPred[k]-qRamo[k])<=tol,
        "initial-v endpoint-side hypothesis differs from RKF endpoints");
      NA60Plane::Require(std::abs(qDirect[k]-qRamo[k])<=tol,
        "integrated sample waveform differs from endpoint Ramo charge");
    }
  }

  std::ofstream summary(prefix+"_summary.csv");
  NA60Plane::Require(bool(summary),"could not create result summary");
  summary<<std::setprecision(16);
  summary<<"crossing_id,channel,station,det_id,archive_mode,response_mode,"
            "true_ion_births,retained_ion_rows,dt_ns,window_us,max_ion_arrival_ns,"
            "plus_endpoint_side_estimate_fC,minus_endpoint_side_estimate_fC,"
            "plus_RKF_window_fC,minus_RKF_window_fC,"
            "plus_RKF_ramo_fC,minus_RKF_ramo_fC\n";
  summary<<input.id<<","<<input.source<<","<<input.station<<","<<input.det
         <<","<<input.mode<<","<<mode<<","<<input.nTrue<<","<<input.nSaved
         <<","<<dt<<","<<windowUs<<",";
  if(mode=="direct")summary<<lastArrivalNs;else summary<<"NA";
  summary<<","<<qPred[0]<<","<<qPred[1]<<",";
  if(mode=="direct")summary<<qDirect[0]<<","<<qDirect[1]<<","
                           <<qRamo[0]<<","<<qRamo[1]<<"\n";
  else summary<<"NA,NA,NA,NA\n";
  summary.close();

  if(mode=="direct") {
    std::ofstream wave(prefix+"_waveforms.csv");
    NA60Plane::Require(bool(wave),"could not create waveform CSV");
    wave<<std::setprecision(16)
        <<"time_center_ns,ion_plus_fC_per_ns,ion_minus_fC_per_ns,"
          "cumulative_plus_fC,cumulative_minus_fC\n";
    double qp=0,qm=0;
    for(int ib=0;ib<bins;++ib) {
      qp+=iPlus[ib]*dt;qm+=iMinus[ib]*dt;
      wave<<(ib+0.5)*dt<<","<<iPlus[ib]<<","<<iMinus[ib]
          <<","<<qp<<","<<qm<<"\n";
    }
    wave.close();
    const double tol=1.e-8*std::max(1.,std::abs(qp));
    NA60Plane::Require(std::abs(qp-qDirect[0])<=tol &&
                       std::abs(qm-qDirect[1])<=tol,
                       "waveform accumulation lost current");
  }
  TFile root((prefix+"_response.root").c_str(),"RECREATE",
             "NA60+ 2100V WHOLE PLANE pilot",4);
  NA60Plane::Require(!root.IsZombie(),"ROOT output creation failed");
  TNamed note("limitations",
    "Ion ONLY, two WHOLE cathodes, provisional Ar+ in Ar mobility; "
    "no electron signal/real strips/electronics; FULL charge-mode assumes cathode sign(v).");
  note.Write();
  if(mode=="direct") {
    TH1D hp("plus_ion_current","NA60+ Simulation;time (ns);I_{+} (fC/ns)",
            bins,0.,bins*dt);
    TH1D hm("minus_ion_current","NA60+ Simulation;time (ns);I_{-} (fC/ns)",
            bins,0.,bins*dt);
    for(int ib=0;ib<bins;++ib) {
      hp.SetBinContent(ib+1,iPlus[ib]);
      hm.SetBinContent(ib+1,iMinus[ib]);
    }
    hp.Write();hm.Write();
  }
  root.Close();
  std::cout<<std::setprecision(12)<<"PASS plane response input="<<source
           <<" source_mode="<<input.mode<<" response_mode="<<mode
           <<" ions="<<input.nTrue<<" saved="<<input.nSaved<<"\n"
           <<"plus_predicted_Q_fC="<<qPred[0]
           <<" minus_predicted_Q_fC="<<qPred[1]<<"\n";
  if(mode=="direct")
    std::cout<<"plus_RKF_Q_fC="<<qDirect[0]
             <<" minus_RKF_Q_fC="<<qDirect[1]
             <<" max_end_time_ns="<<lastArrivalNs<<"\n";
  std::cout<<"No electronics, no physical stereo strips, no calibrated ion mobility.\n";
  return 0;
 } catch(const std::exception& e) {
  std::cerr<<"PLANE-RESPONSE PILOT FAILED: "<<e.what()<<"\n";return 2;
 }
}
