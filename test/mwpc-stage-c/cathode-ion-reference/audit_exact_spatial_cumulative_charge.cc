// Exact-position cumulative induced-charge response pilot.
// ORIGINAL already examined muon0301 weighted 4096 ion births: no new
// Heed, AvalancheMicroscopic, ion reservoir, FEM, strip or electronics.
// Recompute direct 10ns whole-plane reference and compare to a response
// constructed solely from the t=0 RKF drift points and electrode weighting
// potentials Q_k(tau)=e[phi_k(X(tau))-phi_k(X(0))].
//
// Q(tau) is linearly interpolated on the ADAPTIVE Garfield trajectory sample
// times. This is an auditable *test*, not proof of exact arbitrary-time response:
// RKF point spacing might be insufficient for the sharp earliest current.
// Every birth is shifted at exact continuous t_birth by Q(end)-Q(begin),
// never by averaging pre-binned current samples.
#include "Plane2100Archive.hh"

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
constexpr int NB = 10000;
constexpr double DT = 10.;
using B=NA60Plane::Archive::Birth;
using Sample=std::array<double,3>; // relative time ns; Qplus fC; Qminus fC
using W=std::array<std::vector<double>,2>;
std::string Arg(int n,char** v,const std::string& option) {
  for(int j=1;j+1<n;++j) if(option==v[j]) return v[j+1];
  return {};
}
W NewWave() {
  W v;
  for(auto& x:v)x.assign(NB,0.);
  return v;
}
void AddDirect(W& wave,const B& b, NA60Plane::Transport& transport) {
  for(int c=0;c<2;++c) for(int k=0;k<NB;++k)
    wave[c][k]+=b.weight*transport.Current(c,k);
}
std::vector<Sample> ReadTrajectoryCharge(NA60Plane::Transport& tr,
                const B& b, const NA60Plane::Direct& d,
                double& maxEarlyPathGapNs) {
  NA60Plane::Require(tr.DriftPathSize()>1,
      "RKF path has too few stored points for continuous Q(t) test");
  const double p0=tr.WeightingPotentialAt(b.u,b.v,b.w,0);
  const double m0=tr.WeightingPotentialAt(b.u,b.v,b.w,1);
  std::vector<Sample> curve;
  curve.reserve(tr.DriftPathSize()+2);
  curve.push_back({0.,0.,0.});
  double lastT=0.;
  for(std::size_t j=0;j<tr.DriftPathSize();++j) {
    double x=0,y=0,z=0,t=0;
    tr.DriftPathPoint(j,x,y,z,t);
    NA60Plane::Require(std::isfinite(x)&&std::isfinite(y)&&
                       std::isfinite(z)&&std::isfinite(t),
                       "Nonfinite RKF path point");
    NA60Plane::Require(t>=-1.e-9 &&
                       t<=d.endT+1.e-5*std::max(1.,d.endT),
                       "RKF path point outside drift endpoint duration");
    const double qp=NA60Plane::eFc*(tr.WeightingPotentialAt(x,y,z,0)-p0);
    const double qm=NA60Plane::eFc*(tr.WeightingPotentialAt(x,y,z,1)-m0);
    NA60Plane::Require(std::isfinite(qp)&&std::isfinite(qm),
                       "Nonfinite Q at RKF path point");
    if(t<lastT-1e-8)
      throw std::runtime_error("RKF drift path times are not monotonic");
    if(t<=lastT+1e-12) {
      if(j==0 && std::abs(t)<1e-12) continue; // explicit anchor already set
      NA60Plane::Require(std::abs(qp-curve.back()[1])<1e-8 &&
                         std::abs(qm-curve.back()[2])<1e-8,
                         "Same-time RKF path points have distinct Q");
      continue;
    }
    if(lastT<200.)
      maxEarlyPathGapNs=std::max(maxEarlyPathGapNs,
                                 std::min(t,200.)-lastT);
    curve.push_back({t,qp,qm});
    lastT=t;
  }
  // Endpoint Ramo values are recorded by the SAME RKF run, and their
  // readout convention already passed the direct ion transport audit.
  if(d.endT>curve.back()[0]+1e-10) {
    curve.push_back({d.endT,d.endpointQ[0],d.endpointQ[1]});
  } else {
    NA60Plane::Require(std::abs(d.endT-curve.back()[0])<
      1e-5*std::max(1.,d.endT),
      "Last stored RKF trajectory point and endpoint time disagree");
    curve.back()={d.endT,d.endpointQ[0],d.endpointQ[1]};
  }
  NA60Plane::Require(curve.size()>2,"No finite-time RKF points");
  NA60Plane::Require(curve.back()[0]>0. &&
      std::abs(curve.back()[1]-d.endpointQ[0])<1e-12 &&
      std::abs(curve.back()[2]-d.endpointQ[1])<1e-12,
      "Continuous trajectory endpoint does not obey Ramo charge");
  return curve;
}
// Mathematically exact *endpoint differencing* once Q(t) is supplied;
// interpolation of Q(t) between RKF trajectory points remains approximate.
std::array<double,2> At(const std::vector<Sample>& x,double relativeNs) {
  if(relativeNs<=0.)return {};
  if(relativeNs>=x.back()[0])return {x.back()[1],x.back()[2]};
  const auto upper=std::upper_bound(x.begin(),x.end(),relativeNs,
      [](double t,const Sample& s){return t<s[0];});
  NA60Plane::Require(upper!=x.begin()&&upper!=x.end(),
                     "Bad adaptive Q curve interpolation range");
  const auto& b=*upper;
  const auto& a=*(upper-1);
  const double f=(relativeNs-a[0])/(b[0]-a[0]);
  NA60Plane::Require(f>=-1e-10&&f<=1.+1e-10,
                     "Numerically invalid Q path interpolation fraction");
  return {a[1]+f*(b[1]-a[1]),a[2]+f*(b[2]-a[2])};
}
void DepositRamo(W& wave,const B& b,const std::vector<Sample>& q) {
  // Output is the same Garfield Sensor 10ns bin average CURRENT fC/ns.
  const int begin=std::max(0,(int)std::floor(b.t/DT));
  const int end=std::min(NB,2+(int)std::ceil((b.t+q.back()[0])/DT));
  for(int j=begin;j<end;++j) {
    const auto l=At(q,j*DT-b.t);
    const auto h=At(q,(j+1)*DT-b.t);
    for(int k=0;k<2;++k)
      wave[k][j]+=b.weight*(h[k]-l[k])/DT;
  }
}
W ReadReference(const std::string& path){
  std::ifstream file(path);
  NA60Plane::Require(bool(file),"Cannot open frozen reference "+path);
  std::string line;
  NA60Plane::Require(bool(std::getline(file,line))&&
     line=="time_center_ns,ion_plus_fC_per_ns,ion_minus_fC_per_ns,"
           "cumulative_plus_fC,cumulative_minus_fC",
     "Unexpected original direct-RKF reference header");
  W result=NewWave();
  for(int j=0;j<NB;++j){
    NA60Plane::Require(bool(std::getline(file,line)),
                       "Direct-RKF reference truncated before 100us");
    for(char& ch:line)if(ch==',')ch=' ';
    std::istringstream s(line);
    double t,p,m,qp,qm;
    NA60Plane::Require(bool(s>>t>>p>>m>>qp>>qm) &&
                       std::isfinite(p)&&std::isfinite(m)&&
                       std::abs(t-(j+.5)*DT)<1e-8,
                       "Malformed direct-RKF reference row/time");
    result[0][j]=p;result[1][j]=m;
  }
  while(std::getline(file,line))
    NA60Plane::Require(line.empty(),"Extra direct-RKF CSV rows");
  return result;
}
void ValidateSourceSummary(const std::string&fileName,
                           const NA60Plane::Archive& ion) {
  std::ifstream file(fileName);
  std::string header,line;
  NA60Plane::Require(bool(file)&&bool(std::getline(file,header))&&
                     bool(std::getline(file,line)),
                     "Missing original RKF summary");
  std::replace(line.begin(),line.end(),',',' ');
  std::istringstream s(line);
  std::vector<std::string> f;
  std::string word;
  while(s>>word)f.push_back(word);
  NA60Plane::Require(f.size()>=10 && f[0]==ion.id &&
                     f[2]==std::to_string(ion.station) &&
                     f[3]==std::to_string(ion.det) &&
                     f[4]=="weighted" && f[5]=="direct" &&
                     std::stoll(f[6])==ion.nTrue &&
                     std::stoll(f[7])==ion.nSaved &&
                     std::abs(std::stod(f[8])-DT)<1e-8 &&
                     std::abs(std::stod(f[9])-100.)<1e-8,
                     "Original RKF summary incompatible with K256 source");
}
struct Metric {
  double L2=0, qp=0,qm=0,absRms=0;
};
Metric Compare(const W& calculated,const W& direct,int lo,int hi) {
  double se=0,sr=0,qp=0,qm=0;
  for(int j=lo;j<hi;++j)for(int k=0;k<2;++k){
    const double d=calculated[k][j]-direct[k][j];
    se+=d*d;sr+=direct[k][j]*direct[k][j];
    if(k==0)qp+=DT*d;
    else qm+=DT*d;
  }
  return {std::sqrt(se/std::max(sr,1e-50)),qp,qm,
          std::sqrt(se/(2.*(hi-lo)))};
}
void SaveWaveforms(const W& predicted,const W& direct,
                   const W& frozen,const std::string&filename) {
  std::ofstream csv(filename);
  NA60Plane::Require(bool(csv),"cannot create new waveform CSV "+filename);
  csv<<std::setprecision(17);
  csv<<"time_center_ns,path_Q_plus_fC_per_ns,path_Q_minus_fC_per_ns,"
        "recomputed_direct_plus_fC_per_ns,recomputed_direct_minus_fC_per_ns,"
        "frozen_direct_plus_fC_per_ns,frozen_direct_minus_fC_per_ns\n";
  for(int i=0;i<NB;++i)
    csv<<(i+.5)*DT<<","<<predicted[0][i]<<","<<predicted[1][i]
       <<","<<direct[0][i]<<","<<direct[1][i]
       <<","<<frozen[0][i]<<","<<frozen[1][i]<<"\n";
}
} // namespace

int main(int argc,char**argv){
  try {
    const std::string source=Arg(argc,argv,"--ion-file");
    const std::string frozen=Arg(argc,argv,"--reference-waveforms");
    const std::string frozenSummary=Arg(argc,argv,"--reference-summary");
    const std::string mobility=Arg(argc,argv,"--mobility");
    const std::string out=Arg(argc,argv,"--out-prefix");
    NA60Plane::Require(!source.empty()&&!frozen.empty()&&
                       !frozenSummary.empty()&&!mobility.empty()&&!out.empty(),
      "usage: --ion-file ORIGINAL_WEIGHTED_ROOT "
      "--reference-waveforms FROZEN_DIRECT_CSV --reference-summary SUMMARY_CSV "
      "--mobility IONMOBILITY --out-prefix FRESH_OUTPUT");
    NA60Plane::Require(!std::filesystem::exists(out+"_summary.csv")&&
                       !std::filesystem::exists(out+"_comparison.csv"),
                       "Refusing to overwrite an existing cumulative-Q audit");
    NA60Plane::Archive ions(source);
    NA60Plane::Require(ions.mode=="weighted" && ions.Entries()>0 &&
                       ions.Entries()<=20000,
                       "Requires original weighted reservoir <= 20000 rows");
    ValidateSourceSummary(frozenSummary,ions);
    const W original=ReadReference(frozen);
    NA60Plane::Transport drift(mobility,DT,NB,.05,true);
    W direct=NewWave(),continuousQ=NewWave();
    double maxEarlySegmentNs=0;
    long long totalPathPoints=0,late=0;
    std::array<double,2> weightedEndpointQ{};
    for(long long i=0;i<ions.Entries();++i) {
      const B b=ions.Get(i);
      drift.CheckStart(b.u,b.v,b.w,b.t);
      const auto actual=drift.Drift(b.u,b.v,b.w,b.t);
      AddDirect(direct,b,drift);
      const auto zero=drift.Drift(b.u,b.v,b.w,0.);
      if(zero.late)++late;
      weightedEndpointQ[0]+=b.weight*zero.endpointQ[0];
      weightedEndpointQ[1]+=b.weight*zero.endpointQ[1];
      auto q=ReadTrajectoryCharge(drift,b,zero,maxEarlySegmentNs);
      totalPathPoints+=(long long)q.size();
      DepositRamo(continuousQ,b,q);
      if((i+1)%1000==0)
        std::cout<<"CUMULATIVE-Q rows="<<i+1<<"/"<<ions.Entries()<<"\n"
                 <<std::flush;
    }
    ions.CheckComplete();
    const auto directIdentity=Compare(direct,original,0,NB);
    std::cout<<std::setprecision(17)
             <<"INTEGRITY recomputed_direct_vs_frozen_L2="
             <<directIdentity.L2<<"\n";
    NA60Plane::Require(directIdentity.L2<1e-5,
                       "Original real-birth-time RKF reference mismatch; stop");
    const auto whole=Compare(continuousQ,direct,0,NB);
    const auto early=Compare(continuousQ,direct,0,10);
    const auto mid=Compare(continuousQ,direct,10,100);
    const auto mid2=Compare(continuousQ,direct,100,1000);
    const auto lateMetric=Compare(continuousQ,direct,1000,NB);
    SaveWaveforms(continuousQ,direct,original,out+"_comparison.csv");
    std::ofstream csv(out+"_summary.csv");
    NA60Plane::Require(bool(csv),"Cannot open cumulative Q summary");
    csv<<std::setprecision(17);
    csv<<"crossing_id,original_weighted_macroions,original_exact_ion_births,"
          "recomputed_direct_frozen_L2,continuous_Q_full_L2,"
          "continuous_Q_0to100ns_L2,continuous_Q_100ns_to1us_L2,"
          "continuous_Q_1us_to10us_L2,continuous_Q_10us_to100us_L2,"
          "continuous_Q_first100plus_error_fC,"
          "continuous_Q_first100minus_error_fC,"
          "continuous_Q_100us_plus_error_fC,"
          "continuous_Q_100us_minus_error_fC,"
          "late_ion_rows,drift_path_point_count,largest_first200ns_path_gap_ns,"
          "weighted_full_endpoint_plus_fC,weighted_full_endpoint_minus_fC\n";
    csv<<ions.id<<","<<ions.Entries()<<","<<ions.nTrue<<","
       <<directIdentity.L2<<","<<whole.L2<<","<<early.L2<<","
       <<mid.L2<<","<<mid2.L2<<","<<lateMetric.L2<<","
       <<early.qp<<","<<early.qm<<","<<whole.qp<<","<<whole.qm<<","
       <<late<<","<<totalPathPoints<<","<<maxEarlySegmentNs<<","
       <<weightedEndpointQ[0]<<","<<weightedEndpointQ[1]<<"\n";
    csv.close();
    std::cout<<"PASS exact-spatial cumulative Ramo charge endpoint test\n"
             <<"relative_L2_full="<<whole.L2
             <<" first100ns="<<early.L2
             <<" 100ns_to_1us="<<mid.L2
             <<" 1us_to_10us="<<mid2.L2
             <<" 10us_to_100us="<<lateMetric.L2<<"\n"
             <<"100ns_errors_fC="<<early.qp<<","<<early.qm
             <<" 100us_errors_fC="<<whole.qp<<","<<whole.qm<<"\n"
             <<"RKF_points="<<totalPathPoints
             <<" worst_initial_200ns_path_gap_ns="<<maxEarlySegmentNs
             <<" late_ion_rows="<<late<<"\n"
             <<"NOTE: Q(t) between RKF drift points is time-linearly interpolated;"
                " errors are empirical, NOT mathematically guaranteed.\n";
    return 0;
  } catch(const std::exception& error) {
    std::cerr<<"EXACT-SPATIAL CUMULATIVE Q FAILED: "<<error.what()<<"\n";
    return 2;
  }
}
