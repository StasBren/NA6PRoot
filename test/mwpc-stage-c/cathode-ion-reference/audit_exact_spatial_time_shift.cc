// Exact-spatial ion replay identity test: same real archived K256 ion starts,
// Garfield++ DriftLineRKF direct at original t_birth vs at t=0 followed by
// fractional 10ns linear time translation (the kernel-replay approximation).
//
// This diagnoses the *temporal discretization/replay floor* without spatial
// interpolation. No new Heed/AvalancheMicroscopic, FEM/U-V, or electronics.
// In the 2100V symmetric P3 whole-plane ion model only.
#include "Plane2100Archive.hh"
#include <array>
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
constexpr int NB=10000;
constexpr double DT=10.0;
using W=std::array<std::vector<double>,2>;
using Birth=NA60Plane::Archive::Birth;

std::string Arg(int n,char**v,const std::string&key) {
  for(int j=1;j+1<n;++j)if(key==v[j])return v[j+1];
  return {};
}
W NewWave() {
  W s;
  for(auto& c:s)c.assign(NB,0.);
  return s;
}
int IntArg(int n,char**v,const std::string&key,int defaultValue) {
  const auto t=Arg(n,v,key);
  return t.empty()?defaultValue:std::stoi(t);
}
void AccumulateDirect(W& total, NA60Plane::Transport& tr, const Birth& b) {
  for(int channel=0;channel<2;++channel)
    for(int k=0;k<NB;++k)
      total[channel][k]+=b.weight*tr.Current(channel,k);
}
// Deposit a signal computed by DriftIon(b.u,b.v,b.w,LOCAL_PHASE_NS)
 // into the integer-shifted original 10ns window.
 // This is direct Garfield fractional-PHASE interpolation, not interpolation
 // of the 10ns histogram samples of a t=0 current.
void AccumulateIntegerShift(W& total,NA60Plane::Transport& tr,
                            const Birth& b,int offset,double factor) {
  NA60Plane::Require(offset>=0 && offset<=NB &&
                     std::isfinite(factor) && factor>=0. &&
                     factor<=1.+1e-10,
                     "Invalid phase interpolation coefficient/shift");
  if(factor==0)return;
  for(int channel=0;channel<2;++channel)
    for(int j=0;j+offset<NB;++j)
      total[channel][j+offset]+=factor*b.weight*tr.Current(channel,j);
}
void AccumulatePhaseInterpolated(W& total,NA60Plane::Transport& tr,
                                 const Birth& b,int nPhase) {
  NA60Plane::Require(nPhase>=1 && nPhase<=16,
                     "phase subdivision must be 1..16");
  const double shift=b.t/DT;
  const int offset=(int)std::floor(shift);
  const double f=shift-offset;
  NA60Plane::Require(offset>=0 && offset<NB && f>=0 && f<1,
                     "Original ion birth phase outside window");
  const double x=f*nPhase;
  const int k=std::min(nPhase-1,(int)std::floor(x));
  const double part=x-k;
  NA60Plane::Require(part>=0 && part<1.+1e-10,
                     "Phase bin has invalid fraction");

  // Required calling order: tr currently contains the t=0 signal.
  // First consume all phase endpoints supported by t=0 before any
  // additional DriftIon call clears its Sensor history.
  if(k==0)AccumulateIntegerShift(total,tr,b,offset,1.-part);
  if(k+1==nPhase)AccumulateIntegerShift(total,tr,b,offset+1,part);

  if(k>0 && part<1.) {
    const double phase=DT*k/nPhase;
    tr.Drift(b.u,b.v,b.w,phase);
    AccumulateIntegerShift(total,tr,b,offset,1.-part);
  }
  if(k+1<nPhase && part>0.) {
    const double phase=DT*(k+1)/nPhase;
    tr.Drift(b.u,b.v,b.w,phase);
    AccumulateIntegerShift(total,tr,b,offset,part);
  }
}
void AccumulateShifted(W& total,NA60Plane::Transport& tr,
                       const Birth& b) {
  const double shift=b.t/DT;
  const int offset=(int)std::floor(shift);
  const double f=shift-offset;
  NA60Plane::Require(std::isfinite(shift)&&offset>=0&&offset<NB&&
                     f>=0&&f<1,"invalid archive ion birth time for shift");
  for(int channel=0;channel<2;++channel) {
    for(int j=0;j<NB;++j) {
      const double c=b.weight*tr.Current(channel,j);
      if(c==0)continue;
      if(j+offset<NB)total[channel][j+offset]+=(1.-f)*c;
      if(j+offset+1<NB)total[channel][j+offset+1]+=f*c;
    }
  }
}
double RelativeL2(const W& x,const W& y,int lo,int hi) {
  double numer=0,denom=0;
  for(int c=0;c<2;++c)for(int t=lo;t<hi;++t) {
    numer+=(x[c][t]-y[c][t])*(x[c][t]-y[c][t]);
    denom+=y[c][t]*y[c][t];
  }
  return std::sqrt(numer/std::max(denom,1.e-50));
}
std::array<double,2> ChargeError(const W& x,const W& y,int lo,int hi){
  std::array<double,2> q{};
  for(int c=0;c<2;++c)for(int t=lo;t<hi;++t)
    q[c]+=DT*(x[c][t]-y[c][t]);
  return q;
}
W ReadOriginalReference(const std::string& path) {
  std::ifstream f(path);
  NA60Plane::Require(bool(f),"cannot open ORIGINAL reference CSV "+path);
  std::string line;
  NA60Plane::Require(bool(std::getline(f,line))&&
    line=="time_center_ns,ion_plus_fC_per_ns,ion_minus_fC_per_ns,"
          "cumulative_plus_fC,cumulative_minus_fC",
    "original direct CSV header unexpected");
  W out=NewWave();
  for(int j=0;j<NB;++j) {
    NA60Plane::Require(bool(std::getline(f,line)),"original direct CSV too short");
    for(char& ch:line)if(ch==',')ch=' ';
    std::istringstream s(line);
    double t,plus,minus,qplus,qminus;
    NA60Plane::Require(bool(s>>t>>plus>>minus>>qplus>>qminus),
                       "bad original CSV row");
    NA60Plane::Require(std::abs(t-(j+.5)*DT)<1e-7 &&
                       std::isfinite(plus)&&std::isfinite(minus),
                       "unexpected original time grid or current");
    out[0][j]=plus;
    out[1][j]=minus;
  }
  NA60Plane::Require(!std::getline(f,line) || line.empty(),
                     "extra original waveform lines");
  return out;
}
void ValidateOriginalSummary(const std::string&path,
                             const NA60Plane::Archive& a) {
  std::ifstream f(path);
  NA60Plane::Require(bool(f),"cannot open reference summary");
  std::string h,line;
  NA60Plane::Require(bool(std::getline(f,h)) && bool(std::getline(f,line)),
                     "reference summary empty");
  std::istringstream s(line);
  std::vector<std::string> v;
  std::string cell;
  while(std::getline(s,cell,','))v.push_back(cell);
  NA60Plane::Require(v.size()>=10 && v[0]==a.id &&
                       v[2]==std::to_string(a.station) &&
                       v[3]==std::to_string(a.det) &&
                       v[4]=="weighted"&&v[5]=="direct"&&
                       std::stoll(v[6])==a.nTrue&&std::stoll(v[7])==a.nSaved &&
                       std::abs(std::stod(v[8])-DT)<1e-6&&
                       std::abs(std::stod(v[9])-100.)<1e-6,
                       "original direct summary does not match archive/time grid");
}
} // namespace

int main(int argc,char**argv) {
  try {
    const auto source=Arg(argc,argv,"--ion-file");
    const auto original=Arg(argc,argv,"--reference-waveforms");
    const auto originalSummary=Arg(argc,argv,"--reference-summary");
    const auto mobility=Arg(argc,argv,"--mobility");
    const auto prefix=Arg(argc,argv,"--out-prefix");
    const int phaseSubdivisions=IntArg(argc,argv,"--phase-subdivisions",1);
    NA60Plane::Require(phaseSubdivisions==1 || phaseSubdivisions==2 ||
                       phaseSubdivisions==4 || phaseSubdivisions==8,
                       "phase subdivisions must be one of 1,2,4,8");
    NA60Plane::Require(!source.empty()&&!original.empty()&&
                       !originalSummary.empty()&&!mobility.empty()&&
                       !prefix.empty(),
      "usage: --ion-file ROOT --reference-waveforms CSV "
      "--reference-summary CSV --mobility FILE --out-prefix NEW_PATH");
    NA60Plane::Require(!std::filesystem::exists(prefix+"_comparison.csv") &&
                       !std::filesystem::exists(prefix+"_summary.csv"),
                       "refuse to overwrite exact-spatial replay output");
    NA60Plane::Archive ions(source);
    NA60Plane::Require(ions.mode=="weighted"&&ions.Entries()>0&&
                       ions.Entries()<=20000,
                       "requires complete weighted <=20000 original macroion rows");
    ValidateOriginalSummary(originalSummary,ions);
    const W ref=ReadOriginalReference(original);
    NA60Plane::Transport tr(mobility,DT,NB,.05,true);
    W direct=NewWave(),shifted=NewWave(),phaseCorrected=NewWave();
    int late0=0,lateActual=0;
    double sumWeights=0;
    for(long long i=0;i<ions.Entries();++i) {
      const auto b=ions.Get(i);
      tr.CheckStart(b.u,b.v,b.w,b.t);
      const auto referenceDrift=tr.Drift(b.u,b.v,b.w,b.t);
      if(referenceDrift.late)++lateActual;
      AccumulateDirect(direct,tr,b);
      const auto zeroDrift=tr.Drift(b.u,b.v,b.w,0.);
      if(zeroDrift.late)++late0;
      AccumulateShifted(shifted,tr,b);
      // The zero-time ion waveform is still in the sensor.
      AccumulatePhaseInterpolated(phaseCorrected,tr,b,phaseSubdivisions);
      NA60Plane::Require(std::abs((referenceDrift.endT-b.t)-
                                   zeroDrift.endT)<2.,
                        "Unexpected time-dependence in static-field RKF drift");
      sumWeights+=b.weight;
      if((i+1)%1000==0)
        std::cout<<"EXACT-SPATIAL rows="<<i+1<<"/"<<ions.Entries()<<"\n";
    }
    ions.CheckComplete();
    NA60Plane::Require(std::abs(sumWeights-ions.nTrue)<
                       1e-6*std::max<long long>(1,ions.nTrue),
                       "weighted census lost");
    const double orig=RelativeL2(direct,ref,0,NB);
    std::cout<<std::setprecision(16)
             <<"INTEGRITY direct_recomputed_vs_frozen L2="<<orig<<"\n";
    NA60Plane::Require(orig<1e-5,
                       "recomputed same-ion direct RKF differs from frozen benchmark; "
                       "STOP: check original mobility, RKF step and time grid");
    const double full=RelativeL2(shifted,direct,0,NB);
    const double phaseFull=RelativeL2(phaseCorrected,direct,0,NB);
    const auto phaseEarly=ChargeError(phaseCorrected,direct,0,10);
    const auto phaseWindow=ChargeError(phaseCorrected,direct,0,NB);
    const double phaseWindows[4]={
        RelativeL2(phaseCorrected,direct,0,10),
        RelativeL2(phaseCorrected,direct,10,100),
        RelativeL2(phaseCorrected,direct,100,1000),
        RelativeL2(phaseCorrected,direct,1000,NB)
    };
    if(phaseSubdivisions==1) {
      NA60Plane::Require(RelativeL2(phaseCorrected,shifted,0,NB)<1e-11,
        "nPhase=1 does not reproduce original fractional t=0 time shift");
    }
    const auto early=ChargeError(shifted,direct,0,10);
    const auto all=ChargeError(shifted,direct,0,NB);
    const double windows[4]={
        RelativeL2(shifted,direct,0,10),
        RelativeL2(shifted,direct,10,100),
        RelativeL2(shifted,direct,100,1000),
        RelativeL2(shifted,direct,1000,NB)
    };
    std::ofstream wave(prefix+"_comparison.csv");
    NA60Plane::Require(bool(wave),"cannot open output comparison CSV");
    wave<<std::setprecision(17);
    wave<<"time_center_ns,recomputed_direct_plus_fC_per_ns,"
           "recomputed_direct_minus_fC_per_ns,"
           "perfect_spatial_shifted_plus_fC_per_ns,"
           "perfect_spatial_shifted_minus_fC_per_ns,"
           "phase_corrected_plus_fC_per_ns,phase_corrected_minus_fC_per_ns,"
           "frozen_direct_plus_fC_per_ns,frozen_direct_minus_fC_per_ns\n";
    for(int j=0;j<NB;++j)
      wave<<(j+.5)*DT<<","<<direct[0][j]<<","<<direct[1][j]<<","
          <<shifted[0][j]<<","<<shifted[1][j]<<","
          <<phaseCorrected[0][j]<<","<<phaseCorrected[1][j]<<","
          <<ref[0][j]<<","<<ref[1][j]<<"\n";
    wave.close();
    std::ofstream summary(prefix+"_summary.csv");
    NA60Plane::Require(bool(summary),"cannot open output summary CSV");
    summary<<std::setprecision(16);
    summary<<"crossing_id,original_macroions,original_exact_ions,"
              "recomputed_direct_vs_frozen_L2,perfect_spatial_shifted_vs_direct_L2,"
              "first100ns_L2,from100ns_to1us_L2,from1us_to10us_L2,"
              "from10us_to100us_L2,"
              "early100plus_error_fC,early100minus_error_fC,"
              "window100us_plus_error_fC,window100us_minus_error_fC,"
              "n_late_direct,n_late_zero_start,"
              "phase_subdivisions,phase_corrected_full_L2,"
              "phase_corrected_first100ns_L2,"
              "phase_corrected_100ns_to1us_L2,"
              "phase_corrected_1us_to10us_L2,"
              "phase_corrected_10us_to100us_L2,"
              "phase_corrected_early100plus_error_fC,"
              "phase_corrected_early100minus_error_fC,"
              "phase_corrected_100usplus_error_fC,"
              "phase_corrected_100usminus_error_fC\n";
    summary<<ions.id<<","<<ions.Entries()<<","<<ions.nTrue<<","
           <<orig<<","<<full<<","<<windows[0]<<","<<windows[1]<<","
           <<windows[2]<<","<<windows[3]<<","
           <<early[0]<<","<<early[1]<<","<<all[0]<<","<<all[1]<<","
           <<lateActual<<","<<late0<<","
           <<phaseSubdivisions<<","<<phaseFull<<","
           <<phaseWindows[0]<<","<<phaseWindows[1]<<","
           <<phaseWindows[2]<<","<<phaseWindows[3]<<","
           <<phaseEarly[0]<<","<<phaseEarly[1]<<","
           <<phaseWindow[0]<<","<<phaseWindow[1]<<"\n";
    summary.close();
    std::cout<<"PASS perfect-spatial 10ns shift control; same original K256 reservoir\n"
             <<"combined_time_shift_L2="<<full
             <<" first100ns_L2="<<windows[0]
             <<" 100ns-1us_L2="<<windows[1]
             <<" 1us-10us_L2="<<windows[2]
             <<" 10us-100us_L2="<<windows[3]<<"\n"
             <<"early100_err_fC="<<early[0]<<","<<early[1]
             <<" window100us_err_fC="<<all[0]<<","<<all[1]<<"\n"
             <<"PHASE_CONTROL n_subdivisions="<<phaseSubdivisions
             <<" phase_corrected_full_L2="<<phaseFull
             <<" phase_0to100ns_L2="<<phaseWindows[0]
             <<" phase_100ns_to1us_L2="<<phaseWindows[1]
             <<" phase_1us_to10us_L2="<<phaseWindows[2]
             <<" phase_10us_to100us_L2="<<phaseWindows[3]<<"\n"
             <<"PHASE early100_error_fC="<<phaseEarly[0]<<","<<phaseEarly[1]
             <<" window100us_error_fC="<<phaseWindow[0]<<","<<phaseWindow[1]<<"\n"
             <<"This is the exact-SPATIAL replay discrepancy, NOT spatial interpolation\n";
    return 0;
  } catch(const std::exception& ex) {
    std::cerr<<"EXACT-SPATIAL CONTROL FAILED: "<<ex.what()<<"\n";
    return 2;
  }
}
