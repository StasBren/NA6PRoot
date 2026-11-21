// Read-only developmental archived-ion trajectory sample for spectral U/V.
// Export fixed positions at ion ages, not a full muon waveform.
#include "Plane2100Archive.hh"
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>
namespace {
using NA60Plane::Require;
struct Point {double u,v,w,t;}; // cm, ns
class PhysicalIonPaths {
  Garfield::MediumMagboltz gas_;
  Garfield::ComponentAnalyticField field_;
  Garfield::Sensor sensor_;
  Garfield::DriftLineRKF drift_;
 public:
  explicit PhysicalIonPaths(const std::string& mobility,
                            double stepMm) : drift_(&sensor_) {
    Require(stepMm>0. && stepMm<=.05, "step-mm must be in (0,0.05]");
    gas_.SetComposition("ar",70.,"co2",30.);
    gas_.SetTemperature(293.15);
    gas_.SetPressure(760.);
    gas_.Initialise(false);
    Require(gas_.LoadIonMobility(mobility), "ion mobility load failed");
    field_.SetMedium(&gas_);
    for (int k=-4;k<=4;++k)
      field_.AddWire(k*NA60Plane::pitch,0.,.003,2100.,"anode");
    field_.AddPlaneY(-NA60Plane::gap,0.,"cathode_minus");
    field_.AddPlaneY(+NA60Plane::gap,0.,"cathode_plus");
    field_.SetMagneticField(0.,0.,0.);
    sensor_.AddComponent(&field_);
    sensor_.SetArea(-1.35,-NA60Plane::gap,-5.,
                     1.35,+NA60Plane::gap,+5.);
    drift_.EnableSignalCalculation(false);
    // Exactly the same limit and analytic field model as Plane2100Response.hh.
    drift_.SetMaximumStepSize(0.1*stepMm);
  }
  struct Result {
    std::vector<Point> points;
    Point end;
    int status=0;
  };
  Result Run(const NA60Plane::Archive::Birth& b) {
    Require(std::isfinite(b.u)&&std::isfinite(b.v)&&std::isfinite(b.w)&&
            std::isfinite(b.t)&&b.t>=0.,"invalid ion birth coordinates/time");
    const double wire=std::round(b.u/NA60Plane::pitch)*NA60Plane::pitch;
    Require(std::hypot(b.u-wire,b.v)>NA60Plane::wireRadius &&
            std::abs(b.u)<1.35 && std::abs(b.v)<NA60Plane::gap &&
            std::abs(b.w)<5 && b.v!=0.,"ion birth outside physical gas");
    Require(drift_.DriftIon(b.u,b.v,b.w,b.t), "ion RKF drift failed");
    Result out;
    drift_.GetEndPoint(out.end.u,out.end.v,out.end.w,
                       out.end.t,out.status);
    Require(std::isfinite(out.end.t)&&out.end.t>b.t &&
            std::abs(std::abs(out.end.v)-NA60Plane::gap)<.001 &&
            (out.end.v>0.)==(b.v>0.),
            "ion endpoint is not the expected physical cathode");
    const auto n=drift_.GetNumberOfDriftLinePoints();
    Require(n>=1 && n<10000000, "invalid RKF drift-line point count");
    out.points.reserve(n+2);
    for (std::size_t j=0;j<n;++j) {
      Point p{};
      drift_.GetDriftLinePoint(j,p.u,p.v,p.w,p.t);
      Require(std::isfinite(p.u)&&std::isfinite(p.v)&&
              std::isfinite(p.w)&&std::isfinite(p.t),
              "nonfinite RKF path vertex");
      if (!out.points.empty()) {
        Require(p.t>=out.points.back().t,"nonmonotonic RKF times");
        if (p.t==out.points.back().t) continue; // terminal duplicate
      }
      out.points.push_back(p);
    }
    const Point start{b.u,b.v,b.w,b.t};
    Require(!out.points.empty() && out.points.front().t>=b.t-1.e-6,
            "RKF path begins before ion birth");
    if (out.points.front().t>b.t+1.e-7) out.points.insert(out.points.begin(),start);
    else out.points.front()=start;
    Require(out.points.back().t<=out.end.t+1.e-6,
            "RKF path extends beyond reported endpoint");
    if (out.points.back().t<out.end.t-1.e-7)
      out.points.push_back(out.end);
    else out.points.back()=out.end;
    Require(out.points.size()>=2 && out.points.back().t>out.points.front().t,
            "ion trajectory too short");
    return out;
  }
};

// Two frames coexist in frozen ROOT: local Garfield patch coordinates and
// physical chamber coordinates. Both MUST be retained for stereo-strip phase.
class ChamberFrame {
  std::unique_ptr<TFile> file_;
  TTree* ions_ = nullptr;
  double uChamber_ = 0., wChamber_ = 0.;
 public:
  double patchShiftCm=0.,uEntranceCm=0.,wEntranceCm=0.;
  int uCentralIndex=0,vCentralIndex=0;
  explicit ChamberFrame(const std::string& filename) {
    file_.reset(TFile::Open(filename.c_str(),"READ"));
    Require(file_&&!file_->IsZombie(),"cannot reopen physical position ROOT metadata");
    auto* mu=dynamic_cast<TTree*>(file_->Get("muons"));
    ions_=dynamic_cast<TTree*>(file_->Get("ions"));
    Require(mu&&ions_&&mu->GetEntries()==1,
            "missing geometry metadata/ions in frozen ROOT archive");
    for(const char* s:{"patch_shift_cm","u_in_mm","w_in_mm"})
      Require(mu->GetBranch(s)!=nullptr,
              "missing required muon physical-frame branch "+std::string(s));
    for(const char* s:{"u_chamber_cm","w_chamber_cm"})
      Require(ions_->GetBranch(s)!=nullptr,
              "missing required ion physical-frame branch "+std::string(s));
    double uEntranceMm=0.,wEntranceMm=0.;
    mu->SetBranchAddress("patch_shift_cm",&patchShiftCm);
    mu->SetBranchAddress("u_in_mm",&uEntranceMm);
    mu->SetBranchAddress("w_in_mm",&wEntranceMm);
    mu->GetEntry(0);
    uEntranceCm=.1*uEntranceMm;
    wEntranceCm=.1*wEntranceMm;
    Require(std::isfinite(patchShiftCm)&&std::isfinite(uEntranceCm)&&
            std::isfinite(wEntranceCm),"nonfinite physical-frame metadata");
    uCentralIndex=int(std::llround((wEntranceCm-.1*uEntranceCm)/.17));
    vCentralIndex=int(std::llround((wEntranceCm+.1*uEntranceCm)/.17));
    ions_->SetBranchAddress("u_chamber_cm",&uChamber_);
    ions_->SetBranchAddress("w_chamber_cm",&wChamber_);
  }
  void CheckBirth(long long j,const NA60Plane::Archive::Birth& b) {
    Require(ions_->GetEntry(j)>0,"physical ion location row is absent");
    Require(std::isfinite(uChamber_)&&std::isfinite(wChamber_) &&
            std::abs(uChamber_-(b.u+patchShiftCm))<1.e-8 &&
            std::abs(wChamber_-(b.w+wEntranceCm))<1.e-8,
            "local/physical ion-coordinate mapping failed");
  }
};


Point AtAge(const PhysicalIonPaths::Result& path,double t) {
  if(t>=path.end.t) {
    auto p=path.end;
    p.v=std::copysign(NA60Plane::gap,p.v); // exact ideal cathode
    p.t=t;
    return p;
  }
  const auto& p=path.points;
  auto it=std::lower_bound(p.begin(),p.end(),t,
    [](const Point& q,double a){return q.t<a;});
  Require(it!=p.end(),"cannot bracket age");
  if(it==p.begin()||it->t==t) return *it;
  const auto& l=*(it-1); const auto& r=*it;
  double f=(t-l.t)/(r.t-l.t);
  return {l.u+f*(r.u-l.u),l.v+f*(r.v-l.v),l.w+f*(r.w-l.w),t};
}
}
int main(int argc,char** argv) {
 try {
  Require(argc==4,"usage: PROGRAM frozen_ion_births.root MOBILITY new.csv");
  Require(!std::filesystem::exists(argv[3]),"output exists");
  NA60Plane::Archive archive(argv[1]);
  ChamberFrame frame(argv[1]);
  Require(archive.mode=="weighted"&&archive.Entries()>=16,"need weighted archive >=16 rows");
  PhysicalIonPaths transport(argv[2],.005);
  std::vector<double> ages;
  for(int i=0;i<=10;++i) ages.push_back(10.*i);
  ages.insert(ages.end(),{1000.,10000.,100000.});
  std::ofstream out(argv[3]);
  Require(bool(out),"cannot open output");
  out<<std::setprecision(17)
     <<"row,primary,weight,birth_ns,age_ns,u_mm,v_mm,w_mm,upper_w0_mm,lower_w0_mm,"
       "upper_index,lower_index,end_age_ns,at_cathode\n";
  for(int i=0;i<16;++i) {
    long long j=(2*i+1)*archive.Entries()/32;
    auto b=archive.Get(j);
    frame.CheckBirth(j,b);
    auto path=transport.Run(b);
    for(double age:ages) {
      auto p=AtAge(path,b.t+age);
      double wu=10.*(p.w+frame.wEntranceCm-.1*frame.patchShiftCm-.17*frame.uCentralIndex);
      double wv=10.*(p.w+frame.wEntranceCm+.1*frame.patchShiftCm-.17*frame.vCentralIndex);
      out<<j<<","<<b.primary<<","<<b.weight<<","<<b.t<<","<<age<<","
         <<10*p.u<<","<<10*p.v<<","<<10*p.w<<","<<wu<<","<<wv<<","
         <<frame.uCentralIndex<<","<<frame.vCentralIndex<<","<<path.end.t-b.t<<","
         <<int(b.t+age>=path.end.t)<<"\n";
    }
  }
  out.close(); Require(bool(out),"CSV write failed");
  std::cout<<"PASS: 16 deterministic archive rows, 14 ages each; "
           <<"not a complete archive or absolute-time waveform\n";
  return 0;
 }catch(const std::exception& e){std::cerr<<e.what()<<"\n";return 2;}
}
