// Full read-only RKF vertex export, with exact statistical census verification.
// Native little-endian float64 vertices: absolute ns, u mm, v mm, U w0 mm, V w0 mm.
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



}
int main(int argc,char** argv) {
 try {
  Require(argc==4,"usage: PROGRAM frozen.root MOBILITY NEW_DIR");
  namespace fs=std::filesystem;
  const fs::path out=argv[3];
  Require(!fs::exists(out),"output exists");
  NA60Plane::Archive archive(argv[1]);
  ChamberFrame frame(argv[1]);
  Require(archive.mode=="weighted"&&archive.Entries()>0,"need complete weighted archive");
  Require(fs::create_directories(out),"cannot create output directory");
  PhysicalIonPaths transport(argv[2],.005);
  std::ofstream data(out/"vertices.f64",std::ios::binary);
  std::ofstream meta(out/"ions.csv");
  Require(bool(data)&&bool(meta),"cannot create output files");
  static_assert(sizeof(double)==8,"requires binary64");
  meta<<std::setprecision(17)
      <<"row,primary,weight,birth_ns,end_ns,first_vertex,vertices,upper_index,lower_index\n";
  long long offset=0;
  double sum=0.;
  for(long long j=0;j<archive.Entries();++j) {
    const auto b=archive.Get(j);
    frame.CheckBirth(j,b);
    auto path=transport.Run(b);
    Require(path.points.size()>=2,"invalid path");
    // Snap only the final normal coordinate to the ideal electrode.
    path.points.back().v=std::copysign(NA60Plane::gap,path.end.v);
    meta<<j<<","<<b.primary<<","<<b.weight<<","<<b.t<<","<<path.end.t<<","
        <<offset<<","<<path.points.size()<<","<<frame.uCentralIndex<<","
        <<frame.vCentralIndex<<"\n";
    for(const auto& p:path.points) {
      Require(std::abs(p.w-b.w)<1e-8,"unexpected longitudinal motion");
      double wu=10.*(p.w+frame.wEntranceCm-.1*frame.patchShiftCm-.17*frame.uCentralIndex);
      double wv=10.*(p.w+frame.wEntranceCm+.1*frame.patchShiftCm-.17*frame.vCentralIndex);
      const double row[5]={p.t,10*p.u,10*p.v,wu,wv};
      data.write(reinterpret_cast<const char*>(row),sizeof(row));
    }
    offset+=path.points.size(); sum+=b.weight;
    if((j+1)%512==0) std::cout<<"RKF rows "<<j+1<<"/"<<archive.Entries()<<std::endl;
  }
  archive.CheckComplete();
  data.close();meta.close();
  Require(bool(data)&&bool(meta),"output I/O failed");
  Require(fs::file_size(out/"vertices.f64")==static_cast<unsigned long long>(offset)*40,
          "binary size mismatch");
  std::ofstream seal(out/"complete.txt");
  seal<<std::setprecision(17)<<"rows="<<archive.Entries()<<"\n"
      <<"exact_ions="<<archive.nTrue<<"\nweight_sum="<<sum
      <<"\nvertices="<<offset<<"\nstep_mm=0.005\n";
  seal.close();Require(bool(seal),"seal write failed");
  std::cout<<"PASS archive census rows="<<archive.Entries()<<" exact_ions="
           <<archive.nTrue<<" vertices="<<offset<<"\n";
  return 0;
 }catch(const std::exception& e){std::cerr<<e.what()<<"\n";return 2;}
}
