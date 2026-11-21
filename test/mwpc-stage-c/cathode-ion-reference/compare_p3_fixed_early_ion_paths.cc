// Three-wire V5/V6 weighting-only cumulative-Ramo comparison.
// Fixed analytic 9-wire trajectories, identical coordinates for both maps.
// No archived ions, avalanches, signals from electrons or production claims.
#include "Plane2100Archive.hh"
#include "Garfield/ComponentElmer.hh"
#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
namespace {
namespace fs=std::filesystem;
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


struct Trial {
  double radiusUm, theta, wMm;
  std::array<Point,11> samples;
};
Point Sample(const std::vector<Point>& p,double t) {
  Require(t>=p.front().t && t<=p.back().t,"time outside RKF path");
  auto hi=std::lower_bound(p.begin(),p.end(),t,
     [](const Point& q,double x){return q.t<x;});
  Require(hi!=p.end(),"missing path bracket");
  if(hi==p.begin() || hi->t==t) return *hi;
  const auto &a=*(hi-1), &b=*hi;
  const double f=(t-a.t)/(b.t-a.t);
  return {a.u+f*(b.u-a.u),a.v+f*(b.v-a.v),a.w+f*(b.w-a.w),t};
}
void Evaluate(const fs::path& dir,const fs::path& out,
              const std::vector<Trial>& trials) {
  Require(!fs::exists(out),"output exists: "+out.string());
  for(const char* f:{"chamber/mesh.header","chamber/mesh.elements",
       "chamber/mesh.nodes","dielectrics.dat","chamber/drift.result",
       "chamber/weight_plus_0.result","chamber/weight_minus_0.result"})
    Require(fs::is_regular_file(dir/f)&&fs::file_size(dir/f)>0,
            "missing FEM file "+(dir/f).string());
  Garfield::ComponentElmer map((dir/"chamber/mesh.header").string(),
          (dir/"chamber/mesh.elements").string(),(dir/"chamber/mesh.nodes").string(),
          (dir/"dielectrics.dat").string(),(dir/"chamber/drift.result").string(),"mm");
  Require(map.SetWeightingField((dir/"chamber/weight_plus_0.result").string(),"U"),
          "U import failed");
  Require(map.SetWeightingField((dir/"chamber/weight_minus_0.result").string(),"V"),
          "V import failed");
  std::ofstream csv(out);
  Require(bool(csv),"cannot write output");
  csv<<std::setprecision(16)
     <<"trial,radius_um,theta_rad,w_mm,side,t_ns,u_cm,v_cm,w_cm,phi,Q_over_e,bin_Q_over_e\n";
  for(std::size_t j=0;j<trials.size();++j) {
    const auto& trial=trials[j];
    for(const char* label:{"U","V"}) {
      double initial=0.,prev=0.;
      for(std::size_t i=0;i<trial.samples.size();++i) {
        const auto &q=trial.samples[i];
        Require(std::abs(q.u)<.5 && std::abs(q.v)<.3 && std::abs(q.w)<.4 &&
           std::hypot(q.u-std::round(q.u/.3)*.3,q.v)>NA60Plane::wireRadius,
           "trial coordinate outside three-wire gas");
        const double phi=map.WeightingPotential(q.u,q.v,q.w,label);
        // For these early central-wire test points, an exact zero is suspicious
        // (possible out-of-map return). Fail rather than silently accept.
        Require(std::isfinite(phi)&&phi>-.05&&phi<1.05&&phi!=0.,
                "zero/nonfinite/out-of-range trial weighting potential");
        if(i==0) initial=prev=phi;
        csv<<j<<","<<trial.radiusUm<<","<<trial.theta<<","<<trial.wMm<<","
           <<label<<","<<q.t<<","<<q.u<<","<<q.v<<","<<q.w<<","
           <<phi<<","<<phi-initial<<","<<phi-prev<<"\n";
        prev=phi;
      }
    }
  }
  csv.close();
  Require(bool(csv),"failed writing trial CSV");
  std::cout<<"Wrote "<<trials.size()<<" trials, both U/V: "<<out<<"\n";
}
} // namespace
int main(int argc,char** argv) {
 try {
  Require(argc==5,"usage: PROGRAM MAP_V5 MAP_V6 MOBILITY NEW_OUTPUT_DIR");
  const fs::path out=argv[4];
  Require(!fs::exists(out),"output directory exists");
  Require(fs::create_directories(out),"cannot create new output directory");
  PhysicalIonPaths transport(argv[3],.005);
  constexpr double pi=3.14159265358979323846;
  std::vector<Trial> trials;
  for(double r:{15.1,16.5,19.,23.}) for(int a=0;a<8;++a)
    for(double w:{0.,.85}) {
      const double theta=(a+.5)*pi/4.;
      const double rad=r*1.e-4; // um -> cm
      NA60Plane::Archive::Birth b{rad*std::cos(theta),rad*std::sin(theta),
                                   .1*w,0.,1.,0};
      const auto path=transport.Run(b);
      Require(path.end.t>100.,"trajectory ended before requested early window");
      Trial tr{r,theta,w,{}};
      for(std::size_t i=0;i<tr.samples.size();++i)
        tr.samples[i]=Sample(path.points,10.*i);
      trials.push_back(tr);
    }
  Require(trials.size()==64,"wrong canonical trajectory count");
  std::cout<<"Generated 64 fixed central-wire analytic paths (0..100ns sampling)\n";
  // Map objects are destroyed between evaluations to limit RAM.
  Evaluate(argv[1],out/"v5.csv",trials);
  Evaluate(argv[2],out/"v6.csv",trials);
  std::cout<<"PASS FIXED_PATH_SAMPLING; comparison is NOT mesh convergence certification\n";
  return 0;
 } catch(const std::exception& e) {
  std::cerr<<"FIXED_PATH_FAILURE: "<<e.what()<<"\n"; return 2;
 }
}
