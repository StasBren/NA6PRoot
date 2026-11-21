// NA60+ P3/2100 V: genuine frozen K256 ion births -> provisional U/V
// central-strip signal. No avalanche regeneration, no electronics.
// Uses the original analytic nine-wire physical drift and numerically computed
// Elmer single-strip weighting POTENTIALS. Currents are charge per 10ns bin.
//
// IMPORTANT: archived baseline FEM wire surfaces are faceted; until a curved
// mesh and convergence are checked, every calculated signal is PROVISIONAL.
// Noncentral translated channels are opt-in and NOT validated by this program.
#include "Plane2100Archive.hh"
#include "TimeBinnedRamoCharge.hh"
#include "Garfield/ComponentElmer.hh"
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
namespace fs = std::filesystem;
using NA60Plane::Require;
struct Point { double u, v, w, t; };  // cm, cm, cm, ns

std::string Arg(int argc, char** argv, const std::string& key) {
  for (int i=1; i+1<argc; ++i) if (key==argv[i]) return argv[i+1];
  return "";
}
bool Flag(int argc, char** argv, const std::string& key) {
  for (int i=1; i<argc; ++i) if (key==argv[i]) return true;
  return false;
}
double Number(int argc, char** argv, const std::string& key, double def) {
  const auto s=Arg(argc,argv,key);
  return s.empty() ? def : std::stod(s);
}
void NeedFile(const fs::path& p) {
  Require(fs::is_regular_file(p) && fs::file_size(p)>0,
          "missing map file "+p.string());
}

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

struct Channel {
  std::string name, label;
  int k=0;  // reference physical strip index (zero phase assumed)
  int relativeK=0;
  bool plus=false;
  std::unique_ptr<NA60Stereo::TimeBinnedRamoCharge> current;
  double analyticEndpointFc=0.;
};

class StereoMap {
  Garfield::ComponentElmer map_;
  static constexpr double kPitchCm = 0.17;  // physical 1.7mm along wires
  static constexpr double kSlope = 0.1;
  static constexpr double kHalfWidthCm = 0.085;
 public:
  explicit StereoMap(const fs::path& dir)
      : map_((dir/"chamber/mesh.header").string(),
             (dir/"chamber/mesh.elements").string(),
             (dir/"chamber/mesh.nodes").string(),
             (dir/"dielectrics.dat").string(),
             (dir/"chamber/drift.result").string(),"mm") {
    for (const auto* f: {"chamber/mesh.header","chamber/mesh.elements",
         "chamber/mesh.nodes","dielectrics.dat","chamber/drift.result",
         "chamber/weight_plus_0.result","chamber/weight_minus_0.result"})
      NeedFile(dir/f);
    Require(map_.SetWeightingField(
      (dir/"chamber/weight_plus_0.result").string(),"plus_strip_0"),
      "cannot import plus FEM weighting map");
    Require(map_.SetWeightingField(
      (dir/"chamber/weight_minus_0.result").string(),"minus_strip_0"),
      "cannot import minus FEM weighting map");
    const double nearPlus=map_.WeightingPotential(0.,.295,0.,"plus_strip_0");
    const double farPlus=map_.WeightingPotential(0.,-.295,0.,"plus_strip_0");
    const double nearMinus=map_.WeightingPotential(0.,-.295,0.,"minus_strip_0");
    const double farMinus=map_.WeightingPotential(0.,.295,0.,"minus_strip_0");
    Require(nearPlus>.5 && nearMinus>.5 &&
            std::abs(farPlus)<.1 && std::abs(farMinus)<.1,
            "FEM channel sign/weighting boundary sanity failed");
  }
  double Phi(const Point& p, const Channel& c, const ChamberFrame& frame,
             bool endpoint=false) {
    const double slope=c.plus?kSlope:-kSlope;
    // Convert to the station-local strip phase while querying the same FEM
    // central-strip model in patch-local coordinates.
    const double w=p.w+frame.wEntranceCm-slope*frame.patchShiftCm-
                   kPitchCm*c.k;
    Require(std::abs(p.u)<1.55-1.e-4 && std::abs(w)<1.1-1.e-4 &&
            std::abs(p.v)<=.300001,
            "FEM map coordinate out of patch; refusing silent zero lookup");
    if (endpoint) {
      // At a physical cathode, the exact ideal contiguous-strip Dirichlet
      // boundary is either 0 or 1. Never query outside-gas endpoints as 0.
      const bool onPlane=(p.v>0.)==c.plus;
      if (!onPlane) return 0.;
      const double phase=p.w+frame.wEntranceCm-
                         slope*(p.u+frame.patchShiftCm)-kPitchCm*c.k;
      Require(std::abs(std::abs(phase)-kHalfWidthCm)>1.e-7,
              "ion landed at ambiguous ideal strip boundary");
      return std::abs(phase)<kHalfWidthCm ? 1. : 0.;
    }
    const double rho=std::hypot(p.u-std::round(p.u/.3)*.3,p.v);
    Require(rho>NA60Plane::wireRadius,
            "weighting query starts within wire metal");
    const double phi=map_.WeightingPotential(p.u,p.v,w,c.label);
    Require(std::isfinite(phi) && phi>-.05 && phi<1.05,
            "nonfinite/out-of-range FEM strip potential");
    return phi;
  }
};

} // namespace

int main(int argc,char** argv) {
 try {
  const auto ionFile=Arg(argc,argv,"--ion-file");
  const auto mapDir=Arg(argc,argv,"--map-dir");
  const auto mobility=Arg(argc,argv,"--mobility");
  const auto prefix=Arg(argc,argv,"--out-prefix");
  const double dt=Number(argc,argv,"--dt-ns",10.);
  const double windowUs=Number(argc,argv,"--window-us",100.);
  const double stepMm=Number(argc,argv,"--step-mm",.05);
  const int halfStrips=int(Number(argc,argv,"--half-strips",0.));
  const long long limit=(long long)Number(argc,argv,"--ion-limit",0.);
  const bool allowTranslation=Flag(argc,argv,"--allow-translated-strips");
  Require(!ionFile.empty()&&!mapDir.empty()&&!mobility.empty()&&!prefix.empty(),
    "usage: --ion-file ION_ROOT --map-dir EXISTING_FEM_MAP "
    "--mobility FILE --out-prefix NEW_PREFIX [--ion-limit 1] "
    "[--half-strips K --allow-translated-strips]");
  Require(dt==10. && windowUs==100., 
    "first validation stage is fixed at 10ns x 100us");
  Require(halfStrips>=0&&halfStrips<=5&&limit>=0&&limit<=20000,
    "invalid half-strips or ion limit");
  Require(halfStrips==0||allowTranslation,
    "noncentral strips require explicit --allow-translated-strips; not FEM validated");
  for(const char* suffix:{"_summary.csv","_channels.csv","_waveforms.csv",
                          "_pilot.txt","_complete.txt"}) {
    Require(!fs::exists(prefix+suffix),"output exists; choose fresh prefix");
  }
  NA60Plane::Archive archive(ionFile);
  ChamberFrame frame(ionFile);
  Require(archive.mode=="weighted" && archive.Entries()<=20000,
          "requires a complete frozen weighted K256 archive <=20k rows");
  const long long nToRun=(limit==0?archive.Entries():
                           std::min(limit,archive.Entries()));
  Require(nToRun>0,"no ion births found");
  const bool complete=(nToRun==archive.Entries());
  const int nBins=int(windowUs*1000./dt);
  StereoMap map(mapDir);
  PhysicalIonPaths transport(mobility,stepMm);
  std::vector<Channel> channels;
  for (int k=-halfStrips;k<=halfStrips;++k)
    for (int side=0;side<2;++side) {
      Channel c;
      c.plus=(side==0);
      c.relativeK=k;
      c.k=(c.plus?frame.uCentralIndex:frame.vCentralIndex)+k;
      c.label=c.plus?"plus_strip_0":"minus_strip_0";
      c.name=(c.plus?"U_refindex_":"V_refindex_")+std::to_string(c.k);
      c.current=std::make_unique<NA60Stereo::TimeBinnedRamoCharge>(dt,nBins);
      channels.push_back(std::move(c));
    }
  long long nLate=0;
  double weightedLate=0.,maxEndNs=0.;
  for(long long j=0;j<nToRun;++j) {
    const auto b=archive.Get(j);
    try {
      frame.CheckBirth(j,b);
      auto track=transport.Run(b);
      maxEndNs=std::max(maxEndNs,track.end.t);
      if(track.end.t>=dt*nBins) {++nLate;weightedLate+=b.weight;}
      for(auto& ch:channels) {
        double previous=map.Phi(track.points.front(),ch,frame);
        const double first=previous;
        for(std::size_t p=1;p<track.points.size();++p) {
          const bool last=p+1==track.points.size();
          const double next=map.Phi(track.points[p],ch,frame,last);
          if(track.points[p].t>track.points[p-1].t)
            ch.current->Add(track.points[p-1].t,track.points[p].t,
                 b.weight*NA60Plane::eFc*(next-previous));
          previous=next;
        }
        ch.analyticEndpointFc+=b.weight*NA60Plane::eFc*(previous-first);
      }
    } catch (const std::exception& err) {
      std::ostringstream msg;
      msg<<"ION ROW FAILED: crossing="<<archive.id
         <<" index="<<j<<" primary="<<b.primary
         <<" start_cm=("<<b.u<<","<<b.v<<","<<b.w<<")"
         <<" birth_ns="<<b.t<<" detail={"<<err.what()<<"}";
      throw std::runtime_error(msg.str());
    }
    if((j+1)%500==0) std::cout<<"PROGRESS ion_birth_rows="<<(j+1)
         <<"/"<<nToRun<<"\n";
  }
  if(complete) archive.CheckComplete();
  for(const auto& c:channels) {
    const double d=std::abs(c.analyticEndpointFc-c.current->FullEndpointFc());
    Require(d <= 1.e-8*std::max(1.,std::abs(c.analyticEndpointFc)),
            "full endpoint Ramo charge failed per channel");
    Require(c.current->CheckClosure(1.e-7),
            "10ns bin closure failed per channel");
  }
  // Write only after all computations and validations completed.
  std::ofstream chout(prefix+"_channels.csv");
  Require(bool(chout),"cannot create channel CSV");
  chout<<std::setprecision(16)
    <<"channel,strip_index,window_charge_fC,full_endpoint_charge_fC,"
      "unrecorded_tail_fC,first100ns_fC\n";
  std::vector<std::vector<double>> wave;
  for(const auto& c:channels) {
    const auto bins=c.current->BinChargeFc();
    double early=0.;
    for(int i=0;i<10;++i) early+=bins[i];
    chout<<c.name<<","<<c.k<<","<<c.relativeK<<","
         <<c.current->WindowChargeFc()<<","
         <<c.current->FullEndpointFc()<<","
         <<c.current->TailFc()<<","<<early<<"\n";
    wave.push_back(bins);
  }
  chout.close();
  std::ofstream csv(prefix+"_waveforms.csv");
  Require(bool(csv),"cannot create waveform CSV");
  csv<<std::setprecision(16)<<"time_center_ns";
  for(const auto& c:channels) csv<<","<<c.name<<"_ion_fC_per_ns";
  csv<<"\n";
  for(int k=0;k<nBins;++k) {
    csv<<(k+.5)*dt;
    for(const auto& b:wave) csv<<","<<b[k]/dt;
    csv<<"\n";
  }
  csv.close();
  std::ofstream summary(prefix+"_summary.csv");
  Require(bool(summary),"cannot create summary CSV");
  summary<<"crossing_id,source,station,det_id,weighting_map,"
           "archive_ion_rows,processed_ion_rows,exact_ion_births,"
           "dt_ns,window_us,step_mm,half_strips,"
           "late_ion_rows,late_weighted_ions,max_end_ns,"
           "archive_complete,field_accuracy_validated,electron_current_included\n";
  summary<<archive.id<<","<<archive.source<<","<<archive.station<<","
         <<archive.det<<","<<mapDir<<","<<archive.Entries()<<","<<nToRun
         <<","<<archive.nTrue<<","<<dt<<","<<windowUs<<","<<stepMm
         <<","<<halfStrips<<","<<nLate<<","<<weightedLate<<","
         <<maxEndNs<<","<<(complete?1:0)<<",0,0\n";
  summary.close();
  std::ofstream seal(prefix+(complete?"_complete.txt":"_pilot.txt"));
  Require(bool(seal),"cannot write completion note");
  seal<<(complete?"COMPLETE_WEIGHTED_ION_CENSUS":"INCOMPLETE_SINGLE_ION_PILOT")<<"\n";
  seal<<"PROVISIONAL_FEM_FIELD=1\n"
      <<"ELECTRON_CURRENTS_INCLUDED=0\n"
      <<"NONCENTRAL_TRANSLATED_STRIPS_NOT_VALIDATED="<<(halfStrips>0?1:0)<<"\n"
      <<"GARFIELD_CONVENTION=positive_q_delta_phi\n";
  seal.close();
  std::cout<<std::setprecision(10)
           <<"PASS_RKFPATH_Ramo "<<(complete?"FULL_ARCHIVE":"EXPLICIT_PILOT")
           <<" ion_rows="<<nToRun<<" channels="<<channels.size()
           <<" late_rows="<<nLate<<" max_end_ns="<<maxEndNs<<"\n"
           <<"PROVISIONAL FEM field. Ion only. No electron currents/electronics.\n"
           <<"CSV prefix: "<<prefix<<"\n";
  return 0;
 } catch (const std::exception& err) {
   std::cerr<<"STEREO_ION_RESPONSE_FAILED: "<<err.what()<<"\n";
   return 2;
 }
}
