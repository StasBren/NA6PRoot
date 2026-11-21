// Deterministic polar single-ion RKF bank for the 2100 V symmetric P3 model.
// NO new Geant4/TrackHeed/AvalancheMicroscopic. Whole cathode planes ONLY.
// Explicit grid bounds, finite-window tails, and no overwrite of frozen data.
// Uses exactly Plane2100Response.hh field, mobility and 10 ns sensor settings.
#include "Plane2100Response.hh"

#include <TFile.h>
#include <TTree.h>
#include <TNamed.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
namespace fs = std::filesystem;
std::string Arg(int n, char** argv, const std::string& key) {
  for (int i = 1; i+1 < n; ++i) if (key == argv[i]) return argv[i+1];
  return {};
}
int IntArg(int n, char** argv, const std::string& key, int fallback) {
  const auto s = Arg(n, argv, key);
  return s.empty() ? fallback : std::stoi(s);
}
double DoubleArg(int n, char** argv, const std::string& key, double fallback) {
  const auto s = Arg(n, argv, key);
  return s.empty() ? fallback : std::stod(s);
}
constexpr double PI=3.14159265358979323846264338327950288;
constexpr double CM_TO_UM=1.e4;
}  // namespace

int main(int argc, char** argv) {
  try {
    const auto target = Arg(argc, argv, "--out-root");
    const auto mobility = Arg(argc, argv, "--mobility");
    const int nr = IntArg(argc, argv, "--n-radius", 6);
    const int nt = IntArg(argc, argv, "--n-angle", 12);
    const int nPhase = IntArg(argc, argv, "--phase-subdivisions", 1);
    const double clearanceUm = DoubleArg(argc, argv, "--clearance-min-um", 0.00001);
    const double radiusMaxUm = DoubleArg(argc, argv, "--radius-max-um", 80.);
    const double thetaMinDeg = DoubleArg(argc, argv, "--theta-min-deg", 0.001);
    const double dt = DoubleArg(argc, argv, "--dt-ns", 10.);
    const double windowUs = DoubleArg(argc, argv, "--window-us", 100.);
    const double stepMm = DoubleArg(argc, argv, "--max-step-mm", 0.05);
    const int wireIndex = IntArg(argc, argv, "--wire-index", 0);
    NA60Plane::Require(!target.empty() && !mobility.empty(),
                       "Usage: --out-root FILE --mobility ION_MOBILITY "
                       "[--n-radius 6 --n-angle 12]");
    NA60Plane::Require(nPhase == 1 || nPhase == 2 ||
                       nPhase == 4 || nPhase == 8,
                       "phase subdivisions must be 1,2,4,8");
    NA60Plane::Require(nr >= 3 && nr <= 64 && nt >= 6 && nt <= 128 &&
                       nr * nt * nPhase <= 2048,
                       "bounded grid supports at most 2048 RKF trajectory starts");
    NA60Plane::Require(wireIndex == 0,
        "first study is explicitly limited to the physical central wire");
    NA60Plane::Require(std::isfinite(clearanceUm) && clearanceUm > 0. &&
                       clearanceUm < .1 &&
                       std::isfinite(radiusMaxUm) &&
                       radiusMaxUm > 30. && radiusMaxUm <= 100.,
                       "invalid radial bounds");
    NA60Plane::Require(std::isfinite(thetaMinDeg) && thetaMinDeg > 0. &&
                       thetaMinDeg < 1.,
                       "must stay off exact v=0 separatrix");
    NA60Plane::Require(dt == 10. && windowUs == 100. &&
                       stepMm > 0. && stepMm <= .1,
                       "freeze dt=10ns, window=100us, RKF max step<=0.1mm");
    NA60Plane::Require(!fs::exists(target),
                       "refusing to overwrite any existing output ROOT: "+target);
    const fs::path destination = target;
    NA60Plane::Require(fs::exists(destination.parent_path()),
                       "parent output directory must already exist");

    const int nbins = 10000;
    const double rMinUm = NA60Plane::wireRadius * CM_TO_UM + clearanceUm;
    const double minLog = std::log(rMinUm);
    const double logSpan = std::log(radiusMaxUm / rMinUm);
    const double tSpan = 180. - 2. * thetaMinDeg;
    NA60Plane::Transport drift(mobility, dt, nbins, stepMm, true);
    TFile root(target.c_str(), "CREATE", "Deterministic polar 2100V full-plane RKF", 4);
    NA60Plane::Require(!root.IsZombie(), "cannot create ROOT bank");

    int iRadius=0, iAngle=0, iPhase=0, status=0;
    double rUm=0., thetaDeg=0., u=0., v=0., w=0., endT=0., phaseNs=0.;
    double recordedPlus=0., recordedMinus=0., unrecordedPlus=0., unrecordedMinus=0.;
    bool late=false;
    std::vector<double> waveform;
    TTree bank("polar_bank", "Central-wire deterministic positive-v ion impulses");
    bank.Branch("i_radius", &iRadius);
    bank.Branch("i_angle", &iAngle);
    bank.Branch("i_phase", &iPhase);
    bank.Branch("phase_ns", &phaseNs);
    bank.Branch("r_um", &rUm);
    bank.Branch("theta_deg", &thetaDeg);
    bank.Branch("u_cm", &u);
    bank.Branch("v_cm", &v);
    bank.Branch("w_cm", &w);
    bank.Branch("end_time_ns", &endT);
    bank.Branch("endpoint_status", &status);
    bank.Branch("late_after_100us", &late);
    bank.Branch("recorded_plus_fC", &recordedPlus);
    bank.Branch("recorded_minus_fC", &recordedMinus);
    bank.Branch("unrecorded_plus_fC", &unrecordedPlus);
    bank.Branch("unrecorded_minus_fC", &unrecordedMinus);
    bank.Branch("waveform_plus_then_minus", &waveform);

    TNamed scope("scope",
      "Deterministic 2100V nine-wire P3 noFW; one central wire, positive-v "
      "single ion birth, dt=10ns, 100us window, constant w=0; TRUE Garfield++ "
      "DriftLineRKF (allow late cathode arrivals); provisional Ar+ in Ar ion "
      "mobility for Ar/CO2. Two ideal whole cathode planes, NO actual stereo "
      "strips or electronics. Existing muon ROOTs are READ ONLY. No new "
      "microscopic avalanche. Grid interior excludes v=0 exact separatrix.");
    scope.Write();
    TNamed radialRule("radial_nodes",
       "Uniform in natural log of r_um, from 15um + clearance to radius_max_um");
    radialRule.Write();
    TNamed angularRule("angular_nodes",
       "Theta = theta_min + (180-2 theta_min)*(1-cos(pi*j/(n_angle-1)))/2; "
       "includes endpoints near (not at) v=0; theta is atan2(v,u_local) in degrees.");
    angularRule.Write();
    TNamed phaseRule("phase_nodes",
      "One phase per ion starting point for nPhase=1; otherwise start-time "
      "phases p * 10ns/nPhase for p in 0..nPhase-1; 10ns endpoint is phase0 "
      "template shifted one integer bin. Positive-v central wire only.");
    phaseRule.Write();
    TNamed mobilityName("mobility_file", mobility.c_str());
    mobilityName.Write();

    const auto start=std::chrono::steady_clock::now();
    int lateCount=0;
    double maxEnd=0.;
    for (iRadius=0; iRadius<nr; ++iRadius) {
      rUm = std::exp(minLog + (logSpan * iRadius) / (nr - 1));
      for (iAngle=0; iAngle<nt; ++iAngle) {
        const double weight = (1. - std::cos(PI*iAngle/(nt-1))) / 2.;
        thetaDeg = thetaMinDeg + tSpan * weight;
        const double theta=thetaDeg*PI/180.;
        u = rUm * std::cos(theta) / CM_TO_UM;
        v = rUm * std::sin(theta) / CM_TO_UM;
        w = 0.;
        for (iPhase=0; iPhase<nPhase; ++iPhase) {
          phaseNs = dt * iPhase / nPhase;
          drift.CheckStart(u, v, w, phaseNs);
          const auto d = drift.Drift(u, v, w, phaseNs);
          endT=d.endT;
          status=d.status;
          late=d.late;
          recordedPlus=d.integralQ[0]; recordedMinus=d.integralQ[1];
          unrecordedPlus=d.unrecordedTailQ[0]; unrecordedMinus=d.unrecordedTailQ[1];
          maxEnd=std::max(maxEnd,endT);
          if(late) ++lateCount;
          waveform.assign(2*nbins, 0.);
          for(int bin=0;bin<nbins;++bin) {
            waveform[bin] = drift.Current(0,bin);
            waveform[nbins+bin] = drift.Current(1,bin);
            NA60Plane::Require(std::isfinite(waveform[bin]) &&
                              std::isfinite(waveform[nbins+bin]),
                              "nonfinite deterministic single-ion current");
          }
          bank.Fill();
          std::cout<<std::setprecision(12)<<"GRID "
            <<((iRadius*nt+iAngle)*nPhase+iPhase+1)<<"/"<<(nr*nt*nPhase)
            <<" r_um="<<rUm<<" theta_deg="<<thetaDeg
            <<" phase_ns="<<phaseNs
            <<" end_time_ns="<<endT
            <<" late="<<late<<"\n"<<std::flush;
        }
      }
    }
    root.cd(); bank.Write(); root.Close();
    const double sec=std::chrono::duration<double>(
        std::chrono::steady_clock::now()-start).count();
    std::cout<<"PASS deterministic central-wire polar RKF grid"
             <<" rows="<<nr*nt*nPhase<<" radius_nodes="<<nr<<" angle_nodes="<<nt
             <<" phase_nodes="<<nPhase
             <<" late_rows="<<lateCount<<" max_end_ns="<<maxEnd
             <<" elapsed_seconds="<<sec
             <<"\nOutput="<<target
             <<"\nNOTE: successful transport is not a verified interpolation error. "
             <<"No exact v=0 boundary or full-muon reference generated.\n";
    return 0;
  } catch(const std::exception& error) {
    std::cerr<<"POLAR RKF BANK FAILED: "<<error.what()<<"\n";
    return 2;
  }
}
