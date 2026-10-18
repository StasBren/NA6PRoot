#include "NA6PMuonSpecHit.h"

#include <TFile.h>
#include <TGeoManager.h>
#include <TGeoMatrix.h>
#include <TTree.h>
#include <TParticle.h>

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace {

struct SensorFrame {
  int station = -1;
  TGeoHMatrix localToMaster;
};

bool buildSensorFrame(int detID, SensorFrame& out)
{
  if (!gGeoManager) return false;

  for (int station = 0; station < 16; ++station) {
    std::ostringstream name;
    name << "MS_Lr" << station << "_Ch" << detID << "_Gas";
    auto* entry = gGeoManager->GetAlignableEntry(name.str().c_str());
    if (!entry) continue;

    const char* path = entry->GetPath();
    if (!path || !path[0]) continue;
    if (!gGeoManager->cd(path)) continue;

    auto* matrix = gGeoManager->GetCurrentMatrix();
    if (!matrix) continue;

    out.station = station;
    out.localToMaster = *matrix;
    return true;
  }
  return false;
}

void masterPointToGarfield(const TGeoHMatrix& m,
                           double x, double y, double z,
                           double& u_mm, double& v_mm, double& w_mm)
{
  const double master[3] = {x, y, z};
  double local[3] = {0., 0., 0.};
  m.MasterToLocal(master, local);

  // Installed +90 degree chamber convention used throughout Stage C:
  //   chamber-local x -> detector Y -> Garfield u (across wires)
  //   chamber-local z -> detector Z -> Garfield v (chamber normal)
  //   chamber-local y -> -detector X -> -Garfield w
  u_mm = 10. * local[0];
  v_mm = 10. * local[2];
  w_mm = -10. * local[1];
}

void masterVectorToGarfield(const TGeoHMatrix& m,
                            double px, double py, double pz,
                            double& pu, double& pv, double& pw)
{
  const double master[3] = {px, py, pz};
  double local[3] = {0., 0., 0.};
  m.MasterToLocalVect(master, local);
  pu = local[0];
  pv = local[2];
  pw = -local[1];
}

} // namespace

void export_c_geant_muon_crossings(
    const char* hitsFile = "HitsMuonSpec.root",
    const char* geometryFile = "geometry.root",
    const char* outputCsv = "c_geant_gas_crossings.csv",
    const char* kineFile = "")
{
  TFile geomIn(geometryFile, "READ");
  if (geomIn.IsZombie()) {
    std::cerr << "Could not open geometry file " << geometryFile << "\n";
    return;
  }

  auto* geo = dynamic_cast<TGeoManager*>(geomIn.Get("Geometry"));
  if (!geo) {
    // geometry.root written by TGeoManager::Export is usually imported this way.
    geo = TGeoManager::Import(geometryFile);
  }
  if (!geo) {
    std::cerr << "Could not load TGeo geometry from " << geometryFile << "\n";
    return;
  }
  gGeoManager = geo;

  TFile in(hitsFile, "READ");
  if (in.IsZombie()) {
    std::cerr << "Could not open " << hitsFile << "\n";
    return;
  }

  auto* tree = dynamic_cast<TTree*>(in.Get("hitsMuonSpec"));
  if (!tree) {
    std::cerr << "Tree hitsMuonSpec not found in " << hitsFile << "\n";
    return;
  }

  std::vector<NA6PMuonSpecHit>* hits = nullptr;
  if (tree->SetBranchAddress("MuonSpec", &hits) < 0) {
    std::cerr << "Branch MuonSpec not found/readable in hitsMuonSpec\n";
    return;
  }

  // Hits are REMAPPED to MCKine kept-track indices. Use their matching
  // merged chunk MCKine tree (NOT gasSteps original-stack track IDs).
  // Truth connection is optional for older 3-argument callers.
  std::unique_ptr<TFile> kineIn;
  TTree* ktree = nullptr;
  std::vector<TParticle>* tracks = nullptr;
  if (kineFile && kineFile[0]) {
    kineIn.reset(TFile::Open(kineFile, "READ"));
    if (!kineIn || kineIn->IsZombie()) {
      std::cerr << "Missing input kine file " << kineFile << "\n";
      return;
    }
    ktree = dynamic_cast<TTree*>(kineIn->Get("mckine"));
    if (!ktree || ktree->SetBranchAddress("tracks", &tracks) < 0 ||
        ktree->GetEntries() != tree->GetEntries()) {
      std::cerr << "Cannot align MCKine tracks and MuonSpec events\n";
      return;
    }
  }

  std::ofstream out(outputCsv);
  if (!out) {
    std::cerr << "Could not create " << outputCsv << "\n";
    return;
  }
  out << std::setprecision(12);
  out
      << "event,hit_index,track_id,det_id,station,"
      << "x_in_cm,y_in_cm,z_in_cm,x_out_cm,y_out_cm,z_out_cm,"
      << "u_in_mm,v_in_mm,w_in_mm,u_out_mm,v_out_mm,w_out_mm,"
      << "px_in_GeV,py_in_GeV,pz_in_GeV,"
      << "pu_in_GeV,pv_in_GeV,pw_in_GeV,"
      << "path_length_mm,edep_GeV,edep_eV,time_ns,"
      << "px_out_GeV,py_out_GeV,pz_out_GeV,"
      << "pu_out_GeV,pv_out_GeV,pw_out_GeV,"
      << "theta_u_in_deg,theta_w_in_deg,"
      << "theta_u_out_deg,theta_w_out_deg,"
      << "delta_theta_u_deg,delta_theta_w_deg,"
      << "status_start,status_end,pdg,mother_track_id\n";

  std::map<int, SensorFrame> frames;
  Long64_t rows = 0;
  double sumEdep = 0.;

  const Long64_t nEvents = tree->GetEntries();
  for (Long64_t event = 0; event < nEvents; ++event) {
    tree->GetEntry(event);
    if (ktree) {
      ktree->GetEntry(event);
      if (!tracks) {
        std::cerr << "MCKine missing tracks at event " << event << "\n";
        return;
      }
    }
    if (!hits) continue;

    for (size_t ih = 0; ih < hits->size(); ++ih) {
      const auto& h = hits->at(ih);
      const int detID = h.getDetectorID();
      int pdg = 0, mother = -1;
      if (ktree) {
        const int tid = h.getTrackID();
        if (tid < 0 || static_cast<size_t>(tid) >= tracks->size()) {
          std::cerr << "MuonSpec trackID=" << tid
                    << " not in retained MCKine for event " << event
                    << ", size=" << tracks->size() << "\n";
          return;
        }
        pdg = tracks->at(tid).GetPdgCode();
        mother = tracks->at(tid).GetFirstMother();
      }

      auto it = frames.find(detID);
      if (it == frames.end()) {
        SensorFrame frame;
        if (!buildSensorFrame(detID, frame)) {
          std::cerr
              << "Could not resolve alignable gas frame for detID=" << detID
              << ". Expected symbolic name MS_Lr*_Ch" << detID << "_Gas\n";
          return;
        }
        it = frames.emplace(detID, frame).first;
      }
      const auto& frame = it->second;

      double u0, v0, w0, u1, v1, w1;
      masterPointToGarfield(
          frame.localToMaster,
          h.getXIn(), h.getYIn(), h.getZIn(),
          u0, v0, w0);
      masterPointToGarfield(
          frame.localToMaster,
          h.getXOut(), h.getYOut(), h.getZOut(),
          u1, v1, w1);

      double pu, pv, pw;
      masterVectorToGarfield(
          frame.localToMaster,
          h.getPXIn(), h.getPYIn(), h.getPZIn(),
          pu, pv, pw);

      double pu1, pv1, pw1;
      masterVectorToGarfield(
          frame.localToMaster,
          h.getPXOut(), h.getPYOut(), h.getPZOut(),
          pu1, pv1, pw1);

      constexpr double deg = 180. / 3.14159265358979323846;
      const double thetaUIn = std::atan2(pu, pv) * deg;
      const double thetaWIn = std::atan2(pw, pv) * deg;
      const double thetaUOut = std::atan2(pu1, pv1) * deg;
      const double thetaWOut = std::atan2(pw1, pv1) * deg;

      const double du = u1 - u0;
      const double dv = v1 - v0;
      const double dw = w1 - w0;
      const double length = std::sqrt(du * du + dv * dv + dw * dw);
      const double edep = h.getHitValue();

      out
          << event << "," << ih << ","
          << h.getTrackID() << "," << detID << "," << frame.station << ","
          << h.getXIn() << "," << h.getYIn() << "," << h.getZIn() << ","
          << h.getXOut() << "," << h.getYOut() << "," << h.getZOut() << ","
          << u0 << "," << v0 << "," << w0 << ","
          << u1 << "," << v1 << "," << w1 << ","
          << h.getPXIn() << "," << h.getPYIn() << "," << h.getPZIn() << ","
          << pu << "," << pv << "," << pw << ","
          << length << "," << edep << "," << 1.e9 * edep << ","
          << h.getTime() << ","
      << h.getPXOut() << "," << h.getPYOut() << "," << h.getPZOut() << ","
      << pu1 << "," << pv1 << "," << pw1 << ","
      << thetaUIn << "," << thetaWIn << ","
      << thetaUOut << "," << thetaWOut << ","
      << (thetaUOut - thetaUIn) << "," << (thetaWOut - thetaWIn) << ","
      << static_cast<int>(h.getStatusStart()) << ","
      << static_cast<int>(h.getStatusEnd()) << ","
      << pdg << "," << mother << "\n";

      ++rows;
      sumEdep += edep;
    }
  }

  std::cout
      << "GEANT MWPC crossing export COMPLETE\n"
      << "  events              : " << nEvents << "\n"
      << "  hit/crossing rows   : " << rows << "\n"
      << "  detector frames     : " << frames.size() << "\n"
      << "  total Edep [eV]     : " << 1.e9 * sumEdep << "\n"
      << "  MCKine PDG join     : " << (ktree ? "enabled" : "not requested") << "\n"
      << "  wrote               : " << outputCsv << "\n";
}
