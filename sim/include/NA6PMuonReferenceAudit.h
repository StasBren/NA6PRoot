// NA6PCCopyright
//
// Optional non-invasive virtual-Z-plane scoring of Geant4 muon tracks.
// Not a detector, not a stepping limiter, and does not alter transport.
// Samples are linearly INTERPOLATED within existing Geant4 steps, so they
// do not claim exact Geant4 transport intersection states inside a long step.
// Stations can be re-projected to new MODULE LAYOUTS if passive material and
// field geometry stay fixed (otherwise this is only an approximation).
//
// Original stack track ID is retained, BEFORE MCKine selection/remapping.
#ifndef NA6P_MUON_REFERENCE_AUDIT_H
#define NA6P_MUON_REFERENCE_AUDIT_H

#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TVirtualMC.h>
#include <TSystem.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

class NA6PMuonReferenceAudit final
{
 public:
  explicit NA6PMuonReferenceAudit(const std::string& directory)
  {
    // Full working z range, including potential future station relocations.
    // 10-cm spacing is a compromise; points 365/430/495 resolve magnet
    // boundary/reference planes from the historical physics analysis.
    for (int z = 280; z <= 900; z += 10) {
      mZ.push_back(static_cast<double>(z));
    }
    for (double z : {365., 430., 495.}) {
      mZ.push_back(z);
    }
    std::sort(mZ.begin(), mZ.end());
    mZ.erase(std::unique(mZ.begin(), mZ.end()), mZ.end());

    const auto filename = (std::filesystem::path(directory) /
                           "MuonReferencePlanes.root").string();
    mFile.reset(TFile::Open(filename.c_str(), "RECREATE"));
    if (!mFile || mFile->IsZombie()) {
      throw std::runtime_error(
          "Could not create geometry-independent muon reference file: " +
          filename);
    }
    mFile->cd();
    mTree = new TTree(
        "Planes",
        "Geant4 transported muon state interpolated at virtual Z slices; "
        "unremapped trackID; cm, GeV/c");
    mTree->Branch("event", &mRow.event, "event/I");
    mTree->Branch("trackID", &mRow.trackID, "trackID/I");
    mTree->Branch("pdg", &mRow.pdg, "pdg/I");
    mTree->Branch("planeID", &mRow.planeID, "planeID/I");
    mTree->Branch("crossingIndex", &mRow.crossingIndex, "crossingIndex/I");
    mTree->Branch("forward", &mRow.forward, "forward/I");
    mTree->Branch("z_cm", &mRow.z_cm, "z_cm/D");
    mTree->Branch("x_cm", &mRow.x_cm, "x_cm/D");
    mTree->Branch("y_cm", &mRow.y_cm, "y_cm/D");
    mTree->Branch("px_GeV", &mRow.px_GeV, "px_GeV/D");
    mTree->Branch("py_GeV", &mRow.py_GeV, "py_GeV/D");
    mTree->Branch("pz_GeV", &mRow.pz_GeV, "pz_GeV/D");
    mTree->Branch("z_step_start_cm", &mRow.z_step_start_cm, "z_step_start_cm/D");
    mTree->Branch("z_step_end_cm", &mRow.z_step_end_cm, "z_step_end_cm/D");
    mTree->Branch("interpolation_fraction", &mRow.interpolation_fraction,
                  "interpolation_fraction/D");
  }

  NA6PMuonReferenceAudit(const NA6PMuonReferenceAudit&) = delete;
  NA6PMuonReferenceAudit& operator=(const NA6PMuonReferenceAudit&) = delete;

  ~NA6PMuonReferenceAudit()
  {
    if (mFile && !mFile->IsZombie()) {
      mFile->cd();
      if (mTree) {
        mTree->Write();
      }
      mFile->Close();
    }
  }

  void beginEvent()
  {
    mHasPrevious = false;
    mCurrentTrack = -1;
    mCrossingIndex = 0;
  }

  void preTrack(int event, int trackID, int pdg)
  {
    mCurrentTrack = trackID;
    mCrossingIndex = 0;
    mHasPrevious = std::abs(pdg) == 13;
    if (!mHasPrevious) {
      return;
    }
    auto* mc = TVirtualMC::GetMC();
    TLorentzVector pos, mom;
    mc->TrackPosition(pos);
    mc->TrackMomentum(mom);
    mPrevious = {pos.X(), pos.Y(), pos.Z(), mom.Px(), mom.Py(), mom.Pz()};
    mEvent = event;
    mPDG = pdg;
  }

  void step(int event, int trackID, int pdg)
  {
    if (std::abs(pdg) != 13) {
      return;
    }
    auto* mc = TVirtualMC::GetMC();
    TLorentzVector pos, mom;
    mc->TrackPosition(pos);
    mc->TrackMomentum(mom);
    const State current =
        {pos.X(), pos.Y(), pos.Z(), mom.Px(), mom.Py(), mom.Pz()};
    if (!mHasPrevious || mCurrentTrack != trackID ||
        mEvent != event || mPDG != pdg) {
      // The normal VMC lifecycle calls PreTrack before Stepping. If it did
      // not, establish an initial state but NEVER invent a first segment.
      mPrevious = current;
      mCurrentTrack = trackID;
      mEvent = event;
      mPDG = pdg;
      mHasPrevious = true;
      mCrossingIndex = 0;
      return;
    }

    const double dz = current.z - mPrevious.z;
    if (std::abs(dz) > 1.e-12) {
      const double low = std::min(current.z, mPrevious.z);
      const double high = std::max(current.z, mPrevious.z);
      // Do not duplicate the same plane at two touching successive steps:
      // low is excluded; upper crossing is included. Reverse trajectories
      // are kept and labelled separately.
      for (size_t i = 0; i < mZ.size(); ++i) {
        const double z = mZ[i];
        if (!(low < z && z <= high)) {
          continue;
        }
        const double fraction = (z - mPrevious.z) / dz;
        auto lerp = [fraction](double a, double b) {
          return a + fraction * (b - a);
        };
        mRow.event = event;
        mRow.trackID = trackID;
        mRow.pdg = pdg;
        mRow.planeID = static_cast<int>(i);
        mRow.crossingIndex = mCrossingIndex++;
        mRow.forward = dz > 0. ? 1 : 0;
        mRow.z_cm = z;
        mRow.x_cm = lerp(mPrevious.x, current.x);
        mRow.y_cm = lerp(mPrevious.y, current.y);
        mRow.px_GeV = lerp(mPrevious.px, current.px);
        mRow.py_GeV = lerp(mPrevious.py, current.py);
        mRow.pz_GeV = lerp(mPrevious.pz, current.pz);
        mRow.z_step_start_cm = mPrevious.z;
        mRow.z_step_end_cm = current.z;
        mRow.interpolation_fraction = fraction;
        mTree->Fill();
      }
    }
    mPrevious = current;
  }

 private:
  struct State {
    double x = 0., y = 0., z = 0.;
    double px = 0., py = 0., pz = 0.;
  };
  struct Row {
    int event = 0, trackID = -1, pdg = 0;
    int planeID = 0, crossingIndex = 0, forward = 0;
    double z_cm = 0., x_cm = 0., y_cm = 0.;
    double px_GeV = 0., py_GeV = 0., pz_GeV = 0.;
    double z_step_start_cm = 0., z_step_end_cm = 0.;
    double interpolation_fraction = 0.;
  };

  std::unique_ptr<TFile> mFile;
  TTree* mTree = nullptr;
  std::vector<double> mZ;
  Row mRow;
  State mPrevious;
  bool mHasPrevious = false;
  int mCurrentTrack = -1, mEvent = -1, mPDG = 0, mCrossingIndex = 0;
};

#endif // NA6P_MUON_REFERENCE_AUDIT_H
