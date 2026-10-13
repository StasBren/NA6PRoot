// NA6PCCopyright

#ifndef NA6P_MUONSPEC_H_
#define NA6P_MUONSPEC_H_

#include "NA6PModule.h"
#include "NA6PMuonSpecHit.h"

class TFile;
class TTree;
class TGeoVolume;

class NA6PMuonSpec : public NA6PModule
{
 public:
  NA6PMuonSpec() : NA6PModule("MuonSpec") { setActiveID(1); }
  ~NA6PMuonSpec() override = default;
  void createMaterials() override;
  void createGeometry(TGeoVolume* world) override;
  bool stepManager(int volID) override;
  size_t getNHits() const override { return mHits.size(); }
  void createHitsOutput(const std::string& outDir) override;
  void closeHitsOutput() override;
  void writeHits(const std::vector<int>& remapping) override;
  void setAlignableEntries() override;

  void clearHits() override { mHits.clear(); }

  const auto& getHits() const { return mHits; }

  NA6PMuonSpecHit* addHit(int trackID, int detID, const TVector3& startPos, const TVector3& endPos, const TVector3& startMom, const TVector3& endMom,
                                 float endTime, float eLoss, unsigned char startStatus, unsigned char endStatus);

 private:
  std::vector<NA6PMuonSpecHit> mHits, *hHitsPtr = &mHits;
  TFile* mHitsFile = nullptr;
  TTree* mHitsTree = nullptr;

  // Optional diagnostic GEANT step dump in the sensitive MWPC gas.
  // Enabled only when NA6P_MWPC_STORE_GAS_STEPS is set to a non-zero value.
  // It stores non-zero per-step energy deposits, not microscopic ionisation
  // clusters.  The latter can be sampled later from the deposited energy.
  bool mStoreGasSteps = false;
  TTree* mGasStepTree = nullptr;
  int mGasStepEvent = 0;
  int mGasStepTrackID = 0;
  int mGasStepDetID = 0;
  double mGasStepX = 0.;
  double mGasStepY = 0.;
  double mGasStepZ = 0.;
  double mGasStepU = 0.;
  double mGasStepV = 0.;
  double mGasStepW = 0.;
  double mGasStepPx = 0.;
  double mGasStepPy = 0.;
  double mGasStepPz = 0.;
  double mGasStepEdep = 0.;
};

#endif
