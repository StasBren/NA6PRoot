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
};

#endif
