// NA6PCCopyright

#include "NA6PMuonSpec.h"
#include "NA6PMWPCChamber.h"
#include "NA6PMWPCParam.h"
#include "NA6PDetector.h"
#include "NA6PLayoutParam.h"
#include "NA6PMCStack.h"

#include <TVirtualMC.h>
#include <TGeoManager.h>
#include <TGeoMatrix.h>
#include <TGeoNode.h>
#include <TGeoVolume.h>
#include <fairlogger/Logger.h>
#include <TFile.h>
#include <TTree.h>

#include <cmath>
#include <stdexcept>
#include <string>

namespace
{
constexpr const char* LayerName[4] = {"A", "B", "C", "D"};
}

void NA6PMuonSpec::createMaterials()
{
  const auto& p = NA6PMWPCParam::Instance();
  const NA6PMWPCChamber::Materials materials = {
    addName("FR4"),
    addName("Copper"),
    addName("Honeycomb"),
    addName("ArCO2")};
  NA6PMWPCChamber chamber(*this, materials);
  chamber.createMaterials();
}

void NA6PMuonSpec::createGeometry(TGeoVolume* world)
{
  if (!world) {
    throw std::runtime_error("MuonSpec requires a world volume");
  }

  const auto& layout = NA6PLayoutParam::Instance();
  const auto& p = NA6PMWPCParam::Instance();
  if (layout.nMSPlanes < 0 || layout.nMSPlanes > NA6PMWPCParam::MaxStations) {
    throw std::runtime_error("nMSPlanes exceeds the MWPC station-layout capacity");
  }

  // The detailed chamber model keeps the prototype-local mechanical axes
  // (short side = local x, long side = local y).  Installation in the detector
  // rotates the chamber around z so that detector X is horizontal/along the
  // long side (wires) and detector Y is vertical/across the wires.
  const double quarterTurnsReal = static_cast<double>(NA6PMWPCChamber::ChamberRotationZDeg) / 90.;
  const long quarterTurns = std::lround(quarterTurnsReal);
  if (!std::isfinite(NA6PMWPCChamber::ChamberRotationZDeg) ||
      std::abs(quarterTurnsReal - static_cast<double>(quarterTurns)) > 1.e-6) {
    throw std::runtime_error("MWPC chamberRotationZDeg must be a multiple of 90 degrees");
  }
  const bool swapXY = (std::abs(quarterTurns) % 2) == 1;

  const NA6PMWPCChamber::Materials materials = {
    addName("FR4"),
    addName("Copper"),
    addName("Honeycomb"),
    addName("ArCO2")};

  const NA6PMWPCChamber standardChamber(*this, materials);
  // detector Y is chamber-local x after the default +/-90 degree installation
  // rotation, so the requested MS0 vertical sensitive height overrides local x.
  const double ms0BodyShort = static_cast<double>(p.ms0GasY) + 2. * NA6PMWPCChamber::InnerFrameWidth;
  const NA6PMWPCChamber narrowMS0Chamber(*this, materials, ms0BodyShort, p.bodyY);

  standardChamber.createMaterials();

  int requestedChambers = 0;
  for (int ist = 0; ist < layout.nMSPlanes; ++ist) {
    const int nx = p.stationGridNX[ist];
    const int ny = p.stationGridNY[ist];
    if (nx <= 0 || ny <= 0) {
      throw std::runtime_error(fmt::format("MS{} has invalid MWPC grid {}x{}", ist, nx, ny));
    }
    if (nx * ny >= MaxNonSensID) {
      throw std::runtime_error(fmt::format("MS{} has too many chambers for local non-sensitive copy IDs", ist));
    }
    requestedChambers += nx * ny;
  }
  if (requestedChambers > MaxVolID - MaxNonSensID) {
    throw std::runtime_error(fmt::format(
      "MWPC layout requests {} sensitive chambers but the module supports at most {}",
      requestedChambers, MaxVolID - MaxNonSensID));
  }

  int chamberID = 0;
  for (int ist = 0; ist < layout.nMSPlanes; ++ist) {
    const int nx = p.stationGridNX[ist];
    const int ny = p.stationGridNY[ist];
    const bool narrow = ist == 0 && p.useNarrowMS0;
    const NA6PMWPCChamber& chamber = narrow ? narrowMS0Chamber : standardChamber;

    const auto localGasSize = chamber.gasFullSize();
    const auto gasCentre = chamber.gasCentre();
    const double gasX = swapXY ? localGasSize[1] : localGasSize[0];
    const double gasY = swapXY ? localGasSize[0] : localGasSize[1];
    const double overlapX = p.activeOverlapX;
    const double overlapY = p.activeOverlapY;
    const double pitchX = gasX - overlapX;
    const double pitchY = gasY - overlapY;

    if (!std::isfinite(overlapX) || !std::isfinite(overlapY) ||
        overlapX < 0. || overlapY < 0. || pitchX <= 0. || pitchY <= 0.) {
      throw std::runtime_error(fmt::format(
        "MS{} has invalid active overlap Ox={} Oy={} for detector gas {}x{} cm",
        ist, overlapX, overlapY, gasX, gasY));
    }
    if (!std::isfinite(p.staggerZStep) || p.staggerZStep <= 0.) {
      throw std::runtime_error("MWPC staggerZStep must be positive");
    }

    const double coverageX = gasX + (nx - 1) * pitchX;
    const double coverageY = gasY + (ny - 1) * pitchY;
    if (p.stationWorkingAreaX[ist] > 0.f && coverageX + 1.e-6 < p.stationWorkingAreaX[ist]) {
      throw std::runtime_error(fmt::format(
        "MS{} grid covers only {:.3f} cm in X, below requested working area {:.3f} cm",
        ist, coverageX, p.stationWorkingAreaX[ist]));
    }
    if (p.stationWorkingAreaY[ist] > 0.f && coverageY + 1.e-6 < p.stationWorkingAreaY[ist]) {
      throw std::runtime_error(fmt::format(
        "MS{} grid covers only {:.3f} cm in Y, below requested working area {:.3f} cm",
        ist, coverageY, p.stationWorkingAreaY[ist]));
    }

    const std::string stationName = fmt::format("MS{}", ist);
    auto* station = new TGeoVolumeAssembly(stationName.c_str());
    const int firstChamberID = chamberID;
    int localCopyID = 0;

    for (int row = 0; row < ny; ++row) {
      for (int col = 0; col < nx; ++col) {
        // Checkerboard parity in detector X/Y:
        //   A B A B ...
        //   C D C D ...
        //   A B A B ...
        // and the four parity classes are placed at four successive z levels.
        const int q = 2 * (row & 1) + (col & 1);
        const double x = (static_cast<double>(col) - 0.5 * (nx - 1)) * pitchX;
        const double y = (static_cast<double>(row) - 0.5 * (ny - 1)) * pitchY;
        const double desiredGasZ = (static_cast<double>(q) - 1.5) * p.staggerZStep;

        // addTo() first creates the chamber with its prototype-local axes.
        // Replace the node transform with a rotation+translation so that the
        // detailed asymmetric mechanical model is rotated as one rigid object;
        // this preserves the electronics/divider geometry rather than swapping
        // only the body dimensions.
        const double chamberZ = desiredGasZ - gasCentre[2];
        const NA6PMWPCChamber::Placement placement = {x, y, chamberZ};
        chamber.addTo(station, chamberID, localCopyID, placement);

        auto* chamberNode = dynamic_cast<TGeoNodeMatrix*>(station->GetNode(station->GetNdaughters() - 1));
        if (!chamberNode) {
          throw std::runtime_error(fmt::format("MS{} chamber {} node does not carry a placement matrix", ist, chamberID));
        }
        auto* rotation = new TGeoRotation();
        rotation->RotateZ(NA6PMWPCChamber::ChamberRotationZDeg);
        chamberNode->SetMatrix(new TGeoCombiTrans(x, y, chamberZ, rotation));

        LOGP(debug,
             "MS{} row {} col {} layer {} chamberID {} gas centre=({:.3f},{:.3f},{:.3f}) cm",
             ist, row, col, LayerName[q], chamberID, x, y, desiredGasZ);
        ++chamberID;
        ++localCopyID;
      }
    }

    const int lastChamberID = chamberID - 1;
    world->AddNode(
      station,
      composeNonSensorVolID(ist),
      new TGeoTranslation(layout.shiftMS[0] + layout.posMSPlaneX[ist],
                          layout.shiftMS[1] + layout.posMSPlaneY[ist],
                          layout.shiftMS[2] + layout.posMSPlaneZ[ist]));

    LOGP(info,
         "Created MS{} MWPC station: grid={}x{} N={} detectorGas={:.3f}x{:.3f} cm "
         "pitch={:.3f}x{:.3f} cm coverage={:.3f}x{:.3f} cm target={:.3f}x{:.3f} cm "
         "overlap={:.3f}x{:.3f} cm rotZ={:.1f} deg "
         "z(A,B,C,D)=({:.3f},{:.3f},{:.3f},{:.3f}) cm chamberID=[{},{}]{}",
         ist, nx, ny, nx * ny, gasX, gasY, pitchX, pitchY,
         coverageX, coverageY,
         p.stationWorkingAreaX[ist], p.stationWorkingAreaY[ist],
         overlapX, overlapY, NA6PMWPCChamber::ChamberRotationZDeg,
         -1.5 * p.staggerZStep, -0.5 * p.staggerZStep,
         +0.5 * p.staggerZStep, +1.5 * p.staggerZStep,
         firstChamberID, lastChamberID, narrow ? " [narrow MS0]" : "");
  }

  LOGP(info, "Created {} MWPC chambers in {} Muon Spectrometer stations",
       chamberID, layout.nMSPlanes);
}

void NA6PMuonSpec::setAlignableEntries()
{
  const auto& layout = NA6PLayoutParam::Instance();
  if (!gGeoManager || !gGeoManager->GetTopNode() || !gGeoManager->GetTopVolume()) {
    LOGP(error, "Cannot define MWPC alignable entries without a complete TGeo geometry");
    return;
  }

  const std::string topNodeName = gGeoManager->GetTopNode()->GetName();
  TGeoVolume* topVolume = gGeoManager->GetTopVolume();
  int nAlignable = 0;

  for (int ist = 0; ist < layout.nMSPlanes; ++ist) {
    const std::string stationName = fmt::format("MS{}", ist);
    TGeoNode* stationNode = nullptr;
    for (int inode = 0; inode < topVolume->GetNdaughters(); ++inode) {
      TGeoNode* candidate = topVolume->GetNode(inode);
      if (candidate && candidate->GetVolume() &&
          stationName == candidate->GetVolume()->GetName()) {
        stationNode = candidate;
        break;
      }
    }
    if (!stationNode) {
      LOGP(error, "Could not find station node {} below the top volume", stationName);
      continue;
    }

    TGeoVolume* stationVolume = stationNode->GetVolume();
    for (int ich = 0; ich < stationVolume->GetNdaughters(); ++ich) {
      TGeoNode* chamberNode = stationVolume->GetNode(ich);
      if (!chamberNode || !chamberNode->GetVolume()) {
        continue;
      }
      TGeoVolume* chamberVolume = chamberNode->GetVolume();
      TGeoNode* sensorNode = nullptr;
      for (int ipart = 0; ipart < chamberVolume->GetNdaughters(); ++ipart) {
        TGeoNode* partNode = chamberVolume->GetNode(ipart);
        if (partNode && NA6PModule::isSensor(partNode->GetNumber())) {
          if (sensorNode) {
            LOGP(error, "Chamber {} contains more than one sensitive node", chamberNode->GetName());
          }
          sensorNode = partNode;
        }
      }
      if (!sensorNode) {
        LOGP(error, "Chamber {} contains no sensitive gas node", chamberNode->GetName());
        continue;
      }

      const int sensorID = NA6PModule::volID2SensID(sensorNode->GetNumber());
      const int alignableID = getActiveID() * 1000 + sensorID;
      const std::string symbolicName = fmt::format("MS_Lr{}_Ch{}_Gas", ist, sensorID);
      const std::string path = fmt::format("/{}/{}/{}/{}",
                                           topNodeName,
                                           stationNode->GetName(),
                                           chamberNode->GetName(),
                                           sensorNode->GetName());
      TGeoPNEntry* entry = gGeoManager->SetAlignableEntry(
        symbolicName.c_str(), path.c_str(), alignableID);
      if (entry) {
        LOGP(debug, "Added alignable MWPC sensor {} path={} id={}",
             symbolicName, path, alignableID);
        ++nAlignable;
      } else {
        LOGP(error, "FAILED to add alignable MWPC sensor {} path={}", symbolicName, path);
      }
    }
  }
  LOGP(info, "Defined {} alignable MWPC sensitive-gas volumes", nAlignable);
}

bool NA6PMuonSpec::stepManager(int volID)
{
  int sensID = NA6PModule::volID2SensID(volID);
  if (sensID < 0) {
    LOGP(fatal, "Non-sensor volID={} was provided to stepManager of {}", volID, getName());
  }
  auto mc = TVirtualMC::GetMC();
  bool startHit = false, stopHit = false;
  unsigned char status = 0;
  if (mc->IsTrackEntering()) {
    status |= NA6PBaseHit::kTrackEntering;
  }
  if (mc->IsTrackInside()) {
    status |= NA6PBaseHit::kTrackInside;
  }
  if (mc->IsTrackExiting()) {
    status |= NA6PBaseHit::kTrackExiting;
  }
  if (mc->IsTrackOut()) {
    status |= NA6PBaseHit::kTrackOut;
  }
  if (mc->IsTrackStop()) {
    status |= NA6PBaseHit::kTrackStopped;
  }
  if (mc->IsTrackAlive()) {
    status |= NA6PBaseHit::kTrackAlive;
  }

  // track is entering or created in the volume
  if ((status & NA6PBaseHit::kTrackEntering) || (status & NA6PBaseHit::kTrackInside && !mTrackData.mHitStarted)) {
    startHit = true;
  } else if ((status & (NA6PBaseHit::kTrackExiting | NA6PBaseHit::kTrackOut | NA6PBaseHit::kTrackStopped))) {
    stopHit = true;
  }
  // increment energy loss at all steps except entrance
  if (!startHit) {
    mTrackData.mEnergyLoss += mc->Edep();
  }
  if (!(startHit | stopHit)) {
    return false; // do noting
  }
  auto stack = (NA6PMCStack*)mc->GetStack();
  if (startHit) {
    mTrackData.mEnergyLoss = 0.;
    mc->TrackMomentum(mTrackData.mMomentumStart);
    mc->TrackPosition(mTrackData.mPositionStart);
    mTrackData.mTrkStatusStart = status;
    mTrackData.mHitStarted = true;
  }
  if (stopHit) {
    TLorentzVector positionStop, momentumStop;
    mc->TrackMomentum(momentumStop);
    mc->TrackPosition(positionStop);
    // sensID is the global MWPC chamber ID assigned by createGeometry().
    auto* p = addHit(stack->GetCurrentTrackNumber(), sensID, mTrackData.mPositionStart.Vect(), positionStop.Vect(),
                     mTrackData.mMomentumStart.Vect(), momentumStop.Vect(), positionStop.T(),
                     mTrackData.mEnergyLoss, mTrackData.mTrkStatusStart, status);
    if (mVerbosity > 0) {
      LOGP(info, "{} Tr{} {}", getName(), stack->GetCurrentTrackNumber(), p->asString());
    }
    // register det points in TParticle
    stack->addHit(getActiveIDBit());
    mTrackData.mHitStarted = false;
    return true;
  }
  return false;
}

NA6PMuonSpecHit* NA6PMuonSpec::addHit(int trackID, int detID, const TVector3& startPos, const TVector3& endPos, const TVector3& startMom, const TVector3& endMom,
                                                    float endTime, float eLoss, unsigned char startStatus, unsigned char endStatus)
{
  mHits.emplace_back(trackID, detID, startPos, endPos, startMom, endMom, endTime, eLoss, startStatus, endStatus);
  return &(mHits.back());
}

void NA6PMuonSpec::createHitsOutput(const std::string& outDir)
{
  auto nm = fmt::format("{}Hits{}.root", outDir, getName());
  mHitsFile = TFile::Open(nm.c_str(), "recreate");
  mHitsTree = new TTree(fmt::format("hits{}", getName()).c_str(), fmt::format("{} Hits", getName()).c_str());
  mHitsTree->Branch(getName().c_str(), &hHitsPtr);
  LOGP(info, "Will store {} hits in {}", getName(), nm);
}

void NA6PMuonSpec::closeHitsOutput()
{
  if (mHitsTree && mHitsFile) {
    mHitsFile->cd();
    mHitsTree->Write();
    delete mHitsTree;
    mHitsTree = 0;
    mHitsFile->Close();
    delete mHitsFile;
    mHitsFile = 0;
  }
}

void NA6PMuonSpec::writeHits(const std::vector<int>& remapping)
{
  int nh = mHits.size();
  for (int i = 0; i < nh; i++) {
    auto& h = mHits[i];
    if (remapping[h.getTrackID()] < 0) {
      LOGP(error, "Track {} hit {} in {} was not remapped!", h.getTrackID(), i, getName());
    }
    h.setTrackID(remapping[h.getTrackID()]);
  }
  if (mHitsTree) {
    mHitsTree->Fill();
  }
}
