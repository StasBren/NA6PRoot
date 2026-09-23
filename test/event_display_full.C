#include <TEveManager.h>
#include <TEveEventManager.h>
#include <TEveLine.h>
#include <TEveGeoNode.h>
#include <TEvePointSet.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TParticle.h>
#include <TGeoManager.h>
#include <iostream>
#include <vector>
#include <TEveSelection.h>
#include <iostream>
#include <TQObject.h>
#include "NA6PBaseCluster.h"
#include "NA6PMuonSpecHit.h"
#include "NA6PMuonSpecCluster.h"
#include "NA6PVerTelHit.h"
#include "NA6PVerTelCluster.h"
#include "NA6PTrack.h"
#include "NA6PMatch.h"
#include "NA6PMCTruthContainer.h"
#include "NA6PMCComposedLabel.h"
#include "NA6PFastTrackFitter.h"
#include "Propagator.h"
#include "NA6PLayoutParam.h"
#include "MagneticField.h"

void ResetClusterHighlight(int);
void HighlightClusterIndex(int, int);

// ---------- GLOBALS ----------
// cluster index --> EVE object
std::vector<std::vector<TEvePointSet*>> gClusterEve;
std::vector<std::vector<Color_t>> gClusterOrigColor;
std::vector<std::vector<Float_t>> gClusterOrigSize;

std::vector<std::unique_ptr<NA6PMuonSpecHit>> gAllHitsM;
std::vector<std::unique_ptr<NA6PVerTelHit>> gAllHitsV;
std::vector<std::unique_ptr<NA6PTrack>> gAllTracksM;
std::vector<std::unique_ptr<NA6PTrack>> gAllTracksV;
std::vector<std::unique_ptr<NA6PMatch>> gAllTracksT;
std::vector<std::unique_ptr<NA6PBaseCluster>> gAllClustersM;
std::vector<std::unique_ptr<NA6PBaseCluster>> gAllClustersV;
std::vector<std::vector<const NA6PBaseCluster*>> gEventClustersM;
std::vector<std::vector<const NA6PBaseCluster*>> gEventClustersV;
std::vector<std::vector<std::unique_ptr<TParticle>>> gMCParticles;

// ------------------------------------------------------
// Define struct for various elements
// ------------------------------------------------------

enum class DetectorType { Muon,
                          Vertex,
                          Matched };

enum class EveObjType { Hit,
                        Cluster,
                        TrackM,
                        TrackV,
                        TrackT };

struct EveUserData : public TObject {
  EveObjType type;
  const void* ptr;
  int event;
  int idPart;
  std::string mcString;

  EveUserData(EveObjType t, const void* p, int e, int jp, std::string lab = "")
    : type(t), ptr(p), event(e), idPart(jp), mcString(lab) {}

  ClassDef(EveUserData, 1);
};

struct DetectorConfig {
  const char* label;

  int layerMin;
  int layerMax;

  DetectorType det;

  std::vector<std::unique_ptr<NA6PBaseCluster>>* allClusters;
  std::vector<std::vector<const NA6PBaseCluster*>>* eventClusters;

  Color_t hitColor;
  Color_t clusterColor;
};

DetectorConfig kMuonSpec{
  "MS TRACK",
  5, 11,
  DetectorType::Muon,
  &gAllClustersM,
  &gEventClustersM,
  kRed,
  kGreen + 2};

DetectorConfig kVerTel{
  "VT TRACK",
  0, 5,
  DetectorType::Vertex,
  &gAllClustersV,
  &gEventClustersV,
  kMagenta,
  kCyan + 2};

DetectorConfig kMatched{
  "MATCHED TRACK",
  0, 11,
  DetectorType::Matched,
  nullptr,
  nullptr,
  kOrange,
  kViolet + 2};

class EvePicker : public TObject
{
 public:
  template <typename TRK>
  void DumpTrack(const TRK* tr,
                 int ev,
                 std::string mcstr,
                 int mcpart,
                 const DetectorConfig& cfg)
  {
    auto xyz = tr->template getXYZ<double>();
    auto pxyz = tr->template getPXYZ<double>();

    double p = tr->getP();
    double pt = TMath::Sqrt(pxyz[0] * pxyz[0] + pxyz[1] * pxyz[1]);
    double eta = -TMath::Log(pt / (p + pxyz[2]));
    int nHits = -1;
    if (std::is_same<TRK, NA6PTrack>::value) {
      nHits = tr->getNHits();
    }
    std::cout << "[" << cfg.label << "]\n"
              << "  event " << ev << "\n"
              << "  n_clusters = " << nHits << "\n"
              << "  track origin X = " << xyz[0]
              << " Y = " << xyz[1]
              << " Z = " << xyz[2] << "\n"
              << "  p = " << p
              << " pt = " << pt
              << " eta = " << eta << "\n"
              << "  Particle ID = " << mcstr.c_str()
              << std::endl;

    ResetClusterHighlight(ev);

    std::cout << "  Associated clusters:\n";

    if (std::is_same<TRK, NA6PTrack>::value) {
      for (int lr = cfg.layerMin; lr < cfg.layerMax; ++lr) {
        int ic = tr->getClusterIndex(lr);
        if (ic < 0)
          continue;

        HighlightClusterIndex(ev, ic);

        const std::vector<const NA6PBaseCluster*>* evClustersPtr = nullptr;

        if (cfg.det == DetectorType::Muon) {
          evClustersPtr = &gEventClustersM[ev];
        } else if (cfg.det == DetectorType::Vertex) {
          evClustersPtr = &gEventClustersV[ev];
        } else if (cfg.det == DetectorType::Matched) {

          // layer-based routing
          if (lr < 5)
            evClustersPtr = &gEventClustersV[ev];
          else
            evClustersPtr = &gEventClustersM[ev];
        }

        if (!evClustersPtr)
          continue;

        const auto& evClusters = *evClustersPtr;

        if (ic >= (int)evClusters.size())
          continue;

        const auto* clu = evClusters[ic];

        std::cout
          << "    layer " << clu->getLayer()
          << "  size " << clu->getClusterSize()
          << "  X " << clu->getX()
          << "  Y " << clu->getY()
          << "  Z " << clu->getZ()
          << "\n";
      }
    }
    // ---- MC printout
    int n = TMath::Abs(mcpart);

    if (ev >= 0 &&
        ev < (int)gMCParticles.size() &&
        n >= 0 &&
        n < (int)gMCParticles[ev].size()) {

      const TParticle* p = gMCParticles[ev][n].get();

      if (mcpart < 0)
        std::cout << "WARNING: cluster mismatch\n";

      std::cout
        << "[MC PARTICLE (index = " << n << ")]\n"
        << "  PDG = " << p->GetPdgCode() << "\n"
        << "  p = (" << p->Px()
        << ", " << p->Py()
        << ", " << p->Pz() << ")\n"
        << "  E = " << p->Energy() << "\n"
        << "  vtx = ("
        << p->Vx() << ", "
        << p->Vy() << ", "
        << p->Vz() << ")"
        << std::endl;
    } else {
      std::cout << "[MC PARTICLE] invalid ParticleID = "
                << n << std::endl;
    }

    std::cout << std::endl;
    gEve->Redraw3D(kFALSE);
  }

  void OnSelectionAdded(TEveElement* el)
  {
    if (!el)
      return;

    auto* ud = static_cast<EveUserData*>(el->GetUserData());
    if (!ud)
      return;

    switch (ud->type) {

      case EveObjType::Hit: {

        auto* hit =
          static_cast<const NA6PBaseHit*>(ud->ptr);

        std::cout
          << "[HIT]\n"
          << "  X = " << hit->getX()
          << "  Y = " << hit->getY()
          << "  Z = " << hit->getZ()
          << std::endl;

        break;
      }
      case EveObjType::Cluster: {
        auto* clu =
          static_cast<const NA6PBaseCluster*>(ud->ptr);
        std::cout
          << "[CLUSTER]\n"
          << "  layer = " << clu->getLayer()
          << "  size  = " << clu->getClusterSize()
          << "\n"
          << "  X = " << clu->getX()
          << "  Y = " << clu->getY()
          << "  Z = " << clu->getZ()
          << std::endl;
        break;
      }

      case EveObjType::TrackM: {
        auto* tr = static_cast<const NA6PTrack*>(ud->ptr);
        DumpTrack(tr, ud->event, ud->mcString, ud->idPart, kMuonSpec);
        break;
      }

      case EveObjType::TrackV: {
        auto* tr = static_cast<const NA6PTrack*>(ud->ptr);
        DumpTrack(tr, ud->event, ud->mcString, ud->idPart, kVerTel);
        break;
      }

      case EveObjType::TrackT: {
        auto* tr = static_cast<const NA6PTrack*>(ud->ptr);
        DumpTrack(tr, ud->event, ud->mcString, ud->idPart, kMatched);
        break;
      }
    }
  }
  ClassDef(EvePicker, 0);
};

template <typename HitType>
void DrawHits(const std::vector<HitType>& hits,
              std::vector<std::unique_ptr<HitType>>& storage,
              int evIdx,
              Color_t color)
{
  for (const auto& h : hits) {

    storage.emplace_back(
      std::make_unique<HitType>(h));

    auto* hitPtr = storage.back().get();

    auto* hitEve = new TEvePointSet(1);
    hitEve->SetNextPoint(h.getX(), h.getY(), h.getZ());
    hitEve->SetMarkerStyle(4);
    hitEve->SetMarkerSize(1.);
    hitEve->SetMarkerColor(color);
    hitEve->SetPickable(kTRUE);

    hitEve->SetUserData(
      new EveUserData(
        EveObjType::Hit,
        static_cast<NA6PBaseHit*>(hitPtr),
        evIdx,
        h.getTrackID(),
        ""));

    gEve->AddElement(hitEve);
  }
}

template <typename ClusterType>
void DrawClusters(const std::vector<ClusterType>& clusters,
                  const NA6PMCTruthContainer& cluMCLabels,
                  DetectorConfig& cfg,
                  int evIdx)
{
  cfg.eventClusters->emplace_back();

  int nClusters = clusters.size();
  for (int j = 0; j < nClusters; ++j) {
    auto& c = clusters[j];
    std::span labels = cluMCLabels.getLabels(j);
    std::string cluLabels = "";
    int nLabels = labels.size();
    for (int jLab = 0; jLab < nLabels; jLab++) {
      if (nLabels > 1)
        cluLabels.append(Form("(%d)", jLab));
      NA6PMCComposedLabel lbl = labels[jLab];
      cluLabels.append(lbl.asString());
      if (nLabels > 1)
        cluLabels.append(",");
    }

    cfg.allClusters->emplace_back(
      std::make_unique<ClusterType>(c));

    auto* cluPtr = cfg.allClusters->back().get();

    auto* eveClu = new TEvePointSet(1);
    eveClu->SetNextPoint(
      cluPtr->getX(),
      cluPtr->getY(),
      cluPtr->getZ());

    eveClu->SetMarkerStyle(20);
    eveClu->SetMarkerSize(1.5);
    eveClu->SetMarkerColor(cfg.clusterColor);
    eveClu->SetPickable(kTRUE);

    (*cfg.eventClusters)[evIdx].push_back(cluPtr);

    eveClu->SetUserData(
      new EveUserData(EveObjType::Cluster, cluPtr, evIdx, 9999999, cluLabels));

    gClusterEve[evIdx].push_back(eveClu);
    gClusterOrigColor[evIdx].push_back(eveClu->GetMarkerColor());
    gClusterOrigSize[evIdx].push_back(eveClu->GetMarkerSize());

    gEve->AddElement(eveClu);
  }
}

void event_display_full(int firstEv = 0, int nEv = 1,
                        const char* fgeo = "geometry.root",
                        const char* fini = "na6pLayout.ini",
                        const char* fkine = "MCKine.root",
                        const char* fhitsM = "HitsMuonSpec.root",
                        const char* fclustersM = "ClustersMuonSpec.root",
                        const char* ftracksM = "TracksMuonSpec.root",
                        const char* fhitsV = "HitsVerTel.root",
                        const char* fclustersV = "ClustersVerTel.root",
                        const char* ftracksV = "TracksVerTel.root",
                        const char* ftracksT = "TracksMatching.root")
{
  // ------------------------------------------------------
  // EVE initialization (only once)
  // ------------------------------------------------------
  if (!gEve) {
    TEveManager::Create();
  }

  // ------------------------------------------------------
  // Geometry (only once)
  // ------------------------------------------------------
  static bool geomLoaded = false;
  if (!geomLoaded) {
    if (!Propagator::loadGeometry(fgeo) || !Propagator::loadField()) {
      return;
    }

    TEveGeoTopNode* geom =
      new TEveGeoTopNode(gGeoManager, gGeoManager->GetTopNode());
    geom->SetPickableRecursively(kFALSE);
    gEve->AddGlobalElement(geom);

    geomLoaded = true;
  }

  auto& param = NA6PLayoutParam::Instance();
  if (fini)
    param.updateFromFile(fini, "", true);
  const int nVTstations = param.nVerTelPlanes;
  const int nMSstations = param.nMSPlanes;

  // ------------------------------------------------------
  // Clear previous event
  // ------------------------------------------------------
  if (gEve->GetCurrentEvent()) {
    gEve->GetCurrentEvent()->DestroyElements();
  }

  // ------------------------------------------------------
  // Create NEW event manager
  // ------------------------------------------------------
  TEveEventManager* ev =
    new TEveEventManager("Event", "NA6 event");
  gEve->AddEvent(ev);
  gEve->SetCurrentEvent(ev);

  // ------------------------------------------------------
  // Clear your persistent storage
  // ------------------------------------------------------
  gAllHitsM.clear();
  gAllTracksM.clear();
  gAllClustersM.clear();
  gEventClustersM.clear();
  gAllHitsV.clear();
  gAllTracksV.clear();
  gAllClustersV.clear();
  gEventClustersV.clear();
  gAllTracksT.clear();

  gClusterEve.clear();
  gClusterOrigColor.clear();
  gClusterOrigSize.clear();

  gMCParticles.clear();

  static EvePicker* gPicker = new EvePicker;

  gEve->GetSelection()->Disconnect(
    "SelectionAdded(TEveElement*)",
    gPicker,
    "OnSelectionAdded(TEveElement*)");
  gEve->GetSelection()->Connect(
    "SelectionAdded(TEveElement*)",
    "EvePicker",
    gPicker,
    "OnSelectionAdded(TEveElement*)");

  // Enable selection
  gEve->GetSelection()->SetPickToSelect(kTRUE);

  // Load kine
  TFile* fk = new TFile(fkine);
  TTree* tk = (TTree*)fk->Get("mckine");

  std::vector<TParticle>* mcArr = nullptr;
  tk->SetBranchAddress("tracks", &mcArr);

  // Load hits
  TFile* fM = TFile::Open(fhitsM);
  TTree* tM = (TTree*)fM->Get("hitsMuonSpec");
  TFile* fV = TFile::Open(fhitsV);
  TTree* tV = (TTree*)fV->Get("hitsVerTel");

  std::vector<NA6PMuonSpecHit> msHits;
  std::vector<NA6PMuonSpecHit>* msHitsPtr = &msHits;
  tM->SetBranchAddress("MuonSpec", &msHitsPtr);
  std::vector<NA6PVerTelHit> vtHits;
  std::vector<NA6PVerTelHit>* vtHitsPtr = &vtHits;
  tV->SetBranchAddress("VerTel", &vtHitsPtr);

  // Load clusters
  TFile* fcM = TFile::Open(fclustersM);
  TTree* tcM = (TTree*)fcM->Get("clustersMuonSpec");
  TFile* fcV = TFile::Open(fclustersV);
  TTree* tcV = (TTree*)fcV->Get("clustersVerTel");

  std::vector<NA6PMuonSpecCluster> clustersM;
  std::vector<NA6PMuonSpecCluster>* clustersPtrM = &clustersM;
  NA6PMCTruthContainer cluMCLabelsM, *cluMCLabelsPtrM = &cluMCLabelsM;
  tcM->SetBranchAddress("MuonSpec", &clustersPtrM);
  tcM->SetBranchAddress("MuonSpecMCTruth", &cluMCLabelsPtrM);
  std::vector<NA6PVerTelCluster> clustersV;
  std::vector<NA6PVerTelCluster>* clustersPtrV = &clustersV;
  NA6PMCTruthContainer cluMCLabelsV, *cluMCLabelsPtrV = &cluMCLabelsV;
  tcV->SetBranchAddress("VerTel", &clustersPtrV);
  tcV->SetBranchAddress("VerTelMCTruth", &cluMCLabelsPtrV);

  // Load tracks

  TFile* ftM = new TFile(ftracksM);
  TTree* ttM = (TTree*)ftM->Get("tracksMuonSpec");
  TFile* ftV = new TFile(ftracksV);
  TTree* ttV = (TTree*)ftV->Get("tracksVerTel");
  TFile* ftT = new TFile(ftracksT);
  TTree* ttT = (TTree*)ftT->Get("tracksMatching");

  std::vector<NA6PTrack> msTracks, *msTracksPtr = &msTracks;
  std::vector<NA6PMCComposedLabel> msTrackLabs, *msTrackLabPtr = &msTrackLabs;
  ttM->SetBranchAddress("MuonSpec", &msTracksPtr);
  ttM->SetBranchAddress("MuonSpecMCTruth", &msTrackLabPtr);
  std::vector<NA6PTrack> vtTracks, *vtTracksPtr = &vtTracks;
  std::vector<NA6PMCComposedLabel> vtTrackLabs, *vtTrackLabPtr = &vtTrackLabs;
  ttV->SetBranchAddress("VerTel", &vtTracksPtr);
  ttV->SetBranchAddress("VerTelMCTruth", &vtTrackLabPtr);
  std::vector<NA6PMatch> maTracks, *maTracksPtr = &maTracks;
  std::vector<NA6PMCComposedLabel> maTrackLabs, *maTrackLabPtr = &maTrackLabs;
  ttT->SetBranchAddress("Matching", &maTracksPtr);
  ttT->SetBranchAddress("MatchingMCTruth", &maTrackLabPtr);

  for (int iEv = firstEv; iEv < firstEv + nEv; iEv++) {

    int evIdx = iEv - firstEv;

    tM->GetEntry(iEv);
    tcM->GetEntry(iEv);
    ttM->GetEntry(iEv);
    tV->GetEntry(iEv);
    tcV->GetEntry(iEv);
    ttV->GetEntry(iEv);
    ttT->GetEntry(iEv);
    tk->GetEntry(iEv);

    gClusterEve.emplace_back();
    gClusterOrigColor.emplace_back();
    gClusterOrigSize.emplace_back();
    gEventClustersM.emplace_back();
    gEventClustersV.emplace_back();

    gMCParticles.emplace_back();

    auto& mcEv = gMCParticles.back();

    for (const auto& p : *mcArr) {
      mcEv.emplace_back(
        std::make_unique<TParticle>(p));
    }

    DrawHits(msHits, gAllHitsM, evIdx, kMuonSpec.hitColor);
    DrawHits(vtHits, gAllHitsV, evIdx, kVerTel.hitColor);

    DrawClusters(clustersM, cluMCLabelsM, kMuonSpec, evIdx);
    DrawClusters(clustersV, cluMCLabelsV, kVerTel, evIdx);

    std::cout << " Number of reconstructed tracks MS " << msTracks.size() << std::endl;

    int nMSTracks = msTracks.size();
    std::cout << " Number of reconstructed MS tracks " << nMSTracks << std::endl;

    for (int j = 0; j < nMSTracks; ++j) {
      auto& track = msTracks[j];
      NA6PMCComposedLabel lbl = msTrackLabs[j];

      gAllTracksM.emplace_back(std::make_unique<NA6PTrack>(track));
      auto* trPtr = gAllTracksM.back().get();

      TEveLine* eveTrack = new TEveLine();
      eveTrack->SetLineColor(kBlue);
      eveTrack->SetLineWidth(2);
      eveTrack->SetPickable(kTRUE);

      for (int np = (int)param.posMSPlaneZ[0] - 10; np < (int)param.posMSPlaneZ[nMSstations - 1] + 10; np++) {
        Propagator::Instance()->propagateToZ(track, 0. + np);
        eveTrack->SetNextPoint(track.getX(), track.getY(), track.getZ());
      }

      eveTrack->SetUserData(
        new EveUserData(EveObjType::TrackM, trPtr, evIdx, lbl.getTrackIDSigned(), lbl.asString()));

      gEve->AddElement(eveTrack);
    }

    std::cout << " Number of reconstructed tracks VT " << vtTracks.size() << std::endl;

    int nVTTracks = vtTracks.size();
    std::cout << " Number of reconstructed VT tracks " << nVTTracks << std::endl;

    for (int j = 0; j < nVTTracks; ++j) {
      auto& track = vtTracks[j];
      NA6PMCComposedLabel lbl = vtTrackLabs[j];

      gAllTracksV.emplace_back(std::make_unique<NA6PTrack>(track));
      auto* trPtr = gAllTracksV.back().get();

      TEveLine* eveTrack = new TEveLine();
      eveTrack->SetLineColor(kRed);
      eveTrack->SetLineWidth(1);
      eveTrack->SetPickable(kTRUE);

      for (int np = (int)param.posVerTelPlaneZ[0] - 1; np < (int)param.posVerTelPlaneZ[nVTstations - 1] + 2; np++) {
        Propagator::Instance()->propagateToZ(track, 0. + np);
        eveTrack->SetNextPoint(track.getX(), track.getY(), track.getZ());
      }

      eveTrack->SetUserData(
        new EveUserData(EveObjType::TrackV, trPtr, evIdx, lbl.getTrackIDSigned(), lbl.asString()));

      gEve->AddElement(eveTrack);
    }

    int nMATracks = maTracks.size();
    std::cout << " Number of matched tracks " << nMATracks << std::endl;

    for (int j = 0; j < nMATracks; ++j) {
      auto& track = maTracks[j];
      NA6PMCComposedLabel lbl = maTrackLabs[j];

      gAllTracksT.emplace_back(std::make_unique<NA6PMatch>(track));
      auto* trPtr = gAllTracksT.back().get();

      TEveLine* eveTrack = new TEveLine();
      eveTrack->SetLineColor(kOrange);
      eveTrack->SetLineWidth(3);
      eveTrack->SetPickable(kTRUE);

      for (int np = 0; np < (int)param.posMSPlaneZ[nMSstations - 1] + 20; np++) {
        Propagator::Instance()->propagateToZ(track, 0. + np);
        eveTrack->SetNextPoint(track.getX(), track.getY(), track.getZ());
      }

      eveTrack->SetUserData(
        new EveUserData(EveObjType::TrackT, trPtr, evIdx, lbl.getTrackIDSigned(), lbl.asString()));

      gEve->AddElement(eveTrack);
    }

    gEve->Redraw3D(kTRUE);
  }
}

void ResetClusterHighlight(int ev)
{
  for (size_t i = 0; i < gClusterEve[ev].size(); ++i) {
    gClusterEve[ev][i]->SetMarkerColor(gClusterOrigColor[ev][i]);
    gClusterEve[ev][i]->SetMarkerSize(gClusterOrigSize[ev][i]);
  }
}

void HighlightClusterIndex(int ev, int idx)
{
  if (ev < 0 || ev >= (int)gClusterEve[ev].size())
    return;
  if (idx < 0 || idx >= (int)gClusterEve[ev].size())
    return;

  gClusterEve[ev][idx]->SetMarkerColor(kYellow);
  gClusterEve[ev][idx]->SetMarkerSize(2.5);
}
