// NA6PCCopyright

#include "ConfigurableParam.h"
#include <fairlogger/Logger.h>
#include <TApplication.h>
#include <TMCManager.h>
#include <TGeoGlobalMagField.h>
#include <boost/program_options.hpp>
#include <filesystem>
#include <algorithm>
#include <TTree.h>
#include <TFile.h>
#include <TParticle.h>
#include <TStopwatch.h>
#include <utility>
#include "NA6PVerTelHit.h"
#include "NA6PMuonSpecHit.h"
#include "MagneticField.h"
#include "StringUtils.h"
#include "NA6PVerTelReconstruction.h"
#include "NA6PMuonSpecReconstruction.h"
#include "NA6PMatching.h"
#include "NA6PMCTruthContainer.h"

class TreeFromFile
{
 public:
  TreeFromFile(const char* fileName, const char* treeName)
  {
    file = TFile::Open(fileName);
    if (!file || file->IsZombie()) {
      LOGP(fatal, "Failed to open file {:p}", fileName);
    }
    tree = (TTree*)file->Get(treeName);
    if (!tree) {
      LOGP(fatal, "Failed to get tree {:p} from file {:p}", treeName, fileName);
    }
  }

  ~TreeFromFile()
  {
    delete tree;
    file->Close();
    delete file;
  }
  TTree* getTree() { return tree; }
  TFile* getFile() { return file; }

 private:
  TFile* file = nullptr;
  TTree* tree = nullptr;
};

double getPrimaryVertexZ(TTree* mcTree, std::vector<TParticle>* mcArr, int eventID)
{
  if (!mcTree || !mcArr) {
    return 0.0;
  }
  mcTree->GetEvent(eventID);
  for (const auto& part : *mcArr) {
    if (part.IsPrimary()) {
      return part.Vz();
    }
  }
  return 0.0;
}

int main(int argc, char** argv)
{
  const std::string paramsIni = "na6pRecoParam.ini";
  fair::Logger::OnFatal([]() { throw std::runtime_error("Fatal error"); });

  namespace bpo = boost::program_options;
  bpo::variables_map vm;
  bpo::options_description opt_general("Usage:\n  " + std::string(argv[0]));
  bpo::options_description opt_hidden("");
  bpo::options_description opt_all;
  bpo::positional_options_description opt_pos;

  try {
    auto add_option = opt_general.add_options();
    add_option("help,h", "Print this help message");
    add_option("verbosity,v", bpo::value<int>()->default_value(0), "verbosity level [0 = no output]");
    add_option("load-ini", bpo::value<std::string>()->default_value(""), "load configurables from ini file (if defined), overridden by configKeyValues");
    add_option("load-recoparam", bpo::value<std::string>()->default_value(""), "load reco parameters from ini file (if defined), overridden by configKeyValues");
    add_option("configKeyValues", bpo::value<std::string>()->default_value(""), "comma-separated configKeyValues");
    add_option("disable-write-ini", bpo::value<bool>()->default_value(false)->implicit_value(true), "do not write reco parameters ini file");
    add_option("geometry,g", bpo::value<std::string>()->default_value("geometry.root"), "geometry file name");
    add_option("firstevent,f", bpo::value<int32_t>()->default_value(0), "first event");
    add_option("lastevent,l", bpo::value<int32_t>()->default_value(-1), "last event (inclusive; -1 means end of file)");
    add_option("readMC", bpo::value<bool>()->default_value(true)->implicit_value(true), "read MC truth info");
    add_option("doHitsToRecPoints,hitcl", bpo::value<bool>()->default_value(true), "run hits->clusters");
    add_option("doDigitsToRecPoints,cl", bpo::value<bool>()->default_value(false), "run digits->clusters");
    add_option("doTrackletVertex,vert", bpo::value<bool>()->default_value(true), "run tracklet vertexer");
    add_option("doVTTracking,vt", bpo::value<bool>()->default_value(true), "run VT tracker");
    add_option("allowHolesInVTTracks,holesVT", bpo::value<bool>()->default_value(false), "allow holes in VT tracks");
    add_option("doMSTracking,ms", bpo::value<bool>()->default_value(true), "run MS tracker");
    add_option("doMatching,mt", bpo::value<bool>()->default_value(true), "run matching between VT and MS");
    opt_all.add(opt_general).add(opt_hidden);
    bpo::store(bpo::command_line_parser(argc, argv).options(opt_all).positional(opt_pos).run(), vm);

    if (vm.count("help")) {
      std::cout << opt_general << std::endl;
      exit(0);
    }

    bpo::notify(vm);
  } catch (bpo::error& e) {
    std::cerr << "ERROR: " << e.what() << std::endl
              << std::endl;
    std::cerr << opt_general << std::endl;
    exit(1);
  } catch (std::exception& e) {
    std::cerr << e.what() << ", application will now exit" << std::endl;
    exit(2);
  }

  auto flini = vm["load-ini"].as<std::string>();
  try {
    if (!flini.empty()) {
      na6p::conf::ConfigurableParam::updateFromFile(flini, "", true);
    }
    auto flrp = vm["load-recoparam"].as<std::string>();
    if (!flrp.empty()) {
      na6p::conf::ConfigurableParam::updateFromFile(flrp, "", true);
    }
    na6p::conf::ConfigurableParam::updateFromString(vm["configKeyValues"].as<std::string>()); // highest priority
  } catch (const std::exception& error) {
    std::cerr << "ERROR: Failed to load configuration: " << error.what() << std::endl;
    return 2;
  } catch (...) {
    std::cerr << "ERROR: Failed to load configuration: unknown error" << std::endl;
    return 2;
  }
  LOGP(info, "Printing all configs");
  na6p::conf::ConfigurableParam::printAllKeyValuePairs();

  const bool doHitsToRecPoints = vm["doHitsToRecPoints"].as<bool>();
  const bool doDigitsToRecPoints = vm["doDigitsToRecPoints"].as<bool>();
  const bool doTrackletVertex = vm["doTrackletVertex"].as<bool>();
  const bool doVTTracking = vm["doVTTracking"].as<bool>();
  const bool doMSTracking = vm["doMSTracking"].as<bool>();
  const bool doMatching = vm["doMatching"].as<bool>();
  const bool readMC = vm["readMC"].as<bool>();
  const bool skipVTlays = vm["allowHolesInVTTracks"].as<bool>();

  int firstEv = vm["firstevent"].as<int32_t>();
  int lastEv = vm["lastevent"].as<int32_t>();

  const auto getEventRange = [firstEv, lastEv](int nEvents) {
    const int first = std::clamp(firstEv, 0, nEvents);
    const int last = lastEv < 0 ? nEvents - 1 : std::min(lastEv, nEvents - 1);
    return std::pair<int, int>{first, last};
  };

  if (!Propagator::loadField() || !Propagator::loadGeometry(vm["geometry"].as<std::string>())) {
    return -1;
  }

  if (doHitsToRecPoints && doDigitsToRecPoints) {
    LOGP(info, "Will do digits to recpoints for VT instead of hits to recpoints");
  }

  std::unique_ptr<NA6PVerTelReconstruction> vtrec;
  if (doHitsToRecPoints || doDigitsToRecPoints || doTrackletVertex || doVTTracking) {
    vtrec = std::make_unique<NA6PVerTelReconstruction>();
    vtrec->setReadMCTruth(readMC);
    vtrec->setEnableTracksWithHoles(skipVTlays);
  }
  std::unique_ptr<NA6PMuonSpecReconstruction> msrec;
  if (doHitsToRecPoints || doMSTracking) {
    msrec = std::make_unique<NA6PMuonSpecReconstruction>();
    msrec->setReadMCTruth(readMC);
  }
  std::unique_ptr<NA6PMatching> matching = std::make_unique<NA6PMatching>();

  if (doHitsToRecPoints) {
    { // VTHits
      if (!doDigitsToRecPoints) {
        TreeFromFile tfVT("HitsVerTel.root", "hitsVerTel");
        std::vector<NA6PVerTelHit> vtHits, *vtHitsPtr = &vtHits;
        tfVT.getTree()->SetBranchAddress("VerTel", &vtHitsPtr);
        int nEvVT = tfVT.getTree()->GetEntriesFast();
        const auto [first, last] = getEventRange(nEvVT);

        vtrec->createClustersOutput();
        for (int jEv = first; jEv <= last; jEv++) {
          tfVT.getTree()->GetEvent(jEv);
          int nHits = vtHits.size();
          LOGP(info, "VerTel Event {} nHits= {}", jEv, nHits);
          vtrec->clearClusters();
          vtrec->hitsToRecPoints(vtHits, jEv);
          vtrec->writeClusters();
        }
        vtrec->closeClustersOutput();
      }
    }
    { // Muon Spectrometer hits -> clusters
      TreeFromFile tfMS("HitsMuonSpec.root", "hitsMuonSpec");
      std::vector<NA6PMuonSpecHit> msHits, *msHitsPtr = &msHits;
      tfMS.getTree()->SetBranchAddress("MuonSpec", &msHitsPtr);
      int nEvMS = tfMS.getTree()->GetEntriesFast();
      const auto [first, last] = getEventRange(nEvMS);

      msrec->createClustersOutput();
      for (int jEv = first; jEv <= last; jEv++) {
        tfMS.getTree()->GetEvent(jEv);
        int nHits = msHits.size();
        LOGP(info, "MuonSpec Event {} nHits= {}", jEv, nHits);
        msrec->clearClusters();
        msrec->hitsToRecPoints(msHits, jEv);
        msrec->writeClusters();
      }
      msrec->closeClustersOutput();
    }
  } else {
    LOGP(info, "Hits -> Recpoints disabled from input options");
  }
  if (doDigitsToRecPoints) {
    TreeFromFile tfVT("DigitsVerTel.root", "digitsVerTel");
    std::vector<NA6PVerTelDigit> vtDigits, *vtDigitsPtr = &vtDigits;
    NA6PMCTruthContainer vtDigMCLabels, *vtDigMCLabelsPtr = &vtDigMCLabels;
    tfVT.getTree()->SetBranchAddress("VerTel", &vtDigitsPtr);
    tfVT.getTree()->SetBranchAddress("VerTelMCTruth", &vtDigMCLabelsPtr);
    int nEvVT = tfVT.getTree()->GetEntriesFast();
    const auto [first, last] = getEventRange(nEvVT);
    vtrec->createClustersOutput();
    for (int jEv = first; jEv <= last; jEv++) {
      tfVT.getTree()->GetEvent(jEv);
      int nDigits = vtDigits.size();
      int nMClabels = vtDigMCLabels.getNElements();
      LOGP(info, "VerTel Event {} nDigits = {} nDigMClabels = {}", jEv, nDigits, nMClabels);
      vtrec->clearClusters();
      vtrec->digitsToRecPoints(vtDigits, vtDigMCLabels);
      vtrec->writeClusters();
    }
    vtrec->closeClustersOutput();
  }

  const bool needsMCKine = doTrackletVertex || doVTTracking || doMSTracking || doMatching;
  std::unique_ptr<TreeFromFile> tfKine;
  std::vector<TParticle>* mcArr = nullptr;
  if (needsMCKine) {
    tfKine = std::make_unique<TreeFromFile>("MCKine.root", "mckine");
    tfKine->getTree()->SetBranchAddress("tracks", &mcArr);
  }

  std::unique_ptr<TreeFromFile> tfCVT, tfCMS;
  std::vector<NA6PVerTelCluster> vtClus, *vtClusPtr = &vtClus;
  NA6PMCTruthContainer vtCluMCLabels, *vtCluMCLabelsPtr = &vtCluMCLabels;
  if (doTrackletVertex || doVTTracking) {
    tfCVT = std::make_unique<TreeFromFile>("ClustersVerTel.root", "clustersVerTel");
    tfCVT->getTree()->SetBranchAddress("VerTel", &vtClusPtr);
    tfCVT->getTree()->SetBranchAddress("VerTelMCTruth", &vtCluMCLabelsPtr);
  }
  std::vector<NA6PMuonSpecCluster> msClus, *msClusPtr = &msClus;
  NA6PMCTruthContainer msCluMCLabels, *msCluMCLabelsPtr = &msCluMCLabels;
  if (doMSTracking) {
    tfCMS = std::make_unique<TreeFromFile>("ClustersMuonSpec.root", "clustersMuonSpec");
    tfCMS->getTree()->SetBranchAddress("MuonSpec", &msClusPtr);
    tfCMS->getTree()->SetBranchAddress("MuonSpecMCTruth", &msCluMCLabelsPtr);
  }

  if (doTrackletVertex || doVTTracking || doMSTracking) {
    int nEv = -1;
    if (tfCVT && tfCVT->getTree()) {
      nEv = tfCVT->getTree()->GetEntries();
    }
    if (tfCMS && tfCMS->getTree()) {
      int nEvMS = tfCMS->getTree()->GetEntries();
      nEv = (nEv < 0) ? nEvMS : std::min(nEv, nEvMS);
    }
    if (nEv < 0) {
      LOGP(error, "No input trees for tracking stage");
      return -1;
    }
    const auto [first, last] = getEventRange(nEv);

    if (doTrackletVertex)
      vtrec->initVertexer();
    if (doVTTracking)
      vtrec->initTracker();
    if (doMSTracking)
      msrec->initTracker();

    TStopwatch timer;
    timer.Start();

    NA6PVertex pvert;
    for (int jEv = first; jEv <= last; jEv++) {
      LOGP(info, "Process event {}", jEv);
      const double zvert = getPrimaryVertexZ(tfKine->getTree(), mcArr, jEv);
      pvert.setXYZ(0.f, 0.f, zvert);

      if (doTrackletVertex || doVTTracking) {
        LOGP(info, "VT reconstruction");
        vtrec->clearEvent();
        tfCVT->getTree()->GetEvent(jEv);
        vtrec->setPrimaryVertex(&pvert);
        vtrec->setClusters(vtClus);
        vtrec->setClustersMCLabels(vtCluMCLabels);
        if (doTrackletVertex)
          vtrec->runVertexerTracklets();
        if (doVTTracking)
          vtrec->runTracking();
      }

      if (doMSTracking) {
        LOGP(info, "MS reconstruction");
        msrec->clearEvent();
        tfCMS->getTree()->GetEvent(jEv);
        msrec->setPrimaryVertex(&pvert);
        msrec->setClusters(msClus);
        msrec->setClustersMCLabels(msCluMCLabels);
        msrec->runTracking();
      }
    }

    if (doTrackletVertex || doVTTracking) {
      vtrec->closeVerticesOutput();
      vtrec->closeTracksOutput();
    }
    if (doMSTracking) {
      msrec->closeTracksOutput();
    }

    timer.Stop();
    LOGP(info, "------ Tracking time ------");
    timer.Print();
    LOGP(info, "------------------------------");
  }

  if (doMatching) {
    matching->initMatching();

    TreeFromFile tfTVT("TracksVerTel.root", "tracksVerTel");
    std::vector<NA6PTrack> vtTracks, *vtTracksPtr = &vtTracks;
    std::vector<NA6PMCComposedLabel> vtTrMCLabs, *vtTrMCLabsPtr = &vtTrMCLabs;
    tfTVT.getTree()->SetBranchAddress("VerTel", &vtTracksPtr);
    tfTVT.getTree()->SetBranchAddress("VerTelMCTruth", &vtTrMCLabsPtr);

    TreeFromFile tfTMS("TracksMuonSpec.root", "tracksMuonSpec");
    std::vector<NA6PTrack> msTracks, *msTracksPtr = &msTracks;
    std::vector<NA6PMCComposedLabel> msTrMCLabs, *msTrMCLabsPtr = &msTrMCLabs;
    tfTMS.getTree()->SetBranchAddress("MuonSpec", &msTracksPtr);
    tfTMS.getTree()->SetBranchAddress("MuonSpecMCTruth", &msTrMCLabsPtr);

    TreeFromFile tfCVT("ClustersVerTel.root", "clustersVerTel");
    std::vector<NA6PVerTelCluster> vtClusMatch, *vtClusMatchPtr = &vtClusMatch;
    tfCVT.getTree()->SetBranchAddress("VerTel", &vtClusMatchPtr);

    TreeFromFile tfCMS("ClustersMuonSpec.root", "clustersMuonSpec");
    std::vector<NA6PMuonSpecCluster> msClusMatch, *msClusMatchPtr = &msClusMatch;
    tfCMS.getTree()->SetBranchAddress("MuonSpec", &msClusMatchPtr);

    int nEvVT = tfTVT.getTree()->GetEntries();
    int nEvMS = tfTMS.getTree()->GetEntries();
    int nEv = std::min(nEvVT, nEvMS);
    const auto [first, last] = getEventRange(nEv);

    TStopwatch timer;
    timer.Start();

    NA6PVertex pvert;
    for (int jEv = first; jEv <= last; jEv++) {
      LOGP(info, "Process event {}", jEv);
      const double zvert = getPrimaryVertexZ(tfKine->getTree(), mcArr, jEv);
      pvert.setXYZ(0.f, 0.f, zvert);

      tfCVT.getTree()->GetEvent(jEv);
      tfCMS.getTree()->GetEvent(jEv);
      tfTVT.getTree()->GetEvent(jEv);
      tfTMS.getTree()->GetEvent(jEv);

      matching->setPrimaryVertex(&pvert);
      matching->setVerTelClusters(vtClusMatch);
      matching->setMuonSpecClusters(msClusMatch);
      matching->setVerTelTracks(vtTracks);
      matching->setVerTelTrackMCLabels(vtTrMCLabs);
      matching->setMuonSpecTracks(msTracks);
      matching->setMuonSpecTrackMCLabels(msTrMCLabs);

      matching->runMatching();
    }
    matching->closeTracksOutput();

    timer.Stop();
    LOGP(info, "------ Matching time ------");
    timer.Print();
    LOGP(info, "------------------------------");
  }

  if (!vm["disable-write-ini"].as<bool>()) {
    try {
      std::string pth = na6p::utils::Str::rectifyDirectory(na6p::conf::ConfigurableParam::getOutputDir());
      if (pth.empty()) {
        pth = std::filesystem::current_path();
      }
      na6p::conf::ConfigurableParam::writeINI(paramsIni, "reco");
      LOGP(info, "Stored configurable params to {}/{}", pth, paramsIni);
    } catch (std::exception e) {
      LOGP(error, "Failed to store configurable params, the reason is {}", e.what());
    }
  }

  return 0;
}
