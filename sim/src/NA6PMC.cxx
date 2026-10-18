// NA6PCCopyright

#include "NA6PMC.h"
#include "NA6PSimMisc.h"
#include "NA6PDetector.h"
#include "NA6PModule.h"
#include "NA6PGenerator.h"
#include "NA6PMCStack.h"
#include "NA6PMuonReferenceAudit.h"
#include "StringUtils.h"
#include "MiscUtils.h"
#include "KeyValParam.h"
#include "NA6PLayoutParam.h"
#include "NA6PTGeoHelper.h"
#include "NA6PGenCutParam.h"
#include <TSystem.h>
#include <TGeoManager.h>
#include <TVirtualMC.h>
#include <TMCManager.h>
#include <TVirtualMCStack.h>
#include <TGeoVolume.h>
#include <TMCProcess.h>
#include <TInterpreter.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TMethodCall.h>
#include <fairlogger/Logger.h>
#include <filesystem>
#include <functional>
#include <regex>
#include <unistd.h>
#include <cstdlib>

using Str = na6p::utils::Str;

NA6PMC::NA6PMC(const char* name, const char* title) : TVirtualMCApplication(name, title)
{
  NA6PDetector::instance();
  mStack = std::make_unique<NA6PMCStack>(200);
  LOGP(info, "NA6PMC application initialized.");
}

NA6PMC::~NA6PMC()
{
  // Complete the optional reference-track ROOT output before shutdown.
  mMuonReferenceAudit.reset();
  closeKineOutput();
  NA6PDetector::instance().closeHitsOutput();
  LOGP(info, "NA6PMC application terminated.");
}

void NA6PMC::ConstructGeometry()
{
  LOGP(info, "Constructing geometry");
  NA6PDetector::instance().createGeometry();
  TVirtualMC::GetMC()->SetRootGeometry();
}

void NA6PMC::InitGeometry()
{
  const auto& param = NA6PLayoutParam::Instance();
  if (!param.materialsCutsFile.empty()) {
    NA6PTGeoHelper::instance().loadCutsAndProcessesFromFile(gSystem->ExpandPathName(param.materialsCutsFile.c_str()));
  }
}

void NA6PMC::ConstructOpGeometry()
{
  LOGP(info, "Constructing optional geometry...");
}

int NA6PMC::callUserHook(int hookID, bool inout)
{
  int ret = 0;
  if (mUserHooksMethod) {
    const void* args[2] = {&hookID, &inout};
    mUserHooksMethod->Execute(nullptr, args, 2, &ret);
    if (ret < 0) {
      LOGP(fatal, "Execution of user hook method {}({},{}) failed with {}", mUserHookName, hookID, inout, ret);
    }
  }
  return ret;
}

void NA6PMC::AddParticles()
{
  if (mVerbosity > 0) {
    LOGP(info, "Adding particles to the simulation");
  }
  callUserHook(UserHook::ADDParticles, true); // call at entry
  //
  addSpecialParticles(); // copied from O2
  forceCharmHadronicDecays();
  forceJpsiDecays();
  callUserHook(UserHook::ADDParticles, false); // call at entry
}

void NA6PMC::forceCharmHadronicDecays()
{
  int D0_decay_mode[6][3] = {{0}};
  float D0_BR[6] = {100.f, 0.f, 0.f, 0.f, 0.f, 0.f};
  D0_decay_mode[0][0] = -321; // K-
  D0_decay_mode[0][1] = 211;  // pi+
  TVirtualMC::GetMC()->SetDecayMode(421, D0_BR, D0_decay_mode);
  int D0bar_decay_mode[6][3] = {{0}};
  D0bar_decay_mode[0][0] = 321;  // K+
  D0bar_decay_mode[0][1] = -211; // pi-
  TVirtualMC::GetMC()->SetDecayMode(-421, D0_BR, D0bar_decay_mode);
  //
  int Dp_decay_mode[6][3] = {{0}};
  float Dp_BR[6] = {100.f, 0.f, 0.f, 0.f, 0.f, 0.f};
  Dp_decay_mode[0][0] = -321; // K-
  Dp_decay_mode[0][1] = 211;  // pi+
  Dp_decay_mode[0][2] = 211;  // pi+
  TVirtualMC::GetMC()->SetDecayMode(411, Dp_BR, Dp_decay_mode);
  int Dm_decay_mode[6][3] = {{0}};
  Dm_decay_mode[0][0] = 321;  // K+
  Dm_decay_mode[0][1] = -211; // pi-
  Dm_decay_mode[0][2] = -211; // pi-
  TVirtualMC::GetMC()->SetDecayMode(-411, Dp_BR, Dm_decay_mode);
  //
  int Dsp_decay_mode[6][3] = {{0}};
  float Ds_BR[6] = {90.f, 10.f, 0.f, 0.f, 0.f, 0.f};
  Dsp_decay_mode[0][0] = 333;  // phi
  Dsp_decay_mode[0][1] = 211;  // pi+
  Dsp_decay_mode[1][0] = -313; // K*0bar
  Dsp_decay_mode[1][1] = 321;  // K+
  TVirtualMC::GetMC()->SetDecayMode(431, Ds_BR, Dsp_decay_mode);
  int Dsm_decay_mode[6][3] = {{0}};
  Dsm_decay_mode[0][0] = 333;  // phi
  Dsm_decay_mode[0][1] = -211; // pi-
  Dsm_decay_mode[1][0] = 313;  // K*0
  Dsm_decay_mode[1][1] = -321; // K-
  TVirtualMC::GetMC()->SetDecayMode(-431, Ds_BR, Dsm_decay_mode);
  int Kstar_decay_mode[6][3] = {{0}};
  float Kstar_BR[6] = {100.f, 0.f, 0.f, 0.f, 0.f, 0.f};
  Kstar_decay_mode[0][0] = 321;  // K+
  Kstar_decay_mode[0][1] = -211; // pi-
  TVirtualMC::GetMC()->SetDecayMode(313, Kstar_BR, Kstar_decay_mode);
  int Kstarbar_decay_mode[6][3] = {{0}};
  Kstarbar_decay_mode[0][0] = -321; // K-
  Kstarbar_decay_mode[0][1] = 211;  // pi+
  TVirtualMC::GetMC()->SetDecayMode(-313, Kstar_BR, Kstarbar_decay_mode);
  int Phi_decay_mode[6][3] = {{0}};
  float Phi_BR[6] = {100.f, 0.f, 0.f, 0.f, 0.f, 0.f};
  Phi_decay_mode[0][0] = 321;  // K+
  Phi_decay_mode[0][1] = -321; // K-
  TVirtualMC::GetMC()->SetDecayMode(333, Phi_BR, Phi_decay_mode);
  //
  int Lcp_decay_mode[6][3] = {{0}};
  float Lc_BR[6] = {100.f, 0.f, 0.f, 0.f, 0.f, 0.f};
  Lcp_decay_mode[0][0] = -321; // K-
  Lcp_decay_mode[0][1] = 2212; // p
  Lcp_decay_mode[0][2] = 211;  // pi+
  TVirtualMC::GetMC()->SetDecayMode(4122, Lc_BR, Lcp_decay_mode);
  int Lcm_decay_mode[6][3] = {{0}};
  Lcm_decay_mode[0][0] = 321;   // K+
  Lcm_decay_mode[0][1] = -2212; // pbar
  Lcm_decay_mode[0][2] = -211;  // pi-
  TVirtualMC::GetMC()->SetDecayMode(-4122, Lc_BR, Lcm_decay_mode);
}

void NA6PMC::forceJpsiDecays()
{
  int Jpsi_decay_mode[6][3] = {{0}};
  float Jpsi_BR[6] = {100.f, 0.f, 0.f, 0.f, 0.f, 0.f};
  Jpsi_decay_mode[0][0] = 13;  // mu-
  Jpsi_decay_mode[0][1] = -13; // mu+
  TVirtualMC::GetMC()->SetDecayMode(443, Jpsi_BR, Jpsi_decay_mode);

}

long NA6PMC::canGenerateMaxEvents() const
{
  return mGenerator ? mGenerator->canGenerateMaxEvents() : 0;
}

void NA6PMC::GeneratePrimaries()
{
  // Implement primary particle generation here
  mGenerator->generate();
  if (mVerbosity > 0) {
    LOGP(info, "{} primaries, {} tracks, generator: {}", mStack->GetNprimary(), mStack->GetNtrack(), mGenerator->getName());
    if (mVerbosity > 2) {
      mStack->Print();
    }
  }
}

void NA6PMC::BeginEvent()
{
  if (mVerbosity > 0) {
    LOGP(info, "Beginning event {}", mEvCount);
  }
  mStack->clear();
  clearHits();
  mMCTracks.clear();
  if (mMuonReferenceAudit) {
    mMuonReferenceAudit->beginEvent();
  }
  if (mGenerator) {
    mGenerator->clear();
  }
}

void NA6PMC::FinishEvent()
{
  if (mVerbosity > 0) {
    std::string rep = fmt::format(" Tracked {} particles ({} primaries) | Hits:", mStack->GetNtrack(), mStack->GetNprimary());
    for (int im = 0; im < NA6PDetector::instance().getNActiveModules(); im++) {
      const auto* det = NA6PDetector::instance().getActiveModule(im);
      rep += fmt::format(" {}:{}", det->getName(), det->getNHits());
    }
    LOGP(info, "Finishing event {} | {}", mEvCount, rep);
  }
  selectTracksToSave();
  writeKine();
  NA6PDetector::instance().writeHits(mRemap);
  mEvCount++;
}

void NA6PMC::BeginPrimary()
{
  if (mVerbosity > 1) {
    LOGP(info, "Beginning primary particle tracking");
  }
}

void NA6PMC::FinishPrimary()
{
  if (mVerbosity > 1) {
    LOGP(info, "Finishing primary particle tracking");
  }
}

void NA6PMC::PreTrack()
{
  // Optional geometry-independent virtual plane monitor. Never changes
  // the transported particle or introduces artificial detector material.
  if (mMuonReferenceAudit) {
    auto* mc = TVirtualMC::GetMC();
    mMuonReferenceAudit->preTrack(
        static_cast<int>(mEvCount + mEventOffset),
        mStack->GetCurrentTrackNumber(), mc->TrackPid());
  }
}

void NA6PMC::Stepping()
{
  // The reference-plane audit observes ALL muon steps, even outside
  // sensitive volumes. It interpolates crossed Z planes within steps.
  if (mMuonReferenceAudit) {
    auto* mc = TVirtualMC::GetMC();
    mMuonReferenceAudit->step(
        static_cast<int>(mEvCount + mEventOffset),
        mStack->GetCurrentTrackNumber(), mc->TrackPid());
  }
  if (!TVirtualMC::GetMC()->TrackCharge()) {
    return;
  }
  int volCopy;
  TVirtualMC::GetMC()->CurrentVolID(volCopy);
  if (NA6PModule::isSensor(volCopy)) {
    NA6PDetector::instance().getModule(NA6PModule::volID2ModuleID(volCopy))->stepManager(volCopy);
  } else {
    // non-sensitive volumes treatment
  }
  const auto& cuts = NA6PGenCutParam::Instance();
  if (cuts.transportMaxAbsZ > 0 || cuts.transportMaxR) {
    float x, y, z;
    TVirtualMC::GetMC()->TrackPosition(x, y, z);
    if ((cuts.transportMaxAbsZ > 0 && std::abs(z) > cuts.transportMaxAbsZ) ||
        (cuts.transportMaxR > 0 && x * x + y * y > cuts.transportMaxR * cuts.transportMaxR)) {
      TVirtualMC::GetMC()->StopTrack();
    }
  }
}

bool NA6PMC::setupGenerator(const std::string& s)
{
  if (s.empty()) {
    return false;
  }
  TInterpreter::EErrorCode errCode = TInterpreter::EErrorCode::kNoError;
  NA6PGenerator* gen = nullptr;
  auto parseGen = [&]() {
    const std::regex pattern(R"(^(\S+)((?:\.C|.C\+|\.C\+\+|\.C\+g|\.C\+\+g|\.cxx\+|\.cxx\+\+|\.cxx\+g|\.cxx\+\+g))(?:\(([^)]*)\))?$)");
    std::array<std::string, 3> res{};
    std::smatch match;
    if (std::regex_match(s, match, pattern)) {
      res[0] = match[1];
      res[1] = match[2];
      if (match[3].matched) {
        res[2] = match[3];
      }
    }
    return res;
  };
  auto genData = parseGen();
  if (genData[0].empty() || genData[1].empty()) {
    LOGP(fatal, "Failed to parse generator input string {}", s);
  }
  std::string cmd = genData[0] + genData[1];
  LOGP(info, "Loading generator macro {}", cmd);
  TInterpreter::Instance()->LoadMacro(cmd.c_str(), &errCode);
  if (errCode != TInterpreter::EErrorCode::kNoError) {
    LOGP(fatal, "Failed to load generator macro {}, error TInterpreter::EErrorCode={}", cmd.c_str(), int(errCode));
  }
  auto stemStart = genData[0].find_last_of('/');
  std::string invStr = fmt::format("{}({})", stemStart >= genData[0].length() ? genData[0] : genData[0].substr(stemStart + 1), genData[2]);
  LOGP(info, "Invoking generator initialization as {}", invStr);
  gen = dynamic_cast<NA6PGenerator*>((TObject*)TInterpreter::Instance()->ProcessLineSynch(invStr.c_str(), &errCode));
  if (errCode != TInterpreter::EErrorCode::kNoError || !gen) {
    LOGP(fatal, "Execution of generator function {} from macro {} leads to error TInterpreter::EErrorCode={} and returned pointer {}",
         invStr, s, int(errCode), (void*)gen);
  }
  gen->setStack(mStack.get());
  if (mUserVertexMethod) {
    gen->setUserVertexMethod(mUserVertexMethod.get());
  }
  if (mRandomSeed) {
    gen->setRandomSeed(mRandomSeed);
  }
  gen->init();
  mGenerator.reset(gen);
  LOGP(info, "Initialized generator {} via macro {}", mGenerator->getName(), s);
  return mGenerator != nullptr;
}

void NA6PMC::setupUserHooks(const std::string& s)
{
  LOGP(info, "Setting up user hooks macro {}", s);
  TInterpreter::EErrorCode errCode = TInterpreter::EErrorCode::kNoError;
  std::string sfname{}, sCmd{};
  if (Str::endsWith(s, ".C+") || Str::endsWith(s, ".C++") || Str::endsWith(s, ".C+g") || Str::endsWith(s, ".C++g")) {
    static const std::regex pattern(R"((\+\+g|\+g|\+\+|\+)$)");
    sfname = std::regex_replace(s, pattern, "");
    sCmd = s;
  } else if (Str::endsWith(s, ".C")) {
    sfname = s;
    sCmd = s + "+"; // always compile
  } else {
    LOGP(fatal, "User hooks macro {} extention differs from .C, .C+, .C++, .C+g or .C++g", s);
  }
  TInterpreter::Instance()->LoadMacro(sCmd.c_str(), &errCode);
  if (errCode != TInterpreter::EErrorCode::kNoError) {
    LOGP(fatal, "Failed to compile and load user hooks macro {}, error TInterpreter::EErrorCode={}", sCmd, int(errCode));
  }
  std::filesystem::path mpth(sfname);
  mUserHookName = mpth.stem().string();
  mUserHooksMethod = std::make_unique<TMethodCall>();
  mUserHooksMethod->InitWithPrototype(mUserHookName.c_str(), "int, bool");
}

void NA6PMC::setupUserVertex(const std::string& s)
{
  LOGP(info, "Setting up user vertex generator macro {}", s);
  TInterpreter::EErrorCode errCode = TInterpreter::EErrorCode::kNoError;
  std::string sfname{}, sCmd{};
  if (Str::endsWith(s, ".C+") || Str::endsWith(s, ".C++") || Str::endsWith(s, ".C+g") || Str::endsWith(s, ".C++g")) {
    static const std::regex pattern(R"((\+\+g|\+g|\+\+|\+)$)");
    sfname = std::regex_replace(s, pattern, "");
    sCmd = s;
  } else if (Str::endsWith(s, ".C")) {
    sfname = s;
    sCmd = s + "+"; // always compile
  } else {
    LOGP(fatal, "User vertex generator macro {} extention differs from .C, .C+, .C++, .C+g or .C++g", s);
  }
  TInterpreter::Instance()->LoadMacro(sCmd.c_str(), &errCode);
  if (errCode != TInterpreter::EErrorCode::kNoError) {
    LOGP(fatal, "Failed to compile and load user vertex generator macro {}, error TInterpreter::EErrorCode={}", sCmd, int(errCode));
  }
  std::filesystem::path mpth(sfname);
  mUserVertexMacroName = mpth.stem().string();
  mUserVertexMethod = std::make_unique<TMethodCall>();
  mUserVertexMethod->InitWithPrototype(mUserVertexMacroName.c_str(), "float&, float&, float&");
}

void NA6PMC::setRandomSeed(Long64_t r)
{
  if (r < 0) {
    auto t = std::chrono::system_clock::now().time_since_epoch();
    auto timeSeed = std::chrono::duration_cast<std::chrono::nanoseconds>(t).count();
    auto pidSeed = static_cast<ULong64_t>(getpid());
    mRandomSeed = std::hash<std::string>{}(fmt::format("{}:{}", timeSeed, pidSeed));
  } else {
    mRandomSeed = r;
  }
  if (mRandomSeed) {
    gRandom->SetSeed(mRandomSeed);
    LOGP(info, "Set gRandom seed to {}", mRandomSeed);
  }
}

void NA6PMC::init()
{
  if (mVerbosity < 2) {
    MiscUtils::silenceStdOut("for BuildPhysics()");
  }

  TVirtualMC::GetMC()->Init();

  TVirtualMC::GetMC()->BuildPhysics();
  if (mVerbosity < 2) {
    MiscUtils::reviveStdOut();
  }

  TVirtualMC::GetMC()->SetStack(mStack.get());

  NA6PDetector::instance().setVerbosity(mVerbosity);
  mStack->setVerbosity(mVerbosity);
  if (mGenerator) {
    mGenerator->setVerbosity(mVerbosity);
  }

  auto dir = na6p::utils::Str::rectifyDirectory(na6p::conf::KeyValParam::Instance().output_dir);
  if (dir.empty()) {
    dir = na6p::utils::Str::rectifyDirectory(".");
  }
  createKineOutput(dir);
  NA6PDetector::instance().createHitsOutput(dir);

  const char* refAuditFlag = std::getenv("NA6P_MWPC_STORE_REFERENCE_PLANES");
  if (refAuditFlag && refAuditFlag[0] != '0') {
    mMuonReferenceAudit = std::make_unique<NA6PMuonReferenceAudit>(dir);
    LOGP(info, "MWPC layout-independent reference plane muon states enabled: "
               "MuonReferencePlanes.root (10-cm grid, 365/430/495 cm)");
  }
}

void NA6PMC::clearHits()
{
  if (mVerbosity > 0) {
    LOGP(info, "clearHits");
  }
  for (int i = 0; i < NA6PDetector::instance().getNModules(); i++) {
    auto det = NA6PDetector::instance().getModule(i);
    det->clearHits();
  }
}

void NA6PMC::createKineOutput(const std::string& outDir)
{
  auto nm = fmt::format("{}MCKine.root", outDir);
  mKineFile = TFile::Open(nm.c_str(), "recreate");
  mKineTree = new TTree("mckine", "MC kinematics");
  static NA6PMCEventHeader* mchPtr = mStack->getEventHeader();
  mKineTree->Branch("header", &mchPtr);
  mKineTree->Branch("tracks", &mMCTracksPtr);
  LOGP(info, "Will store MC Kinematics in {}", nm);
}

void NA6PMC::closeKineOutput()
{
  if (mKineTree && mKineFile) {
    mKineFile->cd();
    mKineTree->Write();
    delete mKineTree;
    mKineTree = 0;
    mKineFile->Close();
    delete mKineFile;
    mKineFile = 0;
  }
}

void NA6PMC::writeKine()
{
  if (mKineTree) {
    mKineTree->Fill();
  }
}

void NA6PMC::selectTracksToSave()
{
  // store track if primary or if it has hits (then store all ancestors)

  auto resetMD = [](TParticle& p) {
    p.SetFirstDaughter(-1);
    p.SetLastDaughter(-1);
    p.SetFirstMother(-1);
    p.SetLastMother(-1);
  };

  auto resetMothers = [](TParticle& p) {
    p.SetFirstMother(-1);
    p.SetLastMother(-1);
  };
  auto resetDaughters = [](TParticle& p) {
    p.SetFirstDaughter(-1);
    p.SetLastDaughter(-1);
  };
  callUserHook(UserHook::SelectParticles, true); // call at entry
  int ntrIni = mStack->GetNtrack(), nPrimIni = mStack->GetNprimary();
  mRemap.clear();
  mSavID.clear();
  mDtList.clear();
  mRemap.resize(mStack->GetNtrack(), -1);
  int nkeep = 0;
  auto mcHeader = mStack->getEventHeader();
  auto& genHeaders = mcHeader->getGenHeaders();
  int nGenHeaders = genHeaders.size();
  std::vector<int> family;
  if (nGenHeaders) {
    for (int ig = 0; ig < (int)genHeaders.size(); ig++) {
      auto& genh = genHeaders[ig];
      if (genh.getPrimariesOffset() != nkeep || genh.getPrimariesOffset() + genh.getNPrimaries() > nPrimIni) {
        LOGP(fatal, "Processing primaries of generator {} expected in the [{}:{}) range, but we are at offset {} of total {} primaries, check generators",
             ig, genh.getPrimariesOffset(), genh.getPrimariesOffset() + genh.getNPrimaries(), nkeep, nPrimIni);
      }
      genh.setPrimariesOffset(nkeep);
      for (int i = 0; i < genh.getNPrimaries(); i++) {
        mSavID.push_back(nkeep);
        mRemap[nkeep++] = ig; // register generator id at the moment
      }
      genh.setSecondariesOffset(-1);
    }
  } else {
    for (int i = 0; i < nPrimIni; i++) { // generator primaries are stored as is
      mSavID.push_back(i);
      mRemap[i] = 0;
    }
  }
  for (int i = nPrimIni; i < ntrIni; i++) {
    auto* part = mStack->GetParticle(i);
    bool isFromHFDecay = false;
    int idMoth = part->GetFirstMother();
    int mothPdg = -1;
    while (idMoth >= 0) {
      auto* currMoth = mStack->GetParticle(idMoth);
      int absPdg = std::abs(currMoth->GetPdgCode());
      if (absPdg == 11 || absPdg == 22 || absPdg == 211 || absPdg == 130 || absPdg == 321 || absPdg == 2212 || absPdg == 2112) {
        // stop if a "stable" particle is found in the ancestors
        isFromHFDecay = false;
        break;
      }
      if ((absPdg > 400 && absPdg < 600) || (absPdg > 4000 && absPdg < 6000) || absPdg == 100443 || absPdg == 20443) {
        isFromHFDecay = true;
        mothPdg = absPdg;
        break;
      }
      idMoth = currMoth->GetFirstMother();
    }
    if (NA6PModule::testActiveIDOrKeepBits(*part) || isFromHFDecay) { // has hits
      if (mRemap[i] >= 0) {                    // was already accounted
        continue;
      }
      family.clear();
      family.push_back(i); // we store temporarily IDs of particles to store
      int mothID, daughID = i;
      LOGP(debug, "will save contributor {}(pdg{}), mother {}({}), Z={} | Dt: {} {}", i, part->GetPdgCode(), part->GetFirstMother(), part->GetFirstMother() >= 0 ? mRemap[part->GetFirstMother()] : -1, part->Vz(), part->GetFirstDaughter(), part->GetLastDaughter());
      while ((mothID = part->GetFirstMother()) >= 0) {
        // store tmp daughter index
        part = mStack->GetParticle(mothID);
        int entry0 = mDtList.size(), idd = part->GetFirstDaughter();
        auto& ref = mDtList.emplace_back(daughID, idd); // if idd>=0 then some daughters are already added and this was the entry of the last accounted daughter
        part->SetFirstDaughter(entry0);
        part->SetLastDaughter(idd < 0 ? 1 : part->GetLastDaughter() + 1); // here we temporarily accumulate number of daughters
        if (mRemap[mothID] >= 0) {
          LOGP(debug, "stop while on mothID={}", mothID);
          break;
        }
        family.push_back(mothID);
        LOGP(debug, "will save mother {}(pdg{}), its mother {}({}) Z={} | Dt: {} {}", mothID, part->GetPdgCode(), part->GetFirstMother(), part->GetFirstMother() >= 0 ? mRemap[part->GetFirstMother()] : -1, part->Vz(), part->GetFirstDaughter(), part->GetLastDaughter());
        daughID = mothID;
      }
      // we reached mother or already accounted secondary
      if (mothID < 0 || mRemap[mothID] < 0) {
        LOGP(fatal, "While processing a secondary {} encountered unregistered secondary {}", i, mStack->getParticleIndex(part));
      }
      for (int j = (int)family.size(); j--;) {
        mSavID.push_back(family[j]);
        mRemap[family[j]] = mRemap[mothID]; // generator ID
      }
    }
  }
  // update headers
  mcHeader->setNTracks(mSavID.size());
  mcHeader->setNPrimaries(nPrimIni);
  mcHeader->setNColl(mStack->getNColl());
  mcHeader->setNPart(mStack->getNPart());
  mcHeader->setImpPar(mStack->getImpPar());
  mcHeader->setEventID(mEventOffset + mEvCount);
  // register secondaries offsets and store primaries
  nkeep = 0;
  for (int i : mSavID) {
    if (i >= nPrimIni) {
      auto& gh = genHeaders[mRemap[i]];
      if (gh.getSecondariesOffset() == -1) {
        gh.setSecondariesOffset(nkeep);
        gh.setNSecondaries(0);
      }
      gh.incNSecondaries();
      mRemap[i] = -1; // reset remapping
    } else {
      auto* part = mStack->GetParticle(i); // store primary
      mRemap[i] = mMCTracks.size();
      mMCTracks.push_back(*part);
      //      resetMD(mMCTracks.back());
    }
    nkeep++;
    /*
    // check M->D
    auto* part = mStack->GetParticle(i);
    int de = part->GetFirstDaughter(), nd = part->GetLastDaughter();
    if (de>=0) {
      LOGP(info, "Daughters of {} pdg:{} Vz:{:.2f} Nd={} | de0 entry: {}", i, part->GetPdgCode(), part->Vz(), nd, de);
      int idn = 0;
      while (de>=0) {
  auto ref = mDtList[de];
  part = mStack->GetParticle(ref.first);
  LOGP(info, " -> Dt {}/{} pdg:{} Vz:{:.2f} | {} / {}", idn++,nd, part->GetPdgCode(), part->Vz(), ref.first, ref.second);
  de = ref.second;
      }
    }
    */
  }
  for (int i : mSavID) {
    int mothID = i;
    if (mRemap[i] < 0) {                   // need to store
      auto* part = mStack->GetParticle(i); // store primary
      mothID = mMCTracks.size();
      mMCTracks.push_back(*part);
      resetMothers(mMCTracks.back());
      mRemap[i] = mothID;
    } else {
      mothID = mRemap[i]; // particle was already added, just store its direct daughters and fix duaghter indices
    }
    auto& mother = mMCTracks[mothID];
    int entry = mother.GetFirstDaughter(); // collect IDs of all direct daughters for saving
    family.clear();
    while (entry >= 0) {
      const auto& dtentry = mDtList[entry];
      family.push_back(dtentry.first); // by construction, we accumulate from last to 1st daughter
      entry = dtentry.second;
    }
    if (family.size()) { // there are daughters to register
      mother.SetFirstDaughter(mMCTracks.size());
      mother.SetLastDaughter(mMCTracks.size() + family.size() - 1);
    } else {
      resetDaughters(mother);
    }
    for (unsigned int idd = family.size(); idd--;) { // MUST go in reverse order
      if (mRemap[family[idd]] >= 0) {
        LOGP(fatal, "Daughter {} was already stored as {} !!!", family[idd], mRemap[family[idd]]);
      }
      auto* daughter = mStack->GetParticle(family[idd]);
      mRemap[family[idd]] = mMCTracks.size();
      mMCTracks.push_back(*daughter);
      mMCTracks.back().SetFirstMother(mothID);
      mMCTracks.back().SetLastMother(family[idd]); // orig particle ID
    }
  }
  callUserHook(UserHook::SelectParticles, false); // call at exit
  LOGP(info, "Will save {} tracks out of {}", mMCTracks.size(), mRemap.size());
}

void NA6PMC::addSpecialParticles()
{
  // This is a copy/paste from O2 https://github.com/AliceO2Group/AliceO2/Steer/src/O2MCApplication.cxx
  //
  // Add particles needed for ALICE (not present in Geant3 or Geant4)
  // Code ported 1-1 from AliRoot
  //

  LOG(info) << "Adding custom particles to VMC";

  // Hypertriton
  TVirtualMC::GetMC()->DefineParticle(1010010030, "HyperTriton", kPTHadron, 2.991134, 1.0, 2.632e-10, "Ion", 0.0, 0, 1, 0, 0, 0, 0, 0, 3, kFALSE);
  // Anti-Hypertriton
  TVirtualMC::GetMC()->DefineParticle(-1010010030, "AntiHyperTriton", kPTHadron, 2.991134, 1.0, 2.632e-10, "Ion", 0.0, 0, 1, 0, 0, 0, 0, 0, 3, kFALSE);

  // Hyper hydrogen 4 ground state
  TVirtualMC::GetMC()->DefineParticle(1010010040, "Hyperhydrog4", kPTHadron, 3.922434, 1.0, 2.08e-10, "Ion", 0.0, 0, 1, 0, 0, 0, 0, 0, 4, kFALSE);
  // Anti-Hyper hydrogen 4 ground state
  TVirtualMC::GetMC()->DefineParticle(-1010010040, "AntiHyperhydrog4", kPTHadron, 3.922434, 1.0, 2.08e-10, "Ion", 0.0, 0, 1, 0, 0, 0, 0, 0, 4, kFALSE);

  // Hyper helium 4 ground state
  TVirtualMC::GetMC()->DefineParticle(1010020040, "Hyperhelium4", kPTHadron, 3.921728, 2.0, 2.50e-10, "Ion", 0.0, 0, 1, 0, 0, 0, 0, 0, 4, kFALSE);
  // Anti-Hyper helium 4 ground state
  TVirtualMC::GetMC()->DefineParticle(-1010020040, "AntiHyperhelium4", kPTHadron, 3.921728, 2.0, 2.50e-10, "Ion", 0.0, 0, 1, 0, 0, 0, 0, 0, 4, kFALSE);

  // Lithium 4 ground state
  TVirtualMC::GetMC()->DefineParticle(1000030040, "Lithium4", kPTHadron, 3.7513, 3.0, 9.1e-23, "Ion", 0.003, 0, 1, 0, 0, 0, 0, 0, 4, kFALSE);
  // Anti-Lithium 4 ground state
  TVirtualMC::GetMC()->DefineParticle(-1000030040, "AntiLithium4", kPTHadron, 3.7513, 3.0, 9.1e-23, "Ion", 0.003, 0, 1, 0, 0, 0, 0, 0, 4, kFALSE);

  // Hyper helium 5
  TVirtualMC::GetMC()->DefineParticle(1010020050, "Hyperhelium5", kPTHadron, 4.839961, 2.0, 2.74e-10, "Ion", 0.0, 0, 1, 0, 0, 0, 0, 0, 5, kFALSE);
  // Anti-Hyper helium 5
  TVirtualMC::GetMC()->DefineParticle(-1010020050, "AntiHyperhelium5", kPTHadron, 4.839961, 2.0, 2.74e-10, "Ion", 0.0, 0, 1, 0, 0, 0, 0, 0, 5, kFALSE);

  // Double Hyper hydrogen 4
  TVirtualMC::GetMC()->DefineParticle(1020010040, "DoubleHyperhydrogen4", kPTHadron, 4.106, 1.0, 2.632e-10, "Ion", 0.0, 0, 1, 0, 0, 0, 0, 0, 4, kFALSE);
  // Double Anti-Hyper hydrogen 4
  TVirtualMC::GetMC()->DefineParticle(-1020010040, "DoubleAntiHyperhydrogen4", kPTHadron, 4.106, 1.0, 2.632e-10, "Ion", 0.0, 0, 1, 0, 0, 0, 0, 0, 4, kFALSE);

  // 4Xi(-)H
  TVirtualMC::GetMC()->DefineParticle(1120010040, "4XiH", kPTHadron, 4.128, 1.0, 1.639e-10, "Ion", 0.0, 0, 1, 0, 0, 0, 0, 0, 4, kFALSE);
  // Anti-4Xi(-)H
  TVirtualMC::GetMC()->DefineParticle(-1120010040, "Anti4XiH", kPTHadron, 4.128, 1.0, 1.639e-10, "Ion", 0.0, 0, 1, 0, 0, 0, 0, 0, 4, kFALSE);
  // 4Xi(-)He
  TVirtualMC::GetMC()->DefineParticle(1120020040, "4XiHe", kPTHadron, 4.128, 1.0, 1.639e-10, "Ion", 0.0, 0, 1, 0, 0, 0, 0, 0, 4, kFALSE);
  // Anti-4Xi(-)He
  TVirtualMC::GetMC()->DefineParticle(-1120020040, "Anti4XiHe", kPTHadron, 4.128, 1.0, 1.639e-10, "Ion", 0.0, 0, 1, 0, 0, 0, 0, 0, 4, kFALSE);

  // Hyper helium 4 sigma
  TVirtualMC::GetMC()->DefineParticle(1110020040, "Hyperhelium4sigma", kPTHadron, 3.995, 2.0, 8.018e-11, "Ion", 0.0, 0, 1, 0, 0, 0, 0, 0, 4, kFALSE);
  // Anti-Hyper helium 4 sigma
  TVirtualMC::GetMC()->DefineParticle(-1110020040, "AntiHyperhelium4sigma", kPTHadron, 3.995, 2.0, 8.018e-11, "Ion", 0.0, 0, 1, 0, 0, 0, 0, 0, 4, kFALSE);

  // Lambda-Neutron
  TVirtualMC::GetMC()->DefineParticle(1010000020, "LambdaNeutron", kPTNeutron, 2.054, 0.0, 2.632e-10, "Hadron", 0.0, 0, 1, 0, 0, 0, 0, 0, 2, kFALSE);

  // Anti-Lambda-Neutron
  TVirtualMC::GetMC()->DefineParticle(-1010000020, "AntiLambdaNeutron", kPTNeutron, 2.054, 0.0, 2.632e-10, "Hadron", 0.0, 0, 1, 0, 0, 0, 0, 0, 2, kFALSE);

  // H-Dibaryon
  TVirtualMC::GetMC()->DefineParticle(1020000020, "Hdibaryon", kPTNeutron, 2.23, 0.0, 2.632e-10, "Hadron", 0.0, 0, 1, 0, 0, 0, 0, 0, 2, kFALSE);

  // Anti-H-Dibaryon
  TVirtualMC::GetMC()->DefineParticle(-1020000020, "AntiHdibaryon", kPTNeutron, 2.23, 0.0, 2.632e-10, "Hadron", 0.0, 0, 1, 0, 0, 0, 0, 0, 2, kFALSE);

  // Xi-Proton
  TVirtualMC::GetMC()->DefineParticle(1020010020, "Xi0Proton", kPTHadron, 2.248, 1.0, 1.333e-10, "Ion", 0.0, 0, 1, 0, 0, 0, 0, 0, 2, kFALSE);

  // Anti-Xi-Proton
  TVirtualMC::GetMC()->DefineParticle(-1020010020, "AntiXi0Proton", kPTHadron, 2.248, 1.0, 1.333e-10, "Ion", 0.0, 0, 1, 0, 0, 0, 0, 0, 2, kFALSE);

  // Lambda-Neutron-Neutron
  TVirtualMC::GetMC()->DefineParticle(1010000030, "LambdaNeutronNeutron", kPTNeutron, 2.99, 0.0, 2.632e-10, "Hadron", 0.0, 0, 1, 0, 0, 0, 0, 0, 3, kFALSE);

  // Anti-Lambda-Neutron-Neutron
  TVirtualMC::GetMC()->DefineParticle(-1010000030, "AntiLambdaNeutronNeutron", kPTNeutron, 2.99, 0.0, 2.632e-10, "Hadron", 0.0, 0, 1, 0, 0, 0, 0, 0, 3, kFALSE);

  // Omega-Proton
  TVirtualMC::GetMC()->DefineParticle(1030000020, "OmegaProton", kPTNeutron, 2.592, 0.0, 2.632e-10, "Hadron", 0.0, 2, 1, 0, 0, 0, 0, 0, 2, kFALSE);

  // Anti-Omega-Proton
  TVirtualMC::GetMC()->DefineParticle(-1030000020, "AntiOmegaProton", kPTNeutron, 2.592, 0.0, 2.632e-10, "Hadron", 0.0, 2, 1, 0, 0, 0, 0, 0, 2, kFALSE);

  // Omega-Neutron
  TVirtualMC::GetMC()->DefineParticle(1030010020, "OmegaNeutron", kPTHadron, 2.472, 1.0, 2.190e-22, "Hadron", 0.0, 2, 1, 0, 0, 0, 0, 0, 2, kFALSE);

  // Anti-Omega-Neutron
  TVirtualMC::GetMC()->DefineParticle(-1030010020, "AntiOmegaNeutron", kPTHadron, 2.472, 1.0, 2.190e-22, "Hadron", 0.0, 2, 1, 0, 0, 0, 0, 0, 2, kFALSE);

  // Omega-Omega
  TVirtualMC::GetMC()->DefineParticle(1060020020, "OmegaOmega", kPTHadron, 3.229, 2.0, 2.632e-10, "Hadron", 0.0, 0, 1, 0, 0, 0, 0, 0, 2, kFALSE);

  // Anti-Omega-Omega
  TVirtualMC::GetMC()->DefineParticle(-1060020020, "AntiOmegaOmega", kPTHadron, 3.229, 2.0, 2.632e-10, "Hadron", 0.0, 0, 1, 0, 0, 0, 0, 0, 2, kFALSE);

  // Lambda(1405)-Proton
  TVirtualMC::GetMC()->DefineParticle(1010010021, "Lambda1405Proton", kPTHadron, 2.295, 1.0, 1.316e-23, "Hadron", 0.0, 0, 1, 0, 0, 0, 0, 0, 2, kFALSE);

  // Anti-Lambda(1405)-Proton
  TVirtualMC::GetMC()->DefineParticle(-1010010021, "AntiLambda1405Proton", kPTHadron, 2.295, 1.0, 1.316e-23, "Hadron", 0.0, 0, 1, 0, 0, 0, 0, 0, 2, kFALSE);

  // Lambda(1405)-Lambda(1405)
  TVirtualMC::GetMC()->DefineParticle(1020000021, "Lambda1405Lambda1405", kPTNeutron, 2.693, 0.0, 1.316e-23, "Hadron", 0.0, 0, 1, 0, 0, 0, 0, 0, 2, kFALSE);

  // Anti-Lambda(1405)-Lambda(1405)
  TVirtualMC::GetMC()->DefineParticle(-1020000021, "AntiLambda1405Lambda1405", kPTNeutron, 2.693, 0.0, 1.316e-23, "Hadron", 0.0, 0, 1, 0, 0, 0, 0, 0, 2, kFALSE);

  // c-deuteron
  TVirtualMC::GetMC()->DefineParticle(2010010020, "CDeuteron", kPTHadron, 3.226, 1.0, 2.0e-13, "Ion", 0.0, 0, 1, 0, 0, 0, 0, 0, 3, kFALSE);
  // Anti-c-deuteron
  TVirtualMC::GetMC()->DefineParticle(-2010010020, "AntiCDeuteron", kPTHadron, 3.226, 1.0, 2.0e-13, "Ion", 0.0, 0, 1, 0, 0, 0, 0, 0, 3, kFALSE);

  // c-triton
  TVirtualMC::GetMC()->DefineParticle(2010010030, "CTriton", kPTHadron, 4.162, 1.0, 2.0e-13, "Ion", 0.0, 0, 1, 0, 0, 0, 0, 0, 2, kFALSE);
  // Anti-c-Hypertriton
  TVirtualMC::GetMC()->DefineParticle(-2010010030, "AntiCTriton", kPTHadron, 4.162, 1.0, 2.0e-13, "Ion", 0.0, 0, 1, 0, 0, 0, 0, 0, 2, kFALSE);

  // Resonances not in Generators
  //  f0(980) assume 70 MeV as width (PDG: 40 to 100 MeV)
  TVirtualMC::GetMC()->DefineParticle(9010221, "f0_980", kPTNeutron, 0.98, 0.0, 9.403e-24, "Hadron", 7e-2, 0, 1, 1, 0, 0, 1, 0, 0, kTRUE);

  // f2(1270) (PDG: width = 185 MeV)
  TVirtualMC::GetMC()->DefineParticle(225, "f2_1270", kPTNeutron, 1.275, 0.0, 3.558e-24, "Hadron", 0.185, 4, 1, 1, 0, 0, 1, 0, 0, kTRUE);

  // f1(1285) (PDG: width = 24.20 MeV) Spin/Parity might not be correct
  TVirtualMC::GetMC()->DefineParticle(20223, "f1_1285", kPTNeutron, 1.28210, 0.0, 1e-24, "Hadron", 0.02420, 3, 1, 0, 0, 0, 0, 0, 1, kTRUE);
  // f1(1420) (PDG: width = 54 MeV) Spin/Parity might not be correct
  TVirtualMC::GetMC()->DefineParticle(20333, "f1_1420", kPTNeutron, 1.42640, 0.0, 1e-24, "Hadron", 0.05490, 3, 1, 0, 0, 0, 0, 0, 1, kTRUE);

  // Glueball hunting family
  // Their life times are not known, so we set them to 1e-24
  // f0(1370) (PDG: width = 200-500 MeV) Spin/Parity might not be correct
  TVirtualMC::GetMC()->DefineParticle(10221, "f0_1370", kPTNeutron, 1.37, 0.0, 1e-24, "Hadron", 0.2, 1, 1, 1, 0, 0, 1, 0, 0, kTRUE);
  // a2(1320) (PDG: width = 107.8 MeV) (Spin/Parity might not be correct)
  TVirtualMC::GetMC()->DefineParticle(115, "a2_1320", kPTNeutron, 1.3182, 0.0, 1e-24, "Hadron", 0.1078, 1, 1, 1, 1, 0, 1, 0, 0, kTRUE);
  // f0(1500) (PDG: width = 112 MeV) Spin/Parity might not be correct
  TVirtualMC::GetMC()->DefineParticle(9030221, "f0_1500", kPTNeutron, 1.506, 0.0, 1e-24, "Hadron", 0.112, 0, 1, 1, 0, 0, 1, 0, 0, kTRUE);
  // f0(1710) (PDG: width = 139 MeV) Spin/Parity might not be correct
  TVirtualMC::GetMC()->DefineParticle(10331, "f0_1710", kPTNeutron, 1.71, 0.0, 1e-24, "Hadron", 0.139, 1, 1, 1, 0, 0, 1, 0, 0, kTRUE);
  // f2(1525) (PDG: width = 73 MeV) Spin/Parity might not be correct
  TVirtualMC::GetMC()->DefineParticle(335, "f2_1525", kPTNeutron, 1.525, 0.0, 1e-24, "Hadron", 0.073, 5, 1, 1, 0, 0, 1, 0, 0, kTRUE);

  // Xi_0(1820)
  TVirtualMC::GetMC()->DefineParticle(123324, "Xi_0_1820", kPTNeutron, 1.8234, 0.0, 2.742550e-23, "Hadron", 0.024, 3, -1, 0, 1, 1, 0, 0, 1, kTRUE);
  TVirtualMC::GetMC()->DefineParticle(-123324, "Xi_0_Bar_1820", kPTNeutron, 1.8234, 0.0, 2.742550e-23, "Hadron", 0.024, 3, -1, 0, 1, -1, 0, 0, -1, kTRUE);

  int xi_0_1820_mode[6][3] = {{0}};
  float xi_0_1820_ratio[6] = {100.f, 0.f, 0.f, 0.f, 0.f, 0.f};
  xi_0_1820_mode[0][0] = 3122; // Lambda
  xi_0_1820_mode[0][1] = 310;  // K0s
  TVirtualMC::GetMC()->SetDecayMode(123324, xi_0_1820_ratio, xi_0_1820_mode);
  xi_0_1820_mode[0][0] = -3122; // Lambda-bar
  TVirtualMC::GetMC()->SetDecayMode(-123324, xi_0_1820_ratio, xi_0_1820_mode);

  // Xi-+(1820)
  TVirtualMC::GetMC()->DefineParticle(123314, "Xi_Minus_1820", kPTHadron, 1.8234, -1.0, 2.742550e-23, "Hadron", 0.024, 3, -1, 0, 1, -1, 0, 0, 1, kTRUE);
  TVirtualMC::GetMC()->DefineParticle(-123314, "Xi_Plus_1820", kPTHadron, 1.8234, 1.0, 2.742550e-23, "Hadron", 0.024, 3, -1, 0, 1, 1, 0, 0, -1, kTRUE);

  int xi_charged_1820_mode[6][3] = {{0}};
  float xi_charged_1820_ratio[6] = {100.f, 0.f, 0.f, 0.f, 0.f, 0.f};
  xi_charged_1820_mode[0][0] = 3122; // Lambda
  xi_charged_1820_mode[0][1] = -321; // K-
  TVirtualMC::GetMC()->SetDecayMode(123314, xi_charged_1820_ratio, xi_charged_1820_mode);
  xi_charged_1820_mode[0][0] = -3122; // Lambda-bar
  xi_charged_1820_mode[0][1] = 321;   // K+
  TVirtualMC::GetMC()->SetDecayMode(-123314, xi_charged_1820_ratio, xi_charged_1820_mode);

  // Ps - hidden strange (s-sbar) pentaquarks
  TVirtualMC::GetMC()->DefineParticle(9322134, "Ps_2100", kPTHadron, 2.1, 1.0, 1.6455e-23, "Hadron", 4.e-2, 3, -1, 0, 0, 0, 0, 0, 1, kTRUE);
  TVirtualMC::GetMC()->DefineParticle(-9322134, "AntiPs_2100", kPTHadron, 2.1, -1.0, 1.6455e-23, "Hadron", 4.e-2, 3, -1, 0, 0, 0, 0, 0, -1, kTRUE);
  TVirtualMC::GetMC()->DefineParticle(9322136, "Ps_2500", kPTHadron, 2.5, 1.0, 1.6455e-23, "Hadron", 4.e-2, 5, 1, 0, 0, 0, 0, 0, 1, kTRUE);
  TVirtualMC::GetMC()->DefineParticle(-9322136, "AntiPs_2500", kPTHadron, 2.5, -1.0, 1.6455e-23, "Hadron", 4.e-2, 5, 1, 0, 0, 0, 0, 0, -1, kTRUE);

  Int_t psmode[6][3] = {0};
  Float_t psratio[6] = {0.f};
  psratio[0] = 100.;

  psmode[0][0] = 333;  // phi
  psmode[0][1] = 2212; // proton
  TVirtualMC::GetMC()->SetDecayMode(9322134, psratio, psmode);
  TVirtualMC::GetMC()->SetDecayMode(9322136, psratio, psmode);

  psmode[0][1] = -2212; // anti-proton
  TVirtualMC::GetMC()->SetDecayMode(-9322134, psratio, psmode);
  TVirtualMC::GetMC()->SetDecayMode(-9322136, psratio, psmode);

  // Omega(2012)
  for (int j = 1; j < 6; j++) {
    psmode[j][0] = psmode[j][1] = 0;
    psratio[j] = 0.;
  }

  TVirtualMC::GetMC()->DefineParticle(3335, "Omega2012", kPTHadron, 2.012, -1.0, 1.0285e-22, "Hadron", 0.0064, 3, -1, 0, 0, 0, 0, 0, 1, kTRUE);
  psmode[0][0] = 3312; // Xi-
  psmode[0][1] = 310;  // K0S
  psratio[0] = 100.;
  TVirtualMC::GetMC()->SetDecayMode(3335, psratio, psmode);

  TVirtualMC::GetMC()->DefineParticle(-3335, "AntiOmega2012", kPTHadron, 2.012, 1.0, 1.0285e-22, "Hadron", 0.0064, 3, 1, 0, 0, 0, 0, 0, -1, kTRUE);
  psmode[0][0] = -3312; // anti-Xi+
  psmode[0][1] = 310;   // K0S
  psratio[0] = 100.;
  TVirtualMC::GetMC()->SetDecayMode(-3335, psratio, psmode);

  // d*(2380) - dibaryon resonance
  TVirtualMC::GetMC()->DefineParticle(900010020, "d*_2380", kPTHadron, 2.38, 1.0, 0.94e-23, "Ion", 0.07, 6, 1, 0, 0, 0, 0, 0, 2, kTRUE);
  TVirtualMC::GetMC()->DefineParticle(-900010020, "d*_2380_bar", kPTHadron, 2.38, -1.0, 0.94e-23, "Ion", 0.07, 6, 1, 0, 0, 0, 0, 0, -2, kTRUE);

  Int_t dstmode[6][3] = {0};
  Float_t dstratio[6] = {0.f};
  dstratio[0] = 100; // For now we implement only the mode of interest
  // d* -> d pi+ pi-
  dstmode[0][0] = 1000010020; // deuteron
  dstmode[0][1] = -211;       // negative pion
  dstmode[0][2] = 211;        // positive pion
  TVirtualMC::GetMC()->SetDecayMode(900010020, dstratio, dstmode);

  dstmode[0][0] = -1000010020; // anti-deuteron
  TVirtualMC::GetMC()->SetDecayMode(-900010020, dstratio, dstmode);

  // Heavy vector mesons
  // D*+
  TVirtualMC::GetMC()->DefineParticle(413, "D*+", kPTHadron, 2.0103, 1.0, 0.0, "Hadron", 0.0, 1, -1, 0, 0, 0, 0, 0, 0, kTRUE);
  // D*-
  TVirtualMC::GetMC()->DefineParticle(-413, "D*-", kPTHadron, 2.0103, -1.0, 0.0, "Hadron", 0.0, 1, -1, 0, 0, 0, 0, 0, 0, kTRUE);
  // D*0
  TVirtualMC::GetMC()->DefineParticle(423, "D*0", kPTHadron, 2.0007, 0.0, 0.0, "Hadron", 0.0, 1, -1, 0, 0, 0, 0, 0, 0, kTRUE);
  // D*0bar
  TVirtualMC::GetMC()->DefineParticle(-423, "D*0bar", kPTHadron, 2.0007, 0.0, 0.0, "Hadron", 0.0, 1, -1, 0, 0, 0, 0, 0, 0, kTRUE);
  // D*_s+
  TVirtualMC::GetMC()->DefineParticle(433, "D*_s+", kPTHadron, 2.1123, 1.0, 0.0, "Hadron", 0.0, 1, -1, 0, 0, 0, 0, 0, 0, kTRUE);
  // D*_s-
  TVirtualMC::GetMC()->DefineParticle(-433, "D*_s-", kPTHadron, 2.1123, -1.0, 0.0, "Hadron", 0.0, 1, -1, 0, 0, 0, 0, 0, 0, kTRUE);
  // B*0
  TVirtualMC::GetMC()->DefineParticle(513, "B*0", kPTHadron, 5.3251, 0.0, 0.0, "Hadron", 0.0, 1, -1, 0, 0, 0, 0, 0, 0, kTRUE);
  // B*0bar
  TVirtualMC::GetMC()->DefineParticle(-513, "B*0bar", kPTHadron, 5.3251, 0.0, 0.0, "Hadron", 0.0, 1, -1, 0, 0, 0, 0, 0, 0, kTRUE);
  // B*+
  TVirtualMC::GetMC()->DefineParticle(523, "B*+", kPTHadron, 5.3251, 1.0, 0.0, "Hadron", 0.0, 1, -1, 0, 0, 0, 0, 0, 0, kTRUE);
  // B*-
  TVirtualMC::GetMC()->DefineParticle(-523, "B*-", kPTHadron, 5.3251, -1.0, 0.0, "Hadron", 0.0, 1, -1, 0, 0, 0, 0, 0, 0, kTRUE);
  // B*_s0
  TVirtualMC::GetMC()->DefineParticle(533, "B*_s0", kPTHadron, 5.4128, 0.0, 0.0, "Hadron", 0.0, 1, -1, 0, 0, 0, 0, 0, 0, kTRUE);
  // B*_s0bar
  TVirtualMC::GetMC()->DefineParticle(-533, "B*_s0bar", kPTHadron, 5.4128, 0.0, 0.0, "Hadron", 0.0, 1, -1, 0, 0, 0, 0, 0, 0, kTRUE);
  // B*_c+
  TVirtualMC::GetMC()->DefineParticle(543, "B*_c+", kPTHadron, 6.6020, 1.0, 0.0, "Hadron", 0.0, 1, -1, 0, 0, 0, 0, 0, 0, kTRUE);
  // B*_c-
  TVirtualMC::GetMC()->DefineParticle(-543, "B*_c-", kPTHadron, 6.6020, -1.0, 0.0, "Hadron", 0.0, 1, -1, 0, 0, 0, 0, 0, 0, kTRUE);

  // Charm pentaquarks
  // Theta_c: isospin singlet with J=1/2+ (see https://arxiv.org/abs/hep-ph/0409121)
  TVirtualMC::GetMC()->DefineParticle(9422111, "Anti-Theta_c_3100", kPTHadron, 3.099, 0., 6.9e-21, "Hadron", 83.e-6, 1, 1, 0, 0, 0, 0, 0, -1, kTRUE);
  TVirtualMC::GetMC()->DefineParticle(-9422111, "Theta_c_3100", kPTHadron, 3.099, 0., 6.9e-21, "Hadron", 83.e-6, 1, 1, 0, 0, 0, 0, 0, 1, kTRUE);

  for (int j = 1; j < 6; j++) {
    psmode[j][0] = psmode[j][1] = 0;
    psratio[j] = 0.;
  }
  psmode[0][0] = 413;   // D*+
  psmode[0][1] = -2212; // anti-p
  psratio[0] = 100.;
  TVirtualMC::GetMC()->SetDecayMode(9422111, psratio, psmode);
  psmode[0][0] = -413; // D*-
  psmode[0][1] = 2212; // p
  TVirtualMC::GetMC()->SetDecayMode(-9422111, psratio, psmode);

  // Define the 2- and 3-body phase space decay for the Hyper-Triton
  Int_t mode[6][3];
  Float_t bratio[6];

  for (Int_t kz = 0; kz < 6; kz++) {
    bratio[kz] = 0.;
    mode[kz][0] = 0;
    mode[kz][1] = 0;
    mode[kz][2] = 0;
  }
  bratio[0] = 50.;
  mode[0][0] = 1000020030; // Helium3
  mode[0][1] = -211;       // negative pion

  bratio[1] = 50.;
  mode[1][0] = 1000010020; // deuteron
  mode[1][1] = 2212;       // proton
  mode[1][2] = -211;       // negative pion

  TVirtualMC::GetMC()->SetDecayMode(1010010030, bratio, mode);

  // Define the 2- and 3-body phase space decay for the Anti-Hyper-Triton
  Int_t amode[6][3];
  Float_t abratio[6];

  for (Int_t kz = 0; kz < 6; kz++) {
    abratio[kz] = 0.;
    amode[kz][0] = 0;
    amode[kz][1] = 0;
    amode[kz][2] = 0;
  }
  abratio[0] = 50.;
  amode[0][0] = -1000020030; // anti- Helium3
  amode[0][1] = 211;         // positive pion
  abratio[1] = 50.;
  amode[1][0] = -1000010020; // anti-deuteron
  amode[1][1] = -2212;       // anti-proton
  amode[1][2] = 211;         // positive pion

  TVirtualMC::GetMC()->SetDecayMode(-1010010030, abratio, amode);

  ////// ----------Hypernuclei with Mass=4 ----------- //////////

  // Define the 2- and 3-body phase space decay for the Hyper Hydrogen 4

  Int_t mode3[6][3];
  Float_t bratio3[6];

  for (Int_t kz = 0; kz < 6; kz++) {
    bratio3[kz] = 0.;
    mode3[kz][0] = 0;
    mode3[kz][1] = 0;
    mode3[kz][2] = 0;
  }
  bratio3[0] = 50.;
  mode3[0][0] = 1000020040; // Helium4
  mode3[0][1] = -211;       // negative pion

  bratio3[1] = 50.;
  mode3[1][0] = 1000010030; // tritium
  mode3[1][1] = 2212;       // proton
  mode3[1][2] = -211;       // negative pion

  TVirtualMC::GetMC()->SetDecayMode(1010010040, bratio3, mode3);

  // Define the 2- and 3-body phase space decay for the Hyper Hydrogen 4
  Int_t amode3[6][3];
  Float_t abratio3[6];

  for (Int_t kz = 0; kz < 6; kz++) {
    abratio3[kz] = 0.;
    amode3[kz][0] = 0;
    amode3[kz][1] = 0;
    amode3[kz][2] = 0;
  }
  abratio3[0] = 50.;
  amode3[0][0] = -1000020040; // anti- Helium4
  amode3[0][1] = 211;         // positive pion
  abratio3[1] = 50.;
  amode3[1][0] = -1000010030; // anti-tritium
  amode3[1][1] = -2212;       // anti-proton
  amode3[1][2] = 211;         // positive pion

  TVirtualMC::GetMC()->SetDecayMode(-1010010040, abratio3, amode3);

  // Define the 3-body phase space decay for the Hyper Helium 4
  Int_t mode4[6][3];
  Float_t bratio4[6];

  for (Int_t kz = 0; kz < 6; kz++) {
    bratio4[kz] = 0.;
    mode4[kz][0] = 0;
    mode4[kz][1] = 0;
    mode4[kz][2] = 0;
  }
  bratio4[0] = 50.;
  mode4[0][0] = 1000020030; // Helium3
  mode4[0][1] = -211;       // negative pion
  mode4[0][2] = 2212;       // proton

  bratio4[1] = 50.;
  mode4[1][0] = 1000030040; // lithium-4
  mode4[1][1] = -211;       // negative pion

  TVirtualMC::GetMC()->SetDecayMode(1010020040, bratio4, mode4);

  // Define the 2-body phase space decay for the Anti-Hyper Helium 4
  Int_t amode4[6][3];
  Float_t abratio4[6];

  for (Int_t kz = 0; kz < 6; kz++) {
    abratio4[kz] = 0.;
    amode4[kz][0] = 0;
    amode4[kz][1] = 0;
    amode4[kz][2] = 0;
  }
  abratio4[0] = 50.;
  amode4[0][0] = -1000020030; // anti-Helium 3
  amode4[0][1] = 211;         // positive pion
  amode4[0][2] = -2212;       // anti proton

  abratio4[1] = 50.;
  amode4[1][0] = -1000030040; // antilithium-4
  amode4[1][1] = 211;         // positive pion

  TVirtualMC::GetMC()->SetDecayMode(-1010020040, abratio4, amode4);

  // Define the 2-body phase space decay for the Lithium 4
  Int_t model4[6][3];
  Float_t bratiol4[6];

  for (Int_t kz = 0; kz < 6; kz++) {
    bratiol4[kz] = 0.;
    model4[kz][0] = 0;
    model4[kz][1] = 0;
    model4[kz][2] = 0;
  }
  bratiol4[0] = 100.;
  model4[0][0] = 1000020030; // Helium3
  model4[0][1] = 2212;       // proton

  TVirtualMC::GetMC()->SetDecayMode(1000030040, bratiol4, model4);

  // Define the 2-body phase space decay for the Anti-Lithium 4
  Int_t amodel4[6][3];
  Float_t abratiol4[6];

  for (Int_t kz = 0; kz < 6; kz++) {
    abratiol4[kz] = 0.;
    amodel4[kz][0] = 0;
    amodel4[kz][1] = 0;
    amodel4[kz][2] = 0;
  }
  abratiol4[0] = 100.;
  amodel4[0][0] = -1000020030; // Anti-Helium3
  amodel4[0][1] = -2212;       // Anti-proton

  TVirtualMC::GetMC()->SetDecayMode(-1000030040, abratiol4, amodel4);

  // Define the 3-body phase space decay for the Hyper Helium 5
  Int_t mode41[6][3];
  Float_t bratio41[6];

  for (Int_t kz = 0; kz < 6; kz++) {
    bratio41[kz] = 0.;
    mode41[kz][0] = 0;
    mode41[kz][1] = 0;
    mode41[kz][2] = 0;
  }
  bratio41[0] = 50.;
  mode41[0][0] = 1000020040; // Helium4
  mode41[0][1] = -211;       // negative pion
  mode41[0][2] = 2212;       // proton
  bratio41[1] = 50.;
  mode41[1][0] = 1000020030; // Helium3
  mode41[1][1] = -211;       // negative pion
  mode41[1][2] = 1000010020; // Deuteron

  TVirtualMC::GetMC()->SetDecayMode(1010020050, bratio41, mode41);

  // Define the 2-body phase space decay for the Anti-Hyper Helium 5
  Int_t amode41[6][3];
  Float_t abratio41[6];

  for (Int_t kz = 0; kz < 6; kz++) {
    abratio41[kz] = 0.;
    amode41[kz][0] = 0;
    amode41[kz][1] = 0;
    amode41[kz][2] = 0;
  }
  abratio41[0] = 50.;
  amode41[0][0] = -1000020040; // anti-Helium 4
  amode41[0][1] = 211;         // positive pion
  amode41[0][2] = -2212;       // anti proton
  abratio41[1] = 50.;
  amode41[1][0] = -1000020030; // anti-Helium 3
  amode41[1][1] = 211;         // positive pion
  amode41[1][2] = -1000010020; // anti deuteron

  TVirtualMC::GetMC()->SetDecayMode(-1010020050, abratio41, amode41);

  // Define the 3-body phase space decay for the Double Hyper Hydrogen 4
  Int_t mode42[6][3];
  Float_t bratio42[6];

  for (Int_t kz = 0; kz < 6; kz++) {
    bratio42[kz] = 0.;
    mode42[kz][0] = 0;
    mode42[kz][1] = 0;
    mode42[kz][2] = 0;
  }
  bratio42[0] = 50.;
  mode42[0][0] = 1010020040; // Hyper-Helium4
  mode42[0][1] = -211;       // negative pion

  bratio42[1] = 50.;
  mode42[1][0] = 1010010030; // Hypertriton
  mode42[1][1] = 2212;       // proton
  mode42[1][2] = -211;       // negative pion

  TVirtualMC::GetMC()->SetDecayMode(1020010040, bratio42, mode42);

  // Define the 2-body phase space decay for the Anti Double Hyper Hydrogen 4
  Int_t amode42[6][3];
  Float_t abratio42[6];

  for (Int_t kz = 0; kz < 6; kz++) {
    abratio42[kz] = 0.;
    amode42[kz][0] = 0;
    amode42[kz][1] = 0;
    amode42[kz][2] = 0;
  }
  abratio42[0] = 50.;
  amode42[0][0] = -1010020040; // anti-Hyper-Helium 4
  amode42[0][1] = 211;         // positive pion

  abratio42[1] = 50.;
  amode42[1][0] = -1010010030; // anti-Hypertriton
  amode42[1][1] = -2212;       // antiproton
  amode42[1][2] = 211;         // positive pion

  TVirtualMC::GetMC()->SetDecayMode(-1020010040, abratio42, amode42);

  // Define the decay for the 4Xi(-)He
  Int_t mode4XiHe[6][3];
  Float_t bratio4XiHe[6];

  for (Int_t kz = 0; kz < 6; kz++) {
    bratio4XiHe[kz] = 0.;
    mode4XiHe[kz][0] = 0;
    mode4XiHe[kz][1] = 0;
    mode4XiHe[kz][2] = 0;
  }
  bratio4XiHe[0] = 33.;
  mode4XiHe[0][0] = 1010020040; // HyperHelium4
  mode4XiHe[0][1] = -211;       // negative pion

  bratio4XiHe[1] = 33.;
  mode4XiHe[1][0] = 3122;       // lambda
  mode4XiHe[1][1] = 1000020030; // helium-3
  mode4XiHe[1][2] = -211;       // negative pion

  bratio4XiHe[2] = 33.;
  mode4XiHe[2][0] = 1000030040; // lithium-4
  mode4XiHe[2][1] = -211;       // negative pion
  mode4XiHe[2][2] = -211;       // negative pion

  TVirtualMC::GetMC()->SetDecayMode(1120020040, bratio4XiHe, mode4XiHe);

  // Define the decay for the Anti-4Xi(-)He
  Int_t amode4XiHe[6][3];
  Float_t abratio4XiHe[6];

  for (Int_t kz = 0; kz < 6; kz++) {
    abratio4XiHe[kz] = 0.;
    amode4XiHe[kz][0] = 0;
    amode4XiHe[kz][1] = 0;
    amode4XiHe[kz][2] = 0;
  }
  abratio4XiHe[0] = 33.;
  amode4XiHe[0][0] = -1010020040; // antiHyperHelium-4
  amode4XiHe[0][1] = 211;         // positive pion

  abratio4XiHe[1] = 33.;
  amode4XiHe[1][0] = -3122;       // antilambda
  amode4XiHe[1][1] = -1000020030; // antihelium-3
  amode4XiHe[1][2] = 211;         // positive pion

  abratio4XiHe[2] = 33.;
  amode4XiHe[2][0] = -1000030040; // antilithium-4
  amode4XiHe[2][1] = 211;         // positive pion
  amode4XiHe[2][2] = 211;         // positive pion

  TVirtualMC::GetMC()->SetDecayMode(-1120020040, abratio4XiHe, amode4XiHe);

  // Define the decay for the 4Xi(-)H
  Int_t mode4XiH[6][3];
  Float_t bratio4XiH[6];

  for (Int_t kz = 0; kz < 6; kz++) {
    bratio4XiH[kz] = 0.;
    mode4XiH[kz][0] = 0;
    mode4XiH[kz][1] = 0;
    mode4XiH[kz][2] = 0;
  }
  bratio4XiH[0] = 33.;
  mode4XiH[0][0] = 1010010040; // HyperHydrogen4
  mode4XiH[0][1] = -211;       // negative pion

  bratio4XiH[1] = 33.;
  mode4XiH[1][0] = 3122;       // lambda
  mode4XiH[1][1] = 1000010030; // triton
  mode4XiH[1][2] = -211;       // negative pion

  bratio4XiH[2] = 33.;
  mode4XiH[2][0] = 1000020040; // alpha
  mode4XiH[2][1] = -211;       // negative pion
  mode4XiH[2][2] = -211;       // negative pion

  TVirtualMC::GetMC()->SetDecayMode(1120010040, bratio4XiH, mode4XiH);

  // Define the decay for the Anti-4Xi(-)H
  Int_t amode4XiH[6][3];
  Float_t abratio4XiH[6];

  for (Int_t kz = 0; kz < 6; kz++) {
    abratio4XiH[kz] = 0.;
    amode4XiH[kz][0] = 0;
    amode4XiH[kz][1] = 0;
    amode4XiH[kz][2] = 0;
  }
  abratio4XiH[0] = 33.;
  amode4XiH[0][0] = -1010010040; // antiHyperHydrogen-4
  amode4XiH[0][1] = 211;         // positive pion

  abratio4XiH[1] = 33.;
  amode4XiH[1][0] = -3122;       // antilambda
  amode4XiH[1][1] = -1000010030; // antitriton
  amode4XiH[1][2] = 211;         // positive pion

  abratio4XiH[2] = 33.;
  amode4XiH[2][0] = -1000020040; // antialpha
  amode4XiH[2][1] = 211;         // positive pion
  amode4XiH[2][2] = 211;         // positive pion

  TVirtualMC::GetMC()->SetDecayMode(-1120010040, abratio4XiH, amode4XiH);

  // Define the 2- and 3-body phase space decay for the Hyper Helium 4 sigma
  Int_t mode4s[6][3];
  Float_t bratio4s[6];

  for (Int_t kz = 0; kz < 6; kz++) {
    bratio4s[kz] = 0.;
    mode4s[kz][0] = 0;
    mode4s[kz][1] = 0;
    mode4s[kz][2] = 0;
  }
  bratio4s[0] = 20.;
  mode4s[0][0] = 1000020040; // Helium4
  mode4s[0][1] = 111;        // pion0
  bratio4s[1] = 40.;
  mode4s[1][0] = 1000010030; // tritium
  mode4s[1][2] = 2212;       // proton
  mode4s[1][1] = 111;        // pion0
  bratio4s[2] = 40.;
  mode4s[2][0] = 1000010030; // tritium
  mode4s[2][2] = 2212;       // pion+
  mode4s[2][1] = 2112;       // neutron

  TVirtualMC::GetMC()->SetDecayMode(1110020040, bratio4s, mode4s);

  // Define the 2- and 3-body phase space decay for the Anti Hyper Helium 4 sigma
  Int_t amode4s[6][3];
  Float_t abratio4s[6];

  for (Int_t kz = 0; kz < 6; kz++) {
    abratio4s[kz] = 0.;
    amode4s[kz][0] = 0;
    amode4s[kz][1] = 0;
    amode4s[kz][2] = 0;
  }
  abratio4s[0] = 50.;
  amode4s[0][0] = -1000020040; // anti-Helium4
  amode4s[0][1] = 111;         // pion0
  abratio4s[1] = 50.;
  amode4s[1][0] = -1000010030; // anti-tritium
  amode4s[1][2] = -2212;       // anti-proton
  amode4s[1][1] = 111;         // pion0
  abratio4s[2] = 50.;
  amode4s[2][0] = -1000010030; // anti-tritium
  amode4s[2][2] = -211;        // pion-
  amode4s[2][1] = -2112;       // anti-neutron

  TVirtualMC::GetMC()->SetDecayMode(-1110020040, abratio4s, amode4s);

  // Define the 2-body phase space decay for the Lambda-neutron boundstate
  Int_t mode1[6][3];
  Float_t bratio1[6];

  for (Int_t kz = 0; kz < 6; kz++) {
    bratio1[kz] = 0.;
    mode1[kz][0] = 0;
    mode1[kz][1] = 0;
    mode1[kz][2] = 0;
  }
  bratio1[0] = 100.;
  mode1[0][0] = 1000010020; // deuteron
  mode1[0][1] = -211;       // negative pion

  TVirtualMC::GetMC()->SetDecayMode(1010000020, bratio1, mode1);

  // Define the 2-body phase space decay for the Anti-Lambda-neutron boundstate
  Int_t amode1[6][3];
  Float_t abratio1[6];

  for (Int_t kz = 0; kz < 6; kz++) {
    abratio1[kz] = 0.;
    amode1[kz][0] = 0;
    amode1[kz][1] = 0;
    amode1[kz][2] = 0;
  }
  abratio1[0] = 100.;
  amode1[0][0] = -1000010020; // anti-deuteron
  amode1[0][1] = 211;         // positive pion

  TVirtualMC::GetMC()->SetDecayMode(-1010000020, abratio1, amode1);

  // Define the 2-body phase space decay for the H-Dibaryon
  Int_t mode2[6][3];
  Float_t bratio2[6];

  for (Int_t kz = 0; kz < 6; kz++) {
    bratio2[kz] = 0.;
    mode2[kz][0] = 0;
    mode2[kz][1] = 0;
    mode2[kz][2] = 0;
  }
  bratio2[0] = 100.;
  mode2[0][0] = 3122; // Lambda
  mode2[0][1] = 2212; // proton
  mode2[0][2] = -211; // negative pion

  TVirtualMC::GetMC()->SetDecayMode(1020000020, bratio2, mode2);

  // Define the 2-body phase space decay for the Anti-H-Dibaryon
  Int_t amode2[6][3];
  Float_t abratio2[6];

  for (Int_t kz = 0; kz < 6; kz++) {
    abratio2[kz] = 0.;
    amode2[kz][0] = 0;
    amode2[kz][1] = 0;
    amode2[kz][2] = 0;
  }
  abratio2[0] = 100.;
  amode2[0][0] = -3122; // anti-deuteron
  amode2[0][1] = -2212; // anti-proton
  amode2[0][2] = 211;   // positive pion

  TVirtualMC::GetMC()->SetDecayMode(-1020000020, abratio2, amode2);

  // Define the 2-body phase space decay for the Xi0P
  Int_t mode5[6][3];
  Float_t bratio5[6];

  for (Int_t kz = 0; kz < 6; kz++) {
    bratio5[kz] = 0.;
    mode5[kz][0] = 0;
    mode5[kz][1] = 0;
    mode5[kz][2] = 0;
  }
  bratio5[0] = 100.;
  mode5[0][0] = 3122; // Lambda
  mode5[0][1] = 2212; // proton

  TVirtualMC::GetMC()->SetDecayMode(1020010020, bratio5, mode5);

  // Define the 2-body phase space decay for the Anti-Xi0P
  Int_t amode5[6][3];
  Float_t abratio5[6];

  for (Int_t kz = 0; kz < 6; kz++) {
    abratio5[kz] = 0.;
    amode5[kz][0] = 0;
    amode5[kz][1] = 0;
    amode5[kz][2] = 0;
  }
  abratio5[0] = 100.;
  amode5[0][0] = -3122; // anti-Lambda
  amode5[0][1] = -2212; // anti-proton

  TVirtualMC::GetMC()->SetDecayMode(-1020010020, abratio5, amode5);

  // Define the 2-body phase space decay for the Lambda-Neutron-Neutron
  Int_t mode6[6][3];
  Float_t bratio6[6];

  for (Int_t kz = 0; kz < 6; kz++) {
    bratio6[kz] = 0.;
    mode6[kz][0] = 0;
    mode6[kz][1] = 0;
    mode6[kz][2] = 0;
  }
  bratio6[0] = 100.;
  mode6[0][0] = 1000010030; // triton
  mode6[0][1] = -211;       // pion

  TVirtualMC::GetMC()->SetDecayMode(1010000030, bratio6, mode6);

  // Define the 2-body phase space decay for the Anti-Lambda-Neutron-Neutron
  Int_t amode6[6][3];
  Float_t abratio6[6];

  for (Int_t kz = 0; kz < 6; kz++) {
    abratio6[kz] = 0.;
    amode6[kz][0] = 0;
    amode6[kz][1] = 0;
    amode6[kz][2] = 0;
  }
  abratio6[0] = 100.;
  amode6[0][0] = -1000010030; // anti-triton
  amode6[0][1] = 211;         // pion

  TVirtualMC::GetMC()->SetDecayMode(-1010000030, abratio6, amode6);

  // Define the 3-body phase space decay for the Omega-Proton
  Int_t mode7[6][3];
  Float_t bratio7[6];

  for (Int_t kz = 0; kz < 6; kz++) {
    bratio7[kz] = 0.;
    mode7[kz][0] = 0;
    mode7[kz][1] = 0;
    mode7[kz][2] = 0;
  }
  bratio7[0] = 100.;
  mode7[0][0] = 3122; // Lambda
  mode7[0][1] = -321; // negative Kaon
  mode7[0][2] = 2212; // proton

  TVirtualMC::GetMC()->SetDecayMode(1030000020, bratio7, mode7);

  // Define the 3-body phase space decay for the Anti-Omega-Proton
  Int_t amode7[6][3];
  Float_t abratio7[6];

  for (Int_t kz = 0; kz < 6; kz++) {
    abratio7[kz] = 0.;
    amode7[kz][0] = 0;
    amode7[kz][1] = 0;
    amode7[kz][2] = 0;
  }
  abratio7[0] = 100.;
  amode7[0][0] = -3122; // anti-Lambda
  amode7[0][1] = 321;   // positive kaon
  amode7[0][2] = -2212; // anti-proton

  TVirtualMC::GetMC()->SetDecayMode(-1030000020, abratio7, amode7);

  // Define the 2-body phase space decay for the Omega-Neutron
  Int_t mode8[6][3];
  Float_t bratio8[6];

  for (Int_t kz = 0; kz < 6; kz++) {
    bratio8[kz] = 0.;
    mode8[kz][0] = 0;
    mode8[kz][1] = 0;
    mode8[kz][2] = 0;
  }
  bratio8[0] = 100.;
  mode8[0][0] = 3122; // Lambda
  mode8[0][1] = 3312; // negative Xi

  TVirtualMC::GetMC()->SetDecayMode(1030010020, bratio8, mode8);

  // Define the 2-body phase space decay for the Anti-Omega-Neutron
  Int_t amode8[6][3];
  Float_t abratio8[6];

  for (Int_t kz = 0; kz < 6; kz++) {
    abratio8[kz] = 0.;
    amode8[kz][0] = 0;
    amode8[kz][1] = 0;
    amode8[kz][2] = 0;
  }
  abratio8[0] = 100.;
  amode8[0][0] = -3122; // anti-Lambda
  amode8[0][1] = -3312; // positive Xi

  TVirtualMC::GetMC()->SetDecayMode(-1030010020, abratio8, amode8);

  // Define the 3-body phase space decay for the Omega-Omega
  Int_t mode9[6][3];
  Float_t bratio9[6];

  for (Int_t kz = 0; kz < 6; kz++) {
    bratio9[kz] = 0.;
    mode9[kz][0] = 0;
    mode9[kz][1] = 0;
    mode9[kz][2] = 0;
  }
  bratio9[0] = 100.;
  mode9[0][0] = 3334; // negative Omega
  mode9[0][1] = 3312; // negative Xi

  TVirtualMC::GetMC()->SetDecayMode(1060020020, bratio9, mode9);

  // Define the 3-body phase space decay for the Anti-Omega-Omega
  Int_t amode9[6][3];
  Float_t abratio9[6];

  for (Int_t kz = 0; kz < 6; kz++) {
    abratio9[kz] = 0.;
    amode9[kz][0] = 0;
    amode9[kz][1] = 0;
    amode9[kz][2] = 0;
  }
  abratio9[0] = 100.;
  amode9[0][0] = -3334; // positive Omega
  amode9[0][1] = -3312; // positive Xi

  TVirtualMC::GetMC()->SetDecayMode(-1060020020, abratio9, amode9);

  // Define the 2- and 3-body phase space decay for the Lambda(1405)-Proton
  Int_t mode10[6][3];
  Float_t bratio10[6];

  for (Int_t kz = 0; kz < 6; kz++) {
    bratio10[kz] = 0.;
    mode10[kz][0] = 0;
    mode10[kz][1] = 0;
    mode10[kz][2] = 0;
  }
  bratio10[0] = 50.;
  mode10[0][0] = 3122; // Lambda
  mode10[0][1] = 2212; // proton
  bratio10[1] = 50.;
  mode10[1][0] = 2212; // proton
  mode10[1][1] = -321; // negative kaon
  mode10[1][2] = 2212; // proton

  TVirtualMC::GetMC()->SetDecayMode(1010010021, bratio10, mode10);

  // Define the 2- and 3-body phase space decay for the Anti-Lambda(1405)-Proton
  Int_t amode10[6][3];
  Float_t abratio10[6];

  for (Int_t kz = 0; kz < 6; kz++) {
    abratio10[kz] = 0.;
    amode10[kz][0] = 0;
    amode10[kz][1] = 0;
    amode10[kz][2] = 0;
  }
  abratio10[0] = 50.;
  amode10[0][0] = -3122; // anti-Lambda
  amode10[0][1] = -2212; // anti-proton
  abratio10[1] = 50.;
  amode10[1][0] = -2212; // anti-proton
  amode10[1][1] = 321;   // positive kaon
  amode10[1][2] = -2212; // anti-proton

  TVirtualMC::GetMC()->SetDecayMode(-1010010021, abratio10, amode10);

  // Define the 3-body phase space decay for the Lambda(1405)-Lambda(1405)
  Int_t mode11[6][3];
  Float_t bratio11[6];

  for (Int_t kz = 0; kz < 6; kz++) {
    bratio11[kz] = 0.;
    mode11[kz][0] = 0;
    mode11[kz][1] = 0;
    mode11[kz][2] = 0;
  }
  bratio11[0] = 50.;
  mode11[0][0] = 3122; // Lambda
  mode11[0][1] = 3122; // Lambda
  bratio11[1] = 50.;
  mode11[1][0] = 3122; // Lambda
  mode11[1][1] = 2212; // proton
  mode11[1][2] = -211; // negative pion

  TVirtualMC::GetMC()->SetDecayMode(1020000021, bratio11, mode11);

  // Define the 3-body phase space decay for the Anti-Lambda(1405)-Lambda(1405)
  Int_t amode11[6][3];
  Float_t abratio11[6];

  for (Int_t kz = 0; kz < 6; kz++) {
    abratio11[kz] = 0.;
    amode11[kz][0] = 0;
    amode11[kz][1] = 0;
    amode11[kz][2] = 0;
  }
  abratio11[0] = 50.;
  amode11[0][0] = -3122; // anti-Lambda
  amode11[0][1] = -3122; // anti-Lambda
  abratio11[1] = 50.;
  amode11[1][0] = -3122; // anti-Lambda
  amode11[1][1] = -2212; // anti-proton
  amode11[1][2] = 211;   // positive pion

  TVirtualMC::GetMC()->SetDecayMode(-1020000021, abratio11, amode11);

  // Define the decays for the c-triton
  Int_t ctmode[6][3];
  Float_t ctbratio[6];

  for (Int_t kz = 0; kz < 6; kz++) {
    ctbratio[kz] = 0.;
    ctmode[kz][0] = 0;
    ctmode[kz][1] = 0;
    ctmode[kz][2] = 0;
  }
  ctbratio[0] = 50.;
  ctmode[0][0] = 1000020030; // Helium3
  ctmode[0][1] = 310;        // K0s

  ctbratio[1] = 50.;
  ctmode[1][0] = 1000020030; // Helium3
  ctmode[1][1] = -321;       // negative kaon
  ctmode[1][2] = 211;        // positive pion

  TVirtualMC::GetMC()->SetDecayMode(2010010030, ctbratio, ctmode);

  // Define the decays for the anti-c-triton
  Int_t actmode[6][3];
  Float_t actbratio[6];

  for (Int_t kz = 0; kz < 6; kz++) {
    actbratio[kz] = 0.;
    actmode[kz][0] = 0;
    actmode[kz][1] = 0;
    actmode[kz][2] = 0;
  }
  actbratio[0] = 50.;
  actmode[0][0] = -1000020030; // Helium3
  actmode[0][1] = 310;         // K0s

  actbratio[1] = 50.;
  actmode[1][0] = -1000020030; // Helium3
  actmode[1][1] = 321;         // negative kaon
  actmode[1][2] = -211;        // positive pion

  TVirtualMC::GetMC()->SetDecayMode(-2010010030, actbratio, actmode);

  // Define the decays for the c-deuteron
  Int_t cdmode[6][3];
  Float_t cdbratio[6];

  for (Int_t kz = 0; kz < 6; kz++) {
    cdbratio[kz] = 0.;
    cdmode[kz][0] = 0;
    cdmode[kz][1] = 0;
    cdmode[kz][2] = 0;
  }
  cdbratio[0] = 50.;
  cdmode[0][0] = 1000010020; // deuteron
  cdmode[0][1] = -321;       // negative kaon
  cdmode[0][2] = 211;        // positive pion

  cdbratio[1] = 50.;
  cdmode[1][0] = 1000010020; // deuteron
  cdmode[1][1] = 310;        // K0s

  TVirtualMC::GetMC()->SetDecayMode(2010010020, cdbratio, cdmode);

  // Define the decays for the anti-c-deuteron
  Int_t acdmode[6][3];
  Float_t acdbratio[6];

  for (Int_t kz = 0; kz < 6; kz++) {
    acdbratio[kz] = 0.;
    acdmode[kz][0] = 0;
    acdmode[kz][1] = 0;
    acdmode[kz][2] = 0;
  }
  acdbratio[0] = 50.;
  acdmode[0][0] = -1000010020; // deuteron
  acdmode[0][1] = 321;         // negative kaon
  acdmode[0][2] = -211;        // positive pion

  acdbratio[1] = 50.;
  acdmode[1][0] = -1000010020; // deuteron
  acdmode[1][1] = 310;         // K0s

  TVirtualMC::GetMC()->SetDecayMode(-2010010020, acdbratio, acdmode);

  ///////////////////////////////////////////////////////////////////

  // Define the 2-body phase space decay for the f0(980)
  //  Int_t mode[6][3];
  //  Float_t bratio[6];

  for (Int_t kz = 0; kz < 6; kz++) {
    bratio[kz] = 0.;
    mode[kz][0] = 0;
    mode[kz][1] = 0;
    mode[kz][2] = 0;
  }
  bratio[0] = 100.;
  mode[0][0] = 211;  // pion
  mode[0][1] = -211; // pion

  TVirtualMC::GetMC()->SetDecayMode(9010221, bratio, mode);

  // Define the 2-body phase space decay for the f2(1270)
  //  Int_t mode[6][3];
  //  Float_t bratio[6];

  for (Int_t kz = 0; kz < 6; kz++) {
    bratio[kz] = 0.;
    mode[kz][0] = 0;
    mode[kz][1] = 0;
    mode[kz][2] = 0;
  }
  bratio[0] = 100.;
  mode[0][0] = 211;  // pion
  mode[0][1] = -211; // pion

  TVirtualMC::GetMC()->SetDecayMode(225, bratio, mode);

  // Define the 2-body phase space decay for the resonances: f0(1500), f2(1525), f0(1710
  for (Int_t kz = 0; kz < 6; kz++) {
    bratio[kz] = 0.;
    mode[kz][0] = 0;
    mode[kz][1] = 0;
    mode[kz][2] = 0;
  }
  bratio[0] = 100.;
  mode[0][0] = 310; // K0s
  mode[0][1] = 310; // K0s

  TVirtualMC::GetMC()->SetDecayMode(9030221, bratio, mode); // f0(1500)
  TVirtualMC::GetMC()->SetDecayMode(335, bratio, mode);     // f2(1525)
  TVirtualMC::GetMC()->SetDecayMode(10331, bratio, mode);   // f0(1710)
  TVirtualMC::GetMC()->SetDecayMode(10221, bratio, mode);   // f0(1370)
  TVirtualMC::GetMC()->SetDecayMode(115, bratio, mode);     // a2(1320)

  // Define the 3-body phase space decay for the resonances: f1(1285), f1(1420)
  for (Int_t kz = 0; kz < 6; kz++) {
    bratio[kz] = 0.;
    mode[kz][0] = 0;
    mode[kz][1] = 0;
    mode[kz][2] = 0;
  }

  bratio2[0] = 50.;
  mode[0][0] = 310;  // K0s
  mode[0][1] = -321; // anti-K
  mode[0][2] = 211;  // pion+

  bratio2[1] = 50.;
  mode[1][0] = 310;  // K0s
  mode[1][1] = 321;  // K
  mode[1][2] = -211; // pion-

  TVirtualMC::GetMC()->SetDecayMode(20223, bratio2, mode); // f1(1285)
  TVirtualMC::GetMC()->SetDecayMode(20333, bratio2, mode); // f1(1420)

  // Lambda1520/Lambda1520bar

  TVirtualMC::GetMC()->DefineParticle(102134, "Lambda1520", kPTNeutron, 1.5195, 0.0, 4.22e-23, "Hadron", 0.0156, 3, -1, 0, 0, 0, 0, 0, 1, kTRUE);
  TVirtualMC::GetMC()->DefineParticle(-102134, "Lambda1520bar", kPTNeutron, 1.5195, 0.0, 4.22e-23, "Hadron", 0.0156, 3, -1, 0, 0, 0, 0, 0, -1, kTRUE);

  // Lambda1520 decay modes
  Int_t lmode[9][3];
  Float_t lbratio[9];
  for (Int_t kz = 0; kz < 9; kz++) {
    lbratio[kz] = 0.;
    lmode[kz][0] = 0;
    lmode[kz][1] = 0;
    lmode[kz][2] = 0;
  }

  // L(1520) -> p K-
  lbratio[0] = 0.229944;
  lmode[0][0] = 2212;
  lmode[0][1] = -321;

  // L(1520) -> n K0
  lbratio[1] = 0.229944;
  lmode[1][0] = 2112;
  lmode[1][1] = -311;

  // L(1520) -> Sigma+ pi-
  lbratio[2] = 0.143076;
  lmode[2][0] = 3222;
  lmode[2][1] = -211;

  // L(1520) -> Sigma0 pi0
  lbratio[3] = 0.143076;
  lmode[3][0] = 3212;
  lmode[3][1] = 111;

  // L(1520) -> Sigma- pi+
  lbratio[4] = 0.143076;
  lmode[4][0] = 3112;
  lmode[4][1] = 211;

  // L(1520) -> Sigma*- pi+
  lbratio[5] = 0.034066;
  lmode[5][0] = 3114;
  lmode[5][1] = 211;

  // L(1520) -> Sigma*0 pi0
  lbratio[6] = 0.034066;
  lmode[6][0] = 3214;
  lmode[6][1] = 111;

  // L(1520) -> Sigma*+ pi-
  lbratio[7] = 0.034066;
  lmode[7][0] = 3224;
  lmode[7][1] = -211;

  // L(1520) -> Lambda gamma
  lbratio[8] = 0.008687;
  lmode[8][0] = 3122;
  lmode[8][1] = 22;

  TVirtualMC::GetMC()->SetDecayMode(102134, lbratio, lmode);

  // Lambda1520bar decay modes

  // L(1520)bar -> p- K+
  lbratio[0] = 0.229944;
  lmode[0][0] = -2212;
  lmode[0][1] = 321;

  // L(1520)bar -> nbar K0bar
  lbratio[1] = 0.229944;
  lmode[1][0] = -2112;
  lmode[1][1] = 311;

  // L(1520)bar -> Sigmabar- pi+
  lbratio[2] = 0.143076;
  lmode[2][0] = -3222;
  lmode[2][1] = 211;

  // L(1520)bar -> Sigma0bar pi0
  lbratio[3] = 0.143076;
  lmode[3][0] = -3212;
  lmode[3][1] = 111;

  // L(1520)bar -> Sigmabar+ pi-
  lbratio[4] = 0.143076;
  lmode[4][0] = -3112;
  lmode[4][1] = -211;

  // L(1520)bar -> anti-Sigma*- pi-
  lbratio[5] = 0.034066;
  lmode[5][0] = -3114;
  lmode[5][1] = -211;

  // L(1520)bar -> anti-Sigma*0 pi0
  lbratio[6] = 0.034066;
  lmode[6][0] = -3214;
  lmode[6][1] = 111;

  // L(1520)bar -> anti-Sigma*+ pi+
  lbratio[7] = 0.034066;
  lmode[7][0] = -3224;
  lmode[7][1] = 211;

  // L(1520)bar -> Anti-Lambda gamma
  lbratio[8] = 0.008687;
  lmode[8][0] = -3122;
  lmode[8][1] = 22;

  TVirtualMC::GetMC()->SetDecayMode(-102134, lbratio, lmode);

  // --------------------------------------------------------------------

  // Sexaquark (uuddss): compact, neutral and stable hypothetical bound state (arxiv.org/abs/1708.08951)
  TVirtualMC::GetMC()->DefineParticle(900000020, "Sexaquark", kPTUndefined, 2.0, 0.0, 4.35e+17, "Hadron", 0.0, 0, 1, 0, 0, 0, 0, 0, 2, kTRUE);
  TVirtualMC::GetMC()->DefineParticle(-900000020, "AntiSexaquark", kPTUndefined, 2.0, 0.0, 4.35e+17, "Hadron", 0.0, 0, 1, 0, 0, 0, 0, 0, -2, kTRUE);

  // --------------------------------------------------------------------
  {
    int pdgID = 23;
    if (!TVirtualMC::GetMC()->DefineParticle(pdgID, "Z0", kPTUndefined, 91.1876, 0., 2.638e-25, "Boson", 1., 0, 0, 0, 0, 0, 0, 0, 0, false, true)) {
      LOGP(warn, "Failed to add PDG={} to TVirtualMC::GetMC() particles", pdgID);
    }
    Float_t br[6] = {100.f, 0.f, 0.f, 0.f, 0.f, 0.f};
    Int_t mode[6][3] = {};
    mode[0][0] = 13;
    mode[0][1] = -13;
    if (!TVirtualMC::GetMC()->SetDecayMode(23, br, mode)) {
      LOGP(fatal, "Failed to set 100% dimuon branching ratio for PDG={}", pdgID);
    }
    LOGP(info, "Added Z0 as VirtualPhoton with PDG={} and 100% dimuon branching ratio", pdgID);
  }
}
