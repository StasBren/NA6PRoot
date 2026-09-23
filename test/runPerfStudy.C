#if !defined(__CLING__) || defined(__ROOTCLING__)
#include <TTree.h>
#include <TFile.h>
#include <TParticle.h>
#include <TStopwatch.h>
#include "NA6PVerTelCluster.h"
#include "NA6PMuonSpecCluster.h"
#include "NA6PMuonSpecHit.h"
#include "NA6PMCComposedLabel.h"
#include "NA6PMCTruthContainer.h"
#include "Propagator.h"
#include "NA6PMatch.h"
#include "NA6PTrack.h"
#include "NA6PTreeStreamRedirector.h"
#include "HistoManager.h"
#include <unordered_map>
#include <TDatabasePDG.h>
#include <TParticlePDG.h>
#include <TLorentzVector.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TStyle.h>
#include <TF1.h>
#include <TH1F.h>
#include <TH2F.h>
#endif

//----- CUTS
float cutRapidityMin = 1.8f;
float cutRapidityMax = 4.f;
float cutPtMin = 0.1f;

int minLrVT = 5;
int minLrMS = 6;
//-----

struct MCTrackInfo {
  enum TrStatus : uint8_t { New,
                            Discard,
                            Accept };
  TrStatus status = New;
  int8_t charge = 0;
  int8_t nvtLr = 0;      // number of layers with clusters
  int8_t nmsLr = 0;      // number of layers with clusters
  uint8_t nVTTracks = 0; // number of tracks bearing the label of this MC track
  uint8_t nMSTracks = 0; // number of tracks bearing the label of this MC track
  uint8_t nMatches = 0;  // number of tracks bearing the label of this MC track
  int bestVTTrack = -1;
  int bestMSTrack = -1;
  int bestMatch = -1;
  NA6PMCComposedLabel parent;
  std::vector<int> vtClID;
  std::vector<int> msClID;
  std::array<float, 3> momMS{}, posMS{};
  bool isReconstructableVT() const { return status == TrStatus::Accept && nvtLr >= minLrVT; }
  bool isReconstructableMS() const { return status == TrStatus::Accept && nmsLr >= minLrMS; }
};

const char* prefHisto[3] = {"VT", "MS", "Match"};
const char* prefAxis[2] = {"1/p", "rapidity"};
const char* prefAxisS[2] = {"pInv", "rapidity"};
const char* parName[5] = {"X", "Y", "Tx", "Ty", "q/Pxz"};
const float parMaxDiff[5] = {0.02, 0.02, 0.005, 0.005, 0.1};
const float maxMassDiff = 0.4;

const float hrangeParBias[5][2] = {{-0.015f, 0.015f}, {-0.015f, 0.015f}, {-1e-3f, 1e-3f}, {-1e-3f, 1e-3f}, {-1e-3f, 1e-3f}};
const float hrangeParSigm[5][2] = {{0.f, 0.006}, {0.0f, 0.006f}, {0.f, 0.003f}, {0.f, 0.003f}, {0., 0.02}};

const float hrangeEff[4][2] = {{0.8f, 1.f}, {0.2f, 1.f}, {0.2f, 1.f}, {0.2, 1.f}};
const float hrangePurity[4][2] = {{0.85f, 1.f}, {0.6f, 1.f}, {0.6f, 1.f}, {0.2f, 1.f}};

const std::array<int, 4> HMOffsets = {0, 1000, 2000, 3000};
const int nBinsPtI = 10;
const int nBinsRapidity = 10;
const int nBinsDiffKin = 50;
const int nBinsMassRes = 50;

na6p::HistoManager* hm = nullptr;

std::unordered_map<NA6PMCComposedLabel, MCTrackInfo> mcTruth;
std::unique_ptr<NA6PTreeStreamRedirector> outPerf;

TTree *mcTree = nullptr, *vtClTree = nullptr, *vtTrTree = nullptr, *msClTree = nullptr, *msTrTree = nullptr, *matchTree = nullptr;
TTree* msHitsTree = nullptr;
std::vector<TParticle>* mcArr = nullptr;
std::vector<NA6PVerTelCluster>* vtClVec = nullptr;
std::vector<NA6PTrack>* vtTrVec = nullptr;
std::vector<NA6PMuonSpecCluster>* msClVec = nullptr;
NA6PMCTruthContainer* vtCluMCLabels = nullptr;
NA6PMCTruthContainer* msCluMCLabels = nullptr;
std::vector<NA6PMCComposedLabel>* vtTrkLabels = nullptr;
std::vector<NA6PMCComposedLabel>* msTrkLabels = nullptr;
std::vector<NA6PMCComposedLabel>* matchTrkLabels = nullptr;
std::vector<NA6PTrack>* msTrVec = nullptr;
std::vector<NA6PMatch>* matchVec = nullptr;
std::vector<NA6PMuonSpecHit>* msHitsVec = nullptr;
Long64_t nEntries = 0;

void setupInputs(bool checkVTTracks, bool checkMSTracks, bool checkMatches, const char* dirSimu, const char* dirRec);
void loadEventData(int i);
void buildReconstuctableMCInfo();
void checkMSTracking();
void checkMatching();
void checkVTTracking();
void createOutput(bool checkVTTracks, bool checkMSTracks, bool checkMatches, bool checkDimu);
void bookHistos(bool checkVTTracks, bool checkMSTracks, bool checkMatches, bool checkDimu);
void finalizeHistos(na6p::HistoManager* hman);
void drawHistos(na6p::HistoManager* hman, int copy = 0, bool reproc = false);
void createPDF(const char* outName = "perfrep.pdf");

void runPerfStudy(int firstEv = 0,
                  int lastEv = 99999,
                  bool checkVTTracks = true,
                  bool checkMSTracks = true,
                  bool checkMatches = true,
                  bool checkDimu = true,
                  const char* dirSimu = ".",
                  const char* dirRec = ".")
{
  if (!Propagator::loadField() || !Propagator::loadGeometry(Form("%s/geometry.root", dirSimu))) {
    return;
  }
  setupInputs(checkVTTracks, checkMSTracks, checkMatches, dirSimu, dirRec);
  bookHistos(checkVTTracks, checkMSTracks, checkMatches, checkDimu);
  if (lastEv < 1) {
    lastEv = nEntries;
  } else if (lastEv > (int)nEntries) {
    lastEv = nEntries;
  }
  outPerf = std::make_unique<NA6PTreeStreamRedirector>("recperf.root");

  for (int ev = firstEv; ev < lastEv; ev++) {
    loadEventData(ev);
    buildReconstuctableMCInfo();
    if (checkVTTracks) {
      checkVTTracking();
    }
    if (checkMSTracks) {
      checkMSTracking();
    }
    if (checkMatches) {
      checkMatching();
    }

    createOutput(checkVTTracks, checkMSTracks, checkMatches, checkDimu);
  }

  outPerf->Close();
  outPerf.reset();
  finalizeHistos(hm);
  hm->write();
  drawHistos(hm);
  createPDF();
}

void bookHistos(bool checkVTTracks, bool checkMSTracks, bool checkMatches, bool checkDimu)
{
  hm = new na6p::HistoManager("", "recperf_histos.root");
  bool doHisto[3] = {checkVTTracks, checkMSTracks, checkMatches};
  float cutAxisMin[3][2] = {{0.f, cutRapidityMin}, {0.f, cutRapidityMin}, {0.f, cutRapidityMin}};   // min 1/p and rapidity for VT, MS and Matches
  float cutAxisMax[3][2] = {{1.f, cutRapidityMax}, {0.4f, cutRapidityMax}, {0.4f, cutRapidityMax}}; // max 1/p and rapidity for VT, MS and Matches
  int nbinsAxis[2] = {nBinsPtI, nBinsRapidity};
  for (int itp = 0; itp < 3; itp++) {
    if (!doHisto[itp]) {
      continue;
    }
    int offset = HMOffsets[itp];
    for (int iax = 0; iax < 2; iax++) {
      // efficiency and purity
      auto h1mc = new TH1F(Form("reconstructable_%s_%s", prefHisto[itp], prefAxis[iax]), Form("reconstructable %s-tracks vs %s;%s", prefHisto[itp], prefAxis[iax], prefAxis[iax]), nbinsAxis[iax], cutAxisMin[itp][iax], cutAxisMax[itp][iax]);
      auto h1rec = new TH1F(Form("reconstructed_%s_%s", prefHisto[itp], prefAxis[iax]), Form("reconstructed %s-tracks vs %s;%s", prefHisto[itp], prefAxis[iax], prefAxis[iax]), nbinsAxis[iax], cutAxisMin[itp][iax], cutAxisMax[itp][iax]);
      auto h1recCorr = new TH1F(Form("correct_%s_%s", prefHisto[itp], prefAxis[iax]), Form("reconstructed correctly %s-tracks vs %s;%s", prefHisto[itp], prefAxis[iax], prefAxis[iax]), nbinsAxis[iax], cutAxisMin[itp][iax], cutAxisMax[itp][iax]);
      auto h1recEff = new TH1F(Form("receff_%s_%s", prefHisto[itp], prefAxis[iax]), Form("rec.eff %s-tracks vs %s;%s", prefHisto[itp], prefAxis[iax], prefAxis[iax]), nbinsAxis[iax], cutAxisMin[itp][iax], cutAxisMax[itp][iax]);
      auto h1recPur = new TH1F(Form("purity_%s_%s", prefHisto[itp], prefAxis[iax]), Form("purity, %s-tracks vs %s;%s", prefHisto[itp], prefAxis[iax], prefAxis[iax]), nbinsAxis[iax], cutAxisMin[itp][iax], cutAxisMax[itp][iax]);
      hm->addHisto(h1mc, offset + iax * 100 + 0);
      hm->addHisto(h1rec, offset + iax * 100 + 1);
      hm->addHisto(h1recCorr, offset + iax * 100 + 2);
      hm->addHisto(h1recEff, offset + iax * 100 + 3);
      hm->addHisto(h1recPur, offset + iax * 100 + 4);

      //
      // resolutions
      for (int ip = 0; ip < 5; ip++) {
        auto h2 = new TH2F(Form("difRecMC_%s_vs_%s_%s", parName[ip], prefAxis[iax], prefHisto[itp]),
                           Form("%s Rec-MC %s vs %s;%s;%s data-MC", prefHisto[itp], parName[ip], prefAxis[iax], prefAxis[iax], parName[ip]),
                           nbinsAxis[iax], cutAxisMin[itp][iax], cutAxisMax[itp][iax], nBinsDiffKin, -parMaxDiff[ip], parMaxDiff[ip]);
        hm->addHisto(h2, offset + iax * 100 + 500 + ip * 10);
      }
    }
  }
  if (checkDimu) {
    int offset = HMOffsets[3];
    auto hgMS = new TH1F("reconstructableDimuMS_rapidity", "MS reconstructable Dimu vs rapidity;rapidity", nbinsAxis[1], cutAxisMin[2][1], cutAxisMax[2][1]);
    hm->addHisto(hgMS, offset + 0);
    auto hrMS = new TH1F("reconstructedDimuMS_rapidity", "MS reconstructed Dimu vs rapidity;rapidity", nbinsAxis[1], cutAxisMin[2][1], cutAxisMax[2][1]);
    hm->addHisto(hrMS, offset + 1);
    auto hrcMS = new TH1F("reconstructedCorrDimuMS_rapidity", "MS correctly reconstructed Dimu vs rapidity;rapidity", nbinsAxis[1], cutAxisMin[2][1], cutAxisMax[2][1]);
    hm->addHisto(hrcMS, offset + 2);
    auto hrcMSeff = new TH1F("receffDimuMS_rapidity", "MS reconstruction eff. Dimu vs rapidity;rapidity;eff.", nbinsAxis[1], cutAxisMin[2][1], cutAxisMax[2][1]);
    hm->addHisto(hrcMSeff, offset + 3);
    auto hrcMSpur = new TH1F("recpurDimuMS_rapidity", "MS reconstruction purity. Dimu vs rapidity;rapidity;purity", nbinsAxis[1], cutAxisMin[2][1], cutAxisMax[2][1]);
    hm->addHisto(hrcMSpur, offset + 4);

    auto hgMT = new TH1F("matchableDimu_rapidity", "Matchable Dimu vs rapidity;rapidity", nbinsAxis[1], cutAxisMin[2][1], cutAxisMax[2][1]);
    hm->addHisto(hgMT, offset + 10);
    auto hrMT = new TH1F("matchedDimu_rapidity", "Matched Dimu vs rapidity;rapidity", nbinsAxis[1], cutAxisMin[2][1], cutAxisMax[2][1]);
    hm->addHisto(hrMT, offset + 11);
    auto hrcMT = new TH1F("matchedCorrDimu_rapidity", "Matched correctly Dimu vs rapidity;rapidity", nbinsAxis[1], cutAxisMin[2][1], cutAxisMax[2][1]);
    hm->addHisto(hrcMT, offset + 12);
    auto hrcMTeff = new TH1F("matcheffDimuMS_rapidity", "Matching eff. Dimu vs rapidity;rapidity;eff.", nbinsAxis[1], cutAxisMin[2][1], cutAxisMax[2][1]);
    hm->addHisto(hrcMTeff, offset + 13);
    auto hrcMTpur = new TH1F("matchpurDimuMS_rapidity", "Mathing purity. Dimu vs rapidity;rapidity;purity", nbinsAxis[1], cutAxisMin[2][1], cutAxisMax[2][1]);
    hm->addHisto(hrcMTpur, offset + 14);

    auto massResMS = new TH1F("massresDimuMS", "MS reconstructed mass res;M_{rec}-M_{gen}", nBinsMassRes, -maxMassDiff, maxMassDiff);
    hm->addHisto(massResMS, offset + 100);
    auto massResCorrMS = new TH1F("massresCorrDimuMS", "MS correctly reconstructed mass res;M_{rec}-M_{gen}", nBinsMassRes, -maxMassDiff, maxMassDiff);
    hm->addHisto(massResCorrMS, offset + 101);

    auto massResMatch = new TH1F("massresDimuMatch", "Matched mass res;M_{rec}-M_{gen}", nBinsMassRes, -maxMassDiff, maxMassDiff);
    hm->addHisto(massResMatch, offset + 110);
    auto massResCorrMatch = new TH1F("massresCorrDimuMatch", "Correctly matched mass res;M_{rec}-M_{gen}", nBinsMassRes, -maxMassDiff, maxMassDiff);
    hm->addHisto(massResCorrMatch, offset + 111);
  }
  hm->sumw2();
}

void finalizeHistos(na6p::HistoManager* hman)
{
  std::unique_ptr<TF1> gs = std::make_unique<TF1>("gs", "gaus", -10., 10.);
  TObjArray harr;
  harr.SetOwner(true);
  for (int itp = 0; itp < 3; itp++) {
    int offset = HMOffsets[itp];
    if (!hman->getHisto(offset + 0 * 100 + 0)) {
      continue;
    }
    for (int iax = 0; iax < 2; iax++) {
      hman->getHisto(offset + iax * 100 + 3)->Divide(hman->getHisto(offset + iax * 100 + 1), hman->getHisto(offset + iax * 100 + 0), 1., 1., "B");
      hman->getHisto(offset + iax * 100 + 4)->Divide(hman->getHisto(offset + iax * 100 + 2), hman->getHisto(offset + iax * 100 + 1), 1., 1., "B");
      //
      // resolutions
      for (int ip = 0; ip < 5; ip++) {
        int hid = offset + iax * 100 + 500 + ip * 10;
        auto h2 = hman->getHisto2F(hid);
        h2->FitSlicesY(gs.get(), 0, -1, 0, "QNR", &harr);
        harr.SetOwner(true);
        TH1* hmean = (TH1*)harr.RemoveAt(1);
        if (hmean) {
          hmean->SetTitle(Form("<%s>", h2->GetTitle()));
          hman->addHisto(hmean, hid + 1);
        }
        TH1* hsig = (TH1*)harr.RemoveAt(2);
        if (hsig) {
          hsig->SetTitle(Form("#sigma(%s)", h2->GetTitle()));
          hman->addHisto(hsig, hid + 2);
        }
        harr.Delete();
      }
    }
  }
  // dimuons
  {
    int offset = HMOffsets[3];
    if (hman->getHisto(offset + 0)) {
      hman->getHisto(offset + 3)->Divide(hman->getHisto(offset + 1), hman->getHisto(offset + 0), 1., 1., "B"); // eff dimu MS
      hman->getHisto(offset + 4)->Divide(hman->getHisto(offset + 2), hman->getHisto(offset + 1), 1., 1., "B"); // purity dimu MS

      hman->getHisto(offset + 13)->Divide(hman->getHisto(offset + 11), hman->getHisto(offset + 10), 1., 1., "B"); // eff dimu Match
      hman->getHisto(offset + 14)->Divide(hman->getHisto(offset + 12), hman->getHisto(offset + 11), 1., 1., "B"); // purity dimu Match

      if (hman->getHisto(offset + 100)->GetEntries() > 50) {
        hman->getHisto(offset + 100)->Fit(gs.get(), "q", "");
      }
      if (hman->getHisto(offset + 101)->GetEntries() > 50) {
        hman->getHisto(offset + 101)->Fit(gs.get(), "q", "");
      }
      if (hman->getHisto(offset + 110)->GetEntries() > 50) {
        hman->getHisto(offset + 110)->Fit(gs.get(), "q", "");
      }
      if (hman->getHisto(offset + 111)->GetEntries() > 50) {
        hman->getHisto(offset + 111)->Fit(gs.get(), "q", "");
      }
    }
  }
}

TCanvas* cnvEff[3] = {};
TCanvas* cnvKin[3][2] = {{0, 0}, {0, 0}, {0, 0}};
TCanvas* cnvDimuEff = nullptr;
TCanvas* cnvDimuRes = nullptr;

void drawHistos(na6p::HistoManager* hman, int copy, bool reproc)
{
  auto drawHistoA = [copy, hman](int id, float mn = 1e6, float mx = -1e6, float mrgH = 0.15, float mrgL = 0.15) {
    auto h = hman->getHisto(id);
    gPad->SetBottomMargin(0.15);
    gPad->SetLeftMargin(0.15);
    h->Draw(copy > 0 ? "same" : "");
    h->GetXaxis()->SetLabelSize(0.04);
    h->GetYaxis()->SetLabelSize(0.04);
    h->GetXaxis()->SetTitleSize(0.05);
    h->GetYaxis()->SetTitleSize(0.05);
    gPad->SetGrid();
    gPad->Modified();
    if (mn > mx) {
      na6p::HistoManager::getHistosMinMaxRange(gPad, mn, mx); // estimate min max since the external one was not provided
    }
    na6p::HistoManager::setHistosMinMaxRange(gPad, mn, mx, mrgH, mrgL);
    auto stpad = na6p::HistoManager::getStatPad(h);
    if (stpad) {
      na6p::HistoManager::setStatPad(h, 0.7, 0.99, 0.77 - 0.15 * copy, 0.92 - 0.15 * copy);
    }
  };

  gStyle->SetTitleW(0.5);
  gStyle->SetOptStat(0);

  if (reproc) {
    finalizeHistos(hman);
  }

  const char* opt = copy == 0 ? "" : "same";
  for (int itp = 0; itp < 3; itp++) {
    int offset = HMOffsets[itp];
    if (!hman->getHisto(offset + 0 * 100 + 0)) {
      continue;
    }
    //--------------- Eff
    if (!cnvEff[itp]) {
      cnvEff[itp] = new TCanvas(Form("cnv%s", prefHisto[itp]), Form("cnv%s", prefHisto[itp]), 600, 800);
      cnvEff[itp]->Divide(2, 2);
    }
    for (int iax = 0; iax < 2; iax++) {
      cnvEff[itp]->cd(1 + 2 * iax);
      drawHistoA(offset + iax * 100 + 3, hrangeEff[itp][0], hrangeEff[itp][1]);
      cnvEff[itp]->cd(2 + 2 * iax);
      drawHistoA(offset + iax * 100 + 4, hrangePurity[itp][0], hrangePurity[itp][1]);
    }
    //--------------- Res

    for (int iax = 0; iax < 2; iax++) {
      if (!cnvKin[itp][iax]) {
        cnvKin[itp][iax] = new TCanvas(Form("cnv%s%s", prefHisto[itp], prefAxisS[iax]), Form("cnv%s%s", prefHisto[itp], prefAxisS[iax]), 600, 800);
        cnvKin[itp][iax]->Divide(2, 5);
      }
      for (int ip = 0; ip < 5; ip++) {
        int hid = offset + iax * 100 + 500 + ip * 10;
        cnvKin[itp][iax]->cd(1 + ip * 2);
        drawHistoA(hid + 1, hrangeParBias[ip][0], hrangeParBias[ip][1]);
        cnvKin[itp][iax]->cd(2 + ip * 2);
        drawHistoA(hid + 2, hrangeParSigm[ip][0], hrangeParSigm[ip][1], 0.15f, 0.f);
      }
    }
  }
  if (hman->getHisto(HMOffsets[3] + 0)) {
    int offset = HMOffsets[3];
    if (!cnvDimuEff) {
      cnvDimuEff = new TCanvas("cnvDimuEff", "cnvDimuEff", 600, 800);
      cnvDimuEff->Divide(2, 2);
    }
    for (int i = 0; i < 2; i++) {
      cnvDimuEff->cd(1 + i * 2);
      drawHistoA(offset + i * 10 + 3, hrangeEff[3][0], hrangeEff[3][1]);
      cnvDimuEff->cd(1 + i * 2 + 1);
      drawHistoA(offset + i * 10 + 4, hrangePurity[3][0], hrangePurity[3][1]);
    }
    if (!cnvDimuRes) {
      cnvDimuRes = new TCanvas("cnvDimuRes", "cnvDimuRes", 600, 800);
      cnvDimuRes->Divide(2, 2);
    }
    for (int i = 0; i < 2; i++) {
      cnvDimuRes->cd(1 + i * 2);
      drawHistoA(offset + 100 + i * 10 + 0, 1.f, -1.f, 0.15, 0);

      cnvDimuRes->cd(1 + i * 2 + 1);
      drawHistoA(offset + 100 + i * 10 + 1, 1.f, -1.f, 0.15, 0);
    }
  }
}

void createOutput(bool checkVTTracks, bool checkMSTracks, bool checkMatches, bool checkDimu)
{
  NA6PTrackParCov dummy;
  dummy.invalidate();
  auto prop = Propagator::Instance();
  for (auto& mcInfoPair : mcTruth) {
    const auto& mcTrackInfo = mcInfoPair.second;
    if (mcTrackInfo.status != MCTrackInfo::Accept || !(mcTrackInfo.isReconstructableVT() || mcTrackInfo.isReconstructableMS()))
      continue;

    const auto& partMC = (*mcArr)[mcInfoPair.first.getTrackID()];
    // convert to TrackParam
    float mom[3] = {(float)partMC.Px(), (float)partMC.Py(), (float)partMC.Pz()}, pos[3] = {(float)partMC.Vx(), (float)partMC.Vy(), (float)partMC.Vz()};
    NA6PTrackPar trMC(pos, mom, mcTrackInfo.charge);
    float pInvMC = 1.f / trMC.getP(), rapidityMC = trMC.getRapidity();

    if (checkVTTracks && mcTrackInfo.isReconstructableVT()) { // output generated info at the vertex
      hm->getHisto(HMOffsets[0] + 0 * 100 + 0)->Fill(pInvMC); // reconstructable
      hm->getHisto(HMOffsets[0] + 1 * 100 + 0)->Fill(rapidityMC);
      NA6PTrackParCov trCopy;
      int nVTHits = 0, nClones = 0;
      NA6PMCComposedLabel lblTr;
      float chi2 = -1.f;
      if (mcTrackInfo.bestVTTrack >= 0) {
        const auto& vtTr = (*vtTrVec)[mcTrackInfo.bestVTTrack];
        trCopy = vtTr;
        nVTHits = vtTr.getNHits();
        nClones = mcTrackInfo.nVTTracks;
        lblTr = (*vtTrkLabels)[mcTrackInfo.bestVTTrack];
        chi2 = vtTr.getChi2();
        if (!prop->propagateToZ(trCopy, pos[2])) {
          trCopy.invalidate();
        }
        hm->getHisto(HMOffsets[0] + 0 * 100 + 1)->Fill(pInvMC); // reconstructed
        hm->getHisto(HMOffsets[0] + 1 * 100 + 1)->Fill(rapidityMC);
        //
        if (!lblTr.isFake()) {
          hm->getHisto(HMOffsets[0] + 0 * 100 + 2)->Fill(pInvMC); // reconstructed correctly
          hm->getHisto(HMOffsets[0] + 1 * 100 + 2)->Fill(rapidityMC);
        }
        if (trCopy.isValid()) { // resolutions
          for (int ip = 0; ip < 5; ip++) {
            hm->getHisto2F(HMOffsets[0] + 0 * 100 + 500 + ip * 10)->Fill(pInvMC, trCopy.getParam(ip) - trMC.getParam(ip));     // reconstructed
            hm->getHisto2F(HMOffsets[0] + 1 * 100 + 500 + ip * 10)->Fill(rapidityMC, trCopy.getParam(ip) - trMC.getParam(ip)); // reconstructed
          }
        }
      } else {
        trCopy = dummy;
      }
      (*outPerf) << "perfVT" << "evID=" << int(mcTree->GetReadEntry()) << "mcTrID=" << mcInfoPair.first << "pdg=" << partMC.GetPdgCode() << "mcTr=" << trMC
                 << "nvtLr=" << mcTrackInfo.nvtLr << "nmsLr=" << mcTrackInfo.nmsLr
                 << "vtTr=" << trCopy << "nVTCl=" << nVTHits << "chi2vt=" << chi2 << "nClones=" << nClones << "mcID=" << lblTr << "\n";
    }
    //-----------------------

    if (checkMSTracks && mcTrackInfo.isReconstructableMS()) { // output generated info at the vertex
      NA6PTrackParCov trCopy;
      NA6PTrackPar msMCTr(mcTrackInfo.posMS, mcTrackInfo.momMS, mcTrackInfo.charge);
      hm->getHisto(HMOffsets[1] + 0 * 100 + 0)->Fill(pInvMC); // reconstructable
      hm->getHisto(HMOffsets[1] + 1 * 100 + 0)->Fill(rapidityMC);

      int nMSHits = 0, nClones = 0;
      NA6PMCComposedLabel lblTr;
      float chi2 = -1.f;
      if (mcTrackInfo.bestMSTrack >= 0) {
        const auto& msTr = (*msTrVec)[mcTrackInfo.bestMSTrack];
        trCopy = msTr;
        nMSHits = msTr.getNHits();
        nClones = mcTrackInfo.nMSTracks;
        lblTr = (*msTrkLabels)[mcTrackInfo.bestMSTrack];
        chi2 = msTr.getChi2();
        if (!prop->propagateToZ(trCopy, msMCTr.getZ())) {
          trCopy.invalidate();
        }
        hm->getHisto(HMOffsets[1] + 0 * 100 + 1)->Fill(pInvMC); // reconstructed
        hm->getHisto(HMOffsets[1] + 1 * 100 + 1)->Fill(rapidityMC);
        //
        if (!lblTr.isFake()) {
          hm->getHisto(HMOffsets[1] + 0 * 100 + 2)->Fill(pInvMC); // reconstructed correctly
          hm->getHisto(HMOffsets[1] + 1 * 100 + 2)->Fill(rapidityMC);
        }
        if (trCopy.isValid()) { // resolutions
          for (int ip = 0; ip < 5; ip++) {
            hm->getHisto2F(HMOffsets[1] + 0 * 100 + 500 + ip * 10)->Fill(pInvMC, trCopy.getParam(ip) - msMCTr.getParam(ip));     // reconstructed
            hm->getHisto2F(HMOffsets[1] + 1 * 100 + 500 + ip * 10)->Fill(rapidityMC, trCopy.getParam(ip) - msMCTr.getParam(ip)); // reconstructed
          }
        }
      } else {
        trCopy = dummy;
      }
      (*outPerf) << "perfMS" << "evID=" << int(mcTree->GetReadEntry()) << "mcTrID=" << mcInfoPair.first << "pdg=" << partMC.GetPdgCode()
                 << "mcTrMS0=" << msMCTr << "nvtLr=" << mcTrackInfo.nvtLr << "nmsLr=" << mcTrackInfo.nmsLr
                 << "msTr=" << trCopy << "nMSCl=" << nMSHits << "chi2ms=" << chi2 << "nClones=" << nClones << "mcID=" << lblTr << "\n";
    }
    //-----------------------
    if (checkMatches && mcTrackInfo.isReconstructableVT() && mcTrackInfo.isReconstructableMS()) { // output generated info at the vertex
      NA6PTrackParCov trCopy;
      int nClones = 0;
      NA6PMCComposedLabel lblTr;
      float chi2refit = -1.f, chi2match = -1.f;
      hm->getHisto(HMOffsets[2] + 0 * 100 + 0)->Fill(pInvMC); // matchable
      hm->getHisto(HMOffsets[2] + 1 * 100 + 0)->Fill(rapidityMC);
      if (mcTrackInfo.bestMatch >= 0) {
        const auto& mtcTr = (*matchVec)[mcTrackInfo.bestMatch];
        trCopy = mtcTr;
        nClones = mcTrackInfo.nMatches;
        lblTr = (*matchTrkLabels)[mcTrackInfo.bestMatch];
        chi2refit = mtcTr.getChi2Refit();
        chi2match = mtcTr.getChi2Match();
        if (!prop->propagateToZ(trCopy, pos[2])) {
          trCopy.invalidate();
        }
        hm->getHisto(HMOffsets[2] + 0 * 100 + 1)->Fill(pInvMC); // matched
        hm->getHisto(HMOffsets[2] + 1 * 100 + 1)->Fill(rapidityMC);
        //
        if (!lblTr.isFake()) {
          hm->getHisto(HMOffsets[2] + 0 * 100 + 2)->Fill(pInvMC); // matched correctly
          hm->getHisto(HMOffsets[2] + 1 * 100 + 2)->Fill(rapidityMC);
        }
        if (trCopy.isValid()) { // resolutions
          for (int ip = 0; ip < 5; ip++) {
            hm->getHisto2F(HMOffsets[2] + 0 * 100 + 500 + ip * 10)->Fill(pInvMC, trCopy.getParam(ip) - trMC.getParam(ip));
            hm->getHisto2F(HMOffsets[2] + 1 * 100 + 500 + ip * 10)->Fill(rapidityMC, trCopy.getParam(ip) - trMC.getParam(ip));
          }
        }
      } else {
        trCopy = dummy;
      }
      (*outPerf) << "perfMatch" << "evID=" << int(mcTree->GetReadEntry()) << "mcTrID=" << mcInfoPair.first << "pdg=" << partMC.GetPdgCode() << "mcTr=" << trMC
                 << "mtcTr=" << trCopy << "chi2refit=" << chi2refit << "chi2match=" << chi2match << "nClones=" << nClones << "mcID=" << lblTr << "\n";
    }
    //-----------------------
  }

  // dimuons
  if (checkDimu) {
    auto it0 = mcTruth.begin();
    while (it0 != mcTruth.end()) {
      const auto& mcTrackInfo0 = it0->second;
      if (!mcTrackInfo0.isReconstructableMS()) {
        it0++;
        continue;
      }
      const auto& partMC0 = (*mcArr)[it0->first.getTrackID()];
      TLorentzVector muG0;
      float mom0[3] = {(float)partMC0.Px(), (float)partMC0.Py(), (float)partMC0.Pz()}, pos0[3] = {(float)partMC0.Vx(), (float)partMC0.Vy(), (float)partMC0.Vz()};
      NA6PTrackPar muG0Tr(pos0, mom0, mcTrackInfo0.charge);
      muG0Tr.setPID(PID::Muon);
      muG0.SetXYZM(partMC0.Px(), partMC0.Py(), partMC0.Pz(), 0.105658369);
      TLorentzVector muMS0, muMatch0;
      NA6PMCComposedLabel muMSRec0ID, muMatchRec0ID;
      if (std::abs(partMC0.GetPdgCode()) != 13 || !mcTrackInfo0.parent.isSet()) {
        it0++;
        continue;
      }
      NA6PTrackPar muMSRec0Par;
      if (mcTrackInfo0.bestMSTrack >= 0) {
        muMSRec0Par = (*msTrVec)[mcTrackInfo0.bestMSTrack];
        if (prop->propagateToZ(muMSRec0Par, partMC0.Vz())) {
          muMS0.SetXYZM(muMSRec0Par.getPx(), muMSRec0Par.getPy(), muMSRec0Par.getPz(), 0.105658369);
          muMSRec0ID = (*msTrkLabels)[mcTrackInfo0.bestMSTrack];
        }
      }
      bool matchable0 = false;
      NA6PTrackPar muMatchRec0Par;
      if (mcTrackInfo0.isReconstructableVT()) {
        matchable0 = true;
        if (mcTrackInfo0.bestMatch >= 0) {
          muMatchRec0Par = (*matchVec)[mcTrackInfo0.bestMatch];
          if (prop->propagateToZ(muMatchRec0Par, partMC0.Vz())) {
            muMatch0.SetXYZM(muMatchRec0Par.getPx(), muMatchRec0Par.getPy(), muMatchRec0Par.getPz(), 0.105658369);
            muMatchRec0ID = (*matchTrkLabels)[mcTrackInfo0.bestMatch];
          }
        }
      }

      auto it1 = it0;
      while (++it1 != mcTruth.end()) {
        const auto& mcTrackInfo1 = it1->second;
        if (!mcTrackInfo1.isReconstructableMS()) {
          continue;
        }
        if (mcTrackInfo0.parent != mcTrackInfo1.parent || !mcTrackInfo0.parent.isSet()) {
          continue;
        }
        const auto& partMC1 = (*mcArr)[it1->first.getTrackID()];
        if (std::abs(partMC1.GetPdgCode()) != 13 || partMC0.GetPdgCode() * partMC1.GetPdgCode() > 0) {
          continue;
        }
        TLorentzVector muG1;
        float mom1[3] = {(float)partMC1.Px(), (float)partMC1.Py(), (float)partMC1.Pz()}, pos1[3] = {(float)partMC1.Vx(), (float)partMC1.Vy(), (float)partMC1.Vz()};
        NA6PTrackPar muG1Tr(pos1, mom1, mcTrackInfo1.charge);
        muG1Tr.setPID(PID::Muon);
        muG1.SetXYZM(partMC1.Px(), partMC1.Py(), partMC1.Pz(), 0.105658369);
        TLorentzVector muMS1, muMatch1;
        NA6PMCComposedLabel muMSRec1ID, muMatchRec1ID;
        NA6PTrackPar muMSRec1Par;
        if (mcTrackInfo1.bestMSTrack >= 0) {
          muMSRec1Par = (*msTrVec)[mcTrackInfo1.bestMSTrack];
          if (prop->propagateToZ(muMSRec1Par, partMC1.Vz())) {
            muMS1.SetXYZM(muMSRec1Par.getPx(), muMSRec1Par.getPy(), muMSRec1Par.getPz(), 0.105658369);
            muMSRec1ID = (*msTrkLabels)[mcTrackInfo1.bestMSTrack];
          }
        }
        bool matchable1 = false;
        NA6PTrackPar muMatchRec1Par;
        if (mcTrackInfo1.isReconstructableVT()) {
          matchable1 = true;
          if (mcTrackInfo1.bestMatch >= 0) {
            muMatchRec1Par = (*matchVec)[mcTrackInfo1.bestMatch];
            if (prop->propagateToZ(muMatchRec1Par, partMC1.Vz())) {
              muMatch1.SetXYZM(muMatchRec1Par.getPx(), muMatchRec1Par.getPy(), muMatchRec1Par.getPz(), 0.105658369);
              muMatchRec1ID = (*matchTrkLabels)[mcTrackInfo1.bestMatch];
            }
          }
        }

        TLorentzVector muParentG = muG0, muParentMSRec, muParentMatchRec;
        muParentG += muG1;
        if (muMSRec0ID.isSet() && muMSRec1ID.isSet()) {
          muParentMSRec = muMS0;
          muParentMSRec += muMS1;
        }
        if (matchable0 && matchable1 && muMatchRec0ID.isSet() && muMatchRec1ID.isSet()) {
          muParentMatchRec = muMatch0;
          muParentMatchRec += muMatch1;
        }
        //------------------- fill
        int offset = HMOffsets[3];
        auto dimuRapidityMC = muParentG.Rapidity(), dimuMMC = muParentG.M();
        hm->getHisto(offset + 0)->Fill(dimuRapidityMC);
        if (muMSRec0ID.isSet() && muMSRec1ID.isSet()) {
          hm->getHisto(offset + 1)->Fill(dimuRapidityMC);
          hm->getHisto(offset + 100)->Fill(muParentMSRec.M() - dimuMMC);
          if (muMSRec0ID == it0->first && muMSRec1ID == it1->first) {
            hm->getHisto(offset + 2)->Fill(dimuRapidityMC);
            hm->getHisto(offset + 101)->Fill(muParentMSRec.M() - dimuMMC);
          }
        }

        if (matchable0 && matchable1) {
          hm->getHisto(offset + 10)->Fill(dimuRapidityMC);
          if (muMatchRec0ID.isSet() && muMatchRec1ID.isSet()) {
            hm->getHisto(offset + 11)->Fill(dimuRapidityMC);
            hm->getHisto(offset + 110)->Fill(muParentMatchRec.M() - dimuMMC);
            if (muMatchRec0ID == it0->first && muMatchRec1ID == it1->first) {
              hm->getHisto(offset + 12)->Fill(dimuRapidityMC);
              hm->getHisto(offset + 111)->Fill(muParentMatchRec.M() - dimuMMC);
            }
          }
        }

        (*outPerf) << "perfDimu" << "evID=" << int(mcTree->GetReadEntry()) << "genDimu=" << muParentG << "matchable=" << (matchable0 && matchable1)
                   << "recDimuMS=" << muParentMSRec << "recMuMS0ID=" << muMSRec0ID << "recMuMS1ID=" << muMSRec1ID
                   << "recDimuMatch=" << muParentMatchRec << "recMatch0ID=" << muMatchRec0ID << "recMatch1ID=" << muMatchRec1ID
                   << "genMu0=" << muG0 << "genMu1=" << muG1 << "genMu0ID=" << it0->first << "genMu1ID=" << it1->first
                   << "recMuMS0=" << muMS0 << "recMuMS1=" << muMS1
                   << "recMuMatch0=" << muMatch0 << "recMuMatch1=" << muMatch1
                   << "matchable0=" << matchable0 << "matchable1=" << matchable1
                   << "genMu0Par=" << muG0Tr << "genMu1Par=" << muG1Tr
                   << "recMuMS0Par=" << muMSRec0Par << "recMuMS1Par=" << muMSRec1Par
                   << "recMuMatch0Par=" << muMatchRec0Par << "recMuMatch1Par=" << muMatchRec1Par
                   << "\n";
      }
      it0++;
    }
  }
}

void checkVTTracking()
{
  // match tracks in the vtTrVec to MC info
  int ntr = vtTrVec->size();
  for (int itr = 0; itr < ntr; itr++) {
    const auto& tr = (*vtTrVec)[itr];
    auto trMCID = (*vtTrkLabels)[itr];
    auto& mcInfo = mcTruth[trMCID];
    mcInfo.nVTTracks++;
    if (mcInfo.bestVTTrack < 0) {
      mcInfo.bestVTTrack = itr;
    } else {
      const auto& trPrev = (*vtTrVec)[mcInfo.bestVTTrack];
      if (trPrev.getNHits() < tr.getNHits() || (trPrev.getNHits() == tr.getNHits() && trPrev.getChi2() > tr.getChi2())) {
        mcInfo.bestVTTrack = itr;
      }
    }
  }
}

void checkMSTracking()
{
  // match tracks in the msTrVec to MC info
  int ntr = msTrVec->size();
  for (int itr = 0; itr < ntr; itr++) {
    const auto& tr = (*msTrVec)[itr];
    auto trMCID = (*msTrkLabels)[itr];
    auto& mcInfo = mcTruth[trMCID];
    mcInfo.nMSTracks++;
    if (mcInfo.bestMSTrack < 0) {
      mcInfo.bestMSTrack = itr;
    } else {
      const auto& trPrev = (*msTrVec)[mcInfo.bestMSTrack];
      if (trPrev.getNHits() < tr.getNHits() || (trPrev.getNHits() == tr.getNHits() && trPrev.getChi2() > tr.getChi2())) {
        mcInfo.bestMSTrack = itr;
      }
    }
  }
}

void checkMatching()
{
  // match tracks in the matchVec to MC info
  int ntr = matchVec->size();
  for (int itr = 0; itr < ntr; itr++) {
    const auto& tr = (*matchVec)[itr];
    auto trMCID = (*matchTrkLabels)[itr];
    auto& mcInfo = mcTruth[trMCID];
    mcInfo.nMatches++;
    if (mcInfo.bestMatch < 0) {
      mcInfo.bestMatch = itr;
    } else {
      const auto& trPrev = (*matchVec)[mcInfo.bestMatch];
      if (trPrev.getChi2Match() > tr.getChi2Match()) {
        mcInfo.bestMatch = itr;
      }
    }
  }
}

void setupInputs(bool checkVTTracks, bool checkMSTracks, bool checkMatches, const char* dirSimu, const char* dirRec)
{
  // kine file needed to get the primary vertex position (temporary)
  auto fetchTree = [&](const char* fname, const char* tname, TTree*& tree) {
    auto f = TFile::Open(Form("%s/%s", dirSimu, fname));
    if (!f || f->IsZombie() || !(tree = (TTree*)f->Get(tname))) {
      exit(1);
    }
    nEntries = nEntries > 0 ? std::min(tree->GetEntries(), nEntries) : tree->GetEntries();
  };

  fetchTree(Form("%s/MCKine.root", dirSimu), "mckine", mcTree);
  mcTree->SetBranchAddress("tracks", &mcArr);

  if (checkVTTracks || checkMatches) {
    fetchTree(Form("%s/ClustersVerTel.root", dirSimu), "clustersVerTel", vtClTree);
    vtClTree->SetBranchAddress("VerTel", &vtClVec);
    vtClTree->SetBranchAddress("VerTelMCTruth", &vtCluMCLabels);
    fetchTree(Form("%s/TracksVerTel.root", dirSimu), "tracksVerTel", vtTrTree);
    vtTrTree->SetBranchAddress("VerTel", &vtTrVec);
    vtTrTree->SetBranchAddress("VerTelMCTruth", &vtTrkLabels);
  }

  if (checkMSTracks || checkMatches) {
    fetchTree(Form("%s/ClustersMuonSpec.root", dirSimu), "clustersMuonSpec", msClTree);
    msClTree->SetBranchAddress("MuonSpec", &msClVec);
    msClTree->SetBranchAddress("MuonSpecMCTruth", &msCluMCLabels);

    fetchTree(Form("%s/TracksMuonSpec.root", dirSimu), "tracksMuonSpec", msTrTree);
    msTrTree->SetBranchAddress("MuonSpec", &msTrVec);
    msTrTree->SetBranchAddress("MuonSpecMCTruth", &msTrkLabels);

    fetchTree(Form("%s/HitsMuonSpec.root", dirSimu), "hitsMuonSpec", msHitsTree);
    msHitsTree->SetBranchAddress("MuonSpec", &msHitsVec);
  }

  if (checkMatches) {
    fetchTree(Form("%s/TracksMatching.root", dirSimu), "tracksMatching", matchTree);
    matchTree->SetBranchAddress("Matching", &matchVec);
    matchTree->SetBranchAddress("MatchingMCTruth", &matchTrkLabels);
  }
}

void loadEventData(int i)
{
  if (mcTree)
    mcTree->GetEntry(i);
  if (vtClTree)
    vtClTree->GetEntry(i);
  if (msClTree)
    msClTree->GetEntry(i);
  if (vtTrTree)
    vtTrTree->GetEntry(i);
  if (msTrTree)
    msTrTree->GetEntry(i);
  if (msHitsTree)
    msHitsTree->GetEntry(i);
  if (matchTree)
    matchTree->GetEntry(i);
}

void buildReconstuctableMCInfo()
{
  auto consideMCPart = [&](NA6PMCComposedLabel lbl) {
    auto ent0 = mcTruth.emplace(lbl, MCTrackInfo{});
    auto& ent = ent0.first->second;
    if (ent.status == MCTrackInfo::New) {
      const auto& part = (*mcArr)[lbl.getTrackID()]; // RS: here we assume that the lbl.getEventID is already loaded
      while (true) {
        ent.status = MCTrackInfo::Discard;
        auto ppdg = TDatabasePDG::Instance()->GetParticle(part.GetPdgCode());
        if (!ppdg || !ppdg->Charge()) {
          break;
        }
        auto rapidity = part.Y();
        if (rapidity < cutRapidityMin || rapidity > cutRapidityMax) {
          break;
        }
        auto pt = part.Pt();
        if (pt < cutPtMin) {
          break;
        }
        // all checks passed
        ent.status = MCTrackInfo::Accept;
        ent.charge = ppdg->Charge() / 3;
        if (part.GetFirstMother() >= 0) {
          ent.parent = NA6PMCComposedLabel(part.GetFirstMother(), mcTree->GetReadEntry(), 0);
        }
        break;
      }
    }
    return ent.status == MCTrackInfo::Discard ? mcTruth.end() : ent0.first; // return iterator
  };

  mcTruth.clear();
  if (vtClVec) {
    for (int i = 0; i < (int)vtClVec->size(); i++) {
      auto labels = vtCluMCLabels->getLabels(i);
      for (const auto currLbl : labels) {
        auto ent = consideMCPart(currLbl);
        if (ent != mcTruth.end()) {
          ent->second.vtClID.push_back(i);
        }
      }
    }
  }
  if (msClVec) {
    for (int i = 0; i < (int)msClVec->size(); i++) {
      auto labels = msCluMCLabels->getLabels(i);
      for (const auto currLbl : labels) {
        auto ent = consideMCPart(currLbl);
        if (ent != mcTruth.end()) {
          ent->second.msClID.push_back(i);
        }
      }
    }
  }
  // sort clusters in layer increasing order, count layers
  for (auto& mci : mcTruth) {
    {
      auto& vtc = mci.second.vtClID;
      std::sort(vtc.begin(), vtc.end(), [](int i, int j) { return (*vtClVec)[i].getLayer() < (*vtClVec)[j].getLayer(); });
      int prevL = -1;
      for (auto i : vtc) {
        if ((*vtClVec)[i].getLayer() != prevL) {
          mci.second.nvtLr++;
          prevL = (*vtClVec)[i].getLayer();
        }
      }
    }
    {
      auto& msc = mci.second.msClID;
      std::sort(msc.begin(), msc.end(), [](int i, int j) { return (*msClVec)[i].getLayer() < (*msClVec)[j].getLayer(); });
      int prevL = -1;
      for (auto i : msc) {
        if ((*msClVec)[i].getLayer() != prevL) {
          mci.second.nmsLr++;
          prevL = (*msClVec)[i].getLayer();
        }
      }
      if (!msc.empty()) {
        const auto& cl = (*msClVec)[msc.front()];
        auto hit = (*msHitsVec)[cl.getHitID()];
        mci.second.momMS = std::array<float, 3>{(float)hit.getMomIn()[0], (float)hit.getMomIn()[1], (float)hit.getMomIn()[2]};
        mci.second.posMS = std::array<float, 3>{(float)hit.getPosIn()[0], (float)hit.getPosIn()[1], (float)hit.getPosIn()[2]};
      }
    }
  }
}

void createPDF(const char* outName)
{
  TCanvas dummy("dum", "dum", 600, 800);
  dummy.Print(Form("%s[", outName));
  for (int itp = 0; itp < 3; itp++) {
    if (cnvEff[itp]) {
      cnvEff[itp]->Print(outName);
    }
    for (int iax = 0; iax < 2; iax++) {
      if (cnvKin[itp][iax]) {
        cnvKin[itp][iax]->Print(outName);
      }
    }
  }
  if (cnvDimuEff) {
    cnvDimuEff->Print(outName);
  }
  if (cnvDimuRes) {
    cnvDimuRes->Print(outName);
  }
  dummy.Print(Form("%s]", outName));
}
