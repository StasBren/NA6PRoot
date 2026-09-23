/*
  CA Track Finder Macro

  Usage:
    root -l
    root [0] .L setup_macros.C     // Load NA6PRoot libraries and setup paths
    root [1] .L runMSTrackFinderCA.C+  // Compile the macro with ACLiC
    root [2] runMSTrackFinderCA()      // Run with default parameters

  Or with custom parameters:
    root [2] runMSTrackFinderCA(0, 100, "../fullgeo", 5)
*/

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TTree.h>
#include <TFile.h>
#include <TH1F.h>
#include <TVector3.h>
#include <TParticle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include "NA6PMuonSpecCluster.h"
#include "NA6PTrack.h"
#include "NA6PVertex.h"
#include "NA6PTrackerCA.h"
#include "MagneticField.h"
#include "NA6PMuonSpecHit.h"
#include "NA6PRecoParam.h"
#include <TDatabasePDG.h>
#include <TParticlePDG.h>
#include "NA6PMCTruthContainer.h"
#include "NA6PReconstruction.h"
#endif

const int maxIterationsCA = NA6PTrackerCA::kMaxIterationsCA;

void runMSTrackFinderCA(int firstEv = 0,
                        int lastEv = 9999,
                        const char* dirSimu = ".",
                        int nLayers = 6)
{
  NA6PRecoParam::setValue("reco.msNLayers", std::to_string(nLayers));
  if (!Propagator::loadField() || !Propagator::loadGeometry(Form("%s/geometry.root", dirSimu))) {
    return;
  }
  auto dbPDG = TDatabasePDG::Instance();
  int nMomBins = 40;
  float pMax = 30.0;
  TH1F* hMomGen = new TH1F("hMomGen", ";p (GeV/c);counts", nMomBins, 0., pMax);
  TH1F* hEtaGen = new TH1F("hEtaGen", ";#eta;counts", 20, 1., 5.);
  TH1F* hMomReco = new TH1F("hMomReco", ";p (GeV/c);counts", nMomBins, 0., pMax);
  TH1F* hEtaReco = new TH1F("hEtaReco", ";#eta;counts", 20, 1., 5.);
  TH1F* hMomGoodReco = new TH1F("hMomGoodReco", ";p (GeV/c);counts", nMomBins, 0., pMax);
  TH1F* hEtaGoodReco = new TH1F("hEtaGoodReco", ";#eta;counts", 20, 1., 5.);
  TH1F** hMomRecoIterCA = new TH1F*[maxIterationsCA];
  TH1F** hEtaRecoIterCA = new TH1F*[maxIterationsCA];
  int colors[maxIterationsCA] = {kMagenta + 1, kBlue, kGreen + 1, kOrange + 1, kRed + 1, kRed, kRed - 9, kGray + 2, kGray + 1, kGray};
  for (int jIteration = 0; jIteration < maxIterationsCA; ++jIteration) {
    hMomRecoIterCA[jIteration] = new TH1F(Form("hMomRecoIterCA%d", jIteration), ";p (GeV/c);counts", nMomBins, 0., pMax);
    hEtaRecoIterCA[jIteration] = new TH1F(Form("hEtaRecoIterCA%d", jIteration), ";#eta;counts", 20, 1., 5.);
  }

  std::unique_ptr<NA6PTrackerCA> tracker = std::make_unique<NA6PTrackerCA>();
  tracker->configureFromRecoParamMS();
  tracker->setVerbosity(false);
  tracker->printConfiguration();

  TFile* fk = new TFile(Form("%s/MCKine.root", dirSimu));
  TTree* mcTree = (TTree*)fk->Get("mckine");
  std::vector<TParticle>* mcArr = nullptr;
  mcTree->SetBranchAddress("tracks", &mcArr);

  TFile* fh = new TFile(Form("%s/HitsMuonSpec.root", dirSimu));
  TTree* th = (TTree*)fh->Get("hitsMuonSpec");
  std::vector<NA6PMuonSpecHit> msHits, *msHitsPtr = &msHits;
  th->SetBranchAddress("MuonSpec", &msHitsPtr);

  TFile* fc = new TFile(Form("%s/ClustersMuonSpec.root", dirSimu));
  printf("Open cluster file: %s\n", fc->GetName());
  TTree* tc = (TTree*)fc->Get("clustersMuonSpec");
  std::vector<NA6PMuonSpecCluster> msClus, *msClusPtr = &msClus;
  NA6PMCTruthContainer msCluMCLabels, *msCluMCLabelsPtr = &msCluMCLabels;
  tc->SetBranchAddress("MuonSpec", &msClusPtr);
  tc->SetBranchAddress("MuonSpecMCTruth", &msCluMCLabelsPtr);
  int nEv = tc->GetEntries();
  if (lastEv > nEv || lastEv < 0)
    lastEv = nEv;
  if (firstEv < 0)
    firstEv = 0;
  NA6PVertex primVert;
  NA6PReconstruction rec("dummy");
  std::vector<NA6PMCComposedLabel> trkMCLabs;

  int nIterationsCA = tracker->getNIterations();
  std::vector<int> hitSelID;
  for (int jEv = firstEv; jEv < lastEv; jEv++) {
    //    printf("\nEvent %d *************\n", jEv);
    mcTree->GetEvent(jEv);
    th->GetEvent(jEv);
    int nPart = mcArr->size();
    double zvert = 0;
    // get primary vertex position from the Kine Tree
    for (int jp = 0; jp < nPart; jp++) {
      auto curPart = mcArr->at(jp);
      if (curPart.IsPrimary()) {
        zvert = curPart.Vz();
        break;
      }
    }
    primVert.setXYZ(0., 0., zvert);
    uint nHits = msHits.size();
    for (int jp = 0; jp < nPart; jp++) {
      auto curPart = mcArr->at(jp);
      int maskHits = 0;
      int counter = 0;
      hitSelID.clear();
      for (size_t jHit = 0; jHit < nHits; ++jHit) {
        const auto& hit = msHits.at(jHit);
        int idPart = hit.getTrackID();
        if (idPart == jp) {
          counter++;
          hitSelID.push_back(jHit);
        }
      }
      if (nLayers == counter) {
        double pxPart = curPart.Px();
        double pyPart = curPart.Py();
        double pzPart = curPart.Pz();
        double momPart = curPart.P();
        double phiPart = curPart.Phi();
        double thetaPart = std::acos(pzPart / momPart);
        double etaPart = -std::log(std::tan(thetaPart / 2.));
        hMomGen->Fill(momPart);
        hEtaGen->Fill(etaPart);
        /* // uncomment this to print the state of reconstructable tracks
        for (int jHit = 0; jHit < counter; ++jHit) {
          const auto& hit = msHits.at(hitSelID[jHit]);
          const auto mom = hit.getMomIn(), pos = hit.getPosIn();
          int pdg = curPart.GetPdgCode();
          auto ppdg = dbPDG->GetParticle(pdg);
          int idPart = hit.getTrackID();
          NA6PTrackPar trMC(pos, mom, ppdg->Charge()/3.);
          printf("MC state Part%d Pdg=%d at Hit%zu: %s\n", idPart, pdg, jHit, trMC.asString().c_str());
        }
        */
      }
    }
    tc->GetEvent(jEv);
    for (size_t i = 0; i < msClus.size(); ++i) {
      msClus[i].setClusterIndex(static_cast<int>(i));
    }
    tracker->findTracks(msClus, &primVert);
    std::vector<NA6PTrack> trks = tracker->getTracks();
    rec.assignMCLabels(trks, trkMCLabs, msCluMCLabels);
    int nTrks = trks.size();
    for (int jT = 0; jT < nTrks; jT++) {
      const NA6PTrack& tr = trks[jT];
      NA6PMCComposedLabel mcCompLabel = trkMCLabs[jT];
      int idPartTrack = mcCompLabel.getTrackID();
      int jIteration = tr.getCAIteration();
      if (tr.getNHits() == 6) {
        auto curPart = mcArr->at(idPartTrack);
        double pxPart = curPart.Px();
        double pyPart = curPart.Py();
        double pzPart = curPart.Pz();
        double momPart = curPart.P();
        double phiPart = curPart.Phi();
        double thetaPart = std::acos(pzPart / momPart);
        double etaPart = -std::log(std::tan(thetaPart / 2.));
        hMomReco->Fill(momPart);
        hEtaReco->Fill(etaPart);
        // RS printf("found track#%d iter:%d partID:%d %s\n", jT, jIteration, idPartTrack, ((NA6PTrackPar&)tr).asString().c_str());
        if (jIteration >= 0 && jIteration < maxIterationsCA) {
          hMomRecoIterCA[jIteration]->Fill(momPart);
          hEtaRecoIterCA[jIteration]->Fill(etaPart);
        }
        if (!mcCompLabel.isFake()) {
          hMomGoodReco->Fill(momPart);
          hEtaGoodReco->Fill(etaPart);
        }
      }
    }
  }

  TH1F* hPurityMom = (TH1F*)hMomGoodReco->Clone("hPurityMom");
  hPurityMom->GetYaxis()->SetTitle("purity");
  for (int iBin = 1; iBin <= hMomGoodReco->GetNbinsX(); iBin++) {
    double cg = hMomGoodReco->GetBinContent(iBin);
    double ct = hMomReco->GetBinContent(iBin);
    if (ct == 0) {
      hPurityMom->SetBinContent(iBin, 0.);
      hPurityMom->SetBinError(iBin, 0.);
    } else {
      double p = cg / ct;
      double ep = std::sqrt(p * (1 - p) / ct);
      hPurityMom->SetBinContent(iBin, p);
      hPurityMom->SetBinError(iBin, ep);
    }
  }
  TH1F* hPurityEta = (TH1F*)hEtaGoodReco->Clone("hPurityEta");
  hPurityEta->GetYaxis()->SetTitle("purity");
  for (int iBin = 1; iBin <= hEtaGoodReco->GetNbinsX(); iBin++) {
    double cg = hEtaGoodReco->GetBinContent(iBin);
    double ct = hEtaReco->GetBinContent(iBin);
    if (ct == 0) {
      hPurityEta->SetBinContent(iBin, 0.);
      hPurityEta->SetBinError(iBin, 0.);
    } else {
      double p = cg / ct;
      double ep = std::sqrt(p * (1 - p) / ct);
      hPurityEta->SetBinContent(iBin, p);
      hPurityEta->SetBinError(iBin, ep);
    }
  }

  TCanvas* cef = new TCanvas("cef", "", 1400, 800);
  cef->Divide(2, 2);
  cef->cd(1);
  hMomGen->SetLineColor(kGray + 1);
  hMomGen->SetLineWidth(3);
  hMomGen->Draw();
  hMomReco->SetLineWidth(2);
  hMomReco->Draw("same");
  TLegend* leg = new TLegend(0.5, 0.6, 0.89, 0.8);
  leg->AddEntry(hMomGen, "Generated, 6 hits in MS");
  for (int jIteration = 0; jIteration < nIterationsCA; ++jIteration) {
    if (hMomRecoIterCA[jIteration]->GetEntries() > 0) {
      hMomRecoIterCA[jIteration]->SetLineColor(colors[jIteration]);
      hMomRecoIterCA[jIteration]->SetLineWidth(2);
      hMomRecoIterCA[jIteration]->Draw("same");
      leg->AddEntry(hMomRecoIterCA[jIteration], Form("Reconstructed Interation %d, 6 hits in MS", jIteration));
    }
  }
  leg->Draw();
  gPad->SetGrid();
  cef->cd(2);
  TH1F* hEffMom = (TH1F*)hMomReco->Clone("hEffMom");
  hEffMom->Divide(hMomReco, hMomGen, 1., 1., "B");
  hEffMom->GetYaxis()->SetTitle("Efficiency");
  hEffMom->SetStats(0);
  hEffMom->Draw();
  gPad->SetGrid();
  cef->cd(3);
  hEtaGen->SetLineColor(kGray + 1);
  hEtaGen->SetLineWidth(3);
  hEtaGen->Draw();
  hEtaReco->SetLineWidth(2);
  hEtaReco->Draw("same");
  for (int jIteration = 0; jIteration < nIterationsCA; ++jIteration) {
    if (hEtaRecoIterCA[jIteration]->GetEntries() > 0) {
      hEtaRecoIterCA[jIteration]->SetLineColor(colors[jIteration]);
      hEtaRecoIterCA[jIteration]->SetLineWidth(2);
      hEtaRecoIterCA[jIteration]->Draw("same");
    }
  }
  gPad->SetGrid();
  cef->cd(4);
  TH1F* hEffEta = (TH1F*)hEtaReco->Clone("hEffEta");
  hEffEta->Divide(hEtaReco, hEtaGen, 1., 1., "B");
  hEffEta->GetYaxis()->SetTitle("Efficiency");
  hEffEta->SetStats(0);
  hEffEta->Draw();
  gPad->SetGrid();
  cef->Print("msEff.png");

  TCanvas* cpu = new TCanvas("cpu", "", 1400, 800);
  cpu->Divide(2, 2);
  cpu->cd(1);
  hMomReco->Draw();
  hMomGoodReco->SetLineColor(kGreen + 1);
  ;
  hMomGoodReco->Draw("same");
  gPad->SetGrid();
  cpu->cd(2);
  hPurityMom->GetYaxis()->SetTitle("Purity");
  hPurityMom->SetMinimum(0);
  hPurityMom->SetStats(0);
  hPurityMom->Draw();
  gPad->SetGrid();
  cpu->cd(3);
  hEtaReco->Draw();
  hEtaGoodReco->SetLineColor(kGreen + 1);
  ;
  hEtaGoodReco->Draw("same");
  gPad->SetGrid();
  cpu->cd(4);
  hPurityEta->GetYaxis()->SetTitle("Purity");
  hPurityEta->SetMinimum(0);
  hPurityEta->SetStats(0);
  hPurityEta->Draw();
  gPad->SetGrid();
  cpu->Print("msPur.png");
}
