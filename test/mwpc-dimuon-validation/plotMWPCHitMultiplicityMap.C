#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <memory>
#include <set>
#include <string>
#include <vector>

#include <TCanvas.h>
#include <TFile.h>
#include <TH2D.h>
#include <TLine.h>
#include <TParticle.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

#include "NA6PLayoutParam.h"
#include "NA6PMWPCChamber.h"
#include "NA6PMWPCParam.h"
#include "NA6PMuonSpecHit.h"

namespace
{
constexpr int kNStations = 6;

int parentPDG(const std::string& channel)
{
  if (channel == "Jpsi") return 443;
  if (channel == "Omega") return 223;
  if (channel == "Phi") return 333;
  return 0;
}

std::string displayName(const std::string& channel)
{
  if (channel == "Jpsi") return "J/#psi";
  if (channel == "Omega") return "#omega";
  if (channel == "Phi") return "#phi";
  return channel;
}

int stationFromDetectorID(int detectorID, const NA6PMWPCParam& p)
{
  int first = 0;
  for (int station = 0; station < kNStations; ++station) {
    const int n = p.stationGridNX[station] * p.stationGridNY[station];
    if (detectorID >= first && detectorID < first + n) return station;
    first += n;
  }
  return -1;
}

std::array<double, 2> detectorGasSize(int station, const NA6PMWPCParam& p)
{
  const double gasX = p.bodyY - 2. * NA6PMWPCChamber::InnerFrameWidth;
  const double gasY = (station == 0 && p.useNarrowMS0)
                        ? p.ms0GasY
                        : p.bodyX - 2. * NA6PMWPCChamber::InnerFrameWidth;
  return {gasX, gasY};
}

void drawActiveGrid(int station, double xmin, double xmax, double ymin, double ymax,
                    std::vector<std::unique_ptr<TLine>>& lines)
{
  const auto& p = NA6PMWPCParam::Instance();
  const int nx = p.stationGridNX[station];
  const int ny = p.stationGridNY[station];
  const auto gas = detectorGasSize(station, p);
  const double gasX = gas[0];
  const double gasY = gas[1];
  const double pitchX = gasX - p.activeOverlapX;
  const double pitchY = gasY - p.activeOverlapY;

  auto addV = [&](double x) {
    auto l = std::make_unique<TLine>(x, ymin, x, ymax);
    l->SetLineColorAlpha(kGray + 2, 0.35);
    l->SetLineWidth(1);
    l->Draw("same");
    lines.emplace_back(std::move(l));
  };
  auto addH = [&](double y) {
    auto l = std::make_unique<TLine>(xmin, y, xmax, y);
    l->SetLineColorAlpha(kGray + 2, 0.35);
    l->SetLineWidth(1);
    l->Draw("same");
    lines.emplace_back(std::move(l));
  };

  for (int col = 0; col < nx; ++col) {
    const double xc = (static_cast<double>(col) - 0.5 * (nx - 1)) * pitchX;
    addV(xc - 0.5 * gasX);
    addV(xc + 0.5 * gasX);
  }
  for (int row = 0; row < ny; ++row) {
    const double yc = (static_cast<double>(row) - 0.5 * (ny - 1)) * pitchY;
    addH(yc - 0.5 * gasY);
    addH(yc + 0.5 * gasY);
  }
}
} // namespace

void plotMWPCHitMultiplicityMap(const char* runDir = "test_runs/mwpc_dimuon/Phi",
                                const char* channelName = "Phi",
                                double binSizeCm = 10.,
                                int minMuonEntries = 5)
{
  const std::string channel = channelName;
  const int wantedParentPDG = parentPDG(channel);
  if (!wantedParentPDG) {
    std::cerr << "Unknown channel '" << channel << "'. Use Jpsi, Omega or Phi.\n";
    return;
  }

  std::unique_ptr<TFile> kineFile(TFile::Open((std::string(runDir) + "/MCKine.root").c_str(), "READ"));
  std::unique_ptr<TFile> hitsFile(TFile::Open((std::string(runDir) + "/HitsMuonSpec.root").c_str(), "READ"));
  if (!kineFile || kineFile->IsZombie() || !hitsFile || hitsFile->IsZombie()) {
    std::cerr << "Cannot open input ROOT files in " << runDir << '\n';
    return;
  }

  auto* kineTree = dynamic_cast<TTree*>(kineFile->Get("mckine"));
  auto* hitsTree = dynamic_cast<TTree*>(hitsFile->Get("hitsMuonSpec"));
  if (!kineTree || !hitsTree || kineTree->GetEntries() != hitsTree->GetEntries()) {
    std::cerr << "Missing trees or event-count mismatch.\n";
    return;
  }

  std::vector<TParticle>* tracks = nullptr;
  std::vector<NA6PMuonSpecHit>* hits = nullptr;
  kineTree->SetBranchAddress("tracks", &tracks);
  hitsTree->SetBranchAddress("MuonSpec", &hits);

  const auto& p = NA6PMWPCParam::Instance();
  const auto& layout = NA6PLayoutParam::Instance();

  std::array<std::unique_ptr<TH2D>, kNStations> hHitSum;
  std::array<std::unique_ptr<TH2D>, kNStations> hMuonCount;
  std::array<std::unique_ptr<TH2D>, kNStations> hMeanHits;
  std::array<double, kNStations> xMin{}, xMax{}, yMin{}, yMax{};

  for (int station = 0; station < kNStations; ++station) {
    const int nx = p.stationGridNX[station];
    const int ny = p.stationGridNY[station];
    const auto gas = detectorGasSize(station, p);
    const double gasX = gas[0];
    const double gasY = gas[1];
    const double pitchX = gasX - p.activeOverlapX;
    const double pitchY = gasY - p.activeOverlapY;
    const double activeX = gasX + (nx - 1) * pitchX;
    const double activeY = gasY + (ny - 1) * pitchY;
    const double margin = 10.;

    xMin[station] = -0.5 * activeX - margin;
    xMax[station] = +0.5 * activeX + margin;
    yMin[station] = -0.5 * activeY - margin;
    yMax[station] = +0.5 * activeY + margin;

    const int nbinX = std::max(10, static_cast<int>(std::ceil((xMax[station] - xMin[station]) / binSizeCm)));
    const int nbinY = std::max(10, static_cast<int>(std::ceil((yMax[station] - yMin[station]) / binSizeCm)));

    hHitSum[station] = std::make_unique<TH2D>(Form("hHitSum_MS%d", station), "",
      nbinX, xMin[station], xMax[station], nbinY, yMin[station], yMax[station]);
    hMuonCount[station] = std::make_unique<TH2D>(Form("hMuonCount_MS%d", station), "",
      nbinX, xMin[station], xMax[station], nbinY, yMin[station], yMax[station]);
    hHitSum[station]->SetDirectory(nullptr);
    hMuonCount[station]->SetDirectory(nullptr);
  }

  long nGoodPairs = 0;
  long nDaughterMuons = 0;

  const Long64_t nEvents = kineTree->GetEntries();
  for (Long64_t iev = 0; iev < nEvents; ++iev) {
    kineTree->GetEntry(iev);
    hitsTree->GetEntry(iev);
    if (!tracks || !hits) continue;

    int parentIndex = -1;
    for (std::size_t i = 0; i < tracks->size(); ++i) {
      const auto& tr = (*tracks)[i];
      if (tr.GetPdgCode() == wantedParentPDG && tr.GetFirstMother() < 0) {
        parentIndex = static_cast<int>(i);
        break;
      }
    }
    if (parentIndex < 0) continue;

    std::vector<int> muonIndices;
    for (std::size_t i = 0; i < tracks->size(); ++i) {
      const auto& tr = (*tracks)[i];
      if (std::abs(tr.GetPdgCode()) == 13 && tr.GetFirstMother() == parentIndex) {
        muonIndices.push_back(static_cast<int>(i));
      }
    }
    if (muonIndices.size() != 2) continue;
    ++nGoodPairs;

    for (int muonIndex : muonIndices) {
      ++nDaughterMuons;
      std::array<std::set<int>, kNStations> binsTouched;

      for (const auto& hit : *hits) {
        if (hit.getTrackID() != muonIndex) continue;
        const int station = stationFromDetectorID(hit.getDetectorID(), p);
        if (station < 0 || station >= kNStations) continue;

        const double xCentre = layout.shiftMS[0] + layout.posMSPlaneX[station];
        const double yCentre = layout.shiftMS[1] + layout.posMSPlaneY[station];
        const double x = hit.getXIn() - xCentre;
        const double y = hit.getYIn() - yCentre;

        hHitSum[station]->Fill(x, y);
        const int globalBin = hMuonCount[station]->FindBin(x, y);
        binsTouched[station].insert(globalBin);
      }

      for (int station = 0; station < kNStations; ++station) {
        for (int bin : binsTouched[station]) {
          hMuonCount[station]->AddBinContent(bin, 1.0);
        }
      }
    }
  }

  const std::string plotDir = std::string(runDir) + "/plots_hit_multiplicity";
  gSystem->mkdir(plotDir.c_str(), true);
  gStyle->SetOptStat(0);
  gStyle->SetNumberContours(100);

  for (int station = 0; station < kNStations; ++station) {
    hMeanHits[station].reset(static_cast<TH2D*>(hHitSum[station]->Clone(Form("hMeanHits_MS%d", station))));
    hMeanHits[station]->Reset("ICES");
    hMeanHits[station]->SetDirectory(nullptr);

    for (int ix = 1; ix <= hMeanHits[station]->GetNbinsX(); ++ix) {
      for (int iy = 1; iy <= hMeanHits[station]->GetNbinsY(); ++iy) {
        const int bin = hMeanHits[station]->GetBin(ix, iy);
        const double nMu = hMuonCount[station]->GetBinContent(bin);
        if (nMu >= minMuonEntries) {
          hMeanHits[station]->SetBinContent(bin, hHitSum[station]->GetBinContent(bin) / nMu);
        } else {
          hMeanHits[station]->SetBinContent(bin, -1.);
        }
      }
    }

    hMeanHits[station]->SetTitle(Form(
      "MS%d - %s - direct Geant4 chamber-hit multiplicity;"
      "detector X at sensitive-gas crossing [cm];detector Y at sensitive-gas crossing [cm];"
      "mean chamber hits per contributing muon",
      station, displayName(channel).c_str()));
    hMeanHits[station]->SetMinimum(0.);
    hMeanHits[station]->SetMaximum(4.);

    TCanvas c(Form("cMeanHits_MS%d", station), "", 1100, 900);
    c.SetRightMargin(0.20);
    hMeanHits[station]->Draw("COLZ");
    std::vector<std::unique_ptr<TLine>> gridLines;
    drawActiveGrid(station, xMin[station], xMax[station], yMin[station], yMax[station], gridLines);
    c.SaveAs(Form("%s/mean_hits_MS%d_%s.png", plotDir.c_str(), station, channel.c_str()));
  }

  TFile out(Form("%s/hit_multiplicity_%s.root", plotDir.c_str(), channel.c_str()), "RECREATE");
  for (int station = 0; station < kNStations; ++station) {
    hHitSum[station]->Write();
    hMuonCount[station]->Write();
    hMeanHits[station]->Write();
  }
  out.Close();

  std::cout << "\n========== MWPC HIT MULTIPLICITY MAP: " << channel << " ==========\n";
  std::cout << "events            = " << nEvents << '\n';
  std::cout << "good dimuon pairs = " << nGoodPairs << '\n';
  std::cout << "daughter muons    = " << nDaughterMuons << '\n';
  std::cout << "plots: " << plotDir << '\n';
  std::cout << "========================================================\n";
}
