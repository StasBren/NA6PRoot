// Reconstruct the old generated J/psi, phi and omega momentum/rapidity plots
// from the 6 saved high-statistics MNP33 bore production MCKine.root files.
//
// IMPORTANT: this is GENERATED parent/daughter truth, not accepted-muon
// kinematics after absorber / magnet or the current checkerboard projection.
// No new simulation. A smoke run can cap events per chunk.
//
// Example (normal users run the accompanying bash wrapper):
// root -l -b -q 'test/mwpc-stage-c/recover_bore_dimuon_kinematics.C+("...", "...", 200)'

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TParticle.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>
#include <TVirtualPad.h>

namespace BoreKinematics
{
struct Channel {
  const char* channel;
  const char* directory;
  int parentPDG;
  const char* title;
};

double momentum(const TParticle& p)
{
  return std::hypot(std::hypot(p.Px(), p.Py()), p.Pz());
}

double pt(const TParticle& p)
{
  return std::hypot(p.Px(), p.Py());
}

double rapidity(const TParticle& p)
{
  // Robust for physically forward-going particles; invalid values stay NaN,
  // never silently mapped onto the y=0 bin.
  const double num = p.Energy() + p.Pz();
  const double den = p.Energy() - p.Pz();
  if (num <= 0. || den <= 0.) return std::numeric_limits<double>::quiet_NaN();
  return .5 * std::log(num / den);
}

struct Histograms {
  TH1D parentP;
  TH1D parentPt;
  TH1D parentY;
  TH1D daughterP;
  TH1D daughterPt;
  TH1D daughterY;
  TH1D daughterThetaU;
  TH1D daughterThetaW;

  explicit Histograms(const std::string& tag)
    : parentP((tag + "_parentP").c_str(), "Parent |p|;|p| [GeV/c];parents", 160, 0., 100.),
      parentPt((tag + "_parentPt").c_str(), "Parent p_{T};p_{T} [GeV/c];parents", 120, 0., 10.),
      parentY((tag + "_parentY").c_str(), "Parent rapidity;rapidity y;parents", 120, 0., 6.),
      daughterP((tag + "_daughterP").c_str(), "Muon |p| at production;|p| [GeV/c];muons", 160, 0., 100.),
      daughterPt((tag + "_daughterPt").c_str(), "Muon p_{T} at production;p_{T} [GeV/c];muons", 120, 0., 10.),
      daughterY((tag + "_daughterY").c_str(), "Muon rapidity at production;rapidity y;muons", 160, -1., 7.),
      daughterThetaU((tag + "_daughterThetaU").c_str(),
                      "Muon #theta_{u} at production;#theta_{u} [deg];muons",
                      240, -60., 60.),
      daughterThetaW((tag + "_daughterThetaW").c_str(),
                      "Muon #theta_{w} at production;#theta_{w} [deg];muons",
                      240, -60., 60.)
  {
    for (TH1D* h : all()) h->SetDirectory(nullptr);
  }

  std::vector<TH1D*> all()
  {
    return {&parentP, &parentPt, &parentY, &daughterP, &daughterPt,
            &daughterY, &daughterThetaU, &daughterThetaW};
  }
};

struct Counts {
  Long64_t events = 0;
  Long64_t parents = 0;
  Long64_t missingParent = 0;
  Long64_t validPairs = 0;
  Long64_t missingOrExtraMuons = 0;
  Long64_t daughters = 0;
  double maxMomentumClosure = 0.;
  double maxEnergyClosure = 0.;
};

void processChunk(const std::string& fileName, int parentPDG,
                  Long64_t eventLimit, Histograms& h, Counts& cnt)
{
  std::unique_ptr<TFile> file(TFile::Open(fileName.c_str(), "READ"));
  if (!file || file->IsZombie()) {
    throw std::runtime_error("Cannot open input MCKine.root: " + fileName);
  }
  auto* tree = dynamic_cast<TTree*>(file->Get("mckine"));
  if (!tree || !tree->GetBranch("tracks")) {
    throw std::runtime_error("Missing mckine tree / tracks branch: " + fileName);
  }
  std::vector<TParticle>* tracks = nullptr;
  if (tree->SetBranchAddress("tracks", &tracks) < 0) {
    throw std::runtime_error("Cannot bind mckine/tracks branch: " + fileName);
  }
  Long64_t n = tree->GetEntries();
  if (eventLimit > 0) n = std::min(n, eventLimit);
  for (Long64_t event = 0; event < n; ++event) {
    if (tree->GetEntry(event) <= 0 || !tracks) {
      throw std::runtime_error(
          "Cannot read event " + std::to_string(event) + " in " + fileName);
    }
    ++cnt.events;
    int parentIndex = -1;
    for (size_t i = 0; i < tracks->size(); ++i) {
      const auto& track = (*tracks)[i];
      if (track.GetPdgCode() == parentPDG && track.GetFirstMother() < 0) {
        parentIndex = static_cast<int>(i);
        break;
      }
    }
    if (parentIndex < 0) {
      ++cnt.missingParent;
      continue;
    }
    const auto& parent = (*tracks)[parentIndex];
    ++cnt.parents;
    h.parentP.Fill(momentum(parent));
    h.parentPt.Fill(pt(parent));
    if (std::isfinite(rapidity(parent))) h.parentY.Fill(rapidity(parent));

    std::vector<int> muonIds;
    for (size_t i = 0; i < tracks->size(); ++i) {
      const auto& track = (*tracks)[i];
      if (std::abs(track.GetPdgCode()) == 13 &&
          track.GetFirstMother() == parentIndex) {
        muonIds.push_back(static_cast<int>(i));
      }
    }
    if (muonIds.size() != 2) {
      ++cnt.missingOrExtraMuons;
      continue;
    }
    ++cnt.validPairs;

    double px = 0., py = 0., pz = 0., energy = 0.;
    for (int index : muonIds) {
      const auto& mu = (*tracks)[index];
      ++cnt.daughters;
      h.daughterP.Fill(momentum(mu));
      h.daughterPt.Fill(pt(mu));
      if (std::isfinite(rapidity(mu))) h.daughterY.Fill(rapidity(mu));
      // Same angular convention used in Stage-C muon station studies.
      if (mu.Pz() > 0.) {
        h.daughterThetaU.Fill(std::atan2(mu.Py(), mu.Pz()) * (180. / std::acos(-1.)));
        h.daughterThetaW.Fill(std::atan2(mu.Px(), mu.Pz()) * (180. / std::acos(-1.)));
      }
      px += mu.Px();
      py += mu.Py();
      pz += mu.Pz();
      energy += mu.Energy();
    }
    cnt.maxMomentumClosure = std::max(
      cnt.maxMomentumClosure,
      std::hypot(std::hypot(px - parent.Px(), py - parent.Py()),
                 pz - parent.Pz()));
    cnt.maxEnergyClosure = std::max(
      cnt.maxEnergyClosure, std::abs(energy - parent.Energy()));
  }
  tree->ResetBranchAddresses();
  std::cout << "  " << fileName << "  processed events=" << n << "\n";
}

void drawOne(const Channel& channel, Histograms& h, const Counts& c,
             const std::string& output)
{
  gStyle->SetOptStat(0);
  TCanvas canvas(("c_" + std::string(channel.channel)).c_str(), "",
                 1500, 1250);
  canvas.Divide(3, 3, .004, .004);
  auto all = h.all();
  for (size_t i = 0; i < all.size(); ++i) {
    canvas.cd(static_cast<int>(i + 1));
    gPad->SetTicks(1, 1);
    if (i == 0 || i == 3) {
      gPad->SetLogy();
      all[i]->SetMinimum(.5);
    }
    all[i]->SetLineWidth(2);
    all[i]->SetLineColor(kBlue + 1);
    all[i]->Draw("HIST");
  }
  canvas.cd(9);
  TLatex label;
  label.SetNDC();
  label.SetTextSize(.056);
  double y = .91;
  for (const auto& s : std::vector<std::string>{
      std::string(channel.title) + ", generated truth",
      "No chamber/bore selection",
      "Events: " + std::to_string(c.events),
      "Parents: " + std::to_string(c.parents),
      "Valid #mu^{+}#mu^{-} pairs: " + std::to_string(c.validPairs),
      "Missing parent: " + std::to_string(c.missingParent),
      "Invalid pair: " + std::to_string(c.missingOrExtraMuons),
      "Compare with Stage-C MS0 angles",
    }) {
    label.DrawLatex(.07, y, s.c_str());
    y -= .108;
  }
  canvas.SaveAs(output.c_str());
}

void checkOverflow(Histograms& h, const Channel& channel)
{
  // Warn about hidden axis overflow rather than silently plotting truncated
  // tails. No physics result is declared verified when substantial overflow
  // occurs.
  const auto all = h.all();
  for (auto* hist : all) {
    const double overflow = hist->GetBinContent(hist->GetNbinsX() + 1);
    const double underflow = hist->GetBinContent(0);
    const double allEntries = hist->GetEntries();
    const double fraction = allEntries > 0 ?
      (overflow + underflow) / allEntries : 0.;
    if (fraction > 0.001) {
      std::cerr << "WARNING: " << channel.channel << " "
                << hist->GetName() << " under/overflow fraction "
                << 100. * fraction << "%; check plotted axis limits\n";
    }
  }
}
} // namespace BoreKinematics

void recover_bore_dimuon_kinematics(
    const char* base = "",
    const char* out = "test_runs/mwpc_stage_c/recovered_dimuon_kinematics",
    Long64_t maxEventsPerChunk = 0)
{
  using namespace BoreKinematics;
  if (!base || !*base || !out || !*out || maxEventsPerChunk < 0) {
    throw std::runtime_error("Supply base directory, output directory, and "
                             "non-negative maxEventsPerChunk");
  }
  gSystem->mkdir(out, true);
  const std::vector<Channel> channels = {
    {"Jpsi", "jpsi40-bore", 443, "J/#psi"},
    {"Phi", "phi40-bore", 333, "#phi"},
    {"Omega", "omega40-bore", 223, "#omega"},
  };
  std::ofstream summary(std::string(out) + "/kinematics_truth_summary.csv");
  if (!summary) throw std::runtime_error("Cannot write kinematics summary");
  summary << "channel,events,parents,missing_parent,valid_pairs,"
             "invalid_pair_events,daughter_muons,max_momentum_closure_GeV,"
             "max_energy_closure_GeV\n";
  for (const auto& channel : channels) {
    std::cout << "=== " << channel.channel << " generated truth ===\n";
    Histograms h(channel.channel);
    Counts count;
    for (int chunk = 1; chunk <= 2; ++chunk) {
      TString filename = Form(
          "%s/%s-chunk%03d/MCKine.root",
          base, channel.directory, chunk);
      processChunk(filename.Data(), channel.parentPDG, maxEventsPerChunk,
                   h, count);
    }
    if (count.parents == 0) {
      throw std::runtime_error(
        std::string("No generated parent found for channel ") +
        channel.channel + " in existing files");
    }
    if (count.validPairs == 0) {
      std::cerr << "WARNING: no direct muon pairs retained for "
                << channel.channel << ". Parent distributions still valid; "
                << "daughter plots will be empty. Inspect MCKine genealogy.\n";
    }
    const auto plot = std::string(out) + "/kinematics_" +
                      channel.channel + ".png";
    drawOne(channel, h, count, plot);
    checkOverflow(h, channel);
    summary << channel.channel << ',' << count.events << ','
            << count.parents << ',' << count.missingParent << ','
            << count.validPairs << ','
            << count.missingOrExtraMuons << ',' << count.daughters << ','
            << count.maxMomentumClosure << ','
            << count.maxEnergyClosure << '\n';
    std::cout << "  parents=" << count.parents
              << ", muon pairs=" << count.validPairs
              << ", invalid pair events=" << count.missingOrExtraMuons
              << ", momentum closure max=" << count.maxMomentumClosure
              << " GeV\n";
  }
  std::cout << "Finished. Three PNG plots and one CSV summary in "
            << out << "\n";
}
