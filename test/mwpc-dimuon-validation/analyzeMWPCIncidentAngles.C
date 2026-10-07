#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TParticle.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

#include "NA6PLayoutParam.h"
#include "NA6PMWPCParam.h"
#include "NA6PMuonSpecHit.h"

namespace
{
constexpr int kNStations = 6;
constexpr int kNStudyStations = 4;
constexpr int kNSelections = 3;
constexpr double kRadToDeg = 57.2957795130823208768;

// Current analytic NA6PDipoleMS geometry.
constexpr double kMNP33LengthCm = 130.;
constexpr double kGeometryHalfXcm = 160.; // 320 cm analytic aperture
constexpr double kGeometryHalfYcm = 120.; // 240 cm analytic aperture

// Conservative horizontal fiducial from the quoted MNP33 useful aperture
// (245 cm); the vertical useful opening is kept at 240 cm.
constexpr double kUsefulHalfXcm = 122.5;
constexpr double kUsefulHalfYcm = 120.;

enum Selection : int {
  kAll = 0,
  kGeometryFiducial = 1,
  kUsefulFiducial = 2
};

const char* selectionName(int selection)
{
  switch (selection) {
    case kAll: return "all";
    case kGeometryFiducial: return "geometry_fiducial";
    case kUsefulFiducial: return "useful_fiducial";
    default: return "unknown";
  }
}

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

struct AnglePoint {
  double thetaX = 0.;
  double thetaY = 0.;
  double theta = 0.;
  double p = 0.;
};

AnglePoint anglePoint(const NA6PMuonSpecHit& hit)
{
  const double px = hit.getPXIn();
  const double py = hit.getPYIn();
  const double pz = hit.getPZIn();
  AnglePoint a;
  a.thetaX = std::atan2(px, pz) * kRadToDeg;
  a.thetaY = std::atan2(py, pz) * kRadToDeg;
  a.theta = std::atan2(std::hypot(px, py), pz) * kRadToDeg;
  a.p = std::sqrt(px * px + py * py + pz * pz);
  return a;
}

struct Samples {
  std::vector<double> thetaX;
  std::vector<double> thetaY;
  std::vector<double> theta;
  std::vector<double> p;

  void fill(const NA6PMuonSpecHit& hit)
  {
    const auto a = anglePoint(hit);
    thetaX.push_back(a.thetaX);
    thetaY.push_back(a.thetaY);
    theta.push_back(a.theta);
    p.push_back(a.p);
  }

  std::size_t size() const { return theta.size(); }
};

struct Stats {
  std::size_t n = 0;
  double mean = 0.;
  double rms = 0.;
  double q05 = 0.;
  double q16 = 0.;
  double q50 = 0.;
  double q84 = 0.;
  double q95 = 0.;
  double q99 = 0.;
};

double quantileSorted(const std::vector<double>& sorted, double q)
{
  if (sorted.empty()) return 0.;
  if (sorted.size() == 1) return sorted.front();

  const double pos = q * static_cast<double>(sorted.size() - 1);
  const auto lo = static_cast<std::size_t>(std::floor(pos));
  const auto hi = static_cast<std::size_t>(std::ceil(pos));
  const double f = pos - static_cast<double>(lo);
  return sorted[lo] * (1. - f) + sorted[hi] * f;
}

Stats calculateStats(const std::vector<double>& values, bool absoluteValue)
{
  Stats s;
  s.n = values.size();
  if (values.empty()) return s;

  std::vector<double> work;
  work.reserve(values.size());

  double sum = 0.;
  double sum2 = 0.;
  for (double v : values) {
    const double x = absoluteValue ? std::abs(v) : v;
    work.push_back(x);
    sum += x;
    sum2 += x * x;
  }

  s.mean = sum / static_cast<double>(s.n);
  const double mean2 = sum2 / static_cast<double>(s.n);
  s.rms = std::sqrt(std::max(0., mean2 - s.mean * s.mean));

  std::sort(work.begin(), work.end());
  s.q05 = quantileSorted(work, 0.05);
  s.q16 = quantileSorted(work, 0.16);
  s.q50 = quantileSorted(work, 0.50);
  s.q84 = quantileSorted(work, 0.84);
  s.q95 = quantileSorted(work, 0.95);
  s.q99 = quantileSorted(work, 0.99);
  return s;
}

double maxAbs(const std::vector<double>& values)
{
  double m = 0.;
  for (double v : values) m = std::max(m, std::abs(v));
  return m;
}

double maxValue(const std::vector<double>& values)
{
  if (values.empty()) return 0.;
  return *std::max_element(values.begin(), values.end());
}

const NA6PMuonSpecHit* representativeHit(
  const std::vector<const NA6PMuonSpecHit*>& stationHits,
  double nominalZ)
{
  if (stationHits.empty()) return nullptr;
  return *std::min_element(
    stationHits.begin(), stationHits.end(),
    [nominalZ](const auto* a, const auto* b) {
      return std::abs(a->getZIn() - nominalZ) < std::abs(b->getZIn() - nominalZ);
    });
}

bool extrapolateToZ(const NA6PMuonSpecHit& hit, double zTarget,
                    double& x, double& y)
{
  const double pz = hit.getPZIn();
  if (std::abs(pz) < 1.e-12) return false;
  const double dz = zTarget - hit.getZIn();
  x = hit.getXIn() + (hit.getPXIn() / pz) * dz;
  y = hit.getYIn() + (hit.getPYIn() / pz) * dz;
  return true;
}

struct FiducialResult {
  bool geometry = false;
  bool useful = false;
  double xEntrance = 0.;
  double yEntrance = 0.;
  double xExit = 0.;
  double yExit = 0.;
};

FiducialResult classifyFiducial(const NA6PMuonSpecHit& ms1,
                                const NA6PMuonSpecHit& ms2,
                                const NA6PLayoutParam& layout)
{
  FiducialResult f;

  const double zEntrance = layout.posDipMS[2] - 0.5 * kMNP33LengthCm;
  const double zExit = layout.posDipMS[2] + 0.5 * kMNP33LengthCm;

  if (!extrapolateToZ(ms1, zEntrance, f.xEntrance, f.yEntrance)) return f;
  if (!extrapolateToZ(ms2, zExit, f.xExit, f.yExit)) return f;

  // Compare to the aperture about the actual configured magnet centre.
  f.xEntrance -= layout.posDipMS[0];
  f.yEntrance -= layout.posDipMS[1];
  f.xExit -= layout.posDipMS[0];
  f.yExit -= layout.posDipMS[1];

  f.geometry =
    std::abs(f.xEntrance) < kGeometryHalfXcm &&
    std::abs(f.yEntrance) < kGeometryHalfYcm &&
    std::abs(f.xExit) < kGeometryHalfXcm &&
    std::abs(f.yExit) < kGeometryHalfYcm;

  f.useful =
    std::abs(f.xEntrance) < kUsefulHalfXcm &&
    std::abs(f.yEntrance) < kUsefulHalfYcm &&
    std::abs(f.xExit) < kUsefulHalfXcm &&
    std::abs(f.yExit) < kUsefulHalfYcm;

  return f;
}

void writeStatsRow(std::ofstream& csv, const std::string& channel,
                   int station, const char* selection,
                   const char* variable, const Stats& s)
{
  csv << channel << ',' << station << ',' << selection << ',' << variable << ','
      << s.n << ','
      << s.mean << ',' << s.rms << ','
      << s.q05 << ',' << s.q16 << ',' << s.q50 << ',' << s.q84 << ','
      << s.q95 << ',' << s.q99 << '\n';
}

void fillHistogram(TH1D& h, const std::vector<double>& values)
{
  for (double v : values) h.Fill(v);
}

void writePlots(const Samples& s, int station, int selection,
                const std::string& weighting, const std::string& channel,
                const std::string& plotDir, TFile& out)
{
  const std::string tag =
    weighting + "_MS" + std::to_string(station) + "_" +
    selectionName(selection) + "_" + channel;
  const std::string titlePrefix =
    "MS" + std::to_string(station) + " - " + displayName(channel) +
    " - " + selectionName(selection) + " - " + weighting;

  const double signedMax = std::max(5., 1.08 * std::max(maxAbs(s.thetaX), maxAbs(s.thetaY)));
  const double thetaMax = std::max(5., 1.08 * maxValue(s.theta));
  const double pMax = std::max(1., 1.08 * maxValue(s.p));

  auto hThetaX = std::make_unique<TH1D>(
    ("hThetaX_" + tag).c_str(),
    (titlePrefix + ";#theta_{x} [deg];entries").c_str(),
    120, -signedMax, signedMax);
  auto hThetaY = std::make_unique<TH1D>(
    ("hThetaY_" + tag).c_str(),
    (titlePrefix + ";#theta_{y} [deg];entries").c_str(),
    120, -signedMax, signedMax);
  auto hTheta = std::make_unique<TH1D>(
    ("hTheta_" + tag).c_str(),
    (titlePrefix + ";#theta [deg];entries").c_str(),
    120, 0., thetaMax);
  auto hAbsThetaX = std::make_unique<TH1D>(
    ("hAbsThetaX_" + tag).c_str(),
    (titlePrefix + ";|#theta_{x}| [deg];entries").c_str(),
    120, 0., signedMax);
  auto hAbsThetaY = std::make_unique<TH1D>(
    ("hAbsThetaY_" + tag).c_str(),
    (titlePrefix + ";|#theta_{y}| [deg];entries").c_str(),
    120, 0., signedMax);
  auto hThetaXY = std::make_unique<TH2D>(
    ("hThetaXY_" + tag).c_str(),
    (titlePrefix + ";#theta_{x} [deg];#theta_{y} [deg]").c_str(),
    100, -signedMax, signedMax, 100, -signedMax, signedMax);
  auto hPTheta = std::make_unique<TH2D>(
    ("hPTheta_" + tag).c_str(),
    (titlePrefix + ";p [GeV/c];#theta [deg]").c_str(),
    100, 0., pMax, 100, 0., thetaMax);

  hThetaX->SetDirectory(nullptr);
  hThetaY->SetDirectory(nullptr);
  hTheta->SetDirectory(nullptr);
  hAbsThetaX->SetDirectory(nullptr);
  hAbsThetaY->SetDirectory(nullptr);
  hThetaXY->SetDirectory(nullptr);
  hPTheta->SetDirectory(nullptr);

  fillHistogram(*hThetaX, s.thetaX);
  fillHistogram(*hThetaY, s.thetaY);
  fillHistogram(*hTheta, s.theta);
  for (double v : s.thetaX) hAbsThetaX->Fill(std::abs(v));
  for (double v : s.thetaY) hAbsThetaY->Fill(std::abs(v));
  for (std::size_t i = 0; i < s.size(); ++i) {
    hThetaXY->Fill(s.thetaX[i], s.thetaY[i]);
    hPTheta->Fill(s.p[i], s.theta[i]);
  }

  out.cd();
  hThetaX->Write();
  hThetaY->Write();
  hTheta->Write();
  hAbsThetaX->Write();
  hAbsThetaY->Write();
  hThetaXY->Write();
  hPTheta->Write();

  TCanvas c(("c_" + tag).c_str(), "", 1800, 950);
  c.Divide(4, 2);
  c.cd(1); hThetaX->Draw("HIST");
  c.cd(2); hThetaY->Draw("HIST");
  c.cd(3); hTheta->Draw("HIST");
  c.cd(4); hAbsThetaX->Draw("HIST");
  c.cd(5); hAbsThetaY->Draw("HIST");
  c.cd(6); hThetaXY->Draw("COLZ");
  c.cd(7); hPTheta->Draw("COLZ");
  c.cd(8);
  auto hCount = std::make_unique<TH1D>(
    ("hCount_" + tag).c_str(),
    Form("%s;sample;entries", titlePrefix.c_str()), 1, 0., 1.);
  hCount->SetDirectory(nullptr);
  hCount->SetBinContent(1, static_cast<double>(s.size()));
  hCount->GetXaxis()->SetBinLabel(1, weighting.c_str());
  hCount->Draw("HIST TEXT0");
  out.cd();
  hCount->Write();

  c.SaveAs(Form("%s/angles_%s.png", plotDir.c_str(), tag.c_str()));
}

} // namespace

void analyzeMWPCIncidentAngles(
  const char* runDir = "test_runs/mwpc_dimuon/Jpsi",
  const char* channelName = "Jpsi")
{
  const std::string channel = channelName;
  const int wantedParentPDG = parentPDG(channel);
  if (!wantedParentPDG) {
    std::cerr << "Unknown channel '" << channel << "'. Use Jpsi, Omega or Phi.\n";
    return;
  }

  const std::string kineName = std::string(runDir) + "/MCKine.root";
  const std::string hitsName = std::string(runDir) + "/HitsMuonSpec.root";

  std::unique_ptr<TFile> kineFile(TFile::Open(kineName.c_str(), "READ"));
  std::unique_ptr<TFile> hitsFile(TFile::Open(hitsName.c_str(), "READ"));
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

  const auto& mwpc = NA6PMWPCParam::Instance();
  const auto& layout = NA6PLayoutParam::Instance();

  using StationSelectionSamples =
    std::array<std::array<Samples, kNSelections>, kNStudyStations>;
  StationSelectionSamples representativeSamples;
  StationSelectionSamples crossingSamples;

  long nGoodPairs = 0;
  long nDaughterMuons = 0;
  long nWithMS1 = 0;
  long nWithMS2 = 0;
  long nWithBoth = 0;
  long nGeometryFiducial = 0;
  long nUsefulFiducial = 0;

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

      std::array<std::vector<const NA6PMuonSpecHit*>, kNStations> stationHits;
      for (const auto& hit : *hits) {
        if (hit.getTrackID() != muonIndex) continue;
        const int station = stationFromDetectorID(hit.getDetectorID(), mwpc);
        if (station < 0 || station >= kNStations) continue;
        stationHits[station].push_back(&hit);
      }

      std::array<const NA6PMuonSpecHit*, kNStations> representative{};
      for (int station = 0; station < kNStations; ++station) {
        const double nominalZ =
          layout.shiftMS[2] + layout.posMSPlaneZ[station];
        representative[station] =
          representativeHit(stationHits[station], nominalZ);
      }

      if (representative[1]) ++nWithMS1;
      if (representative[2]) ++nWithMS2;

      FiducialResult fid;
      if (representative[1] && representative[2]) {
        ++nWithBoth;
        fid = classifyFiducial(*representative[1], *representative[2], layout);
        if (fid.geometry) ++nGeometryFiducial;
        if (fid.useful) ++nUsefulFiducial;
      }

      const std::array<bool, kNSelections> selected = {
        true, fid.geometry, fid.useful
      };

      for (int station = 0; station < kNStudyStations; ++station) {
        if (representative[station]) {
          for (int sel = 0; sel < kNSelections; ++sel) {
            if (selected[sel]) {
              representativeSamples[station][sel].fill(*representative[station]);
            }
          }
        }

        for (const auto* hit : stationHits[station]) {
          for (int sel = 0; sel < kNSelections; ++sel) {
            if (selected[sel]) crossingSamples[station][sel].fill(*hit);
          }
        }
      }
    }
  }

  const std::string plotDir = std::string(runDir) + "/plots_angles";
  gSystem->mkdir(plotDir.c_str(), true);
  gStyle->SetOptStat(0);

  const std::string rootOutName =
    plotDir + "/incident_angles_" + channel + ".root";
  TFile out(rootOutName.c_str(), "RECREATE");

  for (int station = 0; station < kNStudyStations; ++station) {
    for (int sel = 0; sel < kNSelections; ++sel) {
      writePlots(representativeSamples[station][sel], station, sel,
                 "muon", channel, plotDir, out);
      writePlots(crossingSamples[station][sel], station, sel,
                 "crossing", channel, plotDir, out);
    }
  }
  out.Close();

  const std::string csvName =
    plotDir + "/incident_angle_summary_" + channel + ".csv";
  std::ofstream csv(csvName);
  if (!csv) {
    std::cerr << "Cannot create " << csvName << '\n';
    return;
  }

  const double zEntrance = layout.posDipMS[2] - 0.5 * kMNP33LengthCm;
  const double zExit = layout.posDipMS[2] + 0.5 * kMNP33LengthCm;

  csv << std::setprecision(10);
  csv << "# Stage C0 MWPC incidence-angle study\n";
  csv << "# summary weighting,one representative hit per muon per station\n";
  csv << "# representative hit,minimum |zIn - nominal station z|\n";
  csv << "# MNP33 centre z cm," << layout.posDipMS[2] << '\n';
  csv << "# MNP33 entrance z cm," << zEntrance << '\n';
  csv << "# MNP33 exit z cm," << zExit << '\n';
  csv << "# geometry aperture half-widths cm," << kGeometryHalfXcm
      << ',' << kGeometryHalfYcm << '\n';
  csv << "# useful aperture half-widths cm," << kUsefulHalfXcm
      << ',' << kUsefulHalfYcm << '\n';
  csv << "# generated daughter muons," << nDaughterMuons << '\n';
  csv << "# with MS1 hit," << nWithMS1 << '\n';
  csv << "# with MS2 hit," << nWithMS2 << '\n';
  csv << "# with both MS1 and MS2," << nWithBoth << '\n';
  csv << "# geometry fiducial," << nGeometryFiducial << '\n';
  csv << "# useful fiducial," << nUsefulFiducial << '\n';
  csv << "channel,station,selection,variable,N,mean_deg,rms_deg,"
         "q05_deg,q16_deg,q50_deg,q84_deg,q95_deg,q99_deg\n";

  std::cout << "\n========== MWPC INCIDENT ANGLES: " << channel << " ==========\n";
  std::cout << "events                   = " << nEvents << '\n';
  std::cout << "good dimuon pairs        = " << nGoodPairs << '\n';
  std::cout << "daughter muons           = " << nDaughterMuons << '\n';
  std::cout << "with MS1 hit             = " << nWithMS1 << '\n';
  std::cout << "with MS2 hit             = " << nWithMS2 << '\n';
  std::cout << "with both MS1 and MS2    = " << nWithBoth << '\n';
  std::cout << "geometry fiducial        = " << nGeometryFiducial << '\n';
  std::cout << "useful fiducial          = " << nUsefulFiducial << '\n';
  std::cout << "MNP33 faces z            = " << zEntrance
            << ", " << zExit << " cm\n";
  std::cout << "\nRepresentative-hit quantiles [deg]:\n";
  std::cout << "station  selection             N       q95|tx|   q99|tx|"
               "   q95|ty|   q99|ty|   q95(theta) q99(theta)\n";

  for (int station = 0; station < kNStudyStations; ++station) {
    for (int sel = 0; sel < kNSelections; ++sel) {
      const auto& s = representativeSamples[station][sel];
      const auto sx = calculateStats(s.thetaX, true);
      const auto sy = calculateStats(s.thetaY, true);
      const auto st = calculateStats(s.theta, false);

      writeStatsRow(csv, channel, station, selectionName(sel),
                    "abs_thetaX", sx);
      writeStatsRow(csv, channel, station, selectionName(sel),
                    "abs_thetaY", sy);
      writeStatsRow(csv, channel, station, selectionName(sel),
                    "theta", st);

      std::cout << std::setw(7) << ("MS" + std::to_string(station)) << "  "
                << std::setw(21) << selectionName(sel) << "  "
                << std::setw(7) << s.size() << "  "
                << std::setw(9) << std::fixed << std::setprecision(3) << sx.q95
                << std::setw(10) << sx.q99
                << std::setw(10) << sy.q95
                << std::setw(10) << sy.q99
                << std::setw(12) << st.q95
                << std::setw(11) << st.q99 << '\n';
    }
  }

  csv.close();

  std::cout << "\nROOT output: " << rootOutName << '\n';
  std::cout << "CSV summary: " << csvName << '\n';
  std::cout << "PNG plots:   " << plotDir << "/angles_*.png\n";
  std::cout << "=========================================================\n";
}
