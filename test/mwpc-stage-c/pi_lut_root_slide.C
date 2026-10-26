// One-slide, publication-style ROOT figure from the measured 20-vs-500 CSV.
// Two vertically stacked TH1D histograms, identical scales and 1-sigma
// Garfield-LUT statistical error bars. Writes PDF, PNG and a ROOT file.
#include <TCanvas.h>
#include <TH1D.h>
#include <TFile.h>
#include <TStyle.h>
#include <TPad.h>
#include <TLatex.h>
#include <TColor.h>
#include <TString.h>
#include <TROOT.h>
#include <TAxis.h>
#include <TSystem.h>

#include <array>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

void pi_lut_root_slide(
    const char* csv_name,
    const char* output_dir,
    const char* geometry = "p3p0_field75p0um") {
  const std::string geom(geometry);
  if (geom != "p3p0_field75p0um" && geom != "p3p0_anode_only") {
    throw std::runtime_error("Unknown p=3 mm geometry");
  }
  std::ifstream input(csv_name);
  if (!input) throw std::runtime_error("Cannot open original measured CSV");
  std::string line;
  if (!std::getline(input, line) ||
      line.find("phi_lower_deg") == std::string::npos ||
      line.find("se500") == std::string::npos) {
    throw std::runtime_error("Unexpected 20-vs-500 CSV header");
  }
  TH1D old_lut("phi_lut20", "LUT20", 36, -180., 180.);
  TH1D new_lut("phi_lut500", "LUT500", 36, -180., 180.);
  std::array<bool, 36> used{};
  int n = 0;
  while (std::getline(input, line)) {
    if (line.empty()) continue;
    std::stringstream stream(line);
    std::vector<std::string> cells;
    std::string value;
    while (std::getline(stream, value, ',')) cells.push_back(value);
    if (cells.size() < 9) {
      throw std::runtime_error("Malformed data row");
    }
    if (cells[0] != geom) continue;
    const double lo = std::stod(cells[1]);
    const double hi = std::stod(cells[2]);
    const double p20 = std::stod(cells[3]);
    const double e20 = std::stod(cells[4]);
    const double p500 = std::stod(cells[5]);
    const double e500 = std::stod(cells[6]);
    const int bin = static_cast<int>(std::lround((lo + 180.) / 10.));
    if (bin < 0 || bin >= 36 || used[bin] ||
        std::abs(lo - (-180. + 10. * bin)) > 1.e-7 ||
        std::abs(hi - lo - 10.) > 1.e-7 ||
        !std::isfinite(p20) || !std::isfinite(p500) ||
        !std::isfinite(e20) || !std::isfinite(e500) ||
        p20 < 0. || p500 < 0. || e20 < 0. || e500 < 0.) {
      throw std::runtime_error("Invalid or repeated angular bin");
    }
    old_lut.SetBinContent(bin + 1, p20);
    old_lut.SetBinError(bin + 1, e20);
    new_lut.SetBinContent(bin + 1, p500);
    new_lut.SetBinError(bin + 1, e500);
    used[bin] = true;
    ++n;
  }
  if (n != 36 ||
      std::abs(old_lut.Integral() - 1.) > 1.e-4 ||
      std::abs(new_lut.Integral() - 1.) > 1.e-4) {
    throw std::runtime_error("Expected two normalized 36-bin PDFs");
  }

  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  gStyle->SetEndErrorSize(3);
  gStyle->SetTitleFont(42, "XYZ");
  gStyle->SetLabelFont(42, "XYZ");

  const double yMax = 1.16 * std::max(
    old_lut.GetMaximum() + 2. * old_lut.GetBinError(old_lut.GetMaximumBin()),
    new_lut.GetMaximum() + 2. * new_lut.GetBinError(new_lut.GetMaximumBin()));

  TCanvas canvas("c_lut_20_vs_500",
                 "Garfield LUT convergence | one slide", 1600, 900);
  canvas.SetFillColor(kWhite);

  TLatex header;
  header.SetNDC();
  header.SetTextFont(42);
  header.SetTextSize(.036);
  header.DrawLatex(.085, .945,
                   "Garfield LUT statistics: electron-hit azimuth");
  header.SetTextSize(.023);
  header.SetTextColor(kGray + 2);
  const char* config_label =
      geom == "p3p0_field75p0um" ? "field wire 75 #mum"
                                : "no field wire";
  header.DrawLatex(.086, .907,
                   TString::Format("J/#psi  |  MS0  |  pitch 3.0 mm  |  %s  |  B = 0",
                                   config_label));

  auto style_and_draw = [&](TH1D& hist, const char* padname,
                            double yLow, double yHigh,
                            int color, const char* label) {
    TPad* pad = new TPad(padname, "", .075, yLow, .945, yHigh);
    pad->SetFillColor(kWhite);
    pad->SetLeftMargin(.115);
    pad->SetRightMargin(.025);
    pad->SetTopMargin(.08);
    pad->SetBottomMargin(.17);
    pad->SetTicks(1, 1);
    pad->Draw();
    pad->cd();
    hist.SetTitle("");
    hist.SetLineColor(color);
    hist.SetMarkerColor(color);
    hist.SetMarkerStyle(20);
    hist.SetMarkerSize(.72);
    hist.SetLineWidth(2);
    hist.GetYaxis()->SetRangeUser(0., yMax);
    hist.GetXaxis()->SetTitle("#phi_{hit} [deg]");
    hist.GetYaxis()->SetTitle("Probability / 10^{#circ} bin");
    hist.GetXaxis()->SetTitleFont(42);
    hist.GetYaxis()->SetTitleFont(42);
    hist.GetXaxis()->SetLabelFont(42);
    hist.GetYaxis()->SetLabelFont(42);
    hist.GetXaxis()->SetTitleSize(.057);
    hist.GetYaxis()->SetTitleSize(.050);
    hist.GetXaxis()->SetLabelSize(.044);
    hist.GetYaxis()->SetLabelSize(.042);
    hist.GetXaxis()->SetTitleOffset(1.06);
    hist.GetYaxis()->SetTitleOffset(1.05);
    hist.GetXaxis()->SetNdivisions(509);
    hist.GetYaxis()->SetNdivisions(505);
    hist.Draw("HIST");
    hist.Draw("E1 SAME");
    TLatex section;
    section.SetNDC();
    section.SetTextFont(42);
    section.SetTextSize(.058);
    section.SetTextColor(kBlack);
    section.DrawLatex(.15, .89, label);
    canvas.cd();
  };

  style_and_draw(old_lut, "pad_old", .52, .865,
                 kGray + 3, "Original: 20 trajectories / spatial cell");
  style_and_draw(new_lut, "pad_new", .135, .48,
                 kAzure + 2, "Updated: 500 trajectories / spatial cell");

  TLatex footer;
  footer.SetNDC();
  footer.SetTextFont(42);
  footer.SetTextColor(kBlack);
  footer.SetTextSize(.026);
  footer.DrawLatex(.09, .088,
                   "#sigma_{LUT} #propto 1/#sqrt{K}:  20 #rightarrow 500 trajectories/cell  #Rightarrow  5#times lower expected statistical uncertainty");
  footer.SetTextSize(.019);
  footer.SetTextColor(kGray + 2);
  footer.DrawLatex(.09, .046,
                   "Same Geant4 muon crossings and analytically expected Heed spatial occupancy; errors = Garfield-LUT 1#sigma only");

  canvas.Modified();
  canvas.Update();

  const TString out(output_dir);
  gSystem->mkdir(output_dir, kTRUE);
  canvas.SaveAs(out + "/" + geom.c_str() + "_slide.png");
  canvas.SaveAs(out + "/" + geom.c_str() + "_slide.pdf");
  TFile file(out + "/" + geom.c_str() + "_slide.root", "RECREATE");
  if (file.IsZombie()) throw std::runtime_error("Failed to create ROOT file");
  old_lut.Write();
  new_lut.Write();
  canvas.Write();
  file.Close();
  std::cout << "PASS: ROOT TH1D upper/lower 20-vs-500 slide: "
            << output_dir << std::endl;
}
