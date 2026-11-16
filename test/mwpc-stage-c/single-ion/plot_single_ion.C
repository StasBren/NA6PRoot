// ROOT-only ATLAS-style visualization of Garfield++ single-ion results.
// Usage: root -l -b -q 'plot_single_ion.C("/absolute/output/path")'
#include <algorithm>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TMultiGraph.h"
#include "TPad.h"
#include "TStyle.h"

namespace {
struct Row {
  int step = 0, dtIndex = 0;
  double dt = 0., qMinus = 0., qPlus = 0.;
};
struct Points {
  std::vector<double> t, iMinus, iPlus, qMinus, qPlus;
};
std::vector<std::string> Split(const std::string& s) {
  std::vector<std::string> v;
  std::istringstream in(s);
  std::string x;
  while (std::getline(in,x,',')) v.push_back(x);
  return v;
}
std::vector<Row> ReadSummary(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("Cannot read "+path);
  std::string line;
  std::getline(in,line);
  auto headings=Split(line);
  const auto col=[&](const std::string& name) {
    const auto it=std::find(headings.begin(),headings.end(),name);
    if (it==headings.end()) throw std::runtime_error("Missing CSV column "+name);
    return std::size_t(it-headings.begin());
  };
  const auto kStep=col("step_index"), kIdx=col("dt_index"),
             kDt=col("dt_ns"), kMinus=col("minus_Q_ramo_fC"),
             kPlus=col("plus_Q_ramo_fC");
  std::vector<Row> rows;
  while (std::getline(in,line)) {
    if(line.empty())continue;
    auto t=Split(line);
    if (t.size()!=headings.size())throw std::runtime_error("Malformed summary row");
    rows.push_back({std::stoi(t[kStep]),std::stoi(t[kIdx]),std::stod(t[kDt]),
                    std::stod(t[kMinus]),std::stod(t[kPlus])});
  }
  if(rows.empty())throw std::runtime_error("Empty summary");
  return rows;
}
Points ReadWave(const std::string& path) {
  std::ifstream in(path);
  if(!in)throw std::runtime_error("Missing waveform "+path);
  std::string line;
  std::getline(in,line);
  Points out;
  // Reading every bin into ROOT graphics is wasteful. Analysis and current
  // integration were performed by Garfield/C++ on the full-resolution grid.
  // Here the display is downsampled ONLY for plotting.
  std::size_t j=0;
  while(std::getline(in,line)) {
    if(line.empty())continue;
    auto t=Split(line);
    if(t.size()!=5)throw std::runtime_error("Malformed waveform line "+path);
    if(j%40==0) {
      out.t.push_back(0.001*std::stod(t[0]));
      out.iMinus.push_back(std::stod(t[1]));
      out.iPlus.push_back(std::stod(t[2]));
      out.qMinus.push_back(std::stod(t[3]));
      out.qPlus.push_back(std::stod(t[4]));
    }
    ++j;
  }
  if(out.t.empty())throw std::runtime_error("Empty waveform "+path);
  return out;
}
TGraph* MakeGraph(const std::vector<double>& x,
                  const std::vector<double>& y, int color) {
  auto* g=new TGraph(static_cast<int>(x.size()),x.data(),y.data());
  g->SetLineColor(color);
  g->SetLineWidth(2);
  g->SetMarkerColor(color);
  return g;
}
void DrawLabel() {
  TLatex t;
  t.SetNDC();
  t.SetTextFont(72);
  t.SetTextSize(0.043);
  t.DrawLatex(0.17,0.86,"NA60+");
  t.SetTextFont(42);
  t.SetTextSize(0.037);
  t.DrawLatex(0.34,0.86,"Simulation");
  t.SetTextSize(0.030);
  t.DrawLatex(0.17,0.80,"P3 noFW | one positive ion | provisional mobility");
}
void DrawSide(const std::string& dir, const std::vector<Row>& rows,
              const int maxStep, bool plus) {
  const std::string s=plus?"plus":"minus";
  auto* canvas=new TCanvas(("c_"+s).c_str(),
                           ("NA60+ single-ion "+s).c_str(),1020,960);
  canvas->Divide(1,2,0.,0.);
  TMultiGraph* currents=new TMultiGraph();
  TMultiGraph* charges=new TMultiGraph();
  auto* leg=new TLegend(0.77,0.67,0.95,0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  const int colors[5]={kBlack,kBlue+1,kRed+1,kGreen+2,kMagenta+2};
  std::vector<TGraph*> hold;
  int j=0;
  double reference=std::numeric_limits<double>::quiet_NaN();
  for (const auto& row : rows) {
    if(row.step!=maxStep)continue;
    const std::string file=dir+"/wave_"+std::to_string(row.step)+"_"+
                           std::to_string(row.dtIndex)+".csv";
    const auto p=ReadWave(file);
    const int color=colors[j%5];
    auto* current=MakeGraph(p.t,plus?p.iPlus:p.iMinus,color);
    auto* charge=MakeGraph(p.t,plus?p.qPlus:p.qMinus,color);
    currents->Add(current,"L"); charges->Add(charge,"L");
    leg->AddEntry(current,Form("#Deltat = %.3g ns",row.dt),"l");
    hold.push_back(current);hold.push_back(charge);
    reference=plus?row.qPlus:row.qMinus;
    ++j;
  }
  if(j==0)throw std::runtime_error("Missing last drift-step waveforms");
  canvas->cd(1);
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.04);
  gPad->SetTopMargin(0.09);
  gPad->SetBottomMargin(0.13);
  currents->Draw("AL");
  currents->GetXaxis()->SetTitle("Time [#mus]");
  currents->GetYaxis()->SetTitle(("Cathode "+s+" ion current [fC/ns]").c_str());
  currents->GetYaxis()->SetTitleOffset(1.5);
  DrawLabel();
  leg->Draw();
  canvas->cd(2);
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.04);
  gPad->SetTopMargin(0.09);
  gPad->SetBottomMargin(0.13);
  charges->Draw("AL");
  charges->GetXaxis()->SetTitle("Time [#mus]");
  charges->GetYaxis()->SetTitle(("Cathode "+s+" accumulated charge [fC]").c_str());
  charges->GetYaxis()->SetTitleOffset(1.5);
  // Endpoint Ramo reference: signed Q = e [phi_end-phi_start].
  const double xMin=charges->GetXaxis()->GetXmin();
  const double xMax=charges->GetXaxis()->GetXmax();
  auto* ramo=new TLine(xMin,reference,xMax,reference);
  ramo->SetLineStyle(2);
  ramo->SetLineWidth(2);
  ramo->SetLineColor(kGray+2);
  ramo->Draw();
  TLatex t;
  t.SetNDC();
  t.SetTextFont(42);
  t.SetTextSize(0.033);
  t.DrawLatex(0.17,0.87,"Dashed: endpoint Shockley-Ramo integral");
  canvas->SaveAs((dir+"/single_ion_"+s+".png").c_str());
  canvas->SaveAs((dir+"/single_ion_"+s+".pdf").c_str());
}
} // namespace

void plot_single_ion(const char* directory) {
  gStyle->SetOptStat(0);
  gStyle->SetTitleFont(42,"XYZ");
  gStyle->SetLabelFont(42,"XYZ");
  gStyle->SetTitleSize(0.048,"XYZ");
  gStyle->SetLabelSize(0.043,"XYZ");
  try {
    const std::string dir=directory;
    const auto rows=ReadSummary(dir+"/single_ion_summary.csv");
    int maxStep=0;
    for(const auto& r:rows)maxStep=std::max(maxStep,r.step);
    DrawSide(dir,rows,maxStep,false);
    DrawSide(dir,rows,maxStep,true);
    std::cout<<"ROOT plots written to "<<dir<<"\n";
  }catch(const std::exception& e){
    std::cerr<<"ROOT PLOT FAILURE: "<<e.what()<<"\n";
    throw;
  }
}
