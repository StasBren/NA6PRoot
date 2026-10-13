#include <TFile.h>
#include <TTree.h>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

void export_geant_gas_steps(
    const char* inputFile = "HitsMuonSpec.root",
    const char* outputCsv = "geant_mwpc_gas_steps.csv")
{
  TFile in(inputFile, "READ");
  if (in.IsZombie()) {
    std::cerr << "Could not open " << inputFile << "\n";
    return;
  }

  auto* tree = dynamic_cast<TTree*>(in.Get("gasStepsMuonSpec"));
  if (!tree) {
    std::cerr
        << "Tree gasStepsMuonSpec not found in " << inputFile << ".\n"
        << "Run NA6P with NA6P_MWPC_STORE_GAS_STEPS=1.\n";
    return;
  }

  int event = 0;
  int trackID = 0;
  int detID = 0;
  double x = 0., y = 0., z = 0.;
  double u = 0., v = 0., w = 0.;
  double px = 0., py = 0., pz = 0.;
  double edep = 0.;

  tree->SetBranchAddress("event", &event);
  tree->SetBranchAddress("trackID", &trackID);
  tree->SetBranchAddress("detID", &detID);
  tree->SetBranchAddress("x_cm", &x);
  tree->SetBranchAddress("y_cm", &y);
  tree->SetBranchAddress("z_cm", &z);
  tree->SetBranchAddress("u_mm", &u);
  tree->SetBranchAddress("v_mm", &v);
  tree->SetBranchAddress("w_mm", &w);
  tree->SetBranchAddress("px_GeV", &px);
  tree->SetBranchAddress("py_GeV", &py);
  tree->SetBranchAddress("pz_GeV", &pz);
  tree->SetBranchAddress("edep_GeV", &edep);

  std::ofstream out(outputCsv);
  if (!out) {
    std::cerr << "Could not create " << outputCsv << "\n";
    return;
  }

  out << std::setprecision(12);
  out
      << "event,track_id,det_id,"
      << "x_cm,y_cm,z_cm,u_mm,v_mm,w_mm,"
      << "px_GeV,py_GeV,pz_GeV,edep_GeV,edep_eV\n";

  const Long64_t n = tree->GetEntries();
  double sumEdepGeV = 0.;
  for (Long64_t i = 0; i < n; ++i) {
    tree->GetEntry(i);
    sumEdepGeV += edep;
    out
        << event << "," << trackID << "," << detID << ","
        << x << "," << y << "," << z << ","
        << u << "," << v << "," << w << ","
        << px << "," << py << "," << pz << ","
        << edep << "," << 1.e9 * edep << "\n";
  }

  std::cout
      << "Exported " << n << " non-zero MWPC gas energy-deposit steps\n"
      << "Total Edep = " << 1.e9 * sumEdepGeV << " eV\n"
      << "Wrote      = " << outputCsv << "\n";
}
