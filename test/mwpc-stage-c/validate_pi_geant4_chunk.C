// Validate one native NA6PRoot Geant4 production chunk WITHOUT reading
// object vectors: no detector dictionaries are needed merely to inspect
// TTree metadata, so this check is robust to ROOT class loading order.
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <iostream>
#include <memory>
#include <string>

void validate_pi_geant4_chunk(const char* directory, int expectedEvents)
{
  const std::string dir(directory ? directory : "");
  auto fail = [&](const std::string& msg) {
    std::cerr << "PI_GEANT_CHUNK_FAIL: " << msg << " | " << dir << '\n';
    gSystem->Exit(2);
  };
  if (dir.empty() || expectedEvents <= 0) {
    fail("Invalid output directory or expected event count");
    return;
  }

  std::unique_ptr<TFile> kine(TFile::Open((dir + "/MCKine.root").c_str()));
  std::unique_ptr<TFile> hits(TFile::Open((dir + "/HitsMuonSpec.root").c_str()));
  std::unique_ptr<TFile> geometry(TFile::Open((dir + "/geometry.root").c_str()));
  std::unique_ptr<TFile> reference(
      TFile::Open((dir + "/MuonReferencePlanes.root").c_str()));
  if (!kine || kine->IsZombie() || !hits || hits->IsZombie() ||
      !geometry || geometry->IsZombie() ||
      !reference || reference->IsZombie()) {
    fail("Missing/unreadable mandatory MCKine, HitsMuonSpec, geometry or MuonReferencePlanes ROOT file");
    return;
  }
  auto* mt = dynamic_cast<TTree*>(kine->Get("mckine"));
  auto* ht = dynamic_cast<TTree*>(hits->Get("hitsMuonSpec"));
  auto* st = dynamic_cast<TTree*>(hits->Get("gasStepsMuonSpec"));
  auto* rt = dynamic_cast<TTree*>(reference->Get("Planes"));
  if (!mt || !ht || !st || !rt) {
    fail("Missing mckine, hitsMuonSpec, gasStepsMuonSpec or Planes tree. "
         "The last two require rebuilt Stage-C simulator and BOTH "
         "NA6P_MWPC_STORE_GAS_STEPS / NA6P_MWPC_STORE_REFERENCE_PLANES.");
    return;
  }
  if (!mt->GetBranch("tracks") || !ht->GetBranch("MuonSpec")) {
    fail("Missing tracks or MuonSpec object branch");
    return;
  }
  if (mt->GetEntriesFast() != expectedEvents ||
      ht->GetEntriesFast() != expectedEvents) {
    fail("Event counts differ from requested count");
    return;
  }
  for (const char* col : {"event", "trackID", "detID", "u_mm",
                          "v_mm", "w_mm", "edep_GeV"}) {
    if (!st->GetBranch(col)) {
      fail(std::string("Missing gas-step branch ") + col);
      return;
    }
  }
  for (const char* col : {"event", "trackID", "pdg", "planeID",
                          "forward", "z_cm", "px_GeV", "py_GeV", "pz_GeV",
                          "z_step_start_cm", "z_step_end_cm"}) {
    if (!rt->GetBranch(col)) {
      fail(std::string("Missing portable muon reference branch ") + col);
      return;
    }
  }
  std::cout << "PI_GEANT_CHUNK_PASS events=" << expectedEvents
            << " virtualPlaneStates=" << rt->GetEntriesFast()
            << " gasSteps=" << st->GetEntriesFast()
            << " trees OK | " << dir << std::endl;
  gSystem->Exit(0);
}
