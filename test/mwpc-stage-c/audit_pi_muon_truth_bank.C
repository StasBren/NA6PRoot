// Read-only audit of existing FULL-geometry 50k-parent Geant4 ROOT chunks.
// Native ROOT remains the LOSSLESS transport/event bank for later prototypes.
// No change to frozen Geant4 physics, generator, or original exporter macro.
#include "NA6PMuonSpecHit.h"
#include <TFile.h>
#include <TTree.h>
#include <TParticle.h>
#include <TSystem.h>
#include <TBranch.h>
#include <TObjArray.h>
#include <cmath>
#include <memory>
#include <iostream>
#include <vector>
#include <string>

void audit_pi_muon_truth_bank(const char* directory, int expected_events = 50000) {
  const std::string dir(directory ? directory : "");
  auto fail=[&](const std::string& reason) {
    std::cerr << "MUON_TRUTH_BANK_AUDIT_FAIL: " << reason << " | " << dir << '\n';
    gSystem->Exit(4);
  };
  if (dir.empty() || expected_events!=50000) {
    fail("50k-parent chunk and its directory required");
    return;
  }
  std::unique_ptr<TFile> kine(TFile::Open((dir+"/MCKine.root").c_str()));
  std::unique_ptr<TFile> hits(TFile::Open((dir+"/HitsMuonSpec.root").c_str()));
  std::unique_ptr<TFile> planes(TFile::Open((dir+"/MuonReferencePlanes.root").c_str()));
  std::unique_ptr<TFile> geo(TFile::Open((dir+"/geometry.root").c_str()));
  if (!kine || kine->IsZombie() || !hits || hits->IsZombie() ||
      !planes || planes->IsZombie() || !geo || geo->IsZombie()) {
    fail("Missing/corrupt MCKine/HitsMuonSpec/MuonReferencePlanes/geometry");
    return;
  }
  auto* kt=dynamic_cast<TTree*>(kine->Get("mckine"));
  auto* ht=dynamic_cast<TTree*>(hits->Get("hitsMuonSpec"));
  auto* gt=dynamic_cast<TTree*>(hits->Get("gasStepsMuonSpec"));
  auto* pt=dynamic_cast<TTree*>(planes->Get("Planes"));
  if (!kt || !ht || !gt || !pt || kt->GetEntries()!=expected_events ||
      ht->GetEntries()!=expected_events) {
    fail("Missing trees or mismatched 50k event count");return;
  }
  std::vector<TParticle>* tracks=nullptr;
  std::vector<NA6PMuonSpecHit>* hitvec=nullptr;
  if (kt->SetBranchAddress("tracks",&tracks)<0 ||
      ht->SetBranchAddress("MuonSpec",&hitvec)<0) {
    fail("Cannot read event-index aligned TParticle / NA6PMuonSpecHit");return;
  }
  for (const auto* key : {"event","trackID","detID","u_mm","v_mm",
                          "w_mm","edep_GeV"}) {
    if (!gt->GetBranch(key)) {fail(std::string("Gas-step branch absent: ")+key);return;}
  }
  for (const auto* key : {"event","trackID","pdg","planeID","forward","z_cm",
                          "px_GeV","py_GeV","pz_GeV","z_step_start_cm",
                          "z_step_end_cm"}) {
    if (!pt->GetBranch(key)) {fail(std::string("Plane branch absent: ")+key);return;}
  }
  std::cout<<"MUON_TRUTH_BANK_AUDIT: "<<dir<<"\n";
  std::cout<<"  exact generated event entries="<<ht->GetEntries()
           <<", Geant4 positive-Edep gas steps="<<gt->GetEntries()
           <<", interpolated reference plane states="<<pt->GetEntries()<<"\n";
  std::cout<<"  full per-hit position/exit momentum/time/Edep/status: "
              "NA6PMuonSpecHit object in original ROOT; per-event mother "
              "and kinematics: matching TParticle vector\n";
  std::cout<<"  caution: gasStepsMuonSpec & Planes original trackIDs "
              "MUST NOT be blindly joined to remapped MCKine trackIDs\n";
  auto branchlist=[&](const char* title,TTree* tree) {
    std::cout<<"  "<<title<<" branches:";
    for (int j=0;j<tree->GetListOfBranches()->GetEntries();++j) {
      auto* b=dynamic_cast<TBranch*>(tree->GetListOfBranches()->At(j));
      if (b)std::cout<<" "<<b->GetName();
    }
    std::cout<<"\n";
  };
  branchlist("gasSteps",gt);
  branchlist("referencePlanes",pt);
  long nmuon=0,nhits=0,ntracks=0,nmother=0;
  // Sample events throughout original 50k sequence, not merely at first event.
  // Full 50k integrity is handled by existing validate_pi_geant4_chunk.C;
  // this macro independently checks object-level ancestry and hit truth.
  constexpr int nprobe=200;
  for (int k=0;k<nprobe;++k) {
    Long64_t ev=(static_cast<Long64_t>(k)*expected_events)/nprobe;
    kt->GetEntry(ev);ht->GetEntry(ev);
    if (!tracks || !hitvec){fail("Null tracks/hits on sampled event");return;}
    ntracks+=tracks->size();
    for (const auto& hit:*hitvec) {
      const int id=hit.getTrackID();
      if (id<0 || static_cast<size_t>(id)>=tracks->size()) {
        fail("Remapped MuonSpec trackID out of MCKine bounds");return;
      }
      const auto& particle=tracks->at(id);
      const int mother=particle.GetFirstMother();
      if (mother>=0 && mother<static_cast<int>(tracks->size()))++nmother;
      // NA6PMuonSpecHit getters are float; TParticle getters are double.
      // Force an explicit double-valued list instead of letting C++ deduce
      // conflicting initializer_list element types (ACLiC compilation error).
      const std::vector<double> measured_values{
          hit.getXIn(),hit.getYIn(),hit.getZIn(),
          hit.getXOut(),hit.getYOut(),hit.getZOut(),
          hit.getPXIn(),hit.getPYIn(),hit.getPZIn(),
          hit.getPXOut(),hit.getPYOut(),hit.getPZOut(),
          hit.getTime(),hit.getHitValue(),
          particle.Px(),particle.Py(),particle.Pz(),
          particle.Vx(),particle.Vy(),particle.Vz()};
      for (const double x: measured_values) {
        if (!std::isfinite(x)) {
          fail("Nonfinite native particle vertex, momentum, hit or time");return;
        }
      }
      if (hit.getDetectorID()<0){fail("Negative detector ID");return;}
      ++nhits;
      if (std::abs(particle.GetPdgCode())==13)++nmuon;
    }
  }
  if (nmuon==0 || nhits==0 || ntracks==0) {
    fail("No usable sampled muon gas crossings and truth tracks");return;
  }
  std::cout<<"  checked "<<nprobe<<" distributed EVENTS: "<<ntracks
           <<" stored truth tracks, "<<nhits<<" chamber crossings ("
           <<nmuon<<" muon), "<<nmother
           <<" joined kept-mother links\n";
  std::cout<<"MUON_TRUTH_BANK_AUDIT_PASS\n";
}
