#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TTree.h>
#include <TFile.h>
#include <TParticle.h>
#include <TParticlePDG.h>
#include <TVector3.h>
#include <TRandom3.h>
#include "NA6PVerTelHit.h"
#include "NA6PMuonSpecHit.h"
#include "ConfigurableParam.h"
#include <fairlogger/Logger.h>
#include "NA6PVerTelReconstruction.h"
#include "NA6PMuonSpecHit.h"
#include "NA6PMuonSpecReconstruction.h"
#endif

void convertHitsToRecPoints(bool simulateDeadVTRegions = true)
{
  // Process VerTel hits
  TFile* fhVT = TFile::Open("HitsVerTel.root");
  if (!fhVT || fhVT->IsZombie()) {
    LOGP(error, "Cannot open HitsVerTel.root");
    if (fhVT)
      delete fhVT;
  } else {
    TTree* thVT = (TTree*)fhVT->Get("hitsVerTel");
    if (!thVT) {
      LOGP(error, "Cannot find tree 'hitsVerTel' in HitsVerTel.root");
    } else {
      std::vector<NA6PVerTelHit> vtHits, *vtHitsPtr = &vtHits;
      thVT->SetBranchAddress("VerTel", &vtHitsPtr);
      int nEvVT = thVT->GetEntriesFast();

      NA6PVerTelReconstruction* recVT = new NA6PVerTelReconstruction();
      recVT->createClustersOutput();
      for (int jEv = 0; jEv < nEvVT; jEv++) {
        thVT->GetEvent(jEv);
        int nHits = vtHits.size();
        printf("VerTel Event %d nHits=%d\n", jEv, nHits);
        recVT->clearClusters();
        recVT->hitsToRecPoints(vtHits, jEv, simulateDeadVTRegions);
        recVT->writeClusters();
      }
      recVT->closeClustersOutput();
      delete recVT;
    }
    fhVT->Close();
    delete fhVT;
  }

  // Process MuonSpec hits
  TFile* fhMS = TFile::Open("HitsMuonSpec.root");
  if (!fhMS || fhMS->IsZombie()) {
    LOGP(error, "Cannot open HitsMuonSpec.root");
    if (fhMS)
      delete fhMS;
  } else {
    TTree* thMS = (TTree*)fhMS->Get("hitsMuonSpec");
    if (!thMS) {
      LOGP(error, "Cannot find tree 'hitsMuonSpec' in HitsMuonSpec.root");
    } else {
      std::vector<NA6PMuonSpecHit> msHits, *msHitsPtr = &msHits;
      thMS->SetBranchAddress("MuonSpec", &msHitsPtr);
      int nEvMS = thMS->GetEntriesFast();

      NA6PMuonSpecReconstruction* recMS = new NA6PMuonSpecReconstruction();
      recMS->createClustersOutput();
      for (int jEv = 0; jEv < nEvMS; jEv++) {
        thMS->GetEvent(jEv);
        int nHits = msHits.size();
        printf("MuonSpec Event %d nHits=%d\n", jEv, nHits);
        recMS->clearClusters();
        recMS->hitsToRecPoints(msHits);
        recMS->writeClusters();
      }
      recMS->closeClustersOutput();
      delete recMS;
    }
    fhMS->Close();
    delete fhMS;
  }
}
