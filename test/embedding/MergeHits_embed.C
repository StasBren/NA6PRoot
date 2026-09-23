#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TParticle.h>
#include "NA6PMCEventHeader.h"
#include "NA6PMCGenHeader.h"
#include "NA6PMuonSpecHit.h"
#include "NA6PVerTelHit.h"


void MergeHits_embed(int bckEvent,
                     const char *fKS = "MCKine.root",
                     const char *fKB = "MCKine_bck.root",
                     const char *fHSM = "HitsMuonSpec.root",
                     const char *fHBM = "HitsMuonSpec_bck.root",
                     const char *fHSV = "HitsVerTel.root",
                     const char *fHBV = "HitsVerTel_bck.root"
		     )
{
  // Background event under study
  const char* evStr = gSystem->Getenv("BKG_EVENT");
  bckEvent = evStr ? atoi(evStr) : 0;

  // ===============================
  // Files
  // ===============================
  TFile* fKSig = TFile::Open(fKS);
  TFile* fKBkg = TFile::Open(fKB);
  TFile* fKOut = TFile::Open(Form("MCKine_mix_%d.root", bckEvent),"RECREATE");
  
  TFile* fHSigM = TFile::Open(fHSM);
  TFile* fHBkgM = TFile::Open(fHBM);
  TFile* fHOutM = TFile::Open(Form("HitsMuonSpec_mix_%d.root", bckEvent), "RECREATE");

  TFile* fHSigV = TFile::Open(fHSV);
  TFile* fHBkgV = TFile::Open(fHBV);
  TFile* fHOutV = TFile::Open(Form("HitsVerTel_mix_%d.root", bckEvent), "RECREATE");

  // ===============================
  // Trees
  // ===============================
  TTree* tSigKin  = (TTree*)fKSig->Get("mckine");
  TTree* tSigHitM = (TTree*)fHSigM->Get("hitsMuonSpec");
  TTree* tSigHitV = (TTree*)fHSigV->Get("hitsVerTel");

  TTree* tBkgKin  = (TTree*)fKBkg->Get("mckine");
  TTree* tBkgHitM = (TTree*)fHBkgM->Get("hitsMuonSpec");
  TTree* tBkgHitV = (TTree*)fHBkgV->Get("hitsVerTel");

  // Clone structure ONLY
  fKOut->cd();
  TTree* tOutKin = new TTree("mckine", "mixed kinematics");
  std::vector<TParticle>* mixmcArr = new std::vector<TParticle>();
  tOutKin->Branch("tracks", &mixmcArr);
  NA6PMCEventHeader* mixmcHead = nullptr;
  tOutKin->Branch("header", &mixmcHead);

  fHOutM->cd();
  TTree* tOutHitM = tSigHitM->CloneTree(0); 
  std::vector<NA6PMuonSpecHit>* mixHitsM = new std::vector<NA6PMuonSpecHit>();
  tOutHitM->SetBranchAddress("MuonSpec", &mixHitsM);

  fHOutV->cd();
  TTree* tOutHitV = tSigHitV->CloneTree(0);
  std::vector<NA6PVerTelHit>* mixHitsV = new std::vector<NA6PVerTelHit>();
  tOutHitV->SetBranchAddress("VerTel", &mixHitsV);

  // ===============================
  // Input branches
  // ===============================

  std::vector<TParticle>* sigmcArr = nullptr;
  std::vector<TParticle>* bkgmcArr = nullptr;

  NA6PMCEventHeader* sigmcHead = nullptr;
  NA6PMCEventHeader* bkgmcHead = nullptr;

  tSigKin->SetBranchAddress("header", &sigmcHead);
  tSigKin->SetBranchAddress("tracks", &sigmcArr);
  tBkgKin->SetBranchAddress("header", &bkgmcHead);
  tBkgKin->SetBranchAddress("tracks", &bkgmcArr);

  std::vector<NA6PMuonSpecHit>* sigHitsM = nullptr;
  std::vector<NA6PMuonSpecHit>* bkgHitsM = nullptr;
  std::vector<NA6PVerTelHit>* sigHitsV = nullptr;
  std::vector<NA6PVerTelHit>* bkgHitsV = nullptr;

  tSigHitM->SetBranchAddress("MuonSpec", &sigHitsM);
  tBkgHitM->SetBranchAddress("MuonSpec", &bkgHitsM);
  tSigHitV->SetBranchAddress("VerTel", &sigHitsV);
  tBkgHitV->SetBranchAddress("VerTel", &bkgHitsV);

  // ===============================
  // Event loop
  // ===============================
  Long64_t nEvents = tSigKin->GetEntries();
  
  tBkgKin->GetEntry(bckEvent);
  tBkgHitM->GetEntry(bckEvent);
  tBkgHitV->GetEntry(bckEvent);
    
  for (Long64_t iev = 0; iev < nEvents; ++iev) {

    tSigKin->GetEntry(iev);
    tSigHitM->GetEntry(iev);
    tSigHitV->GetEntry(iev);

// Kinematics HEADER
    
    mixmcHead->clear();    
      
// Get event header info for signal and background
    int eventIDs = sigmcHead->getEventID();
    int eventIDb = bkgmcHead->getEventID();
    int RunNumbers = sigmcHead->getRunNumber();
    int RunNumberb = bkgmcHead->getRunNumber();
    float VXs = sigmcHead->getVX();
    float VYs = sigmcHead->getVY();
    float VZs = sigmcHead->getVZ();
    float VXb = bkgmcHead->getVX();
    float VYb = bkgmcHead->getVY();
    float VZb = bkgmcHead->getVZ();
    int NTrackss = sigmcHead->getNTracks();
    int NTracksb = bkgmcHead->getNTracks();
    int NPrimariess = sigmcHead->getNPrimaries();
    int NPrimariesb = bkgmcHead->getNPrimaries();
    cout << endl;
    cout << "Mixing signal eventID " <<  eventIDs << " + background eventID " <<  eventIDb << endl;
    cout << "Sig vertex x " << VXs << " y " << VYs << " z " << VZs << " Ntracks " << NTrackss << " Nprimaries " << NPrimariess << endl; 
    cout << "Bck vertex x " << VXb << " y " << VYb << " z " << VZb << " Ntracks " << NTracksb << " Nprimaries " << NPrimariesb << endl; 

// building the header of the mixed event 
    mixmcHead->setEventID(nEvents*bckEvent+iev); // build progressive event number
    mixmcHead->setRunNumber(RunNumberb); // use RunNumber from background file (looks not in use until now) 
    mixmcHead->setVX(VXs); // vertex position of signal and background events are the same by construction
    mixmcHead->setVY(VYs);
    mixmcHead->setVZ(VZs);
    mixmcHead->setNTracks(NTrackss+NTracksb); // total number of tracks for mixed event
    mixmcHead->setNPrimaries(NPrimariess+NPrimariesb); // total number of primaries for mixed event

// add generator headers for signal and background
    for(int ighs = 0; ighs < (int)sigmcHead->getNGenHeaders(); ighs++){
      const NA6PMCGenHeader& genHeadsig = sigmcHead->getGenHeader(ighs); // 
      cout << " Sig gen header " << ighs << " primary offset signal " <<  genHeadsig.getPrimariesOffset() << " secondary offset signal " << genHeadsig.getSecondariesOffset()  << " info " << genHeadsig.getInfo() << endl;
      if (genHeadsig.getNPrimaries() !=0 && genHeadsig.getNSecondaries() !=0) mixmcHead->addGenHeader(genHeadsig);
   }
    for(int ighb = 0; ighb < (int)bkgmcHead->getNGenHeaders(); ighb++){
      const NA6PMCGenHeader& genHeadbck = bkgmcHead->getGenHeader(ighb); // 
      cout << " Bck gen header " << ighb << " primary offset signal " <<  genHeadbck.getPrimariesOffset() << " secondary offset signal " << genHeadbck.getSecondariesOffset()  << " info " << genHeadbck.getInfo() << endl;
      if (genHeadbck.getNPrimaries() !=0 && genHeadbck.getNSecondaries() !=0) mixmcHead->addGenHeader(genHeadbck);
    }
    
// Check
    cout << " Mixed event n. " << mixmcHead->getEventID() << " Total tracks " << mixmcHead->getNTracks() << " Total primaries " << mixmcHead->getNPrimaries() << endl;
    for(int ighm = 0; ighm < (int)mixmcHead->getNGenHeaders(); ighm++){
    const NA6PMCGenHeader& genHeadmix = mixmcHead->getGenHeader(ighm); // 
    cout << " gen header " <<  genHeadmix.getInfo() << endl;
    cout << " n. primaries " << genHeadmix.getNPrimaries() << " n. secondaries " << genHeadmix.getNSecondaries() << endl; 
    cout << " primaries offset " << genHeadmix.getPrimariesOffset() << " secondaries offset " << genHeadmix.getSecondariesOffset() << endl; 
   }


    // -----------------------------
    // 1) Kinematics Signal 
    // -----------------------------

    mixmcArr->clear();

    int nsigkin=0;
    
    for (const auto& trk : *sigmcArr) {
      mixmcArr->push_back(trk);
//    cout << " Signal event " << iev << " Sig track " << nsigkin << "  is a " << trk.GetName() << " first mother is track n. " << trk.GetFirstMother() << " first daughter is track n. " << trk.GetFirstDaughter() << endl;
      nsigkin++;
    }

    int offset = mixmcArr->size(); // to be added to background tracks

    // -----------------------------
    // 2) Kinematics Background 
    // -----------------------------
    
    int nbckkin=0;
    
    for (auto trk : *bkgmcArr) {

//    cout << "  Background event " << bckEvent << " Bck track " << nbckkin << "  is a " << trk.GetName() << " first mother is track n. " << trk.GetFirstMother() << " first daughter is track n. " << trk.GetFirstDaughter() << endl;

      nbckkin++;
      int fm = trk.GetFirstMother();
      if (fm >= 0) trk.SetFirstMother(fm + offset);
      int fd = trk.GetFirstDaughter();
      int ld = trk.GetLastDaughter();
      if (fd >= 0) trk.SetFirstDaughter(fd + offset);
      if (ld >= 0) trk.SetLastDaughter(ld + offset);

      mixmcArr->push_back(trk);
    } 

// Check
    int nmixkin=0;
    for (auto trk : *mixmcArr) {
//    cout << " Signal + Background event " << iev << " Mixed track " << nmixkin << "  is a " << trk.GetName() << " first mother is track n. " << trk.GetFirstMother() << " first daughter is track n. " << trk.GetFirstDaughter()<< endl;
      nmixkin++;
    }

    // -----------------------------
    // 3) Signal hits
    // -----------------------------
    
    mixHitsM->clear();
    mixHitsV->clear();

    int nsighitM=0;
    int nsighitV=0;

    // Muon spectrometer
    for (const auto& h : *sigHitsM){
      mixHitsM->push_back(h);
//    cout << " Signal event " << iev << " Signal Hit Muon " << nsighitM << "  x = " << h.getX() << "  y = " << h.getY() <<"  z = " << h.getZ() <<" ref. mixed track n. " << h.getTrackID() << endl;
      nsighitM++;      
    }

    // Vertex spectrometer    
    for (const auto& h : *sigHitsV){
      mixHitsV->push_back(h);
//    cout << " Signal event " << iev << " Signal Hit VT " << nsighitV << "  x = " << h.getX() << "  y = " << h.getY() <<"  z = " << h.getZ() <<" ref. mixed track n. " << h.getTrackID() << endl;
      nsighitV++;      
    }

    // -----------------------------
    // 4) Background hits
    // -----------------------------
    
    int nbckhitM=0;
    int nbckhitV=0;
    
    for (const auto& hit : *bkgHitsM) {
      mixHitsM->emplace_back(hit).setTrackID(hit.getTrackID() + offset); // adjust reference to track number
   //      cout << " Background event " << bckEvent << " Bck Hit Muon " << nbckhitM << "  x = " << hit.getX() << "  y = " << h.getY() <<"  z = " << h.getZ() << " ref. mixed track n. " << hit.getTrackID() << endl;
      nbckhitM++;
    }

    for (const auto& hit : *bkgHitsV) {
      mixHitsV->emplace_back(hit).setTrackID(hit.getTrackID() + offset); // adjust reference to track number
//      cout << " Background event " << bckEvent << " Bck Hit VT " << nbckhitV << "  x = " << hit.getX() << "  y = " << h.getY() <<"  z = " << h.getZ() << " ref. mixed track n. " << hit.getTrackID() << endl;
      nbckhitV++;
    }

    // -----------------------------
    // Fill mixed event trees
    // -----------------------------
    fKOut->cd();
    tOutKin->Fill();
    fHOutM->cd();
    tOutHitM->Fill();
    fHOutV->cd();
    tOutHitV->Fill();
    
  } // End event loop

  // ===============================
  // Write output
  // ===============================
  fKOut->cd();
  tOutKin->Write();
  fKOut->Close();

  fHOutM->cd();
  tOutHitM->Write();
  fHOutM->Close();

  fHOutV->cd();
  tOutHitV->Write();
  fHOutV->Close();

  cout << "Mixed " << nEvents << " events successfully " << endl;
}
