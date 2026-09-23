#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TString.h>
#include <TMath.h>
#include <TTree.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TH1F.h>
#include <TF1.h>
#include <TH2F.h>
#include <TFile.h>
#include <TGeoManager.h>
#include <TGeoNavigator.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TPaveStats.h>
#include <TParticle.h>
#include <TParticlePDG.h>
#include <TVector3.h>
#include <TRandom3.h>
#include "NA6PVerTelHit.h"
#include "NA6PMuonSpecHit.h"
#include "NA6PTrack.h"
#include "NA6PBaseCluster.h"
#include "NA6PFastTrackFitter.h"
#include "MagneticField.h"
#include "NA6PTreeStreamRedirector.h"
#include "ConfigurableParam.h"
#include "NA6PLayoutParam.h"
#endif

// simple macro with an example of how to use NA6PFastFitter to fit tracks starting from their hits in the VT
// uses the MC truth information to group the hits produced by the same particle

void SuperposHistos(TH1* h1, TH1* h2, TH1* h3)
{
  h1->SetLineColor(1);
  h1->SetLineWidth(2);
  if (h2->GetMaximum() > h1->GetMaximum())
    h1->SetMaximum(1.02 * h2->GetMaximum());
  if (h3->GetMaximum() > h1->GetMaximum())
    h1->SetMaximum(1.02 * h3->GetMaximum());
  h1->SetMinimum(0);
  h1->Draw();
  h2->SetLineColor(2);
  h2->SetLineWidth(2);
  h2->Draw("sames");
  h3->SetLineColor(kGreen + 1);
  h3->SetLineWidth(2);
  h3->Draw("sames");
  gPad->Update();
  TPaveStats* st1 = (TPaveStats*)h1->GetListOfFunctions()->FindObject("stats");
  if (st1) {
    st1->SetY1NDC(0.72);
    st1->SetY2NDC(0.92);
    st1->SetTextColor(1);
  }
  TPaveStats* st2 = (TPaveStats*)h2->GetListOfFunctions()->FindObject("stats");
  if (st2) {
    st2->SetY1NDC(0.51);
    st2->SetY2NDC(0.71);
    st2->SetTextColor(2);
  }
  TPaveStats* st3 = (TPaveStats*)h3->GetListOfFunctions()->FindObject("stats");
  if (st3) {
    st3->SetY1NDC(0.30);
    st3->SetY2NDC(0.50);
    st3->SetTextColor(kGreen + 1);
  }
  gPad->Modified();
}

void FillMeanAndRms(TH2F* hImpParVsP, TH1F* hImpParMean, TH1F* hImpParRms, TH1F* hImpParSig)
{
  if (!hImpParVsP || !hImpParMean || !hImpParRms || !hImpParSig) {
    printf("One of pointers is not set: hImpParVsP=%p hImpParMean=%p hImpParRms=%p hImpParSig=%p\n", hImpParVsP, hImpParMean, hImpParRms, hImpParSig);
    return;
  }
  int nMomBins = hImpParVsP->GetXaxis()->GetNbins();
  hImpParMean->SetStats(0);
  hImpParRms->SetStats(0);
  hImpParSig->SetStats(0);
  for (int jp = 1; jp <= nMomBins; jp++) {
    TH1D* htmp1 = hImpParVsP->ProjectionY("htmp1", jp, jp);
    float mean = htmp1->GetMean();
    float emean = htmp1->GetMeanError();
    float rms = htmp1->GetRMS();
    float erms = htmp1->GetRMSError();
    hImpParMean->SetBinContent(jp, mean);
    hImpParMean->SetBinError(jp, emean);
    hImpParRms->SetBinContent(jp, rms);
    hImpParRms->SetBinError(jp, erms);
    float gw = rms;
    float egw = erms;
    if (htmp1->GetEntries() > 20.) {
      int fitStatus = htmp1->Fit("gaus", "Q", "", -3. * rms, 3. * rms);
      if (fitStatus == 0) {
        TF1* fg = (TF1*)htmp1->GetListOfFunctions()->FindObject("gaus");
        if (fg) {
          gw = fg->GetParameter(2);
          egw = fg->GetParError(2);
        }
      }
    }
    hImpParSig->SetBinContent(jp, gw);
    hImpParSig->SetBinError(jp, egw);
    delete htmp1;
  }
}

// Template function to run fast fit on hits - works with both NA6PVerTelHit and NA6PMuonSpecHit
template <typename HitType>
void runFastFitOnHitsTemplate(int firstEv = 0,
                              int lastEv = 99999999,
                              float cluresx = 5.e-4,
                              float cluresy = 5.e-4,
                              const char* dirSimu = ".",
                              const char* hitFileName = "HitsVerTel.root",
                              const char* hitTreeName = "hitsVerTel",
                              const char* hitBranchName = "VerTel",
                              int nLayers = 5,
                              int minHits = 5,
                              bool verbose = false,
                              bool useVT = true)
{

  // na6p::conf::ConfigurableParam::updateFromFile(Form("%s/na6pLayout.ini",dirSimu),"",true);
  // const auto& param = NA6PLayoutParam::Instance();

  TFile* fk = new TFile(Form("%s/MCKine.root", dirSimu));
  TTree* mcTree = (TTree*)fk->Get("mckine");
  int nEv = mcTree->GetEntries();
  std::vector<TParticle>* mcArr = nullptr;
  mcTree->SetBranchAddress("tracks", &mcArr);

  TFile* fh = new TFile(Form("%s/%s", dirSimu, hitFileName));
  TTree* th = (TTree*)fh->Get(hitTreeName);
  std::vector<HitType> hits, *hitsPtr = &hits;
  th->SetBranchAddress(hitBranchName, &hitsPtr);

  int nMomBins = 20;
  float maxP = 10.;
  TH1F* hEtaGen = new TH1F("hEtaGen", ";#eta;counts", 20, 0., 5.);
  TH1F* hEtaTrackable = new TH1F("hEtaTrackable", ";#eta;counts", 20, 0., 5.);
  TH1F* hEtaFitted = new TH1F("hEtaFitted", ";#eta;counts", 20, 0., 5.);
  TH1F* hEtaReco = new TH1F("hEtaReco", ";#eta;counts", 20, 0., 5.);
  TH1F* hMomGen = new TH1F("hMomGen", ";p (GeV/C);counts", nMomBins, 0., maxP);
  TH1F* hMomTrackable = new TH1F("hMomTrackable", ";p (GeV/C);counts", nMomBins, 0., maxP);
  TH1F* hMomFitted = new TH1F("hMomFitted", ";p (GeV/C);counts", nMomBins, 0., maxP);
  TH1F* hDeltaPx = new TH1F("hDeltaPx", ";p_{x}^{rec}-p_{x}^{gen} (GeV/c);counts", 100, -0.2, 0.2);
  TH1F* hDeltaPy = new TH1F("hDeltaPy", ";p_{y}^{rec}-p_{y}^{gen} (GeV/c);counts", 100, -0.2, 0.2);
  TH1F* hDeltaPz = new TH1F("hDeltaPz", ";p_{z}^{rec}-p_{z}^{gen} (GeV/c);counts", 100, -0.5, 0.5);
  TH1F* hDeltaP = new TH1F("hDeltaP", ";p^{rec}-p^{gen} (GeV/c);counts", 100, -0.5, 0.5);
  TH1F* hDeltaPhi = new TH1F("hDeltaPhi", ";#varphi^{rec}-#varphi^{gen};counts", 100, -M_PI / 4., M_PI / 4.);
  float impMax = 60000.; // in microns
  float deltaMax = 1;    // in GeV/c
  if (useVT) {
    impMax = 500.;  // in microns
    deltaMax = 0.5; // in GeV/c
  }
  TH1F* hDeltaEta = new TH1F("hDeltaEta", ";#eta^{rec}-#eta^{gen};counts", 100, -0.5, 0.5);
  TH1F* hImpParX = new TH1F("hImpParX", ";Track Imp. Par. X (#mum)};counts", 100, -impMax, impMax);
  TH1F* hImpParY = new TH1F("hImpParY", ";Track Imp. Par. Y (#mum)};counts", 100, -impMax, impMax);
  TH2F* hImpParXVsP = new TH2F("hImpParXVsP", ";p (GeV/c);Track Imp. Par. X (#mum)};counts", nMomBins, 0., maxP, 100, -impMax, impMax);
  TH2F* hImpParYVsP = new TH2F("hImpParYVsP", ";p (GeV/c);Track Imp. Par. Y (#mum)};counts", nMomBins, 0., maxP, 100, -impMax, impMax);
  TH2F* hDeltaPxVsP = new TH2F("hDeltaPxVsP", ";p (GeV/c);p_{x}^{rec}-p_{x}^{gen} (GeV/c);counts", nMomBins, 0., maxP, 500, -deltaMax, deltaMax);
  TH2F* hDeltaPyVsP = new TH2F("hDeltaPyVsP", ";p (GeV/c);p_{y}^{rec}-p_{y}^{gen} (GeV/c);counts", nMomBins, 0., maxP, 500, -deltaMax, deltaMax);
  TH2F* hDeltaPzVsP = new TH2F("hDeltaPzVsP", ";p (GeV/c);p_{z}^{rec}-p_{z}^{gen} (GeV/c);counts", nMomBins, 0., maxP, 500, -deltaMax, deltaMax);
  TH2F* hDeltaPVsP = new TH2F("hDeltaPVsP", ";p (GeV/c);p^{rec}-p^{gen} (GeV/c);counts", nMomBins, 0., maxP, 100, -deltaMax, deltaMax);
  TH2F* hRelDeltaPxVsP = new TH2F("hRelDeltaPxVsP", ";p (GeV/c);(p_{x}^{rec}-p_{x}^{gen})/p_{x}^{gen};counts", nMomBins, 0., maxP, 100, -deltaMax, deltaMax);
  TH2F* hRelDeltaPyVsP = new TH2F("hRelDeltaPyVsP", ";p (GeV/c);(p_{y}^{rec}-p_{y}^{gen})/p_{y}^{gen};counts", nMomBins, 0., maxP, 100, -deltaMax, deltaMax);
  TH2F* hRelDeltaPzVsP = new TH2F("hRelDeltaPzVsP", ";p (GeV/c);(p_{z}^{rec}-p_{z}^{gen})/p_{z}^{gen};counts", nMomBins, 0., maxP, 100, -deltaMax, deltaMax);
  TH2F* hRelDeltaPVsP = new TH2F("hRelDeltaPVsP", ";p (GeV/c);(p^{rec}-p^{gen})/p_{gen};counts", nMomBins, 0., maxP, 100, -deltaMax, deltaMax);
  TH1F* hIsGood = new TH1F("hIsGood", "", 2, -0.5, 1.5);
  hIsGood->GetXaxis()->SetBinLabel(1, "bad");
  hIsGood->GetXaxis()->SetBinLabel(2, "good");
  TH1F* hNclu = new TH1F("hNclu", ";n_{ITSclus};counts", 7, -0.5, 6.5);
  TH1F* hChi2NDF = new TH1F("hChi2NDF", ";#chi^{2}/nDOF;counts", 100, 0., 10.);

  if (lastEv > nEv || lastEv < 0)
    lastEv = nEv;
  if (firstEv < 0)
    firstEv = 0;

  NA6PTreeStreamRedirector outTr("fitCheck.root");
  if (!Propagator::loadField() || !Propagator::loadGeometry(Form("%s/geometry.root", dirSimu))) {
    return;
  }
  NA6PFastTrackFitter* fitter = new NA6PFastTrackFitter();
  fitter->setMaxChi2Cl(100.);
  fitter->setPropagateToPrimaryVertex(true);
  // fitter->setSeedFromTwoOutermostHits();
  // fitter->setCharge(2);
  // fitter->setParticleHypothesisPDG(2212);
  // fitter->disableMaterialCorrections();
  const auto& layout = NA6PLayoutParam::Instance();

  for (int jEv = firstEv; jEv < lastEv; jEv++) {
    mcTree->GetEvent(jEv);
    th->GetEvent(jEv);
    int nPart = mcArr->size();
    int nHits = hits.size();
    float xvert = 0;
    float yvert = 0;
    float zvert = 0;
    // get primary vertex position from the Kine Tree
    for (int jp = 0; jp < nPart; jp++) {
      auto curPart = mcArr->at(jp);
      if (curPart.IsPrimary()) {
        xvert = curPart.Vx();
        yvert = curPart.Vy();
        zvert = curPart.Vz();
        break;
      }
    }
    fitter->setPrimaryVertexZ(zvert);
    for (int jp = 0; jp < nPart; jp++) {
      const auto& curPart = mcArr->at(jp);
      float pxPart = curPart.Px();
      float pyPart = curPart.Py();
      float pzPart = curPart.Pz();
      float momPart = curPart.P();
      float phiPart = curPart.Phi();
      float thetaPart = std::acos(pzPart / momPart);
      float etaPart = -std::log(std::tan(thetaPart / 2.));
      std::array<std::unique_ptr<NA6PBaseCluster>, 6> clusters{};
      float xclu[6], yclu[6], zclu[6];
      if (curPart.IsPrimary()) {
        hMomGen->Fill(momPart);
        hEtaGen->Fill(etaPart);
        int maskHits = 0;
        fitter->cleanupAndStartFit();
        for (const auto& hit : hits) {
          int idPart = hit.getTrackID();
          if (idPart == jp) {
            int nDet = hit.getDetectorID();
            int nLay = -1; // initialize to invalid value
            if (useVT) {
              nLay = nDet / 4; // Get layer from detector ID for VT
            } else {
              // Lookup layer from layout parameters for MS
              float hitZ = (hit.getZIn() + hit.getZOut()) / 2.f; // use average Z
              const float dzWindow = 20.f;                       // half-width window
              for (int i = 0; i < layout.nMSPlanes; ++i) {
                float zPlane = layout.posMSPlaneZ[i];
                if (hitZ > (zPlane - dzWindow) && hitZ < (zPlane + dzWindow)) {
                  nLay = i;
                  break;
                }
              }
              if (nLay < 0)
                continue; // skip hits not in expected z ranges
            }
            if (verbose)
              printf("nLay %d nLayers %d\n", nLay, nLayers);
            if (nLay < 0 || nLay >= nLayers)
              continue;              // safety check: skip if layer index out of bounds
            maskHits |= (1 << nLay); // use bitmask instead of counter
            if (cluresx > 0 && cluresy > 0) {
              xclu[nLay] = gRandom->Gaus(hit.getX(), cluresx);
              yclu[nLay] = gRandom->Gaus(hit.getY(), cluresy);
              zclu[nLay] = hit.getZ();
            } else {
              xclu[nLay] = hit.getX();
              yclu[nLay] = hit.getY();
              zclu[nLay] = hit.getZ();
            }
            clusters[nLay] = std::make_unique<NA6PBaseCluster>(xclu[nLay], yclu[nLay], zclu[nLay], idPart, nLay);
            clusters[nLay]->setDetectorID(nDet);
            if (cluresx > 0 && cluresy > 0)
              clusters[nLay]->setErr(cluresx * cluresx, 0., cluresy * cluresy);
            else
              clusters[nLay]->setErr(1.e-6, 0., 1e-6);
          }
        }
        int expectedMask = (1 << nLayers) - 1; // e.g., for 6 layers: 63 = 0b111111

        if (maskHits != expectedMask) {
          continue;
        }

        for (int nLay = 0; nLay < (int)clusters.size(); ++nLay) {
          if (clusters[nLay]) {
            if (verbose)
              std::cout << "Adding cluster at layer " << nLay << " for track " << jp << "\n";
            fitter->addCluster(*clusters[nLay]);
          }
        }

        // Uncomment the next lines to use the MC truth as seed for the track
        // float xyz0[3]={xclu[nLayers-1],yclu[nLayers-1],zclu[nLayers-1]};
        // float pxyz0[3]={curPart.Px(),curPart.Py(),curPart.Pz()};
        // printf("Initialize seed at %f %f %f p = %f %f %f\n",xyz0[0],xyz0[1],xyz0[2],pxyz0[0],pxyz0[1],pxyz0[2]);
        // int sign = curPart.GetPdgCode() > 0 ? 1 : -1;
        // fitter->setSeed(xyz0,pxyz0,sign);
        hMomTrackable->Fill(momPart);
        hEtaTrackable->Fill(etaPart);
        TParticle trFit;
        NA6PTrack currTr;
        float chiFit = 0.f;
        if (fitter->computeSeed(-1, &currTr) && (chiFit = fitter->fitSeed(currTr)) >= 0.f) {
          //  if (fitter->fitTrackPoints(currTr)) {
          if (verbose)
            std::cout << "Track fit done.\n";
          fitter->addClustersToTrack(currTr);
          currTr.setChi2(chiFit);
          // fitter->propagateToZ(&currTr,zvert);
          hIsGood->Fill(1);
          int nClusters = currTr.getNHits();
          hNclu->Fill(nClusters);
          if (nClusters != minHits) {
            if (verbose)
              std::cout << "Skipping track with only " << nClusters << " hits (minHits=" << minHits << ")\n";
            continue;
          } else {
            if (verbose)
              std::cout << "Fitted track with " << nClusters << " hits\n";
          }
          chiFit = currTr.getChi2();
          if (verbose)
            printf("chi2 = %f \n", chiFit);
          auto pxyz = currTr.getPXYZ();
          float momtr = currTr.getP();
          float phitr = currTr.getPhi();
          float thetatr = currTr.getTheta();
          float etatr = currTr.getEta();
          float impparX = currTr.getX() - xvert;
          float impparY = currTr.getY() - yvert;
          hMomFitted->Fill(momPart);
          hEtaFitted->Fill(etaPart);
          hEtaReco->Fill(etatr);
          hDeltaPx->Fill(pxyz[0] - pxPart);
          hDeltaPy->Fill(pxyz[1] - pyPart);
          hDeltaPz->Fill(pxyz[2] - pzPart);
          hDeltaP->Fill(momtr - momPart);
          hDeltaPxVsP->Fill(momtr, pxyz[0] - pxPart);
          hDeltaPyVsP->Fill(momtr, pxyz[1] - pyPart);
          hDeltaPzVsP->Fill(momtr, pxyz[2] - pzPart);
          hDeltaPVsP->Fill(momtr, momtr - momPart);
          hRelDeltaPxVsP->Fill(momtr, (pxyz[0] - pxPart) / pxPart);
          hRelDeltaPyVsP->Fill(momtr, (pxyz[1] - pyPart) / pyPart);
          hRelDeltaPzVsP->Fill(momtr, (pxyz[2] - pzPart) / pzPart);
          hRelDeltaPVsP->Fill(momtr, (momtr - momPart) / momPart);
          hDeltaPhi->Fill(phitr - phiPart);
          hDeltaEta->Fill(etatr - etaPart);
          hImpParX->Fill(impparX * 1e4);
          hImpParY->Fill(impparY * 1e4);
          hImpParXVsP->Fill(momtr, impparX * 1e4);
          hImpParYVsP->Fill(momtr, impparY * 1e4);
          hChi2NDF->Fill(chiFit);
          TLorentzVector v;
          v.SetXYZM(pxyz[0], pxyz[1], pxyz[2], 0.14);
          trFit.SetMomentum(v);
          trFit.SetProductionVertex(currTr.getX(), currTr.getY(), currTr.getZ(), 0);
        } else {
          hIsGood->Fill(0);
        }
        outTr << "out"
              << "mcTr=" << curPart << "fitTr=" << trFit << "fitChi2=" << chiFit << "\n";
      }
    }
  }
  outTr.Close();

  TH1F* hImpParXMean = new TH1F("hImpParXMean", ";p (GeV/c);<Imp Par X> (#mum)", nMomBins, 0., maxP);
  TH1F* hImpParYMean = new TH1F("hImpParYMean", ";p (GeV/c);<Imp Par Y> (#mum)", nMomBins, 0., maxP);
  TH1F* hImpParXRms = new TH1F("hImpParXRms", ";p (GeV/c);rms (Imp Par X) (#mum)", nMomBins, 0., maxP);
  TH1F* hImpParYRms = new TH1F("hImpParYRms", ";p (GeV/c);rms (Imp Par Y) (#mum)", nMomBins, 0., maxP);
  TH1F* hImpParXSig = new TH1F("hImpParXSig", ";p (GeV/c);#sigma(Imp Par X) (#mum)", nMomBins, 0., maxP);
  TH1F* hImpParYSig = new TH1F("hImpParYSig", ";p (GeV/c);#sigma(Imp Par Y) (#mum)", nMomBins, 0., maxP);
  FillMeanAndRms(hImpParXVsP, hImpParXMean, hImpParXRms, hImpParXSig);
  FillMeanAndRms(hImpParYVsP, hImpParYMean, hImpParYRms, hImpParYSig);

  TH1F* hDeltaPMean = new TH1F("hDeltaPMean", ";p (GeV/c);<p_{rec}-p_{gen}> (GeV/c)", nMomBins, 0., maxP);
  TH1F* hDeltaPRms = new TH1F("hDeltaPRms", ";p (GeV/c);rms (p_{rec}-p_{gen}) (GeV/c)", nMomBins, 0., maxP);
  TH1F* hDeltaPSig = new TH1F("hDeltaPSig", ";p (GeV/c);#sigma(p_{rec}-p_{gen}) (GeV/c)", nMomBins, 0., maxP);
  TH1F* hRelDeltaPMean = new TH1F("hRelDeltaPMean", ";p (GeV/c);<(p_{rec}-p_{gen})/p_{gen}>", nMomBins, 0., maxP);
  TH1F* hRelDeltaPRms = new TH1F("hRelDeltaPRms", ";p (GeV/c);rms (p_{rec}-p_{gen})/p_{gen}", nMomBins, 0., maxP);
  TH1F* hRelDeltaPSig = new TH1F("hRelDeltaPSig", ";p (GeV/c);#sigma(p_{rec}-p_{gen})/p_{gen} (GeV/c)", nMomBins, 0., maxP);
  FillMeanAndRms(hDeltaPVsP, hDeltaPMean, hDeltaPRms, hDeltaPSig);
  FillMeanAndRms(hRelDeltaPVsP, hRelDeltaPMean, hRelDeltaPRms, hRelDeltaPSig);

  TH1F* hDeltaPxMean = new TH1F("hDeltaPxMean", ";p (GeV/c);<p_{x}^{rec}-p_{x}^{gen}> (GeV/c)", nMomBins, 0., maxP);
  TH1F* hDeltaPxRms = new TH1F("hDeltaPxRms", ";p (GeV/c);rms (p_{x}^{rec}-p_{x}^{gen}) (GeV/c)", nMomBins, 0., maxP);
  TH1F* hDeltaPxSig = new TH1F("hDeltaPxSig", ";p (GeV/c);#sigma(p_{x}^{rec}-p_{x}^{gen}) (GeV/c)", nMomBins, 0., maxP);
  TH1F* hRelDeltaPxMean = new TH1F("hRelDeltaPxMean", ";p (GeV/c);<(p_{x}^{rec}-p_{x}^{gen})/p_{x}^{gen}>", nMomBins, 0., maxP);
  TH1F* hRelDeltaPxRms = new TH1F("hRelDeltaPxRms", ";p (GeV/c);rms (p_{x}^{rec}-p_{x}^{gen})/p_{x}^{gen}", nMomBins, 0., maxP);
  TH1F* hRelDeltaPxSig = new TH1F("hRelDeltaPxSig", ";p (GeV/c);#sigma(p_{x}^{rec}-p_{x}^{gen})/p_{x}^{gen} (GeV/c)", nMomBins, 0., maxP);
  FillMeanAndRms(hDeltaPxVsP, hDeltaPxMean, hDeltaPxRms, hDeltaPxSig);
  FillMeanAndRms(hRelDeltaPxVsP, hRelDeltaPxMean, hRelDeltaPxRms, hRelDeltaPxSig);

  TH1F* hDeltaPyMean = new TH1F("hDeltaPyMean", ";p (GeV/c);<p_{y}^{rec}-p_{x}^{gen}> (GeV/c)", nMomBins, 0., maxP);
  TH1F* hDeltaPyRms = new TH1F("hDeltaPyRms", ";p (GeV/c);rms (p_{y}^{rec}-p_{x}^{gen}) (GeV/c)", nMomBins, 0., maxP);
  TH1F* hDeltaPySig = new TH1F("hDeltaPySig", ";p (GeV/c);#sigma(p_{y}^{rec}-p_{x}^{gen}) (GeV/c)", nMomBins, 0., maxP);
  TH1F* hRelDeltaPyMean = new TH1F("hRelDeltaPyMean", ";p (GeV/c);<(p_{y}^{rec}-p_{y}^{gen})/p_{y}^{gen}>", nMomBins, 0., maxP);
  TH1F* hRelDeltaPyRms = new TH1F("hRelDeltaPyRms", ";p (GeV/c);rms (p_{y}^{rec}-p_{y}^{gen})/p_{y}^{gen}", nMomBins, 0., maxP);
  TH1F* hRelDeltaPySig = new TH1F("hRelDeltaPySig", ";p (GeV/c);#sigma(p_{y}^{rec}-p_{y}^{gen})/p_{y}^{gen} (GeV/c)", nMomBins, 0., maxP);
  FillMeanAndRms(hDeltaPyVsP, hDeltaPyMean, hDeltaPyRms, hDeltaPySig);
  FillMeanAndRms(hRelDeltaPyVsP, hRelDeltaPyMean, hRelDeltaPyRms, hRelDeltaPySig);

  TH1F* hDeltaPzMean = new TH1F("hDeltaPzMean", ";p (GeV/c);<p_{z}^{rec}-p_{z}^{gen}> (GeV/c)", nMomBins, 0., maxP);
  TH1F* hDeltaPzRms = new TH1F("hDeltaPzRms", ";p (GeV/c);rms (p_{z}^{rec}-p_{z}^{gen}) (GeV/c)", nMomBins, 0., maxP);
  TH1F* hDeltaPzSig = new TH1F("hDeltaPzSig", ";p (GeV/c);#sigma(p_{z}^{rec}-p_{z}^{gen}) (GeV/c)", nMomBins, 0., maxP);
  TH1F* hRelDeltaPzMean = new TH1F("hRelDeltaPzMean", ";p (GeV/c);<(p_{z}^{rec}-p_{z}^{gen})/p_{z}^{gen}>", nMomBins, 0., maxP);
  TH1F* hRelDeltaPzRms = new TH1F("hRelDeltaPzRms", ";p (GeV/c);rms (p_{z}^{rec}-p_{z}^{gen})/p_{z}^{gen}", nMomBins, 0., maxP);
  TH1F* hRelDeltaPzSig = new TH1F("hRelDeltaPzSig", ";p (GeV/c);#sigma(p_{z}^{rec}-p_{z}^{gen})/p_{z}^{gen} (GeV/c)", nMomBins, 0., maxP);
  FillMeanAndRms(hDeltaPzVsP, hDeltaPzMean, hDeltaPzRms, hDeltaPzSig);
  FillMeanAndRms(hRelDeltaPzVsP, hRelDeltaPzMean, hRelDeltaPzRms, hRelDeltaPzSig);

  TCanvas* ct = new TCanvas("ct", "", 1000, 800);
  ct->Divide(2, 2);
  ct->cd(1);
  hIsGood->Draw();
  ct->cd(2);
  hNclu->Draw();
  ct->cd(3);
  hChi2NDF->Draw();

  float deltaMaxImp = 60000.;
  if (useVT)
    deltaMaxImp = 500.;
  TCanvas* cip = new TCanvas("cip", "", 1200, 800);
  cip->Divide(2, 2);
  cip->cd(1);
  gPad->SetTickx();
  gPad->SetTicky();
  hImpParXMean->SetMarkerStyle(25);
  hImpParXMean->SetMarkerColor(1);
  hImpParXMean->SetLineColor(1);
  hImpParXMean->SetMinimum(-deltaMaxImp);
  hImpParXMean->SetMaximum(deltaMaxImp);
  hImpParXMean->Draw("P");
  cip->cd(2);
  gPad->SetTickx();
  gPad->SetTicky();
  hImpParXSig->SetMarkerStyle(25);
  hImpParXSig->SetMarkerColor(1);
  hImpParXSig->SetLineColor(1);
  hImpParXSig->SetMinimum(0);
  hImpParXSig->SetMaximum(deltaMaxImp);
  hImpParXSig->Draw("P");
  cip->cd(3);
  gPad->SetTickx();
  gPad->SetTicky();
  hImpParYMean->SetMarkerStyle(25);
  hImpParYMean->SetMarkerColor(1);
  hImpParYMean->SetLineColor(1);
  hImpParYMean->SetMinimum(-deltaMaxImp);
  hImpParYMean->SetMaximum(deltaMaxImp);
  hImpParYMean->Draw("P");
  cip->cd(4);
  gPad->SetTickx();
  gPad->SetTicky();
  hImpParYSig->SetMarkerStyle(25);
  hImpParYSig->SetMarkerColor(1);
  hImpParYSig->SetLineColor(1);
  hImpParYSig->SetMinimum(0);
  hImpParYSig->SetMaximum(deltaMaxImp);
  hImpParYSig->Draw("P");

  float maxDelP = 0.75;
  if (useVT)
    maxDelP = 0.05;
  TCanvas* cmom = new TCanvas("cmom", "", 1400, 900);
  cmom->Divide(3, 3);
  cmom->cd(1);
  gPad->SetTickx();
  gPad->SetTicky();
  hDeltaPxMean->SetMarkerStyle(25);
  hDeltaPxMean->SetMarkerColor(1);
  hDeltaPxMean->SetLineColor(1);
  hDeltaPxMean->SetMinimum(-maxDelP);
  hDeltaPxMean->SetMaximum(maxDelP);
  hDeltaPxMean->Draw("P");
  cmom->cd(2);
  gPad->SetTickx();
  gPad->SetTicky();
  hDeltaPxSig->SetMarkerStyle(25);
  hDeltaPxSig->SetMarkerColor(1);
  hDeltaPxSig->SetLineColor(1);
  hDeltaPxSig->SetMinimum(0);
  hDeltaPxSig->SetMaximum(maxDelP);
  hDeltaPxSig->Draw("P");
  cmom->cd(3);
  gPad->SetTickx();
  gPad->SetTicky();
  hRelDeltaPxSig->SetMarkerStyle(25);
  hRelDeltaPxSig->SetMarkerColor(1);
  hRelDeltaPxSig->SetLineColor(1);
  hRelDeltaPxSig->SetMinimum(0);
  hRelDeltaPxSig->SetMaximum(maxDelP * 2.);
  hRelDeltaPxSig->Draw("P");
  cmom->cd(4);
  gPad->SetTickx();
  gPad->SetTicky();
  hDeltaPyMean->SetMarkerStyle(25);
  hDeltaPyMean->SetMarkerColor(1);
  hDeltaPyMean->SetLineColor(1);
  hDeltaPyMean->SetMinimum(-maxDelP);
  hDeltaPyMean->SetMaximum(maxDelP);
  hDeltaPyMean->Draw("P");
  cmom->cd(5);
  gPad->SetTickx();
  gPad->SetTicky();
  hDeltaPySig->SetMarkerStyle(25);
  hDeltaPySig->SetMarkerColor(1);
  hDeltaPySig->SetLineColor(1);
  hDeltaPySig->SetMinimum(0);
  hDeltaPySig->SetMaximum(maxDelP);
  hDeltaPySig->Draw("P");
  cmom->cd(6);
  gPad->SetTickx();
  gPad->SetTicky();
  hRelDeltaPySig->SetMarkerStyle(25);
  hRelDeltaPySig->SetMarkerColor(1);
  hRelDeltaPySig->SetLineColor(1);
  hRelDeltaPySig->SetMinimum(0);
  hRelDeltaPySig->SetMaximum(maxDelP * 2.);
  hRelDeltaPySig->Draw("P");
  cmom->cd(7);
  gPad->SetTickx();
  gPad->SetTicky();
  hDeltaPzMean->SetMarkerStyle(25);
  hDeltaPzMean->SetMarkerColor(1);
  hDeltaPzMean->SetLineColor(1);
  hDeltaPzMean->SetMinimum(-maxDelP);
  hDeltaPzMean->SetMaximum(maxDelP);
  hDeltaPzMean->Draw("P");
  cmom->cd(8);
  gPad->SetTickx();
  gPad->SetTicky();
  hDeltaPzSig->SetMarkerStyle(25);
  hDeltaPzSig->SetMarkerColor(1);
  hDeltaPzSig->SetLineColor(1);
  hDeltaPzSig->SetMinimum(0);
  hDeltaPzSig->SetMaximum(maxDelP * 4.);
  hDeltaPzSig->Draw("P");
  cmom->cd(9);
  hRelDeltaPzSig->SetMarkerStyle(25);
  hRelDeltaPzSig->SetMarkerColor(1);
  hRelDeltaPzSig->SetLineColor(1);
  hRelDeltaPzSig->SetMinimum(0);
  hRelDeltaPzSig->SetMaximum(maxDelP);
  hRelDeltaPzSig->Draw("P");

  TCanvas* ce = new TCanvas("ce", "", 1400, 800);
  ce->Divide(2, 2);
  ce->cd(1);
  hMomGen->SetLineColor(kGray);
  hMomGen->SetLineWidth(2);
  hMomGen->Draw();
  hMomTrackable->SetLineColor(1);
  hMomTrackable->SetLineWidth(2);
  hMomTrackable->Draw("sames");
  hMomFitted->SetLineColor(kGreen + 1);
  hMomFitted->SetLineWidth(2);
  hMomFitted->Draw("sames");
  ce->cd(2);
  TH1F* hEffMom = (TH1F*)hMomFitted->Clone("hEffMom");
  hEffMom->Divide(hMomFitted, hMomTrackable, 1., 1., "B");
  hEffMom->SetMaximum(1.05);
  hEffMom->SetStats(0);
  hEffMom->Draw("");
  ce->cd(3);
  hEtaGen->SetLineColor(kGray);
  hEtaGen->SetLineWidth(2);
  hEtaGen->Draw();
  hEtaTrackable->SetLineColor(1);
  hEtaTrackable->SetLineWidth(2);
  hEtaTrackable->Draw("sames");
  hEtaFitted->SetLineColor(kGreen + 1);
  hEtaFitted->SetLineWidth(2);
  hEtaFitted->Draw("sames");
  ce->cd(4);
  TH1F* hEffEta = (TH1F*)hEtaFitted->Clone("hEffEta");
  hEffEta->Divide(hEtaFitted, hEtaTrackable, 1., 1., "B");
  hEffEta->SetStats(0);
  hEffEta->SetMaximum(1.05);
  hEffEta->Draw("");

  TFile* outFFit = new TFile("TrackingFastFit.root", "recreate");
  hImpParXVsP->Write();
  hImpParYVsP->Write();
  hDeltaPxVsP->Write();
  hDeltaPyVsP->Write();
  hDeltaPzVsP->Write();
  hDeltaPVsP->Write();
  hRelDeltaPxVsP->Write();
  hRelDeltaPyVsP->Write();
  hRelDeltaPzVsP->Write();
  hRelDeltaPVsP->Write();
  hEtaGen->Write();
  hEtaReco->Write();
  outFFit->Close();
}

// Wrapper functions for backward compatibility and ease of use

// Original function signature for Vertex Telescope hits
void runFastFitOnVerTelHits(int firstEv = 0,
                            int lastEv = 99999999,
                            float clures = 5.e-4,
                            const char* dirSimu = ".",
                            int minHits = 5,
                            bool verbose = false)
{
  runFastFitOnHitsTemplate<NA6PVerTelHit>(firstEv, lastEv, clures, clures, dirSimu,
                                          "HitsVerTel.root", "hitsVerTel", "VerTel",
                                          5, minHits, verbose, true);
}

// New function for Muon Spectrometer Modular hits
void runFastFitOnMuonSpecHits(int firstEv = 0,
                              int lastEv = 99999999,
                              float cluresx = 500.e-4,
                              float cluresy = 1000.e-4,
                              const char* dirSimu = ".",
                              int minHits = 6,
                              bool verbose = false)
{
  runFastFitOnHitsTemplate<NA6PMuonSpecHit>(firstEv, lastEv, cluresx, cluresy, dirSimu,
                                                   "HitsMuonSpec.root", "hitsMuonSpec", "MuonSpec",
                                                   6, minHits, verbose, false);
}
