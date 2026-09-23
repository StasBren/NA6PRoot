#if !defined(__CINT__) || defined(__MAKECINT__)
#include "NA6PVerTelHit.h"
#include "NA6PMuonSpecHit.h"
#include "NA6PMCGenHeader.h"
#include "NA6PMCEventHeader.h"
#include "ConfigurableParam.h"
#include "NA6PLayoutParam.h"
#include <TParticle.h>
#include <TParticlePDG.h>
#include <TTree.h>
#include <TChain.h>
#include <TFile.h>
#include <TKey.h>
#include <fairlogger/Logger.h>
#include <cstdio>
#include <bitset>

#include <iostream>
#include <fstream>
#include <nlohmann/json.hpp>
#endif

using json = nlohmann::json;

// Structure to hold the extracted surface data
struct SurfaceData {
  uint64_t geo_id = 0;
  float z = 0.0f;
  float range_x[2] = {0.0f, 0.0f};
  float range_y[2] = {0.0f, 0.0f};
};

// Function to extract specific fields from surface entries
std::vector<SurfaceData> extractSurfaceData(const std::string& filename)
{
  std::ifstream file(filename);

  // Check if file opened successfully
  if (!file.is_open()) {
    throw std::runtime_error("Could not open file: " + filename);
  }

  std::vector<SurfaceData> surfaces;
  json jsonData;

  try {
    file >> jsonData;

    if (!jsonData.contains("Surfaces") || !jsonData["Surfaces"].contains("entries")) {
      std::cerr << "No Surfaces->entries found in JSON" << std::endl;
      return surfaces;
    }

    for (const auto& entry : jsonData["Surfaces"]["entries"]) {
      // Extract sensitive (optional)
      if (entry.contains("sensitive")) {
        SurfaceData surface;
        // Extract from nested "value" object
        if (entry.contains("value")) {
          const auto& value = entry["value"];
          // Extract geo_id
          if (value.contains("geo_id")) {
            surface.geo_id = value["geo_id"];
          }

          // Surface position and ranges
          if (value.contains("transform")) {
            surface.z = value["transform"]["translation"][2];
            surface.range_x[0] = static_cast<float>(value["material"]["binUtility"]["binningdata"][0]["min"]) + static_cast<float>(value["transform"]["translation"][0]);
            surface.range_x[1] = static_cast<float>(value["material"]["binUtility"]["binningdata"][0]["max"]) + static_cast<float>(value["transform"]["translation"][0]);
            surface.range_y[0] = static_cast<float>(value["material"]["binUtility"]["binningdata"][1]["min"]) + static_cast<float>(value["transform"]["translation"][1]);
            surface.range_y[1] = static_cast<float>(value["material"]["binUtility"]["binningdata"][1]["max"]) + static_cast<float>(value["transform"]["translation"][1]);
            surfaces.push_back(surface);
          }
        }
      }
    }
  } catch (const json::parse_error& e) {
    throw std::runtime_error("JSON parse error: " + std::string(e.what()));
  }

  file.close();
  return surfaces;
}

uint64_t getACTSgeomID(float x, float y, float z, const std::vector<SurfaceData>& surfaces)
{
  auto getGeoID = [](float x, float y, float z, const std::vector<SurfaceData>& surfaces) -> uint64_t {
    for (const auto& surface : surfaces) {
      if (x >= surface.range_x[0] && x <= surface.range_x[1] &&
          y >= surface.range_y[0] && y <= surface.range_y[1] &&
          abs(z - surface.z) < 0.1f) {
        return surface.geo_id;
      }
    }
    return 0; // Return 0 if no matching surface is found
  };
  return getGeoID(x, y, z, surfaces);
}

TChain* loadUserChain(const char* inpData, const char* chName, bool recursive = false);
void addRootFileDirectory(TChain* chain, const char* flName, const char* chName, bool recursive = false);

// Template function per convertire gli hit generici
template <typename MSHitType>
void convert2ACTSImpl(int iev,
                      const std::vector<TParticle>& mcParts,
                      const std::vector<NA6PVerTelHit>& vtHits,
                      const std::vector<MSHitType>& msHits,
                      const std::vector<SurfaceData>& actsSurfaces);

// Wrapper functions per entrambi i tipi
void convert2ACTS(int iev,
                  const std::vector<TParticle>& mcParts,
                  const std::vector<NA6PVerTelHit>& vtHits,
                  const std::vector<NA6PMuonSpecHit>& msHits,
                  const std::vector<SurfaceData>& actsSurfaces)
{
  convert2ACTSImpl(iev, mcParts, vtHits, msHits, actsSurfaces);
}

void convert2ACTS(int iev,
                  const std::vector<TParticle>& mcParts,
                  const std::vector<NA6PVerTelHit>& vtHits,
                  const std::vector<NA6PMuonSpecHit>& msHits,
                  const std::vector<SurfaceData>& actsSurfaces)
{
  convert2ACTSImpl(iev, mcParts, vtHits, msHits, actsSurfaces);
}

void convert2ACTS(const std::string& dirname = "./",
                  const std::string& fnameMCKin = "MCKine.root",
                  const std::string& fnameVTHits = "HitsVerTel.root",
                  const std::string& fnameMSHits = "HitsMuonSpec.root",
                  const std::string& geometryFile = "geometry-map.json",
                  const std::string& confOpts = "",
                  bool useModularSpec = true)
{
  auto actsSurfaces = extractSurfaceData(geometryFile);

  std::string dirnameL = dirname;
  if (dirnameL.empty()) {
    dirnameL = ".";
  } else if (dirnameL.back() != '/') {
    dirnameL += '/';
  }
  if (!confOpts.empty()) {
    na6p::conf::ConfigurableParam::updateFromString(confOpts);
  }

  NA6PMCEventHeader mcHeader, *mcHeaderPtr = &mcHeader;
  std::vector<TParticle> mcParts, *mcPartsPtr = &mcParts;
  std::vector<NA6PVerTelHit> vtHits, *vtHitsPtr = &vtHits;

  TChain* treeKin = loadUserChain(fmt::format("{}{}", dirnameL, fnameMCKin).c_str(), "mckine");
  TChain* treeVTH = loadUserChain(fmt::format("{}{}", dirnameL, fnameVTHits).c_str(), "hitsVerTel");

  // Determina quale branch name usare
  const char* msBranchName = useModularSpec ? "hitsMuonSpec" : "hitsMuonSpec";
  const char* msTreeName = useModularSpec ? "MuonSpec" : "MuonSpec";

  TChain* treeMSH = loadUserChain(fmt::format("{}{}", dirnameL, fnameMSHits).c_str(), msBranchName);

  int nent = 0;
  if (!treeKin || !treeVTH || !treeMSH || !(nent = treeKin->GetEntries()) ||
      treeVTH->GetEntries() != nent || treeMSH->GetEntries() != nent) {
    LOGP(fatal, "Problem in opening input file or empty input");
  }

  treeKin->SetBranchAddress("header", &mcHeaderPtr);
  treeKin->SetBranchAddress("tracks", &mcPartsPtr);
  treeVTH->SetBranchAddress("VerTel", &vtHitsPtr);

  if (useModularSpec) {
    std::vector<NA6PMuonSpecHit> msHits, *msHitsPtr = &msHits;
    treeMSH->SetBranchAddress(msTreeName, &msHitsPtr);

    for (int iev = 0; iev < nent; iev++) {
      treeKin->GetEntry(iev);
      treeVTH->GetEntry(iev);
      treeMSH->GetEntry(iev);
      LOGP(info, "Event#{} Vtx:[{:.2f},{:.2f},{:.2f}] {} Tracks ({} primaries), Hits: VT: {} MS: {}",
           iev, mcHeader.getVX(), mcHeader.getVY(), mcHeader.getVZ(),
           mcHeader.getNTracks(), mcHeader.getNPrimaries(), vtHits.size(), msHits.size());
      convert2ACTS(iev, mcParts, vtHits, msHits, actsSurfaces);
    }
  } else {
    std::vector<NA6PMuonSpecHit> msHits, *msHitsPtr = &msHits;
    treeMSH->SetBranchAddress(msTreeName, &msHitsPtr);

    for (int iev = 0; iev < nent; iev++) {
      treeKin->GetEntry(iev);
      treeVTH->GetEntry(iev);
      treeMSH->GetEntry(iev);
      LOGP(info, "Event#{} Vtx:[{:.2f},{:.2f},{:.2f}] {} Tracks ({} primaries), Hits: VT: {} MS: {}",
           iev, mcHeader.getVX(), mcHeader.getVY(), mcHeader.getVZ(),
           mcHeader.getNTracks(), mcHeader.getNPrimaries(), vtHits.size(), msHits.size());
      convert2ACTS(iev, mcParts, vtHits, msHits, actsSurfaces);
    }
  }
}

TChain* loadUserChain(const char* inpData, const char* chName, bool recursive)
{
  if (!inpData) {
    return nullptr;
  }
  TChain* chain = new TChain(chName);
  TString inpDtStr = inpData;
  if (inpDtStr.EndsWith(".root")) {
    addRootFileDirectory(chain, inpData, chName, recursive);
  } else {
    //
    ifstream inpf(inpData);
    if (!inpf.good()) {
      LOGP(error, "Failed on input filename {}", inpData);
      return nullptr;
    }
    //
    TString flName;
    flName.ReadLine(inpf);
    while (!flName.IsNull()) {
      flName = flName.Strip(TString::kBoth, ' ');
      if (flName.BeginsWith("//") || flName.BeginsWith("#")) {
        flName.ReadLine(inpf);
        continue;
      }
      flName = flName.Strip(TString::kBoth, ',');
      flName = flName.Strip(TString::kBoth, '"');
      addRootFileDirectory(chain, flName.Data(), chName, recursive);
      flName.ReadLine(inpf);
    }
  }
  //
  int n = chain->GetEntries();
  if (n < 1) {
    LOGP(error, "Obtained chain is empty");
    return nullptr;
  } else
    LOGP(info, "Opened {} chain with {} entries", chName, n);
  return chain;
}

void addRootFileDirectory(TChain* chain, const char* flName, const char* chName, bool recursive)
{
  if (!recursive) {
    printf("Adding %s\n", flName);
    chain->AddFile(flName);
  } else {
    TFile* f = TFile::Open(flName);
    if (!f || f->IsZombie()) {
      LOGP(error, "Failed to open {}", flName);
      delete f;
      return;
    }
    auto* lstD = f->GetListOfKeys();
    TIter iter(lstD->MakeIterator());
    while (TObject* obj = iter()) {
      TKey* theKey = (TKey*)obj;
      TString kname = theKey->GetClassName();
      if (!kname.BeginsWith("TDirectory")) {
        continue;
      }
      if (f->Get(Form("%s/%s", theKey->GetName(), chName))) {
        LOGP(info, "Adding {}/{}", flName, theKey->GetName());
        chain->AddFile(flName, TTree::kMaxEntries, Form("%s/%s", theKey->GetName(), chName));
      }
    }
    delete f;
  }
}

template <typename MSHitType>
void convert2ACTSImpl(int iev,
                      const std::vector<TParticle>& mcParts,
                      const std::vector<NA6PVerTelHit>& vtHits,
                      const std::vector<MSHitType>& msHits,
                      const std::vector<SurfaceData>& actsSurfaces)
{
  std::string eventFileName = fmt::format("event{:09}-particles.csv", iev);
  std::string primFileName = fmt::format("event{:09}-primaries.csv", iev);
  std::string hitsFileName = fmt::format("event{:09}-hits.csv", iev);
  std::FILE* fpcsvKin = std::fopen(eventFileName.c_str(), "w");
  std::FILE* fpcsvPrim = std::fopen(primFileName.c_str(), "w");

  fprintf(fpcsvKin, "particle_id_pv,particle_id_sv,particle_id_part,particle_id_gen,particle_id_subpart,particle_type,process,vx,vy,vz,vt,px,py,pz,m,q,nhits,nhits_vs,nhits_ms\n");
  fprintf(fpcsvPrim, "particle_id_pv,particle_id_sv,particle_id_part,particle_id_gen,particle_id_subpart,particle_type,process,vx,vy,vz,vt,px,py,pz,m,q,nhits,nhits_vs,nhits_ms\n");
  int decCount = 0;
  std::vector<uint64_t> doneBC;
  doneBC.resize(mcParts.size());
  std::vector<uint64_t> donePDG;
  donePDG.resize(mcParts.size());
  std::vector<uint64_t> secV;
  secV.resize(mcParts.size());
  printf("Loop on %d particles\n", (int)mcParts.size());
  for (int ipartTop = 0; ipartTop < (int)mcParts.size(); ipartTop++) {
    // clang-format off
    auto getBarCode = [&iev, &ipartTop](int gen, int sub, int isv)
    {
      uint64_t b =
	( (uint64_t(iev+1)    & ((0x1<<12)-1)) << ( 12 + 16 + 8 + 16 ) ) |
	( (uint64_t(isv)      & ((0x1<<12)-1)) << (      16 + 8 + 16 ) ) |
	( (uint64_t(ipartTop) & ((0x1<<16)-1)) << (           8 + 16 ) ) |
	( (uint64_t(gen)      & ((0x1<<8)-1))  << (               16 ) ) |
	uint64_t(sub);
      return b;
    };
    // clang-format on

    auto storePart = [iev, ipartTop, fpcsvKin, fpcsvPrim, &doneBC, &donePDG, &secV, &decCount, &mcParts, &vtHits, &msHits, getBarCode](int ip, int gen, int sub, int imoth, const auto& self) -> void {
      if (doneBC[ip]) {
        //	printf("Particle %d already stored\n",ip);
        return;
      }
      const auto& part = mcParts[ip];
      const auto pdgPart = part.GetPDG();
      double fcmtomm = 10.f;

      if (!pdgPart) {
        LOGP(info, "Skipping unknown particle {}", part.GetPdgCode());
        return;
      }
      doneBC[ip] = getBarCode(gen, sub, imoth);
      donePDG[ip] = std::abs(part.GetPdgCode());
      // Count hits for this particle
      int vtHitCount = 0;
      int msHitCount = 0;
      int totalHitCount = 0;

      //  Count vertex telescope hits
      for (const auto& hit : vtHits) {
        if (hit.getTrackID() == ip) {
          vtHitCount++;
        }
      }

      // Count muon spectrometer hits
      for (const auto& hit : msHits) {
        if (hit.getTrackID() == ip) {
          msHitCount++;
        }
      }

      totalHitCount = vtHitCount + msHitCount;

      // Modified fprintf to include hit counts - adding totalHitCount, vtHitCount, msHitCount
      fprintf(fpcsvKin, "%i,%i,%i,%i,%i,%d,%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%d,%d,%d\n",
              iev + 1, imoth, ipartTop, gen, sub, std::abs(part.GetPdgCode()), 0,
              part.Vx() * fcmtomm, part.Vy() * fcmtomm, part.Vz() * fcmtomm, part.T(),
              part.Px(), part.Py(), part.Pz(), part.GetCalcMass(), pdgPart->Charge() / 3,
              totalHitCount, vtHitCount, msHitCount);

      if (part.IsPrimary()) {
        fprintf(fpcsvPrim, "%i,%i,%i,%i,%i,%d,%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%d,%d,%d\n",
                iev + 1, imoth, ipartTop, gen, sub, std::abs(part.GetPdgCode()), 0,
                part.Vx() * fcmtomm, part.Vy() * fcmtomm, part.Vz() * fcmtomm, part.T(),
                part.Px(), part.Py(), part.Pz(), part.GetCalcMass(), pdgPart->Charge() / 3,
                totalHitCount, vtHitCount, msHitCount);
      }

      int dtF = part.GetFirstDaughter();
      int dtL = part.GetLastDaughter();
      if (dtF > -1) {
        gen++;
        decCount++;
        secV[ip] = decCount;
        int subG = 0;
        for (int i = dtF; i <= dtL; i++) {
          self(i, gen, subG++, secV[ip], self);
        }
      }
    };
    storePart(ipartTop, 0, 0, 0, storePart);
  }
  std::fclose(fpcsvKin);
  std::fclose(fpcsvPrim);

  std::FILE* fpcsvHits = std::fopen(hitsFileName.c_str(), "w");
  fprintf(fpcsvHits, "particle_id_pv,particle_id_sv,particle_id_part,particle_id_gen,particle_id_subpart,geometry_id,tx,ty,tz,tt,tpx,tpy,tpz,te,deltapx,deltapy,deltapz,deltae,index\n");
  std::map<uint64_t, int> hitPerPart;
  double fcmtomm = 10.f;
  auto decodeBarCode = [](uint64_t b) {
    uint64_t sub = b & ((1ULL << 16) - 1);
    uint64_t gen = (b >> 16) & ((1ULL << 8) - 1);
    uint64_t ipartTop = (b >> 24) & ((1ULL << 16) - 1);
    uint64_t isv = (b >> 40) & ((1ULL << 12) - 1);
    uint64_t iev1 = (b >> 52) & ((1ULL << 12) - 1); // this is (iev + 1)

    return std::make_tuple(iev1, isv, ipartTop, gen, sub);
  };
  for (const auto& hit : vtHits) {
    double m = mcParts[hit.getTrackID()].GetCalcMass();
    double en = std::sqrt(hit.getPXIn() * hit.getPXIn() + hit.getPYIn() * hit.getPYIn() + hit.getPZIn() * hit.getPZIn() + m * m);
    auto particle_id = doneBC[hit.getTrackID()];
    if (particle_id == 0) {
      int idP = hit.getTrackID();
      const auto& part = mcParts[idP];
      const auto pdgPart = part.GetPDG();
      printf("Hit w/o particle: %lu %d pdgPart=%x, pdgcode=%d\n", particle_id, hit.getTrackID(), pdgPart, part.GetPdgCode());
    }
    uint64_t geometry_id = getACTSgeomID(hit.getX() * fcmtomm, hit.getY() * fcmtomm, hit.getZ() * fcmtomm, actsSurfaces);
    if (geometry_id) {

      if (hitPerPart.find(particle_id) == hitPerPart.end()) {
        hitPerPart[particle_id] = 0;
      } else {
        hitPerPart[particle_id]++;
      }
      auto [iev1, isv, ipartTop, gen, sub] = decodeBarCode(particle_id);

      fprintf(fpcsvHits, "%lu,%lu,%lu,%lu,%lu,%lu,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%i\n", iev1, isv, ipartTop, gen, sub, geometry_id,
              hit.getX() * fcmtomm, hit.getY() * fcmtomm, hit.getZ() * fcmtomm, hit.getTime(),
              hit.getPXIn(), hit.getPYIn(), hit.getPZIn(), en,
              hit.getPXIn() - hit.getPXOut(), hit.getPYIn() - hit.getPYOut(), hit.getPZIn() - hit.getPZOut(), hit.getHitValue(), hitPerPart[particle_id]);
    }
  }
  for (const auto& hit : msHits) {
    double m = mcParts[hit.getTrackID()].GetCalcMass();
    double en = std::sqrt(hit.getPXIn() * hit.getPXIn() + hit.getPYIn() * hit.getPYIn() + hit.getPZIn() * hit.getPZIn() + m * m);
    auto particle_id = doneBC[hit.getTrackID()];

    uint64_t geometry_id = getACTSgeomID(hit.getX() * fcmtomm, hit.getY() * fcmtomm, hit.getZ() * fcmtomm, actsSurfaces);
    if (geometry_id) {
      if (hitPerPart.find(particle_id) == hitPerPart.end()) {
        hitPerPart[particle_id] = 0;
      } else {
        hitPerPart[particle_id]++;
      }
      auto [iev1, isv, ipartTop, gen, sub] = decodeBarCode(particle_id);

      fprintf(fpcsvHits, "%lu,%lu,%lu,%lu,%lu,%lu,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%i\n", iev1, isv, ipartTop, gen, sub, geometry_id,
              hit.getX() * fcmtomm, hit.getY() * fcmtomm, hit.getZ() * fcmtomm, hit.getTime(),
              hit.getPXIn(), hit.getPYIn(), hit.getPZIn(), en,
              hit.getPXIn() - hit.getPXOut(), hit.getPYIn() - hit.getPYOut(), hit.getPZIn() - hit.getPZOut(), hit.getHitValue(), hitPerPart[particle_id]);
    }
  }
  std::fclose(fpcsvHits);
}