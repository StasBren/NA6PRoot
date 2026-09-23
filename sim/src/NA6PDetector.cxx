// NA6PCCopyright
#include "NA6PDetector.h"
#include "NA6PModule.h"
#include "NA6PTGeoHelper.h"
#include "NA6PDipoleVT.h"
#include "NA6PDipoleMS.h"
#include "NA6PMagnetsGDML.h"
#include "NA6PTarget.h"
#include "NA6PVerTel.h"
#include "NA6PAbsorber.h"
#include "NA6PMuonSpec.h"
#include "NA6PLayoutParam.h"
#include "ConfigurableParam.h"
#include "StringUtils.h"
#include <TGeoManager.h>
#include <TColor.h>
#include <TSystem.h>
#include <fairlogger/Logger.h>

NA6PDetector::NA6PDetector()
{
  const auto& param = NA6PLayoutParam::Instance();

  // Declare modules
  if (!param.use_gdml_magnets) {
    // Analytic magnets from C++ (default behaviour)
    addModule(new NA6PDipoleVT());
    addModule(new NA6PDipoleMS());
  } else {
    // GDML-based magnets: add the GDML magnet module
    addModule(new NA6PMagnetsGDML());
  }

  addModule(new NA6PTarget());
  addModule(new NA6PVerTel());
  addModule(new NA6PAbsorber());
  addModule(new NA6PMuonSpec());
}

void NA6PDetector::createCommonMaterials()
{
  auto& matPool = NA6PTGeoHelper::instance().getMatPool();
  if (matPool.find("Air") == matPool.end()) {
    auto mixt = new TGeoMixture("Air", 2, 0.001);
    mixt->AddElement(new TGeoElement("N", "Nitrogen", 7, 14.01), 0.78);
    mixt->AddElement(new TGeoElement("O", "Oxygen", 8, 16.00), 0.22);
    matPool["Air"] = mixt;
    NA6PTGeoHelper::instance().addMedium("Air");
  }
}

void NA6PDetector::addModule(NA6PModule* m)
{
  if (!m) {
    LOGP(fatal, "null pointer is provided for a module");
  }
  if (getModule(m->getName())) {
    LOGP(fatal, "module {} was already defined", m->getName());
  }
  if (m->isActive()) {
    for (const auto mo : mModulesVec) {
      if (mo->getActiveID() == m->getActiveID()) {
        LOGP(fatal, "Redefinition of ActiveModuleID:{} {}, previosly defined module: {}", m->getActiveID(), m->getName(), mo->getName());
      }
    }
  }
  mModulesMap[m->getName()] = m;
  m->setID((int)mModulesVec.size() + 1); // ID should be counted from 1!!!
  mModulesVec.push_back(m);
  if (m->isActive()) {
    mActiveModulesVec.push_back(m);
  }
  LOGP(info, "Adding module:{} {} ActiveID:{}", m->getID(), m->getName(), m->getActiveID());
}

void NA6PDetector::createGeometry(const std::string& name)
{
  const auto& param = NA6PLayoutParam::Instance();
  TGeoManager* geom = nullptr;
  TGeoVolume* world = nullptr;

  if (!param.use_gdml_magnets) {
    // === Analytic geometry (original): build from C++ code ===
    geom = new TGeoManager(name.c_str(), (name + " TGeometry").c_str());
    createCommonMaterials();
    world = geom->MakeBox("World",
                          NA6PTGeoHelper::instance().getMedium("Air"),
                          1000., 1000., 1000.);
    geom->SetTopVolume(world);

    // Build modules (analytic magnets and other detectors)
    for (auto m : mModulesVec) {
      m->createMaterials();
      m->createGeometry(world);
    }
  } else {
    // === GDML-based magnets: import BothMagnets.gdml ===
    // Get GDML module (it was added in constructor)
    NA6PMagnetsGDML* gdmlMagnetModule = nullptr;
    for (auto m : mModulesVec) {
      if (auto* gdmlMod = dynamic_cast<NA6PMagnetsGDML*>(m)) {
        gdmlMagnetModule = gdmlMod;
        break;
      }
    }
    if (!gdmlMagnetModule) {
      LOGP(fatal, "GDML magnet mode is enabled, but NA6PMagnetsGDML module not found");
    }

    std::string gdmlPath = gSystem->ExpandPathName(param.gdmlMagPath.c_str());
    LOGP(info, "Importing geometry from GDML file: {}", gdmlPath);
    TGeoManager::Import(gdmlPath.c_str());
    geom = gGeoManager;
    if (!geom) {
      LOGP(fatal, "Failed to create TGeoManager from GDML file {}", gdmlPath);
    }
    // Rename geometry to match the analytic case
    geom->SetName(name.c_str());
    geom->SetTitle((name + " TGeometry").c_str());

    createCommonMaterials();
    world = geom->GetTopVolume();
    if (!world) {
      LOGP(fatal, "No top volume found after importing GDML file {}", gdmlPath);
    }

    // Create materials for all modules BEFORE assigning them to GDML volumes
    for (auto m : mModulesVec) {
      m->createMaterials();
    }

    // Let the GDML module do magnet-specific setup
    gdmlMagnetModule->assignMediaToGDMLVolumes();
    gdmlMagnetModule->alignMagnetsToLayout(world);
    gdmlMagnetModule->harmoniseVolumeNames();

    // Build geometry for non-magnet modules
    for (auto m : mModulesVec) {
      m->createGeometry(world);
    }
    world->SetName("World");    
    auto* topNode = geom->GetTopNode();
    topNode->SetName(
      fmt::format("{}_{}", world->GetName(), topNode->GetNumber()).c_str());
  }

  // Apply consistent coloring: Iron/Fe → blue, Copper/Cu → yellow
  {
    auto* vols = geom->GetListOfVolumes();
    if (vols) {
      for (int i = 0; i < vols->GetEntriesFast(); ++i) {
        auto* v = static_cast<TGeoVolume*>(vols->At(i));
        if (!v) {
          continue;
        }
        auto* med = v->GetMedium();
        if (!med) {
          continue;
        }
        const TGeoMaterial* mat = med->GetMaterial();
        const char* mname = nullptr;
        if (mat && mat->GetName()) {
          mname = mat->GetName();
        } else if (med->GetName()) {
          mname = med->GetName();
        }
        if (!mname) {
          continue;
        }
        std::string nm = mname;
        if (nm.find("Iron") != std::string::npos || nm.find("Fe") != std::string::npos) {
          v->SetLineColor(kBlue);
        } else if (nm.find("Copper") != std::string::npos || nm.find("Cu") != std::string::npos) {
          v->SetLineColor(kYellow);
        }
      }
    }
  }
  // Finalize
  geom->GetTopVolume()->Voxelize("");
  geom->CloseGeometry();
  for (auto m : mModulesVec) {
    m->setAlignableEntries();
  }
  auto outDir = na6p::utils::Str::rectifyDirectory(na6p::conf::ConfigurableParam::getOutputDir());
  geom->Export((outDir + "geometry.root").c_str());
}

void NA6PDetector::setVerbosity(int v)
{
  mVerbosity = v;
  for (auto m : mModulesVec) {
    m->setVerbosity(v);
  }
}

void NA6PDetector::createHitsOutput(const std::string& dir)
{
  for (auto* det : mActiveModulesVec) {
    det->createHitsOutput(dir);
  }
}

void NA6PDetector::writeHits(const std::vector<int>& remapping)
{
  for (auto* det : mActiveModulesVec) {
    det->writeHits(remapping);
  }
}

void NA6PDetector::closeHitsOutput()
{
  for (auto* det : mActiveModulesVec) {
    det->closeHitsOutput();
  }
}
