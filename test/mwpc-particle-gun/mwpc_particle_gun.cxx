// NA6PCCopyright

#include "MWPCTestModule.h"

#include "NA6PMC.h"
#include "NA6PDetector.h"
#include "NA6PModule.h"
#include "NA6PTGeoHelper.h"

#include <TG4RunConfiguration.h>
#include <TGeant4.h>
#include <TGeoManager.h>
#include <TGeoOverlap.h>
#include <TGeoVolume.h>
#include <TVirtualMC.h>

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>

#ifndef MWPC_TEST_SOURCE_DIR
#define MWPC_TEST_SOURCE_DIR "."
#endif

namespace
{
int findFreeActiveID()
{
  auto& detector = NA6PDetector::instance();
  for (int activeID = 0; activeID < NA6PModule::MaxActiveID; ++activeID) {
    bool used = false;
    for (int i = 0; i < detector.getNActiveModules(); ++i) {
      used = used || detector.getActiveModule(i)->getActiveID() == activeID;
    }
    if (!used) {
      return activeID;
    }
  }
  throw std::runtime_error("No unused NA6PRoot active-module ID is available for the MWPC test");
}

class MWPCParticleGunApplication final : public NA6PMC
{
 public:
  explicit MWPCParticleGunApplication(MWPCTestModule& module)
    : NA6PMC("NA6PMWPCParticleGun", "Isolated detailed MWPC particle-gun test"),
      mModule(module),
      mEvents("events.csv"),
      mGasSteps("gas_steps.csv")
  {
    mEvents << "event,gas_steps,gas_entries,gas_exits,primary_steps,primary_entries,energy_deposit_keV\n";
    mGasSteps
      << "event,track_id,x_cm,y_cm,z_cm,u_mm,v_mm,w_mm,"
      << "px_GeV,py_GeV,pz_GeV,edep_GeV,edep_eV\n";
  }

  void ConstructGeometry() override
  {
    auto* geometry = gGeoManager;
    if (!geometry) {
      geometry = new TGeoManager("NA6PMWPCParticleGun", "One detailed MWPC chamber");
    }
    if (geometry->GetTopVolume()) {
      throw std::runtime_error("Particle-gun test geometry was already constructed");
    }

    auto& detector = NA6PDetector::instance();
    detector.createCommonMaterials();
    auto* world = geometry->MakeBox("World", NA6PTGeoHelper::instance().getMedium("Air"), 100., 100., 100.);
    geometry->SetTopVolume(world);
    world->SetVisibility(false);

    mModule.createMaterials();
    mModule.createGeometry(world);

    geometry->CloseGeometry();
    geometry->CheckOverlaps(1.e-5);

    const int overlapCount = geometry->GetListOfOverlaps() ? geometry->GetListOfOverlaps()->GetEntriesFast() : 0;
    std::ofstream overlapReport("geometry_checks.txt");
    overlapReport << "ROOT overlap reports=" << overlapCount << '\n';
    if (geometry->GetListOfOverlaps()) {
      for (auto* object : *geometry->GetListOfOverlaps()) {
        auto* overlap = dynamic_cast<TGeoOverlap*>(object);
        if (!overlap) {
          continue;
        }
        overlapReport << overlap->GetFirstVolume()->GetName() << " / "
                      << (overlap->GetSecondVolume() ? overlap->GetSecondVolume()->GetName() : "extrusion")
                      << " : " << overlap->GetOverlap() << " cm\n";
      }
    }
    if (overlapCount != 0) {
      throw std::runtime_error("ROOT found overlaps in the isolated MWPC geometry; see geometry_checks.txt");
    }

    geometry->Export("geometry.root");
    TVirtualMC::GetMC()->SetRootGeometry();
  }

  void FinishEvent() override
  {
    const auto& s = mModule.stats();
    mEvents << mEventCount << ','
            << s.gasSteps << ','
            << s.gasEntries << ','
            << s.gasExits << ','
            << s.primarySteps << ','
            << s.primaryEntries << ','
            << std::setprecision(10) << s.energyDepositGeV * 1.e6 << '\n';
    mEvents.flush();

    for (const auto& step : mModule.gasSteps()) {
      mGasSteps
          << mEventCount << ","
          << step.trackID << ","
          << step.xCm << "," << step.yCm << "," << step.zCm << ","
          << step.uMm << "," << step.vMm << "," << step.wMm << ","
          << step.pxGeV << "," << step.pyGeV << "," << step.pzGeV << ","
          << step.eDepGeV << "," << step.eDepGeV * 1.e9 << "\n";
    }
    mGasSteps.flush();

    ++mEventCount;
    if (s.primarySteps > 0) {
      ++mEventsWithPrimaryGasSteps;
    }
    if (s.primaryEntries > 0) {
      ++mEventsWithPrimaryEntry;
    }
    mTotalGasSteps += s.gasSteps;
    mTotalEnergyDepositGeV += s.energyDepositGeV;

    NA6PMC::FinishEvent();
  }

  int eventCount() const { return mEventCount; }
  int eventsWithPrimaryGasSteps() const { return mEventsWithPrimaryGasSteps; }
  int eventsWithPrimaryEntry() const { return mEventsWithPrimaryEntry; }
  long long totalGasSteps() const { return mTotalGasSteps; }
  double totalEnergyDepositGeV() const { return mTotalEnergyDepositGeV; }

 private:
  MWPCTestModule& mModule;
  std::ofstream mEvents;
  std::ofstream mGasSteps;
  int mEventCount = 0;
  int mEventsWithPrimaryGasSteps = 0;
  int mEventsWithPrimaryEntry = 0;
  long long mTotalGasSteps = 0;
  double mTotalEnergyDepositGeV = 0.;
};
} // namespace

int main(int argc, char** argv)
{
  if (argc < 4 || argc > 6) {
    std::cerr << "Usage: mwpc_particle_gun OUTPUT_DIR EVENTS MOMENTUM_GEV [PDG=13] [SEED=20260914]\n";
    return 2;
  }

  try {
    const std::filesystem::path output = std::filesystem::absolute(argv[1]);
    const int events = std::stoi(argv[2]);
    const double momentumGeV = std::stod(argv[3]);
    const int pdg = argc >= 5 ? std::stoi(argv[4]) : 13;
    const int seed = argc >= 6 ? std::stoi(argv[5]) : 20260914;
    const bool gridMode = std::getenv("MWPC_GUN_GRID_N") != nullptr;

    if (events < 1 || events > 100000) {
      throw std::runtime_error("EVENTS must be in the range 1..100000");
    }
    if (momentumGeV <= 0.) {
      throw std::runtime_error("MOMENTUM_GEV must be positive");
    }
    if (pdg != 13 && pdg != -13) {
      throw std::runtime_error("This smoke test accepts PDG 13 (mu-) or -13 (mu+) only");
    }
    if (seed <= 0) {
      throw std::runtime_error("SEED must be positive");
    }

    std::filesystem::create_directories(output);
    if (!std::filesystem::is_empty(output)) {
      throw std::runtime_error("Use a fresh empty output directory");
    }
    std::filesystem::current_path(output);

    TGeoManager::SetDefaultUnits(TGeoManager::kRootUnits);

    auto& detector = NA6PDetector::instance();
    auto module = std::make_unique<MWPCTestModule>(findFreeActiveID());
    detector.addModule(module.get());

    auto* application = new MWPCParticleGunApplication(*module);
    application->setVerbosity(1);
    application->setRandomSeed(seed);

    auto* runConfiguration = new TG4RunConfiguration("geomRoot", "FTFP_BERT", "stepLimiter");
    runConfiguration->SetMTApplication(false);
    auto* geant4 = new TGeant4("TGeant4", "NA6PRoot isolated MWPC particle gun", runConfiguration);

    const std::filesystem::path gunPath = std::filesystem::path(MWPC_TEST_SOURCE_DIR) / "mwpc_gun.C";
    const std::string generator = gunPath.string() + "+(" + std::to_string(momentumGeV) + "," + std::to_string(pdg) + ")";
    if (!application->setupGenerator(generator)) {
      throw std::runtime_error("Failed to configure the MWPC particle gun");
    }

    application->init();
    geant4->ProcessGeantCommand(("/random/setSeeds " + std::to_string(seed) + " 31").c_str());
    geant4->ProcessRun(events);

    // The central smoke test expects every primary to enter the gas. In grid mode
    // many tracks are intentionally outside the active rectangle; the dedicated
    // area-scan analysis checks the expected/observed gas-ID map point by point.
    const bool passed = application->eventCount() == events &&
                        (gridMode ||
                         (application->eventsWithPrimaryGasSteps() == events &&
                          application->eventsWithPrimaryEntry() == events));

    std::ofstream summary("summary.txt");
    summary << "mode=" << (gridMode ? "grid" : "central") << '\n'
            << "events=" << application->eventCount() << '\n'
            << "events_with_primary_gas_steps=" << application->eventsWithPrimaryGasSteps() << '\n'
            << "events_with_primary_entry=" << application->eventsWithPrimaryEntry() << '\n'
            << "total_gas_steps=" << application->totalGasSteps() << '\n'
            << "total_energy_deposit_keV=" << std::setprecision(10)
            << application->totalEnergyDepositGeV() * 1.e6 << '\n'
            << "momentum_GeV=" << momentumGeV << '\n'
            << "pdg=" << pdg << '\n'
            << "seed=" << seed << '\n'
            << "result=" << (passed ? "PASS" : "FAIL") << '\n';
    summary.close();

    std::cout << "\nMWPC particle-gun transport result: " << (passed ? "PASS" : "FAIL") << '\n'
              << "  primary entered gas: " << application->eventsWithPrimaryEntry() << '/' << events << " events\n"
              << "  primary had gas steps: " << application->eventsWithPrimaryGasSteps() << '/' << events << " events\n"
              << "  total gas steps: " << application->totalGasSteps() << '\n'
              << "  total gas energy deposit: " << application->totalEnergyDepositGeV() * 1.e6 << " keV\n"
              << "  outputs: " << output << '\n';

    delete geant4;
    delete application;
    return passed ? 0 : 1;
  } catch (const std::exception& error) {
    std::cerr << "MWPC particle-gun test failed: " << error.what() << '\n';
    return 1;
  }
}
