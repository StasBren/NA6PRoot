// NA6PCCopyright

#include "MWPCTestModule.h"

#include "NA6PMCStack.h"
#include "NA6PMWPCChamber.h"

#include <TVirtualMC.h>
#include <TGeoManager.h>
#include <TGeoVolume.h>
#include <TLorentzVector.h>

#include <fstream>
#include <stdexcept>

MWPCTestModule::MWPCTestModule(int activeID) : NA6PModule("MWPCTest")
{
  setActiveID(activeID);
}

void MWPCTestModule::createMaterials()
{
  const NA6PMWPCChamber::Materials materials = {
    addName("FR4"),
    addName("Copper"),
    addName("Honeycomb"),
    addName("ArCO2")};
  NA6PMWPCChamber chamber(*this, materials);
  chamber.createMaterials();
}

void MWPCTestModule::createGeometry(TGeoVolume* world)
{
  if (!world) {
    throw std::runtime_error("MWPC particle-gun test requires a world volume");
  }

  const NA6PMWPCChamber::Materials materials = {
    addName("FR4"),
    addName("Copper"),
    addName("Honeycomb"),
    addName("ArCO2")};
  NA6PMWPCChamber chamber(*this, materials);

  // One complete chamber, body centred at the origin. chamberID=0 is the
  // sensitive gas ID; localCopyID=0 is the passive chamber assembly ID.
  chamber.addTo(world, 0, 0, {0., 0., 0.});

  const auto gas = chamber.gasFullSize();
  const auto gasCentre = chamber.gasCentre();
  const auto envelope = chamber.envelopeFullSize();
  std::ofstream out("geometry_summary.txt");
  out << "gas_cm=" << gas[0] << ',' << gas[1] << ',' << gas[2] << '\n'
      << "gas_centre_cm=" << gasCentre[0] << ',' << gasCentre[1] << ',' << gasCentre[2] << '\n'
      << "envelope_cm=" << envelope[0] << ',' << envelope[1] << ',' << envelope[2] << '\n';
}

bool MWPCTestModule::stepManager(int volID)
{
  const int sensorID = NA6PModule::volID2SensID(volID);
  if (sensorID != 0) {
    throw std::runtime_error("MWPC particle-gun test received an unexpected sensitive-volume ID");
  }

  auto* mc = TVirtualMC::GetMC();
  if (!mc) {
    throw std::runtime_error("MWPC particle-gun test has no active VMC engine");
  }

  ++mStats.gasSteps;
  const double eDep = mc->Edep();
  mStats.energyDepositGeV += eDep;

  if (eDep > 0.) {
    TLorentzVector pos;
    TLorentzVector mom;
    mc->TrackPosition(pos);
    mc->TrackMomentum(mom);

    double master[3] = {pos.X(), pos.Y(), pos.Z()};
    double local[3] = {0., 0., 0.};
    if (gGeoManager) {
      gGeoManager->MasterToLocal(master, local);
    }

    auto* stackNow = static_cast<NA6PMCStack*>(mc->GetStack());
    GasStepRecord rec;
    rec.trackID = stackNow ? stackNow->GetCurrentTrackNumber() : -1;
    rec.xCm = pos.X();
    rec.yCm = pos.Y();
    rec.zCm = pos.Z();

    // Isolated-chamber local axes before installation rotation:
    //   local x = across wires (Garfield u)
    //   local z = chamber normal (Garfield v)
    //   local y = along wires (Garfield w)
    rec.uMm = 10. * local[0];
    rec.vMm = 10. * local[2];
    rec.wMm = 10. * local[1];

    rec.pxGeV = mom.Px();
    rec.pyGeV = mom.Py();
    rec.pzGeV = mom.Pz();
    rec.eDepGeV = eDep;
    mGasSteps.push_back(rec);
  }
  if (mc->IsTrackEntering()) {
    ++mStats.gasEntries;
  }
  if (mc->IsTrackExiting() || mc->IsTrackOut() || mc->IsTrackStop()) {
    ++mStats.gasExits;
  }

  auto* stack = static_cast<NA6PMCStack*>(mc->GetStack());
  if (stack) {
    if (stack->GetCurrentTrackNumber() == 0) {
      ++mStats.primarySteps;
      if (mc->IsTrackEntering()) {
        ++mStats.primaryEntries;
      }
    }
    stack->addHit(getActiveIDBit());
  }
  return true;
}
