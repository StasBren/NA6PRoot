// NA6PCCopyright

#include "NA6PMWPCChamber.h"

#include "NA6PMWPCParam.h"
#include "NA6PModule.h"
#include "NA6PTGeoHelper.h"

#include <TColor.h>
#include <TGeoManager.h>
#include <TGeoMaterial.h>
#include <TGeoMatrix.h>
#include <TGeoVolume.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>
#include <utility>

namespace
{
constexpr double kMolarMassH = 1.00794;
constexpr double kMolarMassC = 12.0107;
constexpr double kMolarMassN = 14.0067;
constexpr double kMolarMassO = 15.9994;
constexpr double kMolarMassSi = 28.0855;
constexpr double kMolarMassAr = 39.948;
constexpr double kMolarMassCu = 63.546;

constexpr double kHoneycombEdgeWall = 0.1;
constexpr double kReadoutEndMargin = 0.78;
constexpr double kElectronicsPlaneOverhang = 6.8;
constexpr double kDividerPlaneLengthMargin = 0.2;
constexpr double kDividerPlaneOverhang = 4.0;
constexpr double kDividerPlaneThickness = 0.2;
constexpr double kOuterSkinFront = 0.015;
constexpr double kHoneycombFront = 0.6;
constexpr double kReadoutFR4 = 0.08;
constexpr double kReadoutCopper = 0.0017;
constexpr double kGasGap = 0.6;
constexpr double kCoverFR4 = 0.015;
constexpr double kCoverCopper = 0.0017;
constexpr double kHoneycombBack = 0.6;
constexpr double kOuterSkinBack = 0.015;
constexpr bool kCoverCopperOnGasFace = false;
constexpr bool kReadoutCopperOnElectronicsExtension = true;
constexpr bool kIncludeDividerPlane = true;
constexpr double kArgonMoleFraction = 0.7;
constexpr double kGasTemperatureK = 293.15;
constexpr double kGasPressurePa = 101325.0;
}

NA6PMWPCChamber::NA6PMWPCChamber(const NA6PModule& module, Materials materials,
                                 double bodyXOverride, double bodyYOverride)
  : mModule(module), mMaterials(std::move(materials)),
    mBodyXOverride(bodyXOverride), mBodyYOverride(bodyYOverride)
{
}

void NA6PMWPCChamber::requirePositive(double value, const char* name, bool zeroAllowed)
{
  if (!std::isfinite(value) || value < 0. || (!zeroAllowed && value == 0.)) {
    throw std::runtime_error(std::string("Invalid MWPC parameter: ") + name);
  }
}

bool NA6PMWPCChamber::overlaps(const Part& a, const Part& b)
{
  constexpr double tolerance = 1.e-8;
  for (int axis = 0; axis < 3; ++axis) {
    const double ahi = a.centre[axis] + a.fullSize[axis] / 2.;
    const double alo = a.centre[axis] - a.fullSize[axis] / 2.;
    const double bhi = b.centre[axis] + b.fullSize[axis] / 2.;
    const double blo = b.centre[axis] - b.fullSize[axis] / 2.;
    if (std::min(ahi, bhi) - std::max(alo, blo) <= tolerance) {
      return false;
    }
  }
  return true;
}

void NA6PMWPCChamber::createMaterials() const
{
  auto& helper = NA6PTGeoHelper::instance();
  auto& matPool = helper.getMatPool();

  // Materials are kept in a shared pool. If a material is already there we
  // simply reuse it; otherwise we define it here.
  if (!matPool.count(mMaterials.copper)) {
    matPool[mMaterials.copper] = new TGeoMaterial(mMaterials.copper.c_str(), kMolarMassCu, 29., 8.96);
  }

  if (!matPool.count(mMaterials.fr4)) {
    auto* fr4 = new TGeoMixture(mMaterials.fr4.c_str(), 4, 1.86);
    fr4->DefineElement(0, kMolarMassSi, 14., 0.2804);
    fr4->DefineElement(1, kMolarMassO, 8., 0.3874);
    fr4->DefineElement(2, kMolarMassC, 6., 0.3052);
    fr4->DefineElement(3, kMolarMassH, 1., 0.0270);
    matPool[mMaterials.fr4] = fr4;
  }

  if (!matPool.count(mMaterials.honeycomb)) {
    const double molarMass = 14. * kMolarMassC + 10. * kMolarMassH +
                             2. * kMolarMassN + 2. * kMolarMassO;
    auto* hc = new TGeoMixture(mMaterials.honeycomb.c_str(), 4, 0.048);
    hc->DefineElement(0, kMolarMassC, 6., 14. * kMolarMassC / molarMass);
    hc->DefineElement(1, kMolarMassH, 1., 10. * kMolarMassH / molarMass);
    hc->DefineElement(2, kMolarMassN, 7., 2. * kMolarMassN / molarMass);
    hc->DefineElement(3, kMolarMassO, 8., 2. * kMolarMassO / molarMass);
    matPool[mMaterials.honeycomb] = hc;
  }

  if (!matPool.count(mMaterials.gas)) {
    requirePositive(kArgonMoleFraction, "argonMoleFraction");
    requirePositive(kGasTemperatureK, "gasTemperatureK");
    requirePositive(kGasPressurePa, "gasPressurePa");
    if (kArgonMoleFraction >= 1.f) {
      throw std::runtime_error("MWPC argonMoleFraction must be between 0 and 1");
    }

    // The chamber gas is a binary Ar/CO2 mixture. We configure the Ar mole
    // fraction and obtain the CO2 fraction from the remainder.
    const double xAr = kArgonMoleFraction;
    const double xCO2 = 1. - xAr;

    const double mCO2 = kMolarMassC + 2. * kMolarMassO;
    const double meanMolarMass = xAr * kMolarMassAr + xCO2 * mCO2;

    // Ideal-gas density. meanMolarMass is in g/mol, so the final factor
    // converts from g/m^3 to g/cm^3, which is what TGeo expects.
    const double density = kGasPressurePa * meanMolarMass /
                           (8.31446261815324 * kGasTemperatureK) * 1.e-6;

    // TGeoMixture is defined through elemental mass fractions rather than
    // molecular fractions. Convert 70/30 Ar/CO2 (or any configured ratio)
    // into the corresponding Ar, C and O mass fractions.
    const double wAr = xAr * kMolarMassAr / meanMolarMass;
    const double wC = xCO2 * kMolarMassC / meanMolarMass;
    const double wO = xCO2 * 2. * kMolarMassO / meanMolarMass;

    auto* gas = new TGeoMixture(mMaterials.gas.c_str(), 3, density);
    gas->DefineElement(0, kMolarMassAr, 18., wAr);
    gas->DefineElement(1, kMolarMassC, 6., wC);
    gas->DefineElement(2, kMolarMassO, 8., wO);
    gas->SetState(TGeoMaterial::kMatStateGas);
    gas->SetTemperature(kGasTemperatureK);
    matPool[mMaterials.gas] = gas;
  }

  // Geant/VMC uses TGeoMedium objects for volumes. The material definitions
  // above live in matPool; addMedium() creates the corresponding media only
  // when they do not already exist.
  for (const auto& name : {mMaterials.fr4, mMaterials.copper, mMaterials.honeycomb, mMaterials.gas}) {
    if (!helper.getMedPool().count(name)) {
      helper.addMedium(name);
    }
  }
}

std::vector<NA6PMWPCChamber::Part> NA6PMWPCChamber::buildParts() const
{
  const auto& p = NA6PMWPCParam::Instance();

  const double w = mBodyXOverride > 0. ? mBodyXOverride : p.bodyX;
  const double h = mBodyYOverride > 0. ? mBodyYOverride : p.bodyY;
  requirePositive(w, "bodyX");
  requirePositive(h, "bodyY");
  requirePositive(NA6PMWPCChamber::InnerFrameWidth, "innerFrameWidth");
  requirePositive(kHoneycombEdgeWall, "honeycombEdgeWall");
  requirePositive(kReadoutEndMargin, "readoutEndMargin", true);
  requirePositive(kElectronicsPlaneOverhang, "electronicsPlaneOverhang", true);
  requirePositive(kDividerPlaneLengthMargin, "dividerPlaneLengthMargin", true);
  requirePositive(kDividerPlaneOverhang, "dividerPlaneOverhang", true);
  requirePositive(kDividerPlaneThickness, "dividerPlaneThickness");

  requirePositive(kOuterSkinFront, "outerSkinFront");
  requirePositive(kHoneycombFront, "honeycombFront");
  requirePositive(kReadoutFR4, "readoutFR4");
  requirePositive(kReadoutCopper, "readoutCopper", true);
  requirePositive(kGasGap, "gasGap");
  requirePositive(kCoverFR4, "coverFR4");
  requirePositive(kCoverCopper, "coverCopper", true);
  requirePositive(kHoneycombBack, "honeycombBack");
  requirePositive(kOuterSkinBack, "outerSkinBack");

  const double frame = NA6PMWPCChamber::InnerFrameWidth;
  const double wall = kHoneycombEdgeWall;
  const double readoutW = w + kElectronicsPlaneOverhang;
  const double readoutH = h - 2. * kReadoutEndMargin;
  const double gasW = w - 2. * frame;
  const double gasH = h - 2. * frame;

  if (std::min({gasW, gasH, w - 2. * wall, h - 2. * wall, readoutH,
                w - kDividerPlaneLengthMargin}) <= 0.) {
    throw std::runtime_error("MWPC dimensions leave a non-positive derived box");
  }
  if (readoutH < gasH) {
    throw std::runtime_error("MWPC readout board does not cover the gas opening");
  }

  std::vector<std::pair<std::string, double>> stack = {
    {"outerSkinFront", kOuterSkinFront},
    {"honeycombFront", kHoneycombFront},
    {"readoutFR4", kReadoutFR4},
    {"readoutCopper", kReadoutCopper},
    {"gasGap", kGasGap}};
  if (kCoverCopperOnGasFace) {
    stack.emplace_back("coverCopper", kCoverCopper);
    stack.emplace_back("coverFR4", kCoverFR4);
  } else {
    stack.emplace_back("coverFR4", kCoverFR4);
    stack.emplace_back("coverCopper", kCoverCopper);
  }
  stack.emplace_back("honeycombBack", kHoneycombBack);
  stack.emplace_back("outerSkinBack", kOuterSkinBack);

  double depth = 0.;
  for (const auto& layer : stack) {
    depth += layer.second;
  }
  double cursor = -depth / 2.;
  std::map<std::string, double> z;
  for (const auto& layer : stack) {
    z[layer.first] = cursor + layer.second / 2.;
    cursor += layer.second;
  }

  std::vector<Part> parts;
  auto add = [&](std::string name, MaterialKind material,
                 std::array<double, 3> size, std::array<double, 3> centre,
                 bool sensitive = false) {
    if (std::min({size[0], size[1], size[2]}) <= 0.) {
      throw std::runtime_error("Invalid MWPC solid size: " + name);
    }
    parts.push_back({std::move(name), material, size, centre, sensitive});
  };

  add("OuterSkinFront", MaterialKind::FR4,
      {w, h, kOuterSkinFront}, {0., 0., z.at("outerSkinFront")});
  add("OuterSkinBack", MaterialKind::FR4,
      {w, h, kOuterSkinBack}, {0., 0., z.at("outerSkinBack")});

  const std::array<std::pair<const char*, double>, 2> honeycombLayers{{
    {"Front", kHoneycombFront}, {"Back", kHoneycombBack}}};
  for (const auto& layer : honeycombLayers) {
    const std::string side = layer.first;
    const std::string key = side == "Front" ? "honeycombFront" : "honeycombBack";
    const double t = layer.second;
    const double zz = z.at(key);
    add("Honeycomb" + side, MaterialKind::Honeycomb,
        {w - 2. * wall, h - 2. * wall, t}, {0., 0., zz});
    add("HCWallXMinus" + side, MaterialKind::FR4,
        {wall, h, t}, {-(w - wall) / 2., 0., zz});
    add("HCWallXPlus" + side, MaterialKind::FR4,
        {wall, h, t}, {+(w - wall) / 2., 0., zz});
    add("HCWallYMinus" + side, MaterialKind::FR4,
        {w - 2. * wall, wall, t}, {0., -(h - wall) / 2., zz});
    add("HCWallYPlus" + side, MaterialKind::FR4,
        {w - 2. * wall, wall, t}, {0., +(h - wall) / 2., zz});
  }

  add("ElectronicsPlane", MaterialKind::FR4,
      {readoutW, readoutH, kReadoutFR4},
      {kElectronicsPlaneOverhang / 2., 0., z.at("readoutFR4")});
  if (kReadoutCopper > 0.) {
    const double copperW = kReadoutCopperOnElectronicsExtension ? readoutW : w;
    const double copperX = kReadoutCopperOnElectronicsExtension ? kElectronicsPlaneOverhang / 2. : 0.;
    add("ReadoutCopper", MaterialKind::Copper,
        {copperW, readoutH, kReadoutCopper},
        {copperX, 0., z.at("readoutCopper")});
  }

  add("Gas", MaterialKind::Gas,
      {gasW, gasH, kGasGap}, {0., 0., z.at("gasGap")}, true);
  add("InnerFrameXMinus", MaterialKind::FR4,
      {frame, h, kGasGap}, {-(w - frame) / 2., 0., z.at("gasGap")});
  add("InnerFrameXPlus", MaterialKind::FR4,
      {frame, h, kGasGap}, {+(w - frame) / 2., 0., z.at("gasGap")});
  add("InnerFrameYMinus", MaterialKind::FR4,
      {gasW, frame, kGasGap}, {0., -(h - frame) / 2., z.at("gasGap")});
  add("InnerFrameYPlus", MaterialKind::FR4,
      {gasW, frame, kGasGap}, {0., +(h - frame) / 2., z.at("gasGap")});

  add("ChamberCover", MaterialKind::FR4,
      {w, h, kCoverFR4}, {0., 0., z.at("coverFR4")});
  if (kCoverCopper > 0.) {
    add("CoverCopper", MaterialKind::Copper,
        {w, h, kCoverCopper}, {0., 0., z.at("coverCopper")});
  }

  if (kIncludeDividerPlane && kDividerPlaneOverhang > 0.) {
    add("DividerPlaneApproximation", MaterialKind::FR4,
        {w - kDividerPlaneLengthMargin, kDividerPlaneOverhang, kDividerPlaneThickness},
        {0., -(h + kDividerPlaneOverhang) / 2.,
         z.at("gasGap") - kGasGap / 2. + kDividerPlaneThickness / 2.});
  }

  validate(parts);
  return parts;
}

void NA6PMWPCChamber::validate(const std::vector<Part>& parts) const
{
  int sensitiveCount = 0;
  for (size_t i = 0; i < parts.size(); ++i) {
    const auto& a = parts[i];
    if (a.sensitive) {
      ++sensitiveCount;
    }
    for (int axis = 0; axis < 3; ++axis) {
      if (!std::isfinite(a.fullSize[axis]) || a.fullSize[axis] <= 0. ||
          !std::isfinite(a.centre[axis])) {
        throw std::runtime_error("Invalid MWPC box geometry: " + a.name);
      }
    }
    for (size_t j = 0; j < i; ++j) {
      if (overlaps(a, parts[j])) {
        throw std::runtime_error("Overlapping MWPC parts: " + a.name + " / " + parts[j].name);
      }
    }
  }
  if (sensitiveCount != 1) {
    throw std::runtime_error("MWPC chamber must contain exactly one sensitive gas box");
  }
}

TGeoVolumeAssembly* NA6PMWPCChamber::addTo(TGeoVolume* parent, int chamberID, int localCopyID,
                                            const Placement& placement) const
{
  if (!parent || !gGeoManager) {
    throw std::runtime_error("MWPC chamber requires an existing parent volume and TGeoManager");
  }
  if (chamberID < 0 || chamberID >= NA6PModule::MaxVolID - NA6PModule::MaxNonSensID) {
    throw std::runtime_error("MWPC chamberID is outside the sensitive-volume ID range");
  }
  if (localCopyID < 0 || localCopyID >= NA6PModule::MaxNonSensID) {
    throw std::runtime_error("MWPC localCopyID is outside the non-sensitive ID range");
  }

  createMaterials();
  const auto parts = buildParts();
  auto& helper = NA6PTGeoHelper::instance();

  const std::string assemblyName = "MWPCChamber_" + std::to_string(chamberID);
  auto* assembly = new TGeoVolumeAssembly(assemblyName.c_str());

  int passiveID = 100;
  for (const auto& part : parts) {
    const std::string volumeName = assemblyName + "_" + part.name;
    std::string mediumName;
    int color = kGreen + 2;
    switch (part.material) {
      case MaterialKind::FR4:
        mediumName = mMaterials.fr4;
        break;
      case MaterialKind::Copper:
        mediumName = mMaterials.copper;
        color = kOrange + 7;
        break;
      case MaterialKind::Honeycomb:
        mediumName = mMaterials.honeycomb;
        color = kOrange - 3;
        break;
      case MaterialKind::Gas:
        mediumName = mMaterials.gas;
        color = kCyan + 1;
        break;
    }

    auto* volume = gGeoManager->MakeBox(volumeName.c_str(), helper.getMedium(mediumName),
                                        part.fullSize[0] / 2., part.fullSize[1] / 2., part.fullSize[2] / 2.);
    volume->SetLineColor(color);
    volume->SetTransparency(part.sensitive ? 75 : 15);

    const int copyID = part.sensitive ? mModule.composeSensorVolID(chamberID)
                                      : mModule.composeNonSensorVolID(passiveID++);
    assembly->AddNode(volume, copyID,
                      new TGeoTranslation(part.centre[0], part.centre[1], part.centre[2]));
  }

  parent->AddNode(assembly, mModule.composeNonSensorVolID(localCopyID),
                  new TGeoTranslation(placement.x, placement.y, placement.z));
  return assembly;
}

std::array<double, 3> NA6PMWPCChamber::envelopeFullSize() const
{
  const auto parts = buildParts();
  std::array<double, 3> low{
    std::numeric_limits<double>::infinity(),
    std::numeric_limits<double>::infinity(),
    std::numeric_limits<double>::infinity()};
  std::array<double, 3> high{
    -std::numeric_limits<double>::infinity(),
    -std::numeric_limits<double>::infinity(),
    -std::numeric_limits<double>::infinity()};

  for (const auto& part : parts) {
    for (int axis = 0; axis < 3; ++axis) {
      low[axis] = std::min(low[axis], part.centre[axis] - part.fullSize[axis] / 2.);
      high[axis] = std::max(high[axis], part.centre[axis] + part.fullSize[axis] / 2.);
    }
  }
  return {high[0] - low[0], high[1] - low[1], high[2] - low[2]};
}

std::array<double, 3> NA6PMWPCChamber::gasFullSize() const
{
  for (const auto& part : buildParts()) {
    if (part.sensitive) {
      return part.fullSize;
    }
  }
  throw std::runtime_error("MWPC gas volume is missing");
}

std::array<double, 3> NA6PMWPCChamber::gasCentre() const
{
  for (const auto& part : buildParts()) {
    if (part.sensitive) {
      return part.centre;
    }
  }
  throw std::runtime_error("MWPC gas volume is missing");
}
