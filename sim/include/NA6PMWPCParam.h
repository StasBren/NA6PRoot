// NA6PCCopyright

#ifndef NA6P_MWPC_PARAMS_H_
#define NA6P_MWPC_PARAMS_H_

#include "ConfigurableParam.h"
#include "ConfigurableParamHelper.h"

struct NA6PMWPCParam : public na6p::conf::ConfigurableParamHelper<NA6PMWPCParam> {
  static constexpr int MaxStations = 7;

  // User-facing acceptance/layout parameters only. Detailed chamber-construction
  // constants live in NA6PMWPCChamber so an old .ini file cannot silently
  // override an updated expert detector definition.

  // Maximum standard chamber body size, cm, in chamber-local coordinates.
  float bodyX = 53.2f;
  float bodyY = 68.56f;

  // Approximate required detector working rectangles (global detector X,Y), cm.
  float stationWorkingAreaX[MaxStations] = {220.f, 230.f, 310.f, 320.f, 440.f, 500.f, 0.f};
  float stationWorkingAreaY[MaxStations] = {220.f, 240.f, 310.f, 320.f, 410.f, 440.f, 0.f};

  // Regular chamber grids in detector coordinates. NX is horizontal (X), NY
  // vertical (Y). These are intended to be varied for acceptance studies.
  int stationGridNX[MaxStations] = {4, 4, 6, 6, 8, 9, 0};
  int stationGridNY[MaxStations] = {8, 5, 7, 7, 9, 10, 0};

  // Active-gas overlap between neighbouring modules.
  float activeOverlapX = 3.0f;
  float activeOverlapY = 3.0f;

  // Four-level checkerboard longitudinal offset:
  // zGas-zStation = (q-1.5)*staggerZStep, q=A,B,C,D -> 0,1,2,3.
  float staggerZStep = 4.0f;

  // Optional reduced vertical active size for MS0.
  bool useNarrowMS0 = true;
  float ms0GasY = 30.72f;

  NA6PParamDef(NA6PMWPCParam, "mwpc");
};

#endif
