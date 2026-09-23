// NA6PCCopyright

#ifdef __CLING__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class NA6PMCStack + ;
#pragma link C++ class NA6PGenerator + ;
#pragma link C++ class NA6PGenCutParam + ;
#pragma link C++ class NA6PMWPCParam + ;

#pragma link C++ class na6p::conf::ConfigurableParamHelper < NA6PGenCutParam> + ;
#pragma link C++ class na6p::conf::ConfigurableParamHelper < NA6PMWPCParam> + ;

#pragma link C++ class NA6PBaseHit + ;
#pragma link C++ class std::vector < NA6PBaseHit> + ;

#pragma link C++ class NA6PVerTelHit + ;
#pragma link C++ class std::vector < NA6PVerTelHit> + ;

#pragma link C++ class NA6PMuonSpecHit + ;
#pragma link C++ class std::vector < NA6PMuonSpecHit> + ;

#pragma link C++ struct VTPixID+;
#pragma link C++ class NA6PVerTelDigit + ;
#pragma link C++ class std::vector < NA6PVerTelDigit> + ;

#pragma link C++ class std::unordered_map < std::string, float> + ;

#pragma link C++ class NA6PMCEventHeader + ;
#pragma link C++ class std::vector < NA6PMCEventHeader> + ;

#pragma link C++ class NA6PMCGenHeader + ;
#pragma link C++ class std::vector < NA6PMCGenHeader> + ;

#pragma link C++ class NA6PGenBox + ;
#pragma link C++ class NA6PGenParam + ;
#pragma link C++ class NA6PGenHepMC + ;
#pragma link C++ class NA6PGenHisto + ;
#pragma link C++ class NA6PGenCocktail + ;

#pragma link C++ class NA6PVerTelSegmentation + ;
#pragma link C++ class NA6PVerTelDigitizer + ;

#pragma link C++ class std::vector < TParticle> + ;

#pragma link C++ class std::vector < TLorentzVector> + ;
#pragma link C++ class std::vector < ROOT::Math::XYZTVector> + ;
#pragma link C++ class std::vector < ROOT::Math::XYZTVectorF> + ;
#pragma link C++ class std::vector < ROOT::Math::PtEtaPhiEVector> + ;
#pragma link C++ class std::vector < ROOT::Math::PtEtaPhiMVector > + ;
#pragma link C++ class std::vector < ROOT::Math::PxPyPzMVector > + ;

#pragma link C++ class GenMUONLMR + ;

#pragma link C++ class UserHook + ;

#endif
