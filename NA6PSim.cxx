// NA6PCCopyright

#include "ConfigurableParam.h"
#include <fairlogger/Logger.h>
#include <TApplication.h>
#include <TMCManager.h>
#include <TGeoGlobalMagField.h>
#include <boost/program_options.hpp>
#include <filesystem>
#include <TTree.h>
#include <TFile.h>
#include "MagneticField.h"
#include "StringUtils.h"
#include "NA6PMC.h"
#include "NA6PGenerator.h"
#include "NA6PVerTelHit.h"
#include "NA6PVerTelDigitizer.h"
#include "TG4RunConfiguration.h"
#include "TGeant4.h"

int main(int argc, char** argv)
{
  const std::string layoutIni = "na6pLayout.ini";
  fair::Logger::OnFatal([]() { throw std::runtime_error("Fatal error"); });

  namespace bpo = boost::program_options;
  bpo::variables_map vm;
  bpo::options_description opt_general("Usage:\n  " + std::string(argv[0]));
  bpo::options_description opt_hidden("");
  bpo::options_description opt_all;
  bpo::positional_options_description opt_pos;

  try {
    auto add_option = opt_general.add_options();
    add_option("help,h", "Print this help message");
    add_option("verbosity,v", bpo::value<int>()->default_value(0), "verbosity level [0 = no output]");
    add_option("configKeyValues", bpo::value<std::string>()->default_value(""), "comma-separated configKeyValues");
    add_option("load-ini", bpo::value<std::string>()->default_value(""), "load configurables from ini file (if defined), overridden by configKeyValues");
    add_option("disable-write-ini", bpo::value<bool>()->default_value(false)->implicit_value(true), "do not write ini file");
    add_option("nevents,n", bpo::value<int32_t>()->default_value(1), "number of events to generate");
    add_option("generator,g", bpo::value<std::string>()->default_value(""), "generator defintion root C macro, must return NA6PGenerator pointer");
    add_option("user-hooks,u", bpo::value<std::string>()->default_value(""), "root macro C macro with user hooks for initialization");
    add_option("user-vertex,V", bpo::value<std::string>()->default_value(""), "root macro C macro with user method for vertex generation");
    add_option("rnd-seed,r", bpo::value<int64_t>()->default_value(-1), "random number seed, <0: generate from time + processID");
    add_option("save-bfield,b", bpo::value<bool>()->default_value(false)->implicit_value(true), "Save magnetic field to file");
    add_option("doDigitization,dig", bpo::value<bool>()->default_value(true), "run hits->digits");
    add_option("digitize-only", bpo::value<bool>()->default_value(false)->implicit_value(true), "only run hits->digits using existing output files");
    add_option("event-offset", bpo::value<uint32_t>()->default_value(0), "global event ID offset for this job");
    add_option("skip-events", bpo::value<int64_t>()->default_value(0), "events to skip in the input generator before simulating");
    opt_all.add(opt_general).add(opt_hidden);
    bpo::store(bpo::command_line_parser(argc, argv).options(opt_all).positional(opt_pos).run(), vm);
    if (vm.count("help")) {
      std::cout << opt_general << std::endl;
      exit(0);
    }

    bpo::notify(vm);
  } catch (bpo::error& e) {
    std::cerr << "ERROR: " << e.what() << std::endl
              << std::endl;
    std::cerr << opt_general << std::endl;
    exit(1);
  } catch (std::exception& e) {
    std::cerr << e.what() << ", application will now exit" << std::endl;
    exit(2);
  }
  auto flini = vm["load-ini"].as<std::string>();
  try {
    na6p::conf::ConfigurableParam::updateFromString(vm["configKeyValues"].as<std::string>());
    if (!flini.empty()) {
      na6p::conf::ConfigurableParam::updateFromFile(flini, "", true);
    }
  } catch (const std::exception& error) {
    std::cerr << "ERROR: Failed to load configuration: " << error.what() << std::endl;
    return 2;
  } catch (...) {
    std::cerr << "ERROR: Failed to load configuration: unknown error" << std::endl;
    return 2;
  }
  LOGP(info, "Printing all configs");
  na6p::conf::ConfigurableParam::printAllKeyValuePairs();

  auto runDigitization = [&]() {
    LOGP(info, "Run hits -> digits");
    auto dir = na6p::utils::Str::rectifyDirectory(na6p::conf::ConfigurableParam::getOutputDir());
    NA6PVerTelDigitizer dig;
    dig.init((dir + "geometry.root").c_str());
    TFile* fhVT = TFile::Open((dir + "HitsVerTel.root").c_str());
    if (!fhVT || fhVT->IsZombie()) {
      LOGP(fatal, "Failed to open file {}", dir + "HitsVerTel.root");
    }
    TTree* thVT = (TTree*)fhVT->Get("hitsVerTel");
    if (!thVT) {
      LOGP(fatal, "Failed to get tree hitsVerTel from file {}", dir + "HitsVerTel.root");
    }
    std::vector<NA6PVerTelHit> vtHits, *vtHitsPtr = &vtHits;
    thVT->SetBranchAddress("VerTel", &vtHitsPtr);
    int nEvVT = thVT->GetEntriesFast();
    for (int jEv = 0; jEv < nEvVT; jEv++) {
      thVT->GetEvent(jEv);
      int nHits = vtHits.size();
      LOGP(info, "Digitize VerTel Event {} nHits = {}", jEv, nHits);
      dig.setEventMetaData(vm["event-offset"].as<uint32_t>() + jEv);
      dig.process(vtHits);
    }
    dig.closeDigitsOutput();
    delete thVT;
    fhVT->Close();
    delete fhVT;
  };

  if (vm["digitize-only"].as<bool>()) {
    runDigitization();
    return 0;
  }

  auto mc = new NA6PMC("NA6PMCApp", "NA6P Virtual Monte Carlo Application");
  mc->setVerbosity(vm["verbosity"].as<int>());
  mc->setEventOffset(vm["event-offset"].as<uint32_t>());
  mc->setRandomSeed(vm["rnd-seed"].as<int64_t>());
  if (!vm["user-hooks"].as<std::string>().empty()) {
    mc->setupUserHooks(vm["user-hooks"].as<std::string>());
  }
  if (!vm["user-vertex"].as<std::string>().empty()) {
    mc->setupUserVertex(vm["user-vertex"].as<std::string>());
  }
  // mag field definition
  auto magField = new MagneticField();
  magField->loadField();
  magField->setAsGlobalField();
  if (vm["save-bfield"].as<bool>()) {
    LOGP(info, "Dumping the magnetic field map");
    magField->dumpMagFieldFromConfig();
  }

  auto runConfig = new TG4RunConfiguration("geomRoot", "FTFP_BERT");
  auto geant4 = new TGeant4("TGeant4", "Geant4 Monte Carlo Engine", runConfig, argc, argv);

  TVirtualMC::GetMC()->SetMagField(TGeoGlobalMagField::Instance()->GetField());
  //
  if (!vm["generator"].as<std::string>().empty()) {
    mc->setupGenerator(vm["generator"].as<std::string>());
  } else {
    LOGP(info, "generator is not set");
  }
  mc->init();

  int nEvents = vm["nevents"].as<int32_t>(); // Number of events to simulate
  auto maxEvFromGenerator = mc->canGenerateMaxEvents();
  if (maxEvFromGenerator > 0) {
    LOGP(info, "Generator can generate at most {} events", maxEvFromGenerator);
    if (nEvents < 0 || nEvents > maxEvFromGenerator) {
      nEvents = maxEvFromGenerator;
    }
  }
  if (nEvents && mc->getGenerator()) {
    LOGP(info, "Processing {} events", nEvents);
    auto skipEvents = vm["skip-events"].as<int64_t>();
    if (skipEvents > 0) {
      LOGP(info, "Skipping {} generator events", skipEvents);
      mc->getGenerator()->skipEvents(skipEvents);
    }
    // RS: if we want to silence long geant initialization output, we can add here MiscUtils::silenceStdOut("for Geant4 inits");
    // and add in the NA6PMC::BeginEvent() MiscUtils::reviveStdOut();
    TVirtualMC::GetMC()->ProcessRun(nEvents);
  }
  if (!vm["disable-write-ini"].as<bool>()) {
    try {
      std::string pth = na6p::utils::Str::rectifyDirectory(na6p::conf::ConfigurableParam::getOutputDir());
      if (pth.empty()) {
        pth = std::filesystem::current_path();
      }
      na6p::conf::ConfigurableParam::writeINI(layoutIni);
      LOGP(info, "Stored configurable params to {}/{}", pth, layoutIni);
    } catch (std::exception e) {
      LOGP(error, "Failed to store configurable params, the reason is {}", e.what());
    }
  }
  //  ~/aliroot/sw/SOURCES/GEANT4_VMC/v5-3/v5-3/examples/Gflash
  delete geant4;
  delete mc;

  const bool doHitsToDigits = vm["doDigitization"].as<bool>();
  if (doHitsToDigits) {
    runDigitization();
  }

  return 0;
}
