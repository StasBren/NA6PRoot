// ROOT-only 2100-V working point diagnostics, with independent-muon bootstrap.
// Reuses validated CSV parser, ROOT formatting, statistical definitions from HV scan.
// Do NOT interpret 1 or 5 muons as a precision result.
#include "../hv-high-2100-2250/AnalyzeHVHigh_ROOT.C"

void AnalyzeHV2100_ROOT(const char* outputDir) {
  using namespace HvAna;
  try {
    Style();
    const std::string dir(outputDir);
    gSystem->mkdir((dir+"/analysis").c_str(),kTRUE);
    VSample a=Load(dir,2100);
    Bootstrap(a,500);
    std::ofstream f(dir+"/analysis/2100V_gain_statistics.csv");
    if(!f)throw std::runtime_error("Cannot write 2100V statistics file");
    f << "voltage_V,n_muons,n_primary,n_positive,zero_fraction,"
         "mean_gain_all,mean_gain_all_bootstrap_error,"
         "rms_over_mean_all,bootstrap_error_all,"
         "mean_gain_positive,rms_over_mean_positive,bootstrap_error_positive,"
         "mean_muon_effective_gain,muon_mean_bootstrap_error,"
         "rms_over_mean_muon,bootstrap_error_muon\n";
    f << std::setprecision(12) << 2100 << "," << a.muons.size()
      << "," << a.allGain.size() << "," << a.positiveGain.size()
      << "," << a.zeros/a.allGain.size() << ","
      << a.seedMean << "," << a.eMean << ","
      << a.rAll << "," << a.eAll << ","
      << a.posMean << "," << a.rPos << "," << a.ePos << ","
      << a.gmuMean << "," << a.eGmuMean << ","
      << a.rMuon << "," << a.eMuon << "\n";
    f.close();

    auto gain=MakeHist("positive_gain_2100_",2100,a.positiveGain,true);
    TCanvas cGain("c2100gain","2100 V gain spectrum",900,650);
    gain->SetLineColor(kBlue+1);
    gain->SetTitle(";log_{10}(collected gain per Heed primary);Probability density");
    gain->GetYaxis()->SetTitleOffset(1.25);
    gain->Draw("HIST");
    Label("2100 V, P3 noFW, 30 #mum anode, Penning r=0.57; g>0");
    Save(cGain,dir+"/analysis","positive_gain_2100V");

    auto asym=MakeHist("signed_muon_asym_2100_",2100,a.A,false);
    TCanvas cAsym("c2100asym","2100 V event asymmetry",900,650);
    asym->SetLineColor(kBlack);
    asym->SetTitle(";Signed event ion-birth asymmetry A_{#mu};Probability density");
    asym->GetYaxis()->SetTitleOffset(1.25);
    asym->Draw("HIST");
    Label("2100 V; signed ion-birth asymmetry, not cathode charge");
    Save(cAsym,dir+"/analysis","signed_asymmetry_2100V");

    std::cout << "2100 V: " << a.muons.size() << " muons, "
              << a.allGain.size() << " Heed primaries, mean gain="
              << a.seedMean << ", RMS/mean (all)=" << a.rAll
              << ", RMS/mean (positive)=" << a.rPos
              << ", muon RMS/mean=" << a.rMuon
              << ", zero fraction=" << a.zeros/a.allGain.size() << "\n"
              << "ROOT results: " << dir << "/analysis/\n";
  } catch (const std::exception& e) {
    ::Error("AnalyzeHV2100_ROOT","%s",e.what());
    throw;
  }
}
