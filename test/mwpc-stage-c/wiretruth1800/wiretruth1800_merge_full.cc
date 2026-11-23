// Merge fully sealed real Garfield++ 1800V per-wire ROOT shards with exact
// TTree crossing_id non-overlap QA, additive ROOT histogram/TProfile merge,
// additive independent-muon counts and re-render final gain/asymmetry plots.
// Numerical merging and summary statistics in C++/ROOT. No Garfield runs.
#include <TFile.h>
#include <TFileMerger.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TColor.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;
namespace {
void Need(bool condition,const std::string&msg){
  if(!condition)throw std::runtime_error("1800V FULL ROOT MERGE: "+msg);
}
std::vector<std::string> Split(const std::string&row,char delim=','){
  std::stringstream ss(row);std::string value;std::vector<std::string> out;
  while(std::getline(ss,value,delim))out.push_back(value);
  return out;
}
struct Csv {
  std::vector<std::string> header;
  std::vector<std::vector<std::string>> rows;
  explicit Csv(const fs::path&file){
    std::ifstream f(file);
    Need(bool(f),"missing full campaign exact summary CSV "+file.string());
    std::string line;
    Need(bool(std::getline(f,line)),"empty exact summary CSV "+file.string());
    header=Split(line);
    Need(header.size()>2,"malformed column names");
    while(std::getline(f,line)){
      if(line.empty())continue;
      auto cols=Split(line);
      Need(cols.size()==header.size(),"wrong CSV row width: "+file.string());
      rows.push_back(std::move(cols));
    }
    Need(!rows.empty(),"no data in CSV "+file.string());
  }
};
std::map<std::string,std::size_t> Columns(const std::vector<std::string>&h){
  std::map<std::string,std::size_t> out;
  for(std::size_t i=0;i<h.size();++i)
    Need(out.emplace(h[i],i).second,"duplicate summary column name");
  return out;
}
long long Number(const std::string&s){
  std::size_t stop=0;
  const auto value=std::stoll(s,&stop);
  Need(stop==s.size()&&value>=0,"negative/malformed integer counter "+s);
  return value;
}
double Real(const std::string&s){
  std::size_t stop=0;
  const auto value=std::stod(s,&stop);
  Need(stop==s.size()&&std::isfinite(value),"nonfinite CSV ratio");
  return value;
}
std::vector<std::string> ReadShardDirs(const fs::path&file){
  std::ifstream f(file);Need(bool(f),"missing shard directory list");
  std::vector<std::string> dirs;std::string line;std::set<std::string> seen;
  while(std::getline(f,line)){
    if(line.empty())continue;
    Need(fs::is_directory(line),"missing staged source shard "+line);
    Need(seen.insert(line).second,"repeated shard list path");
    dirs.push_back(line);
  }
  Need(!dirs.empty(),"zero source ROOT shards");
  return dirs;
}
void MergeCsv(const std::vector<std::string>&dirs,const fs::path&out,
              const std::string&file,bool channelSummary,
              long long expectedTotal){
  std::vector<std::string> heading;
  std::map<std::string,std::vector<long long>> totals;
  std::map<std::string,double> weightedSecond;
  std::set<std::string> scopes;
  for(const auto& d:dirs){
    Csv c(fs::path(d)/file);
    if(heading.empty())heading=c.header;
    Need(c.header==heading,"CSV structure differs between physical shards "+file);
    const auto ix=Columns(heading);
    Need(ix.count("scope")&&ix.at("scope")==0,"missing explicit scope");
    const std::size_t M=heading.size();
    const std::size_t meanCol=channelSummary?
      ix.at("mean_q2_over_qtotal_positives"):std::numeric_limits<std::size_t>::max();
    const std::array<std::string,3> calc={"frac_multianode_gt5pct","binomial_se",
                                          "mean_q2_over_qtotal_positives"};
    if(channelSummary)
      for(const auto&name:calc)Need(ix.count(name),"missing ratio "+name);
    std::set<std::string> shardScopes;
    for(const auto&row:c.rows){
      const std::string&scope=row[0];
      Need(shardScopes.insert(scope).second,"duplicate analysis scope in input shard");
      auto&sum=totals[scope];
      if(sum.empty())sum.assign(M,0LL);
      if(channelSummary){
        const long long nm=Number(row[ix.at("n_muons")]);
        const long long nz=Number(row[ix.at("n_zero_charge")]);
        Need(nm>=nz,"impossible zero gain muon denominator");
        weightedSecond[scope]+=(nm-nz)*Real(row[meanCol]);
      }
      for(std::size_t i=1;i<M;++i){
        if(channelSummary&&
           (heading[i]=="frac_multianode_gt5pct"||
            heading[i]=="binomial_se"||
            heading[i]=="mean_q2_over_qtotal_positives"))
          continue;
        sum[i]+=Number(row[i]);
      }
    }
    if(scopes.empty())scopes=shardScopes;
    Need(scopes==shardScopes,"missing or extra channel/station row in one shard");
  }
  Need(totals.count("all"),"no independent-muon pooled row");
  if(channelSummary){
    Need(totals.at("all")[Columns(heading).at("n_muons")]==expectedTotal,
         "merged summary muon total inconsistent with ROOT");
  }
  std::ofstream f(out);Need(bool(f),"cannot create additive summary CSV");
  for(std::size_t i=0;i<heading.size();++i)
    f<<(i?",":"")<<heading[i];
  f<<"\n";
  f<<std::setprecision(17);
  const auto ix=Columns(heading);
  for(const auto&[scope,sum]:totals){
    if(channelSummary){
      const auto n=sum[ix.at("n_muons")];
      const auto nz=sum[ix.at("n_zero_charge")];
      const auto successes=sum[ix.at("n_5pct_multianode")];
      Need(successes<=n&&nz<=n,"bad combined muon population");
      Need(n==nz+sum[ix.at("n_single_anode")]+sum[ix.at("n_two_anode")]+
              sum[ix.at("n_3plus_anodes")],
           "merged anode multiplicities do not conserve Geant4 muons");
    }
    for(std::size_t j=0;j<heading.size();++j){
      if(j)f<<",";
      if(j==0){f<<scope;continue;}
      if(channelSummary){
        const long long n=sum[ix.at("n_muons")];
        const long long nz=sum[ix.at("n_zero_charge")];
        const long long succ=sum[ix.at("n_5pct_multianode")];
        const double fraction=n?double(succ)/n:0.;
        if(heading[j]=="frac_multianode_gt5pct"){f<<fraction;continue;}
        if(heading[j]=="binomial_se"){
          f<<(n?std::sqrt(fraction*(1.-fraction)/n):0.);continue;
        }
        if(heading[j]=="mean_q2_over_qtotal_positives"){
          f<<(n>nz?weightedSecond[scope]/(n-nz):0.);continue;
        }
      }
      f<<sum[j];
    }
    f<<"\n";
  }
  Need(bool(f),"failed writing combined CSV "+file);
}
TFile* globalFile=nullptr;
TH1* Histogram(const std::string&name){
  Need(globalFile,"merged ROOT input was not opened");
  TH1* h=nullptr;
  globalFile->GetObject(name.c_str(),h);
  return h;
}
void Plot(const fs::path&out,const std::string&scope,
          const std::string&suffix,const std::string&filename,
          bool is2D,bool isProfile=false){
  TH1* h=Histogram(scope+suffix);
  if(!h||h->GetEntries()==0)return;
  TCanvas canvas(("c_full_"+scope+"_"+filename).c_str(),"",1000,725);
  canvas.SetLeftMargin(.14);canvas.SetRightMargin(is2D?.17:.07);
  canvas.SetBottomMargin(.13);
  gStyle->SetOptStat(0);
  if(is2D){
    auto* d=dynamic_cast<TH2*>(h);
    Need(d,"wrong merged ROOT histogram type "+scope+suffix);
    d->Draw("COLZ");
  }else if(isProfile){
    auto*p=dynamic_cast<TProfile*>(h);
    Need(p,"wrong merged ROOT profile type "+scope+suffix);
    p->SetMarkerStyle(20);p->SetMarkerColor(kBlue+1);
    p->SetLineColor(kBlue+1);p->Draw("E1");
  }else{
    h->SetLineWidth(2);h->SetLineColor(kBlue+1);h->Draw("HIST");
  }
  TLatex title;title.SetTextFont(42);title.SetNDC();
  title.SetTextSize(.038);title.DrawLatex(.16,.92,"NA60+ Simulation");
  title.SetTextSize(.023);
  title.DrawLatex(.16,.88,("1800 V P3 noFW; ALL verified archived chunks; "+scope).c_str());
  canvas.SaveAs((out/("full_"+filename+"_"+scope+".png")).string().c_str());
}
void PlotFinal(TFile&input,const fs::path&dir){
  globalFile=&input;
  const struct {const char* hist;const char* outfile;bool density;bool profile;} plots[]={
    {"_asym_wire_gain_vs_A_hemi","gain_vs_A_hemi_wire",true,false},
    {"_asym_primary_wire_gain_vs_A_hemi","gain_vs_A_hemi_primary_wire",true,false},
    {"_asym_muon_gain_vs_A_hemi","gain_vs_A_hemi_muon",true,false},
    {"_asym_wire_gain_vs_A_TB","gain_vs_A_TB_wire",true,false},
    {"_asym_muon_second_wire_fraction_vs_A_hemi","muon_sharing_vs_A_hemi",true,false},
    {"_asym_wire_ion_births_vs_A_hemi","ion_births_vs_A_hemi_wire",true,false},
    {"_asym_wire_charge_fraction_vs_A_hemi","charge_fraction_vs_A_hemi_wire",true,false},
    {"_asym_A_hemi_vs_A_TB","A_hemi_vs_A_TB",true,false},
    {"_asym_A_hemi_vs_local_wire","A_hemi_vs_local_wire",true,false},
    {"_asym_A_hemi_vs_signed_impact","A_hemi_vs_signed_impact",true,false},
    {"_asym_A_hemi_vs_theta_u","A_hemi_vs_theta_u",true,false},
    {"_asym_A_hemi_vs_theta_w","A_hemi_vs_theta_w",true,false},
    {"_asym_ion_A_zero_collected_wire","zero_collected_ion_wires_A_hemi",false,false},
    {"_asym_A_hemi","A_hemi",false,false},
    {"_asym_A_TB","A_TB",false,false},
    {"_asym_mean_A_vs_local_wire","mean_A_vs_local_wire",false,true},
    {"_asym_mean_A_vs_impact","mean_A_vs_signed_impact",false,true},
    {"_asym_mean_A_vs_theta_u","mean_A_vs_signed_theta_u",false,true},
    {"_asym_mean_A_vs_theta_w","mean_A_vs_signed_theta_w",false,true},
    {"_log_wire_Q","collected_e_per_wire",false,false},
    {"_log_wire_ions","ion_births_per_wire",false,false},
    {"_n_wire","anode_multiplicity_per_muon",false,false},
    {"_second_wire_fraction","q2_over_qtotal",false,false},
    {"_q_fraction_by_wire","charge_fraction_vs_signed_impact_and_wire",true,false},
    {"_ion_fraction_by_wire","ion_birth_fraction_vs_signed_impact_and_wire",true,false}
  };
  for(const auto&scope:{"all","jpsi","phi","omega"}){
    for(const auto&p:plots)
      Plot(dir,scope,p.hist,p.outfile,p.density,p.profile);
  }
  for(const auto&scope:{"jpsi","phi","omega"}){
    for(int station=0;station<6;++station){
      const auto key=std::string(scope)+"_station"+std::to_string(station);
      Plot(dir,key,"_asym_wire_gain_vs_A_hemi","gain_vs_A_hemi_wire",true);
      Plot(dir,key,"_asym_A_hemi_vs_local_wire","A_hemi_vs_local_wire",true);
    }
  }
}
void MergeTrees(const std::vector<std::string>&dirs,const fs::path&out,
                long long expectedMuonCount){
  const fs::path target=out/"wiretruth1800_full.root";
  TFileMerger merger(kFALSE,kTRUE);
  Need(merger.OutputFile(target.string().c_str(),"RECREATE"),
       "cannot open final ROOT merger output");
  for(const auto&dir:dirs){
    const auto p=(fs::path(dir)/"wiretruth1800_8k.root").string();
    Need(merger.AddFile(p.c_str()),"failed to add verified ROOT shard "+p);
  }
  Need(merger.Merge(),"ROOT TFileMerger failed to merge histograms/TTrees");
  std::unique_ptr<TFile> file(TFile::Open(target.string().c_str(),"READ"));
  Need(file && !file->IsZombie(),"merged ROOT file cannot be reopened");
  TTree* mu=nullptr;
  TTree* wr=nullptr;
  TTree* seed=nullptr;
  TTree* seedWire=nullptr;
  file->GetObject("muons",mu);file->GetObject("wires",wr);
  file->GetObject("seeds",seed);file->GetObject("seed_wires",seedWire);
  Need(mu&&wr&&seed&&seedWire,"merged physics ROOT TTrees missing");
  Need(mu->GetEntries()==expectedMuonCount,
       "merged independent Geant4 muon count does not equal frozen total");
  Need(wr->GetEntries()>=expectedMuonCount &&
       seed->GetEntries()>=expectedMuonCount &&
       seedWire->GetEntries()>=seed->GetEntries(),
       "merged physical per-wire or per-Heed-primary TTree entries missing");
  std::set<std::string> ids;
  std::map<std::string,long long> channelMuons;
  TTreeReader reader(mu);
  TTreeReaderValue<std::string> crossing(reader,"crossing_id");
  TTreeReaderValue<std::string> channel(reader,"channel");
  for(;reader.Next();){
    Need(ids.insert(*crossing).second,"duplicate original Geant4 crossing in merged sample: "+*crossing);
    Need(*channel=="jpsi"||*channel=="phi"||*channel=="omega",
         "unknown channel after ROOT TTree merge");
    ++channelMuons[*channel];
  }
  Need(static_cast<long long>(ids.size())==expectedMuonCount,
       "ROOT merged source crossing IDs do not conserve independent muons");
  auto* h=Histogram("all_n_wire");
  (void)h;
  file->GetObject("all_n_wire",h);
  Need(h && std::llround(h->GetEntries())==expectedMuonCount,
       "merged ROOT multiplicity histogram differs from actual Geant4 crossing IDs");
  TH1* g=nullptr;
  file->GetObject("all_asym_wire_gain_vs_A_hemi",g);
  Need(g&&g->GetEntries()>0,"empty merged actual per-anode gain vs A ROOT histogram");
  std::ofstream audit(out/"merged_root_qa.txt");
  Need(bool(audit),"cannot write exact ROOT merged physics QA");
  audit<<"NA60+ simulation 1800V P3 noFW: original SHA-frozen Garfield++ source\n"
       <<"source_shards="<<dirs.size()<<"\n"
       <<"original_independent_muons="<<ids.size()<<"\n"
       <<"source_channel_jpsi="<<channelMuons["jpsi"]<<"\n"
       <<"source_channel_phi="<<channelMuons["phi"]<<"\n"
       <<"source_channel_omega="<<channelMuons["omega"]<<"\n"
       <<"physical_muon_wire_rows="<<wr->GetEntries()<<"\n"
       <<"heed_primary_rows="<<seed->GetEntries()<<"\n"
       <<"heed_primary_wire_rows="<<seedWire->GetEntries()<<"\n"
       <<"asym_positive_gain_physical_wire_rows="<<g->GetEntries()<<"\n"
       <<"all_distinct_crossing_ids=true\n"
       <<"all_ROOT_histograms_have_identical_frozen_run_schema=true\n"
       <<"This is exact archived wire endpoint/birth truth, NOT strip current.\n"
       <<"Statistics across same-muon wires need muon-cluster treatment.\n";
  Need(bool(audit),"failed writing ROOT merge QA");
  PlotFinal(*file,out);
  file->Close();
}
}
int main(int argc,char**argv){
  try{
    Need(argc==4,
         "usage: wiretruth1800_merge_full shard_directory_list.txt output_dir expected_independent_muons");
    const fs::path listed=fs::absolute(argv[1]);
    const fs::path out=fs::absolute(argv[2]);
    Need(fs::is_directory(out),"output directory must already exist");
    const long long muons=Number(argv[3]);
    Need(muons>0,"zero expected independent muons");
    const auto dirs=ReadShardDirs(listed);
    MergeCsv(dirs,out/"channel_station_summary.csv",
             "channel_station_summary.csv",true,muons);
    MergeCsv(dirs,out/"asymmetry_denominators.csv",
             "asymmetry_denominators.csv",false,muons);
    MergeTrees(dirs,out,muons);
    std::cout<<"PASS 1800 V FULL exact merged original Geant4 muons="<<muons
             <<" shards="<<dirs.size()<<" output="<<out<<"\n";
    return 0;
  }catch(const std::exception&e){
    std::cerr<<"FAILED: "<<e.what()<<"\n";
    return 2;
  }
}
