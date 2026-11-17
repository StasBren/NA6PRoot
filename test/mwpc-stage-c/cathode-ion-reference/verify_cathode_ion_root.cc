// ROOT/C++ physics-integrity audit for one uncapped 2100 V ion-birth archive.
// No statistical information is inferred from unsealed or failed jobs.
// Usage: verify_cathode_ion_root ROOTFILE MICROSCOPIC_PREFIX
#include <TFile.h>
#include <TTree.h>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
namespace {
void Require(bool b,const std::string& why) {
  if(!b)throw std::runtime_error("ROOT PHYSICS AUDIT FAILED: "+why);
}
std::vector<std::string> Split(const std::string& s){
  std::stringstream ss(s);std::string x;std::vector<std::string> out;
  while(std::getline(ss,x,','))out.push_back(x);
  return out;
}
std::map<std::string,std::size_t> Columns(const std::string& header) {
  std::map<std::string,std::size_t> m;
  auto fields=Split(header);
  for(std::size_t i=0;i<fields.size();++i) {
    auto key=fields[i];
    if(!key.empty()&&key.back()=='\r')key.pop_back();
    Require(m.emplace(key,i).second,"duplicate CSV column "+key);
  }
  return m;
}
std::string Col(const std::vector<std::string>& row,
                const std::map<std::string,std::size_t>& names,
                const std::string& key) {
  auto it=names.find(key);
  Require(it!=names.end()&&it->second<row.size(),"missing CSV value "+key);
  return row[it->second];
}
struct Census {
  long long n=0,saved=0,collected=0,endpoints=0,rootIonRows=0;
  bool successful=false;
};
}
int main(int argc,char**argv) {
  try {
    Require(argc==3,"usage verify_cathode_ion_root FILE.root MICROSCOPIC_PREFIX");
    std::unique_ptr<TFile> file(TFile::Open(argv[1],"READ"));
    Require(file&&!file->IsZombie(),"cannot open ROOT input");
    auto tree=[&](const char* name)->TTree*{
      auto* t=dynamic_cast<TTree*>(file->Get(name));
      Require(t!=nullptr,std::string("missing TTree ")+name);
      return t;
    };
    auto* muons=tree("muons");
    auto* primaries=tree("primaries");
    auto* ions=tree("ions");
    auto* ends=tree("electron_endpoints");
    auto* electron=tree("electron_control");
    tree("clusters");
    Require(muons->GetEntries()==1,"single-muon ROOT archive requires exactly one row");
    int st=0,det=0,base=0,reservoir=0;double pitch=0,phase=0,patchShift=0,entryW=0;
    std::string *mode=nullptr,*source=nullptr,*crossingId=nullptr;
    muons->SetBranchAddress("storage_mode",&mode);
    muons->SetBranchAddress("source_channel",&source);
    muons->SetBranchAddress("crossing_id",&crossingId);
    muons->SetBranchAddress("station",&st);
    muons->SetBranchAddress("det_id",&det);
    muons->SetBranchAddress("base_global_wire",&base);
    muons->SetBranchAddress("reservoir_per_primary",&reservoir);
    muons->SetBranchAddress("pitch_cm",&pitch);
    muons->SetBranchAddress("station_phase_mm",&phase);
    muons->SetBranchAddress("patch_shift_cm",&patchShift);
    muons->SetBranchAddress("w_in_mm",&entryW);
    muons->GetEntry(0);
    Require(*mode=="weighted"||*mode=="full","unrecognized archive mode");
    Require(st>=0&&st<6&&det>=0&&pitch>0&&std::abs(pitch-.3)<1e-9,
            "geometry/station not P3");
    Require(std::abs(phase-(st==0?1.5:0.))<1e-9,
            "station-aware physical wire parity mismatch");
    Require(std::abs(patchShift-.1*(phase+base*10*pitch))<1e-8,
            "physical origin of nine-wire patch drifted");
    Require(!source->empty()&&!crossingId->empty(),"lost muon identity");

    int pri=0,succ=0;long long total=0,saved=0,collected=0,endpoints=0;
    primaries->SetBranchAddress("primary_index",&pri);
    primaries->SetBranchAddress("ion_births_exact",&total);
    primaries->SetBranchAddress("ion_births_saved",&saved);
    primaries->SetBranchAddress("collected_electrons",&collected);
    primaries->SetBranchAddress("electron_endpoints",&endpoints);
    primaries->SetBranchAddress("successful",&succ);
    std::map<int,Census> census;
    long long totalIons=0,totalSaved=0,totalCollected=0;
    for(Long64_t i=0;i<primaries->GetEntries();++i){
      primaries->GetEntry(i);
      Require(pri==i&&total>=0&&saved>=0&&saved<=total&&collected>=0,
              "invalid primary order/census");
      Require(succ==1,"failed primary in ROOT archive");
      if(*mode=="full")Require(saved==total,"FULL mode silently sampled ions");
      else Require((reservoir==256||reservoir==512||reservoir==1024) &&
                     saved==std::min<long long>(total,reservoir),
                   "weighted reservoir size disagrees with exact ion census");
      census.emplace(pri,Census{total,saved,collected,endpoints,0,true});
      totalIons+=total;totalSaved+=saved;totalCollected+=collected;
    }
    Require(primaries->GetEntries()>0,"no Heed primaries in detailed event");

    int ipri=0,wire=0,local=0;
    long long ordinal=0;
    double u=0,v=0,w=0,t=0,ug=0,wg=0,wireMm=0,pi=0,weight=0;
    ions->SetBranchAddress("primary_index",&ipri);
    ions->SetBranchAddress("callback_ordinal",&ordinal);
    ions->SetBranchAddress("local_wire_index",&local);
    ions->SetBranchAddress("global_wire_index",&wire);
    ions->SetBranchAddress("u_patch_cm",&u);
    ions->SetBranchAddress("v_absolute_cm",&v);
    ions->SetBranchAddress("w_relative_cm",&w);
    ions->SetBranchAddress("t_birth_ns",&t);
    ions->SetBranchAddress("u_chamber_cm",&ug);
    ions->SetBranchAddress("w_chamber_cm",&wg);
    ions->SetBranchAddress("wire_u_global_mm",&wireMm);
    ions->SetBranchAddress("inclusion_probability",&pi);
    ions->SetBranchAddress("inverse_inclusion_weight",&weight);
    int previousPrimary=-1;long long prevOrdinal=-1;
    std::set<long long> distinctOrdinals;
    long long currentPrimaryCount=0;
    auto checkGroup=[&](){
      if(previousPrimary>=0)
        Require(currentPrimaryCount==census.at(previousPrimary).saved,
                "sample/full ROOT rows differ from primary census");
    };
    for(Long64_t i=0;i<ions->GetEntries();++i){
      ions->GetEntry(i);
      Require(census.count(ipri),"ion row refers to unknown primary");
      if(ipri!=previousPrimary){
        checkGroup();
        Require(ipri>previousPrimary,"ion rows are not ordered by primary");
        previousPrimary=ipri;prevOrdinal=-1;distinctOrdinals.clear();
        currentPrimaryCount=0;
      }
      const auto& c=census.at(ipri);
      Require(ordinal>=0&&ordinal<c.n,"invalid ion callback ordinal");
      if(*mode=="full"){
        Require(ordinal==prevOrdinal+1,"FULL ion callback stream has a gap");
      }else{
        Require(distinctOrdinals.insert(ordinal).second,
                "reservoir repeated a callback ordinal");
      }
      prevOrdinal=ordinal;
      Require(std::isfinite(u)&&std::isfinite(v)&&std::isfinite(w)&&
              std::isfinite(t)&&std::abs(v)<=.3000001 &&
              std::abs(w)<5.&&std::abs(u)<1.5,
              "unphysical ion birth coordinate/time");
      Require(local>=-4&&local<=4&&wire==base+local,
              "ion birth nearest physical anode drift");
      Require(std::abs(ug-(patchShift+u))<1e-9 &&
              std::abs(wg-(entryW*.1+w))<1e-9 &&
              std::abs(wireMm-(phase+wire*10*pitch))<1e-8,
              "physical coordinate transform drift");
      const double expectedPi=double(c.saved)/double(c.n);
      Require(std::abs(pi-expectedPi)<1e-12 &&
              std::abs(weight-double(c.n)/c.saved)<1e-10*
                    std::max(1.,double(c.n)/c.saved),
              "non-inverse-inclusion-probability ion weight");
      census[ipri].rootIonRows++;
      currentPrimaryCount++;
    }
    checkGroup();
    Require(ions->GetEntries()==totalSaved,"total ROOT ion entries != exact saved rows");
    for(const auto& [idx,c]:census)
      Require(c.rootIonRows==c.saved,"primary ion count incomplete or corrupted");
    if(*mode=="full"){
      Require(ends->GetEntries()==census.at(0).endpoints,
              "full first primary electron endpoints incomplete");
      Require(electron->GetEntries()==4000,
              "no direct 0.05ns whole-cathode electron signal reference");
    }else{
      Require(ends->GetEntries()==0&&electron->GetEntries()==0,
              "weighted mode pretends to contain full electron signal");
    }

    // Independent exact seed and wire census against native original CSVs.
    std::ifstream seeds(std::string(argv[2])+"_seed.csv");
    Require(seeds.good(),"missing original microscopic_seed.csv");
    std::string line;std::getline(seeds,line);const auto sc=Columns(line);
    long long csvSeedIons=0,csvSeedElec=0;int csvSeedCount=0;
    while(std::getline(seeds,line)){
      if(line.empty())continue;
      auto row=Split(line);
      int j=std::stoi(Col(row,sc,"primary_index"));
      Require(census.count(j),"CSV primary not in ROOT");
      const auto& c=census.at(j);
      const auto si=std::stoll(Col(row,sc,"ion_births"));
      const auto se=std::stoll(Col(row,sc,"collected_e"));
      Require(si==c.n&&se==c.collected,"ROOT vs microscopic seed charge mismatch");
      Require(Col(row,sc,"capped")=="0" &&
              Col(row,sc,"successful")=="1",
              "capped or failed avalanche in microscopic CSV");
      csvSeedIons+=si;csvSeedElec+=se;++csvSeedCount;
    }
    Require(csvSeedIons==totalIons && csvSeedElec==totalCollected &&
            csvSeedCount==(int)census.size(),"CSV and ROOT seed sums disagree");
    std::ifstream wireSeeds(std::string(argv[2])+"_wire_seed.csv");
    Require(wireSeeds.good(),"missing exact physical wire seed CSV");
    std::getline(wireSeeds,line);const auto wc=Columns(line);
    std::map<int,std::pair<long long,long long>> wireSums;
    while(std::getline(wireSeeds,line)){
      if(line.empty())continue;
      auto row=Split(line);
      int j=std::stoi(Col(row,wc,"primary_index"));
      int global=std::stoi(Col(row,wc,"wire_index_global"));
      Require(global>=base-4&&global<=base+4,"wire outside patch");
      wireSums[j].first+=std::stoll(Col(row,wc,"ion_births"));
      wireSums[j].second+=std::stoll(Col(row,wc,"collected_e"));
    }
    for(const auto& [idx,c]:census)
      Require(wireSums[idx].first==c.n &&
              wireSums[idx].second==c.collected,
              "per-wire electron collection / ion birth census mismatch");
    const double bytes=double(std::filesystem::file_size(argv[1]));
    std::cout<<"PASS CATHODE ROOT PHYSICS AUDIT\n"
             <<"channel="<<*source<<" station="<<st<<" det_id="<<det
             <<" mode="<<*mode<<" Heed_primaries="<<census.size()
             <<" true_ion_births="<<totalIons
             <<" saved_ion_rows="<<totalSaved
             <<" collected_electrons="<<totalCollected
             <<" root_bytes="<<(long long)bytes
             <<" bytes_per_saved_ion="<<(totalSaved?bytes/totalSaved:0.)
             <<"\n";
    return 0;
  }catch(const std::exception& e){
    std::cerr<<e.what()<<"\n";return 2;
  }
}
