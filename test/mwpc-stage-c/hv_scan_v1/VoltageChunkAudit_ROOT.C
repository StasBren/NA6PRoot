// ROOT/C++ physics-integrity gate for one full-microscopic Garfield++ chunk.
// Run with root -l -b -q 'VoltageChunkAudit_ROOT.C("PREFIX","INPUT",20,"PASS_FILE")'
#include <TSystem.h>
#include <TString.h>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <stdexcept>
#include <iostream>
#include <cstdint>
#include <cstdlib>
#include <iomanip>

namespace NA60HV {
using Row=std::map<std::string,std::string>;
std::vector<std::string> Split(const std::string& s,const char delim=','){
  std::vector<std::string> ret; std::string item;std::stringstream ss(s);
  while(std::getline(ss,item,delim))ret.push_back(item);
  if(!s.empty()&&s.back()==delim)ret.emplace_back();
  return ret;
}
std::vector<Row> Load(const std::string& filename){
  std::ifstream f(filename);if(!f)throw std::runtime_error("missing CSV: "+filename);
  std::string line;if(!std::getline(f,line))throw std::runtime_error("empty CSV: "+filename);
  auto cols=Split(line);std::vector<Row> rows;
  while(std::getline(f,line)){
    if(line.empty())continue;
    const auto vals=Split(line);
    if(vals.size()!=cols.size())throw std::runtime_error("malformed CSV "+filename+" line "+std::to_string(rows.size()+2)+" cols "+std::to_string(vals.size())+" != "+std::to_string(cols.size()));
    Row r;for(std::size_t i=0;i<cols.size();++i)r.emplace(cols[i],vals[i]);rows.emplace_back(std::move(r));
  }
  return rows;
}
const std::string& Get(const Row& r,const std::string& key){
  auto it=r.find(key);if(it==r.end())throw std::runtime_error("Missing column "+key);
  return it->second;
}
long long Count(const Row& r,const std::string& key){return std::stoll(Get(r,key));}
}
void VoltageChunkAudit_ROOT(const char* pref,const char* input,int expected,const char* passFile){
  using namespace NA60HV;
  try{
    const auto in=Load(input);
    if((int)in.size()!=expected)throw std::runtime_error("input event count mismatch");
    const auto cr=Load(std::string(pref)+"_crossing.csv");
    const auto seeds=Load(std::string(pref)+"_seed.csv");
    const auto wireSeed=Load(std::string(pref)+"_wire_seed.csv");
    const auto wireMuon=Load(std::string(pref)+"_wire_muon.csv");
    const auto heed=Load(std::string(pref)+"_heed_primaries.csv");
    if((int)cr.size()!=expected)throw std::runtime_error("Garfield crossing count mismatch");
    std::map<std::string,std::pair<long long,long long>> fromSeed,fromMuon,fromWireMuon;
    std::map<std::pair<std::string,long long>,std::pair<long long,long long>> fromWireSeed,fromSeedTruth;
    std::map<std::string,long long> nseed,nheed;
    std::set<std::string> seen;
    long long capped=0,failed=0,ionTotal=0,eTotal=0;
    for(const auto& r:seeds){
      const std::string id=Get(r,"crossing_id");const long long idx=Count(r,"primary_index");
      const long long q=Count(r,"collected_e"),n=Count(r,"ion_births");
      fromSeedTruth[{id,idx}]={q,n};
      if(++nseed[id]>100000)throw std::runtime_error("suspicious Heed seed count");
      fromSeed[id].first+=q;fromSeed[id].second+=n;
      capped+=Count(r,"capped");failed+=Count(r,"successful")==0;
    }
    for(const auto& r:heed)++nheed[Get(r,"crossing_id")];
    for(const auto& r:wireSeed){
      const std::string id=Get(r,"crossing_id");
      auto& p=fromWireSeed[{id,Count(r,"primary_index")}];
      p.first+=Count(r,"collected_e");p.second+=Count(r,"ion_births");
    }
    for(const auto& p:fromSeedTruth){
      if(p.second!=fromWireSeed[p.first])throw std::runtime_error("per-seed wire charge/ion mismatch for "+p.first.first);
    }
    for(const auto& r:wireMuon){
      auto& p=fromWireMuon[Get(r,"crossing_id")];
      p.first+=Count(r,"collected_e");p.second+=Count(r,"ion_births");
    }
    const std::string cidHeader=cr.front().begin()->first;
    // The first column contains 'source|chunk|event|...'; a header string is stable.
    for(const auto& r:cr){
      const auto id=Get(r,"source|chunk|event|hit_index|track_id|station|det_id");
      if(!seen.insert(id).second)throw std::runtime_error("duplicate crossing ID");
      const long long q=Count(r,"collected_avalanche_electrons"),n=Count(r,"ion_births");
      if(fromSeed[id]!=std::make_pair(q,n)||fromWireMuon[id]!=std::make_pair(q,n))
        throw std::runtime_error("per-muon wire/seed sums disagree with Garfield crossing");
      if(nseed[id]!=Count(r,"heed_primaries")||nheed[id]!=nseed[id])
        throw std::runtime_error("Heed primary counts mismatch");
      if(Count(r,"capped_seeds")!=0)throw std::runtime_error("capped seeds in crossing");
      ionTotal+=n;eTotal+=q;
    }
    if(capped)throw std::runtime_error("AVALANCHE SIZE LIMIT: "+std::to_string(capped)+" capped seeds");
    // It is possible for unsuccessful Garfield::AvalancheElectron to have physics outputs.
    // Record those flags; do not silently censor affected events.
    std::ofstream out(passFile);
    if(!out)throw std::runtime_error("Cannot write audit seal");
    out<<"AUDIT_PASS,muons,"<<seen.size()<<",seeds,"<<seeds.size()
       <<",collected_e,"<<eTotal<<",ion_births,"<<ionTotal
       <<",capped,"<<capped<<",unsuccessful,"<<failed<<"\n";
    out.close();
    std::cout<<"PASS ROOT CHUNK AUDIT: "<<seen.size()<<" muons, "
             <<seeds.size()<<" avalanches, "<<eTotal<<" collected e, "
             <<ionTotal<<" ion births, capped="<<capped<<"\n";
  }catch(const std::exception& e){
    std::cerr<<"FATAL ROOT CHUNK AUDIT: "<<e.what()<<"\n";
    gSystem->Exit(12);
  }
}
