// Extend a COPY of the verified 1800 V 5620-muon ROOT with exact Heed records.
// Does not calculate electron drift paths. No Garfield++ rerun.
// Input normalized TSV comes from SHA-audited original Heed/seed archives.
#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TObject.h>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
void Check(bool b,const std::string& s){if(!b)throw std::runtime_error(s);}
std::string Key(const std::string& s,int i) {
  return s+"\x1f"+std::to_string(i);
}
std::vector<std::string> Split(const std::string& line) {
  std::vector<std::string> r;std::string s;std::stringstream st(line);
  while(std::getline(st,s,'\t'))r.push_back(s);
  if(!line.empty()&&line.back()=='\t')r.push_back("");
  return r;
}
int Int(const std::string& s){
  std::size_t n=0;const long v=std::stol(s,&n);
  Check(n==s.size()&&v>=std::numeric_limits<int>::min()&&
        v<=std::numeric_limits<int>::max(),"invalid integer '"+s+"'");
  return static_cast<int>(v);
}
long long Long(const std::string& s) {
  std::size_t n=0;const auto v=std::stoll(s,&n);
  Check(n==s.size(),"invalid long integer '"+s+"'");
  return v;
}
double Real(const std::string&s) {
  std::size_t n=0;const double v=std::stod(s,&n);
  Check(n==s.size()&&std::isfinite(v),"invalid float '"+s+"'");
  return v;
}
std::ifstream OpenTSV(const std::string& file,const std::vector<std::string>& header){
  std::ifstream f(file);
  Check(bool(f),"cannot open "+file);
  std::string line;Check(bool(std::getline(f,line)),"empty TSV "+file);
  Check(Split(line)==header,"wrong TSV columns "+file);
  return f;
}
void TypeBranch(TTree& t,const char* name,int* p){
  t.Branch(name,p,(std::string(name)+"/I").c_str());
}
void TypeBranch(TTree& t,const char* name,Long64_t* p){
  t.Branch(name,p,(std::string(name)+"/L").c_str());
}
void TypeBranch(TTree& t,const char* name,double* p){
  t.Branch(name,p,(std::string(name)+"/D").c_str());
}
}
int main(int argc,char** argv) {
  try {
    Check(argc==5,"usage: enrich_heed_root ROOT_COPY CLUSTERS_TSV PRIMARIES_TSV EXPECTED_MUONS");
    const std::string path=argv[1];
    const long long expected=Long(argv[4]);
    Check(expected==5620,"this frozen-export converter is only for 5620 muons");
    std::unique_ptr<TFile> file(TFile::Open(path.c_str(),"UPDATE"));
    Check(file&&!file->IsZombie(),"cannot update ROOT copy");
    TTree *mu=nullptr,*seed=nullptr;
    file->GetObject("muons",mu);file->GetObject("seeds",seed);
    Check(mu&&seed&&mu->GetEntries()==expected,"existing ROOT muons/seeds absent or mismatched");
    Check(!file->Get("heed_clusters")&&!file->Get("heed_primaries"),
          "output already enriched; refuse replacement");

    std::map<std::string,long long> expectedPrim, seenPrim;
    {
      TTreeReader rdr(mu);
      TTreeReaderValue<std::string> crossing(rdr,"crossing_id");
      TTreeReaderValue<Long64_t> count(rdr,"heed_primaries");
      for(;rdr.Next();) {
        Check(*count>=0&&expectedPrim.emplace(*crossing,*count).second,
              "invalid/duplicate original muon crossing_id");
      }
    }
    std::set<std::string> seedKeys;
    {
      TTreeReader rdr(seed);
      TTreeReaderValue<std::string> crossing(rdr,"crossing_id");
      TTreeReaderValue<Int_t> primary(rdr,"primary_index");
      for(;rdr.Next();) {
        Check(expectedPrim.count(*crossing)&&*primary>=0,
              "original seed without original muon");
        Check(seedKeys.insert(Key(*crossing,*primary)).second,
              "duplicate original ROOT seed");
      }
    }
    const auto expectedSeeds=seedKeys.size();

    std::string crossing;
    int index=0,cluster_index=0,n_original=0,n_accepted=0;
    double u=0,v=0,w=0,t=0;
    TTree clusters("heed_clusters","Exact archived TrackHeed cluster starts");
    clusters.Branch("crossing_id",&crossing);
    TypeBranch(clusters,"cluster_index",&cluster_index);
    TypeBranch(clusters,"u_rel_cm",&u);TypeBranch(clusters,"v_abs_cm",&v);
    TypeBranch(clusters,"w_rel_cm",&w);TypeBranch(clusters,"t_ns",&t);
    TypeBranch(clusters,"n_original",&n_original);
    TypeBranch(clusters,"n_accepted",&n_accepted);
    std::map<std::string,std::pair<int,int>> clusterCounts; // wanted,actual
    std::map<std::string,long long> clustersPerMuon;
    {
      auto in=OpenTSV(argv[2],{"crossing_id","cluster_index",
        "u_rel_cm","v_abs_cm","w_rel_cm","t_ns","n_original","n_accepted"});
      std::string line;
      while(std::getline(in,line)){
        if(line.empty())continue;
        auto x=Split(line);Check(x.size()==8,"malformed cluster row");
        crossing=x[0];cluster_index=Int(x[1]);
        u=Real(x[2]);v=Real(x[3]);w=Real(x[4]);t=Real(x[5]);
        n_original=Int(x[6]);n_accepted=Int(x[7]);
        Check(expectedPrim.count(crossing)&&cluster_index>=0&&
              n_original>=n_accepted&&n_accepted>=0,
              "invalid Heed cluster payload");
        Check(clusterCounts.emplace(Key(crossing,cluster_index),
              std::make_pair(n_accepted,0)).second,"duplicate cluster ID");
        clustersPerMuon[crossing]++;
        clusters.Fill();
      }
    }
    // The source stores Heed's local coordinate frame. Global u includes a
    // patch/entrance phase and should not be reconstructed without input G4.
    double raw_energy=0,raw_dx=0,raw_dy=0,raw_dz=0,
           seed_u=0,seed_v=0,seed_w=0,seed_energy=0;
    int garfield_seed=0,garfield_ne=0,garfield_ni=0,
        successful=0,capped=0;
    Long64_t ions=0,collected=0,attached=0,other=0;
    TTree primaries("heed_primaries",
                    "Exact TrackHeed primary starts and Garfield avalanche counts");
    primaries.Branch("crossing_id",&crossing);
    TypeBranch(primaries,"primary_index",&index);
    TypeBranch(primaries,"cluster_index",&cluster_index);
    TypeBranch(primaries,"u_rel_cm",&u);
    TypeBranch(primaries,"v_abs_cm",&v);
    TypeBranch(primaries,"w_rel_cm",&w);
    TypeBranch(primaries,"t_ns",&t);
    TypeBranch(primaries,"raw_heed_energy_eV",&raw_energy);
    TypeBranch(primaries,"raw_dx",&raw_dx);
    TypeBranch(primaries,"raw_dy",&raw_dy);
    TypeBranch(primaries,"raw_dz",&raw_dz);
    TypeBranch(primaries,"seed_u_mm",&seed_u);
    TypeBranch(primaries,"seed_v_mm",&seed_v);
    TypeBranch(primaries,"seed_w_mm",&seed_w);
    TypeBranch(primaries,"seed_energy_eV",&seed_energy);
    TypeBranch(primaries,"garfield_seed",&garfield_seed);
    TypeBranch(primaries,"garfield_ne",&garfield_ne);
    TypeBranch(primaries,"garfield_ni",&garfield_ni);
    TypeBranch(primaries,"ion_births",&ions);
    TypeBranch(primaries,"collected_e",&collected);
    TypeBranch(primaries,"attached",&attached);
    TypeBranch(primaries,"other",&other);
    TypeBranch(primaries,"successful",&successful);
    TypeBranch(primaries,"capped",&capped);
    long long nPrimary=0;
    {
      auto in=OpenTSV(argv[3],{"crossing_id","primary_index","cluster_index",
        "u_rel_cm","v_abs_cm","w_rel_cm","t_ns","raw_heed_energy_eV",
        "raw_dx","raw_dy","raw_dz","seed_u_mm","seed_v_mm","seed_w_mm",
        "seed_energy_eV","garfield_seed","garfield_ne","garfield_ni",
        "ion_births","collected_e","attached","other","successful","capped"});
      std::string line;
      while(std::getline(in,line)){
        if(line.empty())continue;
        auto x=Split(line);Check(x.size()==24,"malformed primary row");
        crossing=x[0];index=Int(x[1]);cluster_index=Int(x[2]);
        u=Real(x[3]);v=Real(x[4]);w=Real(x[5]);t=Real(x[6]);
        raw_energy=Real(x[7]);raw_dx=Real(x[8]);
        raw_dy=Real(x[9]);raw_dz=Real(x[10]);
        seed_u=Real(x[11]);seed_v=Real(x[12]);
        seed_w=Real(x[13]);seed_energy=Real(x[14]);
        garfield_seed=Int(x[15]);garfield_ne=Int(x[16]);garfield_ni=Int(x[17]);
        ions=Long(x[18]);collected=Long(x[19]);attached=Long(x[20]);
        other=Long(x[21]);successful=Int(x[22]);capped=Int(x[23]);
        Check(expectedPrim.count(crossing)&&index>=0&&cluster_index>=0&&
              ions>=0&&collected>=0&&attached>=0&&other>=0&&
              garfield_seed>=0&&garfield_ne>=0&&garfield_ni>=0&&
              (successful==0||successful==1)&&capped==0,
              "invalid primary source or Garfield counts");
        Check(std::abs(seed_v-10.*v)<1.e-7&&std::abs(seed_w-10.*w)<1.e-7,
              "seed frame not consistent with cached Heed v,w");
        const auto key=Key(crossing,index);
        Check(seedKeys.erase(key)==1,"Heed primary not found in ROOT seeds");
        const auto ck=Key(crossing,cluster_index);
        auto it=clusterCounts.find(ck);
        Check(it!=clusterCounts.end(),"primary has no archived parent cluster");
        it->second.second++;
        seenPrim[crossing]++;
        primaries.Fill();nPrimary++;
      }
    }
    Check(seedKeys.empty()&&nPrimary==static_cast<long long>(expectedSeeds),
          "not all original ROOT seeds have matching Heed rows");
    Check(seenPrim.size()==expectedPrim.size(),"missing muon primary populations");
    for(const auto&[id,n]:expectedPrim)
      Check(seenPrim[id]==n,"per-muon Heed count differs from original ROOT");
    for(const auto&[id,p]:clusterCounts)
      Check(p.first==p.second,"cluster n_accepted differs from actual electrons");
    file->cd();
    Check(clusters.Write()>0&&primaries.Write()>0,"ROOT tree write failed");
    file->Write();file->Close();
    std::cout<<"PASS_HEED_ENRICHED muons="<<expected
             <<" clusters="<<clusterCounts.size()
             <<" primaries="<<nPrimary
             <<" original_root_seeds="<<expectedSeeds
             <<" original_root_objects_preserved=true\n";
    return 0;
  } catch(const std::exception& e) {
    std::cerr<<"FAILED_HEED_ENRICHMENT: "<<e.what()<<"\n";
    return 2;
  }
}
