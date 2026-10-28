// REAL native Geant4 crossing -> one common TrackHeed cluster population
// -> four independent full microscopic Garfield avalanche geometries.
// B1 is an ionization-BIRTH count, not induced cathode signal.
// No old 2+4mm code, original ROOT files or frozen LUTs are modified.
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/ComponentAnalyticField.hh"
#include "Garfield/ComponentConstant.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Random.hh"
#include "Garfield/RandomEngineRoot.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/TrackHeed.hh"

namespace {
constexpr double PI=3.1415926535897932384626433832795;
constexpr double CM_PER_MM=.1;
constexpr double ANGLE_EPS=1.e-12;

struct Crossing {
  std::string source;
  int chunk=0,event=0,hit=0,track=0,station=0,det=0,pdg=0;
  std::array<double,3> a{},b{},p{};
  double thetaU=0,thetaW=0,length=0,edep=0,momentum=0;
  unsigned int heedSeed=0;
};
struct Primary {
  double u=0,v=0,w=0,t=0,energy=0;
  double dx=0,dy=0,dz=0;
};
struct Counts {
  std::uint64_t up=0,down=0,boundary=0,top=0,bottom=0,side=0;
  void Add(double du,double v) {
    const double phi=std::atan2(v,du)*180./PI;
    if (phi>0 && phi<180.)++up;
    else if (phi<0 && phi>-180.)++down;
    else ++boundary;
    if(phi>=45.&&phi<135.)++top;
    else if(phi>=-135.&&phi<-45.)++bottom;
    else ++side;
  }
  void Merge(const Counts& o) {
    up+=o.up;down+=o.down;boundary+=o.boundary;
    top+=o.top;bottom+=o.bottom;side+=o.side;
  }
};
std::string Asym(std::uint64_t a,std::uint64_t b) {
  if(a+b==0)return ""; // Empty denominators are UNDEFINED.
  std::ostringstream s;
  s<<std::setprecision(12)<<(double(a)-double(b))/double(a+b);
  return s.str();
}
struct BirthAccumulator {
  Counts counts;
  double pitch=0;
  std::uint64_t total=0;
  // Independent reservoir RNG; MUST NOT consume Garfield microscopic RNG.
  std::mt19937_64 sampleRng;
  std::vector<std::array<double,4>> sample;
  BirthAccumulator(unsigned int seed,double p)
    :pitch(p),sampleRng(std::uint64_t(seed)+0x9e3779b97f4a7c15ULL){}
  void Add(double u,double v,double w,double t) {
    counts.Add(u-std::round(u/pitch)*pitch,v);
    ++total;
    const std::array<double,4> birth{{u,v,w,t}};
    constexpr std::size_t limit=100;
    if(sample.size()<limit)sample.push_back(birth);
    else {
      std::uniform_int_distribution<std::uint64_t> choose(0,total-1);
      const auto j=choose(sampleRng);
      if(j<limit)sample[static_cast<std::size_t>(j)]=birth;
    }
  }
};
BirthAccumulator* gBirth=nullptr;
// Garfield stores a process-wide RNG pointer. Keep the current engine alive
// across every field/sensor allocation, not just inside a function scope.
std::unique_ptr<Garfield::RandomEngineRoot> gEngine;
void SeedGarfield(unsigned int n){
  gEngine=std::make_unique<Garfield::RandomEngineRoot>(n);
  Garfield::Random::SetEngine(*gEngine);
}
void OnIonisation(const double x,const double y,const double z,
                  const double t,const int,const int,Garfield::Medium*) {
  if(gBirth)gBirth->Add(x,y,z,t);
}
double Numeric(const std::string& s) {
  if(s.empty())throw std::runtime_error("Empty numeric column in native input");
  const double x=std::stod(s);
  if(!std::isfinite(x))throw std::runtime_error("Non-finite native input");
  return x;
}
int Integer(const std::string& s){return std::stoi(s);}
std::vector<std::string> Split(const std::string& line) {
  std::vector<std::string> out;
  std::stringstream ss(line);
  std::string item;
  while(std::getline(ss,item,','))out.push_back(item);
  return out;
}
std::vector<Crossing> ReadCrossings(const std::string& file) {
  std::ifstream in(file);
  if(!in)throw std::runtime_error("Cannot open immutable native selection: "+file);
  const std::string expected=
   "source,chunk,event,hit_index,track_id,station,det_id,pdg,"
   "u_in_mm,v_in_mm,w_in_mm,u_out_mm,v_out_mm,w_out_mm,"
   "theta_u_in_deg,theta_w_in_deg,path_length_mm,edep_eV,"
   "pu_in_GeV,pv_in_GeV,pw_in_GeV,momentum_GeV,heed_seed";
  std::string header;
  if(!std::getline(in,header))
    throw std::runtime_error("Empty native Geant4 muon CSV header");
  // Python's csv.DictWriter uses CRLF by default, even on Linux.
  // std::getline removes LF but retains CR. Accept CRLF and LF.
  if(!header.empty() && header.back()=='\r')header.pop_back();
  if(header!=expected)
    throw std::runtime_error("Wrong real Geant4 muon input column schema; got ["+
                             header+"], expected ["+expected+"]");
  std::vector<Crossing> out;
  std::string line;
  while(std::getline(in,line)){
    if(!line.empty() && line.back()=='\r')line.pop_back();
    if(line.empty())continue;
    const auto f=Split(line);
    if(f.size()!=23)throw std::runtime_error("Malformed native crossing CSV row");
    Crossing h;
    h.source=f[0];h.chunk=Integer(f[1]);h.event=Integer(f[2]);h.hit=Integer(f[3]);
    h.track=Integer(f[4]);h.station=Integer(f[5]);h.det=Integer(f[6]);
    h.pdg=Integer(f[7]);
    for(int j=0;j<3;++j){h.a[j]=Numeric(f[8+j]);h.b[j]=Numeric(f[11+j]);}
    h.thetaU=Numeric(f[14]);h.thetaW=Numeric(f[15]);
    h.length=Numeric(f[16]);h.edep=Numeric(f[17]);
    for(int j=0;j<3;++j)h.p[j]=Numeric(f[18+j]);
    h.momentum=Numeric(f[21]);h.heedSeed=std::stoul(f[22]);
    const auto p=std::hypot(h.p[0],h.p[1],h.p[2]);
    const double chord=std::hypot(h.b[0]-h.a[0],
                                 h.b[1]-h.a[1],h.b[2]-h.a[2]);
    double dot=0.;
    for(int j=0;j<3;++j)dot+=(h.b[j]-h.a[j])*h.p[j];
    if(std::abs(h.pdg)!=13 || h.station<0 || h.station>5 ||
       h.momentum<=.11 || std::abs(h.momentum-p)>1.e-5*p ||
       h.length<=.05 || std::abs(h.length-chord)>.005*chord+.02 ||
       dot/(chord*p)<.999 || h.heedSeed<1 ||
       std::max(std::abs(h.a[1]),std::abs(h.b[1]))>3.05)
      throw std::runtime_error("Unphysical native muon or direction mismatch");
    out.push_back(h);
  }
  if(out.empty()||out.size()>24)
    throw std::runtime_error("B1 accepts 1..24 real verified muon crossings only");
  return out;
}
struct Geometry {std::string label;double pitchMm;double fwUm;};
const std::array<Geometry,4> GEOMS{{{"p3_noFW",3.,0.},{"p3_FW75",3.,75.},
                                    {"p4_noFW",4.,0.},{"p4_FW75",4.,75.}}};
int NAnodes(int st,double pitch) {
  const double width=st==0?307.2:512.;
  return int(std::floor(width/pitch+.5)); // EXACT existing frozen map rule
}
double EntrancePhase(double uMm,int station,double pitch) {
  const int n=NAnodes(station,pitch);
  const double phase=n%2==0?pitch*.5:0.;
  const double nearest=std::round((uMm-phase)/pitch)*pitch+phase;
  return uMm-nearest;
}
// The primary-ionisation source is an independent, gas-filled box.
 // Do NOT model it using two identical-potential ComponentAnalyticField
 // planes: Garfield's CellCheck rejects a cell with no field solution,
 // then TrackHeed cannot find an ionisable medium at the start position.
struct HeedSourceVolume {
  Garfield::ComponentConstant field;
  Garfield::Sensor sensor;
  std::array<double,3> start{}, direction{};
  explicit HeedSourceVolume(Garfield::MediumMagboltz& gas,
                            const Crossing& c) {
    const double du=c.b[0]-c.a[0],dv=c.b[1]-c.a[1],dw=c.b[2]-c.a[2];
    const double maxU=std::max(1.,std::abs(du)*.1+.2);
    const double maxW=std::max(1.,std::abs(dw)*.1+.2);
    // The component AND the sensor require their own active bounds.
    field.SetArea(-maxU,-.3,-maxW,maxU,+.3,+maxW);
    field.SetMedium(&gas);
    field.SetElectricField(0.,1.,0.);
    sensor.AddComponent(&field);
    sensor.SetArea(-maxU,-.3,-maxW,maxU,+.3,+maxW);
    direction={{c.p[0]/c.momentum,c.p[1]/c.momentum,
                c.p[2]/c.momentum}};
    // Define gas in local u around zero so we can use the original
    // muon entry with only a *translation*, never a random angle.
    constexpr double inwardCm=1.e-4; // 1 micron inside gas
    start={{inwardCm*direction[0], .1*c.a[1]+inwardCm*direction[1],
            inwardCm*direction[2]}};
    // Sensor.hh in the user's installed Garfield version declares
    // Medium* GetMedium(double x, double y, double z), not a 4-argument
    // bool-plus-output-pointer API.
    Garfield::Medium* medium=sensor.GetMedium(start[0],start[1],start[2]);
    if(!medium || !medium->IsIonisable()) {
      std::ostringstream details;
      details<<"Heed primary source is not ionisable at muon start: "
             <<"y="<<start[1]<<"cm (cathodes ±0.3cm), "
             <<"native v_in="<<c.a[1]<<"mm, direction v="<<direction[1];
      throw std::runtime_error(details.str());
    }
  }
};
std::vector<Primary> HeedPrimaries(Garfield::MediumMagboltz& gas,
                                  const Crossing& c,int& nClusters,
                                  int& nDiscarded,int& nFailed){
  // A geometry-independent HEED PRIMARY source is generated ONCE per
  // original crossing, inside a valid ionisable gas volume. This weak
  // placeholder uniform field defines medium lookup ONLY; it is NOT
  // used for Garfield avalanche transport, which uses each chamber's
  // actual 1800V analytic wire field.
  HeedSourceVolume volume(gas,c);
  Garfield::TrackHeed heed;
  heed.SetSensor(&volume.sensor);
  heed.SetParticle(c.pdg==13?"mu-":"mu+");
  heed.SetMomentum(c.momentum*1.e9);
  heed.EnableDeltaElectronTransport();

  SeedGarfield(c.heedSeed);
  // Use the actual Geant4 muon momentum and the verified start location.
  const auto& d=volume.direction;
  const auto& start=volume.start;
  if(!heed.NewTrack(start[0],start[1],start[2],0.,d[0],d[1],d[2]))\n    throw std::runtime_error("TrackHeed::NewTrack returned false for verified native muon");

  std::vector<Primary> primaries;
  double xc=0,yc=0,zc=0,tc=0,ec=0,extra=0;
  int nc=0;
  nClusters=nDiscarded=nFailed=0;
  while(heed.GetCluster(xc,yc,zc,tc,nc,ec,extra)){
    // TrackHeed propagates to the SENSOR boundary, not necessarily the
    // original gas-hit OUT point. Clip along actual native Geant4 chord.
    const double along=10.*(xc*d[0]+(yc-.1*c.a[1])*d[1]+zc*d[2]);
    if(along<-.005||along>c.length+.005){++nDiscarded;continue;}
    ++nClusters;
    for(int i=0;i<nc;++i){
      Primary x;
      if(!heed.GetElectron(i,x.u,x.v,x.w,x.t,x.energy,x.dx,x.dy,x.dz)){
        ++nFailed;continue;
      }
      if(x.v<=-.3||x.v>=.3){++nDiscarded;continue;}
      primaries.push_back(x);
    }
  }
  return primaries;
}
struct Response {
  Counts ions,electrons;
  std::uint64_t births=0,ionized=0, eCollected=0, attached=0,other=0;
  int zero=0,capped=0,failed=0;
};
void RunGeometry(Garfield::MediumMagboltz& gas,const Geometry& geom,
                 const Crossing& h,const std::vector<Primary>& primaries,
                 const std::string& crossingLabel,
                 std::ofstream& seedOut,std::ofstream& birthOut,
                 std::ofstream& hitOut,int clusters,int dropped,int heedFailed) {
  const double p=.1*geom.pitchMm;
  const double wireRadius=.0015; // 30 um diameter = 0.0015 cm radius
  // For the finite centered anode grid, translate the Geant4 physical
  // entrance u to a central periodic 9-anode patch without changing the
  // electron's offset from its nearest actual wire.
  const double phaseCm=.1*EntrancePhase(h.a[0],h.station,geom.pitchMm);
  const double fw=geom.fwUm*1.e-4;
  Garfield::ComponentAnalyticField field;
  field.SetMedium(&gas);
  constexpr int half=4;
  for(int i=-half;i<=half;++i)field.AddWire(i*p,0.,.003,1800.,"anode");
  if(fw>0)for(int i=-half;i<half;++i)
    field.AddWire((i+.5)*p,0.,fw,0.,"field_wire");
  field.AddPlaneY(-.3,0.,"cathode_minus");
  field.AddPlaneY(+.3,0.,"cathode_plus");
  field.SetMagneticField(0.,0.,0.);
  Garfield::Sensor sensor;
  sensor.AddComponent(&field);
  sensor.SetArea(-(half+.5)*p,-.3,-5.,(half+.5)*p,.3,5.);
  Response tot;
  for(std::size_t i=0;i<primaries.size();++i) {
    const auto& e=primaries[i];
    const double u=e.u+phaseCm;
    if(std::abs(u)>4.2*p||std::abs(e.w)>4.8)
      throw std::runtime_error("True Heed electron exceeds validated analytic periodic wire patch");
    const std::uint64_t raw=std::uint64_t(h.heedSeed)+
                            std::uint64_t(i+1)*104729u;
    const unsigned int rngSeed=1+unsigned(raw%2'000'000'000ULL);
    SeedGarfield(rngSeed);
    Garfield::AvalancheMicroscopic avalanche(&sensor);
    constexpr unsigned int cap=200000;
    avalanche.EnableAvalancheSizeLimit(cap);
    BirthAccumulator births(rngSeed,p);gBirth=&births;
    avalanche.SetUserHandleIonisation(OnIonisation);
    // With TrackHeed::EnableDeltaElectronTransport(), GetElectron gives
    // low-energy conduction electron POSITIONS and times. Its returned
    // energy and direction are NOT physically meaningful (Garfield guide
    // section 5.1.1). Give these thermalized seeds 0.1eV with isotropic
    // initial direction, consistently across all four candidate fields.
    constexpr double startEnergy=0.1;
    const bool ok=avalanche.AvalancheElectron(u,e.v,e.w,e.t,startEnergy,
                                             0.,0.,0.);
    gBirth=nullptr;
    int ne=0,ni=0;
    avalanche.GetAvalancheSize(ne,ni);
    Counts eCounts;
    std::uint64_t collected=0,attached=0,other=0;
    const auto nEnd=avalanche.GetNumberOfElectronEndpoints();
    for(std::size_t j=0;j<nEnd;++j){
      double x0,y0,z0,t0,E0,x1,y1,z1,t1,E1;int status=0;
      avalanche.GetElectronEndpoint(j,x0,y0,z0,t0,E0,
                                   x1,y1,z1,t1,E1,status);
      const double du=x1-std::round(x1/p)*p;
      if(std::hypot(du,y1)<wireRadius*1.5) {
        ++collected;eCounts.Add(du,y1);
      } else if(status==-7)++attached;
      else ++other;
    }
    const bool capped=ne>=int(cap)||nEnd>=cap;
    const bool zero=collected==0;
    if(ni>0 && births.total==0)
      throw std::runtime_error("Garfield reports avalanche ions but ionization callback recorded none");
    tot.ions.Merge(births.counts);
    tot.electrons.Merge(eCounts);
    tot.births+=births.total;
    tot.ionized+=std::max(ni,0);
    tot.eCollected+=collected;tot.attached+=attached;tot.other+=other;
    tot.zero+=zero;tot.capped+=capped;tot.failed+=!ok;
    seedOut<<crossingLabel<<","<<geom.label<<","<<i<<","<<rngSeed<<","
           <<(10.*u)<<","<<(10.*e.v)<<","<<(10.*e.w)<<","
           <<startEnergy<<","<<ne<<","<<ni<<","<<births.total<<","
           <<collected<<","<<attached<<","<<other<<","<<ok<<","
           <<capped<<","<<births.counts.up<<","<<births.counts.down<<","
           <<births.counts.top<<","<<births.counts.bottom<<","
           <<Asym(births.counts.up,births.counts.down)<<","
           <<Asym(births.counts.top,births.counts.bottom)<<"\n";
    // The full birth census is used for A; this bounded independent-RNG
    // reservoir gives unbiased spatial draws WITHIN each seed avalanche.
    // Weight by N_birth / sample_size for pooling different-sized avalanches.
    const double reservoirWeight=births.sample.empty()?0.:
        double(births.total)/births.sample.size();
    for(const auto& v:births.sample)
      birthOut<<crossingLabel<<","<<geom.label<<","<<i<<","
              <<(10.*v[0])<<","<<(10.*v[1])<<","<<(10.*v[2])<<","
              <<v[3]<<","<<births.total<<","<<reservoirWeight<<"\n";
  }
  hitOut<<crossingLabel<<","<<geom.label<<","<<h.station<<","<<h.det<<","
        <<h.thetaU<<","<<h.thetaW<<","<<h.length<<","<<h.momentum<<","
        <<clusters<<","<<primaries.size()<<","<<dropped<<","<<heedFailed<<","
        <<tot.eCollected<<","<<tot.births<<","<<tot.ionized<<","
        <<tot.zero<<","<<tot.capped<<","<<tot.failed<<","
        <<tot.ions.up<<","<<tot.ions.down<<","<<tot.ions.boundary<<","
        <<tot.ions.top<<","<<tot.ions.bottom<<","<<tot.ions.side<<","
        <<Asym(tot.ions.up,tot.ions.down)<<","
        <<Asym(tot.ions.top,tot.ions.bottom)<<","
        <<tot.electrons.up<<","<<tot.electrons.down<<","
        <<tot.electrons.top<<","<<tot.electrons.bottom<<","
        <<Asym(tot.electrons.up,tot.electrons.down)<<","
        <<Asym(tot.electrons.top,tot.electrons.bottom)<<",";
  if(!primaries.empty())hitOut<<double(tot.eCollected)/primaries.size();
  hitOut<<"\n";
  std::cout<<crossingLabel<<" "<<geom.label
           <<" Heed_e="<<primaries.size()
           <<" ionBirth="<<tot.births
           <<" collected_e="<<tot.eCollected
           <<" capped="<<tot.capped<<"\n";
}
} // namespace

int main(int argc,char**argv) {
  try{
    if(argc==3 && std::string(argv[1])=="--validate-input") {
      const auto crossings=ReadCrossings(argv[2]);
      // The former --validate-input checked only CSV syntax. Also check
      // that the ACTUAL Heed source volume has a valid ionisable medium
      // at every selected original muon's start position.
      Garfield::MediumMagboltz preflightGas;
      preflightGas.SetComposition("ar",70.,"co2",30.);
      preflightGas.SetTemperature(293.15);preflightGas.SetPressure(760.);
      for(const auto& h:crossings){
        HeedSourceVolume preflightVolume(preflightGas,h);
      }
      std::cout<<"PASS: parsed "<<crossings.size()
               <<" original Geant4 muon crossings (CSV CRLF/LF supported)\n"
               <<"PASS: Heed gas source medium is ionisable at all "
               <<crossings.size()<<" muon starts (±3 mm)\n";
      return 0; // No Heed ionisation or Garfield avalanche physics.
    }
    if(argc==3 && std::string(argv[1])=="--validate-heed") {
      const auto crossings=ReadCrossings(argv[2]);
      Garfield::MediumMagboltz gas;
      gas.SetComposition("ar",70.,"co2",30.);
      gas.SetTemperature(293.15);gas.SetPressure(760.);
      int totalClusters=0,totalElectrons=0;
      for(const auto& h:crossings) {
        int clusters=0,discarded=0,failed=0;
        const auto electrons=HeedPrimaries(gas,h,clusters,discarded,failed);
        if(failed)
          throw std::runtime_error("Heed GetElectron failed for native muon "+
                                   std::to_string(h.event)+
                                   "; failed="+std::to_string(failed));
        totalClusters+=clusters;
        totalElectrons+=static_cast<int>(electrons.size());
        std::cout<<"Heed native crossing event="<<h.event
                 <<" station="<<h.station
                 <<" thetaU_deg="<<h.thetaU
                 <<" thetaW_deg="<<h.thetaW
                 <<" p_GeV="<<h.momentum
                 <<" clusters="<<clusters
                 <<" electrons="<<electrons.size()
                 <<" discarded="<<discarded<<"\n";
      }
      if(totalClusters==0||totalElectrons==0)
        throw std::runtime_error("Heed-only native muon preflight yielded no ionisation");
      std::cout<<"PASS: true TrackHeed source generated "<<totalClusters
               <<" clusters and "<<totalElectrons
               <<" low-energy electrons from "<<crossings.size()
               <<" physical Geant4 muon crossings; NO avalanches run\n";
      return 0;
    }
    if(argc!=4)throw std::runtime_error(
      "Usage: mwpc_phase_b_native_muon_3p3 INPUT_NATIVE_CSV OUTPUT_PREFIX MAX_PRIMARY_ELECTRONS | --validate-input INPUT_NATIVE_CSV | --validate-heed INPUT_NATIVE_CSV");
    const std::string input=argv[1],prefix=argv[2];
    const int maxPrimary=std::stoi(argv[3]);
    if(maxPrimary<1||maxPrimary>1000)throw std::runtime_error("max primaries must be 1..1000");
    const auto crossings=ReadCrossings(input);
    Garfield::MediumMagboltz gas;
    gas.SetComposition("ar",70.,"co2",30.);
    gas.SetTemperature(293.15);gas.SetPressure(760.);
    gas.SetMaxElectronEnergy(200.);gas.Initialise(false);
    std::ofstream hitOut(prefix+"_crossing.csv");
    std::ofstream seedOut(prefix+"_seed.csv");
    std::ofstream birthOut(prefix+"_ion_birth_sample.csv");
    if(!hitOut||!seedOut||!birthOut)
      throw std::runtime_error("Cannot open new B1 output CSVs");
    hitOut<<"source|chunk|event|hit_index|track_id|station|det_id,configuration,"
             "station,det_id,theta_u_in_deg,theta_w_in_deg,path_mm,momentum_GeV,"
             "heed_clusters,heed_primaries,heed_discarded,heed_failed,"
             "collected_avalanche_electrons,ion_births,garfield_ion_count,"
             "zero_seeds,capped_seeds,failed_seeds,"
             "ion_up,ion_down,ion_boundary,ion_top,ion_bottom,ion_side,"
             "A_ion_hemi,A_ion_TB,e_up,e_down,e_top,e_bottom,"
             "A_e_hemi,A_e_TB,G_eff_per_Heed_primary\n";
    seedOut<<"crossing_id,configuration,primary_index,garfield_seed,"
             "seed_u_mm,seed_v_mm,seed_w_mm,seed_energy_eV,"
             "garfield_ne,garfield_ni,ion_births,collected_e,"
             "attached,other,successful,capped,"
             "ion_up,ion_down,ion_top,ion_bottom,A_ion_hemi,A_ion_TB\n";
    birthOut<<"crossing_id,configuration,primary_index,u_mm,v_mm,w_mm,t_ns,"
               "seed_ion_births_total,reservoir_weight\n";
    for(const auto& h:crossings){
      const std::string id=h.source+"|"+std::to_string(h.chunk)+"|"+
         std::to_string(h.event)+"|"+std::to_string(h.hit)+"|"+
         std::to_string(h.track)+"|"+std::to_string(h.station)+"|"+
         std::to_string(h.det);
      int nCl=0,discarded=0,failed=0;
      const auto electrons=HeedPrimaries(gas,h,nCl,discarded,failed);
      if(static_cast<int>(electrons.size())>maxPrimary)
        throw std::runtime_error("Heed primary budget exceeded, preserve partial outputs and use a different pilot: "+id);
      // One TRUE Heed electron list used verbatim across all four candidates.
      for(const auto& geom:GEOMS)
        RunGeometry(gas,geom,h,electrons,id,seedOut,birthOut,
                    hitOut,nCl,discarded,failed);
      std::cout<<"PASS genuine native muon "<<id<<" shared Heed primaries="
               <<electrons.size()<<"\n";
    }
    std::cout<<"PASS: muon-driven B1 all crossing/geometry rows written\n";
    return 0;
  }catch(const std::exception& e){
    std::cerr<<"B1 FAILURE: "<<e.what()<<"\n";
    return 2;
  }
}
