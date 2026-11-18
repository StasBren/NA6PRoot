// NA60+ 2100V P3: C++/ROOT direct-single-ion training + held-out basis pilot.
// TRAIN from archived weighted microscopic ion birth POSITIONS (not currents).
// TRUE signals: Garfield++ RKF, t_birth=0; NO electronics, NO physical strips.
// Two hypotheses on independent test ions:
//   A) raw inverse-distance weighted response templates (4 nearest SAME v sign)
//   B) low-rank PCA basis with interpolated coefficients (ranks 4,8,16).
// The PCA "oracle" projects a TRUE held-out waveform and measures compression
// alone; it is NOT an available fast prediction and is never used as one.
// No microscopic avalanche simulation, no full-muon waveform invented here.
#include "Plane2100Archive.hh"
#include <TFile.h>
#include <TTree.h>
#include <TNamed.h>
#include <TH1D.h>
#include <TMatrixDSym.h>
#include <TMatrixDSymEigen.h>
#include <TMatrixD.h>
#include <TVectorD.h>
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace {
using Vec=std::vector<double>;
using Birth=NA60Plane::Archive::Birth;
using Clock=std::chrono::steady_clock;
struct Point {Birth b;Vec wave;Vec coeff;};
struct Metrics {
  double relL2=0, maxCumErrorE=0, finalPlusErrorE=0,
         finalMinusErrorE=0, early1usErrorE=0;
};
std::string Arg(int n,char**v,const std::string& key){
  for(int i=1;i+1<n;++i)if(key==v[i])return v[i+1];return "";
}
int NArg(int n,char**v,const std::string& key,int def) {
  const auto s=Arg(n,v,key);return s.empty()?def:std::stoi(s);
}
double DArg(int n,char**v,const std::string& key,double def) {
  const auto s=Arg(n,v,key);return s.empty()?def:std::stod(s);
}
double Seconds(Clock::time_point a,Clock::time_point b){
  return std::chrono::duration<double>(b-a).count();
}
double Dot(const Vec&a,const Vec&b){
  NA60Plane::Require(a.size()==b.size(),"Dot: vector dimension mismatch");
  double s=0;for(std::size_t i=0;i<a.size();++i)s+=a[i]*b[i];return s;
}
double Distance(const Birth&a,const Birth&b){
  // Near-wire cylindrical geometry. Refuse interpolation ACROSS v=0:
  // ions would drift to opposite cathodes. r/angle are in chamber-local u,v.
  const double ra=std::hypot(a.u,a.v),rb=std::hypot(b.u,b.v);
  const double du=(a.u-b.u)/0.003; // 30 um
  const double dv=(a.v-b.v)/0.003;
  const double dlog=std::log(ra/rb)/0.4;
  return du*du+dv*dv+dlog*dlog;
}
std::vector<std::pair<int,double>> Neighbors(
   const Birth&b,const std::vector<Point>&train,int neighbors) {
  std::vector<std::pair<double,int>> distances;
  for(int i=0;i<(int)train.size();++i) {
    if((train[i].b.v>0)!=(b.v>0))continue;
    distances.emplace_back(Distance(b,train[i].b),i);
  }
  NA60Plane::Require((int)distances.size()>=neighbors,
    "insufficient same-cathode training positions");
  std::partial_sort(distances.begin(),distances.begin()+neighbors,
    distances.end());
  std::vector<std::pair<int,double>> out;
  if(distances[0].first<1.e-16) {
    out.emplace_back(distances[0].second,1.);
    return out;
  }
  double sum=0.;
  for(int j=0;j<neighbors;++j) {
    const double d=std::max(1.e-12,distances[j].first);
    const double weight=1./d;
    out.emplace_back(distances[j].second,weight);sum+=weight;
  }
  for(auto&v:out)v.second/=sum;
  return out;
}
Vec Simulate(NA60Plane::Transport&tr,const Birth&b,int bins){
  // All reference signals start at t=0 to isolate the spatial response K(u,v,t).
  // Original t_birth is preserved in archive but not needed in this kernel test.
  tr.CheckStart(b.u,b.v,b.w,0.);
  const auto d=tr.Drift(b.u,b.v,b.w,0.);
  (void)d;
  Vec wave(2*bins,0.);
  for(int t=0;t<bins;++t){
    wave[t]=tr.Current(0,t);
    wave[bins+t]=tr.Current(1,t);
  }
  return wave;
}
Metrics Evaluate(const Vec&pred,const Vec&ref,int bins,double dt){
  NA60Plane::Require(pred.size()==ref.size()&&pred.size()==2u*(unsigned)bins,
    "waveform shape mismatch");
  Metrics m;
  double sse=0,den=0,cumP=0,cumR=0;
  constexpr double q=NA60Plane::eFc;
  int early=std::min(bins,std::max(1,(int)std::ceil(1000./dt)));
  for(int c=0;c<2;++c) {
    double dFinal=0.;
    for(int t=0;t<bins;++t) {
      const int k=c*bins+t;
      const double d=pred[k]-ref[k];
      sse+=d*d;den+=ref[k]*ref[k];
      dFinal+=dt*d;
      cumP+=0.; // cumulative charge is computed separately per electrode
    }
    if(c==0)m.finalPlusErrorE=std::abs(dFinal)/q;
    else m.finalMinusErrorE=std::abs(dFinal)/q;
    double cumulative=0.;
    for(int t=0;t<bins;++t) {
      cumulative+=dt*(pred[c*bins+t]-ref[c*bins+t]);
      m.maxCumErrorE=std::max(m.maxCumErrorE,std::abs(cumulative)/q);
      if(t==early-1)m.early1usErrorE=std::max(m.early1usErrorE,
                                       std::abs(cumulative)/q);
    }
  }
  m.relL2=std::sqrt(sse/std::max(den,1.e-60));
  return m;
}
Vec RawIDW(const std::vector<Point>&train,
           const std::vector<std::pair<int,double>>&near,int len){
  Vec y(len,0.);
  for(const auto &[j,w]:near)for(int t=0;t<len;++t)y[t]+=w*train[j].wave[t];
  return y;
}
Vec Reconstruct(const Vec&mean,const std::vector<Vec>&basis,
                const Vec&coeff,int rank){
  Vec y=mean;
  for(int p=0;p<rank;++p)
    for(std::size_t k=0;k<y.size();++k)y[k]+=coeff[p]*basis[p][k];
  return y;
}
Vec InterpolatedCoeffs(const std::vector<Point>&train,
                  const std::vector<std::pair<int,double>>&near,int rank){
  Vec coeff(rank,0.);
  for(auto [j,w]:near)for(int p=0;p<rank;++p)coeff[p]+=w*train[j].coeff[p];
  return coeff;
}
Vec OracleCoeff(const Vec&truth,const Vec&mean,
                const std::vector<Vec>&basis,int rank) {
  Vec x=truth;
  for(std::size_t t=0;t<x.size();++t)x[t]-=mean[t];
  Vec c(rank,0.);
  for(int j=0;j<rank;++j)c[j]=Dot(x,basis[j]);
  return c;
}
void PrintSummary(const std::string&label,const std::vector<Metrics>&data,
                  std::ostream&out) {
  NA60Plane::Require(!data.empty(),"no metrics to summarize");
  std::vector<double> errors;
  double sum=0,maxCum=0,maxFinal=0,maxEarly=0;
  for(const auto&m:data){
    errors.push_back(m.relL2);sum+=m.relL2;
    maxCum=std::max(maxCum,m.maxCumErrorE);
    maxFinal=std::max(maxFinal,std::max(m.finalPlusErrorE,m.finalMinusErrorE));
    maxEarly=std::max(maxEarly,m.early1usErrorE);
  }
  std::sort(errors.begin(),errors.end());
  auto percentile=[&](double frac){
    return errors[std::min(errors.size()-1,
                          (std::size_t)std::ceil(frac*errors.size())-1)];
  };
  out<<label<<","<<data.size()<<","<<sum/data.size()<<","
     <<errors[errors.size()/2]<<","<<percentile(0.9)<<","
     <<errors.back()<<","<<maxCum<<","<<maxFinal<<","<<maxEarly<<"\n";
}
} // namespace
int main(int argc,char**argv){
 try{
  const auto file=Arg(argc,argv,"--ion-file");
  const auto prefix=Arg(argc,argv,"--out-prefix");
  const auto mobility=Arg(argc,argv,"--mobility");
  const int nTrain=NArg(argc,argv,"--train",96);
  const int nTest=NArg(argc,argv,"--test",48);
  const int seed=NArg(argc,argv,"--seed",210018);
  const int neighbors=NArg(argc,argv,"--neighbors",4);
  const double dt=DArg(argc,argv,"--dt-ns",100.);
  const double windowUs=DArg(argc,argv,"--window-us",100.);
  const double stepMm=DArg(argc,argv,"--max-step-mm",0.05);
  const int bins=(int)std::ceil(windowUs*1000./dt);
  NA60Plane::Require(!file.empty()&&!prefix.empty()&&!mobility.empty(),
   "usage: --ion-file WEIGHTED.root --mobility FILE --out-prefix PREFIX");
  NA60Plane::Require(nTrain>=32&&nTrain<=256&&nTest>=16&&nTest<=128,
   "pilot uses train=32..256 and test=16..128");
  NA60Plane::Require(neighbors>=1&&neighbors<=8&&neighbors<nTrain/4,
   "bad neighbor count");
  NA60Plane::Require(dt>=20&&windowUs>=50&&windowUs<=150&&
   bins>=500&&bins<=3000,"pilot time grid must be 500..3000 bins, dt>=20 ns");
  NA60Plane::Require(stepMm>0&&stepMm<=0.1,"bad RKF max step");
  for(const char* ext:{"_metrics.csv","_summary.csv","_basis.root"}) {
   NA60Plane::Require(!std::filesystem::exists(prefix+ext),
     "refusing to overwrite existing pilot "+prefix+ext);
  }

  const auto start=Clock::now();
  NA60Plane::Archive input(file);
  NA60Plane::Require(input.mode=="weighted","train on verified weighted ion archive");
  NA60Plane::Require(input.Entries()>=nTrain+nTest,
    "not enough archived ion birth positions");
  std::vector<Birth> source;
  source.reserve(input.Entries());
  for(long long i=0;i<input.Entries();++i)source.push_back(input.Get(i));
  input.CheckComplete();
  std::vector<int> positive,negative;
  for(int i=0;i<(int)source.size();++i) {
    if(source[i].v>0)positive.push_back(i);
    else if(source[i].v<0)negative.push_back(i);
  }
  NA60Plane::Require((int)positive.size()>=nTrain/2+nTest/2+2 &&
                     (int)negative.size()>=nTrain/2+nTest/2+2,
    "not enough births from both cathode hemispheres");
  std::mt19937_64 rng((unsigned long long)seed);
  std::shuffle(positive.begin(),positive.end(),rng);
  std::shuffle(negative.begin(),negative.end(),rng);
  std::vector<int> selectedTrain,selectedTest;
  for(int i=0;i<nTrain;++i) {
    auto &v=(i%2==0?positive:negative);
    selectedTrain.push_back(v[i/2]);
  }
  for(int i=0;i<nTest;++i) {
    auto &v=(i%2==0?positive:negative);
    selectedTest.push_back(v[(nTrain+1)/2+i/2]);
  }
  // All draws are without replacement, and test identities excluded from train.
  std::vector<int> all=selectedTrain;all.insert(all.end(),
     selectedTest.begin(),selectedTest.end());
  std::sort(all.begin(),all.end());
  NA60Plane::Require(std::adjacent_find(all.begin(),all.end())==all.end(),
    "train/test ion positions reused");
  const auto sampled=Clock::now();
  NA60Plane::Transport tr(mobility,dt,bins,stepMm);
  const int len=2*bins;
  std::vector<Point> train(nTrain);
  for(int i=0;i<nTrain;++i) {
    train[i].b=source[selectedTrain[i]];
    train[i].wave=Simulate(tr,train[i].b,bins);
    if((i+1)%32==0)std::cout<<"training RKF "<<i+1<<"/"<<nTrain<<"\n";
  }
  const auto trained=Clock::now();

  // A centered Gram-matrix PCA: the only large dense eigensolve is
  // nTrain x nTrain; never a 2*nBins square matrix.
  Vec mean(len,0.);
  for(const auto&row:train)for(int t=0;t<len;++t)mean[t]+=row.wave[t]/nTrain;
  std::vector<Vec> centered(nTrain,Vec(len,0.));
  for(int i=0;i<nTrain;++i)
   for(int t=0;t<len;++t)centered[i][t]=train[i].wave[t]-mean[t];
  TMatrixDSym gram(nTrain);
  for(int i=0;i<nTrain;++i)for(int j=i;j<nTrain;++j) {
    double dot=Dot(centered[i],centered[j]);
    gram(i,j)=dot;gram(j,i)=dot;
  }
  TMatrixDSymEigen eig(gram);
  TVectorD eigen=eig.GetEigenValues();
  TMatrixD vectors=eig.GetEigenVectors();
  std::vector<int> order(nTrain);
  std::iota(order.begin(),order.end(),0);
  std::sort(order.begin(),order.end(),[&](int a,int b){
    return eigen[a]>eigen[b];
  });
  const int maxRank=std::min(16,nTrain-2);
  std::vector<Vec> basis;
  std::vector<double> lambda;
  const double maxEval=std::max(0.,eigen[order[0]]);
  for(int p=0;p<maxRank;++p){
    int col=order[p];double ev=eigen[col];
    NA60Plane::Require(ev>std::max(maxEval*1.e-12,1.e-30),
      "insufficient meaningful PCA rank, try smaller target rank");
    Vec v(len,0.);
    for(int i=0;i<nTrain;++i)
      for(int t=0;t<len;++t)
        v[t]+=vectors(i,col)*centered[i][t]/std::sqrt(ev);
    // ROOT Gram-matrix eigenbasis is orthonormal in ideal arithmetic.
    double norm=std::sqrt(Dot(v,v));
    NA60Plane::Require(norm>0.99&&norm<1.01,"PCA eigenvector norm invalid");
    for(auto&x:v)x/=norm;
    basis.push_back(std::move(v));lambda.push_back(ev);
  }
  for(auto&trn:train){
    trn.coeff=OracleCoeff(trn.wave,mean,basis,maxRank);
  }
  const auto decomposed=Clock::now();

  std::vector<Point> test(nTest);
  for(int i=0;i<nTest;++i){
    test[i].b=source[selectedTest[i]];
    test[i].wave=Simulate(tr,test[i].b,bins);
    if((i+1)%16==0)std::cout<<"held-out RKF "<<i+1<<"/"<<nTest<<"\n";
  }
  const auto tested=Clock::now();

  std::ofstream row(prefix+"_metrics.csv");
  NA60Plane::Require(bool(row),"cannot create per-ion metrics CSV");
  row<<std::setprecision(12)
     <<"ion_source_row,primary,u_cm,v_cm,w_cm,birth_t_ns,model,rank,"
       "relative_L2,max_cumulative_error_e,final_plus_error_e,"
       "final_minus_error_e,early_1us_error_e\n";
  std::vector<Metrics> metricRaw, metricBasis4,metricBasis8,metricBasis16,
                       oracle4,oracle8,oracle16;
  for(int i=0;i<nTest;++i){
    const auto&b=test[i].b;
    const auto near=Neighbors(b,train,neighbors);
    auto evaluate=[&](const std::string& name,int rank,const Vec&wave,
                      std::vector<Metrics>&list){
      const auto m=Evaluate(wave,test[i].wave,bins,dt);
      list.push_back(m);
      row<<selectedTest[i]<<","<<b.primary<<","<<b.u<<","<<b.v<<","
         <<b.w<<","<<b.t<<","<<name<<","<<rank<<","
         <<m.relL2<<","<<m.maxCumErrorE<<","<<m.finalPlusErrorE<<","
         <<m.finalMinusErrorE<<","<<m.early1usErrorE<<"\n";
    };
    evaluate("raw_spatial_IDW",0,RawIDW(train,near,len),metricRaw);
    for(int rank:{4,8,16}){
      auto coeff=InterpolatedCoeffs(train,near,rank);
      auto interp=Reconstruct(mean,basis,coeff,rank);
      auto oracle=Reconstruct(mean,basis,
                     OracleCoeff(test[i].wave,mean,basis,rank),rank);
      if(rank==4) {
        evaluate("PCA_spatial_IDW",rank,interp,metricBasis4);
        evaluate("PCA_ORACLE_not_fast",rank,oracle,oracle4);
      } else if(rank==8) {
        evaluate("PCA_spatial_IDW",rank,interp,metricBasis8);
        evaluate("PCA_ORACLE_not_fast",rank,oracle,oracle8);
      } else {
        evaluate("PCA_spatial_IDW",rank,interp,metricBasis16);
        evaluate("PCA_ORACLE_not_fast",rank,oracle,oracle16);
      }
    }
  }
  row.close();
  const auto predicted=Clock::now();
  std::ofstream summary(prefix+"_summary.csv");
  NA60Plane::Require(bool(summary),"cannot create pilot summary CSV");
  summary<<std::setprecision(12)
    <<"model,n_test,mean_relative_L2,median_relative_L2,p90_relative_L2,"
      "max_relative_L2,max_cumulative_error_e,max_final_charge_error_e,"
      "max_early_1us_error_e\n";
  PrintSummary("raw_spatial_IDW",metricRaw,summary);
  PrintSummary("PCA_rank4_spatial_IDW",metricBasis4,summary);
  PrintSummary("PCA_rank8_spatial_IDW",metricBasis8,summary);
  PrintSummary("PCA_rank16_spatial_IDW",metricBasis16,summary);
  PrintSummary("PCA_rank4_ORACLE_not_fast",oracle4,summary);
  PrintSummary("PCA_rank8_ORACLE_not_fast",oracle8,summary);
  PrintSummary("PCA_rank16_ORACLE_not_fast",oracle16,summary);
  summary.close();

  TFile root((prefix+"_basis.root").c_str(),"RECREATE",
     "NA60+ Simulation: prototype 2100-V whole-plane PCA, no stereo strips",4);
  NA60Plane::Require(!root.IsZombie(),"ROOT basis output error");
  TNamed limits("limitations",
    "Experimental single-ion K(u,v,t), mean-drift RKF and Ar+ in Ar provisional "
    "mobility for Ar/CO2; independent held-out sampled births; "
    "t_birth zeroed for kernel study. NO full-muon waveform, electronics, "
    "physical stereo strips, diffusion, space charge or calibrated ion species.");
  limits.Write();
  TNamed provenance("input_ion_archive",file.c_str());provenance.Write();
  TH1D variances("pca_variance","PCA Gram eigenvalues;rank;lambda",maxRank,0.5,maxRank+0.5);
  for(int r=0;r<maxRank;++r)variances.SetBinContent(r+1,lambda[r]);
  variances.Write();
  TH1D averagePlus("average_plus","Training mean;time since ion birth (ns);I (fC/ns)",
                   bins,0.,bins*dt);
  TH1D averageMinus("average_minus","Training mean;time since ion birth (ns);I (fC/ns)",
                    bins,0.,bins*dt);
  for(int k=0;k<bins;++k){
    averagePlus.SetBinContent(k+1,mean[k]);
    averageMinus.SetBinContent(k+1,mean[bins+k]);
  }
  averagePlus.Write();averageMinus.Write();
  for(int r=0;r<maxRank;++r)for(int c=0;c<2;++c) {
    const std::string name="basis_"+std::to_string(r+1)+(c==0?"_plus":"_minus");
    TH1D h(name.c_str(),"NA60+ Simulation;time since birth (ns);unit PCA basis",
           bins,0.,bins*dt);
    for(int k=0;k<bins;++k)h.SetBinContent(k+1,basis[r][c*bins+k]);
    h.Write();
  }
  // Save full training bank with coordinates and coefficients, to reuse.
  TTree bank("training_bank","Held-out separated ion-start kernel templates");
  double u=0,v=0,w=0,t=0;int originalRow=0,primary=0;
  Vec savedWave,savedCoeff;
  bank.Branch("original_ion_row",&originalRow);
  bank.Branch("original_primary",&primary);
  bank.Branch("u_cm",&u);bank.Branch("v_cm",&v);
  bank.Branch("w_cm",&w);bank.Branch("original_birth_t_ns",&t);
  bank.Branch("waveform_plus_then_minus",&savedWave);
  bank.Branch("coefficients",&savedCoeff);
  for(int i=0;i<nTrain;++i) {
    const auto&b=train[i].b;originalRow=selectedTrain[i];
    primary=b.primary;u=b.u;v=b.v;w=b.w;t=b.t;
    savedWave=train[i].wave;savedCoeff=train[i].coeff;
    bank.Fill();
  }
  bank.Write();root.Close();
  std::cout<<std::setprecision(6)<<std::fixed
    <<"PASS experimental compact kernel pilot, source="<<input.id<<"\n"
    <<"training_ions="<<nTrain<<" held_out_ions="<<nTest
    <<" waveform_bins="<<bins<<" dt_ns="<<dt<<"\n"
    <<"times_sec read="<<Seconds(start,sampled)
    <<" train_RKF="<<Seconds(sampled,trained)
    <<" PCA="<<Seconds(trained,decomposed)
    <<" test_RKF="<<Seconds(decomposed,tested)
    <<" fast_predictions="<<Seconds(tested,predicted)<<"\n"
    <<"IMPORTANT: PASS means computations completed; accuracy is in CSV, "
      "not automatically accepted.\n";
  std::ifstream view(prefix+"_summary.csv");std::string line;
  while(std::getline(view,line))std::cout<<line<<"\n";
  return 0;
 }catch(const std::exception& e) {
  std::cerr<<"COMPACT PLANE KERNEL PILOT FAILURE: "<<e.what()<<"\n";
  return 2;
 }
}
