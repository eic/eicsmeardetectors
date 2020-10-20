#include "DrcPidFast.h"

// using namespace std;

DrcPidFast::DrcPidFast(){

  fMass[0]=0.000511;
  fMass[1]=0.105658;
  fMass[2]=0.139570;
  fMass[3]=0.49368;
  fMass[4]=0.938272;
  
  // read Cherenkov track resolution map
  ReadMap("ctr_map_p1_0.95.root");
}


void  DrcPidFast::ReadMap(TString name){
  TFile* file = TFile::Open(name);
  fTrrMap = new TH2F();
  file->GetObject("htrr", fTrrMap);  
}

DrcPidInfo DrcPidFast::GetInfo(int pdg,TVector3 mom, double track_err){
  double p = mom.Mag();
  double theta = mom.Theta()*TMath::RadToDeg();
  return GetInfo(pdg,p,theta,track_err);
}

DrcPidInfo DrcPidFast::GetInfo(int pdg, double p, double theta, double track_err){

  const int max = 5;
  DrcPidInfo info;
  int pid = get_pid(pdg);
  
  // set default values
  for(int i=0; i<max; i++){
    info.probability[i]=0.25;
    info.sigma[i]=100;
  }
  info.cangle=0;
  info.cctr=0;
  
  // check range
  if(theta<25 || theta>153) return info;
      
  int bin = fTrrMap->FindBin(theta,(p>10)? 10:p); // ctr map is till 10 GeV/c
  double ctr = fTrrMap->GetBinContent(bin);  // Cherenkov track resolution [mrad]
  double cctr = sqrt(ctr*ctr+track_err*track_err)*0.001; // combined Cherenkov track resolution[rad]
  
  // 1.46907 - fused silica
  double true_cangle = acos(sqrt(p*p + fMass[pid]*fMass[pid])/p/1.46907);  
  true_cangle += fRand.Gaus(0,cctr);

  // return default values if momentum below Cherenkov threshold (true_cangle is NaN)
  if(true_cangle != true_cangle) return info;
  
  double cangle,sum=0,fsum=0;
  double delta[max]={0}, probability[max]={0};

  for(int i=0; i<max; i++){
    cangle = acos(sqrt(p*p + fMass[i]*fMass[i])/p/1.46907);
    if(cangle != cangle) continue;
    delta[i] = fabs(cangle-true_cangle);
    sum += delta[i];    
    info.sigma[i]=(cangle-true_cangle)/cctr;
    if(i==pid) info.cangle = cangle;
  }
  // normalization
  for(int i=0; i<max; i++){
    if(delta[i]>0) info.probability[i] = sum/delta[i];
    fsum += info.probability[i];
  }
  for(int i=0; i<max; i++) info.probability[i] /= fsum;
  info.cctr = cctr;
  
  return info;
}

int DrcPidFast::get_pid(int pdg){
  int pid=0;
  if(pdg==11)   pid=0; //e
  if(pdg==13)   pid=1; //mu
  if(pdg==211)  pid=2; //pi
  if(pdg==321)  pid=3; //K
  if(pdg==2212) pid=4; //p
  return pid;
}
