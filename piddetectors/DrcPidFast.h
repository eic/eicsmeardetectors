// Fast reconstruction for EIC Barrel DIRC
// original author: r.dzhygadlo at gsi.de

#ifndef DrcPidFast_h
#define DrcPidFast_h 1

#include "TFile.h"
#include "TH2F.h"
#include "TVector3.h"
#include "TRandom.h"
#include <iostream>

// probability - normalized to 1 probability for e,mu,pi,k,p
// sigma - deviation of the determined Cherenkov angle from expected in terms of Cherenkov track resolution
// cangle - Cherenkov angle
// cctr -  combined Cherenkov track resolution
struct DrcPidInfo {
  double probability[5];
  double sigma[5];
  double cangle;
  double cctr;
};

class DrcPidFast{
  
 public:
  DrcPidFast();
  ~DrcPidFast(){}

  // read Cherenkov track resolution map from a file
  void ReadMap(TString name);
  
  // pdg - Particle Data Group code of the particle
  // mom - 3-momentum of the particle [GeV/c]
  // track_err - error assosiated with track direction [mrad]
  DrcPidInfo GetInfo(int pdg,TVector3 mom, double track_err=0);

  // p - momentum of the particle [GeV/c]
  // theta - polar angle of the particle [deg]
  DrcPidInfo GetInfo(int pdg, double p, double theta, double track_err=0);
  TH2F *GetTrrMap(){ return fTrrMap; }
  
 private:
  int get_pid(int pdg);
  TH2F *fTrrMap;
  double fMass[5];
  TRandom  fRand;
};

#endif
