#include "barrelDirc.h"

using namespace std;

barrelDirc::barrelDirc(double trackResolution, double timePrecision, int qe, double eL, double eH)
{
  fTrackResolution = trackResolution;
  fTimePrecision = timePrecision;
  etaLow  = eL;
  etaHigh = eH;
  fQe=qe;

  int id[]={27,22};
  char nameit[500];
  sprintf(nameit,"Barrel DIRC TR=%1.1f [mrad] dT=%1.1f ns QE = %d %%",fTrackResolution,fTimePrecision,id[fQe]);
  myName = nameit;

  if(fQe==1) pid.ReadMap("ctr_map_p2_0.95.root"); // map for 22% quantum efficiency
  else pid.ReadMap("ctr_map_p1_0.95.root"); // map for 27% quantum efficiency
}

double barrelDirc::numSigma(double eta, double p, PID::type PID)
{
  if(valid(eta,p)){
    double theta = 2.0*atan(exp(-eta))*TMath::RadToDeg();

    if(PID == pi_k){
      DrcPidInfo info_pi = pid.GetInfo(211,p,theta,fTrackResolution);
      DrcPidInfo info_k = pid.GetInfo(321,p,theta,fTrackResolution);
      return fabs(info_pi.cangle-info_k.cangle)/info_pi.cctr;
    }else if(PID == k_p){
      DrcPidInfo info_k = pid.GetInfo(321,p,theta,fTrackResolution);
      DrcPidInfo info_p = pid.GetInfo(2212,p,theta,fTrackResolution);
      return fabs(info_p.cangle-info_k.cangle)/info_p.cctr;
    }else{
      cout << "barrelDirc.C:  Unrecognized PID type requested." <<endl;
    }
  }else{
    cout << "barrelDirc.C:  Invalid (eta,p) for this detector." <<endl;
  }
  return 0;
}

double barrelDirc::numSigma (double eta, double p, const PID::Species truth, const PID::Species reference)
{
  if(valid(eta,p)){
    double theta = 2.0*atan(exp(-eta))*TMath::RadToDeg();

    // sigh. These enums are silly.
    int truthpdg = 0;
    switch ( truth ){
    case kElectron : truthpdg=11; break;
    case kPion     : truthpdg=211; break;
    case kProton   : truthpdg=2212; break;
    case kKaon     : truthpdg=321; break;
    case kMuon     : truthpdg=13; break;
    }
    int refpdg = 0;
    switch ( reference ){
    case kElectron : refpdg=11; break;
    case kPion     : refpdg=211; break;
    case kProton   : refpdg=2212; break;
    case kKaon     : refpdg=321; break;
    case kMuon     : refpdg=13; break;
    }


    DrcPidInfo info_ref = pid.GetInfo(refpdg,p,theta,fTrackResolution);
    DrcPidInfo info_truth = pid.GetInfo(truthpdg,p,theta,fTrackResolution);
    // cout << info_truth.cangle << "  " << info_ref.cangle << "  " << info_truth.cctr << endl;
    return fabs(info_truth.cangle-info_ref.cangle)/info_truth.cctr;

  }else{
    cout << "barrelDirc.C:  Invalid (eta,p) for this detector." <<endl;
  }
  return 0;

}

double barrelDirc::maxP(double eta, double numSigma, PID::type PID)
{
  double theta = 2.0*atan(exp(-eta))*TMath::RadToDeg();  
  if (valid(eta,1.0)){
    if(PID == pi_k){
      DrcPidInfo info_pi,info_k;
      for(double p=15; p>0.4; p-=0.01){
	info_pi = pid.GetInfo(211,p,theta,fTrackResolution);
	info_k = pid.GetInfo(321,p,theta,fTrackResolution);
	double sep = fabs(info_pi.cangle-info_k.cangle)/info_pi.cctr;
	if(sep>numSigma) return p;
      }
    }else if(PID == k_p){
      DrcPidInfo info_k,info_p;
      for(double p=15; p>0.4; p-=0.01){
	info_k = pid.GetInfo(321,p,theta,fTrackResolution);
	info_p = pid.GetInfo(2212,p,theta,fTrackResolution);
	double sep = fabs(info_p.cangle-info_k.cangle)/info_p.cctr;
	if(sep>numSigma) return p;
      }
    }else{
      cout << "barrelDirc.C:  Unrecognized PID type requested." <<endl;
    }  
  }else{
    cout << "barrelDirc.C:  Invalid (eta) for this detector." <<endl;
  }
  return 0;
}

void barrelDirc::description()
{
  //  Here one should describe what this detector is and what assumptions are 
  //  made in calculating the performance specifications.

  int id[]={27,22};  
  cout << "My name is \"" << myName << "\" and I am described as follows:" <<endl;
  cout << "    Eta coverage =  [" << etaLow << "," << etaHigh<<"]"<<endl;
  cout << "    Assumed time precision = " << fTimePrecision << " ns" <<endl;
  cout << "    Assumed track resolution = "<< fTrackResolution << " mrad" <<endl;
  cout << "    Assumed quantum efficiency of the MCP-PMT = "<< id[fQe] << "%" <<endl;
  cout << endl;
}
