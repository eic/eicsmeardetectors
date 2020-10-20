#ifndef __BARRELDIRC_H__
#define __BARRELDIRC_H__
	
//
//  wrapper class for DrcPidFast (fast PID for for the EIC Barrel DIRC)
//

#include "eicsmear/smear/NumSigmaPid.h"
#include "DrcPidFast.h"
// #include "DrcPidFast.cxx"
#include <string>

class barrelDirc: public PID
{
public:
  barrelDirc(double trackResolution=0.5, double timePrecision=0.1, int qe=0, double etaLow=-1.5, double etaHigh=1.5);
  virtual ~barrelDirc() {}
	
  bool   valid   (double eta, double p) {return (eta>etaLow && eta<etaHigh);}
  double numSigma(double eta, double p,        PID::type PID);
  double maxP    (double eta, double numSigma, PID::type PID);
  double minP    (double eta, double numSigma, PID::type PID) {return 0;}
  std::string name    () {return myName;}
  void   description ();

  double numSigma (double eta, double p, const PID::Species truth, const PID::Species reference);

    
protected:
  std::string myName;

  DrcPidFast pid;
  DrcPidInfo info;

  double fTrackResolution; // resolution of the traker [mrad]
  double fTimePrecision;   // time precision of the MCP-PMT [ns]
  int    fQe;              // id for Quantum efficiency of the MCP-PMT
  double etaLow;
  double etaHigh;
};
	
#endif
