#ifndef __TOFBARREL_H__
#define __TOFBARREL_H__
	
//
//  Hello tofBarrel Fans:
//
//  This is an example class that inherits from the PID base class.
//  We shall implement all the necessary functions as required by the
//  base class.
//
//  The UNIQUE identifiers for this class are radius, eta extent, and 
//  time resolution. 
//
//  Note that in keeping with any well written base class, we shall
//  give default values to all the constructor arguments.
//
//  This routine assumes units of cm for distances and picoSeconds for time.
//
	
#include "eicsmear/smear/NumSigmaPid.h"

class tofBarrel: public PID
{
public:
  tofBarrel(double radius=100, double etaLow=-1.0, double etaHigh=1.0, double sigmaT=10);
  virtual ~tofBarrel() {}
	
  bool   valid   (double eta, double p) {return (eta>etaLow && eta<etaHigh);}
  double numSigma(double eta, double p,        PID::type PID);
  double maxP    (double eta, double numSigma, PID::type PID);
  double minP    (double eta, double numSigma, PID::type PID) {return 0;}
  std::string name    () {return myName;}
  void   description ();
		
protected:
  std::string myName;

  // utility function
  double tof(double L, double p, double m);

  // TOF wall parameters
  double radius;   // cm
  double etaLow;
  double etaHigh;
  double sigmaT;   // picosecond

  // Physical constants (should come from elsewhere!)
  double mPion;    // GeV/c^2
  double mKaon;    // GeV/c^2
  double mProton;  // GeV/c^2
  double c;        // cm/picosecond;
};
	
#endif /* __PID_H__ */
