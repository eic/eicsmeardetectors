#ifndef INCLUDE_DIRCBARREL_H
#define INCLUDE_DIRCBARREL_H

#include "eicsmear/smear/NumSigmaPid.h"
#include "barrelDirc.h"

namespace Smear{
  class DircBarrelSmearer : public Smear::NumSigmaPid {
  public :
    
    /** standard ctor
	trackResolution is pointing resolution in mrad
    */
    DircBarrelSmearer( double trackResolution=0.5, double timePrecision=0.1, int qe=0, double etaLow=-1.5, double etaHigh=1.5){
      ThePidObject = std::make_shared<barrelDirc>(trackResolution, timePrecision, qe, etaLow, etaHigh);
    }; // ctor
  };// class
};// namespace
#endif //  INCLUDE_DIRCBARREL_H
