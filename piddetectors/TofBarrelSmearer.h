#ifndef INCLUDE_TOFBARREL_H
#define INCLUDE_TOFBARREL_H

#include "eicsmear/smear/NumSigmaPid.h"
#include "tofBarrel.h"

namespace Smear{
  class TofBarrelSmearer : public Smear::NumSigmaPid {
  public :
    
    /** standard ctor
     */
    TofBarrelSmearer( double radius=100, double etaLow=-1.0, double etaHigh=1.0, double sigmaT=10 ){
      ThePidObject = std::make_shared<tofBarrel>(radius, etaLow, etaHigh,sigmaT);
    };
  };
};

#endif //  INCLUDE_TOFBARREL_H
