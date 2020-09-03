#ifndef EICSMEARDETECTORS_HH
#define EICSMEARDETECTORS_HH

#include <string>
#include <iostream>
#include <cctype>

#include "eicsmear/smear/Detector.h"
// #include "eicsmear/smear/NumSigmaPid.h"
// #include "piddetectors/TofBarrelSmearer.h"
// #include "piddetectors/tofBarrel.h"

Smear::Detector BuildMatrixDetector_0_1();
Smear::Detector BuildHandBook_1_2();
Smear::Detector BuildPerfectDetector();
Smear::Detector BuildJLEIC_0_1();
Smear::Detector BuildBeAST_0_1();
Smear::Detector BuildBeAST_0_0();
Smear::Detector BuildSTAR_0_0();
Smear::Detector BuildZEUS_0_0();
Smear::Detector BuildeSTAR_0_0();
Smear::Detector BuildePHENIX_0_0(bool multipleScattering=true);

// experimental
Smear::Detector BuildWithTof();
Smear::Detector BuildMatrixDetector_0_1_FF( const int beam_mom_nn=275 );
Smear::Detector BuildMatrixDetector_0_1_TOF();
Smear::Detector BuildMatrixDetector_0_1_JHFCalo(const int radius=4000 );


/** For convenience.
    Not case-sensitive.Should be all upper case.
    Use like this:
    Smear::Detector detector = BuildByName("MATRIX");
    
    Provides multiple aliases
*/
// We could probably pull some tricks with variadic arguments,
// https://en.cppreference.com/w/cpp/utility/variadic
// But it's probably safer and more readable to
// overload below for scripts that need parameters
// Note that if you allow a default parameter, the detector needs
// to appear here here as well
// Notes:
// - If we put it in the Smear namespace, for some reason it doesn't get picked up by the autoloader
// - Tab completion for plain functions isn't supported by root (modules are the future, but that doesn't help)

Smear::Detector BuildByName (std::string dname);

/** Overloaded version of  Smear::Detector BuildByName ( std::string dname )
    for detectors with a bool parameter
*/
Smear::Detector BuildByName ( std::string dname, const bool b);


/** Overloaded version of  Smear::Detector BuildByName ( std::string dname )
    for detectors with an int parameter
*/
Smear::Detector BuildByName ( std::string dname, const int i);
  
#endif //EICSMEARDETECTORS_HH
