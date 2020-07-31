#ifndef EICSMEARDETECTORS_HH
#define EICSMEARDETECTORS_HH

#include <string>
#include <iostream>
#include <cctype>

#include "eicsmear/smear/Detector.h"

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


namespace Smear{
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
  Smear::Detector BuildByName (std::string dname){
    // transform to upper case
    for (auto & c: dname) c = toupper(static_cast<unsigned char>(c));

    // -- Online, OFFICIAL, matrix from https://physdiv.jlab.org/DetectorMatrix
    if ( dname == "MATRIXDETECTOR_0_1"  ||
	 dname == "MATRIX_0_1" ||
	 dname == "MATRIX" ) return BuildMatrixDetector_0_1();

    // -- Handbook matrix from http://www.eicug.org/web/sites/default/files/EIC_HANDBOOK_v1.2.pdf
    if ( dname == "HANDBOOK_1_2" ||
	 dname == "HANDBOOK" ) return BuildHandBook_1_2();

    // -- Perfect detection and PID in |eta|<15
    if ( dname == "PERFECTDETECTOR" ||
         dname == "PERFECT" ) return BuildPerfectDetector();

    // -- Inofficial detector scripts
    // ---- BeAST
    if ( dname == "BEAST_0_1" ||
	 dname == "BEAST" ) return BuildBeAST_0_1();
    // ---- Jleic
    if ( dname == "JLEIC_0_1" ||
         dname == "JLEIC" ) return BuildJLEIC_0_1();
    // -- Older legacy detector scripts; may require adding
    //    det.SetLegacyMode(true);
    // ---- BeAST
    if ( dname == "BEAST_0_0" ||
         dname == "BEAST" ) return  BuildBeAST_0_1();
    // ---- ZEUS
    if ( dname == "ZEUS_0_0" ||
         dname == "ZEUS" ) return  BuildZEUS_0_0();
    // eSTAR
    if ( dname == "ESTAR_0_0" ||
         dname == "ESTAR" ) return  BuildeSTAR_0_0();
    // STAR
    if ( dname == "STAR_0_0" ||
	 dname == "STAR" ) return  BuildSTAR_0_0();

    // Note that BuildePHENIX(bool) appears again below
    if ( dname == "EPHENIX_0_0" ||
	 dname == "EPHENIX" ) return  BuildePHENIX_0_0();

    std::cerr << "Detector sepcified as " << dname
	      << " not recognized or empty." << std::endl;
    throw;
    return Smear::Detector();
  };

  /** Overloaded version of  Smear::Detector BuildByName ( std::string dname )
      for detectors with a bool parameter
  */
  Smear::Detector BuildByName ( std::string dname, const bool b){
    // transform to upper case
    for (auto & c: dname) c = toupper(static_cast<unsigned char>(c));

    if ( dname == "EPHENIX_0_0" ||
	 dname == "EPHENIX" ) return  BuildePHENIX_0_0( b ); // b is multipleScattering

    std::cerr << "Detector sepcified as " << dname
	      << " not recognized or empty." << std::endl;
    throw;
    return Smear::Detector();
  };

}


#endif //EICSMEARDETECTORS_HH
