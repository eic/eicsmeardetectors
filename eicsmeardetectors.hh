#ifndef EICSMEARDETECTORS_HH
#define EICSMEARDETECTORS_HH

#include <map>
#include <string>

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

// For convenience.
// Should be all upper case.
// Use like this:
// Smear::Detector detector = BuildyName["MATRIX"]();
// Or safer,
// auto detfunc = BuildyName["MATRIX"];
// Smear::Detector detector;
// if (detfunc) detector = detfunc();


// You can transform a general string like this:
// string detstring = "BeAsT";
// for (auto & c: detstring) c = toupper(c);

std::map< std::string, Smear::Detector (*)()> BuildyName = {
  // -- Online, OFFICIAL, matrix from https://physdiv.jlab.org/DetectorMatrix
  { "MATRIXDETECTOR_0_1", BuildMatrixDetector_0_1 },
  { "MATRIXDETECTOR", BuildMatrixDetector_0_1 },
  { "MATRIX_0_1", BuildMatrixDetector_0_1 },
  { "MATRIX", BuildMatrixDetector_0_1 },
  // -- Handbook matrix from http://www.eicug.org/web/sites/default/files/EIC_HANDBOOK_v1.2.pdf
  { "HANDBOOK_1_2", BuildHandBook_1_2 },
  { "HANDBOOK", BuildHandBook_1_2 },
  // -- Perfect detection and PID in |eta|<15
  { "PERFECTDETECTOR", BuildPerfectDetector },
  { "PERFECT", BuildPerfectDetector },
  // -- Inofficial detector scripts
  // ---- BeAST
  {"BEAST_0_1", BuildBeAST_0_1},
  {"BEAST", BuildBeAST_0_1},
  // ---- Jleic
  {"JLEIC_0_1", BuildJLEIC_0_1},
  {"JLEIC", BuildJLEIC_0_1},
  // -- Older legacy detector scripts; may require adding
  //    det.SetLegacyMode(true);
  // ---- BeAST
  {"BEAST_0_0", BuildBeAST_0_0},
  {"BEAST", BuildBeAST_0_1},
  // ---- ZEUS
  {"ZEUS_0_0", BuildZEUS_0_0},
  {"ZEUS", BuildZEUS_0_0},
  // eSTAR
  {"ESTAR_0_0", BuildeSTAR_0_0},
  {"ESTAR", BuildeSTAR_0_0},
  // STAR
  {"STAR_0_0", BuildSTAR_0_0},
  {"STAR", BuildSTAR_0_0}
  // Note that BuildePHENIX(bool) has a different signature and can't be delivered the same way
};


#endif //EICSMEARDETECTORS_HH

