#include "eicsmeardetectors.hh"

Smear::Detector BuildByName (std::string dname){
  // transform to upper case
  for (auto & c: dname) c = toupper(static_cast<unsigned char>(c));
  
  // -- Online, OFFICIAL, matrix from https://physdiv.jlab.org/DetectorMatrix
  // From November 21, 2020
  if ( dname == "MATRIXDETECTOR_0_2_B1_5T"  ||
       dname == "MATRIX_0_2_B1_5T" ||
       dname == "MATRIX02B15" ) return BuildMatrixDetector_0_2_B1_5T();

  if ( dname == "MATRIXDETECTOR_0_2_B3T"  ||
       dname == "MATRIX_0_2_B3T" ||
       dname == "MATRIX02B3" ) return BuildMatrixDetector_0_2_B3T();

  // -- Online, OFFICIAL, original matrix from https://physdiv.jlab.org/DetectorMatrix
  // From June 16 2020
  if ( dname == "MATRIXDETECTOR_0_1"  ||
       dname == "MATRIX_0_1" ||
       dname == "MATRIX01" ) return BuildMatrixDetector_0_1();
  
  // -- Handbook matrix from http://www.eicug.org/web/sites/default/files/EIC_HANDBOOK_v1.2.pdf
  if ( dname == "HANDBOOK_1_2" ||
       dname == "HANDBOOK" ) return BuildHandBook_1_2();
  
  // -- Perfect detection and PID in |eta|<15
  if ( dname == "PERFECTDETECTOR" ||
       dname == "PERFECT" ) return BuildPerfectDetector();

  // EXPERIMENTAL
  // if ( dname == "MATRIXTOF" ) return BuildMatrixDetector_0_1_TOF();
  
  if ( dname == "TRACKINGPREVIEW_0_2_B1_5T"  ||
       dname == "TRACKING_0_2_B1_5T" ||
       dname == "TRACKING02B15" ) return BuildTrackingPreview_0_2_B1_5T();

  if ( dname == "TRACKINGPREVIEW_0_2_B3T"  ||
       dname == "TRACKING_0_2_B3" ||
       dname == "TRACKING02B3" ) return BuildTrackingPreview_0_2_B3T();

    

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

Smear::Detector BuildByName ( std::string dname, const int i){
  // transform to upper case
  for (auto & c: dname) c = toupper(static_cast<unsigned char>(c));
  
  // -- UNOFFICIAL addition of far forward detectors to the matrix
  if ( dname == "MATRIXDETECTOR_0_1_FF"  ||
       dname == "MATRIX_0_1_FF" ||
       dname == "MATRIXFF" ) return BuildMatrixDetector_0_1_FF( i ); // i is beam_mom_nn
  
  std::cerr << "Detector sepcified as " << dname
	    << " not recognized or empty." << std::endl;
  throw;
  return Smear::Detector();
};



