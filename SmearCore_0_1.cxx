// Based on
// https://indico.bnl.gov/event/7913/contributions/41704/attachments/30564/47972/200924_Parametrizations.pdf
// https://indico.bnl.gov/event/11053/contributions/46969/attachments/33324/53540/core.pdf

// From March 29, 2021
// 
// By default, use 3T field
// Command line argument allows any other one, but be aware that
// it's purely a linear (affine; pol1) inter-/extrapolation
// between the parameters at 1.4 and 3T
// 
// Important Notes:
// - Where ranges are given, the more conservative number is chosen.
// - Without available specifications, angular resolution is assumed to be perfect.
// - Momentum and energy acceptance specifications are NOT implemented
// -- It is not a priori clear whether cuts should apply to truth or smeared value
// -- A "looper" may not reach the detector.
// -- On the other hand, a calorimetry deposit may fluctuate above threshold
// - PID for pi/K/p is implemented with perfect PID in the acceptance (and momentum range)
// - Electron PID is not implemented. Very high purity is important and not handled well with the above approach

#include "eicsmear/erhic/VirtualParticle.h"
#include "eicsmear/smear/Acceptance.h"
#include "eicsmear/smear/Device.h"
#include "eicsmear/smear/Detector.h"
#include "eicsmear/smear/Smearer.h"
#include "eicsmear/smear/ParticleMCS.h"
#include "eicsmear/smear/PerfectID.h"
#include <eicsmear/smear/Smear.h>
#include <eicsmear/erhic/ParticleMC.h>
#include "Math/Vector4D.h"

// using valarray instead of vector because entry-wise arithmetic is built-in
#include <valarray>

// declare static --> local to this file, won't clash with others
static double ThetaFromEta( const double eta );

// FIXME: Comments
static void AssembleCoreTracker ( Smear::Detector& det,
				  const std::valarray<double>& eta_min, const std::valarray<double>& eta_max,
				  const std::valarray<double>& A,   const std::valarray<double>& B );


/**
   Momentum resolution: sigma_p/p ~ A% * p + B% 
   Here: Inter-/Extrapolate A from two different B fields 
   Linear interpolation works well.
 */
static std::valarray<double> CalcA ( const double Bfield,
				  const double x1, const std::valarray<double> A1,
				  const double x2, const std::valarray<double> A2 );

/**
   Momentum resolution: sigma_p/p ~ A% * p + B% 
   Here: Inter-/Extrapolate B from two different B fields 
   Linear interpolation works well.
 */
static std::valarray<double> CalcB ( const double Bfield,
				  const double x1, const std::valarray<double> B1,
				  const double x2, const std::valarray<double> B2 );

Smear::Detector BuildCore_0_1( const double Bfield ) {
  gSystem->Load("libeicsmear");

  // Create the detector object to hold all devices
  Smear::Detector det;

  // The framework provides implementations of three kinematic calculation methods
  // from smeared values
  // NM - "Null method" uses the scattered lepton.
  // DA - Double Angle method
  // JB - Jacquet-Blondel method
  // Users should be mindful of the limitations and assumptions in
  // the calculation of smeared kinematics. Sophisticated calculations are best redone in user code
  // from the smeared particles directly.
  det.SetEventKinematicsCalculator("NM DA JB");

  // Perfect phi and theta for all particles
  // ---------------------------------------
  // total coverage for tracker, ecal, and hcal is -4.0 < eta < 4.0
  Smear::Acceptance::Zone AngleZoneCommon(ThetaFromEta ( 4.0 ),ThetaFromEta ( -4.0 ));
  Smear::Device SmearThetaCommon(Smear::kTheta, "0.0");
  SmearThetaCommon.Accept.AddZone(AngleZoneCommon);
  SmearThetaCommon.Accept.SetGenre(Smear::kAll);
  det.AddDevice(SmearThetaCommon);

  Smear::Device SmearPhiCommon(Smear::kPhi, "0.0");
  SmearPhiCommon.Accept.AddZone(AngleZoneCommon);
  SmearPhiCommon.Accept.SetGenre(Smear::kAll);
  det.AddDevice(SmearPhiCommon);

  // Tracking  (B = 3 T)
  // --------
  // Note: Smear::kCharged checks pdg charge, so includes muons (good)

  // |eta|
  const std::valarray<double> eta_min = {   0,   0.5,     1,   1.5,     2,   2.5,     3,   3.5 };
  const std::valarray<double> eta_max = { 0.5,     1,   1.5,     2,   2.5,     3,   3.5,     4 };
  // sigma_p/p ~ A% * p + B% 
  const std::valarray<double> A3T   = { 0.018, 0.016, 0.016, 0.012, 0.018, 0.039, 0.103, 0.295 };
  const std::valarray<double> B3T   = { 0.369, 0.428, 0.427, 0.462, 0.719, 1.336, 2.428, 4.552 };
  const std::valarray<double> A1p4T = { 0.038, 0.035, 0.035, 0.026, 0.041, 0.088, 0.217, 0.610 };
  const std::valarray<double> B1p4T = { 0.816, 0.898, 0.921, 0.997, 1.548, 2.830, 5.234, 9.797 };
  
  const std::valarray<double> A = CalcA ( Bfield, 1.4, A1p4T, 3.0, A3T );
  const std::valarray<double> B = CalcB ( Bfield, 1.4, B1p4T, 3.0, B3T );
  AssembleCoreTracker ( det, eta_min, eta_max, A, B );
  
  // PID
  // ---
  // Perfect pi/K/P PID everywhere for now
  // Make sure to not cover more than is covered by the other detectors.
  // No minimum momentum is used.
  // Accept charged hadrons (excludes e, mu)

  // p < 10 GeV
  Smear::Acceptance::Zone PidZone(ThetaFromEta(-4.0),ThetaFromEta(-4.0),
				  0., TMath::TwoPi(), // phi
				  0., TMath::Infinity(), // E
				  0., 10 ); // p
  Smear::PerfectID Pid;
  Pid.Accept.AddZone(PidZone);
  Pid.Accept.SetCharge(Smear::kCharged);
  Pid.Accept.SetGenre(Smear::kHadronic);
  det.AddDevice( Pid );

  // EM Calorimeters
  // ---------------
  // Note: Smear::kElectromagnetic == gamma + e. Does not include muons (good)
  // Calorimeter resolution usually given as sigma_E/E = A% / E + B%/Sqrt{E} + C%
  // EIC Smear needs absolute sigma: sigma_E = Sqrt{A*A + B*B*E + C*C*E*E}

  // E> 20 (50) MeV not used because it suppresses smeared-up particles

  // Inner Back
  // eta = -3.5 -- -2.0
  // 1%/E + 2.5%/sqrtE + 1%
  // A       B           C
  Smear::Acceptance::Zone EmcalInnerBackZone(ThetaFromEta ( -2.0 ),ThetaFromEta ( -3.5 ));
  Smear::Device EmcalInnerBack(Smear::kE, "sqrt( pow ( 0.01,2 ) + pow( 0.025,2)*E + pow ( 0.01*E,2 ) )");
  EmcalInnerBack.Accept.AddZone(EmcalInnerBackZone);
  EmcalInnerBack.Accept.SetGenre(Smear::kElectromagnetic);
  det.AddDevice(EmcalInnerBack);

  // Outer Back / Barrel
  // eta = -2.0 -- 0.0
  // 2%/E + 4%/sqrtE + 2%
  // A       B           C
  Smear::Acceptance::Zone EmcalOuterBackZone(ThetaFromEta ( 0.0 ),ThetaFromEta ( -2.0 ));
  Smear::Device EmcalOuterBack(Smear::kE, "sqrt( pow ( 0.02,2 ) + pow( 0.04,2)*E + pow ( 0.02*E,2 ) )");
  EmcalOuterBack.Accept.AddZone(EmcalOuterBackZone);
  EmcalOuterBack.Accept.SetGenre(Smear::kElectromagnetic);
  det.AddDevice(EmcalOuterBack);


  // Forward barrel - assume Matrix specs (waiting for Craig Woody's talk)
  // eta = 0 -- 1
  // E > 100 MeV  (50 MeV if higher resolution)
  //  not used because it suppresses smeared-up particles
  // 2%/E⊕(12-14)%/√E⊕(2-3)% for 30 cm space
  // choosing
  // 2%/E  + 14%/sqrtE + 3%
  // A       B           C
  Smear::Acceptance::Zone EmcalForwardBarrelZone(ThetaFromEta ( 1 ),ThetaFromEta ( 0 )); 
  Smear::Device EmcalForwardBarrel(Smear::kE, "sqrt( pow ( 0.02,2 ) + pow( 0.14,2)*E + pow ( 0.03*E,2 ) )");
  EmcalForwardBarrel.Accept.AddZone(EmcalForwardBarrelZone);
  EmcalForwardBarrel.Accept.SetGenre(Smear::kElectromagnetic);
  det.AddDevice(EmcalForwardBarrel);

  // Forward outer - assume Matrix specs (waiting for Craig Woody's talk)
  // eta = 1 -- 2
  // E > 50 MeV not used because it suppresses smeared-up particles
  // "2%/E + (4*-12)%/sqrtE  + 2%  Upper limit achievable with 40cm space
  // *Better resolution requires ~65 cm   space allocated"
  // choosing
  // 2%/E  + 12%/sqrtE + 2%  
  // A       B           C
  Smear::Acceptance::Zone EmcalFwdOuterZone(ThetaFromEta ( 2 ),ThetaFromEta ( 1 )); // E
  Smear::Device EmcalFwdOuter(Smear::kE, "sqrt( pow ( 0.02,2 ) + pow( 0.12,2)*E + pow ( 0.02*E,2 ) )");
  EmcalFwdOuter.Accept.AddZone(EmcalFwdOuterZone);
  EmcalFwdOuter.Accept.SetGenre(Smear::kElectromagnetic);
  det.AddDevice(EmcalFwdOuter);

  // > Shashlyk, Pb/Sc 5.5 x5.5, 8% stochastic, 2% constant.

  // Forward inner - Shashlyk guess
  // eta = 2 -- 3.5
  // E > 50 MeV not used because it suppresses smeared-up particles
  // Optimistic Shashlyk
  // 2%/E  + 8%/sqrtE + 2%  
  // A       B           C
  Smear::Acceptance::Zone EmcalFwdInnerZone(ThetaFromEta ( 3.5 ),ThetaFromEta ( 2 )); // E
  Smear::Device EmcalFwdInner(Smear::kE, "sqrt( pow ( 0.02,2 ) + pow( 0.8,2)*E + pow ( 0.02*E,2 ) )");
  EmcalFwdInner.Accept.AddZone(EmcalFwdInnerZone);
  EmcalFwdInner.Accept.SetGenre(Smear::kElectromagnetic);
  det.AddDevice(EmcalFwdInner);


  // Hadronic  Calorimeters
  // Assume marix specs (waiting for Oleg Tsai's talk)
  // Except, aspire to 35 % / sqrt(E)  + 2% forward
  // ----------------------
  // Note: kHadronic == |pdg|>110.

  // Back
  // eta = -3.5 -- -1
  // E>500 MeV  not used because it suppresses smeared-up particles
  // (Better resolution required more space and R&D)
  // stoch. = 50%, const = 10%
  Smear::Acceptance::Zone HcalBackZone(ThetaFromEta ( -1 ),ThetaFromEta ( -3.5 ) );
  Smear::Device HcalBack(Smear::kE, "sqrt(pow( 0.1*E, 2) + pow ( 0.5,2) *E)");
  HcalBack.Accept.AddZone(HcalBackZone);
  HcalBack.Accept.SetGenre(Smear::kHadronic);
  det.AddDevice(HcalBack);

  // Barrel
  // eta = -1 -- 1
  // E>500 MeV  not used because it suppresses smeared-up particles
  // stoch. = 100%, const = 10%
  Smear::Acceptance::Zone HcalBarrelZone(ThetaFromEta ( 1 ),ThetaFromEta ( -1 ) ); // E
  Smear::Device HcalBarrel(Smear::kE, "sqrt( pow( 0.1*E, 2) + pow( 1.0, 2)*E)");
  HcalBarrel.Accept.AddZone(HcalBarrelZone);
  HcalBarrel.Accept.SetGenre(Smear::kHadronic);
  det.AddDevice(HcalBarrel);

  // Forward
  // eta = 1 -- 3.5
  // E>500 MeV not used because it suppresses smeared-up particles
  // stoch. = 35%, const = 2% (aspirational goal)
  Smear::Acceptance::Zone HcalFwdZone(ThetaFromEta ( 3.5 ),ThetaFromEta ( 1 ));
  Smear::Device HcalFwd(Smear::kE, "sqrt(pow( 0.02*E, 2) + pow ( 0.35,2) *E)");
  HcalFwd.Accept.AddZone(HcalFwdZone);
  HcalFwd.Accept.SetGenre(Smear::kHadronic);
  det.AddDevice(HcalFwd);

  return det;
}

// -------------------------------------------------------------------
double ThetaFromEta( const double eta ) {
  if ( !std::isnan(eta) && !std::isinf(eta)   ) {
    return 2.0 * atan( exp( -eta ));
  }
  throw std::runtime_error("ThetaFromEta called with NaN or Inf");
  return -1;
}

// -------------------------------------------------------------------
std::valarray<double> CalcA ( const double Bfield,
			     const double x1, const std::valarray<double> A1,
			     const double x2, const std::valarray<double> A2 ){
  // TODO: Add tests:
  // - vectors have the same length
  // - B>=0
  // - x1<x2
  // - Result is exact at x1 and x2
  // - all vector components are >=0

  // Basic linear interpolation
  // f ( x ) = f(x1) + ( x-x1 ) * ( f(x2) - f(x1) ) / (x2 - x1 )
  return A1 + (Bfield - x1 ) * (A2 - A1) / (x2 - x1 );  
}
// -------------------------------------------------------------------
std::valarray<double> CalcB ( const double Bfield,
			     const double x1, const std::valarray<double> B1,
			     const double x2, const std::valarray<double> B2 ){
  // TODO: Add tests:
  // - vectors have the same length
  // - B>=0
  // - x1<x2
  // - Result is exact at x1 and x2
  // - all vector components are >=0

  // Basic linear interpolation
  // f ( x ) = f(x1) + ( x-x1 ) * ( f(x2) - f(x1) ) / (x2 - x1 )
  return B1 + (Bfield - x1 ) * (B2 - B1) / (x2 - x1 );  
}

// -------------------------------------------------------------------
void AssembleCoreTracker ( Smear::Detector& det,
			   const std::valarray<double>& eta_min, const std::valarray<double>& eta_max,
			   const std::valarray<double>& A, const std::valarray<double>& B ){

  // sigma_p/p ~ A% * p + B % 
  const TString SmearString = "sqrt( pow (  AAA *P*P, 2) + pow ( BBB * P, 2) )";  

  for( auto i = 0; i< A.size(); ++i)  {
    auto CurrentString = SmearString;
    TString AAA = ""; AAA += 0.01*A[i];
    TString BBB = ""; BBB += 0.01*B[i];
    CurrentString.ReplaceAll ( "AAA", AAA );
    CurrentString.ReplaceAll ( "BBB", BBB );

    Smear::Device TrackPos(Smear::kP, CurrentString);
    TrackPos.Accept.SetCharge(Smear::kCharged);
    Smear::Device TrackNeg = TrackPos;
  
    // Two zones, acceptance is symmetrical in eta
    Smear::Acceptance::Zone TrackZonePos( ThetaFromEta (  eta_max[i] ), ThetaFromEta (  eta_min[i] ));
    Smear::Acceptance::Zone TrackZoneNeg( ThetaFromEta ( -eta_min[i] ), ThetaFromEta ( -eta_max[i] ));

    TrackPos.Accept.AddZone(TrackZonePos);
    TrackNeg.Accept.AddZone(TrackZoneNeg);
    det.AddDevice(TrackPos);
    det.AddDevice(TrackNeg);

    // std::cout << i << "  " << CurrentString << std::endl;
  }  

  return;
}

