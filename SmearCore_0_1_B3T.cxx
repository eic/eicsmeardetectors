// Based on
// https://indico.bnl.gov/event/7913/contributions/41704/attachments/30564/47972/200924_Parametrizations.pdf
// https://indico.bnl.gov/event/11053/contributions/46969/attachments/33324/53540/core.pdf

// From March 29, 2021
// 
// Tracking assumes a 3T field
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

// declare static --> local to this file, won't clash with others
static double ThetaFromEta( const double eta );

Smear::Detector BuildCore_0_1_B3T() {
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
  
  // |eta| = 0 -- 0.5
  // sigma_p/p ~ 0.018 % p + 0.369 % 
  Smear::Device TrackPos1(Smear::kP, "sqrt( pow (  0.018 * 0.01 *P*P, 2) + pow ( 0.369 * 0.01 *P, 2) )");
  TrackPos1.Accept.SetCharge(Smear::kCharged);
  Smear::Device TrackNeg1 = TrackPos1;
  Smear::Acceptance::Zone TrackZonePos1(ThetaFromEta ( 0.5 ), ThetaFromEta (  0.0 ));
  Smear::Acceptance::Zone TrackZoneNeg1(ThetaFromEta ( 0.0 ), ThetaFromEta ( -0.5 ));
  TrackPos1.Accept.AddZone(TrackZonePos1);
  TrackNeg1.Accept.AddZone(TrackZoneNeg1);
  det.AddDevice(TrackPos1);
  det.AddDevice(TrackNeg1);

  // |eta| = 0.5 -- 1
  // sigma_p/p ~ 0.016 % p + 0.428 % 
  Smear::Device TrackPos2(Smear::kP, "sqrt( pow (  0.016 * 0.01 *P*P, 2) + pow ( 0.428 * 0.01 *P, 2) )");
  TrackPos2.Accept.SetCharge(Smear::kCharged);
  Smear::Device TrackNeg2 = TrackPos2;
  Smear::Acceptance::Zone TrackZonePos2(ThetaFromEta (  1.0 ), ThetaFromEta (  0.5 ));
  Smear::Acceptance::Zone TrackZoneNeg2(ThetaFromEta ( -0.5 ), ThetaFromEta ( -1.0 ));
  TrackPos2.Accept.AddZone(TrackZonePos2);
  TrackNeg2.Accept.AddZone(TrackZoneNeg2);
  det.AddDevice(TrackPos2);
  det.AddDevice(TrackNeg2);

  // |eta| = 1 -- 1.5
  // sigma_p/p ~ 0.016 % p + 0.427 % 
  Smear::Device TrackPos3(Smear::kP, "sqrt( pow (  0.016 * 0.01 *P*P, 2) + pow ( 0.427 * 0.01 *P, 2) )");
  TrackPos3.Accept.SetCharge(Smear::kCharged);
  Smear::Device TrackNeg3 = TrackPos3;
  Smear::Acceptance::Zone TrackZonePos3(ThetaFromEta (  1.5 ), ThetaFromEta (  1.0 ));
  Smear::Acceptance::Zone TrackZoneNeg3(ThetaFromEta ( -1.0 ), ThetaFromEta ( -1.5 ));
  TrackPos3.Accept.AddZone(TrackZonePos3);
  TrackNeg3.Accept.AddZone(TrackZoneNeg3);
  det.AddDevice(TrackPos3);
  det.AddDevice(TrackNeg3);


  // |eta| = 1.5 -- 2
  // sigma_p/p ~ 0.012 % p + 0.462 % 
  Smear::Device TrackPos4(Smear::kP, "sqrt( pow (  0.012 * 0.01 *P*P, 2) + pow ( 0.462 * 0.01 *P, 2) )");
  TrackPos4.Accept.SetCharge(Smear::kCharged);
  Smear::Device TrackNeg4 = TrackPos4;
  Smear::Acceptance::Zone TrackZonePos4(ThetaFromEta (  2.0 ), ThetaFromEta (  1.5 ));
  Smear::Acceptance::Zone TrackZoneNeg4(ThetaFromEta ( -1.5 ), ThetaFromEta ( -2.0 ));
  TrackPos4.Accept.AddZone(TrackZonePos4);
  TrackNeg4.Accept.AddZone(TrackZoneNeg4);
  det.AddDevice(TrackPos4);
  det.AddDevice(TrackNeg4);

  // |eta| = 2 -- 2.5
  // sigma_p/p ~ 0.018 % p + 0.719 % 
  Smear::Device TrackPos5(Smear::kP, "sqrt( pow (  0.018 * 0.01 *P*P, 2) + pow ( 0.719 * 0.01 *P, 2) )");
  TrackPos5.Accept.SetCharge(Smear::kCharged);
  Smear::Device TrackNeg5 = TrackPos5;
  Smear::Acceptance::Zone TrackZonePos5(ThetaFromEta (  2.5 ), ThetaFromEta (  2.0 ));
  Smear::Acceptance::Zone TrackZoneNeg5(ThetaFromEta ( -2.0 ), ThetaFromEta ( -2.5 ));
  TrackPos5.Accept.AddZone(TrackZonePos5);
  TrackNeg5.Accept.AddZone(TrackZoneNeg5);
  det.AddDevice(TrackPos5);
  det.AddDevice(TrackNeg5);

  // |eta| = 2.5 -- 3
  // sigma_p/p ~ 0.039 % p + 1.336 % 
  Smear::Device TrackPos6(Smear::kP, "sqrt( pow (  0.039 * 0.01 *P*P, 2) + pow ( 1.336 * 0.01 *P, 2) )");
  TrackPos6.Accept.SetCharge(Smear::kCharged);
  Smear::Device TrackNeg6 = TrackPos6;
  Smear::Acceptance::Zone TrackZonePos6(ThetaFromEta (  3.0 ), ThetaFromEta (  2.5 ));
  Smear::Acceptance::Zone TrackZoneNeg6(ThetaFromEta ( -2.5 ), ThetaFromEta ( -3.0 ));
  TrackPos6.Accept.AddZone(TrackZonePos6);
  TrackNeg6.Accept.AddZone(TrackZoneNeg6);
  det.AddDevice(TrackPos6);
  det.AddDevice(TrackNeg6);


  // |eta| = 3 -- 3.5
  // sigma_p/p ~ 0.103 % p + 2.428 % 
  Smear::Device TrackPos7(Smear::kP, "sqrt( pow (  0.103 * 0.01 *P*P, 2) + pow ( 2.428 * 0.01 *P, 2) )");
  TrackPos7.Accept.SetCharge(Smear::kCharged);
  Smear::Device TrackNeg7 = TrackPos7;
  Smear::Acceptance::Zone TrackZonePos7(ThetaFromEta (  3.5 ), ThetaFromEta (  3.0 ));
  Smear::Acceptance::Zone TrackZoneNeg7(ThetaFromEta ( -3.0 ), ThetaFromEta ( -3.5 ));
  TrackPos7.Accept.AddZone(TrackZonePos7);
  TrackNeg7.Accept.AddZone(TrackZoneNeg7);
  det.AddDevice(TrackPos7);
  det.AddDevice(TrackNeg7);

  // |eta| = 3.5 -- 4.0
  // sigma_p/p ~ 0.295 % p + 4.552 % 
  Smear::Device TrackPos8(Smear::kP, "sqrt( pow (  0.295 * 0.01 *P*P, 2) + pow ( 4.552 * 0.01 *P, 2) )");
  TrackPos8.Accept.SetCharge(Smear::kCharged);
  Smear::Device TrackNeg8 = TrackPos8;
  Smear::Acceptance::Zone TrackZonePos8(ThetaFromEta (  4.0 ), ThetaFromEta (  3.5 ));
  Smear::Acceptance::Zone TrackZoneNeg8(ThetaFromEta ( -3.5 ), ThetaFromEta ( -4.0 ));
  TrackPos8.Accept.AddZone(TrackZonePos8);
  TrackNeg8.Accept.AddZone(TrackZoneNeg8);
  det.AddDevice(TrackPos8);
  det.AddDevice(TrackNeg8);

  
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

  // Forward - assume Matrix specs (waiting for Craig Woody's talk)
  // eta = 1 -- 3.5
  // E > 50 MeV not used because it suppresses smeared-up particles
  // "2%/E + (4*-12)%/sqrtE  + 2%  Upper limit achievable with 40cm space
  // *Better resolution requires ~65 cm   space allocated"
  // choosing
  // 2%/E  + 12%/sqrtE + 2%  
  // A       B           C
  Smear::Acceptance::Zone EmcalFwdZone(ThetaFromEta ( 3.5 ),ThetaFromEta ( 1 )); // E
  Smear::Device EmcalFwd(Smear::kE, "sqrt( pow ( 0.02,2 ) + pow( 0.12,2)*E + pow ( 0.02*E,2 ) )");
  EmcalFwd.Accept.AddZone(EmcalFwdZone);
  EmcalFwd.Accept.SetGenre(Smear::kElectromagnetic);
  det.AddDevice(EmcalFwd);

  // Hadronic  Calorimeters
  // Assume marix specs (waiting for Oleg Tsai's talk)
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
  // stoch. = 50%, const = 10%
  // "(35%/sqrtE not achievable)"
  Smear::Acceptance::Zone HcalFwdZone(ThetaFromEta ( 3.5 ),ThetaFromEta ( 1 ));
  Smear::Device HcalFwd(Smear::kE, "sqrt(pow( 0.1*E, 2) + pow ( 0.5,2) *E)");
  HcalFwd.Accept.AddZone(HcalFwdZone);
  HcalFwd.Accept.SetGenre(Smear::kHadronic);
  det.AddDevice(HcalFwd);

  // Not covered:
  // Tracker material budget X/X0 <~5%
  // Tracker pointing resolution
  // Low-Q^2 tagger: -6.9<eta<-5.8: Delta_theta/theta < 1.5%; 10^-6 < Q2 < 10^-2 GeV2
  // Proton spectrometer:  eta>6.2: sigma_intrinsic(|t|)/|t| < 1%; Acceptance: 0.2 < pT < 1.2 GeV/c
  // Barrel vertexing: sigma_xyz ~ 20 microns, d0(z) ~ d0(r phi) ~ (20 microns)/(pT [GeV])  + 5 microns
  // Central muon detection

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
