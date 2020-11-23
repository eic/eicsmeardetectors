// Based on
// https://physdiv.jlab.org/DetectorMatrix/
// From November 21, 2020
// 
// Tracking assumes a 1.5T field
// 
// Important Notes:
// - Where ranges are given, the more conservative number is chosen.
// - Without available specifications, angular resolution is assumed to be perfect.
// - Momentum and energy acceptance specifications are NOT implemented
// -- a) "90% acceptance" is not specific enough to introduce efficiency
// -- b) It is not a priori clear whether cuts should apply to truth or smeared value
// --    A "looper" may not reach the detector.
// --    On the other hand, a calorimetry deposit may fluctuate above threshold
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

Smear::Detector BuildMatrixDetector_0_2_B1_5T() {
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
  // total coverage of the handbook for tracker, ecal, and hcal is -3.5 < eta < 3.5
  Smear::Acceptance::Zone AngleZoneCommon(ThetaFromEta ( 3.5 ),ThetaFromEta ( -3.5 ));
  Smear::Device SmearThetaCommon(Smear::kTheta, "0.0");
  SmearThetaCommon.Accept.AddZone(AngleZoneCommon);
  SmearThetaCommon.Accept.SetGenre(Smear::kAll);
  det.AddDevice(SmearThetaCommon);

  Smear::Device SmearPhiCommon(Smear::kPhi, "0.0");
  SmearPhiCommon.Accept.AddZone(AngleZoneCommon);
  SmearPhiCommon.Accept.SetGenre(Smear::kAll);
  det.AddDevice(SmearPhiCommon);

  // Tracking  (B = 1.5 T)
  // --------
  // Note: Smear::kCharged checks pdg charge, so includes muons (good)
  // NOT implemented:
  // Minimum pT for B = 1.5 T:
  // 100 MeV/c for -3.0 < eta < -2.5
  // 130 MeV/c for -2.5 < eta < -2.0
  // 70 MeV/c for -2.0 < eta < -1.5
  // 150 MeV/c for -1.5 < eta < -1.0
  
  // eta = -3.5 --  -2.5
  // sigma_p/p ~ 0.2% p + 5%
  Smear::Acceptance::Zone TrackBack1Zone(ThetaFromEta ( -2.5 ),ThetaFromEta ( -3.5 ));
  Smear::Device TrackBack1P(Smear::kP, "sqrt( pow ( 0.002*P*P, 2) + pow ( 0.05*P, 2) )");
  TrackBack1P.Accept.AddZone(TrackBack1Zone);
  TrackBack1P.Accept.SetCharge(Smear::kCharged);
  det.AddDevice(TrackBack1P);

  // eta = -2.5 --  -1
  // sigma_p/p ~ 0.04% p + 2%
  Smear::Acceptance::Zone TrackBack2Zone(ThetaFromEta ( -1 ),ThetaFromEta ( -2.5 ));
  Smear::Device TrackBack2P(Smear::kP, "sqrt( pow ( 0.0004*P*P, 2) + pow ( 0.02*P, 2) )");
  TrackBack2P.Accept.AddZone(TrackBack2Zone);
  TrackBack2P.Accept.SetCharge(Smear::kCharged);
  det.AddDevice(TrackBack2P);

  // eta = -1 -- +1
  // sigma_p/p ~ 0.04% p + 1%
  Smear::Acceptance::Zone TrackBarrelZone(ThetaFromEta ( 1 ),ThetaFromEta ( -1 ));
  Smear::Device TrackBarrelP(Smear::kP, "sqrt( pow ( 0.0004*P*P, 2) + pow ( 0.01*P, 2) )");
  TrackBarrelP.Accept.AddZone(TrackBarrelZone);
  TrackBarrelP.Accept.SetCharge(Smear::kCharged);
  det.AddDevice(TrackBarrelP);

  // eta = 1 -- 2.5
  // sigma_p/p ~ 0.04% p+2%
  Smear::Acceptance::Zone TrackFwd2Zone(ThetaFromEta ( 2.5 ),ThetaFromEta ( 1 ));
  Smear::Device TrackFwd2P(Smear::kP, "sqrt( pow ( 0.0004*P*P, 2) + pow ( 0.02*P, 2) )");
  TrackFwd2P.Accept.AddZone(TrackFwd2Zone);
  TrackFwd2P.Accept.SetCharge(Smear::kCharged);
  det.AddDevice(TrackFwd2P);

  // eta = 2.5 -- 3.5
  // sigma_p/p ~ 0.2% p+5%
  Smear::Acceptance::Zone TrackFwd1Zone(ThetaFromEta ( 3.5 ),ThetaFromEta ( 2.5 ));
  Smear::Device TrackFwd1P(Smear::kP, "sqrt( pow ( 0.002*P*P, 2) + pow ( 0.05*P, 2) )");
  TrackFwd1P.Accept.AddZone(TrackFwd1Zone);
  TrackFwd1P.Accept.SetCharge(Smear::kCharged);
  det.AddDevice(TrackFwd1P);

  // PID
  // ---
  // Perfect pi/K/P PID in the region indicating >3 sigma separation
  // Make sure to not cover more than is covered by the other detectors.
  // No minimum momentum is used.
  // Accept charged hadrons (excludes e, mu)

  // Back
  // eta = -3.5 -- -1
  // p < 10 GeV
  Smear::Acceptance::Zone PidBackZone(ThetaFromEta(-1),ThetaFromEta(-3.5),
				      0., TMath::TwoPi(), // phi
				      0., TMath::Infinity(), // E
				      0., 10 ); // p
  Smear::PerfectID PidBack;
  PidBack.Accept.AddZone(PidBackZone);
  PidBack.Accept.SetCharge(Smear::kCharged);
  PidBack.Accept.SetGenre(Smear::kHadronic);
  det.AddDevice( PidBack );

  // Barrel
  // eta = -1 -- 1
  // p < 6 GeV
  Smear::Acceptance::Zone PidBarrelZone(ThetaFromEta(1),ThetaFromEta(-1),
				      0., TMath::TwoPi(), // phi
				      0., TMath::Infinity(), // E
				      0., 6 ); // p
  Smear::PerfectID PidBarrel;
  PidBarrel.Accept.AddZone(PidBarrelZone);
  PidBarrel.Accept.SetCharge(Smear::kCharged);
  PidBarrel.Accept.SetGenre(Smear::kHadronic);
  det.AddDevice( PidBarrel );

  // Forward
  // eta = -1 -- 1
  // p < 50 GeV
  // not sure what to make of "(worse approaching 3.5)"
  Smear::Acceptance::Zone PidFwdZone(ThetaFromEta(3.5),ThetaFromEta(1),
				      0., TMath::TwoPi(), // phi
				      0., TMath::Infinity(), // E
				      0., 50 ); // p
  Smear::PerfectID PidFwd;
  PidFwd.Accept.AddZone(PidFwdZone);
  PidFwd.Accept.SetCharge(Smear::kCharged);
  PidFwd.Accept.SetGenre(Smear::kHadronic);
  det.AddDevice( PidFwd );

  // EM Calorimeters
  // ---------------
  // Note: Smear::kElectromagnetic == gamma + e. Does not include muons (good)
  // Calorimeter resolution usually given as sigma_E/E = A% / E + B%/Sqrt{E} + C%
  // EIC Smear needs absolute sigma: sigma_E = Sqrt{A*A + B*B*E + C*C*E*E}

  // Back
  // eta = -3.5 -- -2
  // E> 20 MeV not used because it suppresses smeared-up particles
  // 1%/E + 2.5%/sqrtE + 1%
  // A       B           C
  Smear::Acceptance::Zone EmcalBackZone(ThetaFromEta ( -2 ),ThetaFromEta ( -3.5 ));
  Smear::Device EmcalBack(Smear::kE, "sqrt( pow ( 0.01,2 ) + pow( 0.025,2)*E + pow ( 0.01*E,2 ) )");
  EmcalBack.Accept.AddZone(EmcalBackZone);
  EmcalBack.Accept.SetGenre(Smear::kElectromagnetic);
  det.AddDevice(EmcalBack);

  // MidBack
  // eta = -2 -- -1
  // E> 50 MeV not used because it suppresses smeared-up particles
  // 2%/E   +  (4-8)%/sqrtE   + 2% (Upper limit achievable with 50 cm space. Better resolution requires ~65 cm   space allocated  )
  // A       B choose 8%        C  
  Smear::Acceptance::Zone EmcalMidBackZone(ThetaFromEta ( -1 ),ThetaFromEta ( -2 ));
  Smear::Device EmcalMidBack(Smear::kE, "sqrt( pow ( 0.02,2 ) + pow( 0.08,2)*E + pow ( 0.02*E,2 ) )");
  EmcalMidBack.Accept.AddZone(EmcalMidBackZone);
  EmcalMidBack.Accept.SetGenre(Smear::kElectromagnetic);
  det.AddDevice(EmcalMidBack);

  // Barrel
  // eta = -1 -- 1
  // E > 100 MeV  (50 MeV if higher resolution)
  //  not used because it suppresses smeared-up particles
  // 2%/E⊕(12-14)%/√E⊕(2-3)% for 30 cm space
  // A better stochastic term can be achieved with more space:
  // 2.5% with crystals    35cm
  // 10% sampling          40cm
  // 4% SciGlass            65cm
  // (*Better resolution requires ~65 cm   space allocated  )"
  // choosing
  // 2%/E  + 14%/sqrtE + 3%
  // A       B           C
  Smear::Acceptance::Zone EmcalBarrelZone(ThetaFromEta ( 1 ),ThetaFromEta ( -1 )); 
  Smear::Device EmcalBarrel(Smear::kE, "sqrt( pow ( 0.02,2 ) + pow( 0.14,2)*E + pow ( 0.03*E,2 ) )");
  EmcalBarrel.Accept.AddZone(EmcalBarrelZone);
  EmcalBarrel.Accept.SetGenre(Smear::kElectromagnetic);
  det.AddDevice(EmcalBarrel);

  // Forward
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
