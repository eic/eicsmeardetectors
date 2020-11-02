// Based on
// https://physdiv.jlab.org/DetectorMatrix/
// Using everything from v 0.1 from June 16 2020 except:
// - EMCal coverage reduced to |eta|<3.5
// - Tracking parameters from B=1.5T Matrix preview at
//   https://indico.bnl.gov/event/9984/contributions/43066/attachments/31173/49186/YR_Detector_Matrix_Tracking_only_10282020.xlsx

// Reminder:
// Acceptance::Zone(double theta = 0., double = TMath::Pi(),
//      double phi = 0., double = TMath::TwoPi(),
//      double E = 0., double = TMath::Infinity(),
//      double p = 0., double = TMath::Infinity(),
//      double pt = 0., double = TMath::Infinity(),
//      double pz = -TMath::Infinity(), double = TMath::Infinity());


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

Smear::Detector BuildTrackingPreview_0_2_B1_5T() {
  gSystem->Load("libeicsmear");

  // Create the detector object to hold all devices
  Smear::Detector det;

  // The framework provides implementations of three kinematic calculation methods
  // from smeared values
  // NM - "Null method" uses the scattered lepton.
  // DA - Double Angle method
  // JB - Jacquet-Blondel method
  // Important Notes:
  // - All methods rely on measured energy and momentum components (lepton for NM, hadrons for DA, JB).
  //   In order to rely on these methods, you therefore _need_ to cover the relevant
  //   regions with smearers for momentum, angles and energy.
  // - This is exacerbated by the fact that in the current implementation a smearing of e.g. P
  //   may fool the framework into assuming theta is measured to be the initialization value 0,
  //   leading to nan or inf or otherwise wrong calculations.
  //   NM catches P=0 and replaces it in calculations with E,
  //   but JB and DA either don't work at all or report nonsenical values.
  // - It may be advantageous to use a P measurement in place of E in a more general case.
  //   In the future, we plan to change to a weighted mean approach.
  // The summary of the above is that the user should be mindful of the limitations and assumptions in
  // the calculation of smeared kinematics. Sophisticated calculations are best redone in user code
  // from the smeared particles directly.
  det.SetEventKinematicsCalculator("NM DA JB");

  // Perfect phi and theta for all particles
  // ---------------------------------------
  // TODO: Better options?
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

  // TODO: What to make of this?
  // "100 MeV/c with 50% acceptance (similar for pi and K)"
  // Add efficiency?
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

  // EM Calorimeters
  // ---------------
  // Note: Smear::kElectromagnetic == gamma + e. Does not include muons (good)

  // Calorimeter resolution usually given as sigma_E/E = const% + stocastic%/Sqrt{E}
  // EIC Smear needs absolute sigma: sigma_E = Sqrt{const*const*E*E + stoc*stoc*E}

  // Back
  // eta = -3.5 -- -2
  // stoch. = 2%
  Smear::Acceptance::Zone EmcalBackZone(ThetaFromEta ( -2 ),ThetaFromEta ( -3.5 ));
  Smear::Device EmcalBack(Smear::kE, "sqrt( pow ( 0.0*E,2 ) + pow( 0.02,2)*E)");
  EmcalBack.Accept.AddZone(EmcalBackZone);
  EmcalBack.Accept.SetGenre(Smear::kElectromagnetic);
  det.AddDevice(EmcalBack);

  // MidBack
  // eta = -2 -- -1
  // stoch. = 7%
  Smear::Acceptance::Zone EmcalMidBackZone(ThetaFromEta ( -1 ),ThetaFromEta ( -2 ));
  Smear::Device EmcalMidBack(Smear::kE, "sqrt( pow ( 0.0*E,2 ) + pow( 0.07,2)*E)");
  EmcalMidBack.Accept.AddZone(EmcalMidBackZone);
  EmcalMidBack.Accept.SetGenre(Smear::kElectromagnetic);
  det.AddDevice(EmcalMidBack);

  // Forward
  // eta = -1 -- 3.5
  // stoch. = 10-12%, use 12%
  Smear::Acceptance::Zone EmcalFwdZone(ThetaFromEta ( 3.5 ),ThetaFromEta ( -1 ));
  Smear::Device EmcalFwd(Smear::kE, "sqrt( pow ( 0.0*E,2 ) + pow( 0.12,2)*E)");
  EmcalFwd.Accept.AddZone(EmcalFwdZone);
  EmcalFwd.Accept.SetGenre(Smear::kElectromagnetic);
  det.AddDevice(EmcalFwd);

  // Could turn on perfect PID
  // Make sure to not cover more than is covered by the other detectors.
  // Using the smallest region here, but we could add a second one for
  // the extended EmCal range
  // Smear::Acceptance::Zone acceptpid(ThetaFromEta(3.5),ThetaFromEta(-3.5));
  // Smear::PerfectID pid;
  // pid.Accept.AddZone(acceptpid);
  // det.AddDevice( pid );

  // Hadronic  Calorimeters
  // ----------------------
  // Note: kHadronic == |pdg|>110.

  // Back
  // eta = -3.5 -- -1
  // stoch. = 50%
  Smear::Acceptance::Zone HcalBackZone(ThetaFromEta ( -1 ),ThetaFromEta ( -3.5 ));
  Smear::Device HcalBack(Smear::kE, "sqrt(pow( 0.0*E, 2) + pow ( 0.5,2) *E)");
  HcalBack.Accept.AddZone(HcalBackZone);
  HcalBack.Accept.SetGenre(Smear::kHadronic);
  det.AddDevice(HcalBack);

  // Barrel
  // eta = -1 -- 1
  // The matrix has nothing. As examples, one could turn to
  // ~CMS
  // Smear::Device HcalBarrel(Smear::kE, "sqrt( pow( 0.07*E, 2) + pow( 0.85, 2)*E)");
  // ~Zeus
  // Smear::Device HcalBarrel(Smear::kE, "sqrt( pow( 0.02*E, 2) + pow( 0.35,2) *E)");

  // Smear::Acceptance::Zone HcalBarrelZone(ThetaFromEta ( 1 ),ThetaFromEta ( -1 ));
  // HcalBarrel.Accept.AddZone(HcalBarrelZone);
  // HcalBarrel.Accept.SetGenre(Smear::kHadronic);
  // det.AddDevice(HcalBarrel);

  // Forward
  // eta = 1 -- 3.5
  // stoch. = 50%
  Smear::Acceptance::Zone HcalFwdZone(ThetaFromEta ( 3.5 ),ThetaFromEta ( 1 ));
  Smear::Device HcalFwd(Smear::kE, "sqrt(pow( 0.0*E, 2) + pow ( 0.5,2) *E)");
  HcalFwd.Accept.AddZone(HcalFwdZone);
  HcalFwd.Accept.SetGenre(Smear::kHadronic);
  det.AddDevice(HcalFwd);

  // Not covered:
  // Tracker material budget X/X0 <~5%
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
