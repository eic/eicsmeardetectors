// Based on
// https://physdiv.jlab.org/DetectorMatrix/
// From June 16 2020

// Here added devices for the far forward region on request by the WG

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

#include <string>
using std::string;

// declare static --> local to this file, won't clash with others
static double ThetaFromEta( const double eta );

/**
   beam_mom_nn: ion beam momentum per nucleon in GeV. Using int to avoid rounding issues in switch
*/

Smear::Detector BuildMatrixDetector_0_1_FF( const int beam_mom_nn  ) {
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

  // IMPORTANT: There are two traps (will be addressed in future releases):
  //            1) If you smear E but don't provide phi, theta smearers, those values will be
  //               set to 0, not to a fault value and not to the truth level
  //            2) If you do provide such a smearer, pt and pz will be changed
  //               by a consistency enforcer in Detector::Smear()


  // Perfect phi and theta for all particles
  // ---------------------------------------
  // total coverage of the handbook for tracker and hcal is -3.5 < eta < 3.5
  Smear::Acceptance::Zone AngleZoneHadronic(ThetaFromEta ( 3.5 ),ThetaFromEta ( -3.5 ));
  Smear::Device SmearThetaHadronic(Smear::kTheta, "0.0");
  SmearThetaHadronic.Accept.AddZone(AngleZoneHadronic);
  SmearThetaHadronic.Accept.SetGenre(Smear::kHadronic);
  det.AddDevice(SmearThetaHadronic);

  Smear::Device SmearPhiHadronic(Smear::kPhi, "0.0");
  SmearPhiHadronic.Accept.AddZone(AngleZoneHadronic);
  SmearPhiHadronic.Accept.SetGenre(Smear::kHadronic);
  det.AddDevice(SmearPhiHadronic);

  // muons are neither hadrons nor electromgnetic
  Smear::Acceptance::Zone AngleZoneMuon(ThetaFromEta ( 3.5 ),ThetaFromEta ( -3.5 ));
  Smear::Device SmearThetaMuon(Smear::kTheta, "0.0");
  SmearThetaMuon.Accept.AddZone(AngleZoneMuon);
  SmearThetaMuon.Accept.AddParticle(13);
  SmearThetaMuon.Accept.AddParticle(-13);
  det.AddDevice(SmearThetaMuon);

  Smear::Device SmearPhiMuon(Smear::kPhi, "0.0");
  SmearPhiMuon.Accept.AddZone(AngleZoneMuon);
  SmearPhiMuon.Accept.AddParticle(13);
  SmearPhiMuon.Accept.AddParticle(-13);
  det.AddDevice(SmearPhiMuon);

  // emcal stretches to -4.5 < eta < 4.5
  Smear::Acceptance::Zone AngleZoneEmcal(ThetaFromEta ( 4.5 ),ThetaFromEta ( -4.5 ));
  Smear::Device SmearThetaEmcal(Smear::kTheta, "0.0");
  SmearThetaEmcal.Accept.AddZone(AngleZoneEmcal);
  SmearThetaEmcal.Accept.SetGenre(Smear::kElectromagnetic);
  det.AddDevice(SmearThetaEmcal);

  Smear::Device SmearPhiEmcal(Smear::kPhi, "0.0");
  SmearPhiEmcal.Accept.AddZone(AngleZoneEmcal);
  SmearPhiEmcal.Accept.SetGenre(Smear::kElectromagnetic);
  det.AddDevice(SmearPhiEmcal);

  // Tracking
  // --------
  // Note: Smear::kCharged checks pdg charge, so includes muons (good)
  // eta = -3.5 --  -2.0
  // sigma_p/p ~ 0.1% p+0.5%
  Smear::Acceptance::Zone TrackBack1Zone(ThetaFromEta ( -2.5 ),ThetaFromEta ( -3.5 ));
  Smear::Device TrackBack1P(Smear::kP, "sqrt( pow ( 0.001*P*P, 2) + pow ( 0.005*P, 2) )");
  TrackBack1P.Accept.AddZone(TrackBack1Zone);
  TrackBack1P.Accept.SetCharge(Smear::kCharged);
  det.AddDevice(TrackBack1P);

  // eta = -2.0 --  -1
  // sigma_p/p ~ 0.05% p+ 0.5%
  Smear::Acceptance::Zone TrackBack2Zone(ThetaFromEta ( -1 ),ThetaFromEta ( -2.5 ));
  Smear::Device TrackBack2P(Smear::kP, "sqrt( pow ( 0.0005*P*P, 2) + pow ( 0.005*P, 2) )");
  TrackBack2P.Accept.AddZone(TrackBack2Zone);
  TrackBack2P.Accept.SetCharge(Smear::kCharged);
  det.AddDevice(TrackBack2P);

  // eta = -1 -- +1
  // sigma_p/p ~ 0.05% p+0.5%
  Smear::Acceptance::Zone TrackBarrelZone(ThetaFromEta ( 1 ),ThetaFromEta ( -1 ));
  Smear::Device TrackBarrelP(Smear::kP, "sqrt( pow ( 0.0005*P*P, 2) + pow ( 0.005*P, 2) )");
  TrackBarrelP.Accept.AddZone(TrackBarrelZone);
  TrackBarrelP.Accept.SetCharge(Smear::kCharged);
  det.AddDevice(TrackBarrelP);

  // eta = 1 -- 2.5
  // sigma_p/p ~ 0.05% p+1.0%
  Smear::Acceptance::Zone TrackFwd2Zone(ThetaFromEta ( 2.5 ),ThetaFromEta ( 1 ));
  Smear::Device TrackFwd2P(Smear::kP, "sqrt( pow ( 0.0005*P*P, 2) + pow ( 0.01*P, 2) )");
  TrackFwd2P.Accept.AddZone(TrackFwd2Zone);
  TrackFwd2P.Accept.SetCharge(Smear::kCharged);
  det.AddDevice(TrackFwd2P);

  // eta = 2.5 -- 3.5
  // sigma_p/p ~ 0.1% p+2.0%
  Smear::Acceptance::Zone TrackFwd1Zone(ThetaFromEta ( 3.5 ),ThetaFromEta ( 2.5 ));
  Smear::Device TrackFwd1P(Smear::kP, "sqrt( pow ( 0.001*P*P, 2) + pow ( 0.02*P, 2) )");
  TrackFwd1P.Accept.AddZone(TrackFwd1Zone);
  TrackFwd1P.Accept.SetCharge(Smear::kCharged);
  det.AddDevice(TrackFwd1P);

  // EM Calorimeters
  // ---------------
  // Note: Smear::kElectromagnetic == gamma + e. Does not include muons (good)

  // Calorimeter resolution usually given as sigma_E/E = const% + stochastic%/Sqrt{E}
  // EIC Smear needs absolute sigma: sigma_E = Sqrt{const*const*E*E + stoc*stoc*E}

  // Back
  // eta = -4.5 -- -2
  // stoch. = 2%
  Smear::Acceptance::Zone EmcalBackZone(ThetaFromEta ( -2 ),ThetaFromEta ( -4.5 ));
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
  // eta = -1 -- 4.5
  // stoch. = 10-12%, use 12%
  Smear::Acceptance::Zone EmcalFwdZone(ThetaFromEta ( 4.5 ),ThetaFromEta ( -1 ));
  Smear::Device EmcalFwd(Smear::kE, "sqrt( pow ( 0.0*E,2 ) + pow( 0.12,2)*E)");
  EmcalFwd.Accept.AddZone(EmcalFwdZone);
  EmcalFwd.Accept.SetGenre(Smear::kElectromagnetic);
  det.AddDevice(EmcalFwd);

  // TODO: Add PID
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

  // Far forward
  // -----------
  //   Neutrons:
  // - Assume uniform acceptancefor 0<theta< 4.5 mrad
  // - Assume an overall energy resolution sigmaE/E = 50% / sqrt (E) oplus 5%
  // - Assume angular resolution of sigmaTheta = 3 mrad/rootE
  // sigma_E/E = stochastic%/Sqrt{E} + const% 
  // EIC Smear needs absolute sigma: sigma_E = Sqrt{const*const*E*E + stoc*stoc*E}
  // Note: In principle,   ZDC{|theta|phi}.Accept.SetCharge(Smear::kNeutral); is the more correct
  // way - which would also accept K0L for example.
  // Anticipating the actual needs of users, we're going to only accept neutrons and gammas for now
  Smear::Acceptance::Zone ZDCzone( 1e-7, 4.5e-3 );
  Smear::Device ZDC(Smear::kE, "sqrt(pow( 0.05*E, 2) + pow ( 0.5,2) *E)");
  ZDC.Accept.AddZone(ZDCzone);
  ZDC.Accept.AddParticle(2112);
  ZDC.Accept.AddParticle(22);
  det.AddDevice(ZDC);
  
  Smear::Device ZDCtheta(Smear::kTheta, "3e-3 / sqrt(E)");
  ZDCtheta.Accept.AddZone(ZDCzone);
  ZDCtheta.Accept.AddParticle(2112);
  ZDCtheta.Accept.AddParticle(22);
  det.AddDevice(ZDCtheta);
  
  Smear::Device ZDCphi(Smear::kPhi, "0");
  ZDCphi.Accept.AddZone(ZDCzone);
  ZDCphi.Accept.AddParticle(2112);
  ZDCphi.Accept.AddParticle(22);
  det.AddDevice(ZDCphi);

  // Protons
  // All detectors: Reasonable to assume sigma_p/p= 0.5% sigmaPt/Pt = 3%
  // CHANGE on September 10, 2020: 5% was a typo in the source, it should be 0.5%
  std::string pformula  = "0.005*P";
  std::string ptformula = "0.03*pT";

  // Assume uniform acceptance for 6<theta <20 mrad – "B0 spectrometer"
  // Note that anti-protons bend the other way
  Smear::Acceptance::Zone B0zone( 6e-3, 20e-3 );
  Smear::Device B0P(Smear::kP, pformula);
  B0P.Accept.AddZone(B0zone);
  B0P.Accept.AddParticle(2212);
  det.AddDevice(B0P);

  Smear::Device B0Pt(Smear::kPt, ptformula);
  B0Pt.Accept.AddZone(B0zone);
  B0Pt.Accept.AddParticle(2212);
  det.AddDevice(B0Pt);
  
  Smear::Device B0phi(Smear::kPhi, "0");
  B0phi.Accept.AddZone(B0zone);
  B0phi.Accept.AddParticle(2212);
  det.AddDevice(B0phi);

  // For protons with p_z/(beam momentum / nucleus )>.6 – "Roman pots"
  auto RP_minpz = 0.6 * beam_mom_nn;
  // 275 GeV -or- 135 GeV/n deuterons: Assume uniform acceptance for .5<theta<5.0 mrad
  // 100 GeV: Assume uniform acceptance for .2<theta<5.0 mrad
  // 41 GeV: Assume uniform acceptance for 1.0<theta<4.5 mrad
  //
  // for protons from nuclear breakup, the TOTAL momentum of the beam must be specified
  // so 41 GeV/n He-3 is 61 GeV total, and 41 GeV/n deuteron is 82 GeV total
  float thetamin = 0;
  float thetamax = 0; 
  switch ( beam_mom_nn  ){ // switch needs an int. add a little to avoid rounding problems
  case 275 : // e+P
  case 135 : // e+D
    thetamin = 0.5e-3;
    thetamax = 5e-3;
    break;
  case 110 :
  case 100 :
    thetamin = 0.2e-3;
    thetamax = 5e-3;
    break;
  case 41 :
    thetamin = 1.0e-3;
    thetamax = 4.5e-3;
    break;
  case 61 :   // 41 GeV/n He-3 beam setting
  case 82 :   // 41 GeV/n deuteron beam setting
  case 165 :  // 110 GeV/n He-3 beam setting
  case 220 :  // 110 GeV/n deuteron beam setting
    thetamin = 1.0e-6;
    thetamax = 5.0e-3;
	break;
  default :
    throw std::runtime_error ( "Unsupported beam momentum for far forward detectors");
  }

  Smear::Acceptance::Zone RPzone( thetamin, thetamax,
				  0, TMath::TwoPi(), // phi
				  0., TMath::Infinity(), // E
				  0., TMath::Infinity(), // p
				  0., TMath::Infinity(), // pt
				  RP_minpz, TMath::Infinity() // pz
				  );
      
  Smear::Device RPP(Smear::kP, pformula);
  RPP.Accept.AddZone(RPzone);
  RPP.Accept.AddParticle(2212);
  det.AddDevice(RPP);

  Smear::Device RPPt(Smear::kPt, ptformula);
  RPPt.Accept.AddZone(RPzone);
  RPPt.Accept.AddParticle(2212);
  det.AddDevice(RPPt);
  
  Smear::Device RPphi(Smear::kPhi, "0");
  RPphi.Accept.AddZone(RPzone);
  RPphi.Accept.AddParticle(2212);
  det.AddDevice(RPphi);

  // For protons with .25<p_z/(beam momentum)<.6 – "Off-momentum Detectors"
  auto OM_minpz = 0.25 * beam_mom_nn;
  auto OM_maxpz = RP_minpz;

  // Assume uniform acceptance for 0.0<theta<2.0 mrad  
  Smear::Acceptance::Zone OMfullzone( 1e-7, 2e-3,
  				      0, TMath::TwoPi(), // phi
  				      0., TMath::Infinity(), // E
  				      0., TMath::Infinity(), // p
  				      0., TMath::Infinity(), // pt
  				      OM_minpz, OM_maxpz // pz
  				      );

    
  Smear::Device OMfullP(Smear::kP, pformula);
  OMfullP.Accept.AddZone(OMfullzone);
  OMfullP.Accept.AddParticle(2212);
  det.AddDevice(OMfullP);

  Smear::Device OMfullPt(Smear::kPt, ptformula);
  OMfullPt.Accept.AddZone(OMfullzone);
  OMfullPt.Accept.AddParticle(2212);
  det.AddDevice(OMfullPt);
  
  Smear::Device OMfullphi(Smear::kPhi, "0");
  OMfullphi.Accept.AddZone(OMfullzone);
  OMfullphi.Accept.AddParticle(2212);
  det.AddDevice(OMfullphi);

  // for 2.0<theta<5.0 mrad, only accepted for |phi|>1 radian
  Smear::Acceptance::Zone OMpartialzone( 2e-3, 5e-3,
  					 1, TMath::TwoPi()-1, // phi
  					 0., TMath::Infinity(), // E
  					 0., TMath::Infinity(), // p
  					 0., TMath::Infinity(), // pt
  					 OM_minpz, OM_maxpz // pz
  					 );

  Smear::Device OMpartialP(Smear::kP, pformula);
  OMpartialP.Accept.AddZone(OMpartialzone);
  OMpartialP.Accept.AddParticle(2212);
  det.AddDevice(OMpartialP);

  Smear::Device OMpartialPt(Smear::kPt, ptformula);
  OMpartialPt.Accept.AddZone(OMpartialzone);
  OMpartialPt.Accept.AddParticle(2212);
  det.AddDevice(OMpartialPt);
  
  Smear::Device OMpartialphi(Smear::kPhi, "0");
  OMpartialphi.Accept.AddZone(OMpartialzone);
  OMpartialphi.Accept.AddParticle(2212);
  det.AddDevice(OMpartialphi);
  
  // Not covered:
  // Tracker material budget X/X0 <~5%
  // Low-Q^2 tagger: -6.9<eta<-5.8: Delta_theta/theta < 1.5%; 10^-6 < Q2 < 10^-2 GeV2
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
