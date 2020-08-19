//-----------------------
//Simple Detector: perfect resolution
//Acceptance from [-15,15] in eta
// Original author: Barak Schmookler
//-----------------------

#include "eicsmear/erhic/VirtualParticle.h"
#include "eicsmear/smear/Acceptance.h"
#include "eicsmear/smear/Device.h"
#include "eicsmear/smear/Detector.h"
#include "eicsmear/smear/Smearer.h"
#include "eicsmear/smear/ParticleMCS.h"
#include "eicsmear/smear/PerfectID.h"
#include "piddetectors/TofBarrelSmearer.h"


//Convert pseudorapidity (eta) to polar angle (theta) in radians.
//Make use of TLorentzVector to do eta-to-theta conversion.
static double etaToTheta(const double eta) {
  TLorentzVector v;
  v.SetPtEtaPhiM(1., eta, 0., 0.);
  return v.Theta();
}

//Build Detector
Smear::Detector BuildWithTof(){

  //Create Devices
  Smear::Device energy(Smear::kE,"0"); //"0":Perfect Resolution
  Smear::Device momentum(Smear::kP,"0");
  Smear::Device theta(Smear::kTheta,"0");
  Smear::Device phi(Smear::kPhi,"0");

  //Detector Acceptance
  Smear::Acceptance::Zone acceptall(etaToTheta(15.),etaToTheta(-15.));

  energy.Accept.AddZone(acceptall);
  momentum.Accept.AddZone(acceptall);
  theta.Accept.AddZone(acceptall);
  phi.Accept.AddZone(acceptall);

  Smear::TofBarrelSmearer tofBarrel(100, -1.0, 1.0, 10);
  tofBarrel.Accept.AddZone( acceptall ); // cuts are done by the detector!
  tofBarrel.Accept.SetCharge(Smear::kCharged);
  
  // //PID performance is unparameterised as of now
  // Smear::PerfectID pid;
  // pid.Accept.AddZone(acceptall);

  //Create the detector and add devices
  Smear::Detector det;
  det.AddDevice(energy);
  det.AddDevice(momentum);
  det.AddDevice(theta);
  det.AddDevice(phi);
  det.AddDevice(tofBarrel);
  // det.AddDevice(pid);
  det.SetEventKinematicsCalculator("NM JB DA");

  return det;

}
