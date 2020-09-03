#include <TMath.h>
#include <TString.h>

#include <iostream>

using namespace std;

static double ThetaFromEta( const double eta );

// Tester for implementation of slide 2 from
// https://indico.bnl.gov/event/8231/contributions/37910/attachments/28335/43607/talk_eic_yr-cal_2020_05.pdf
// For radii, using rough estimates from an ePHENIX  sketch in the white paper
// p. 146 in https://arxiv.org/abs/1212.1701

// Hradius - HCal cylinder radius in mm
// Eradius - EMCCal cylinder radius in mm
// If you need better than integer precision mm,
// you're using the wrong package
int CaloTester( const int Hradius = 1800 /* mm */, const int Eradius = 800 /* mm */){

  // Assumptions:
  // - endcap starts where barrel ends
  // - Granularity is folded into the terms - not used by itself
 
  // Notes
  // - All this assumes single tracks, no overlap
  //   Should warn users to evaluate how far tracks are from each other or
  //   realize that this is too optimistic for overlapping showers.
  // - Don't forget the decimal points in number strings!
  // Also, presumably
  // - hcal is only important for neutral hadrons
  // - ecal is important for photons and for electrons outside the tracker


  // Some options to explore
  // Endcap: is always non-projective and always has
  // X0=10 mm and lambda = 200 mm <- should be pretty safe
  // The table allows for two options in the outer backward ECal
  bool UseWorseOuterBack = true;

  // Barrel:
  // -- Safe: X0=10 mm and lambda = 200 mm
  // -- ePHENIX barrel ecal: X0=7mm
  // Projectivity:
  // -- not at all
  // -- perfect
  // -- imperfect
  // ---- ePhenix:
  //       9 degrees in ECal phi and similar in theta (to be checked)
  //       12 degrees in HCal phi and perfect (for vertex 0) in theta
  //       We know X0 = 7mm for the Ecal
  //
  // We'll provide perfectly projective barrels by default, using
  // conservative X0 (and lambda), and alternatively provide something
  // as close to ePHENIX as we can get it.
  bool MockEPhenixBarrel = false;
  // MockEPhenixBarrel = true;

  bool ProjectiveHCalBarrel = true;
  bool ProjectiveECalBarrel = true;

  TString HCalCylRadius = ""; HCalCylRadius+=Hradius; // cylinder radius in mm
  TString HCalBarrelSigzroot  = "50.";  // in mm, for 1/rootE term
  TString HCalBarrelSigzconst = "30.";  // in mm, for constant term

  // X can stand for lambda or X0
  // lambda  - HCal nuclear interaction length in mm
  // X0      - EMCal radiation length in mm

  // optional, deviation from projectivity
  TString HCalBarrelThetaAngle  = "0.0";  // 12 degrees in radians
  TString HCalBarrelPhiAngle    = "0.0"; 
  TString HCalBarrelX           = "200";      // nuclear interaction length in mm
  if ( MockEPhenixBarrel ) { 
    HCalBarrelThetaAngle  = "0.20944";  // 12 degrees in radians
    HCalBarrelPhiAngle    = "0.0"; 
    HCalBarrelX           = "200";      // nuclear interaction length in mm
  }

  TString HCalEndcapSigRroot   = "50.";  // in mm, for 1/rootE term
  TString HCalEndcapSigRconst  = "30.";  // in mm, for constant term
  TString HCalEndCapX          = "0";  // radiation length in mm

  TString ECalCylRadius = ""; ECalCylRadius+=Eradius; // cylinder radius in mm
  TString ECalBarrelSigzroot  = "3.";  // in mm, for 1/rootE term
  TString ECalBarrelSigzconst = "1.";  // in mm, for constant term

  // optional, deviation from projectivity
  TString ECalBarrelThetaAngle  = "0.0";  // 9 degrees in radians
  TString ECalBarrelPhiAngle    = "0.0";  // 9 degrees in radians
  TString ECalBarrelX           = "10";       // radiation length in mm
  if ( MockEPhenixBarrel ) { 
    ECalBarrelThetaAngle  = "0.15708";  // 9 degrees in radians
    ECalBarrelPhiAngle    = "0.15708";  // 9 degrees in radians
    ECalBarrelX           = "7";       // radiation length in mm
  }
    
  TString ECalEndCapX              = "10";  // radiation length in mm
  // forward
  TString ECalEndcapForwSigRroot   = "3.";  // in mm, for 1/rootE term
  TString ECalEndcapForwSigRconst  = "1.";  // in mm, for constant term

  // inner back
  TString ECalEndcapBack1SigRroot   = "3.";  // in mm, for 1/rootE term
  TString ECalEndcapBack1SigRconst  = "1.";  // in mm, for constant term
  
  // outer back
  TString ECalEndcapBack2SigRroot   = "3.";  // in mm, for 1/rootE term
  if ( UseWorseOuterBack ) ECalEndcapBack2SigRroot   = "6.";  // in mm, for 1/rootE term -- ALTERNATIVE
  TString ECalEndcapBack2SigRconst  = "1.";  // in mm, for constant term

  // =====================
  // Math:
  // =====================
  // R = z * tan (theta), where R = distance from z-axis (usually cylinder radius)
  // <==> theta = atan ( R / z ) == atan2 ( R ,  z ) if needed
  // sigma := sigma_x = sigma_y = sigma_z

  // =====================
  // BARREL: R=const
  // -- sigma(phi) = sigma / R
  // -- sigma(theta) non-projective:
  // sigma(th) = | del theta / del z | sigma
  //           = R / (R^2 + z^2)  * sigma
  //           = sin^2(theta) / R * sigma

  // -- sigma(theta) projective: tower is perpendicular to trajectory, removing one sine:
  // sigma(th) = sin(theta) / R * sigma
  // 
  // 
  // -- for approximately projective cases, replace sigma with sqrt ( sigma^2 + (X*sin( angle ))^2)
  // while this could be folded over into the constant term, the material choices matter quite a bit
  //
  // --> sigma ( phi )   = 1/R            *  sqrt( sigma^2 + (X*sin( angle ))^2)
  // --> sigma ( theta ) = sin(theta) / R *  sqrt( sigma^2 + (X*sin( angle ))^2)

  TString BarrelPhiString   = "1. / cylradius * sigma "; 
  // phi: same for projective, but don't forget to replace sigma in approx. projective cases

  TString BarrelThetaString           = "1. / cylradius * pow( sin(theta ), 2 ) * sigma ";
  TString ProjectiveBarrelThetaString = "1. / cylradius * sin(theta) * sigma ";
  
  // =====================
  // ENDCAP: z=const
  // theta = atan ( sqrt ( x^2 + y^2 ) /z )
  // phi = atan ( y/x )
  // with sigma_x = sigma_y =: sigma

  // -- sigma( phi ) = sqrt ( (del phi / del x )^2 + (del phi / del y )^2 ) * sigma = 1/R * sigma 
  //              = sigma / (z * tan (theta ) ) 

    
  // -- sigma(th) =  sqrt ( (del theta / del x )^2 + (del theta / del y )^2 ) * sigma
  //              = z / (R^2 + z^2)  * sigma
  //              = cos^2(theta) / z * sigma
  // https://www.wolframalpha.com/input/?i=diff+%28+atan+%28+R%2Fz%29%2C+R%29
  // https://www.wolframalpha.com/input/?i=simplify+%28+z+%2F+%28R%5E2%2Bz%5E2%29%2C+R+%3D+z+*+tan%28theta%29%29
  // 
  // However, there is uncertainty in the z direction, for non-projective endcap.
  // again, account for by replacing sigma with sqrt( sigma^2 + (X*sin( angle ))^2)
  // --> sigma(th) = cos^2(theta)/z * sqrt( sigma^2 + (X*sin( angle ))^2)
  // where now angle is theta

  // We don't plan for projective endcaps here.    
  TString EndcapPhiString   = "1. / abs(endcapz) / abs(tan(theta)) * sigma  ";
  TString EndcapThetaString = "1. / endcapz * pow( cos(theta ), 2 ) * sigma ";
  
  // Locations and extents
  // =====================
  double EtaBarrelMin = -1;
  double EtaBarrelMax =  1;

  // convert to theta - inverts min and max
  double ThetaBarrelMin = ThetaFromEta ( EtaBarrelMax );
  double ThetaBarrelMax = ThetaFromEta ( EtaBarrelMin );

  cout << "===============================================================" << endl;
  cout << "Barrel in theta = "
       << ThetaBarrelMin*TMath::RadToDeg() << " -- "
       << ThetaBarrelMax*TMath::RadToDeg() << " degrees" << endl;


  // HCal endcaps are symmetric
  // --------------------------
  double HCalEtaEndcap1Min = -3.5;
  double HCalEtaEndcap1Max = EtaBarrelMin;

  double HCalEtaEndcap2Min = EtaBarrelMax;
  double HCalEtaEndcap2Max = 3.5;

  double HCalThetaEnd1Min = ThetaFromEta ( HCalEtaEndcap1Max );
  double HCalThetaEnd1Max = ThetaFromEta ( HCalEtaEndcap1Min );
  
  double HCalThetaEnd2Min = ThetaFromEta ( HCalEtaEndcap2Max );
  double HCalThetaEnd2Max = ThetaFromEta ( HCalEtaEndcap2Min );

  cout << "===============================================================" << endl;
  cout << "HCAL Endcaps in theta = "
       << HCalThetaEnd2Min*TMath::RadToDeg() << " -- "
       << HCalThetaEnd2Max*TMath::RadToDeg() << " degrees" << endl;
  cout << "            and theta = "
       << HCalThetaEnd1Min*TMath::RadToDeg() << " -- "
       << HCalThetaEnd1Max*TMath::RadToDeg() << " degrees" << endl;

  // calculate endcap distance, assuming R is the same. Use absolute value later.
  double HCalZ1 = Hradius / tan ( HCalThetaEnd1Min );
  double HCalZ2 = Hradius / tan ( HCalThetaEnd2Max );
  cout << "With an HCal cylinder radius of " << Hradius
       << " mm,  endcaps sit at z = " << HCalZ2 << " and  " << HCalZ1 << " mm." << endl;
  cout << "Cross-check: HCal extent is r = "
       << HCalZ2 * tan ( HCalThetaEnd2Min ) << " -- " 
       << HCalZ2 * tan ( HCalThetaEnd2Max ) << " mm " << endl
       << "                        and r = " 
       << HCalZ1 * tan ( HCalThetaEnd1Max ) << " -- " 
       << HCalZ1 * tan ( HCalThetaEnd1Min ) << " mm "
       << endl;

  cout << "===============================================================" << endl;
  cout << endl;
  

  // =============== Repeat for ECal ==================
  
  double ECalEtaBack1Min = -3.5;
  double ECalEtaBack1Max = -2;

  double ECalEtaBack2Min = ECalEtaBack1Max;
  double ECalEtaBack2Max = EtaBarrelMin;
  
  double ECalEtaForwMin = EtaBarrelMax;
  double ECalEtaForwMax = 3.5;

 
  // convert to theta - inverts min and max
  double ECalThetaBack1Min = ThetaFromEta ( ECalEtaBack1Max );
  double ECalThetaBack1Max = ThetaFromEta ( ECalEtaBack1Min );
  
  double ECalThetaBack2Min = ThetaFromEta ( ECalEtaBack2Max );
  double ECalThetaBack2Max = ThetaFromEta ( ECalEtaBack2Min );

  double ECalThetaForwMin = ThetaFromEta ( ECalEtaForwMax );
  double ECalThetaForwMax = ThetaFromEta ( ECalEtaForwMin );

  cout << "===============================================================" << endl;
  cout << "EMCAL Endcaps in theta = "
       << ECalThetaForwMin*TMath::RadToDeg() << " -- "
       << ECalThetaForwMax*TMath::RadToDeg() << " degrees" << endl;
  cout << "             and theta = "    
       << ECalThetaBack2Min*TMath::RadToDeg() << " -- "
       << ECalThetaBack2Max*TMath::RadToDeg() << " degrees" << endl;
  cout << "             and theta = "
       << ECalThetaBack1Min*TMath::RadToDeg() << " -- "
       << ECalThetaBack1Max*TMath::RadToDeg() << " degrees" << endl;
  
  // calculate endcap distances and radii. Use absolute value later.
  double ECalZBack = Eradius / tan ( ECalThetaBack2Min );
  double ECalZForw = Eradius / tan ( ECalThetaForwMax );

  cout << "With an ECal cylinder radius of " << Eradius
       << " mm, endcaps sit at z = " << ECalZForw
       << " and z = " << ECalZBack
       << " mm." << endl;

  cout << "Cross-check: ECal extents are r = "
       << ECalZForw * tan ( ECalThetaForwMin ) << " -- "
       << ECalZForw * tan ( ECalThetaForwMax ) << " mm " << endl
       << "                          and r = " 
       << ECalZBack * tan ( ECalThetaBack1Max ) << " -- "
       << ECalZBack * tan ( ECalThetaBack1Min ) << " mm " << endl
       << "                          and r = " 
       << ECalZBack * tan ( ECalThetaBack2Max ) << " -- "
       << ECalZBack * tan ( ECalThetaBack2Min ) << " mm " << endl;
  cout << "===============================================================" << endl;
  cout << endl;

  // Barrel
  // ------

  // HCal
  // ----
  TString HCalBarrelPhiString = BarrelPhiString;

  if ( ProjectiveHCalBarrel && HCalBarrelPhiAngle.Atof() > 0.0 ){
    HCalBarrelPhiString.ReplaceAll ( "sigma", "sqrt( sigma^2 + (X*sin( angle ))^2)");
    HCalBarrelPhiString.ReplaceAll ("angle", HCalBarrelPhiAngle);
    HCalBarrelPhiString.ReplaceAll ("X", HCalBarrelX);
  }

  HCalBarrelPhiString.ReplaceAll ( "sigma", "(sqrt ( sigzroot*sigzroot / E + sigzconst*sigzconst) )");
  HCalBarrelPhiString.ReplaceAll ("sigzroot", HCalBarrelSigzroot);
  HCalBarrelPhiString.ReplaceAll ("sigzconst", HCalBarrelSigzconst);
  HCalBarrelPhiString.ReplaceAll ("cylradius", HCalCylRadius);

  std::cout << " HCal barrel sigma(phi) described by   " << HCalBarrelPhiString << std::endl;

  TString HCalBarrelThetaString = BarrelThetaString;
  if ( ProjectiveHCalBarrel ) {
    HCalBarrelThetaString = ProjectiveBarrelThetaString;
    if ( HCalBarrelThetaAngle.Atof() > 0.0 ){
      HCalBarrelThetaString.ReplaceAll ( "sigma", "sqrt( sigma^2 + (X*sin( angle ))^2)");
      HCalBarrelThetaString.ReplaceAll ("angle", HCalBarrelThetaAngle);
      HCalBarrelThetaString.ReplaceAll ("X", HCalBarrelX);
    }
  }

  HCalBarrelThetaString.ReplaceAll ( "sigma", "(sqrt ( sigzroot*sigzroot / E + sigzconst*sigzconst) )");
  HCalBarrelThetaString.ReplaceAll ("sigzroot", HCalBarrelSigzroot);
  HCalBarrelThetaString.ReplaceAll ("sigzconst", HCalBarrelSigzconst);
  HCalBarrelThetaString.ReplaceAll ("cylradius", HCalCylRadius);

  std::cout << " HCal barrel sigma(theta) described by " << HCalBarrelThetaString << std::endl;
  std::cout  << std::endl;

  // ECal
  // ----
  TString ECalBarrelPhiString = BarrelPhiString;
  
  if ( ProjectiveECalBarrel && ECalBarrelPhiAngle.Atof() > 0.0 ){
    ECalBarrelPhiString.ReplaceAll ( "sigma", "sqrt( sigma^2 + (X*sin( angle ))^2)");
    ECalBarrelPhiString.ReplaceAll ("angle", ECalBarrelPhiAngle);
    ECalBarrelPhiString.ReplaceAll ("X", ECalBarrelX);
  }

  ECalBarrelPhiString.ReplaceAll ( "sigma", "(sqrt ( sigzroot*sigzroot / E + sigzconst*sigzconst) )");
  ECalBarrelPhiString.ReplaceAll ("sigzroot", ECalBarrelSigzroot);
  ECalBarrelPhiString.ReplaceAll ("sigzconst", ECalBarrelSigzconst);
  ECalBarrelPhiString.ReplaceAll ("cylradius", ECalCylRadius);

  std::cout << " ECal barrel sigma(phi) described by   " << ECalBarrelPhiString << std::endl;
  
  TString ECalBarrelThetaString = BarrelThetaString;
  if ( ProjectiveECalBarrel ) {
    ECalBarrelThetaString = ProjectiveBarrelThetaString;
    if ( ECalBarrelThetaAngle.Atof() > 0.0 ){
      ECalBarrelThetaString.ReplaceAll ( "sigma", "sqrt( sigma^2 + (X*sin( angle ))^2)");
      ECalBarrelThetaString.ReplaceAll ("angle", ECalBarrelThetaAngle);
      ECalBarrelThetaString.ReplaceAll ("X", ECalBarrelX);
    }
  }

  ECalBarrelThetaString.ReplaceAll ( "sigma", "(sqrt ( sigzroot*sigzroot / E + sigzconst*sigzconst) )");
  ECalBarrelThetaString.ReplaceAll ("sigzroot", ECalBarrelSigzroot);
  ECalBarrelThetaString.ReplaceAll ("sigzconst", ECalBarrelSigzconst);
  ECalBarrelThetaString.ReplaceAll ("cylradius", ECalCylRadius);

  ECalBarrelThetaString.ReplaceAll ("cylradius", ECalCylRadius);

  std::cout << " ECal barrel sigma(theta) described by " << ECalBarrelThetaString << std::endl;
  std::cout << std::endl;
  
  cout << "===============================================================" << endl;

  // Endcaps
  // -------
  TString zpos;

  // HCal endcap 1
  zpos =""; zpos += std::abs(HCalZ1);

  TString HCalEndcap1PhiString = EndcapPhiString;
  HCalEndcap1PhiString.ReplaceAll ( "sigma", "sqrt( sigma^2 + (X*sin( theta ))^2)");
  HCalEndcap1PhiString.ReplaceAll ("X", HCalEndCapX);

  HCalEndcap1PhiString.ReplaceAll ( "sigma", "(sqrt ( sigrroot*sigrroot / E + sigrconst*sigrconst) )");

  HCalEndcap1PhiString.ReplaceAll ("sigrroot", HCalEndcapSigRroot );
  HCalEndcap1PhiString.ReplaceAll ("sigrconst", HCalEndcapSigRconst );
  HCalEndcap1PhiString.ReplaceAll ("endcapz", zpos );
  std::cout << " HCal endcap 1 sigma(phi) described by   " << HCalEndcap1PhiString << std::endl;

  TString HCalEndcap1ThetaString = EndcapThetaString;
  HCalEndcap1ThetaString.ReplaceAll ( "sigma", "sqrt( sigma^2 + (X*sin( theta ))^2)");
  HCalEndcap1ThetaString.ReplaceAll ("X", HCalEndCapX);

  HCalEndcap1ThetaString.ReplaceAll ( "sigma", "(sqrt ( sigrroot*sigrroot / E + sigrconst*sigrconst) )");

  HCalEndcap1ThetaString.ReplaceAll ("sigrroot", HCalEndcapSigRroot );
  HCalEndcap1ThetaString.ReplaceAll ("sigrconst", HCalEndcapSigRconst );
  HCalEndcap1ThetaString.ReplaceAll ("endcapz", zpos );
  std::cout << " HCal endcap 1 sigma(theta) described by " << HCalEndcap1ThetaString << std::endl;

  // HCal endcap 2
  zpos ="("; zpos += HCalZ2;  zpos +=")"; 
  TString HCalEndcap2PhiString = EndcapPhiString;
  HCalEndcap2PhiString.ReplaceAll ( "sigma", "sqrt( sigma^2 + (X*sin( theta ))^2)");
  HCalEndcap2PhiString.ReplaceAll ("X", HCalEndCapX);

  HCalEndcap2PhiString.ReplaceAll ( "sigma", "(sqrt ( sigrroot*sigrroot / E + sigrconst*sigrconst) )");
  HCalEndcap2PhiString.ReplaceAll ("sigrroot", HCalEndcapSigRroot );
  HCalEndcap2PhiString.ReplaceAll ("sigrconst", HCalEndcapSigRconst );
  HCalEndcap2PhiString.ReplaceAll ("endcapz", zpos );
  std::cout << " HCal endcap 2 sigma(phi) described by   " << HCalEndcap2PhiString << std::endl;

  TString HCalEndcap2ThetaString = EndcapThetaString;
  HCalEndcap2ThetaString.ReplaceAll ( "sigma", "sqrt( sigma^2 + (X*sin( theta ))^2)");
  HCalEndcap2ThetaString.ReplaceAll ("X", HCalEndCapX);

  HCalEndcap2ThetaString.ReplaceAll ( "sigma", "(sqrt ( sigrroot*sigrroot / E + sigrconst*sigrconst) )");
  HCalEndcap2ThetaString.ReplaceAll ("sigrroot", HCalEndcapSigRroot );
  HCalEndcap2ThetaString.ReplaceAll ("sigrconst", HCalEndcapSigRconst );
  HCalEndcap2ThetaString.ReplaceAll ("endcapz", zpos );
  std::cout << " HCal endcap 2 sigma(theta) described by " << HCalEndcap2ThetaString << std::endl;
  std::cout << std::endl;

  // ----
  
  // ECal forward endcap
  zpos =""; zpos += std::abs(ECalZForw);
    
  TString ECalEndcapForwPhiString = EndcapPhiString;
  ECalEndcapForwPhiString.ReplaceAll ( "sigma", "sqrt( sigma^2 + (X*sin( theta ))^2)");
  ECalEndcapForwPhiString.ReplaceAll ("X", ECalEndCapX);
  ECalEndcapForwPhiString.ReplaceAll ( "sigma", "(sqrt ( sigrroot*sigrroot / E + sigrconst*sigrconst) )");

  ECalEndcapForwPhiString.ReplaceAll ("sigrroot", ECalEndcapForwSigRroot );
  ECalEndcapForwPhiString.ReplaceAll ("sigrconst", ECalEndcapForwSigRconst );
  ECalEndcapForwPhiString.ReplaceAll ("endcapz", zpos );
  std::cout << " ECal Forward endcap sigma(phi) described by   " << ECalEndcapForwPhiString << std::endl;

  TString ECalEndcapForwThetaString = EndcapThetaString;
  ECalEndcapForwThetaString.ReplaceAll ( "sigma", "sqrt( sigma^2 + (X*sin( theta ))^2)");
  ECalEndcapForwThetaString.ReplaceAll ("X", ECalEndCapX);

  ECalEndcapForwThetaString.ReplaceAll ( "sigma", "(sqrt ( sigrroot*sigrroot / E + sigrconst*sigrconst) )");

  ECalEndcapForwThetaString.ReplaceAll ("sigrroot", ECalEndcapForwSigRroot );
  ECalEndcapForwThetaString.ReplaceAll ("sigrconst", ECalEndcapForwSigRconst );
  ECalEndcapForwThetaString.ReplaceAll ("endcapz", zpos );
  std::cout << " ECal Forward endcap sigma(theta) described by " << ECalEndcapForwThetaString << std::endl;
  std::cout << std::endl;

  
  // ECal inner backward endcap
  zpos =""; zpos += std::abs(ECalZBack);
  TString ECalEndcapBack1PhiString = EndcapPhiString;
  ECalEndcapBack1PhiString.ReplaceAll ( "sigma", "sqrt( sigma^2 + (X*sin( theta ))^2)");
  ECalEndcapBack1PhiString.ReplaceAll ("X", ECalEndCapX);
  ECalEndcapBack1PhiString.ReplaceAll ( "sigma", "(sqrt ( sigrroot*sigrroot / E + sigrconst*sigrconst) )");
  ECalEndcapBack1PhiString.ReplaceAll ("sigrroot", ECalEndcapBack1SigRroot );
  ECalEndcapBack1PhiString.ReplaceAll ("sigrconst", ECalEndcapBack1SigRconst );
  ECalEndcapBack1PhiString.ReplaceAll ("endcapz", zpos );
  std::cout << " ECal inner backward endcap sigma(phi) described by   " << ECalEndcapBack1PhiString << std::endl;

  TString ECalEndcapBack1ThetaString = EndcapThetaString;
  ECalEndcapBack1ThetaString.ReplaceAll ( "sigma", "sqrt( sigma^2 + (X*sin( theta ))^2)");
  ECalEndcapBack1ThetaString.ReplaceAll ("X", ECalEndCapX);

  ECalEndcapBack1ThetaString.ReplaceAll ( "sigma", "(sqrt ( sigrroot*sigrroot / E + sigrconst*sigrconst) )");

  ECalEndcapBack1ThetaString.ReplaceAll ("sigrroot", ECalEndcapBack1SigRroot );
  ECalEndcapBack1ThetaString.ReplaceAll ("sigrconst", ECalEndcapBack1SigRconst );
  ECalEndcapBack1ThetaString.ReplaceAll ("endcapz", zpos );
  std::cout << " ECal inner backward endcap sigma(theta) described by " << ECalEndcapBack1ThetaString << std::endl;

  // ECal outer backward endcap
  TString ECalEndcapBack2PhiString = EndcapPhiString;
  ECalEndcapBack2PhiString.ReplaceAll ( "sigma", "sqrt( sigma^2 + (X*sin( theta ))^2)");
  ECalEndcapBack2PhiString.ReplaceAll ("X", ECalEndCapX);
  ECalEndcapBack2PhiString.ReplaceAll ( "sigma", "(sqrt ( sigrroot*sigrroot / E + sigrconst*sigrconst) )");

  ECalEndcapBack2PhiString.ReplaceAll ("sigrroot", ECalEndcapBack2SigRroot );
  ECalEndcapBack2PhiString.ReplaceAll ("sigrconst", ECalEndcapBack2SigRconst );
  ECalEndcapBack2PhiString.ReplaceAll ("endcapz", zpos );
  std::cout << " ECal outer backward endcap sigma(phi) described by   " << ECalEndcapBack2PhiString << std::endl;

  TString ECalEndcapBack2ThetaString = EndcapThetaString;
  ECalEndcapBack2ThetaString.ReplaceAll ( "sigma", "sqrt( sigma^2 + (X*sin( theta ))^2)");
  ECalEndcapBack2ThetaString.ReplaceAll ("X", ECalEndCapX);

  ECalEndcapBack2ThetaString.ReplaceAll ( "sigma", "(sqrt ( sigrroot*sigrroot / E + sigrconst*sigrconst) )");

  ECalEndcapBack2ThetaString.ReplaceAll ("sigrroot", ECalEndcapBack2SigRroot );
  ECalEndcapBack2ThetaString.ReplaceAll ("sigrconst", ECalEndcapBack2SigRconst );
  ECalEndcapBack2ThetaString.ReplaceAll ("endcapz", zpos );
  std::cout << " ECal outer backward endcap sigma(theta) described by " << ECalEndcapBack2ThetaString << std::endl;

  std::cout << std::endl;
  

  // =============== Plot ==================
  gStyle->SetOptStat(0);
  double lw = 5;

  TString formulastring;
  // Guide the eye
  TLine l;
  l.SetLineStyle(7);
  l.SetLineWidth(1);
  l.SetLineColor(kGray+1);
  const double sigref = 0.004;
  vector <TString> energies = {"1", "10", "100" };

  TH1D* dummy = new TH1D("dummy","HCal #theta Resolution;#theta [rad];#sigma(#theta) [rad]", 100, 0, TMath::Pi());
  // Lazily guessing the maximum
  formulastring = HCalBarrelThetaString;  formulastring.ReplaceAll("theta", "x");
  formulastring.ReplaceAll("E", "1");
  TF1* maxguess = new TF1("maxguess",formulastring, ThetaBarrelMin, ThetaBarrelMax);
  dummy->SetAxisRange(0,1.8 * maxguess->Eval(TMath::Pi()/2), "y");

  // Plot for 1, 10, and 100 GeV
  new TCanvas;
  dummy->DrawClone("");
  l.DrawLine(0, sigref, TMath::Pi(), sigref);
  auto leg = new TLegend ( .35, .65, .65, .88, "Neutral hadrons");
  
  for ( auto energy : energies ){
    int color = kRed;
    if ( energy == "10" ) color = kGreen+1;
    if ( energy == "100" ) color = kMagenta+1;
    
    formulastring = HCalBarrelThetaString;  formulastring.ReplaceAll("theta", "x");
    formulastring.ReplaceAll("E", energy);
    TF1* fB = new TF1("fB_" + energy,formulastring, ThetaBarrelMin, ThetaBarrelMax);
    fB->SetLineColor( color );    fB->SetLineWidth( lw );
    fB->Draw("same");
  
    formulastring = HCalEndcap1ThetaString;  formulastring.ReplaceAll("theta", "x");
    formulastring.ReplaceAll("E", energy);
    TF1* f1 = new TF1("f1_" + energy,formulastring, HCalThetaEnd1Min, HCalThetaEnd1Max);
    f1->SetLineColor( color );    f1->SetLineWidth( lw );
    f1->Draw("same");

    formulastring = HCalEndcap2ThetaString;  formulastring.ReplaceAll("theta", "x");
    formulastring.ReplaceAll("E", energy);
    TF1* f2 = new TF1("f2_" + energy,formulastring, HCalThetaEnd2Min, HCalThetaEnd2Max);
    f2->SetLineColor( color );    f2->SetLineWidth( lw );
    f2->Draw("same");

	
    leg->AddEntry( f1, "Energy = " + energy + " GeV");
  }
  leg->Draw();
  if (MockEPhenixBarrel) gPad->SaveAs("HCal_ePhenix.png");
  else gPad->SaveAs("HCal.png");

  // Draw hcal phi resolution
  
  new TCanvas;
  // formulastring = HCalBarrelPhiString;  formulastring.ReplaceAll("theta", "x");
  // formulastring.ReplaceAll("E", "1");
  // TF1* maxguess = new TF1("maxguess",formulastring, ThetaBarrelMin, ThetaBarrelMax);
  // dummy->SetAxisRange(0,1.8 * maxguess->Eval(TMath::Pi()/2), "y");
  TH1D* dummyphi = new TH1D("dummyphi","HCal #phi Resolution ;#theta [rad];#sigma(#phi) [rad]", 100, 0, TMath::Pi());
  dummyphi->SetAxisRange(0,0.5, "y");
  dummyphi->DrawClone("");
  l.DrawLine(0, sigref, TMath::Pi(), sigref); 
  leg = new TLegend ( .35, .65, .65, .88, "Neutral hadrons");
 
  for ( auto energy : energies ){
    int color = kRed;
    if ( energy == "10" ) color = kGreen+1;
    if ( energy == "100" ) color = kMagenta+1;
    
    // phi
    formulastring = HCalBarrelPhiString;  formulastring.ReplaceAll("theta", "x");
    formulastring.ReplaceAll("E", energy);
    TF1* fPhiB = new TF1("fPhiB_" + energy,formulastring, ThetaBarrelMin, ThetaBarrelMax);
    fPhiB->SetLineColor( color );    fPhiB->SetLineWidth( lw );
    fPhiB->Draw("same");

    formulastring = HCalEndcap1PhiString;  formulastring.ReplaceAll("theta", "x");
    formulastring.ReplaceAll("E", energy);
    TF1* fPhi1 = new TF1("fPhi1_" + energy,formulastring, HCalThetaEnd1Min, HCalThetaEnd1Max);
    fPhi1->SetLineColor( color );    fPhi1->SetLineWidth( lw );
    fPhi1->Draw("same");

    formulastring = HCalEndcap2PhiString;  formulastring.ReplaceAll("theta", "x");
    formulastring.ReplaceAll("E", energy);
    TF1* fPhi2 = new TF1("f2Phi_" + energy,formulastring, HCalThetaEnd2Min, HCalThetaEnd2Max);
    fPhi2->SetLineColor( color );    fPhi2->SetLineWidth( lw );
    fPhi2->Draw("same");

	
    leg->AddEntry( fPhiB, "Energy = " + energy + " GeV");
  }
  leg->Draw();
  if (MockEPhenixBarrel) gPad->SaveAs("HCal_Phi_ePhenix.png");
  else gPad->SaveAs("HCal_Phi.png");


  
  // =============== Plot ECal ==================

  TH1D* dummy2 = new TH1D("dummy2","ECal #theta Resolution;#theta [rad];#sigma(#theta) [rad]", 100, 0,TMath::Pi());

  // Lazily guessing the maximum
  formulastring = ECalBarrelThetaString;  formulastring.ReplaceAll("theta", "x");
  formulastring.ReplaceAll("E", "1");
  maxguess = new TF1("maxguess2",formulastring, ThetaBarrelMin, ThetaBarrelMax);
  dummy2->SetAxisRange(0,1.7 * maxguess->Eval(TMath::Pi()/2), "y");
  
  // Plot for 1, 10, and 100 GeV
  new TCanvas;
  dummy2->DrawClone("");
  l.DrawLine(0, sigref, TMath::Pi(), sigref);
  auto leg1 = new TLegend ( .35, .65, .65, .88, "Electrons & #gamma");
  
  for ( auto energy : energies ){
    int color = kRed;
    if ( energy == "10" ) color = kGreen+1;
    if ( energy == "100" ) color = kMagenta+1;
    
    formulastring = ECalBarrelThetaString;  formulastring.ReplaceAll("theta", "x");
    formulastring.ReplaceAll("E", energy);
    TF1* fB = new TF1("fB_" + energy,formulastring, ThetaBarrelMin, ThetaBarrelMax);
    fB->SetLineColor( color );    fB->SetLineWidth( lw );
    fB->Draw("same");
  
    formulastring = ECalEndcapBack2ThetaString;  formulastring.ReplaceAll("theta", "x");
    formulastring.ReplaceAll("E", energy);
    TF1* fB2 = new TF1("fB2_" + energy,formulastring, ECalThetaBack2Min, ECalThetaBack2Max);
    fB2->SetLineColor( color );    fB2->SetLineWidth( lw );
    fB2->Draw("same");

    formulastring = ECalEndcapBack1ThetaString;  formulastring.ReplaceAll("theta", "x");
    formulastring.ReplaceAll("E", energy);
    TF1* fB1 = new TF1("fB1_" + energy,formulastring, ECalThetaBack1Min, ECalThetaBack1Max);
    fB1->SetLineColor( color );    fB1->SetLineWidth( lw );
    fB1->Draw("same");

    formulastring = ECalEndcapForwThetaString;  formulastring.ReplaceAll("theta", "x");
    formulastring.ReplaceAll("E", energy);
    TF1* fF = new TF1("fF_" + energy,formulastring, ECalThetaForwMin, ECalThetaForwMax);
    fF->SetLineColor( color );    fF->SetLineWidth( lw );
    fF->Draw("same");
        
    leg1->AddEntry( fB, "Energy = " + energy + " GeV");
  }
  leg1->Draw();

  if (MockEPhenixBarrel) gPad->SaveAs("ECal_ePhenix.png");
  else gPad->SaveAs("ECal.png");
  
  TH1D* dummy2phi = new TH1D("dummy2phi","ECal #phi Resolution;#theta [rad];#sigma(#phi) [rad]", 100, 0,TMath::Pi());
  dummy2phi->SetAxisRange(0,0.06, "y");
  new TCanvas;
  dummy2phi->DrawClone("");
  l.DrawLine(0, sigref, TMath::Pi(), sigref);
  leg1 = new TLegend ( .35, .65, .65, .88, "Electrons & #gamma");
  
  for ( auto energy : energies ){
    int color = kRed;
    if ( energy == "10" ) color = kGreen+1;
    if ( energy == "100" ) color = kMagenta+1;
        
    formulastring = ECalBarrelPhiString;  formulastring.ReplaceAll("theta", "x");
    formulastring.ReplaceAll("E", energy);
    TF1* fPhiB = new TF1("fPhiB_" + energy,formulastring, ThetaBarrelMin, ThetaBarrelMax);
    fPhiB->SetLineColor( color );    fPhiB->SetLineWidth( lw );
    fPhiB->Draw("same");

    formulastring = ECalEndcapForwPhiString;  formulastring.ReplaceAll("theta", "x");
    formulastring.ReplaceAll("E", energy);
    TF1* fPhiF = new TF1("fPhiF_" + energy,formulastring, ECalThetaForwMin, ECalThetaForwMax);
    fPhiF->SetLineColor( color );    fPhiF->SetLineWidth( lw );
    fPhiF->Draw("same");
    
    formulastring = ECalEndcapBack1PhiString;  formulastring.ReplaceAll("theta", "x");
    formulastring.ReplaceAll("E", energy);
    TF1* fPhiB1 = new TF1("fPhiB1_" + energy,formulastring, ECalThetaBack1Min, ECalThetaBack1Max);
    fPhiB1->SetLineColor( color );    fPhiB1->SetLineWidth( lw );
    fPhiB1->Draw("same");

    formulastring = ECalEndcapBack2PhiString;  formulastring.ReplaceAll("theta", "x");
    formulastring.ReplaceAll("E", energy);
    TF1* fPhiB2 = new TF1("fPhiB2_" + energy,formulastring, ECalThetaBack2Min, ECalThetaBack2Max);
    fPhiB2->SetLineColor( color );    fPhiB2->SetLineWidth( lw );
    fPhiB2->Draw("same");

    leg1->AddEntry( fPhiB, "Energy = " + energy + " GeV");
  }
  leg1->Draw();

  if (MockEPhenixBarrel) gPad->SaveAs("ECal_Phi_ePhenix.png");
  else gPad->SaveAs("ECal_Phi.png");

  return 0;
}
// -------------------------------------------------------------------
double ThetaFromEta( const double eta ) {
  if ( !std::isnan(eta) && !std::isinf(eta)   ) {
    return 2.0 * atan( exp( -eta ));
  }
  throw std::runtime_error("ThetaFromEta called with NaN or Inf");
  return -1;
}
