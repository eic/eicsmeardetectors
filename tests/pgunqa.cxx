// Create similar QA plots as qaplots.cxx but using a particle gun
// This implementation also shows how to smear an individual particle

// reusing the same header since it mostly sets up
// default parameters
#include "qaplots.hh"

#include <eicsmear/functions.h>
#include <eicsmear/smear/functions.h>

#include "eicsmear/erhic/Particle.h"
#include "eicsmear/smear/Smear.h"
#include "eicsmear/smear/ParticleMCS.h"
#include "eicsmear/smear/Detector.h"

#include "eicsmeardetectors.hh"

// Note: The remaining includes are not necessary for eic-smear usage
#include <TSystem.h>
#include <TRandom.h>
#include <TDatabasePDG.h>
#include <TFile.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TProfile.h>
#include <TLegend.h>
#include <TText.h>
#include <TStyle.h>

#include <string>
#include <iomanip>
#include <cctype>
#include <exception>
#include <vector>
#include <map>

// Convenience only
using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::vector;
using std::map;

// some helpers
static void initializepidqabook(const qaparameters& qapars, map<int,pidqacollection>& qabook );
static void FillParticleQA( map<int,pidqacollection>& qabook, const Particle* const inParticle, const Smear::ParticleMCS* const inParticleS  );
static void PlotQA ( const qaparameters& qapars, map<int,pidqacollection>& qabook );

static double EtaFromTheta( const double theta );
static double ThetaFromEta( const double eta );

int main(int argc, char* argv[]){
  
  // Parse arguments
  // ---------------
  // Defaults are set in qaplots.h
  qaparameters qapars = ParseArguments ( argc, argv );

  // First try to instantiate the detector 
  // to avoid pointlessly transforming if that doesn't work
  // ------------------------------------------------------
  Smear::Detector detector;
  if ( qapars.beam_mom_nn < 0 ) detector = BuildByName(qapars.detstring);
  else                          detector = BuildByName(qapars.detstring, qapars.beam_mom_nn);
  
  if ( detector.GetNDevices() == 0 ) {
    cerr << "Detector sepcified as " << qapars.detstring
	 << " not recognized or empty." << endl;
    return -1;
  }
  
  // Open histo root file and book histograms
  // ----------------------------------------
  TFile * outfile = new TFile ( qapars.outfilebase + ".pgun." + qapars.detstring + ".root", "RECREATE");

  // A collection of histos and maybe other info for every pid
  // pidqacollection is defined in the header
  map<int,pidqacollection> qabook;
  // By default, use the standard particles
  if ( qapars.pids.size() == 0 ) qapars.pids = { 11, 211, 321, 2212, 2112 }; // e, pi, K, p, n
  initializepidqabook ( qapars, qabook );

  // the initialization routine determines ranges,
  // use those for range limits
  TH2D* h = qabook[ qapars.pids[0] ].dPhi_p;
  // reject super low momenta
  double pmin = std::max ( h->GetXaxis()->GetXmin(), 0.1);
  double pmax = h->GetXaxis()->GetXmax();
  
  // push theta away a litle bit from 0, pi
  double eps = 1e-5;
  h = qabook[ qapars.pids[0] ].DelP_th;
  double thmin = std::max ( h->GetXaxis()->GetXmin(), eps );
  double thmax = std::min ( h->GetXaxis()->GetXmax(), TMath::Pi() - eps);
  // want to generate in eta though
  double etamin =  EtaFromTheta(thmin);
  double etamax =  EtaFromTheta(thmax);
  
  double phimin = -TMath::Pi();
  double phimax = TMath::Pi();

  // helper for mass, charge
  TDatabasePDG* db = TDatabasePDG::Instance();
    
  // -------------------
  // Loop over Particles
  // -------------------
  Particle* inParticle=nullptr;
  Smear::ParticleMCS* inParticleS=nullptr;
  for(int np=0; np<qapars.nparticles; np++){
    if ( !(np%20000) ) cout << "Generated " << np*qapars.pids.size() << " particles" << endl;

    for ( auto pid : qapars.pids ){
      // helper
      TParticlePDG* pdg_particle = db->GetParticle(pid);

      // random charge
      if ( std::abs(pdg_particle->Charge()) <1e-1){
	if ( gRandom->Integer(2) == 1) pid = -pid;
      }
      // Create not smeared particle
      // flat in phi, eta, pt (!)

      // double eta   = gRandom->Uniform (etamin,etamax);
      double theta   = gRandom->Uniform (thmin,thmax);
      double eta     = EtaFromTheta(theta);

      double phi   = gRandom->Uniform (phimin,phimax);
      // double mom   = gRandom->Uniform (pmin,pmax);
      // double pt    = mom * sin (theta);
      double pt=-1;
      double mom=-1;
      do{
	pt = gRandom->Uniform (0,pmax);
	mom = pt / sin ( ThetaFromEta(eta) );
	// cout << pt << "  " << mom << endl;
      } while ( mom < pmin || mom > pmax);
      
      TLorentzVector input_vect;
      input_vect.SetPtEtaPhiM ( pt, eta , phi, pdg_particle->Mass());

      if ( inParticle ) delete inParticle;
      inParticle = new Particle();
      inParticle->SetId(pid);    // PDG particle code
      inParticle->SetIndex(0);   // Particle index in event
      inParticle->SetStatus(1);  // Particle status code: like in PYTHIA, Beagle, etc
      inParticle->Set4Vector(input_vect);
      
      // Smear the particle
      if ( inParticleS ) delete inParticleS;
      inParticleS = detector.Smear( *inParticle);

      // Particle was not smeared
      if( !inParticleS ) continue; 
      
      // ----------------
      // Particle-wise QA
      // ----------------
      // Following function just contains many statements of the form
      // coll.dEta_p->Fill(inParticle->GetP(), inParticle->GetEta() - inParticleS->GetEta());
      FillParticleQA( qabook, inParticle, inParticleS );
    }
  }
  if ( inParticle ) delete inParticle;
  if ( inParticleS ) delete inParticleS;
  
  // NOTE: The remainder of this long file is tedious and explicit creation and filling of histograms
  // The Fill*QA functions demonstrate how to access values in detail
  PlotQA( qapars, qabook );
  outfile->Write();
  
  return 0;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

void FillParticleQA( map<int,pidqacollection>& qabook, const Particle* const inParticle, const Smear::ParticleMCS* const inParticleS  ){
    
  // If any component is smeared, all others are either smeared or 0 (meaning "not detected").
  
  // Fill histograms
  for ( auto& pidcoll : qabook ){
    auto& pid = pidcoll.first;
    auto& coll = pidcoll.second;
    if ( pid==0 || inParticle->GetPdgCode() == pid ){

      if ( inParticleS->IsPSmeared() ){
	auto delP = (inParticle->GetP() - inParticleS->GetP()) / inParticle->GetP();
	coll.DelP_th->Fill(inParticle->GetTheta(), delP);
	coll.DelP_eta->Fill(inParticle->GetEta(), delP);
      }

      if ( inParticleS->IsESmeared() ){
	auto delE = (inParticle->GetE() - inParticleS->GetE()) / inParticle->GetE();
	coll.DelE_E->Fill(inParticle->GetE(), delE);
	coll.DelE_th->Fill(inParticle->GetTheta(), delE);
	coll.DelE_eta->Fill(inParticle->GetEta(), delE);
      }
	
      if ( inParticleS->IsThetaSmeared() ){
	coll.dTh_p->Fill(inParticle->GetP(), inParticle->GetTheta() - inParticleS->GetTheta());
	coll.dEta_p->Fill(inParticle->GetP(), inParticle->GetEta() - inParticleS->GetEta());
      }

      if ( inParticleS->IsPhiSmeared() ){
	coll.dPhi_p->Fill(inParticle->GetP(), inParticle->GetPhi() - inParticleS->GetPhi());
      }
    }
  }
}

// -----------------------------------------------------------------------------

qaparameters ParseArguments ( int argc, char* argv[] ){
  vector<string> arguments(argv + 1, argv + argc);
  bool argsokay=true;
  qaparameters qapars;
  try{
    for ( auto parg = arguments.begin() ; parg!=arguments.end() ; ++parg){
      string arg=*parg;
      if ( arg == "-h" ){
	argsokay=false;
	break; 
      } else if ( arg == "-o" ){
	if (++parg == arguments.end() ){ argsokay=false; break; }
	qapars.outfilebase=*parg;
      } else if ( arg == "-N" ){      
	if (++parg==arguments.end() ){ argsokay=false; break; }
	qapars.nparticles=std::stoi(parg->data());
      } else if ( arg == "-addpid" ){
	if ( ++parg == arguments.end() ){ argsokay=false; break; }
	qapars.pids.push_back(std::stoi(parg->data()));
      } else if ( arg == "-det" ){
	if (++parg == arguments.end() ){ argsokay=false; break; }
	qapars.detstring=*parg;
	for (auto & c: qapars.detstring) c = toupper(c);
	if ( TString(qapars.detstring).Contains("MATRIX") && TString(qapars.detstring).Contains("FF")){
	  if (++parg == arguments.end() ){ argsokay=false; break; }
	  qapars.beam_mom_nn = std::stoi(parg->data());
	}
	if ( TString(qapars.detstring).Contains("CORE") && !TString(qapars.detstring).Contains("B")){
	  if (++parg == arguments.end() ){ argsokay=false; break; }
	  qapars.beam_mom_nn = std::stod(parg->data());
	}	
      } else {
	argsokay=false;
	break;
      }
    }
  } catch ( const std::exception& e){
    cerr << "Caught exception during argument parsing: "
	 << e.what() << endl;
    argsokay=false; 
  }  
  
  if ( !argsokay ) {
    cerr << "usage: " << argv[0] << endl
	 << " [-o OutFileBase] (extension will be added)"  << endl
      	 << " [-N #particles per pid]" << endl
      	 << " [-addpid pid] (can be called multiple times)" << endl
	 << " [-det detstring] matrix, matrixff [beam_mom_nn], handbook, perfect, beast, ephenix, zeus, jleic (capitalization does not matter.)" << endl
	 << endl;
    throw std::runtime_error("Not a valid list of options");
  }
  for (auto & c: qapars.detstring) c = toupper(c);

  return qapars;
}

// ---------------------------------------------------------------
void initializepidqabook(const qaparameters& qapars, map<int,pidqacollection>& qabook ){
  gStyle->SetHistLineColor(kRed); // for Profiles

  TString s;
  float pmin = 0;
  float pmax = 20;
  int pbins = 80;

  float dpmin = -0.15;
  float dpmax = 0.15;
  int dpbins = 100;
  
  float emin = 0;
  float emax = 20;
  int ebins = 80;

  float demin = -1.5;
  float demax = 1.5;
  int debins = 100;

  float thmin = 0;
  float thmax = TMath::Pi();
  int thbins = 64;
  
  float dthmin = -0.1;
  float dthmax = 0.1;
  int dthbins = 100;

  float etamin = -5;
  float etamax = 5;
  int etabins = 100;
  
  float detamin = -0.1;
  float detamax = 0.1;
  int detabins = 100;

  double phimin = -TMath::Pi();
  double phimax = TMath::Pi();
  int phibins = 64;
  
  float dphimin = -0.1;
  float dphimax = 0.1;
  int dphibins = 100;

  for ( auto pid : qapars.pids ){
    pid = abs ( pid ); // ignoring charge

    s = qapars.detstring + "_DelE_E_"; s += pid;
    qabook[pid].DelE_E = new TH2D( s,s+";E;#Delta E/E", ebins, emin, emax, debins, demin, demax);

    s = qapars.detstring + "_dPhi_p_"; s += pid;
    qabook[pid].dPhi_p = new TH2D( s,s+";p;#Delta#phi", pbins, pmin, pmax, dphibins, dphimin, dphimax );    

    s = qapars.detstring + "_DelP_th_"; s += pid;
    qabook[pid].DelP_th = new TH2D( s,s+";#theta;#Delta p/p", thbins, thmin, thmax, dpbins, dpmin, dpmax);

    s = qapars.detstring + "_DelE_th_"; s += pid;
    qabook[pid].DelE_th = new TH2D( s,s+";#theta;#Delta E/E", thbins, thmin, thmax, debins, demin, demax);
    
    s = qapars.detstring + "_dTh_p_"; s += pid;
    qabook[pid].dTh_p = new TH2D( s,s+";p;#Delta#theta", pbins, pmin, pmax, dthbins, dthmin, dthmax );

    s = qapars.detstring + "_DelP_eta_"; s += pid;
    qabook[pid].DelP_eta = new TH2D( s,s+";#eta;#Delta p/p", etabins, etamin, etamax, dpbins, dpmin, dpmax);

    s = qapars.detstring + "_DelE_eta_"; s += pid;
    qabook[pid].DelE_eta = new TH2D( s,s+";#eta;#Delta E/E", etabins, etamin, etamax, debins, demin, demax);
    
    s = qapars.detstring + "_dEta_p_"; s += pid;
    qabook[pid].dEta_p = new TH2D( s,s+";p;#Delta#eta", pbins, pmin, pmax, detabins, detamin, detamax );

  }

}

// ---------------------------------------------------------------
void PlotQA ( const qaparameters& qapars, map<int,pidqacollection>& qabook ){

  // Stat  position and size
  // -----------------------
  gStyle->SetStatX(0.55);
  gStyle->SetStatW(0.15);
  gStyle->SetStatY(0.9);
  gStyle->SetStatH(0.15);


  // Suppress obnoxious flood of
  // "Current canvas added to pdf file"
  gErrorIgnoreLevel = kWarning;

  // prep a pdf collection
  new TCanvas;
  auto pdfname = qapars.outfilebase + ".pgun." + qapars.detstring + ".pdf";
  gPad->SaveAs( pdfname + "[" );

  // particle QA
  // -----------
  for ( auto& pidcoll : qabook ){
    auto& pid = pidcoll.first;
    auto& coll = pidcoll.second;

    // option "s" in Profile shows rms
    
    coll.DelP_th->Draw("colz");
    gPad->SaveAs( pdfname );

    coll.DelP_eta->Draw("colz");
    coll.DelP_eta->ProfileX("_px",1,-1,"s")->Draw("same");
    gPad->SaveAs( pdfname );

    gStyle->SetStatX(0.25); // reposition stat box
    coll.DelE_th->Draw("colz");
    gPad->SaveAs( pdfname );
    
    coll.DelE_eta->Draw("colz");
    coll.DelE_eta->ProfileX("_px",1,-1,"s")->Draw("same");
    gPad->SaveAs( pdfname );
    gStyle->SetStatX(0.55); // reposition stat box

    coll.DelE_E->Draw("colz");
    coll.DelE_E->ProfileX("_px",1,-1,"s")->Draw("same");
    gPad->SaveAs( pdfname );

    coll.dTh_p->Draw("colz");
    gPad->SaveAs( pdfname );
    
    coll.dEta_p->Draw("colz");
    coll.dEta_p->ProfileX("_px",1,-1,"s")->Draw("same");
    gPad->SaveAs( pdfname );

    coll.dPhi_p->Draw("colz");
    coll.dPhi_p->ProfileX("_px",1,-1,"s")->Draw("same");
    gPad->SaveAs( pdfname );
  }

  // return to standard warning level
  gErrorIgnoreLevel = kInfo;

  // close the pdf collection
  gPad->SaveAs( pdfname + "]" );
}
// ---------------------------------------------------------------
double EtaFromTheta( const double theta ) {
  return log ( tan ( 0.5*theta ));
}; 
// -------------------------------------------------------------------
double ThetaFromEta( const double eta ) {
  if ( !std::isnan(eta) && !std::isinf(eta)   ) {
    return 2.0 * atan( exp( -eta ));
  }
  throw std::runtime_error("ThetaFromEta called with NaN or Inf");
  return -1;
}
