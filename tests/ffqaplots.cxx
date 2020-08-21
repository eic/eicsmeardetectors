// Based on qaplots.cxx with a more focused selection
// of far forward plots
// It's assumed you've already created the input root file
// for example with qaplots. Otherwise uncomment the
// BuildTree command below



#include "ffqaplots.hh"

#include <eicsmear/functions.h>
#include <eicsmear/smear/functions.h>

#include "eicsmear/erhic/EventBase.h"
#include "eicsmear/erhic/EventPythia.h"
#include "eicsmear/erhic/Particle.h"
#include "eicsmear/smear/EventSmear.h"
#include "eicsmear/smear/EventS.h"
#include "eicsmear/smear/Smear.h"
#include "eicsmear/smear/ParticleMCS.h"
#include "eicsmear/smear/Detector.h"

#include "eicsmeardetectors.hh"

// Note: The remaining includes are not necessary for eic-smear usage
#include <TSystem.h>
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
static const TString getrootname(const qaparameters& qapars );
static void initializepidqabook(const qaparameters& qapars, map<int,pidqacollection>& qabook );
static void initializeeventqa(const qaparameters& qapars, eventqacollection& eventqa );
static void FillEventQA( eventqacollection& eventqa , const erhic::EventMC* const inEvent, const Smear::Event* const inEventS );
static void FillParticleQA( map<int,pidqacollection>& qabook,
			    const Particle* const inParticle,
			    const Smear::ParticleMCS* const inParticleS,
			    const qaparameters& qapars);
static void PlotQA ( const qaparameters& qapars, eventqacollection& eventqa, map<int,pidqacollection>& qabook );
  
int main(int argc, char* argv[]){
  
  // Parse arguments
  // ---------------
  // Defaults are set in qaplots.h
  qaparameters qapars = ParseArguments ( argc, argv );

  // Set up output name
  TString rootname = getrootname(qapars);

  // First try to instantiate the detector 
  // to avoid pointlessly transforming if that doesn't work
  // ------------------------------------------------------
  Smear::Detector detector;
  if ( qapars.beam_mom_nn < 0 ) detector = BuildByName(qapars.detstring);
  else                          detector = BuildByName(qapars.detstring, qapars.beam_mom_nn);
  if ( detector.GetNDevices() == 0 ) {
    cerr << "Detector sepcified as " << qapars.detstring
	 << " with beam_mom_nn=" << qapars.beam_mom_nn
         << " not recognized or empty." << endl;
    return -1;
  }
  
  // Convert input file to tree
  // --------------------------
  if ( gSystem->AccessPathName( rootname ) ){ // quoting the ROOT documentation: "Attention, bizarre convention of return value!!"
    cout << " ======================= " << endl;
    cout << " Transforming input file " << endl
	 << qapars.txtfilename << endl
	 << " into root file " << endl
	 << rootname << endl;  
    auto buildresult = BuildTree(qapars.txtfilename.c_str(), qapars.outpath.c_str(), qapars.nevents);
    if ( buildresult !=0 ){
      cerr << "Failed to build a tree from " << qapars.txtfilename << endl;
      return buildresult;
    }
  } else {
    cout << " ======================= " << endl;
    cout << " Reusing existing root file " << rootname << endl;  
  }

  // Smear the tree
  // --------------
  TString smearedname = rootname;
  smearedname.ReplaceAll (".root",".smeared.root" );
  // Can disable warnings here.
  // Many warnings are harmless ( x or y can be smeared to values >1)
  // But it is recommended to leave it on and follow up on "inf", "nan" etc. if you test a new detector
  erhic::DisKinematics::BoundaryWarning=false;
  SmearTree( detector, rootname.Data(), smearedname.Data());  

  // -------------
  // Load the tree
  // -------------
  TChain* inTree = new TChain("EICTree");
  inTree->Add(rootname);
  inTree->AddFriend("Smeared",smearedname);
  
  // Setup Input Event Buffer
  erhic::EventMC* inEvent(NULL);
  Smear::Event* inEventS(NULL);
  inTree->SetBranchAddress("event",&inEvent);
  inTree->SetBranchAddress("eventS",&inEventS);

  // Open histo root file and book histograms
  // ----------------------------------------
  TFile * outfile = new TFile ( qapars.outfilebase + qapars.detstring + "_" + qapars.beam_mom_nn + ".root", "RECREATE");

  // A collection of event-wise qa plots, like kinematics
  // eventqacollection is defined in the header
  eventqacollection eventqa = {}; // this syntax initializes everything to 0
  initializeeventqa ( qapars, eventqa );
  
  // We'll also have a collection of histos and maybe other info for every pid
  // pidqacollection is defined in the header
  map<int,pidqacollection> qabook;
  // By default, use the standard particles
  if ( qapars.pids.size() == 0 ) qapars.pids = { 2212, 2112  }; // p, n
  // if ( qapars.pids.size() == 0 ) qapars.pids = { 2212, 2112, 22, 130, 3122 }; // p, n, gamma, K0L, Lambda0
  initializepidqabook ( qapars, qabook );
    
  // -------------
  // Analysis loop
  // -------------  
  for(long iEvent=0; iEvent<inTree->GetEntries(); iEvent++){
    
    //Read next Event
    if(inTree->GetEntry(iEvent) <=0) break;
    if(iEvent%10000 == 0) cout << "Event " << iEvent << endl;
    
    // -------------
    // event-wise QA
    // -------------
    // Following function just contains many statements of the form
    // if ( inEventS->GetQ2()>0 ) eventqa.Q2_NM->Fill ( std::log10(inEvent->GetQ2()), std::log10(inEventS->GetQ2()));
    FillEventQA( eventqa, inEvent, inEventS );

    // -------------------
    // Loop over Particles
    // -------------------
    for(int j=0; j<inEventS->GetNTracks(); j++){
      // Skip beam
      if ( j<3 ) continue;

      const Smear::ParticleMCS* inParticleS = inEventS->GetTrack(j); // Smeared Particle      
      const Particle* inParticle = inEvent->GetTrack(j); // Unsmeared Particle

      // Skip non-final particles. 
      if ( inParticle->GetStatus() != 1 ) continue;

      // Particle was not smeared
      if(inParticleS == NULL) continue; 

      // ----------------
      // Particle-wise QA
      // ----------------
      // Following function just contains many statements of the form
      // coll.dEta_p->Fill(inParticle->GetP(), inParticle->GetEta() - inParticleS->GetEta());
      FillParticleQA( qabook, inParticle, inParticleS, qapars );
    }
  }
  
  // NOTE: The remainder of this long file is tedious and explicit creation and filling of histograms
  // The Fill*QA functions demonstrate how to access values in detail, but this is the end of the 
  // basic work flow to use eic-smear from start to finish! 

  qapars.usedevents = inTree->GetEntries();
  PlotQA( qapars, eventqa, qabook );  
  outfile->Write();
  
  return 0;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

void FillEventQA( eventqacollection& eventqa , const erhic::EventMC* const inEvent, const Smear::Event* const inEventS ){

    // Q2
    // if ( eventqa.Q2_JB && inEventS->GetQ2JacquetBlondel()>0 ) eventqa.Q2_JB->Fill ( std::log10(inEvent->GetQ2()), std::log10(inEventS->GetQ2JacquetBlondel()));
    // else eventqa.missedQ2_JB++;
    // if ( eventqa.Q2_DA && inEventS->GetQ2DoubleAngle()>0 )    eventqa.Q2_DA->Fill ( std::log10(inEvent->GetQ2()), std::log10(inEventS->GetQ2DoubleAngle()));
    // else eventqa.missedQ2_DA++;
    
    // y
    // if ( eventqa.y_JB && inEventS->GetYJacquetBlondel()>0 ) eventqa.y_JB->Fill ( inEvent->GetY(), inEventS->GetYJacquetBlondel());
    // else eventqa.missedy_JB++;
    // if ( eventqa.y_DA && inEventS->GetYDoubleAngle()>0 )    eventqa.y_DA->Fill ( inEvent->GetY(), inEventS->GetYDoubleAngle());
    // else eventqa.missedy_DA++;
    
    // x
    // if ( eventqa.x_JB && inEventS->GetXJacquetBlondel()>0 ) eventqa.x_JB->Fill ( std::log10(inEvent->GetX()), std::log10(inEventS->GetXJacquetBlondel()));
    // else eventqa.missedx_JB++;
    // if ( eventqa.x_DA && inEventS->GetXDoubleAngle()>0 )    eventqa.x_DA->Fill ( std::log10(inEvent->GetX()), std::log10(inEventS->GetXDoubleAngle()));
    // else eventqa.missedx_DA++;
    
}
// ---------------------------------------------------------------
void FillParticleQA( map<int,pidqacollection>& qabook,
		     const Particle* const inParticle,
		     const Smear::ParticleMCS* const inParticleS,
		     const qaparameters& qapars){  
  // If any component is smeared, all others are either smeared or 0 (meaning "not detected").
  // Could detect the latter case (with a given accuracy):
  // const double epsilon = 1e-9;
  
  // Fill histograms
  for ( auto& pidcoll : qabook ){
    auto& pid = pidcoll.first;
    auto& coll = pidcoll.second;
    if ( pid==0 || inParticle->GetPdgCode() == pid ){
            
      auto th = inParticle->GetTheta();
      
      auto P = inParticle->GetP();
      coll.P_th->Fill(th, P);
      
      auto Pt = inParticle->GetPt();
      coll.Pt_th->Fill(th, Pt);
      
      auto xL = inParticle->GetP() / qapars.beam_mom_nn;
      coll.xL_th->Fill(th, xL);

      coll.Phi_theta->Fill(th, inParticle->GetPhi());
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
      } else if ( arg == "-i" ){
	if (++parg ==arguments.end() ){ argsokay=false; break; }
	qapars.txtfilename=*parg;
      } else if ( arg == "-N" ){      
	if (++parg==arguments.end() ){ argsokay=false; break; }
	qapars.nevents=std::stoi(parg->data());
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
      	 << " [-i txtfilename] (Lund-style file)"  << endl
	 << " [-o OutFileBase] (extension will be added)"  << endl
      	 << " [-N Nevents] (<0 for all)" << endl
      	 << " [-addpid pid] (can be called multiple times)" << endl
	 << " [-det detstring] matrix, matrixff [beam_mom_nn], handbook, perfect, beast, ephenix, zeus, jleic (capitalization does not matter.)" << endl
	 << endl;
    throw std::runtime_error("Not a valid list of options");
  }
  for (auto & c: qapars.detstring) c = toupper(c);

  return qapars;
}

// ---------------------------------------------------------------
const TString getrootname(const qaparameters& qapars ){
  // The root file name is created by replacing the extension by ".root"

  TString rootname = qapars.txtfilename;
  
  // Remove zip extension, if there is one.
  if ( rootname.EndsWith(".gz", TString::kIgnoreCase) ||
       rootname.EndsWith(".zip", TString::kIgnoreCase) )
    rootname.Replace(rootname.Last('.'), rootname.Length(), "");
  
  // Remove the remaining extension, if there is one.
  if (rootname.Last('.') > -1) {
    rootname.Replace(rootname.Last('.'), rootname.Length(), "");
  }  // if
  
  // BuildTree includes event number in partial transformation
  if ( qapars.nevents>=0 ) {
    rootname += ".";
    rootname += qapars.nevents;
    rootname += "event";
  }
  rootname += ".root";
  rootname = gSystem->BaseName( rootname );
  rootname.Prepend( qapars.outpath );
  return rootname;
}
// ---------------------------------------------------------------
void initializepidqabook(const qaparameters& qapars, map<int,pidqacollection>& qabook ){
  gStyle->SetHistLineColor(kRed); // for Profiles

  TString s, t;
  auto pdgdb=TDatabasePDG::Instance() ;
  
  // everything happens between 1e-7 and 20 mrad 
  // will try log bins 
  float thmin = 0;
  float thmax = 0.025; // some room to the side, endcap starts at 35 mrad or 60 mrad - expand here to see this
  int thbins = 100;
  int logthbins = 100;

  float pmin = 0;
  float pmax = qapars.beam_mom_nn*1.1;  // going a bit above so you can see if your beam_mom_nn is wrong
  int pbins = 100;
 
  float ptmin = 0;
  float ptmax = 2; // unlikely to exceed 2.5
  int ptbins = 50;

  float xLmin = 0;
  float xLmax = 1.1; // going a bit above so you can see if your beam_mom_nn is wrong
  int xLbins = 100;
  
  float phimin = 0;
  float phimax = TMath::TwoPi();
  int phibins = 64;

  for ( auto pid : qapars.pids ){
    // pid = abs ( pid ); 
    TString pdgname = pdgdb->GetParticle(pid)->GetName();
    t = qapars.detstring + " " + qapars.beam_mom_nn + " " + pdgname;

    s = qapars.detstring + "_P_th_"; s += pid;
    qabook[pid].P_th = new TH2D( s,t+" P vs. #theta;#theta;P", thbins, thmin, thmax, pbins, pmin, pmax);

    s = qapars.detstring + "_Pt_th_"; s += pid;
    qabook[pid].Pt_th = new TH2D( s,t+" p_{T} vs. #theta;#theta;Pt", thbins, thmin, thmax, ptbins, ptmin, ptmax);

    s = qapars.detstring + "_xL_th_"; s += pid;
    qabook[pid].xL_th = new TH2D( s,t+" x_{L} vs. #theta;#theta;xL", thbins, thmin, thmax, xLbins, xLmin, xLmax);

    s = qapars.detstring + "_phi_theta_"; s += pid;
    qabook[pid].Phi_theta = new TH2D( s,t+" #phi vs. #theta;#theta;#phi", thbins, thmin, thmax, phibins, phimin, phimax );

  }

}

// ---------------------------------------------------------------
void initializeeventqa(const qaparameters& qapars, eventqacollection& eventqa ){
  // recording log of x, y, Q2
  // keep around if we want to reactivate them
  return;
  
  // TString s;

  // float Q2min = 1e-3;
  // float Q2max = 1e5;
  // int logQ2bins = 240;

  // s = qapars.detstring + "_LogQ2_JB";
  // eventqa.Q2_JB = new TH2D( s,s+";log Q^{2};log Q^{2}_{JB}", logQ2bins, std::log10(Q2min), std::log10(Q2max), logQ2bins, std::log10(Q2min), std::log10(Q2max));

  // s = qapars.detstring + "_LogQ2_DA";
  // eventqa.Q2_DA = new TH2D( s,s+";log Q^{2};log Q^{2}_{DA}", logQ2bins, std::log10(Q2min), std::log10(Q2max), logQ2bins, std::log10(Q2min), std::log10(Q2max));

  // float ymin = 0;
  // float ymax = 1.2; // see y>1 as well
  // int ybins = 100;
  
  // s = qapars.detstring + "_y_JB";
  // eventqa.y_JB = new TH2D( s,s+";y;y_{JB}", ybins, ymin, ymax, ybins, ymin, ymax);

  // s = qapars.detstring + "_y_DA";
  // eventqa.y_DA = new TH2D( s,s+";y;y_{DA}", ybins, ymin, ymax, ybins, ymin, ymax);

  // float xmin = 1e-4;
  // float xmax = 1e1; // see x>1 as well
  // int logxbins = 200;

  // s = qapars.detstring + "_Logx_JB";
  // eventqa.x_JB = new TH2D( s,s+";log x;log x_{JB}", logxbins, std::log10(xmin), std::log10(xmax), logxbins, std::log10(xmin), std::log10(xmax));

  // s = qapars.detstring + "_Logx_DA";
  // eventqa.x_DA = new TH2D( s,s+";log x;log x_{DA}", logxbins, std::log10(xmin), std::log10(xmax), logxbins, std::log10(xmin), std::log10(xmax));

}

// ---------------------------------------------------------------
void PlotQA ( const qaparameters& qapars, eventqacollection& eventqa, map<int,pidqacollection>& qabook ){

  // Stat  position and size
  // -----------------------
  gStyle->SetStatX(0.25);
  gStyle->SetStatW(0.15);
  gStyle->SetStatY(0.9);
  gStyle->SetStatH(0.15);

  // Position of the "Missed: " box
  float missx = 0.55;
  float missy = 0.2;
  float missy2 = 0.8;
  TText t;
  t.SetNDC();

  // Suppress obnoxious flood of
  // "Current canvas added to pdf file"
  gErrorIgnoreLevel = kWarning;

  // prep a pdf collection
  new TCanvas;
  auto pdfname = qapars.outfilebase + qapars.detstring + "_" + qapars.beam_mom_nn + ".pdf";
  gPad->SaveAs( pdfname + "[" );

  // event-wise qa

  // // DA
  // if ( eventqa.y_DA ) {
  //   eventqa.y_DA->Draw("colz");
  //   t.DrawText( missx,missy, Form("Missed: %ld / %ld",eventqa.missedy_DA, qapars.usedevents));
  //     gPad->SaveAs( pdfname );
  // }
  // if ( eventqa.x_DA ) {
  //   eventqa.x_DA->Draw("colz");
  //   t.DrawText( missx,missy, Form("Missed: %ld / %ld",eventqa.missedx_DA, qapars.usedevents));
  //   gPad->SaveAs( pdfname );
  // }
  // if ( eventqa.Q2_DA ) {
  //   eventqa.Q2_DA->Draw("colz");
  //   t.DrawText( missx,missy, Form("Missed: %ld / %ld",eventqa.missedQ2_DA, qapars.usedevents));
  //   gPad->SaveAs( pdfname );
  // }
  // // JB
  // if ( eventqa.y_JB ) {
  //   eventqa.y_JB->Draw("colz");
  //   t.DrawText( missx,missy, Form("Missed: %ld / %ld",eventqa.missedy_JB, qapars.usedevents));
  //   gPad->SaveAs( pdfname );
  // }
  // if ( eventqa.x_JB ) {
  //   eventqa.x_JB->Draw("colz");
  //   t.DrawText( missx,missy, Form("Missed: %ld / %ld",eventqa.missedx_JB, qapars.usedevents));
  //   gPad->SaveAs( pdfname );
  // }
  // if ( eventqa.Q2_JB ) {
  //   eventqa.Q2_JB->Draw("colz");
  //   t.DrawText( missx,missy, Form("Missed: %ld / %ld",eventqa.missedQ2_JB, qapars.usedevents));
  //   gPad->SaveAs( pdfname );
  // }

  // particle QA
  // -----------
  // reposition stat box
  gStyle->SetStatX(0.89);
  gStyle->SetStatW(0.15);
  gStyle->SetStatY(0.89);
  gStyle->SetStatH(0.15);

  // iterating over the qabook is numerically sorted by pid
  // We'd rather use the order specified in the pids vector
  // for ( auto& pidcoll : qabook ){
  //   auto& pid = pidcoll.first;
  //   auto& coll = pidcoll.second;

  for ( auto& pid : qapars.pids ){
    auto& coll = qabook[pid];

    coll.Phi_theta->Draw("colz");
    gPad->SaveAs( pdfname );
    
    coll.P_th->Draw("colz");
    gPad->SaveAs( pdfname );

    coll.Pt_th->Draw("colz");
    gPad->SaveAs( pdfname );

    coll.xL_th->Draw("colz");
    gPad->SaveAs( pdfname );

  }

  // return to standard warning level
  gErrorIgnoreLevel = kInfo;

  // close the pdf collection
  gPad->SaveAs( pdfname +"]" );
}
// ---------------------------------------------------------------
