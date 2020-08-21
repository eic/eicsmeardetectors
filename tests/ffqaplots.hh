#ifndef FFQAPLOTS_H
#define FFQAPLOTS_H

// Note: Nothing in here is necessary for eic-smear usage

#include <TString.h>
#include <TH1.h>
#include <TH2.h>

#include <string>
#include <vector>

struct qaparameters{
  std::string txtfilename="./ep_hiQ2.20x250.small.txt.gz";
  TString outfilebase="./ffqaplots";
  std::string outpath="./";
  long nevents=-1;
  std::vector<int> pids = {}; // sign will be ignored. 0 for all. leave empty for n, p. 
  std::string detstring = "MATRIXFF"; // Capitalization does not matter

  long usedevents=-1;// pure convenience so I can access the true number when nevents=-1;

  // ion beam momentum per nucleon in GeV for far forward detectors.
  // Using int to avoid rounding issues in switch
  int beam_mom_nn=100; 
};

struct eventqacollection {
  // // could re-add these to check improved performance 
  // // by more accepted hadrons
  // TH2D* Q2_JB;
  // TH2D* Q2_DA;
  // TH2D* y_JB;
  // TH2D* y_DA;
  // TH2D* x_JB;
  // TH2D* x_DA;

  // long missedQ2_JB;
  // long missedQ2_DA;
  // long missedy_JB;
  // long missedy_DA;
  // long missedx_JB;
  // long missedx_DA;
};

struct pidqacollection {
  TH2D* Phi_theta;
  TH2D* P_th;
  TH2D* Pt_th;
  TH2D* xL_th;
};



qaparameters ParseArguments ( int argc, char* argv[] );

#endif // FFQAPLOTS_H
