#ifndef TTGamma_PIDsep_h
#define TTGamma_PIDsep_h

#include <iostream>
#include <iomanip>
#include "TCutG.h"
#include "TH1.h"
#include "TF1.h"
#include "TCutG.h"
#include "TTree.h"
#include "TEntryList.h"
#include "TChain.h"
#include "TH2.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TList.h"
#include "TTigress.h"
#include "TTip.h"
#include "TChannel.h"
#include "TParserLibrary.h"
#include "TEnv.h"

#include "TCanvas.h"
#include "TApplication.h"

using namespace std;


TApplication *theApp;

TList *tttimeList;

//TIGRESS PID Separated plots
TH1F *tigTT;

PIDGates *gates;

class TTGamma_PIDsep{
	public :

		TTGamma_PIDsep(){;} 
		void SortData(char const *afile, char const *calfile, const int nP, const int nA, const float Elow_start, const float Ehigh_start, const float Elow_stop, const float Ehigh_stop, char const *outfile);
		void Initialise();
};
#endif

void TTGamma_PIDsep::Initialise(){

  printf("Start initialization\n");
  printf("Creating lists\n");

  tttimeList = new TList;

  //Setup TIP PID gates
  gates = new PIDGates;
}

