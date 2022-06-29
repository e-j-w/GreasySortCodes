#ifndef EDopp_PIDsep_mca_h
#define EDopp_PIDsep_mca_h

#include <iostream>
#include <iomanip>
#include "TCutG.h"
#include "TH1.h"
#include "TF1.h"
#include "TCutG.h"
#include "TTree.h"
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

//spectra
float mcaOut[NTIGRING+1][S32K]; //output .fmca data, sp 0 is sum, sp 1-6 are TIGRESS rings
float mcaProjOut[NTIGRING+1][S32K]; //output projection

PIDGates *gates;

class EDopp_PIDsep_mca{
	public :

		EDopp_PIDsep_mca(){;} 
		void SortData(const char*, const char*, const unsigned int, const unsigned int, const int, const int, const int);
		void Initialise();
};
#endif

void EDopp_PIDsep_mca::Initialise(){

  cout << "Start initialization." << endl;

  //Setup TIP PID gates
  gates = new PIDGates;

}

