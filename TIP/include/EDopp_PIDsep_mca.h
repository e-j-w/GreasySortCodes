#ifndef EDopp_PIDsep_mca_h
#define EDopp_PIDsep_mca_h

#include <iostream>
#include <iomanip>
#include "TCutG.h"
#include "TTree.h"
#include "TChain.h"
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
float mcaOut[NTIGRING+NTIGSEGRING+1][S32K]; //output .fmca data, sp 0 is sum, sp 1-6 are TIGRESS rings, sp 7-18 are segment rings

PIDGates *gates;

class EDopp_PIDsep_mca{
	public :

		EDopp_PIDsep_mca(){;}
		void WriteData(const unsigned int, const unsigned int);
		void SortData(const char*, const char*, const unsigned int, const unsigned int, double);

};

#endif

