#ifndef TTipTGamma_PIDsep_h
#define TTipTGamma_PIDsep_h

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

TH2F *tipTtigT_tigE, *tipTtigT_doppE, *EEGamma, *EEDopp;

PIDGates *gates;

class TTipTGamma_PIDsep{
	public :

		TTipTGamma_PIDsep(){;}
    void WriteData(char const *outfile);
		void SortData(char const *afile, char const *calfile, const int nP, const int nA);
		void Initialise(const int nP, const int nA);
};
#endif

