#ifndef ReprocessGain_SMOL_h
#define ReprocessGain_SMOL_h

#include <iostream>
#include <iomanip>
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TList.h"
#include "TRandom.h"
#include "TTigress.h"
#include "TGriffin.h"
#include "TGriffinBgo.h"
#include "TChannel.h"
#include "TParserLibrary.h"
#include "TEnv.h"

#include "TCanvas.h"
#include "TApplication.h"

using namespace std;

#define MAX_TIME_WINDOWS_PER_TREE  1024
#define MAX_INPUT_E                8

TApplication *theApp;

class ReprocessGain_SMOL{
	public :

		ReprocessGain_SMOL(){;} 
		void SortData(const char *sfile, const char *outfile, const double evalWindowSize, const double tWindowSize, const double en[MAX_INPUT_E], const int numEnVals);
};
#endif
