#ifndef SeparatorTiming1Tip2Tig_h
#define SeparatorTiming1Tip2Tig_h

#include <iostream>
#include <iomanip>
#include "TTree.h"
#include "TChain.h"
#include "TH2.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TList.h"
#include "TRandom.h"
#include "TReaction.h"
#include "TSRIM.h"
#include "TTigress.h"
#include "TTip.h"
#include "TSpectrum.h"
#include "TChannel.h"
#include "TParserLibrary.h"
#include "TEnv.h"

#include "TCanvas.h"
#include "TApplication.h"

using namespace std;

TApplication *theApp;

class SeparatorTiming1Tip2Tig{
	public :

		SeparatorTiming1Tip2Tig(){;} 
		void SortData(const char*, const char*, const char*);
};
#endif
