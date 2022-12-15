#ifndef SeparatorTiming_h
#define SeparatorTiming_h

#include <iostream>
#include <iomanip>
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TList.h"
#include "TRandom.h"
#include "TTigress.h"
#include "TTip.h"
#include "TChannel.h"
#include "TParserLibrary.h"
#include "TEnv.h"

#include "TCanvas.h"
#include "TApplication.h"

using namespace std;

TApplication *theApp;

class SeparatorTiming{
	public :

		SeparatorTiming(){;} 
		void SortData(const char*, const uint8_t, const uint8_t, const char*, const char*);
};
#endif
