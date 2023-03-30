#ifndef TIPRingPairs_h
#define TIPRingPairs_h

#include <iostream>
#include <iomanip>
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TList.h"
#include "TTigress.h"
#include "TTip.h"
#include "TParserLibrary.h"
#include "TEnv.h"

#include "TCanvas.h"
#include "TApplication.h"

#include "evt_fmt.h"
#include <stdint.h> //allows uint8_t and similiar types

using namespace std;

class TIPRingPairs {

	public :

		TIPRingPairs(){;} 
		void SortData(const char*);
};
#endif
