#ifndef ReprocessAB_SMOL_h
#define ReprocessAB_SMOL_h

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

TApplication *theApp;

class ReprocessAB_SMOL{
	public :

		ReprocessAB_SMOL(){;} 
		uint64_t SortData();
};
#endif
