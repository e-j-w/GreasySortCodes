#ifndef CheckTimingWindows_h
#define CheckTimingWindows_h

#include <iostream>
#include <iomanip>
#include "TCutG.h"
#include "TH1.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCutG.h"
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
#include "TEmma.h"
#include "TParserLibrary.h"
#include "TEnv.h"

#include "TCanvas.h"
#include "TApplication.h"

using namespace std;

TApplication *theApp;

TList *tiptipList, *tigtigList, *tiptigList; 

//TIP
TH1F *tiptipFitT, *tiptipFitTPassed;

//TIGRESS-TIGRESS
TH1F *addT_addT, *addT_addTPassed;

//TIGRESS-TIP
TH1F *tipT_tigT_diff, *tipT_tigT_diffPassed;
TH2I *tiptig_mult, *tiptig_multSupp;
TH2I *tiptig_multPassed, *tiptig_multSuppPassed;

class CheckTimingWindows{
	public :

		CheckTimingWindows(){;} 
		void SortData(const char*, const char*, const char*);
		void Initialise();
};
#endif

void CheckTimingWindows::Initialise(){

  printf("Start initialization\n");
  printf("Creating lists\n");

  tiptipList = new TList;
  tigtigList = new TList;
  tiptigList = new TList;

  printf("Creating histograms\n");

  //TIP-TIP
  tiptipFitT = new TH1F("TIP-TIP time (fit)","TIP-TIP time (fit)",8192,-4096,4096);
  tiptipList->Add(tiptipFitT);
  tiptipFitTPassed = new TH1F("After gate: TIP-TIP time (fit)","After gate: TIP-TIP time (fit)",8192,-4096,4096);
  tiptipList->Add(tiptipFitTPassed);

  //TIGRESS-TIGRESS
  addT_addT = new TH1F("Tigress-Tigress time","Tigress-Tigress time",4096,-2048,2048); 
  tigtigList->Add(addT_addT);
  addT_addTPassed = new TH1F("After gate: Tigress-Tigress time","After gate: Tigress-Tigress time",4096,-2048,2048); 
  tigtigList->Add(addT_addTPassed);

  //TIGRESS-TIP
  tipT_tigT_diff = new TH1F("TIP fit - Tigress time","TIP fit - Tigress time",4096,-4096,4096);
  tipT_tigT_diff->GetXaxis()->SetTitle("t_{TIP, fit} - t_{TIGRESS} (ns)");
  tiptigList->Add(tipT_tigT_diff);
  tipT_tigT_diffPassed = new TH1F("After gate: TIP fit - Tigress time","After gate: TIP fit - Tigress time",4096,-4096,4096);
  tipT_tigT_diffPassed->GetXaxis()->SetTitle("t_{TIP, fit} - t_{TIGRESS} (ns)");
  tiptigList->Add(tipT_tigT_diffPassed);

  tiptig_mult = new TH2I("TIP-TIGRESS multiplicity","TIP-TIGRESS multiplicity",10,0,10,10,0,10);
  tiptig_mult->GetXaxis()->SetTitle("TIP multiplicity");
  tiptig_mult->GetYaxis()->SetTitle("TIGRESS multiplicity (unsuppressed)");
  tiptigList->Add(tiptig_mult);
  tiptig_multSupp = new TH2I("TIP-TIGRESS multiplicity (suppressed)","TIP-TIGRESS multiplicity (suppressed)",10,0,10,10,0,10);
  tiptig_multSupp->GetXaxis()->SetTitle("TIP multiplicity");
  tiptig_multSupp->GetYaxis()->SetTitle("TIGRESS multiplicity (suppressed)");
  tiptigList->Add(tiptig_multSupp);

  tiptig_multPassed = new TH2I("After gate: TIP-TIGRESS multiplicity","After gate: TIP-TIGRESS multiplicity",10,0,10,10,0,10);
  tiptig_multPassed->GetXaxis()->SetTitle("TIP multiplicity");
  tiptig_multPassed->GetYaxis()->SetTitle("TIGRESS multiplicity (unsuppressed)");
  tiptigList->Add(tiptig_multPassed);
  tiptig_multSuppPassed = new TH2I("After gate: TIP-TIGRESS multiplicity (suppressed)","After gate: TIP-TIGRESS multiplicity (suppressed)",10,0,10,10,0,10);
  tiptig_multSuppPassed->GetXaxis()->SetTitle("TIP multiplicity");
  tiptig_multSuppPassed->GetYaxis()->SetTitle("TIGRESS multiplicity (suppressed)");
  tiptigList->Add(tiptig_multSuppPassed);

}
