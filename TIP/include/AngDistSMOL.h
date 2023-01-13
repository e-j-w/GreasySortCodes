#ifndef AngDistG_h
#define AngDistG_h

#include <iostream>
#include <iomanip>
#include "TH1.h"
#include "TH2.h"
#include "TGraphErrors.h"
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

#define NUM_PAST_HIT 128
#define NUM_HIST_BINS 4096

using namespace std;

TApplication *theApp;

TList *adList;

//Angular correlation plots
TH1F *AngDistRaw, *AngDistBG, *AngDistSource;
TH1F *AngDistCoreRaw, *AngDistCoreBG, *AngDistCoreSource;

//TF1 *corr420, *corr520;

class AngDist {

	public :

		AngDist(){;} 
		void SortData(const char*, const char*, const char*);
		void Initialise();
};
#endif

void AngDist::Initialise() {

  cout << "Creating lists" << endl;

  adList = new TList;

  cout << "Creating histograms" << endl;

  TH1::SetDefaultSumw2(kTRUE);
  AngDistRaw = new TH1F("Raw angular distribution", "Raw angular distribution", NUM_HIST_BINS, -1, 1);
  AngDistRaw->GetXaxis()->SetTitle("cos(#theta)");
  AngDistRaw->GetYaxis()->SetTitle("Counts");
  adList->Add(AngDistRaw);
  AngDistBG = new TH1F("Background angular distribution", "Background angular distribution", NUM_HIST_BINS, -1, 1);
  AngDistBG->GetXaxis()->SetTitle("cos(#theta)");
  AngDistBG->GetYaxis()->SetTitle("Counts");
  adList->Add(AngDistBG);
  AngDistSource = new TH1F("Source angular distribution", "Source angular distribution", NUM_HIST_BINS, -1, 1);
  AngDistSource->GetXaxis()->SetTitle("cos(#theta)");
  AngDistSource->GetYaxis()->SetTitle("Counts");
  adList->Add(AngDistSource);
  AngDistCoreRaw = new TH1F("Raw core angular distribution", "Raw core angular distribution", NUM_HIST_BINS, -1, 1);
  AngDistCoreRaw->GetXaxis()->SetTitle("cos(#theta)");
  AngDistCoreRaw->GetYaxis()->SetTitle("Counts");
  adList->Add(AngDistCoreRaw);
  AngDistCoreBG = new TH1F("Background core angular distribution", "Background core angular distribution", NUM_HIST_BINS, -1, 1);
  AngDistCoreBG->GetXaxis()->SetTitle("cos(#theta)");
  AngDistCoreBG->GetYaxis()->SetTitle("Counts");
  adList->Add(AngDistCoreBG);
  AngDistCoreSource = new TH1F("Source core angular distribution", "Source core angular distribution", NUM_HIST_BINS, -1, 1);
  AngDistCoreSource->GetXaxis()->SetTitle("cos(#theta)");
  AngDistCoreSource->GetYaxis()->SetTitle("Counts");
  adList->Add(AngDistCoreSource);

  /*corr420 = new TF1("corr420","[0]*(1 + 0.1020408163265306*0.5*(3*x*x - 1) + 0.009070294784580489*(1/8.0)*(35*x*x*x*x - 30*x*x + 3))",-1,0.6);
  corr420->SetParameter(0,1.00);
  adList->Add(corr420);

  corr520 = new TF1("corr520","[0]*(1 + 0.17857142857142844*0.5*(3*x*x - 1) - 0.0043290043290043255*(1/8.0)*(35*x*x*x*x - 30*x*x + 3))",-1,0.6);
  corr520->SetParameter(0,1.00);
  adList->Add(corr520);*/
  

}
