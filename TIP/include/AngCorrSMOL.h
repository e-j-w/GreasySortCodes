#ifndef AngCorrG_h
#define AngCorrG_h

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

TList *acList;

//Angular correlation plots
TH1F *angCorrRaw, *angCorrBG, *angCorrEvtMix;
TH1F *angCorrCoreRaw, *angCorrCoreBG, *angCorrCoreEvtMix;
TH1F *evtMixDist;

TF1 *corr420, *corr520;

class AngCorr {

	public :

		AngCorr(){;} 
		void SortData(const char*, const char*);
		void Initialise();
};
#endif

void AngCorr::Initialise() {

  cout << "Creating lists" << endl;

  acList = new TList;

  cout << "Creating histograms" << endl;

  TH1::SetDefaultSumw2(kTRUE);
  angCorrRaw = new TH1F("Raw angular correlation", "Raw angular correlation", NUM_HIST_BINS, -1, 1);
  angCorrRaw->GetXaxis()->SetTitle("cos(#theta)");
  angCorrRaw->GetYaxis()->SetTitle("Counts");
  acList->Add(angCorrRaw);
  angCorrBG = new TH1F("Background angular correlation", "Background angular correlation", NUM_HIST_BINS, -1, 1);
  angCorrBG->GetXaxis()->SetTitle("cos(#theta)");
  angCorrBG->GetYaxis()->SetTitle("Counts");
  acList->Add(angCorrBG);
  angCorrEvtMix = new TH1F("Event mixed angular correlation", "Event Mixed angular correlation", NUM_HIST_BINS, -1, 1);
  angCorrEvtMix->GetXaxis()->SetTitle("cos(#theta)");
  angCorrEvtMix->GetYaxis()->SetTitle("Counts");
  acList->Add(angCorrEvtMix);
  angCorrCoreRaw = new TH1F("Raw core angular correlation", "Raw core angular correlation", NUM_HIST_BINS, -1, 1);
  angCorrCoreRaw->GetXaxis()->SetTitle("cos(#theta)");
  angCorrCoreRaw->GetYaxis()->SetTitle("Counts");
  acList->Add(angCorrCoreRaw);
  angCorrCoreBG = new TH1F("Background core angular correlation", "Background core angular correlation", NUM_HIST_BINS, -1, 1);
  angCorrCoreBG->GetXaxis()->SetTitle("cos(#theta)");
  angCorrCoreBG->GetYaxis()->SetTitle("Counts");
  acList->Add(angCorrCoreBG);
  angCorrCoreEvtMix = new TH1F("Event mixed core angular correlation", "Event Mixed core angular correlation", NUM_HIST_BINS, -1, 1);
  angCorrCoreEvtMix->GetXaxis()->SetTitle("cos(#theta)");
  angCorrCoreEvtMix->GetYaxis()->SetTitle("Counts");
  acList->Add(angCorrCoreEvtMix);
  evtMixDist = new TH1F("Event mixing separation", "Event mixing separation", 1024, 0, 1024);
  evtMixDist->GetXaxis()->SetTitle("Event # difference");
  evtMixDist->GetYaxis()->SetTitle("Counts");
  acList->Add(evtMixDist);

  corr420 = new TF1("corr420","[0]*(1 + 0.1020408163265306*0.5*(3*x*x - 1) + 0.009070294784580489*(1/8.0)*(35*x*x*x*x - 30*x*x + 3))",-1,1);
  corr420->SetParameter(0,1.00);
  acList->Add(corr420);

  corr520 = new TF1("corr520","[0]*(1 + 0.17857142857142844*0.5*(3*x*x - 1) - 0.0043290043290043255*(1/8.0)*(35*x*x*x*x - 30*x*x + 3))",-1,1);
  corr520->SetParameter(0,1.00);
  acList->Add(corr520);
  

}
