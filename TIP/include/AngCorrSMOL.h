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
#define NUM_HIST_BINS 256

using namespace std;

TApplication *theApp;

TList *acList;

//Angular correlation plots
TH1F *angCorrRaw, *angCorrEvtMix;
TH1F *angCorrCoreRaw, *angCorrCoreEvtMix;
TH1F *evtMixDist, *corePosAngCorr, *corePosEvtMix;
TH2F *corePosVsAngCorr, *corePosVsEvtMix;

TF1 *corrfit, *corrfitzeroa4;

PIDGates *gates;

class AngCorr {

	public :

		AngCorr(){;} 
		void SortData(const char*, const char*, const int);
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
  angCorrEvtMix = new TH1F("Event mixed angular correlation", "Event Mixed angular correlation", NUM_HIST_BINS, -1, 1);
  angCorrEvtMix->GetXaxis()->SetTitle("cos(#theta)");
  angCorrEvtMix->GetYaxis()->SetTitle("Counts");
  acList->Add(angCorrEvtMix);
  angCorrCoreRaw = new TH1F("Raw core angular correlation", "Raw core angular correlation", NUM_HIST_BINS, -1, 1);
  angCorrCoreRaw->GetXaxis()->SetTitle("cos(#theta)");
  angCorrCoreRaw->GetYaxis()->SetTitle("Counts");
  acList->Add(angCorrCoreRaw);
  angCorrCoreEvtMix = new TH1F("Event mixed core angular correlation", "Event Mixed core angular correlation", NUM_HIST_BINS, -1, 1);
  angCorrCoreEvtMix->GetXaxis()->SetTitle("cos(#theta)");
  angCorrCoreEvtMix->GetYaxis()->SetTitle("Counts");
  acList->Add(angCorrCoreEvtMix);
  evtMixDist = new TH1F("Event mixing separation", "Event mixing separation", 1024, 0, 1024);
  evtMixDist->GetXaxis()->SetTitle("Event # difference");
  evtMixDist->GetYaxis()->SetTitle("Counts");
  acList->Add(evtMixDist);
  corePosAngCorr = new TH1F("Raw angular correlation position (hit 2)", "Raw angular correlation position (hit 2)", 64, 0, 64);
  corePosAngCorr->GetXaxis()->SetTitle("TIGRESS position");
  corePosAngCorr->GetYaxis()->SetTitle("Counts");
  acList->Add(corePosAngCorr);
  corePosEvtMix = new TH1F("Event mixing position (hit 2)", "Event mixing position (hit 2)", 64, 0, 64);
  corePosEvtMix->GetXaxis()->SetTitle("TIGRESS position");
  corePosEvtMix->GetYaxis()->SetTitle("Counts");
  acList->Add(corePosEvtMix);
  corePosVsAngCorr = new TH2F("Raw angular correlation vs. Raw angular correlation position (hit 2)", "Raw angular correlation vs. Raw angular correlation position (hit 2)", NUM_HIST_BINS, -1, 1, 64, 0, 64);
  corePosVsAngCorr->GetXaxis()->SetTitle("cos(#theta)");
  corePosVsAngCorr->GetYaxis()->SetTitle("TIGRESS position");
  acList->Add(corePosVsAngCorr);
  corePosVsEvtMix = new TH2F("Event mixed core angular correlation vs. Event mixing position (hit 2)", "Event mixed core angular correlation vs. Event mixing position (hit 2)", NUM_HIST_BINS, -1, 1, 64, 0, 64);
  corePosVsEvtMix->GetXaxis()->SetTitle("cos(#theta)");
  corePosVsEvtMix->GetYaxis()->SetTitle("TIGRESS position");
  acList->Add(corePosVsEvtMix);

  corrfit = new TF1("corrfit","[0]*(1 + [1]*0.5*(3*x*x - 1) + [2]*(1/8.0)*(35*x*x*x*x - 30*x*x + 3))",-1,1);
  corrfit->SetParameter(0,1.00);
  corrfit->SetParameter(1,0.00);
  corrfit->SetParameter(2,0.00);
  acList->Add(corrfit);

  corrfitzeroa4 = new TF1("corrfitzeroa4","[0]*(1 + [1]*0.5*(3*x*x - 1))",-1,1);
  corrfitzeroa4->SetParameter(0,1.00);
  corrfitzeroa4->SetParameter(1,0.00);
  acList->Add(corrfitzeroa4);
  

}
