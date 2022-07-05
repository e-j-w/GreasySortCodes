#ifndef PlotTimeSepPID_h
#define PlotTimeSepPID_h

#include <iostream>
#include <iomanip>
#include "TCutG.h"
#include "TH1.h"
#include "TF1.h"
#include "TCutG.h"
#include "TTree.h"
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

TList *tipPIDList, *tipPIDGateList; 

//TIP PID
TH2F *tip_E_PID_Sum, *tip_E_PID_Ring[NTIPRING], *tip_E_PID[NTIP];

PIDGates *gates;

class PlotTimeSepPID{
	public :

		PlotTimeSepPID(){;} 
		void SortData(const char*, const char*, const char*);
		void Initialise();
};
#endif

void PlotTimeSepPID::Initialise(){

  printf("Start initialization\n");
  printf("Creating lists\n");

  tipPIDList = new TList;
  tipPIDGateList = new TList;

  printf("Creating histograms\n");

  //TIP PID
  tip_E_PID_Sum = new TH2F("TIP energy vs PID (Sum)","TIP energy vs PID (Sum)",1024,0,128,512,0,512);
  tip_E_PID_Sum->GetYaxis()->SetTitle("A_{S}/A_{F} x 100");
  tip_E_PID_Sum->GetXaxis()->SetTitle("CsI E");
  tipPIDList->Add(tip_E_PID_Sum);
  for(int i=0; i<NTIPRING; i++){
    tip_E_PID_Ring[i] = new TH2F(Form("TIP energy vs PID (Ring %i)",i),Form("TIP energy vs PID (Ring %i)",i),1024,0,128,512,0,512);
    tip_E_PID_Ring[i]->GetYaxis()->SetTitle("A_{S}/A_{F} x 100");
    tip_E_PID_Ring[i]->GetXaxis()->SetTitle("CsI E");
    tipPIDList->Add(tip_E_PID_Ring[i]);
  }
  for(int i=0; i<NTIP; i++){
    tip_E_PID[i] = new TH2F(Form("TIP energy vs PID (Pos %i)",i+1),Form("TIP energy vs PID (Pos %i)",i+1),1024,0,128,512,0,512);
    tip_E_PID[i]->GetYaxis()->SetTitle("A_{S}/A_{F} x 100");
    tip_E_PID[i]->GetXaxis()->SetTitle("CsI E");
    tipPIDList->Add(tip_E_PID[i]);
  }

  //Setup TIP PID gates
  gates = new PIDGates;
  for(int i=0; i<NTIPRING; i++){
    tipPIDGateList->Add(gates->alphaRingCut[i]);
    tipPIDGateList->Add(gates->protonRingCut[i]);
  }

}
