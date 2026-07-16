#ifndef AngCorrS_h
#define AngCorrS_h

#include <iostream>
#include <iomanip>
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TList.h"

#include "TCanvas.h"
#include "TApplication.h"

#include "evt_fmt.h"
#include <stdint.h> //allows uint8_t and similiar types

using namespace std;

TApplication *theApp;

TList *angcorrList; 

//angular correlation plots
TH1D *angCorrRaw, *angCorrCompton, *angCorrEvtMix, *angCorrEvtMixCompton, *angCorrEvtMixCorr, *angCorrEvtMixCorrCompton;

#define MIN_EVT_MIX_OFFSET 200
#define MAX_EVT_MIX_OFFSET 2000

#endif

void InitialiseHists() {

  cout << "Creating lists" << endl;

  angcorrList = new TList;

  cout << "Creating histograms" << endl;
  //histogram names shouldn't use spaces, to aid GRSISort-based analysis

  //Raw HPGe Spectra
  angCorrRaw = new TH1D("angCorrRaw", "Raw angular correlation", 8192, -1, 1);
  angCorrRaw->GetXaxis()->SetTitle("cos #theta");
  angCorrRaw->GetYaxis()->SetTitle("Normalized counts");
  angcorrList->Add(angCorrRaw);
  angCorrCompton = new TH1D("angCorrCompton", "Compton BG angular correlation", 8192, -1, 1);
  angCorrCompton->GetXaxis()->SetTitle("cos #theta");
  angCorrCompton->GetYaxis()->SetTitle("Normalized counts");
  angcorrList->Add(angCorrCompton);
  angCorrEvtMix = new TH1D("angCorrEvtMix", "Event mixing angular correlation", 8192, -1, 1);
  angCorrEvtMix->GetXaxis()->SetTitle("cos #theta");
  angCorrEvtMix->GetYaxis()->SetTitle("Normalized counts");
  angcorrList->Add(angCorrEvtMix);
  angCorrEvtMixCompton = new TH1D("angCorrEvtMixCompton", "Compton BG event mixing angular correlation", 8192, -1, 1);
  angCorrEvtMixCompton->GetXaxis()->SetTitle("cos #theta");
  angCorrEvtMixCompton->GetYaxis()->SetTitle("Normalized counts");
  angcorrList->Add(angCorrEvtMixCompton);
  angCorrEvtMixCorr = new TH1D("angCorrEvtMixCorr", "Event mixing corrected angular correlation", 8192, -1, 1);
  angCorrEvtMixCorr->GetXaxis()->SetTitle("cos #theta");
  angCorrEvtMixCorr->GetYaxis()->SetTitle("Normalized counts");
  angcorrList->Add(angCorrEvtMixCorr);
  angCorrEvtMixCorrCompton = new TH1D("angCorrEvtMixCorrCompton", "Event mixing and Compton corrected angular correlation", 8192, -1, 1);
  angCorrEvtMixCorrCompton->GetXaxis()->SetTitle("cos #theta");
  angCorrEvtMixCorrCompton->GetYaxis()->SetTitle("Normalized counts");
  angcorrList->Add(angCorrEvtMixCorrCompton);

}
