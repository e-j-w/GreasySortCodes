#ifndef EEDopp_PIDsep_h
#define EEDopp_PIDsep_h

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

TList *tigPIDSepList;

//TIGRESS PID Separated plots
TH2F *tigEE_xayp[MAX_NUM_PARTICLE+1][MAX_NUM_PARTICLE+1];

PIDGates *gates;

class EEDopp_PIDsep{
	public :

		EEDopp_PIDsep(){;} 
		void SortData(const char*, const char*, const char*);
		void Initialise();
};
#endif

void EEDopp_PIDsep::Initialise(){

  printf("Start initialization\n");
  printf("Creating lists\n");

  tigPIDSepList = new TList;

  //Setup TIP PID gates
  gates = new PIDGates;

  printf("Creating histograms\n");

  //PID separated data
  for(int i=0; i<MAX_NUM_PARTICLE+1; i++){
    for(int j=0; j<MAX_NUM_PARTICLE+1; j++){
      if((i+j)<=MAX_NUM_PARTICLE){
        //TIGRESS
        tigEE_xayp[i][j] = new TH2F(Form("TIGRESS Doppler corrected addback energy vs. addback energy (%ip%ia gate)",i,j),Form("TIGRESS Doppler corrected addback energy vs. addback energy (%ip%ia gate)",i,j),8192,0,8192,8192,0,8192);
        tigEE_xayp[i][j]->GetXaxis()->SetTitle("E_{#gamma 1} (keV)");
        tigEE_xayp[i][j]->GetYaxis()->SetTitle("E_{#gamma 2} (keV)");
        tigPIDSepList->Add(tigEE_xayp[i][j]);
      }
    }
  }
}

