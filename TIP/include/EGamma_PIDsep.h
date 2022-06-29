#ifndef EGamma_PIDsep_h
#define EGamma_PIDsep_h

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

TList *tigPIDSepList, *tipPIDSepList;

//TIGRESS PID Separated plots
TH1F *tigE_xayp[MAX_NUM_PARTICLE+1][MAX_NUM_PARTICLE+1];

//TIP PID Separated plots
TH1F *tipE_xayp[MAX_NUM_PARTICLE+1][MAX_NUM_PARTICLE+1];
TH1F *tipESum_xayp[MAX_NUM_PARTICLE+1][MAX_NUM_PARTICLE+1];

PIDGates *gates;

class EGamma_PIDsep{
	public :

		EGamma_PIDsep(){;} 
		void SortData(const char*, const char*, const char*);
		void Initialise();
};
#endif

void EGamma_PIDsep::Initialise(){

  printf("Start initialization\n");
  printf("Creating lists\n");

  tigPIDSepList = new TList;
  tipPIDSepList = new TList;

  //Setup TIP PID gates
  gates = new PIDGates;

  printf("Creating histograms\n");

  //PID separated data
  for(int i=0; i<MAX_NUM_PARTICLE+1; i++){
    for(int j=0; j<MAX_NUM_PARTICLE+1; j++){
      if((i+j)<=MAX_NUM_PARTICLE){
        //TIGRESS
        tigE_xayp[i][j] = new TH1F(Form("TIGRESS addback energy (%ip%ia gate)",i,j),Form("TIGRESS addback energy (%ip%ia gate)",i,j),8192,0,8192);
        tigE_xayp[i][j]->GetXaxis()->SetTitle("E_{#gamma} (keV)");
        tigE_xayp[i][j]->GetYaxis()->SetTitle("Counts");
        tigPIDSepList->Add(tigE_xayp[i][j]);
        //TIP
        tipE_xayp[i][j] = new TH1F(Form("TIP energy (%ip%ia gate)",i,j),Form("TIP energy (%ip%ia gate)",i,j),2048,0,4096);
        tipE_xayp[i][j]->GetXaxis()->SetTitle("E (MeV)");
        tipE_xayp[i][j]->GetYaxis()->SetTitle("Counts");
        tipPIDSepList->Add(tipE_xayp[i][j]);
        tipESum_xayp[i][j] = new TH1F(Form("TIP sum (%ip%ia gate)",i,j),Form("TIP sum (%ip%ia gate)",i,j),2048,0,8192);
        tipESum_xayp[i][j]->GetXaxis()->SetTitle("E (MeV)");
        tipESum_xayp[i][j]->GetYaxis()->SetTitle("Counts");
        tipPIDSepList->Add(tipESum_xayp[i][j]);
      }
    }
  }

}
