#ifndef EDopp_PIDsep_h
#define EDopp_PIDsep_h

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

TList *tigPIDSepList[MAX_NUM_PARTICLE+1][MAX_NUM_PARTICLE+1], *tipPIDSepList;

//TIGRESS PID Separated plots
TH1F *tigE_xayp[MAX_NUM_PARTICLE+1][MAX_NUM_PARTICLE+1];
TH1F *tigE_xayp_ring[MAX_NUM_PARTICLE+1][MAX_NUM_PARTICLE+1][NTIGRING];
TH1F *tigE_xayp_segring[MAX_NUM_PARTICLE+1][MAX_NUM_PARTICLE+1][NTIGSEGRING];

//TIP PID Separated plots
TH1F *tipE_xayp[MAX_NUM_PARTICLE+1][MAX_NUM_PARTICLE+1];
TH1F *tipESum_xayp[MAX_NUM_PARTICLE+1][MAX_NUM_PARTICLE+1];

PIDGates *gates;

class EDopp_PIDsep{
	public :

		EDopp_PIDsep(){;} 
		void SortData(const char*, const char*, const char*);
		void Initialise();
};
#endif

void EDopp_PIDsep::Initialise(){

  printf("Start initialization\n");
  printf("Creating lists\n");

  for(int i=0; i<MAX_NUM_PARTICLE+1; i++){
    for(int j=0; j<MAX_NUM_PARTICLE+1; j++){
      tigPIDSepList[i][j] = new TList;
    }
  }
  tipPIDSepList = new TList;

  //Setup TIP PID gates
  gates = new PIDGates;

  printf("Creating histograms\n");

  //PID separated data
  for(int i=0; i<MAX_NUM_PARTICLE+1; i++){
    for(int j=0; j<MAX_NUM_PARTICLE+1; j++){
      if((i+j)<=MAX_NUM_PARTICLE){
        //TIGRESS
        tigE_xayp[i][j] = new TH1F(Form("TIGRESS Doppler corrected addback energy (%ip%ia gate)",i,j),Form("TIGRESS Doppler corrected addback energy (%ip%ia gate)",i,j),8192,0,8192);
        tigE_xayp[i][j]->GetXaxis()->SetTitle("E_{#gamma} (keV)");
        tigE_xayp[i][j]->GetYaxis()->SetTitle("Counts");
        tigPIDSepList[i][j]->Add(tigE_xayp[i][j]);
        //TIP
        tipE_xayp[i][j] = new TH1F(Form("TIP energy (%ip%ia gate)",i,j),Form("TIP energy (%ip%ia gate)",i,j),2048,0,64);
        tipE_xayp[i][j]->GetXaxis()->SetTitle("E (MeV)");
        tipE_xayp[i][j]->GetYaxis()->SetTitle("Counts");
        tipPIDSepList->Add(tipE_xayp[i][j]);
        tipESum_xayp[i][j] = new TH1F(Form("TIP sum (%ip%ia gate)",i,j),Form("TIP sum (%ip%ia gate)",i,j),2048,0,128);
        tipESum_xayp[i][j]->GetXaxis()->SetTitle("E (MeV)");
        tipESum_xayp[i][j]->GetYaxis()->SetTitle("Counts");
        tipPIDSepList->Add(tipESum_xayp[i][j]);
      }
    }
  }
  for(int i=0; i<MAX_NUM_PARTICLE+1; i++){
    for(int j=0; j<MAX_NUM_PARTICLE+1; j++){
      if((i+j)<=MAX_NUM_PARTICLE){
        for(int k=0; k<NTIGRING; k++){
          //TIGRESS ring spectra
          tigE_xayp_ring[i][j][k] = new TH1F(Form("Ring %i TIGRESS Doppler corrected addback energy (%ip%ia gate)",k+1,i,j),Form("Ring %i TIGRESS Doppler corrected addback energy (%ip%ia gate)",k+1,i,j),8192,0,8192);
          tigE_xayp_ring[i][j][k]->GetXaxis()->SetTitle("E_{#gamma} (keV)");
          tigE_xayp_ring[i][j][k]->GetYaxis()->SetTitle("Counts");
          tigPIDSepList[i][j]->Add(tigE_xayp_ring[i][j][k]);
        }
        for(int k=0; k<NTIGSEGRING; k++){
          //TIGRESS segment ring spectra
          tigE_xayp_segring[i][j][k] = new TH1F(Form("Segment ring %i TIGRESS Doppler correctedaddback energy (%ip%ia gate)",k+1,i,j),Form("Segment ring %i TIGRESS Doppler corrected addback energy (%ip%ia gate)",k+1,i,j),8192,0,8192);
          tigE_xayp_segring[i][j][k]->GetXaxis()->SetTitle("E_{#gamma} (keV)");
          tigE_xayp_segring[i][j][k]->GetYaxis()->SetTitle("Counts");
          tigPIDSepList[i][j]->Add(tigE_xayp_segring[i][j][k]);
        }
      }
    }
  }
}

